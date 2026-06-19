"""
merge.py — GUI-agnostic backend for merging two or more master DataFrames into one.

A "master" is the flat one-row-per-well-per-timepoint frame described by
models.MASTER_COLUMNS. Merging clean masters is, by design, NOT a metric recompute
(replicate means are within-file), so a merge is
concatenation + schema alignment + identity dedup + per-source exclusion-rule
re-scoping. This module is Tkinter-free and unit-testable with plain DataFrames.

Public API
----------
    classify_sources(sources)                   -> MergeReport            (dry run, no mutation)
    merge_masters(sources, resolutions, log_fn) -> (merged_df, MergeReport)

Both entry points run the same internal analysis (`_analyze`): classify_sources returns
just its report (for a merge-preview pane), while merge_masters runs the analysis, refuses
to proceed if any FORBIDDEN issue is unresolved, and otherwise builds the merged frame.

Safety taxonomy:

    FORBIDDEN     LABELING_MISMATCH        (no override)
                  TIME_VECTOR_DIVERGENCE   (no override)
                  FILENAME_DATA_COLLISION  (resolution: filename_collision = rename|cancel)
                  RAW_CHANNELS_REQUIRED    (no override)
    NEEDS_INPUT   FILENAME_EXCLUSION_CONFLICT  (same file, sources disagree on excluded
                                                wells; resolution: dup_keep::<File_Name>
                                                = <src label>, default keep-first)
    INFO/AUTO     SCHEMA_MIGRATED, NEW_CONDITIONS, EXTRA_COLUMNS_DROPPED, FILENAME_DUPLICATE
                  (identical copies deduped), FILENAME_COLLISION_RENAMED.

CONDITION_PARTITION_OVERLAP is reserved for a planned condition-canonicalization pass and is
currently emitted by no path (see _detect_condition_partition_overlap).

merge_masters runs the analysis; if any FORBIDDEN issue is present and its matching
override/resolution is not supplied, it raises MergeForbidden (with the report attached).
Otherwise it builds and returns the merged frame.
"""

import hashlib
import logging
from dataclasses import dataclass, field

import pandas as pd

from models import MASTER_COLUMNS, LEGACY_COLUMN_DEFAULTS, APP_VERSION
from export import ensure_master_csv_schema
from processing import coerce_bool
from restore import scope_blob_for_merge, remap_blob_file_names

logger = logging.getLogger("NCollector")

# --- Severities --------------------------------------------------------------
FORBIDDEN = "forbidden"
NEEDS_INPUT = "needs_input"
INFO = "info"

# --- Issue codes (the safety taxonomy) ---------------------------------------
LABELING_MISMATCH = "LABELING_MISMATCH"
TIME_VECTOR_DIVERGENCE = "TIME_VECTOR_DIVERGENCE"
FILENAME_DATA_COLLISION = "FILENAME_DATA_COLLISION"
RAW_CHANNELS_REQUIRED = "RAW_CHANNELS_REQUIRED"
FILENAME_DUPLICATE = "FILENAME_DUPLICATE"
FILENAME_EXCLUSION_CONFLICT = "FILENAME_EXCLUSION_CONFLICT"
CONDITION_PARTITION_OVERLAP = "CONDITION_PARTITION_OVERLAP"   # planned feature
SCHEMA_MIGRATED = "SCHEMA_MIGRATED"
NEW_CONDITIONS = "NEW_CONDITIONS"
EXTRA_COLUMNS_DROPPED = "EXTRA_COLUMNS_DROPPED"
FILENAME_COLLISION_RENAMED = "FILENAME_COLLISION_RENAMED"


# --------------------------------------------------------------------------- #
# Data carriers
# --------------------------------------------------------------------------- #
@dataclass
class MergeSource:
    """One input to the merge. `df` is a master-shaped DataFrame (a CSV master or a
    freshly compiled folder). `kind` is "master" or "folder" (informational; folder
    sources are enriched by construction). `label` must be unique across sources — it is
    used in collision suffixes and dup-resolution keys."""
    label: str
    df: pd.DataFrame
    kind: str = "master"
    path: str = ""


def Issue(code: str, severity: str, message: str, context: dict | None = None) -> dict:
    """An issue record. severity ∈ {forbidden, needs_input, info}."""
    return {"code": code, "severity": severity, "message": message,
            "context": context or {}}


@dataclass
class MergeReport:
    forbidden: list = field(default_factory=list)
    needs_input: list = field(default_factory=list)
    auto: list = field(default_factory=list)
    summary: dict = field(default_factory=dict)

    def add(self, issue: dict):
        sev = issue["severity"]
        if sev == FORBIDDEN:
            self.forbidden.append(issue)
        elif sev == NEEDS_INPUT:
            self.needs_input.append(issue)
        else:
            self.auto.append(issue)

    def all_issues(self) -> list:
        return self.forbidden + self.needs_input + self.auto


class MergeForbidden(Exception):
    """Custom Error. Raised by merge_masters when an unresolved FORBIDDEN issue
    blocks the merge. Carries the full report (and the specific unresolved issues)."""
    def __init__(self, report: MergeReport, issues: list | None = None):
        self.report = report
        self.issues = issues if issues is not None else list(report.forbidden)
        codes = ", ".join(sorted({i["code"] for i in self.issues}))
        super().__init__(f"Merge forbidden: unresolved issue(s) [{codes}]")


# --------------------------------------------------------------------------- #
# Step 0 — normalization
# --------------------------------------------------------------------------- #
def _normalize(src: MergeSource):
    """Migrate one source to the current schema, coerce Is_Excluded, and reorder to
    MASTER_COLUMNS (adding missing columns at their LEGACY default / NaN, dropping unknown
    extras). Returns (normalized_df, was_modified, dropped_extras, blob)."""
    df = src.df.copy()
    df, was_modified, _ = ensure_master_csv_schema(df)

    if 'Is_Excluded' in df.columns:
        df['Is_Excluded'] = df['Is_Excluded'].map(coerce_bool)

    extras = [c for c in df.columns if c not in MASTER_COLUMNS]
    for col in MASTER_COLUMNS:
        if col not in df.columns:
            df[col] = LEGACY_COLUMN_DEFAULTS.get(col, {}).get("default", float('nan'))
    df = df[MASTER_COLUMNS].copy()

    blob = "None"
    if 'Applied_Exclusions' in df.columns and df['Applied_Exclusions'].notna().any():
        blob = df['Applied_Exclusions'].dropna().iloc[0]
    return df, bool(was_modified), extras, blob


# --------------------------------------------------------------------------- #
# Per-source fact extraction
# --------------------------------------------------------------------------- #
def _has_labeling(df: pd.DataFrame) -> bool:
    return bool((df['Replicate'].astype(str) == "labeling control").any())


def _canonical_time_vector(df: pd.DataFrame) -> tuple:
    """sorted unique Time_(min) rounded to 6 dp (kills float noise)."""
    t = pd.to_numeric(df['Time_(min)'], errors='coerce').dropna().round(6).unique()
    return tuple(sorted(float(x) for x in t))


def _conditions(df: pd.DataFrame) -> set:
    """Distinct (Main_Plasmids, Transfection, Cell_Line, Ligand) tuples."""
    cols = ['Main_Plasmids', 'Transfection', 'Cell_Line', 'Ligand']
    sub = df[cols].astype(str).drop_duplicates()
    return set(map(tuple, sub.itertuples(index=False, name=None)))


def _file_identity_key(file_rows: pd.DataFrame) -> tuple:
    """Content-identity key for one (source, File_Name): (sorted unique Date,
    sorted unique Well_ID, sha1 of sorted (Well_ID, Time_(min)@6dp, Raw_BRET_unexcluded@9dp)
    triples). Two occurrences with the SAME key carry the same underlying measurement
    (same wells, same timepoints, same raw BRET) and are therefore the same file. The key
    deliberately ignores exclusion state and derived columns, so two copies of a file that
    differ only in which wells are switched off still share an identity key (and are routed
    to the duplicate / exclusion-conflict path, not the data-collision path)."""
    dates = tuple(sorted(file_rows['Date'].astype(str).unique()))
    wells = tuple(sorted(file_rows['Well_ID'].astype(str).unique()))
    t = pd.to_numeric(file_rows['Time_(min)'], errors='coerce').round(6)
    rb = pd.to_numeric(file_rows['Raw_BRET_unexcluded'], errors='coerce').round(9)
    triples = sorted(zip(
        file_rows['Well_ID'].astype(str).tolist(),
        [None if pd.isna(x) else float(x) for x in t.tolist()],
        [None if pd.isna(x) else float(x) for x in rb.tolist()],
    ))
    # Build identity key hash for less memory usage than comparing tuples
    h = hashlib.sha1(repr(triples).encode('utf-8')).hexdigest()
    return (dates, wells, h)


def _excluded_wells(file_rows: pd.DataFrame) -> frozenset:
    if 'Is_Excluded' not in file_rows.columns:
        return frozenset()
    mask = file_rows['Is_Excluded'].map(coerce_bool)
    return frozenset(file_rows.loc[mask, 'Well_ID'].astype(str).tolist())


# --------------------------------------------------------------------------- #
# Core analysis (shared logic behind classify_sources and merge_masters)
# --------------------------------------------------------------------------- #
def _detect_condition_partition_overlap(norm: list, report: "MergeReport") -> None:
    """HOOK for planned condition canonicalization."""
    return None


def _analyze(sources: list):
    """Normalize all sources and compute the full issue report + the facts the merge
    builder needs. Does not mutate the input DataFrames and does not concatenate.

    Returns (norm, report, facts) where:
      norm  = [{label, df, kind, path, was_modified, extras, blob, files{name->rows}}...]
      facts = {
        'labeling_ok': bool,
        'time_ok': bool,
        'collisions': {File_Name: [src_idx,...]},   # differing identity keys (forbidden)
        'duplicates': {File_Name: {'keeper': idx, 'others': [idx,...], 'agree': bool}},
        'new_condition_count': int,
      }
    """
    if not sources:
        raise ValueError("merge requires at least one source")

    report = MergeReport()
    norm = []
    for src in sources:
        df, was_mod, extras, blob = _normalize(src)
        norm.append({
            "label": src.label, "df": df, "kind": src.kind, "path": src.path,
            "was_modified": was_mod, "extras": extras, "blob": blob,
        })

    # --- SCHEMA_MIGRATED / EXTRA_COLUMNS_DROPPED (info) ---
    for n in norm:
        if n["was_modified"]:
            report.add(Issue(SCHEMA_MIGRATED, INFO,
                             f"Source '{n['label']}' was migrated to the current master schema.",
                             {"source": n["label"]}))
        if n["extras"]:
            report.add(Issue(EXTRA_COLUMNS_DROPPED, INFO,
                             f"Source '{n['label']}' had non-schema columns dropped: "
                             f"{', '.join(n['extras'])}.",
                             {"source": n["label"], "columns": list(n["extras"])}))

    # --- LABELING_MISMATCH (forbidden, NO override) ---
    labeling_flags = {n["label"]: _has_labeling(n["df"]) for n in norm}
    labeling_ok = len(set(labeling_flags.values())) <= 1
    if not labeling_ok:
        report.add(Issue(
            LABELING_MISMATCH, FORBIDDEN,
            "Sources disagree on labeling-control presence. Labeling-corrected data and "
            "non-labeling data cannot be merged.",
            {"per_source": labeling_flags}))

    # --- TIME_VECTOR_DIVERGENCE (forbidden, NO override) ---
    vectors = {n["label"]: _canonical_time_vector(n["df"]) for n in norm}
    distinct_vectors = set(vectors.values())
    time_ok = len(distinct_vectors) <= 1
    if not time_ok:
        ref_label, ref_vec = next(iter(vectors.items()))
        divergences = []
        for lab, vec in vectors.items():
            if vec == ref_vec:
                continue
            first_div = None
            for a, b in zip(ref_vec, vec):
                if a != b:
                    first_div = {"ref": a, "other": b}
                    break
            if first_div is None and len(vec) != len(ref_vec):
                idx = min(len(vec), len(ref_vec))
                first_div = {"ref": ref_vec[idx] if idx < len(ref_vec) else None,
                             "other": vec[idx] if idx < len(vec) else None}
            divergences.append({"source": lab, "length": len(vec),
                                "first_divergent": first_div})
        report.add(Issue(
            TIME_VECTOR_DIVERGENCE, FORBIDDEN,
            "Sources do not share a common Time_(min) vector. Overlaid kinetic plotting "
            "needs a common x-axis; this cannot be merged.",
            {"reference": {"source": ref_label, "length": len(ref_vec)},
             "divergent": divergences}))

    # --- File-level content-identity keys across sources ---
    # name -> list of (src_idx, identity_key, excluded_wells)
    occurrences = {}
    for i, n in enumerate(norm):
        files = {}
        for fname, fr in n["df"].groupby('File_Name', sort=False):
            files[str(fname)] = fr
            occurrences.setdefault(str(fname), []).append(
                (i, _file_identity_key(fr), _excluded_wells(fr)))
        n["files"] = files

    collisions = {}   # File_Name -> [src_idx,...] of differing-identity occurrences
    duplicates = {}   # File_Name -> {keeper, others, agree}
    for fname, occs in occurrences.items():
        if len(occs) < 2:
            continue
        keeper_idx, keeper_key, keeper_excl = occs[0]
        dup_others, dup_excl_agree = [], True
        coll_idxs = []
        for (idx, ikey, excl) in occs[1:]:
            if ikey == keeper_key:
                dup_others.append(idx)
                if excl != keeper_excl:
                    dup_excl_agree = False
            else:
                coll_idxs.append(idx)

        if dup_others:
            duplicates[fname] = {"keeper": keeper_idx, "others": dup_others,
                                 "agree": dup_excl_agree}
            dup_labels = [norm[keeper_idx]["label"]] + [norm[j]["label"] for j in dup_others]
            if dup_excl_agree:
                report.add(Issue(
                    FILENAME_DUPLICATE, INFO,
                    f"File '{fname}' appears in {len(dup_others) + 1} sources with identical "
                    f"data and exclusion state; duplicates will be dropped (keep first).",
                    {"file_name": fname, "sources": dup_labels}))
            else:
                conflict = {}
                for (idx, ikey, excl) in occs:
                    if ikey == keeper_key:
                        conflict[norm[idx]["label"]] = sorted(excl)
                report.add(Issue(
                    FILENAME_EXCLUSION_CONFLICT, NEEDS_INPUT,
                    f"File '{fname}' is the same measurement in multiple sources but the "
                    f"sources disagree on which wells are excluded. Pick which source's "
                    f"exclusion set the merged file should keep via resolution key "
                    f"'dup_keep::{fname}' = <source label> (default: keep the first source).",
                    {"file_name": fname, "resolution_key": f"dup_keep::{fname}",
                     "excluded_wells_by_source": conflict}))

        if coll_idxs:
            coll_labels = [norm[keeper_idx]["label"]] + [norm[j]["label"] for j in coll_idxs]
            collisions[fname] = [keeper_idx] + coll_idxs
            report.add(Issue(
                FILENAME_DATA_COLLISION, FORBIDDEN,
                f"File '{fname}' carries DIFFERENT raw data across sources ({', '.join(coll_labels)}). "
                f" Resolve via filename_collision = 'rename' or 'cancel'.",
                {"file_name": fname, "sources": coll_labels,
                 "resolution_key": "filename_collision"}))

    # --- RAW_CHANNELS_REQUIRED (forbidden, NO override) ---
    for n in norm:
        df = n["df"]
        offending = {}
        for col in ('Donor_Raw_kinetic', 'Acceptor_Raw_kinetic'):
            nan_mask = pd.to_numeric(df[col], errors='coerce').isna()
            if nan_mask.any():
                for fname, fr in df.loc[nan_mask].groupby('File_Name', sort=False):
                    rec = offending.setdefault(str(fname), {"nan_rows": 0, "wells": set()})
                    rec["nan_rows"] += len(fr)
                    rec["wells"].update(fr['Well_ID'].astype(str).tolist())
        if offending:
            ctx_files = {f: {"nan_rows": v["nan_rows"],
                             "sample_wells": sorted(v["wells"])[:5]}
                         for f, v in offending.items()}
            report.add(Issue(
                RAW_CHANNELS_REQUIRED, FORBIDDEN,
                f"Source '{n['label']}' has missing/non-numeric raw channel values "
                f"(Donor_Raw_kinetic / Acceptor_Raw_kinetic). Only fully enriched masters "
                f"may be merged — enrich this master first via the enrich path before merging.",
                {"source": n["label"], "files": ctx_files}))

    # --- NEW_CONDITIONS (info) ---
    base_conditions = _conditions(norm[0]["df"]) if norm else set()
    later_conditions = set()
    for n in norm[1:]:
        later_conditions |= _conditions(n["df"])
    new_conditions = later_conditions - base_conditions
    if new_conditions:
        report.add(Issue(
            NEW_CONDITIONS, INFO,
            f"{len(new_conditions)} new condition(s) detected.",
            {"count": len(new_conditions),
             "conditions": sorted(new_conditions)}))

    # --- CONDITION_PARTITION_OVERLAP (needs_input) — HOOK ---
    _detect_condition_partition_overlap(norm, report)

    facts = {
        "labeling_ok": labeling_ok,
        "time_ok": time_ok,
        "collisions": collisions,
        "duplicates": duplicates,
        "new_condition_count": len(new_conditions),
    }

    # Projected summary (dry-run view).
    report.summary = {
        "total_sources": len(sources),
        "total_files": len({f for f in occurrences}),
        "total_rows": int(sum(len(n["df"]) for n in norm)),
        "conditions_added": len(new_conditions),
        "collisions_handled": sum(len(v) - 1 for v in collisions.values()),
        "duplicates_dropped": sum(len(v["others"]) for v in duplicates.values()),
    }
    return norm, report, facts


# --------------------------------------------------------------------------- #
# Public API
# --------------------------------------------------------------------------- #
def classify_sources(sources: list) -> MergeReport:
    """Dry run: normalize each source through ensure_master_csv_schema, run every check,
    and return the report WITHOUT mutating inputs or concatenating. Use this to populate a
    merge preview pane."""
    _, report, _ = _analyze(sources)
    return report


def _unresolved_forbidden(report: MergeReport, resolutions: dict) -> list:
    """Return the FORBIDDEN issues whose override/resolution is NOT supplied."""
    out = []
    for iss in report.forbidden:
        code = iss["code"]
        if code == FILENAME_DATA_COLLISION and resolutions.get("filename_collision") == "rename":
            continue
        # LABELING_MISMATCH, TIME_VECTOR_DIVERGENCE and RAW_CHANNELS_REQUIRED have no override
        out.append(iss)
    return out


def merge_masters(sources: list, resolutions: dict | None = None, log_fn=None):
    """Merge two or more master DataFrames into one.

    Runs the shared analysis (`_analyze`); if any FORBIDDEN issue lacks its override/resolution,
    raises MergeForbidden (report attached). Otherwise concatenates the normalized, collision-resolved,
    deduped sources, re-scopes every source's exclusion rules so each re-resolves to exactly its origin
    well-set, writes the unified Applied_Exclusions blob to all rows, stamps NCollector_version, and
    returns (merged_df, report).
    """
    resolutions = resolutions or {}
    norm, report, facts = _analyze(sources)

    blocking = _unresolved_forbidden(report, resolutions)
    if blocking:
        raise MergeForbidden(report, blocking)

    def _log(msg):
        if log_fn:
            log_fn(msg)

    # --- Build working copies tagged by source index ---
    work = []
    for i, n in enumerate(norm):
        d = n["df"].copy()
        d["_src_idx"] = i
        work.append(d)

    # --- Apply collision renames (only reached if resolved to 'rename') ---
    rename_maps = {i: {} for i in range(len(norm))}   # src_idx -> {old: new}
    renamed_count = 0
    if facts["collisions"] and resolutions.get("filename_collision") == "rename":
        for fname, idxs in facts["collisions"].items():
            keeper = idxs[0]      # first occurrence keeps the original name
            for idx in idxs[1:]:
                new_name = f"{fname}#{norm[idx]['label']}"
                m = work[idx]['File_Name'].astype(str) == fname
                work[idx].loc[m, 'File_Name'] = new_name
                rename_maps[idx][fname] = new_name
                renamed_count += 1
                report.add(Issue(
                    FILENAME_COLLISION_RENAMED, INFO,
                    f"Collision on '{fname}': source '{norm[idx]['label']}' occurrence "
                    f"renamed to '{new_name}' to keep both files distinct.",
                    {"old": fname, "new": new_name, "source": norm[idx]["label"]}))
                _log(f"   [MERGE] Renamed colliding file '{fname}' -> '{new_name}'.")

    # --- Apply duplicate drops ---
    dropped_count = 0
    for fname, info in facts["duplicates"].items():
        keeper = info["keeper"]
        if not info["agree"]:
            # Disagreement: pick the source whose exclusions win (default = first/keeper).
            chosen_label = resolutions.get(f"dup_keep::{fname}")
            if chosen_label is not None:
                for j, n in enumerate(norm):
                    if n["label"] == chosen_label and (j == keeper or j in info["others"]):
                        keeper = j
                        break
            _log(f"   [MERGE] Duplicate '{fname}' exclusion conflict; keeping "
                 f"'{norm[keeper]['label']}'.")
        drop_idxs = [j for j in ([info["keeper"]] + info["others"]) if j != keeper]
        for idx in drop_idxs:
            m = work[idx]['File_Name'].astype(str) == fname
            dropped_count += int(work[idx][m]['File_Name'].nunique() > 0)
            work[idx] = work[idx].loc[~m]

    # --- Concatenate surviving rows ---
    merged = pd.concat(work, ignore_index=True)
    merged = merged.reindex(columns=MASTER_COLUMNS + ["_src_idx"])

    # --- Re-scope every source's exclusion rules against the provisional merged frame ---
    all_tokens = []
    for i, n in enumerate(norm):
        origin_mask = merged["_src_idx"] == i
        if not origin_mask.any():
            continue  # entire source deduped away
        remapped_blob = remap_blob_file_names(n["blob"], rename_maps.get(i, {}))
        tokens = scope_blob_for_merge(merged, remapped_blob, origin_mask, log_fn=log_fn)
        all_tokens.extend(tokens)

    # De-duplicate tokens preserving order, then join.
    seen, uniq = set(), []
    for t in all_tokens:
        if t not in seen:
            seen.add(t)
            uniq.append(t)
    unified_blob = " || ".join(uniq) if uniq else "None"

    # --- Finalize: drop temp tag, write blob, stamp version, reorder ---
    merged = merged.drop(columns=["_src_idx"])
    merged['Applied_Exclusions'] = unified_blob
    merged['NCollector_version'] = APP_VERSION
    merged = merged.reindex(columns=MASTER_COLUMNS)

    # Actual post-merge summary.
    report.summary.update({
        "total_sources": len(sources),
        "total_files": int(merged['File_Name'].nunique()),
        "total_rows": int(len(merged)),
        "conditions_added": facts["new_condition_count"],
        "collisions_handled": renamed_count,
        "duplicates_dropped": dropped_count,
    })

    _log(f"   [MERGE] Done: {report.summary['total_files']} files, "
         f"{report.summary['total_rows']} rows, "
         f"{len(uniq)} exclusion rule(s).")
    return merged, report