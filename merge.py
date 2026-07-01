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
                  NO_PLASMID_OVERLAP       (sources share no plasmid token (main_plasmids or transfection
                                            at all — different experiments; no override)
                  MULTIPLE_MAIN_PLASMIDS   (overlapping sources need confirmed selected Main_Plasmids;
                                            resolution: main_plasmids = [selected tokens] —
                                            see Main_Plasmids canonicalization)
    NEEDS_INPUT   FILENAME_EXCLUSION_CONFLICT  (same file, sources disagree on excluded
                                                wells; resolution: dup_keep::<File_Name>
                                                = <src label>, default keep-first)
    INFO/AUTO     SCHEMA_MIGRATED, NEW_CONDITIONS, EXTRA_COLUMNS_DROPPED, FILENAME_DUPLICATE
                  (identical copies deduped), FILENAME_COLLISION_RENAMED,
                  CONDITION_PARTITION_OVERLAP (canonicalization-unresolved conditions dropped
                  as no-metadata).

Main_Plasmids canonicalization (plan_canonicalization / apply_canonicalization) reconciles
the cosmetic Main_Plasmids/Transfection split:
  * Main_Plasmids MATCH across sources -> nothing to reconcile: no prompt, merge as-is.
  * Main_Plasmids DIFFER -> judge overlap over each source's FULL token set
    (Main_Plasmids ∪ Transfection, role-agnostic):
      - share >=1 token -> the user MUST pick a global Main_Plasmids backbone
        (resolutions["main_plasmids"]).
      - share NO token  -> NO_PLASMID_OVERLAP (forbidden, different experiments).
The chosen tokens become the global Main_Plasmids; Transfection is re-derived per condition
as the remaining tokens. Conditions that do not contain selected plasmids are dropped as no-metadata
(CONDITION_PARTITION_OVERLAP). Exclusion state (Is_Excluded) is preserved; the unified blob is
re-scoped + backstopped so it re-resolves to the identical excluded-well set.

merge_masters runs the analysis; if any FORBIDDEN issue is present and its matching
override/resolution is not supplied, it raises MergeForbidden (with the report attached).
Otherwise it builds and returns the merged frame.
"""

import hashlib
import logging
from dataclasses import dataclass, field

import pandas as pd

from models import MASTER_COLUMNS, LEGACY_COLUMN_DEFAULTS, APP_VERSION, MAIN_ONLY_CONDITION
from export import ensure_master_csv_schema
from processing import coerce_bool
from exclusions import (scope_blob_for_merge, remap_blob_file_names,
                        well_token, build_resolve_ctx, parse_exclusion_blob)

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
CONDITION_PARTITION_OVERLAP = "CONDITION_PARTITION_OVERLAP"   # emitted for unresolved canon conditions
MULTIPLE_MAIN_PLASMIDS = "MULTIPLE_MAIN_PLASMIDS"            # forbidden unless main_plasmids resolved
NO_PLASMID_OVERLAP = "NO_PLASMID_OVERLAP"                    # sources share no plasmid token (no override)
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
    # Main_Plasmids canonicalization plan (set by _analyze)
    canon_plan: "CanonPlan | None" = None

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
# Main_Plasmids canonicalization
# --------------------------------------------------------------------------- #
# The Main_Plasmids / Transfection split is cosmetic. Canonicalization lets the user
# pick a global Main_Plasmids backbone and re-derives Transfection per condition.
# Trigger: sources disagree on Main_Plasmids
# Tokenization: split on " + ", strip, drop empties. Case-fold for identity/ordering;
# preserve original casing in emitted strings.

# Unresolved conditions (those that cannot honor the chosen backbone) are relabelled into
# this Transfection sink so the existing "Empty" drop/index/N logic removes them.
_CANON_DROP_SINK = "Empty (no selected plasmid)"


@dataclass
class CanonPlan:
    """Plan describing a Main_Plasmids canonicalization (see plan_canonicalization).

    requires_selection : the sources' Main_Plasmids DIFFER and they still overlap (share >=1
                         plasmid token across Main_Plasmids ∪ Transfection) -> the user must
                         pick plasmid(s) for the Main_Plasmids (the GUI shows the selection dialog).
                         False when the backbones already match (nothing to reconcile).
    candidate_tokens   : [(token, count), ...] every plasmid token across all sources,
                         ordered by descending MEASUREMENT count (distinct File_Name carrying
                         the token), then alphabetically (case-insensitive). Drives the dialog.
    preselected        : default plasmid(s) for the Main_Plasmids = the single top
                         candidate (most measurements, then alphabetical).
    selected_main      : the chosen tokens (None until provided), used by apply_canonicalization.
    rewrites           : concrete per-condition rewrites (only when selected_main given):
                         [{match: (cell, ligand, sorted_tokens),
                           from: (main, transf), to: (main, transf), affected_files: [...]}].
    unresolved         : conditions whose token set does NOT contain all selected-main tokens
                         (dropped as no-metadata): [{cell, ligand, from, tokens, affected_files}].
    """
    requires_selection: bool = False
    candidate_tokens: list = field(default_factory=list)
    preselected: list = field(default_factory=list)
    selected_main: "list | None" = None
    rewrites: list = field(default_factory=list)
    unresolved: list = field(default_factory=list)


def _tokens(value) -> list:
    """Split a Main_Plasmids / Transfection string into plasmid tokens (' + ' separated,
    stripped, empties dropped). Casing is preserved here; callers case-fold for identity."""
    if value is None:
        return []
    return [t.strip() for t in str(value).split(" + ") if t.strip()]


def _join_transfection(remaining_cf, casing) -> str:
    """Join the remaining (case-folded) tokens into a canonical Transfection string.

    The main plasmid-only sentinel (MAIN_ONLY_CONDITION, "-") is meaningful ONLY when it
    stands alone — once canonicalization leaves the well with a construct, drop the sentinel"""
    toks = sorted(remaining_cf)
    sentinel = MAIN_ONLY_CONDITION.casefold()
    if sentinel in toks and len(toks) > 1:
        toks = [t for t in toks if t != sentinel]
    return " + ".join(casing.get(k, k) for k in toks)


def _build_casing(selected_main, *value_iterables) -> dict:
    """Map case-folded token -> canonical original casing (first occurrence wins). The
    selected_main casing is registered first so the emitted backbone uses the user's casing."""
    casing = {}

    def reg(tok):
        k = tok.casefold()
        if k not in casing:
            casing[k] = tok

    for t in (selected_main or []):
        s = str(t).strip()
        if s:
            reg(s)
    for it in value_iterables:
        for v in it:
            for tok in _tokens(v):
                reg(tok)
    return casing


def _full_tokens(df) -> set:
    """Case-folded set of every plasmid token in a source — the union over all rows of
    tokens(Main_Plasmids) ∪ tokens(Transfection). Role-agnostic: a token counts whether it
    appears as a Main or a per-condition (Transfection) plasmid."""
    if df is None or df.empty:
        return set()
    out = set()
    for col in ('Main_Plasmids', 'Transfection'):
        if col in df.columns:
            for v in df[col].astype(str).unique():
                for t in _tokens(v):
                    out.add(t.casefold())
    return out


def _declared_main_tokens(df) -> set:
    """Case-folded token set of a source's declared Main_Plasmids column only."""
    if df is None or df.empty or 'Main_Plasmids' not in df.columns:
        return set()
    out = set()
    for v in df['Main_Plasmids'].astype(str).unique():
        for t in _tokens(v):
            out.add(t.casefold())
    return out


def _mains_match(frames) -> bool:
    """True if every source declares the SAME Main_Plasmids (case-/order-insensitive
    token-set equality). When the Main_Plasmids already agree there is nothing to reconcile,
    so the merge neither prompts nor is forbidden."""
    return len({frozenset(_declared_main_tokens(df)) for df in frames}) <= 1


def _overlap_status(token_sets):
    """Cross-source plasmid-token overlap status over full (Main ∪ Transfection) token sets.

      None  -> fewer than 2 non-empty sources (overlap not applicable).
      True  -> every source shares >=1 token with the union of the others (connected) ->
               the sources are related; the merge must pick ONE Main_Plasmids backbone.
      False -> at least one source is token-disjoint from all others -> different
               experiments; the merge is forbidden (NO_PLASMID_OVERLAP, no override).
    """
    ne = [s for s in token_sets if s]
    if len(ne) < 2:
        return None
    for i, s in enumerate(ne):
        others = set().union(*[o for j, o in enumerate(ne) if j != i])
        if s.isdisjoint(others):
            return False
    return True


def _as_source_frames(obj) -> list:
    """Normalize plan_canonicalization input to a list of DataFrames. Accepts a single
    DataFrame, a list of DataFrames, or a list of MergeSource."""
    if isinstance(obj, pd.DataFrame):
        return [obj]
    frames = []
    for it in obj:
        if isinstance(it, pd.DataFrame):
            frames.append(it)
        elif hasattr(it, "df"):
            frames.append(it.df)
        else:
            raise TypeError("plan_canonicalization expects a DataFrame, a list of "
                            "DataFrames, or a list of MergeSource.")
    return frames


def plan_canonicalization(sources_or_merged_df, selected_main=None) -> CanonPlan:
    """Analyse sources for a Main_Plasmids mismatch and (optionally) plan the rewrite.
    Pure analysis — never mutates the input.

    requires_selection is True when sources disagree on the Main_Plasmids backbone AND that
    disagreement is reconcilable: either a true content collision (same Cell_Line/Ligand/
    full-token-set under different Main strings) OR the declared Main Plasmids overlap (share a
    token), e.g. "b2AR-nLuc" vs "b2AR-nLuc + miniG". Genuinely disjoint sources (no shared
    token) do not require a selection and merge untouched.

      * No reconcilable mismatch -> requires_selection=False, no rewrites (no-op).
      * Mismatch + selected_main is None -> requires_selection=True with candidate_tokens
        (descending occurrence) and preselected (declared-main intersection); NO rewrites.
      * selected_main provided -> concrete rewrites + unresolved using that backbone.
    """
    frames = _as_source_frames(sources_or_merged_df)

    casing = {}

    def reg(tok):
        k = tok.casefold()
        if k not in casing:
            casing[k] = tok

    for t in (selected_main or []):
        s = str(t).strip()
        if s:
            reg(s)

    declared_main_sets = []          # per source: set of case-folded declared-main tokens
    cond_agg = {}                    # (main, transf, cell, ligand) -> {cf:set, files:set}
    token_to_files = {}              # cf token -> set(File_Name), for candidate ordering/count

    full_token_sets = []             # per source: case-folded Main_Plasmids ∪ Transfection

    for df in frames:
        if df is None or df.empty or 'Main_Plasmids' not in df.columns:
            declared_main_sets.append(set())
            full_token_sets.append(set())
            continue
        cols = [c for c in ('Main_Plasmids', 'Transfection', 'Cell_Line', 'Ligand', 'File_Name')
                if c in df.columns]
        sub = df[cols].astype(str)
        dmain = set()
        for mv in sub['Main_Plasmids'].unique():
            for t in _tokens(mv):
                reg(t)
                dmain.add(t.casefold())
        declared_main_sets.append(dmain)
        full_token_sets.append(_full_tokens(df))

        grp = (sub.groupby(['Main_Plasmids', 'Transfection', 'Cell_Line', 'Ligand'])['File_Name']
               .apply(lambda s: set(s)).reset_index())
        for _, r in grp.iterrows():
            main, transf = r['Main_Plasmids'], r['Transfection']
            cell, lig = r['Cell_Line'], r['Ligand']
            toks = _tokens(main) + _tokens(transf)
            for t in toks:
                reg(t)
            cf = frozenset(t.casefold() for t in toks)
            agg = cond_agg.setdefault((main, transf, cell, lig), {'cf': set(cf), 'files': set()})
            agg['files'] |= r['File_Name']
            for t in cf:
                token_to_files.setdefault(t, set()).update(r['File_Name'])

    # Candidate tokens ordered by descending MEASUREMENT count, then alphabetically
    ordered = sorted(token_to_files.keys(),
                     key=lambda k: (-len(token_to_files[k]), k))
    candidate_tokens = [(casing.get(k, k), len(token_to_files[k])) for k in ordered]

    # Default = the top candidate
    preselected = [candidate_tokens[0][0]] if candidate_tokens else []

    # Trigger (per the user-confirmed definition):
    #   * Main_Plasmids already MATCH across sources -> nothing to reconcile: NO prompt.
    #   * Main_Plasmids DIFFER -> look for any shared plasmid token across the FULL
    #     Main_Plasmids ∪ Transfection set (role-agnostic: "b2AR" as Main in one source and
    #     Transfection in another counts):
    #         - overlap -> PROMPT the user to choose the backbone (requires_selection).
    #         - no overlap -> different experiments -> forbidden (NO_PLASMID_OVERLAP, emitted
    #           in _analyze).
    mains_match = len({frozenset(s) for s in declared_main_sets}) <= 1
    requires_selection = (not mains_match) and (_overlap_status(full_token_sets) is True)

    rewrites, unresolved = [], []
    sel_cf = {str(t).casefold() for t in selected_main if str(t).strip()} \
        if selected_main is not None else None
    if sel_cf:
        global_main = " + ".join(casing.get(k, k) for k in sorted(sel_cf))
        for (main, transf, cell, lig), agg in cond_agg.items():
            cf = agg['cf']
            if sel_cf <= cf:
                new_transf = _join_transfection(cf - sel_cf, casing)
                if (global_main, new_transf) != (main, transf):
                    rewrites.append({
                        "match": (cell, lig, tuple(sorted(cf))),
                        "from": (main, transf),
                        "to": (global_main, new_transf),
                        "affected_files": sorted(agg['files'])})
            else:
                unresolved.append({
                    "cell": cell, "ligand": lig,
                    "from": (main, transf),
                    "tokens": sorted(cf),
                    "affected_files": sorted(agg['files'])})

    return CanonPlan(
        requires_selection=requires_selection,
        candidate_tokens=candidate_tokens,
        preselected=preselected,
        selected_main=list(selected_main) if selected_main is not None else None,
        rewrites=rewrites,
        unresolved=unresolved)


def apply_canonicalization(df: pd.DataFrame, plan: CanonPlan, log_fn=None) -> pd.DataFrame:
    """Apply a CanonPlan to a master-shaped frame, rewriting ONLY Main_Plasmids and
    Transfection. Resolved conditions get the global Main_Plasmids + the remaining tokens
    (Transfection); unresolved conditions get the _CANON_DROP_SINK Transfection so downstream
    "Empty" handling drops them. No recompute — every numeric / identity / exclusion column
    is preserved. A no-op (returns an equal copy) when the plan carries no Main_plasmids."""
    if df is None or df.empty or plan is None or not getattr(plan, "selected_main", None):
        return df.copy() if df is not None else df

    sel_cf = {str(t).casefold() for t in plan.selected_main if str(t).strip()}
    if not sel_cf:
        return df.copy()

    casing = _build_casing(plan.selected_main,
                           df['Main_Plasmids'].astype(str).unique(),
                           df['Transfection'].astype(str).unique())
    global_main = " + ".join(casing.get(k, k) for k in sorted(sel_cf))

    out = df.copy()
    combo_map = {}
    for _, r in out[['Main_Plasmids', 'Transfection']].astype(str).drop_duplicates().iterrows():
        main, transf = r['Main_Plasmids'], r['Transfection']
        cf = {t.casefold() for t in (_tokens(main) + _tokens(transf))}
        if sel_cf <= cf:
            combo_map[(main, transf)] = (global_main, _join_transfection(cf - sel_cf, casing))
        else:
            combo_map[(main, transf)] = (None, _CANON_DROP_SINK)   # unresolved -> sink

    keys = list(zip(out['Main_Plasmids'].astype(str), out['Transfection'].astype(str)))
    new_main, new_transf = [], []
    for (m, t) in keys:
        nm, nt = combo_map[(m, t)]
        new_main.append(m if nm is None else nm)   # unresolved keeps its original Main
        new_transf.append(nt)
    out['Main_Plasmids'] = new_main
    out['Transfection'] = new_transf

    if log_fn and plan.rewrites:
        log_fn(f"   [MERGE] Canonicalized Main_Plasmids to '{global_main}' "
               f"({len(plan.rewrites)} condition rewrite(s)).")
    return out


def _canonical_conditions(df: pd.DataFrame, selected_main) -> set:
    """Distinct (Main_Plasmids, Transfection, Cell_Line, Ligand) tuples AFTER canonicalizing
    with selected_main (unresolved conditions excluded). With no selected_main, falls back to
    the raw _conditions(df) — so the NEW_CONDITIONS count is computed on the canonical view."""
    if df is None or df.empty:
        return set()
    if not selected_main:
        return _conditions(df)

    sel_cf = {str(t).casefold() for t in selected_main if str(t).strip()}
    casing = _build_casing(selected_main,
                           df['Main_Plasmids'].astype(str).unique(),
                           df['Transfection'].astype(str).unique())
    global_main = " + ".join(casing.get(k, k) for k in sorted(sel_cf))

    out = set()
    sub = df[['Main_Plasmids', 'Transfection', 'Cell_Line', 'Ligand']].astype(str).drop_duplicates()
    for _, r in sub.iterrows():
        main, transf = r['Main_Plasmids'], r['Transfection']
        cf = {t.casefold() for t in (_tokens(main) + _tokens(transf))}
        if sel_cf <= cf:
            remaining = sorted(cf - sel_cf)
            new_transf = " + ".join(casing.get(k, k) for k in remaining)
            out.add((global_main, new_transf, r['Cell_Line'], r['Ligand']))
        # unresolved conditions are dropped -> excluded from the count
    return out


# --------------------------------------------------------------------------- #
# Core analysis (shared logic behind classify_sources and merge_masters)
# --------------------------------------------------------------------------- #
def _detect_condition_partition_overlap(canon_plan: "CanonPlan | None",
                                        report: "MergeReport") -> None:
    """Emit CONDITION_PARTITION_OVERLAP (info) for canonicalization-unresolved conditions —
    conditions whose token set lacks the selected Main_Plasmids, which are dropped
    as no-metadata (never merged). No-op when there is nothing unresolved."""
    if canon_plan is None or not canon_plan.unresolved:
        return None
    conds = [{"cell": u["cell"], "ligand": u["ligand"],
              "main_plasmids": u["from"][0], "transfection": u["from"][1],
              "files": u["affected_files"]}
             for u in canon_plan.unresolved]
    report.add(Issue(
        CONDITION_PARTITION_OVERLAP, INFO,
        f"{len(conds)} condition(s) lack the selected Main_Plasmids and were "
        f"dropped as no-metadata (not merged).",
        {"count": len(conds), "conditions": conds}))
    return None


def _analyze(sources: list, selected_main=None):
    """Normalize all sources and compute the full issue report + the facts the merge
    builder needs. Does not mutate the input DataFrames and does not concatenate.

    selected_main: the chosen Main_Plasmids (from resolutions["main_plasmids"]). When the
    sources share no plasmid token, a forbidden NO_PLASMID_OVERLAP is emitted. When
    they overlap and selected_main is None, a forbidden MULTIPLE_MAIN_PLASMIDS is emitted (the
    merge is blocked until the user picks a backbone); when provided, the canonicalization
    rewrites are planned and NEW_CONDITIONS / the unresolved-condition hook are computed on the
    canonical view.

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

    # --- Main_Plasmids canonicalization plan ---
    canon_plan = plan_canonicalization([n["df"] for n in norm], selected_main=selected_main)

    # --- NO_PLASMID_OVERLAP (forbidden, NO override) ---
    if not _mains_match([n["df"] for n in norm]):
        overlap = _overlap_status([_full_tokens(n["df"]) for n in norm])
        if overlap is False:
            per_source = {n["label"]: sorted(_full_tokens(n["df"])) for n in norm}
            report.add(Issue(
                NO_PLASMID_OVERLAP, FORBIDDEN,
                "Sources have different Main_Plasmids and share no common plasmid (no "
                "overlapping token across Main_Plasmids / Transfection). They look like "
                "different experiments and cannot be merged.",
                {"tokens_by_source": per_source}))

    if canon_plan.requires_selection and not selected_main:
        # Overlapping sources, no Main_Plasmids chosen yet -> block until the user picks one
        distinct = {}
        for n in norm:
            if 'Main_Plasmids' not in n["df"].columns:
                continue
            for mv in pd.unique(n["df"]['Main_Plasmids'].dropna()):
                distinct.setdefault(str(mv), set()).add(n["label"])
        report.add(Issue(
            MULTIPLE_MAIN_PLASMIDS, FORBIDDEN,
            "Choose a single Main_Plasmids backbone for the merged data (the sources overlap "
            "but you must confirm which plasmid tokens are the shared backbone).",
            {"distinct_main_plasmids": {k: sorted(v) for k, v in distinct.items()},
             "candidate_tokens": canon_plan.candidate_tokens,
             "preselected": canon_plan.preselected,
             "resolution_key": "main_plasmids"}))

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

    # --- NEW_CONDITIONS (info) — computed on the CANONICALIZED view ---
    # Use the chosen backbone (or, pre-selection, the provisional preselected one) so cosmetic
    # split differences no longer inflate the count, and dropped/unresolved rows are excluded.
    eff_selected = selected_main if selected_main else (
        canon_plan.preselected if canon_plan.requires_selection else None)
    base_conditions = _canonical_conditions(norm[0]["df"], eff_selected) if norm else set()
    later_conditions = set()
    for n in norm[1:]:
        later_conditions |= _canonical_conditions(n["df"], eff_selected)
    new_conditions = later_conditions - base_conditions
    if new_conditions:
        report.add(Issue(
            NEW_CONDITIONS, INFO,
            f"{len(new_conditions)} new condition(s) detected.",
            {"count": len(new_conditions),
             "conditions": sorted(new_conditions)}))

    # --- CONDITION_PARTITION_OVERLAP (info) — unresolved canonicalization conditions ---
    _detect_condition_partition_overlap(canon_plan, report)

    report.canon_plan = canon_plan

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
        if code == MULTIPLE_MAIN_PLASMIDS and resolutions.get("main_plasmids"):
            continue   # resolved by the chosen backbone (defensive: _analyze omits it then)
        # LABELING_MISMATCH, TIME_VECTOR_DIVERGENCE and RAW_CHANNELS_REQUIRED have no override
        out.append(iss)
    return out


def merge_masters(sources: list, resolutions: dict | None = None, log_fn=None):
    """Merge two or more master DataFrames into one.

    Runs the shared analysis (`_analyze`, passing resolutions["main_plasmids"] as the backbone);
    if any FORBIDDEN issue lacks its override/resolution, raises MergeForbidden (report attached).
    Otherwise concatenates the normalized, collision-resolved, deduped sources, applies the
    Main_Plasmids canonicalization (dropping unresolved conditions as no-metadata), re-scopes every
    source's exclusion rules + backstops so the unified Applied_Exclusions blob re-resolves to the
    identical excluded-well set, stamps NCollector_version, and returns (merged_df, report).
    """
    resolutions = resolutions or {}
    selected_main = resolutions.get("main_plasmids")
    norm, report, facts = _analyze(sources, selected_main=selected_main)

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

    # --- Canonicalize Main_Plasmids / Transfection BEFORE re-scoping exclusion blobs ---
    # Is_Excluded is untouched by canonicalization, so per-row exclusion truth survives.
    # Doing this BEFORE scope_blob_for_merge means the blob is re-scoped against the
    # canonical columns; the coverage backstop below then guarantees every excluded well
    # still re-resolves (rules whose stale Cond:/Main: no longer match are backfilled).
    canon_dropped = 0
    if selected_main and report.canon_plan is not None and report.canon_plan.selected_main:
        merged = apply_canonicalization(merged, report.canon_plan, log_fn=log_fn)
        before = len(merged)
        merged = merged[~merged['Transfection'].astype(str)
                        .str.contains("Empty", na=False)].copy()
        canon_dropped = before - len(merged)
        if canon_dropped:
            _log(f"   [MERGE] Canonicalization dropped {canon_dropped} row(s) in "
                 f"unresolved condition(s) (no selected backbone).")

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

    # --- Coverage backstop: the unified blob MUST re-resolve to every excluded well ---
    # Is_Excluded is authoritative and untouched, but canonicalization can strand a manual rule,
    # so scope_blob_for_merge would drop it and silently lose those exclusions.
    # Any currently-excluded well not covered  by the scoped tokens gets a per-well File-pinned token
    # built from the LIVE canonical row (carrying canonical Cell/Cond/Date).
    if 'Is_Excluded' in merged.columns and 'Well_ID' in merged.columns:
        ctx = build_resolve_ctx(merged)
        covered = set()
        for t in uniq:
            for e in parse_exclusion_blob(t):
                covered |= ctx.rule_wells(e, only_excluded=True)
        excl_mask = merged['Is_Excluded'].map(coerce_bool)
        all_excl = set(zip(merged.loc[excl_mask, 'File_Name'].astype(str),
                           merged.loc[excl_mask, 'Well_ID'].astype(str)))
        for (f, w) in sorted(all_excl - covered):
            tok = well_token(merged, f, w)
            if tok not in seen:
                seen.add(tok)
                uniq.append(tok)

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


def compute_merge_preview(sources, report, main_plasmids_choice=None, log_fn=None):
    """Produce the merged frame for the summary preview WITHOUT touching app state.

    Mirrors merge_run's conflict handling so the summary reflects the ACTUAL merged
    result: hard-forbidden (no-override) issues -> None (cannot preview); filename
    collisions are auto-resolved by 'rename' for the preview (so both copies show);
    conflicting exclusions are unioned (Is_Excluded only — enough for accurate N, no
    recompute needed for counts). main_plasmids_choice is the user's chosen backbone
    (app.merge_main_plasmids_choice) used when canonicalization requires a selection."""
    log = log_fn or (lambda *a, **k: None)

    hard = [i for i in report.forbidden
            if i["code"] not in (FILENAME_DATA_COLLISION,
                                 MULTIPLE_MAIN_PLASMIDS)]
    if hard:
        codes = ", ".join(sorted({i["code"] for i in hard}))
        log(f"[MERGE] Cannot preview merged summary — unresolved forbidden "
            f"issue(s): {codes}. Resolve these before merging.")
        return None

    resolutions = {}
    if any(i["code"] == FILENAME_DATA_COLLISION for i in report.forbidden):
        resolutions["filename_collision"] = "rename"
        log("[MERGE] Summary preview: same-named files with different data are "
            "shown as separate (renamed) copies.")

    # MULTIPLE_MAIN_PLASMIDS: use selected main plasmids for display summary
    canon = getattr(report, "canon_plan", None)
    if canon is not None and canon.requires_selection:
        selected_mp = main_plasmids_choice or list(canon.preselected)
        if not selected_mp:
            log("[MERGE] Cannot preview merged summary — no Main Plasmids chosen.")
            return None
        resolutions["main_plasmids"] = list(selected_mp)

    try:
        merged_df, _r = merge_masters(sources, resolutions, log_fn=None)
    except MergeForbidden as e:
        codes = ", ".join(sorted({i["code"] for i in e.issues}))
        log(f"[MERGE] Cannot preview merged summary: {codes}.")
        return None

    # Union conflicting exclusions (Is_Excluded only) so fully-excluded columns drop
    # from N exactly as they will after a real merge.
    file_col = merged_df['File_Name'].astype(str)
    well_col = merged_df['Well_ID'].astype(str)
    for iss in report.needs_input:
        if iss["code"] != FILENAME_EXCLUSION_CONFLICT:
            continue
        ctx = iss["context"]
        fname = ctx.get("file_name")
        wells = {str(w) for lst in ctx.get("excluded_wells_by_source", {}).values()
                 for w in lst}
        for well in wells:
            m = (file_col == str(fname)) & (well_col == str(well))
            if m.any():
                merged_df.loc[m, 'Is_Excluded'] = True
    return merged_df