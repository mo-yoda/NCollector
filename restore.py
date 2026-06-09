"""
restore.py — GUI-agnostic backend for REVERSIBLE exclusions ("restore").

Exclusion in N Collector is per-row authoritative via Is_Excluded; the engine
(recompute_master_after_exclusion) derives every downstream column from the PRISTINE
raw (Raw_BRET_unexcluded -> Acceptor/Donor -> NaN) with the currently-excluded wells
NaN-d. Exclude and restore are therefore the SAME operation: change Is_Excluded, then
recompute. This module exposes the pure logic that decides WHICH wells a restore touches
and keeps the human-readable Applied_Exclusions blob in sync.

There is NO per-well "which rule excluded this" record and we do not add one. Which wells
a rule affects is INFERRED from the blob + the master:

  * Rules are separated by " || " (v2.0.5+). The blob is robustly re-tokenized by START
    MARKERS ("Ligand:" for manual rules, "AUTO:" for auto rules), which also recovers
    rules from legacy " | "-collapsed blobs where the intra- and inter-rule separators
    collide.
  * A manual rule carries {Ligand, Date, Cell_Line, Condition, Replicate, Row}; an auto
    rule carries {Ligand, Cell_Line, Condition, Date, well_id} (NO Replicate — it is
    resolved by its well_id, which pins row+column).
  * Each rule resolves to a (File_Name, Well_ID) set by filtering the master, then
    intersecting with Is_Excluded == True so inference can never claim a non-excluded well.

No Tkinter. Everything operates on master_df + a ProcessingConfig.

Public API (consumed by the UI layer):
    is_restorable(master_df, file_name, well_id) -> (bool, reason)
    list_active_exclusions(master_df)            -> list[dict]
    restore_rule(master_df, label, config)       -> dict  (report)
    restore_wells(master_df, label, wells, config) -> dict (report)

Helpers also imported elsewhere:
    parse_exclusion_blob(blob) -> list[dict]
    migrate_blob_separator(blob) -> str        (used by export.ensure_master_csv_schema)
"""

import re
import logging
import pandas as pd

from processing import recompute_master_after_exclusion, coerce_bool

logger = logging.getLogger("NCollector")

# Marker that begins a manual / auto rule token.
_MARKER_RE = re.compile(r'(?:Ligand:|AUTO:)')

LEGACY_LABEL = "Legacy exclusions"

# --- Rule grammar (single source of truth for the parser) --------------------
# These mirror the blob FORMAT emitted by the writers and must stay in step with them:
#   * manual rules: app.apply_exclusions / compile_master_dataframe
#       "Ligand: X | Date: Y | Cell: Z | Cond: W | Rep:R | Row:r"
#   * auto rules:   processing.format_warning_str
#       "AUTO: [TAG]   {ligand} | {cell} | {cond} | {date} | {well_id} - value: {v}"
# _MANUAL_KEY_MAP maps the blob's abbreviated keys to canonical master columns;
# the *_DEFAULTS dicts list the canonical criteria keys with their "match-all" defaults
# ("All"/"" both mean "do not filter on this field").
_MANUAL_KEY_MAP = {"Ligand": "Ligand", "Date": "Date", "Cell": "Cell_Line",
                   "Cond": "Condition", "Rep": "Replicate", "Row": "Row"}
_MANUAL_DEFAULTS = {"Ligand": "All", "Date": "All", "Cell_Line": "All",
                    "Condition": "All", "Replicate": "", "Row": ""}
_AUTO_DEFAULTS = {"Ligand": "All", "Cell_Line": "All", "Condition": "All",
                  "Date": "All", "well_id": ""}


# --------------------------------------------------------------------------- #
# Blob tokenizing / parsing
# --------------------------------------------------------------------------- #

def _marker_tokenize(blob: str) -> list[str]:
    """
    Re-tokenize a blob into rule strings by their START MARKERS, ignoring whatever
    separator (" | " or " || ") sits between them. Returns [] if no markers are found
    (e.g. "None", empty, or pre-v2 free text), which the callers treat as opaque/legacy.
    """
    s = "" if blob is None else str(blob)
    matches = list(_MARKER_RE.finditer(s))
    if not matches:
        return []
    tokens = []
    for i, m in enumerate(matches):
        start = m.start()
        end = matches[i + 1].start() if i + 1 < len(matches) else len(s)
        tok = s[start:end].strip()
        # Trim any trailing separator residue left between this rule and the next.
        tok = tok.rstrip().rstrip('|').rstrip()
        if tok:
            tokens.append(tok)
    return tokens


def migrate_blob_separator(blob) -> str:
    """
    Rewrite a blob to the v2.0.5 " || " separator. Idempotent.
      * No markers (None / empty / opaque free text) -> returned unchanged (it stays a
        single opaque token, surfaced later as the LEGACY_LABEL entry).
      * Otherwise -> markers re-tokenized and re-joined with " || ".
    """
    s = "" if blob is None else str(blob)
    if not s.strip() or s.strip().lower() == "none":
        return s
    tokens = _marker_tokenize(s)
    if not tokens:
        return s
    return " || ".join(tokens)


def _parse_manual_token(tok: str) -> dict:
    """Parse 'Ligand: X | Date: Y | Cell: Z | Cond: W | Rep:R | Row:r' -> criteria dict."""
    crit = dict(_MANUAL_DEFAULTS)
    for field in tok.split(" | "):
        if ":" not in field:
            continue
        k, v = field.split(":", 1)
        k, v = k.strip(), v.strip()
        if k in _MANUAL_KEY_MAP:
            crit[_MANUAL_KEY_MAP[k]] = v
    return crit


def _parse_auto_token(tok: str) -> dict:
    """
    Parse 'AUTO: [LOW LUM]   {ligand} | {cell} | {cond} | {date} | {well} - value: {v}'
    -> {Ligand, Cell_Line, Condition, Date, well_id}. (No Replicate — well_id pins it.)
    """
    crit = dict(_AUTO_DEFAULTS)
    body = tok[len("AUTO:"):].strip()
    # Drop a leading "[TAG]" if present (e.g. [LOW LUM], [VEHICLE WARN], [SPLIT]).
    body = re.sub(r'^\[[^\]]*\]\s*', '', body).strip()
    parts = [p.strip() for p in body.split(" | ")]
    if len(parts) >= 1 and parts[0]:
        crit["Ligand"] = parts[0]
    if len(parts) >= 2:
        crit["Cell_Line"] = parts[1]
    if len(parts) >= 3:
        crit["Condition"] = parts[2]
    if len(parts) >= 4:
        crit["Date"] = parts[3]
    if len(parts) >= 5:
        tail = parts[4]
        # tail looks like "{well_id} - value: {v}"
        well = tail.split(" - value:")[0].strip() if " - value:" in tail else tail.strip()
        crit["well_id"] = well
    return crit


def parse_exclusion_blob(blob) -> list[dict]:
    """
    Parse an Applied_Exclusions blob into a list of rule entries:
        {"label": <token text>, "source": "manual"|"AUTO", "criteria": {...}}
    Returns [] when the blob has no parseable markers (opaque / legacy / "None").
    """
    tokens = _marker_tokenize(blob)
    parsed = []
    for tok in tokens:
        if tok.startswith("AUTO:"):
            parsed.append({"label": tok, "source": "AUTO",
                           "criteria": _parse_auto_token(tok)})
        elif tok.startswith("Ligand:"):
            parsed.append({"label": tok, "source": "manual",
                           "criteria": _parse_manual_token(tok)})
    return parsed


# --------------------------------------------------------------------------- #
# Resolution: rule criteria -> (File_Name, Well_ID) set
# --------------------------------------------------------------------------- #

def _norm_dates(df: pd.DataFrame) -> pd.Series:
    """master_df['Date'] normalised to the '%d.%m.%y' dropdown / warning form."""
    # format='mixed' lets pandas parse heterogeneous date strings per-element without
    # emitting a UserWarning on every call; unparseable values coerce to NaT.
    return pd.to_datetime(df['Date'], errors='coerce', format='mixed').dt.strftime('%d.%m.%y')


def _excluded_mask(df: pd.DataFrame) -> pd.Series:
    if 'Is_Excluded' in df.columns:
        return df['Is_Excluded'].map(coerce_bool)
    return pd.Series(False, index=df.index)


def _resolve_criteria(df: pd.DataFrame, source: str, criteria: dict) -> pd.Series:
    """
    Boolean mask over df for one rule's criteria — the SAME filter semantics
    apply_exclusions uses (manual: filter by criteria then expand Row x matched columns;
    auto: the named well_id within files matching Date/Cell/Cond/Ligand). NOT yet
    intersected with Is_Excluded.
    """
    mask = pd.Series(True, index=df.index)
    norm_dates = _norm_dates(df) if 'Date' in df.columns else None

    def _eq(col, val):
        return df[col].astype(str).str.strip() == str(val).strip()

    if source == "AUTO":
        lig = criteria.get("Ligand", "All")
        if lig not in ("All", "") and 'Ligand' in df.columns:
            mask &= _eq('Ligand', lig)
        cell = criteria.get("Cell_Line", "All")
        if cell not in ("All", "") and 'Cell_Line' in df.columns:
            mask &= _eq('Cell_Line', cell)
        cond = criteria.get("Condition", "All")
        if cond not in ("All", "") and 'Transfection' in df.columns:
            mask &= _eq('Transfection', cond)
        date = criteria.get("Date", "All")
        if date not in ("All", "") and norm_dates is not None:
            mask &= (norm_dates == date)
        well = criteria.get("well_id", "")
        if well and 'Well_ID' in df.columns:
            mask &= _eq('Well_ID', well)
        else:
            # No well_id -> cannot pin; refuse to match anything (safer than over-match).
            mask &= False
        return mask

    # manual
    lig = criteria.get("Ligand", "All")
    if lig not in ("All", "") and 'Ligand' in df.columns:
        mask &= _eq('Ligand', lig)
    date = criteria.get("Date", "All")
    if date not in ("All", "") and norm_dates is not None:
        mask &= (norm_dates == date)
    cell = criteria.get("Cell_Line", "All")
    if cell not in ("All", "") and 'Cell_Line' in df.columns:
        mask &= _eq('Cell_Line', cell)
    cond = criteria.get("Condition", "All")
    if cond not in ("All", "") and 'Transfection' in df.columns:
        mask &= _eq('Transfection', cond)
    rep = criteria.get("Replicate", "")
    if rep not in ("All", "") and 'Replicate' in df.columns:
        mask &= _eq('Replicate', rep)
    row = criteria.get("Row", "")
    if row not in ("All", "") and 'Plate_Row' in df.columns:
        mask &= _eq('Plate_Row', row)
    return mask


def _resolve_rule_wells(df: pd.DataFrame, entry: dict, only_excluded: bool = True) -> set:
    """Resolve one parsed rule entry to a set of (File_Name, Well_ID) tuples."""
    mask = _resolve_criteria(df, entry["source"], entry["criteria"])
    if only_excluded:
        mask &= _excluded_mask(df)
    sub = df.loc[mask]
    if sub.empty:
        return set()
    return set(zip(sub['File_Name'].astype(str), sub['Well_ID'].astype(str)))


def _all_excluded_wells(df: pd.DataFrame) -> set:
    m = _excluded_mask(df)
    sub = df.loc[m]
    return set(zip(sub['File_Name'].astype(str), sub['Well_ID'].astype(str)))


# --------------------------------------------------------------------------- #
# Restore feasibility
# --------------------------------------------------------------------------- #

def is_restorable(master_df: pd.DataFrame, file_name: str, well_id: str):
    """
    (bool, reason). True iff a usable pristine source exists for the well — either
    Raw_BRET_unexcluded non-NaN, OR (Donor non-NaN AND Acceptor non-NaN) at any timepoint.
    """
    rows = master_df[(master_df['File_Name'].astype(str) == str(file_name)) &
                     (master_df['Well_ID'].astype(str) == str(well_id))]
    if rows.empty:
        return False, f"well {file_name}/{well_id} not found in master"

    if 'Raw_BRET_unexcluded' in rows.columns:
        if pd.to_numeric(rows['Raw_BRET_unexcluded'], errors='coerce').notna().any():
            return True, ""

    has_donor = ('Donor_Raw_kinetic' in rows.columns
                 and pd.to_numeric(rows['Donor_Raw_kinetic'], errors='coerce').notna().any())
    has_acceptor = ('Acceptor_Raw_kinetic' in rows.columns
                    and pd.to_numeric(rows['Acceptor_Raw_kinetic'], errors='coerce').notna().any())
    if has_donor and has_acceptor:
        return True, ""

    return False, ("no Raw_BRET_unexcluded and no Donor/Acceptor channels — "
                   "likely a legacy master excluded before v2.0.5")


# --------------------------------------------------------------------------- #
# Listing
# --------------------------------------------------------------------------- #

def _split_restorable(master_df, wells):
    """Split a (File_Name, Well_ID) set into restorable vs non-restorable (with reasons)."""
    restorable, non_restorable = [], []
    for f, w in sorted(wells):
        ok, reason = is_restorable(master_df, f, w)
        if ok:
            restorable.append((f, w))
        else:
            non_restorable.append((f, w, reason))
    return restorable, non_restorable


def list_active_exclusions(master_df: pd.DataFrame) -> list[dict]:
    """
    Describe the currently-active exclusions. Each entry:
        {
          "label": <rule text or "Legacy exclusions">,
          "source": "AUTO" | "manual",
          "wells": [(file, well), ...],                # all currently-excluded wells
          "restorable": [(file, well), ...],
          "non_restorable": [(file, well, reason), ...],
        }
    A blob with no parseable markers yields a single LEGACY_LABEL entry covering ALL
    Is_Excluded == True wells (per-well restorable; just no rule granularity).
    """
    if master_df is None or master_df.empty:
        return []

    blob = master_df['Applied_Exclusions'].dropna().iloc[0] \
        if 'Applied_Exclusions' in master_df.columns and master_df['Applied_Exclusions'].notna().any() \
        else ""
    entries = parse_exclusion_blob(blob)

    if not entries:
        wells = _all_excluded_wells(master_df)
        if not wells:
            return []
        restorable, non_restorable = _split_restorable(master_df, wells)
        return [{
            "label": LEGACY_LABEL,
            "source": "manual",
            "wells": sorted(wells),
            "restorable": restorable,
            "non_restorable": non_restorable,
        }]

    out = []
    for e in entries:
        wells = _resolve_rule_wells(master_df, e, only_excluded=True)
        restorable, non_restorable = _split_restorable(master_df, wells)
        out.append({
            "label": e["label"],
            "source": e["source"],
            "wells": sorted(wells),
            "restorable": restorable,
            "non_restorable": non_restorable,
        })
    return out


# --------------------------------------------------------------------------- #
# Per-well token synthesis (well-pinned, same shape auto rules use)
# --------------------------------------------------------------------------- #

def _well_token(master_df: pd.DataFrame, file_name: str, well_id: str) -> str:
    """
    Build a single well-pinned token (auto form) that re-resolves to exactly this well:
        'AUTO: [RESTORE-SPLIT] {ligand} | {cell} | {cond} | {date} | {well_id} - value: n/a'
    """
    rows = master_df[(master_df['File_Name'].astype(str) == str(file_name)) &
                     (master_df['Well_ID'].astype(str) == str(well_id))]
    if rows.empty:
        return f"AUTO: [RESTORE-SPLIT] All | All | All | All | {well_id} - value: n/a"
    first = rows.iloc[0]
    lig = str(first.get('Ligand', 'All'))
    cell = str(first.get('Cell_Line', 'All'))
    cond = str(first.get('Transfection', 'All'))
    try:
        date = pd.to_datetime(first.get('Date')).strftime('%d.%m.%y')
    except Exception:
        date = str(first.get('Date', 'All'))
    return f"AUTO: [RESTORE-SPLIT] {lig} | {cell} | {cond} | {date} | {well_id} - value: n/a"


# --------------------------------------------------------------------------- #
# Restore core
# --------------------------------------------------------------------------- #

def _write_blob(master_df: pd.DataFrame, tokens: list[str]):
    """Broadcast the global Applied_Exclusions blob to every row."""
    new_blob = " || ".join(tokens) if tokens else "None"
    master_df['Applied_Exclusions'] = new_blob


def _set_excluded(master_df: pd.DataFrame, wells: set, value: bool):
    """Flip Is_Excluded for a set of (File_Name, Well_ID) tuples (in place)."""
    if not wells:
        return
    key = pd.MultiIndex.from_arrays(
        [master_df['File_Name'].astype(str), master_df['Well_ID'].astype(str)])
    target = pd.Series(key.isin(list(wells)), index=master_df.index)
    master_df.loc[target, 'Is_Excluded'] = value


def _recompute_in_place(master_df: pd.DataFrame, affected_files: list, config):
    """Recompute affected files via the engine and copy results back into master_df."""
    if not affected_files:
        return
    recomputed = recompute_master_after_exclusion(master_df, list(affected_files), config)
    # recompute_master_after_exclusion preserves the row index; copy every column back so
    # the caller's master_df object stays valid (Raw_BRET_unexcluded is untouched).
    for col in recomputed.columns:
        master_df[col] = recomputed[col].values


def _restore_core(master_df, label, wells_subset, config):
    """
    Shared implementation for restore_rule (wells_subset=None -> all of the rule's wells)
    and restore_wells (wells_subset = an explicit subset). Mutates master_df in place
    (Is_Excluded + Applied_Exclusions + recomputed derived columns) and returns a report.

    Terminology used below:
      * target rule (target_entry / target_idx): the rule named by `label` — the one being
        restored. target_wells is the set of (File_Name, Well_ID) it currently excludes.
      * other rules (other_entries) / other_cov: every OTHER still-active rule and the wells
        they currently cover. A well also covered by another rule is NOT restored here —
        that rule keeps attributing it (see the per-well decision).
    """
    report = {
        "restored_count": 0,
        "nonrestorable": [],
        "blocked_by_other_rule_count": 0,
        "rule_decomposed": False,
        "affected_files": [],
    }
    if master_df is None or master_df.empty:
        return report

    blob = master_df['Applied_Exclusions'].dropna().iloc[0] \
        if 'Applied_Exclusions' in master_df.columns and master_df['Applied_Exclusions'].notna().any() \
        else ""
    entries = parse_exclusion_blob(blob)

    # --- Identify the target rule and the "other" active rules ---
    legacy = False
    if not entries:
        # Opaque/legacy blob -> the single LEGACY_LABEL entry covers all excluded wells.
        if label != LEGACY_LABEL:
            logger.warning(f"restore: label '{label}' not found (legacy/opaque blob).")
            return report
        legacy = True
        target_entry = None
        other_entries = []
        target_wells = _all_excluded_wells(master_df)
    else:
        target_idx = next((i for i, e in enumerate(entries) if e["label"] == label), None)
        if target_idx is None:
            logger.warning(f"restore: rule label not found: {label!r}")
            return report
        target_entry = entries[target_idx]
        other_entries = [e for j, e in enumerate(entries) if j != target_idx]
        target_wells = _resolve_rule_wells(master_df, target_entry, only_excluded=True)

    # OTHER coverage computed against the PRE-restore Is_Excluded state.
    other_cov = set()
    for e in other_entries:
        other_cov |= _resolve_rule_wells(master_df, e, only_excluded=True)

    # --- Determine the subset to attempt restoring ---
    if wells_subset is None:
        subset = set(target_wells)
    else:
        want = {(str(f), str(w)) for f, w in wells_subset}
        subset = want & set(target_wells)

    # --- Per-well decision ---
    restored, blocked, nonrestorable = set(), set(), []
    for (f, w) in subset:
        if (f, w) in other_cov:
            blocked.add((f, w))
            continue
        ok, reason = is_restorable(master_df, f, w)
        if ok:
            restored.add((f, w))
        else:
            nonrestorable.append((f, w, reason))

    # --- Apply: flip restored wells to not-excluded ---
    _set_excluded(master_df, restored, False)

    # --- Decide this rule's rewritten coverage (orphan-free) ---
    # Wells that remain excluded AND are NOT covered by any other active rule must stay
    # attributed by THIS rule. = (target_wells \ restored) \ other_cov.
    remaining_excluded = set(target_wells) - restored
    must_cover = remaining_excluded - other_cov

    if legacy:
        # Opaque blob: re-express the still-excluded wells as explicit per-well tokens.
        all_still = _all_excluded_wells(master_df)
        if all_still:
            new_tokens = [_well_token(master_df, f, w) for (f, w) in sorted(all_still)]
        else:
            new_tokens = []
        _write_blob(master_df, new_tokens)
        report["rule_decomposed"] = bool(restored) and bool(all_still)
    else:
        tokens = [e["label"] for e in entries]
        keep_original = (must_cover == set(target_wells)) and (len(restored) == 0)
        if keep_original:
            # Nothing this rule uniquely covers changed; leave its token as-is.
            new_tokens = tokens
        elif not must_cover:
            # Rule no longer uniquely covers anything -> drop its token entirely.
            new_tokens = [t for j, t in enumerate(tokens) if j != target_idx]
        else:
            # Replace the rule with well-pinned tokens for exactly its remaining wells.
            replacement = [_well_token(master_df, f, w) for (f, w) in sorted(must_cover)]
            new_tokens = []
            for j, t in enumerate(tokens):
                if j == target_idx:
                    new_tokens.extend(replacement)
                else:
                    new_tokens.append(t)
            report["rule_decomposed"] = True
        _write_blob(master_df, new_tokens)

    # --- Recompute the files whose Is_Excluded actually changed ---
    affected = sorted({f for (f, w) in restored})
    _recompute_in_place(master_df, affected, config)

    report["restored_count"] = len(restored)
    report["nonrestorable"] = sorted(nonrestorable)
    report["blocked_by_other_rule_count"] = len(blocked)
    report["affected_files"] = affected
    return report


def restore_rule(master_df: pd.DataFrame, label: str, config) -> dict:
    """
    Restore ALL of a rule's wells (equivalent to restore_wells over the full set). Per
    well: still matched by ANOTHER active rule -> stays excluded (the other rule keeps
    attributing it); else restorable -> Is_Excluded=False; else (non-restorable, only this
    rule) -> stays excluded. The blob is rewritten so the rule still matches EXACTLY its
    remaining excluded wells (the non-restorable leftovers), decomposed into per-well
    tokens; the token is dropped entirely if zero wells remain. Recomputes affected files.

    Report keys: restored_count, nonrestorable [(file, well, reason)],
    blocked_by_other_rule_count, rule_decomposed (bool), affected_files.
    """
    return _restore_core(master_df, label, None, config)


def restore_wells(master_df: pd.DataFrame, label: str, wells, config) -> dict:
    """
    Per-well refinement over a SUBSET of a rule's wells (mostly for manual criteria-rules;
    auto rules are already single-well). Same per-well decision logic as restore_rule. The
    blob is then rewritten so the rule ALWAYS re-resolves to exactly the wells it still
    excludes: the rule's token is dropped and replaced by well-pinned tokens covering only
    its still-excluded, uniquely-attributed wells (so a broad criteria-rule can never
    re-resolve to a well you just restored). Recomputes affected files.

    `wells` is an iterable of (File_Name, Well_ID).
    Report keys as in restore_rule.
    """
    return _restore_core(master_df, label, list(wells), config)