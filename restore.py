"""
restore.py — GUI-agnostic backend for REVERSIBLE exclusions ("restore").

Exclusion in N Collector is per-row authoritative via Is_Excluded; the engine
(recompute_master_after_exclusion) derives every downstream column from the PRISTINE
raw (Raw_BRET_unexcluded -> Acceptor/Donor -> NaN) with the currently-excluded wells
NaN-d. Exclude and restore are therefore the SAME operation: change Is_Excluded, then
recompute. This module exposes the pure logic that decides WHICH wells a restore touches
and keeps the human-readable Applied_Exclusions blob in sync.

Reverting, end to end:
    decide per well  ->  flip Is_Excluded back to off  ->  fix the blob note  ->  redo the
    downstream math once

There is NO per-well "which rule excluded this" record. Which wells
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
    list_active_exclusions(master_df, ctx=None)  -> list[dict]
    restore_rule(master_df, label, config, recompute=True, ctx=None)        -> dict (report)
    restore_wells(master_df, label, wells, config, recompute=True, ctx=None) -> dict (report)
    build_resolve_ctx(master_df) -> _ResolveCtx

The last two arguments on the restore/list calls are optional speed levers for a BATCH
revert (reverting many rules at once); a single one-off call can ignore them:
    * recompute (restore_* only): True (default) redoes the downstream math inside this
      call. Pass False to SKIP it and recompute once yourself afterwards over the union of
      every report's "affected_files" — so N rules trigger 1 recompute, not N.
    * ctx: a shared resolution "notebook" from build_resolve_ctx. Build it ONCE and pass the
      same one to every restore_*/list_active_exclusions call in the batch, so the slow
      resolving work happens a single time instead of once per (rule x rule).
    (See app.revert_exclusions for the canonical batch pattern: build_resolve_ctx -> loop
    with recompute=False, ctx=ctx -> one recompute_master_after_exclusion at the end.)

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
                   "Cond": "Condition", "Rep": "Replicate", "Row": "Row",
                   "Main": "Main_Plasmids", "File": "File_Name"}
_MANUAL_DEFAULTS = {"Ligand": "All", "Date": "All", "Cell_Line": "All",
                    "Condition": "All", "Replicate": "", "Row": "",
                    "Main_Plasmids": "All", "File_Name": "All"}
_AUTO_DEFAULTS = {"Ligand": "All", "Cell_Line": "All", "Condition": "All",
                  "Date": "All", "well_id": "", "File_Name": "All"}


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
    # Optional trailing keyed segment "| File: {fname}" (additive; absent -> stays "All",
    # so every legacy AUTO token parses exactly as before).
    for p in parts[5:]:
        if ":" in p:
            k, v = p.split(":", 1)
            if k.strip() == "File":
                crit["File_Name"] = v.strip()
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
# Pending-rule resolution (exclude/revert dropdowns) -> targets / wells
# --------------------------------------------------------------------------- #
# These resolve the PENDING rule dicts the GUI builds from its dropdowns directly
# against the live master_df. Distinct from _ResolveCtx below, which resolves
# already-applied blob entries (AUTO/manual tokens). Kept together here so all
# exclusion resolution lives in one module.

# Columns that together identify which master rows an exclusion target refers
# to. File_Name + Well_ID alone could collide if two loaded files happen to
# share a name (e.g. same filename in two folders), so the match is widened to
# the well's biological identity. Both resolve_rule_to_targets (which builds
# the target tuples) and apply_exclusions (which builds the match mask) read
# this list, so the two keys are guaranteed to line up by construction.
EXCLUSION_KEY_COLS = ('File_Name', 'Well_ID', 'Ligand', 'Transfection',
                      'Cell_Line', 'Main_Plasmids')


def exclusion_key_cols(df):
    """The exclusion-key columns actually present in df, in fixed order.
    File_Name + Well_ID are always present; the identity columns are added
    when available (they are part of MASTER_COLUMNS, so normally all six)."""
    return [c for c in EXCLUSION_KEY_COLS if c in df.columns]


def rule_mask(df, rule, norm_dates, ligand_locked=False):
    """
    Boolean mask over df for ONE pending criteria rule — the single source of
    the pending-rule filter semantics, shared by both the exclude path
    (resolve_rule_to_targets) and the revert path (resolve_rule_to_wells) so the two
    can never drift. "All ligands" also covers the single-locked-ligand case
    (ligand_locked, i.e. the GUI's cb_lig disabled). norm_dates is df['Date']
    pre-normalised to '%d.%m.%y'.
    """
    mask = pd.Series(True, index=df.index)

    ligand_is_all = (rule.get('Ligand', 'All') in ("All", "") or ligand_locked)
    if not ligand_is_all:
        mask &= (df['Ligand'].astype(str) == str(rule['Ligand']))
    if rule.get('Date', 'All') != "All":
        mask &= (norm_dates == rule['Date'])
    if rule.get('Cell_Line', 'All') != "All":
        mask &= (df['Cell_Line'].astype(str) == str(rule['Cell_Line']))
    if rule.get('Condition', 'All') != "All":
        # Index vocabulary "Condition" maps to master_df "Transfection".
        mask &= (df['Transfection'].astype(str) == str(rule['Condition']))
    rep = rule.get('Replicate', 'All')
    if rep not in ("All", ""):
        mask &= (df['Replicate'].astype(str) == str(rep))
    row = rule.get('Row', 'All')
    if row not in ("All", ""):
        mask &= (df['Plate_Row'].astype(str) == str(row))
    # Scope narrowing (post-merge masters). Skipped when "All"/empty
    main = rule.get('Main_Plasmids', 'All')
    if main not in ("All", "") and 'Main_Plasmids' in df.columns:
        mask &= (df['Main_Plasmids'].astype(str) == str(main))
    fname = rule.get('File_Name', 'All')
    if fname not in ("All", "") and 'File_Name' in df.columns:
        mask &= (df['File_Name'].astype(str) == str(fname))
    return mask


def resolve_rule_to_targets(df, rule, norm_dates, ligand_locked=False, log_fn=None):
    """
    Resolves one pending exclusion rule to a set of identity tuples directly on
    df (the single source of truth). Same rule semantics as before,
    including the "whole date" branch — which resolves to every matching
    well on that date.

    Each returned tuple is keyed on EXCLUSION_KEY_COLS — i.e. not just
    (File_Name, Well_ID) but also Ligand / Transfection / Cell_Line /
    Main_Plasmids — so that a duplicate file name cannot cause the wrong rows
    to be flagged. apply_exclusions builds its match mask from the same column
    list, keeping the two sides consistent.

    norm_dates is df['Date'] pre-normalised to the '%d.%m.%y' dropdown form.
    """
    mask = rule_mask(df, rule, norm_dates, ligand_locked)

    sub = df.loc[mask]
    if sub.empty:
        if log_fn:
            log_fn(f"   [WARNING] Rule {rule} matched 0 records.")
        return set()
    key_cols = exclusion_key_cols(df)
    return set(zip(*[sub[c].astype(str) for c in key_cols]))


def resolve_rule_to_wells(df, rule, norm_dates, ligand_locked=False):
    """
    Resolve one pending criteria rule to a set of (File_Name, Well_ID) tuples, using
    the SAME mask as resolve_rule_to_targets (via rule_mask). Used by the revert
    path, which needs plain (file, well) pairs for the restore API rather than the
    full exclusion-identity tuples the add-only flag mask uses.
    """
    mask = rule_mask(df, rule, norm_dates, ligand_locked)
    sub = df.loc[mask]
    if sub.empty:
        return set()
    return set(zip(sub['File_Name'].astype(str), sub['Well_ID'].astype(str)))


# --------------------------------------------------------------------------- #
# Resolution: rule criteria -> (File_Name, Well_ID) set
# --------------------------------------------------------------------------- #

def _vectorized_excluded_mask(df: pd.DataFrame) -> pd.Series:
    """
    Vectorized equivalent of df['Is_Excluded'].map(coerce_bool). A real-bool column is
    returned as-is; an object/string/numeric column (e.g. from pd.read_csv) is coerced
    with the same truthiness rules as coerce_bool, but without the per-element Python call.
    """
    if 'Is_Excluded' not in df.columns:
        return pd.Series(False, index=df.index)
    s = df['Is_Excluded']
    if s.dtype == bool:
        return s
    # If master has some manual edited values (also see processing.coerce_bool)
    as_str = s.astype('string').str.strip().str.lower()
    truthy_str = as_str.isin(["true", "1", "1.0", "yes"])
    num = pd.to_numeric(s, errors='coerce')
    truthy_num = num.notna() & (num != 0)
    return (truthy_str | truthy_num).fillna(False)


class _ResolveCtx:
    """
    Helper class serving as per-revert resolution cache. Re-resolving rules was the dominant cost of a many-rule
    revert: every restore call re-derived, over the WHOLE master, the normalised dates
    (pd.to_datetime+strftime), the stripped criteria columns, and every other rule's mask
    — i.e. O(rules**2) full-master scans.

    Everything resolution needs is INVARIANT across a revert except Is_Excluded:
      * normalised Date strings, stripped criteria columns        -> precomputed once
      * each rule's criteria match mask (pre-Is_Excluded)         -> cached per criteria
      * per-well restorability                                    -> precomputed once
      * the Is_Excluded mask                                      -> recomputed (vectorized)
        only after wells are actually flipped (invalidate_excluded()).

    The context must be built on (and passed the SAME object as) the master_df being
    mutated, so excluded_mask() reflects live Is_Excluded flips.
    """
    def __init__(self, df: pd.DataFrame):
        # Once cache is created, collect all data that does not change with Is_Excluded switch
        self.df = df
        self.index = df.index
        self._norm_dates = (pd.to_datetime(df['Date'], errors='coerce', format='mixed')
                            .dt.strftime('%d.%m.%y')) if 'Date' in df.columns else None
        # Stripped string views of the columns rules filter on (computed once)
        self._scol = {}
        for c in ('File_Name', 'Well_ID', 'Ligand', 'Cell_Line', 'Transfection',
                  'Plate_Row', 'Replicate', 'Main_Plasmids'):
            if c in df.columns:
                self._scol[c] = df[c].astype(str).str.strip()
        self._crit_cache = {}          # (source, sorted-criteria-items) -> boolean mask
        self._excl_mask = None
        self._restorable_map = None    # invariant; built lazily

    # --- Is_Excluded (mutates during revert) ---
    def invalidate_excluded(self):
        self._excl_mask = None

    def excluded_mask(self) -> pd.Series:
        if self._excl_mask is None:
            self._excl_mask = _vectorized_excluded_mask(self.df)
        return self._excl_mask

    # --- restorability (invariant during revert) ---
    def restorable_map(self) -> dict:
        if self._restorable_map is None:
            self._restorable_map = _build_restorable_map(self.df)
        return self._restorable_map

    # --- criteria matching (invariant; cached per criteria) ---
    def _eq(self, col, val):
        s = self._scol.get(col)
        if s is None:
            return pd.Series(True, index=self.index)
        return s == str(val).strip()

    def criteria_mask(self, source: str, criteria: dict) -> pd.Series:
        # Deciphers which wells are targeted by rules (criteria)
        key = (source, tuple(sorted(criteria.items())))
        cached = self._crit_cache.get(key)
        if cached is not None:
            return cached

        mask = pd.Series(True, index=self.index)

        if source == "AUTO":
            lig = criteria.get("Ligand", "All")
            if lig not in ("All", "") and 'Ligand' in self._scol:
                mask &= self._eq('Ligand', lig)
            cell = criteria.get("Cell_Line", "All")
            if cell not in ("All", "") and 'Cell_Line' in self._scol:
                mask &= self._eq('Cell_Line', cell)
            cond = criteria.get("Condition", "All")
            if cond not in ("All", "") and 'Transfection' in self._scol:
                mask &= self._eq('Transfection', cond)
            date = criteria.get("Date", "All")
            if date not in ("All", "") and self._norm_dates is not None:
                mask &= (self._norm_dates == date)
            well = criteria.get("well_id", "")
            if well and 'Well_ID' in self._scol:
                mask &= self._eq('Well_ID', well)
            else:
                # No well_id -> cannot pin; refuse to match anything.
                mask &= False
            fname = criteria.get("File_Name", "All")
            if fname not in ("All", "") and 'File_Name' in self._scol:
                mask &= self._eq('File_Name', fname)
        else:
            # manual
            lig = criteria.get("Ligand", "All")
            if lig not in ("All", "") and 'Ligand' in self._scol:
                mask &= self._eq('Ligand', lig)
            date = criteria.get("Date", "All")
            if date not in ("All", "") and self._norm_dates is not None:
                mask &= (self._norm_dates == date)
            cell = criteria.get("Cell_Line", "All")
            if cell not in ("All", "") and 'Cell_Line' in self._scol:
                mask &= self._eq('Cell_Line', cell)
            cond = criteria.get("Condition", "All")
            if cond not in ("All", "") and 'Transfection' in self._scol:
                mask &= self._eq('Transfection', cond)
            rep = criteria.get("Replicate", "")
            if rep not in ("All", "") and 'Replicate' in self._scol:
                mask &= self._eq('Replicate', rep)
            row = criteria.get("Row", "")
            if row not in ("All", "") and 'Plate_Row' in self._scol:
                mask &= self._eq('Plate_Row', row)
            main = criteria.get("Main_Plasmids", "All")
            if main not in ("All", "") and 'Main_Plasmids' in self._scol:
                mask &= self._eq('Main_Plasmids', main)
            fname = criteria.get("File_Name", "All")
            if fname not in ("All", "") and 'File_Name' in self._scol:
                mask &= self._eq('File_Name', fname)

        self._crit_cache[key] = mask
        return mask

    def rule_wells(self, entry: dict, only_excluded: bool = True) -> set:
        # Resolving rules to actual set of wells
        mask = self.criteria_mask(entry["source"], entry["criteria"])
        if only_excluded:
            mask = mask & self.excluded_mask()
        if not mask.any():
            return set()
        fcol = self._scol.get('File_Name')
        wcol = self._scol.get('Well_ID')
        if fcol is None or wcol is None:
            sub = self.df.loc[mask]
            return set(zip(sub['File_Name'].astype(str), sub['Well_ID'].astype(str)))
        return set(zip(fcol[mask], wcol[mask]))

    def all_excluded_wells(self) -> set:
        # Every File x Well that is currently switched off, regardless of which rule (used for legacy exclusions)
        mask = self.excluded_mask()
        fcol = self._scol.get('File_Name')
        wcol = self._scol.get('Well_ID')
        if fcol is None or wcol is None:
            sub = self.df.loc[mask]
            return set(zip(sub['File_Name'].astype(str), sub['Well_ID'].astype(str)))
        return set(zip(fcol[mask], wcol[mask]))


def build_resolve_ctx(master_df: pd.DataFrame) -> _ResolveCtx:
    """
    Build a resolution context (_ResolveCtx) for a batch revert. Construct ONCE and pass it to every
    restore_rule / restore_wells / list_active_exclusions call in the same revert so the
    invariant resolution work (dates, stripped columns, per-rule masks, restorability) is
    done a single time instead of once per (rule x rule). Must be built on the same
    master_df object those calls mutate.
    """
    return _ResolveCtx(master_df)


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

_NON_RESTORABLE_REASON = ("no Raw_BRET_unexcluded and no Donor/Acceptor channels — "
                          "likely a legacy master excluded before v2.0.5")


def _build_restorable_map(master_df: pd.DataFrame) -> dict:
    """
    Vectorized restorability for EVERY (File_Name, Well_ID) in the master, computed in a
    SINGLE pass. Mirrors is_restorable exactly: a well is restorable iff
    Raw_BRET_unexcluded is non-NaN anywhere for that well, OR (Donor non-NaN anywhere AND
    Acceptor non-NaN anywhere).

    Returns {(file_name, well_id): (bool, reason)}.

    This replaces the previous pattern of calling is_restorable() per well inside the
    restore/list loops — each of those calls scanned the whole master twice (.astype(str)
    on two full columns), giving O(n_wells * n_rows) string churn. Here every per-well
    answer is precomputed once via one groupby, then looked up in O(1).
    """
    if master_df is None or master_df.empty:
        return {}

    cols = master_df.columns
    agg = pd.DataFrame({
        '_f': master_df['File_Name'].astype(str),
        '_w': master_df['Well_ID'].astype(str),
    })
    has_unexcl = 'Raw_BRET_unexcluded' in cols
    has_donor = 'Donor_Raw_kinetic' in cols
    has_acceptor = 'Acceptor_Raw_kinetic' in cols
    if has_unexcl:
        agg['unexcl'] = pd.to_numeric(master_df['Raw_BRET_unexcluded'], errors='coerce').notna().values
    if has_donor:
        agg['donor'] = pd.to_numeric(master_df['Donor_Raw_kinetic'], errors='coerce').notna().values
    if has_acceptor:
        agg['acceptor'] = pd.to_numeric(master_df['Acceptor_Raw_kinetic'], errors='coerce').notna().values

    # One pass: per (file, well), does ANY row carry usable pristine source data?
    grouped = agg.groupby(['_f', '_w'], sort=False).any()

    out = {}
    for key, row in grouped.iterrows():
        ok = bool(row['unexcl']) if has_unexcl else False
        if not ok and has_donor and has_acceptor:
            ok = bool(row['donor']) and bool(row['acceptor'])
        out[key] = (True, "") if ok else (False, _NON_RESTORABLE_REASON)
    return out


def _split_restorable(master_df, wells, rmap: dict | None = None):
    """
    Split a (File_Name, Well_ID) set into restorable vs non-restorable (with reasons).
    Pass a precomputed `rmap` as lookup (from _build_restorable_map) to avoid per-well scans; if
    omitted it is built once here.
    """
    if rmap is None:
        rmap = _build_restorable_map(master_df)
    restorable, non_restorable = [], []
    for f, w in sorted(wells):
        ok, reason = rmap.get((str(f), str(w)),
                              (False, f"well {f}/{w} not found in master"))
        if ok:
            restorable.append((f, w))
        else:
            non_restorable.append((f, w, reason))
    return restorable, non_restorable


def list_active_exclusions(master_df: pd.DataFrame, ctx=None) -> list[dict]:
    """
    Describe the currently-active exclusions (for GUI display). Each entry:
        {
          "label": <rule text or "Legacy exclusions">,
          "source": "AUTO" | "manual",
          "wells": [(file, well), ...],                # all currently-excluded wells
          "restorable": [(file, well), ...],
          "non_restorable": [(file, well, reason), ...],
        }
    A blob with no parseable markers yields a single LEGACY_LABEL entry covering ALL
    Is_Excluded == True wells (per-well restorable; just no rule granularity).

    Pass a shared ctx (build_resolve_ctx) to reuse precomputed resolution across a batch.
    """
    if master_df is None or master_df.empty:
        return []

    if ctx is None:
        ctx = _ResolveCtx(master_df)

    blob = master_df['Applied_Exclusions'].dropna().iloc[0] \
        if 'Applied_Exclusions' in master_df.columns and master_df['Applied_Exclusions'].notna().any() \
        else ""
    entries = parse_exclusion_blob(blob)

    # Restorability precomputed once for the whole master (via ctx); reused per entry.
    rmap = ctx.restorable_map()

    if not entries:
        wells = ctx.all_excluded_wells()
        if not wells:
            return []
        restorable, non_restorable = _split_restorable(master_df, wells, rmap)
        return [{
            "label": LEGACY_LABEL,
            "source": "manual",
            "wells": sorted(wells),
            "restorable": restorable,
            "non_restorable": non_restorable,
        }]

    out = []
    for e in entries:
        wells = ctx.rule_wells(e, only_excluded=True)
        restorable, non_restorable = _split_restorable(master_df, wells, rmap)
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

def well_token(master_df: pd.DataFrame, file_name: str, well_id: str) -> str:
    """
    Build a single well-pinned token (auto form) that re-resolves to exactly this well:
        'AUTO: [RESTORE-SPLIT] {ligand} | {cell} | {cond} | {date} | {well_id} - value: n/a | File: {file_name}'
    The trailing "| File: {file_name}" pins the token to one (File_Name, Well_ID) even when
    another source file shares the same Ligand/Cell/Cond/Date/Well_ID.

    Public emitter: used internally for partial-rule reverts and merge-time blob scoping,
    and by the app's merge layer to record union exclusions for conflicting duplicate files
    so the added wells stay shown and revertable.
    """
    rows = master_df[(master_df['File_Name'].astype(str) == str(file_name)) &
                     (master_df['Well_ID'].astype(str) == str(well_id))]
    if rows.empty:
        return (f"AUTO: [RESTORE-SPLIT] All | All | All | All | "
                f"{well_id} - value: n/a | File: {file_name}")
    first = rows.iloc[0]
    lig = str(first.get('Ligand', 'All'))
    cell = str(first.get('Cell_Line', 'All'))
    cond = str(first.get('Transfection', 'All'))
    try:
        date = pd.to_datetime(first.get('Date')).strftime('%d.%m.%y')
    except Exception:
        date = str(first.get('Date', 'All'))
    return (f"AUTO: [RESTORE-SPLIT] {lig} | {cell} | {cond} | {date} | "
            f"{well_id} - value: n/a | File: {file_name}")


# --------------------------------------------------------------------------- #
# Token emitters (inverse of the _parse_*_token parsers) + merge-time scoping
# --------------------------------------------------------------------------- #
# These are used by the merge engine (merge.py) to rewrite an imported source's
# Applied_Exclusions so that, AFTER concatenation, every rule re-resolves to EXACTLY
# the (File_Name, Well_ID) set it covered within its own source — no bleed onto the
# other sources' rows — while keeping criteria-rules readable instead of exploding
# them into per-well tokens. The shapes emitted here round-trip cleanly through
# parse_exclusion_blob / criteria_mask, so the existing restore machinery treats a
# merged blob exactly like any other.

def _emit_manual_token(crit: dict) -> str:
    """Inverse of _parse_manual_token: build a canonical manual rule string from a
    criteria dict (canonical master-column keys). Mirrors app.apply_exclusions' writer:
    the six core fields are always present; Main/File are appended only when meaningful."""
    s = (f"Ligand: {crit.get('Ligand', 'All')} | Date: {crit.get('Date', 'All')} | "
         f"Cell: {crit.get('Cell_Line', 'All')} | Cond: {crit.get('Condition', 'All')} | "
         f"Rep:{crit.get('Replicate', '')} | Row:{crit.get('Row', '')}")
    main = crit.get('Main_Plasmids', 'All')
    if main not in ("All", ""):
        s += f" | Main: {main}"
    fname = crit.get('File_Name', 'All')
    if fname not in ("All", ""):
        s += f" | File: {fname}"
    return s


def _emit_auto_token(crit: dict, tag: str = "MERGE") -> str:
    """Inverse of _parse_auto_token: build a canonical AUTO rule string (well-pinned,
    optionally File-pinned) from a criteria dict. Used to rewrite File references on
    auto rules during a merge remap. The trailing 'value: n/a' is cosmetic."""
    lig = crit.get('Ligand', 'All')
    cell = crit.get('Cell_Line', 'All')
    cond = crit.get('Condition', 'All')
    date = crit.get('Date', 'All')
    well = crit.get('well_id', '')
    s = f"AUTO: [{tag}] {lig} | {cell} | {cond} | {date} | {well} - value: n/a"
    fname = crit.get('File_Name', 'All')
    if fname not in ("All", ""):
        s += f" | File: {fname}"
    return s


def _emit_token(entry: dict) -> str:
    """Re-emit a parsed rule entry ({source, criteria}) to its canonical token string."""
    if entry["source"] == "AUTO":
        return _emit_auto_token(entry["criteria"])
    return _emit_manual_token(entry["criteria"])


def remap_blob_file_names(blob, rename_map: dict) -> str:
    """
    Rewrite File references inside an Applied_Exclusions blob according to rename_map
    ({old_File_Name: new_File_Name}). Only rules whose File criterion is an exact key of
    rename_map are touched; everything else is preserved verbatim (re-emitted in canonical
    form). Returns the rejoined blob (" || ").

    The merge engine calls this on a source's blob BEFORE scope_blob_for_merge whenever that
    source had files collision-renamed, so a File-pinned rule keeps pointing at the renamed
    rows rather than silently resolving to nothing. Idempotent for an empty/identity map.
    """
    if not rename_map:
        return "" if blob is None else str(blob)
    entries = parse_exclusion_blob(blob)
    if not entries:
        # Opaque / legacy / "None": nothing parseable to remap.
        return "" if blob is None else str(blob)
    out = []
    for e in entries:
        crit = dict(e["criteria"])
        old = crit.get("File_Name", "All")
        if old in rename_map:
            crit["File_Name"] = rename_map[old]
            out.append(_emit_token({"source": e["source"], "criteria": crit}))
        else:
            out.append(e["label"])
    return " || ".join(out)


# Columns walked, in readability-preference order, when scoping a manual rule at merge time

# ORDER MATTERS here: the ladder tries these one at a time and stops at the first that cleanly separates
# this source from the others, so the most human-readable discriminators (the biological identity, then Date,
# then the technical Replicate/Row) come first.
# Note that Date is the only column whose values must be compared in the rule's normalized form (%d.%m.%y)
# rather than raw (see _disc_series).
_MERGE_DISCRIMINATORS = ['Main_Plasmids', 'Cell_Line', 'Ligand', 'Transfection',
                         'Date', 'Replicate', 'Plate_Row']
# master column -> the criteria key _emit_manual_token / _parse_manual_token use for it.
# This is keyed by MASTER-COLUMN name (what the discriminator ladder iterates over).
# (different from _MANUAL_KEY_MAP as that map is keyed by the blob's ABBREVIATED field
# names ('Cell', 'Cond', 'Rep', 'Row', 'Main', 'File')). Order is irrelevant — this is a
# plain .get() lookup table, never iterated.
_COL_TO_CRIT_KEY = {'Main_Plasmids': 'Main_Plasmids', 'Cell_Line': 'Cell_Line',
                    'Ligand': 'Ligand', 'Transfection': 'Condition', 'Date': 'Date',
                    'Replicate': 'Replicate', 'Plate_Row': 'Row',
                    'File_Name': 'File_Name'}


def _disc_series(df: pd.DataFrame, col: str) -> pd.Series:
    """Per-row string values of `col`, normalized to match how criteria_mask compares it.
    Every column is matched by raw stripped string EXCEPT Date, which criteria_mask compares
    in '%d.%m.%y' form — so a Date discriminator must be computed (and injected) in that same
    form to round-trip through _emit_manual_token -> _parse_manual_token -> criteria_mask."""
    if col == 'Date':
        return (pd.to_datetime(df['Date'], errors='coerce', format='mixed')
                .dt.strftime('%d.%m.%y').fillna(''))
    return df[col].astype(str).str.strip()


def _wells_of_mask(df: pd.DataFrame, mask: pd.Series) -> set:
    """(File_Name, Well_ID) tuples for the rows under `mask`."""
    if not mask.any():
        return set()
    sub = df.loc[mask]
    return set(zip(sub['File_Name'].astype(str), sub['Well_ID'].astype(str)))


def _scoped_resolves_to(ctx: '_ResolveCtx', tokens: list[str]) -> set:
    """Union of (File_Name, Well_ID) that `tokens` re-resolve to against the merged frame
    (criteria ∩ Is_Excluded). Used to VERIFY a scoping attempt before accepting it."""
    acc = set()
    for t in tokens:
        for e in parse_exclusion_blob(t):
            acc |= ctx.rule_wells(e, only_excluded=True)
    return acc


def _inject_manual(label: str, overrides: dict) -> str:
    """Take a manual rule's text, apply {master_column: value} overrides, re-emit."""
    crit = _parse_manual_token(label)
    for col, val in overrides.items():
        crit[_COL_TO_CRIT_KEY.get(col, col)] = val
    return _emit_manual_token(crit)


def _scope_one_manual(merged_df, ctx, entry, origin_mask, w_mask, target_wells, log_fn):
    """Scope a single MANUAL rule against the provisional merged frame. Walks the
    discriminator ladder (no-bleed -> single column -> minimal pair -> File) and returns
    the FIRST attempt that verifies (re-resolves to exactly target_wells). Returns None if
    nothing verified (caller then falls back to per-well tokens)."""
    crit_mask = ctx.criteria_mask(entry["source"], entry["criteria"])
    bleed_mask = crit_mask & (~origin_mask)

    # --- (A) no bleed: the rule already cannot touch other sources -> keep verbatim ---
    if not bleed_mask.any():
        toks = [entry["label"]]
        if _scoped_resolves_to(ctx, toks) == target_wells:
            return toks

    # --- (B) single readable column ---
    for col in _MERGE_DISCRIMINATORS:
        if col not in merged_df.columns:
            continue
        s = _disc_series(merged_df, col)
        vw = set(s[w_mask])
        vb = set(s[bleed_mask])
        if vw and vw.isdisjoint(vb):
            toks = [_inject_manual(entry["label"], {col: v}) for v in sorted(vw)]
            if _scoped_resolves_to(ctx, toks) == target_wells:
                if log_fn and len(toks) == 1:
                    log_fn(f"   [MERGE SCOPE] '{entry['label']}' scoped by {col}.")
                return toks

    # --- (C) minimal pair of columns ---
    present = [c for c in _MERGE_DISCRIMINATORS if c in merged_df.columns]
    series = {c: _disc_series(merged_df, c) for c in present}
    for i in range(len(present)):
        for j in range(i + 1, len(present)):
            c1, c2 = present[i], present[j]
            s1, s2 = series[c1], series[c2]
            pw = set(zip(s1[w_mask], s2[w_mask]))
            pb = set(zip(s1[bleed_mask], s2[bleed_mask]))
            if pw and pw.isdisjoint(pb):
                toks = [_inject_manual(entry["label"], {c1: v1, c2: v2})
                        for (v1, v2) in sorted(pw)]
                if _scoped_resolves_to(ctx, toks) == target_wells:
                    if log_fn:
                        log_fn(f"   [MERGE SCOPE] '{entry['label']}' scoped by {c1}+{c2}.")
                    return toks

    # --- (D) last resort: File scoping (one rule per origin file) ---
    files = sorted({f for (f, _) in target_wells})
    bleed_files = set(merged_df.loc[bleed_mask, 'File_Name'].astype(str).str.strip())
    if files and set(files).isdisjoint(bleed_files):
        toks = [_inject_manual(entry["label"], {'File_Name': f}) for f in files]
        if _scoped_resolves_to(ctx, toks) == target_wells:
            if log_fn:
                log_fn(f"   [MERGE SCOPE] '{entry['label']}' scoped by File ({len(files)}).")
            return toks

    return None


def scope_blob_for_merge(merged_df: pd.DataFrame, src_blob, origin_mask,
                         log_fn=None) -> list[str]:
    """
    Re-scope ONE source's Applied_Exclusions blob for a merge.

    Args:
        merged_df:   the provisional concatenated frame (all sources, collision-resolved).
        src_blob:    this source's Applied_Exclusions blob (File references already remapped
                     by the caller if any of this source's files were collision-renamed).
        origin_mask: boolean Series over merged_df.index, True on THIS source's rows.
        log_fn:      optional logger callback.

    Returns the list of scoped rule tokens for this source such that, when concatenated with
    the other sources' tokens into the unified blob, each rule re-resolves (criteria ∩
    Is_Excluded) against the WHOLE merged frame to exactly the (File_Name, Well_ID) set it
    covered within this source alone.

    Strategy per rule:
      * AUTO rule  -> re-emit one File-pinned well_token per origin well (single-well; exact).
      * manual rule -> inject the most readable discriminating field(s) (Main_Plasmids,
        Cell_Line, Ligand, Transfection, Replicate, Plate_Row), falling back to File scoping,
        then VERIFY against the merged frame.
      * verification failure (pathological shared-everything overlap) -> decompose THAT rule
        only into per-well File-pinned tokens.

    A blob with no parseable markers (opaque/legacy) is re-expressed as per-well tokens for
    this source's currently-excluded wells.
    """
    if merged_df is None or merged_df.empty:
        return []
    if origin_mask is None:
        origin_mask = pd.Series(True, index=merged_df.index)
    origin_mask = origin_mask.reindex(merged_df.index, fill_value=False).astype(bool)

    ctx = _ResolveCtx(merged_df)
    excl = ctx.excluded_mask()
    entries = parse_exclusion_blob(src_blob)

    # Opaque / legacy blob: no rule structure -> per-well tokens for this source's excluded wells.
    if not entries:
        own = _wells_of_mask(merged_df, origin_mask & excl)
        return [well_token(merged_df, f, w) for (f, w) in sorted(own)]

    out = []
    for e in entries:
        crit_mask = ctx.criteria_mask(e["source"], e["criteria"])
        w_mask = crit_mask & origin_mask & excl
        target_wells = _wells_of_mask(merged_df, w_mask)
        if not target_wells:
            # This rule excludes nothing within this source's surviving rows -> drop it.
            continue

        if e["source"] == "AUTO":
            # Single-well auto rules: re-emit File-pinned per origin well (exact).
            out.extend(well_token(merged_df, f, w) for (f, w) in sorted(target_wells))
            continue

        scoped = _scope_one_manual(merged_df, ctx, e, origin_mask, w_mask,
                                   target_wells, log_fn)
        if scoped is not None:
            out.extend(scoped)
        else:
            # Safety net: pathological overlap (sources identical in every meaningful field
            # AND by file) -> decompose this rule alone into per-well File-pinned tokens.
            if log_fn:
                log_fn(f"   [MERGE SCOPE] '{e['label']}' could not be field-scoped; "
                       f"decomposed into {len(target_wells)} per-well token(s).")
            out.extend(well_token(merged_df, f, w) for (f, w) in sorted(target_wells))

    return out


# --------------------------------------------------------------------------- #
# Restore core
# --------------------------------------------------------------------------- #

def _write_blob(master_df: pd.DataFrame, tokens: list[str]):
    """Broadcast the global Applied_Exclusions blob to every row."""
    new_blob = " || ".join(tokens) if tokens else "None"
    master_df['Applied_Exclusions'] = new_blob


def _set_excluded(master_df: pd.DataFrame, wells: set, value: bool):
    """Flip Is_Excluded for a set of (File_Name, Well_ID) tuples (in place).
    value=True for exclusion, value=False for reverting."""
    if not wells:
        return
    key = pd.MultiIndex.from_arrays(
        [master_df['File_Name'].astype(str), master_df['Well_ID'].astype(str)])
    target = pd.Series(key.isin(list(wells)), index=master_df.index)
    master_df.loc[target, 'Is_Excluded'] = value


def _recompute_in_place(master_df: pd.DataFrame, affected_files: list, config):
    """Recompute affected files via the engine and copy results back into master_df.
    Used for single-call restores. Batch revert does this once at the end (— see the recompute flag
    in _restore_core.)"""
    if not affected_files:
        return
    recomputed = recompute_master_after_exclusion(master_df, list(affected_files), config)
    # recompute_master_after_exclusion preserves the row index; copy every column back so
    # the caller's master_df object stays valid (Raw_BRET_unexcluded is untouched).
    for col in recomputed.columns:
        master_df[col] = recomputed[col].values


def _restore_core(master_df, label, wells_subset, config, recompute: bool = True, ctx=None):
    """
    Shared implementation for restore_rule (wells_subset=None -> all of the rule's wells)
    and restore_wells (wells_subset = an explicit subset). Mutates master_df in place
    (Is_Excluded + Applied_Exclusions + recomputed derived columns) and returns a report.

      1. Resolve current rule and the wells it points at
      2. Find every OTHER active rule's wells ("other_cov"): if a well is also wanted by a
         rule it is NOT reverted
      3. Decide each well: covered by another rule -> stays off ("blocked"); can't be
         rebuilt -> stays off ("nonrestorable"); otherwise -> goes in the "restore" pile.
      4. Flip the restore pile back ON, then rewrite the blob: drop this rule's line, or, if
         only some wells came back, replace it with precise one-well lines ("decompose").
      5. Redo the math — UNLESS recompute=False, in which case the caller does it once later.
    The report it returns is the little scoreboard the caller adds up (see `_tally`).

    recompute:
        When True (default, single-call behaviour) the affected files are recomputed in
        place before returning. When False the Is_Excluded flips and the Applied_Exclusions
        blob rewrite still happen, but the (expensive) derived-column recompute is SKIPPED —
        the caller is then responsible for running recompute_master_after_exclusion ONCE
        over the union of every report's "affected_files". This lets a batch revert of many
        rules recompute each affected file a single time instead of once per rule. Deferring
        is safe because no restore decision depends on the derived columns: feasibility
        (is_restorable) reads only the pristine source columns (Raw_BRET_unexcluded /
        Donor / Acceptor), and rule resolution reads only the criteria columns + Is_Excluded,
        all of which are already correct at decision time.

    ctx:
        Optional: pass a shared _ResolveCtx (build_resolve_ctx) for a batch revert (see app.py);
        if omitted, a private one is built for this single call. The ctx must be built on the
        same master_df object being mutated.

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

    if ctx is None:
        ctx = _ResolveCtx(master_df)

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
        target_wells = ctx.all_excluded_wells()
    else:
        target_idx = next((i for i, e in enumerate(entries) if e["label"] == label), None)
        if target_idx is None:
            logger.warning(f"restore: rule label not found: {label!r}")
            return report
        target_entry = entries[target_idx]
        other_entries = [e for j, e in enumerate(entries) if j != target_idx]
        target_wells = ctx.rule_wells(target_entry, only_excluded=True)

    # OTHER coverage computed against the PRE-restore Is_Excluded state.
    other_cov = set()
    for e in other_entries:
        other_cov |= ctx.rule_wells(e, only_excluded=True)

    # --- Determine the subset to attempt restoring ---
    if wells_subset is None:
        subset = set(target_wells)
    else:
        want = {(str(f), str(w)) for f, w in wells_subset}
        subset = want & set(target_wells)

    # --- Per-well decision ---
    # Restorability precomputed once for the whole master (via ctx), then O(1) per well.
    rmap = ctx.restorable_map()
    restored, blocked, nonrestorable = set(), set(), []
    for (f, w) in subset:
        if (f, w) in other_cov:
            blocked.add((f, w))
            continue
        ok, reason = rmap.get((str(f), str(w)),
                              (False, f"well {f}/{w} not found in master"))
        if ok:
            restored.add((f, w))
        else:
            nonrestorable.append((f, w, reason))

    # --- Apply: flip restored wells to not-excluded ---
    _set_excluded(master_df, restored, False)
    ctx.invalidate_excluded()   # Is_Excluded changed -> stale excluded mask

    # --- Decide this rule's rewritten coverage (orphan-free) ---
    # Wells that remain excluded AND are NOT covered by any other active rule must stay
    # attributed by THIS rule. = (target_wells \ restored) \ other_cov.
    remaining_excluded = set(target_wells) - restored
    must_cover = remaining_excluded - other_cov

    if legacy:
        # Opaque blob: re-express the still-excluded wells as explicit per-well tokens.
        all_still = ctx.all_excluded_wells()
        if all_still:
            new_tokens = [well_token(master_df, f, w) for (f, w) in sorted(all_still)]
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
            replacement = [well_token(master_df, f, w) for (f, w) in sorted(must_cover)]
            new_tokens = []
            for j, t in enumerate(tokens):
                if j == target_idx:
                    new_tokens.extend(replacement)
                else:
                    new_tokens.append(t)
            report["rule_decomposed"] = True
        _write_blob(master_df, new_tokens)

    # --- Recompute the files whose Is_Excluded actually changed ---
    # When recompute=False the caller batches a single recompute over the union of all
    # reports' affected_files (see restore_rule / restore_wells docstrings).
    affected = sorted({f for (f, w) in restored})
    if recompute:
        _recompute_in_place(master_df, affected, config)

    report["restored_count"] = len(restored)
    report["nonrestorable"] = sorted(nonrestorable)
    report["blocked_by_other_rule_count"] = len(blocked)
    report["affected_files"] = affected
    return report


def restore_rule(master_df: pd.DataFrame, label: str, config,
                 recompute: bool = True, ctx=None) -> dict:
    """
    Restore ALL of a rule's wells (equivalent to restore_wells over the full set). Per
    well: still matched by ANOTHER active rule -> stays excluded (the other rule keeps
    attributing it); else restorable -> Is_Excluded=False; else (non-restorable, only this
    rule) -> stays excluded. The blob is rewritten so the rule still matches EXACTLY its
    remaining excluded wells (the non-restorable leftovers), decomposed into per-well
    tokens; the token is dropped entirely if zero wells remain. Recomputes affected files.

    Pass recompute=False to defer the derived-column recompute to the caller (batch revert);
    the returned report still lists "affected_files" so the caller can recompute their union
    exactly once. See _restore_core for why deferring is safe.

    Pass a shared ctx (build_resolve_ctx) when reverting many rules in a batch so resolution
    is not recomputed from scratch on every call.

    Report keys: restored_count, nonrestorable [(file, well, reason)],
    blocked_by_other_rule_count, rule_decomposed (bool), affected_files.
    """
    return _restore_core(master_df, label, None, config, recompute=recompute, ctx=ctx)


def restore_wells(master_df: pd.DataFrame, label: str, wells, config,
                  recompute: bool = True, ctx=None) -> dict:
    """
    Per-well refinement over a SUBSET of a rule's wells (mostly for manual criteria-rules;
    auto rules are already single-well). Same per-well decision logic as restore_rule. The
    blob is then rewritten so the rule ALWAYS re-resolves to exactly the wells it still
    excludes: the rule's token is dropped and replaced by well-pinned tokens covering only
    its still-excluded, uniquely-attributed wells (so a broad criteria-rule can never
    re-resolve to a well you just restored). Recomputes affected files.

    Pass recompute=False to defer the derived-column recompute to the caller (batch revert).
    Pass a shared ctx (build_resolve_ctx) for an efficient batch revert.

    `wells` is an iterable of (File_Name, Well_ID).
    Report keys as in restore_rule.
    """
    return _restore_core(master_df, label, list(wells), config, recompute=recompute, ctx=ctx)