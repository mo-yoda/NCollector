import logging
import re
import pandas as pd
from models import LEGACY_COLUMN_DEFAULTS

logger = logging.getLogger("NCollector")

# Columns that should be set to NA when retroactively applying exclusions
_EXCLUSION_DATA_COLS = [
    "Raw_BRET_kinetic", "Lab_BRET_kinetic", "Bl_Corrected_BRET",
    "Veh_Norm_Kinetic", "Kinetic_Mean", "Raw_BRET_CRC",
    "Lab_LP", "Bl_LP", "Veh_Norm_LP", "LP_Mean",
    "Lab_AUC", "Bl_AUC", "Veh_Norm_AUC", "AUC_Mean",
]


def _fix_legacy_exclusion_bug(df: pd.DataFrame, tag: str, regex_suffix: str,
                              match_col: str, label: str) -> str | None:
    """
    Shared fixer for legacy exclusion bugs where warnings were recorded in
    Applied_Exclusions but the corresponding data was never set to NA.

    Args:
        df:            Master DataFrame (modified in-place).
        tag:           Warning tag in brackets, e.g. "LOW LUM" or "VEHICLE WARN".
        regex_suffix:  Regex capturing the identifier after the date field
                       (e.g. replicate number or well ID) and the trailing value.
        match_col:     DataFrame column to match the 5th capture group against
                       (e.g. "Replicate" or "Well_ID").
        label:         Human-readable label for the log message (e.g. "Lum check").

    Returns:
        A bug-fix log message string if rows were corrected, or None if nothing to fix.
    """
    if 'Applied_Exclusions' not in df.columns:
        return None

    escaped_tag = re.escape(tag)
    if not df['Applied_Exclusions'].astype(str).str.contains(
            f"AUTO: \\[{escaped_tag}\\]", regex=True).any():
        return None

    # Build regex: 4 common groups (ligand, cell_line, condition, date) + type-specific suffix
    pattern = re.compile(
        rf"AUTO:\s*\[{escaped_tag}\]\s+"
        r"(.+?)\s*\|\s*"  # group 1: ligand
        r"(.+?)\s*\|\s*"  # group 2: cell_line
        r"(.+?)\s*\|\s*"  # group 3: condition (= Transfection)
        r"(.+?)\s*\|\s*"  # group 4: date (dd.mm.yy)
        + regex_suffix     # group 5: type-specific identifier
    )

    sample_text = df['Applied_Exclusions'].astype(str).iloc[0]
    parsed_entries = pattern.findall(sample_text)
    if not parsed_entries:
        return None

    available_data_cols = [c for c in _EXCLUSION_DATA_COLS if c in df.columns]
    if not available_data_cols:
        return None

    # Normalize dates once for comparison (warnings use dd.mm.yy format)
    df_date_normalized = (
        pd.to_datetime(df['Date'], errors='coerce').dt.strftime('%d.%m.%y')
        if 'Date' in df.columns else pd.Series()
    )

    total_fixed = 0

    for ligand, cell_line, condition, date_str, identifier in parsed_entries:
        ligand, cell_line, condition, date_str, identifier = (
            ligand.strip(), cell_line.strip(), condition.strip(), date_str.strip(), identifier.strip())

        # Build mask: 4 common columns + the type-specific column
        mask = pd.Series(True, index=df.index)
        for col, val in [('Ligand', ligand), ('Cell_Line', cell_line),
                         ('Transfection', condition)]:
            if col in df.columns:
                mask &= df[col].astype(str).str.strip() == val
        if 'Date' in df.columns and not df_date_normalized.empty:
            mask &= df_date_normalized == date_str
        if match_col in df.columns:
            mask &= df[match_col].astype(str).str.strip() == identifier

        # Only fix rows where data is NOT already NA
        if mask.any():
            already_na = df.loc[mask, available_data_cols[0]].isna()
            needs_fix = mask & ~already_na.reindex(df.index, fill_value=True)

            if needs_fix.any():
                df.loc[needs_fix, available_data_cols] = float('nan')
                if 'Is_Excluded' in df.columns:
                    df.loc[needs_fix, 'Is_Excluded'] = True
                total_fixed += needs_fix.sum()

    if total_fixed > 0:
        return (f"Fixed legacy {label} bug: set {total_fixed} rows to NaN "
                f"for {len(parsed_entries)} {tag} exclusion(s).")
    return None


def ensure_master_csv_schema(df: pd.DataFrame, log_fn=None) -> tuple[pd.DataFrame, bool, bool]:
    """
    Fills in missing columns and cleans legacy data for Master CSVs from older NCollector versions.
    Add new columns to the defaults dict as the schema evolves.
    Returns (df, was_modified, was_fixed).
    Bug fixes are logged.
    """

    modified = []
    bugs_fixed = []

    def _modified(msg):
        modified.append(msg)
        logger.info(msg)

    def _bug_fix(msg):
        modified.append(msg)
        bugs_fixed.append(msg)
        logger.info(msg)

    # --- Fill all missing columns with defaults ---
    missing_cols = set(LEGACY_COLUMN_DEFAULTS.keys()) - set(df.columns)

    for col in missing_cols:
        config = LEGACY_COLUMN_DEFAULTS[col]
        default_val = config["default"]
        df[col] = default_val

        # Only log if values are not reconstructable from available data
        if not config["reconstructable"]:
            _modified(f"'{col}' column missing. Assigned default: '{default_val}'.")

    # --- v1.0.0 legacy cleanup ---
    if "Empty/NoID" in df.get('Transfection', pd.Series()).values:
        _modified("Removing legacy 'Empty/NoID' data...")
        df = df[df['Transfection'] != "Empty/NoID"].copy()

    # --- Reconstruct Bl_LP if missing (pre-v2 beta CSVs only had Bl_Corrected_BRET) ---
    if 'Bl_LP' not in df.columns and 'Bl_Corrected_BRET' in df.columns:
        _modified("'Bl_LP' missing in CSV. Reconstructing from 'Bl_Corrected_BRET'...")
        df_sorted = df.sort_values(by=['File_Name', 'Well_ID', 'Time_(min)'])
        last_3 = df_sorted.groupby(['File_Name', 'Well_ID']).tail(3)
        bl_lp_means = last_3.groupby(['File_Name', 'Well_ID'])['Bl_Corrected_BRET'].mean().reset_index()
        bl_lp_means.rename(columns={'Bl_Corrected_BRET': 'Bl_LP'}, inplace=True)
        df = df.merge(bl_lp_means, on=['File_Name', 'Well_ID'], how='left')

    # --- Fix legacy exclusion bugs: exclusions were listed as applied but data was not set to NA ---
    _exclusion_fix_configs = [
        {
            "tag": "LOW LUM",
            "regex_suffix": r"Replicate\s+(\d+)\s*-\s*value:\s*[\d.]+",
            "match_col": "Replicate",
            "label": "Lum check",
        },
        {
            "tag": "VEHICLE WARN",
            "regex_suffix": r"([A-H]\d{1,2})\s*-\s*value:\s*[\d.]+",
            "match_col": "Well_ID",
            "label": "Vehicle warning",
        },
    ]

    for fix_cfg in _exclusion_fix_configs:
        fixed = _fix_legacy_exclusion_bug(df, fix_cfg["tag"], fix_cfg["regex_suffix"],
                                          fix_cfg["match_col"], fix_cfg["label"])
        if fixed:
            _bug_fix(fixed)

    # --- Reconstruct Is_Vehicle for older CSVs ---
    if 'Is_Vehicle' in missing_cols:
        if 'Plate_Row' in df.columns and 'Ligand_Conc' in df.columns:
            df['Is_Vehicle'] = (df['Plate_Row'] == 'H') & (df['Ligand_Conc'] == 0.0)
            _modified("'Is_Vehicle' column missing. Reconstructed from Plate_Row == 'H' & Ligand_Conc == 0.0.")
        elif 'Plate_Row' in df.columns:
            df['Is_Vehicle'] = df['Plate_Row'] == 'H'
            _modified("'Is_Vehicle' column missing. Reconstructed from Plate_Row == 'H'.")
        else:
            # Fallback to default
            default_val = LEGACY_COLUMN_DEFAULTS['Is_Vehicle']['default']
            _modified(f"'Is_Vehicle' column could not be reconstructed. Defaulting to '{default_val}'.")

    # --- Clean up legacy vehicle concentration (pre-v2 CSVs stored 0.0 for vehicle) ---
    if 'Ligand_Conc' in df.columns:
        legacy_vehicle = df['Is_Vehicle'] & (df['Ligand_Conc'] == 0.0)
        if legacy_vehicle.any():
            df.loc[legacy_vehicle, 'Ligand_Conc'] = float('nan')
            _modified(f"Cleaned {legacy_vehicle.sum()} vehicle rows: Ligand_Conc 0.0 -> NaN.")

    # --- Reconstruct Is_Excluded for pre-v2 CSVs ---
    if 'Is_Excluded' in missing_cols:
        has_exclusion_rules = (
            'Applied_Exclusions' in df.columns
            and df['Applied_Exclusions'].dropna().astype(str).ne("None").any()
        )
        if has_exclusion_rules and 'Raw_BRET_kinetic' in df.columns:
            # Exclusion rules were applied: wells with NaN in raw BRET data were excluded
            df['Is_Excluded'] = df['Raw_BRET_kinetic'].isna()
            n_excluded = df['Is_Excluded'].sum()
            _modified(f"Reconstructed 'Is_Excluded' from NaN in Raw_BRET_kinetic ({n_excluded} excluded rows).")
        else:
            default_val = LEGACY_COLUMN_DEFAULTS['Is_Excluded']['default']
            _modified(f"'Is_Excluded' column missing and no exclusions were applied. Defaulting to '{default_val}'.")

    # --- Fix legacy AUC bug: excluded wells had 0.0 instead of NaN (sum of all-NaN) ---
    if 'Is_Excluded' in df.columns and df['Is_Excluded'].any():
        excluded_mask = df['Is_Excluded'].astype(bool)
        auc_cols = ['Lab_AUC', 'Bl_AUC', 'Veh_Norm_AUC']
        fixed_cols = []
        for col in auc_cols:
            if col in df.columns:
                bad_zeros = excluded_mask & (df[col] == 0.0)
                if bad_zeros.any():
                    df.loc[bad_zeros, col] = float('nan')
                    fixed_cols.append(col)
        if fixed_cols:
            _bug_fix(f"Fixed legacy AUC bug: set 0.0 -> NaN for excluded wells in {', '.join(fixed_cols)}.")

        # Recalculate AUC_Mean from corrected Veh_Norm_AUC (mean included bogus 0.0 values)
        if 'Veh_Norm_AUC' in fixed_cols and 'AUC_Mean' in df.columns:
            group_keys = ['File_Name', 'Transfection', 'Cell_Line', 'Ligand', 'Plate_Row']
            available_keys = [k for k in group_keys if k in df.columns]
            recalc = df.groupby(available_keys)['Veh_Norm_AUC'].transform('mean')
            df['AUC_Mean'] = recalc
            _bug_fix("Recalculated 'AUC_Mean' from corrected 'Veh_Norm_AUC'.")

    was_modified = len(modified) > 0
    was_fixed = len(bugs_fixed) > 0

    if bugs_fixed and log_fn:
        for msg in bugs_fixed:
            log_fn(f"   [CSV BUG FIX] {msg}")

    return df, was_modified, was_fixed


def apply_export_filters(df, config):
    """Filters the dataframe based on config dictionary."""
    df_subset = df.copy()

    # Mapping config keys to DataFrame columns
    filters = {
        'ligands': 'Ligand',
        'cells': 'Cell_Line',
        'transfections': 'Transfection'
    }

    for cfg_key, df_col in filters.items():
        if config.get(cfg_key) and config.get(cfg_key) != 'All':
            # Ensure target is a list
            targets = config[cfg_key] if isinstance(config[cfg_key], list) else [config[cfg_key]]
            df_subset = df_subset[df_subset[df_col].isin(targets)]
    return df_subset

def build_row_info(df):
    """
    Creates a lookup dict mapping display strings to their filter criteria.
    Display: 'Row A: -9.0 log(M) Ligand' or 'Row H: Vehicle Ligand'
    Value: dict with keys 'row', 'ligand', 'is_vehicle', and optionally 'conc'
    Used for GUI population and export filtering.
    """
    if df is None or df.empty:
        return {}

    lookup = {}
    # Get unique combinations of row info
    group_cols = ['Plate_Row', 'Ligand', 'Is_Vehicle', 'Ligand_Conc']
    available = [c for c in group_cols if c in df.columns]
    unique_rows = df[available].drop_duplicates()

    for _, row in unique_rows.iterrows():
        plate_row = str(row['Plate_Row'])
        ligand = str(row['Ligand'])
        is_vehicle = bool(row['Is_Vehicle'])

        if is_vehicle:
            display = f"Row {plate_row}: Vehicle {ligand}"
        else:
            conc = row['Ligand_Conc']
            display = f"Row {plate_row}: {conc} log(M) {ligand}"

        criteria = {
            'row': plate_row,
            'ligand': ligand,
            'is_vehicle': is_vehicle,
        }
        if not is_vehicle:
            criteria['conc'] = str(row['Ligand_Conc'])

        lookup[display] = criteria

    # Return sorted by key
    return dict(sorted(lookup.items()))

def generate_header_key(df, group_by=None, include_conc=False):
    """
    Creates the 'Header_Key' column for exporting data.
    include_conc: If True, appends concentration (or 'Vehicle') to the key (used for kinetic exports).
    """
    mapping = {"Cell Line": "Cell_Line", "Transfection": "Transfection"}
    exclude_col = mapping.get(group_by)

    # List of Series to combine
    parts = []
    if exclude_col != "Transfection": parts.append(df["Transfection"])
    if exclude_col != "Cell_Line": parts.append(df["Cell_Line"])
    # Only add ligand, if there is more than one
    if df['Ligand'].nunique() > 1: parts.append(df["Ligand"])

    if include_conc:
        conc_display = df["Ligand_Conc"].astype(str)
        conc_display = conc_display.where(~df['Is_Vehicle'].astype(bool), "Vehicle")
        parts.append(conc_display)

    if parts:
        # Start with the first column
        df["Header_Key"] = parts[0]
        # Append subsequent columns
        for p in parts[1:]:
            df["Header_Key"] = df["Header_Key"] + " | " + p
    else:
        # Fallback
        df["Header_Key"] = "Data"
    return df

def create_clean_pivot(df_input, index_col, value_col, disregard_well_id, drop_labeling_control_col):
    """
    Pivots the table. If Mean Data: One column per File.
    If Raw Data: One column per Well (Technical Replicates side-by-side)
    """
    if drop_labeling_control_col:
        df = df_input[df_input["Replicate"] != "labeling control"].copy()
    else:
        df = df_input.copy()

    is_kinetic = df["Time_(min)"].nunique() > 1

    # Zero-padding for correct sorting of two-digit numbers
    well_parts = df['Well_ID'].astype(str).str.extract(r'([A-Za-z])(\d+)', expand=True)

    if not well_parts.empty and well_parts.shape[1] == 2:
        df['well_sort_key'] = well_parts[0] + well_parts[1].str.zfill(2)
    else:
        # Fallback if regex fails
        df['well_sort_key'] = df['Well_ID']

    # Assign a replicate number (1, 2, 3...) per Header_Key
    if disregard_well_id:
        # For mean data sort only by file name
        df['sort_key'] = df['File_Name'].astype(str)
    else:
        # Consider File_Name and Well ID or Plate Col for technical replicates
        if is_kinetic:
            df['sort_key'] = df['File_Name'].astype(str) + "_" + df['well_sort_key'].astype(str)
        else:
            # Only get Col number instead of complete well id
            df['sort_key'] = df['File_Name'].astype(str) + "_" + df['well_sort_key'].astype(str).str[1:].str.zfill(2)
    df['Rep_Num'] = df.groupby('Header_Key')['sort_key'].rank(method='dense').astype(int)

    # Pivot
    pivot = df.pivot_table(
        index=index_col,
        columns=["Header_Key", "Rep_Num"],
        values=value_col
    )

    # Find max replicates and create full grid
    max_reps = df['Rep_Num'].max()
    replicate_range = range(1, max_reps + 1)
    headers = sorted(df["Header_Key"].unique())

    full_columns = pd.MultiIndex.from_product([headers, replicate_range], names=["Header", "Rep"])

    # Reindex adds NaN columns for missing replicates
    pivot = pivot.reindex(columns=full_columns)

    # Flatten Header (Drop the Replicate Number)
    pivot.columns = pivot.columns.droplevel(1)
    pivot.reset_index(drop=True)
    return pivot


def filter_by_conc(df_in, criteria_list):
    """Filters a DataFrame to rows matching the given concentration criteria."""
    filtered = []
    for criteria in criteria_list:
        mask = (
                (df_in["Plate_Row"] == criteria['row']) &
                (df_in["Ligand"] == criteria['ligand'])
        )
        if criteria.get('is_vehicle'):
            mask = mask & (df_in["Is_Vehicle"].astype(bool))
        else:
            mask = mask & (df_in["Ligand_Conc"].astype(str) == criteria['conc'])
        filtered.append(df_in[mask].copy())
    if filtered:
        return pd.concat(filtered).drop_duplicates()
    return pd.DataFrame()


def create_bargraph_table(df_input, value_col, group_by=None, drop_labeling_control=True):
    """
    Creates a table for Prism-style bar charts.
    Columns = conditions (from Header_Key), rows = individual data points (replicates/files).
    Each column lists all replicate values vertically — Prism reads this as grouped bar data.

    Expected input: filtered to a single concentration/row.
    """
    if drop_labeling_control:
        df = df_input[df_input["Replicate"] != "labeling control"].copy()
    else:
        df = df_input.copy()

    # Deduplicate: AUC values are repeated per timepoint in the master
    # For mean data: one value per condition per file (biological replicate)
    # For replicate data: one value per well per file (technical replicate)
    is_mean = "Mean" in value_col
    if is_mean:
        dedup_cols = ["File_Name", "Transfection", "Cell_Line", "Ligand"]
    else:
        dedup_cols = ["File_Name", "Transfection", "Cell_Line", "Well_ID", "Ligand"]
    available = [c for c in dedup_cols if c in df.columns]
    df = df.drop_duplicates(subset=available).copy()

    if value_col not in df.columns:
        return pd.DataFrame()

    # For each Header_Key, collect all values into a list
    groups = {}
    for key in sorted(df["Header_Key"].unique()):
        vals = df.loc[df["Header_Key"] == key, value_col].dropna().tolist()
        groups[key] = vals

    # Pad to same length (Prism expects rectangular tables)
    max_len = max((len(v) for v in groups.values()), default=0)
    for key in groups:
        groups[key] += [float("nan")] * (max_len - len(groups[key]))

    return pd.DataFrame(groups)


def create_heatmap_table(df_input, value_col, group_by, drop_labeling_control=True):
    """
    Creates a table for Prism-style heatmaps.
    Columns = group_by factor (e.g. Cell Line values).
    Rows = merged remaining factors (e.g. "Transfection - Ligand").

    Expected input: filtered to a single concentration/row, mean data.
    The group_by parameter must be "Cell Line" or "Transfection" (not "None").
    """
    if drop_labeling_control:
        df = df_input[df_input["Replicate"] != "labeling control"].copy()
    else:
        df = df_input.copy()

    # Deduplicate: AUC values are repeated per timepoint in the master
    # For mean data: one value per condition per file
    # For replicate data: one value per well per file
    is_mean = "Mean" in value_col
    if is_mean:
        dedup_cols = ["File_Name", "Transfection", "Cell_Line", "Ligand"]
    else:
        dedup_cols = ["File_Name", "Transfection", "Cell_Line", "Well_ID", "Ligand"]
    available = [c for c in dedup_cols if c in df.columns]
    df = df.drop_duplicates(subset=available).copy()

    if value_col not in df.columns:
        return pd.DataFrame()

    # Map group_by display name to column
    group_col_map = {"Cell Line": "Cell_Line", "Transfection": "Transfection"}
    group_col = group_col_map.get(group_by)
    if not group_col or group_col not in df.columns:
        return pd.DataFrame()

    # Build the row factor from remaining factors
    row_parts = []
    if group_col != "Transfection" and "Transfection" in df.columns:
        row_parts.append("Transfection")
    if group_col != "Cell_Line" and "Cell_Line" in df.columns:
        row_parts.append("Cell_Line")
    if df["Ligand"].nunique() > 1 and "Ligand" in df.columns:
        row_parts.append("Ligand")

    if row_parts:
        df["Row_Factor"] = df[row_parts[0]].astype(str)
        for p in row_parts[1:]:
            df["Row_Factor"] = df["Row_Factor"] + " - " + df[p].astype(str)
    else:
        df["Row_Factor"] = "Data"

    # Aggregate: mean per (Row_Factor, group_col) across files
    pivot = df.groupby(["Row_Factor", group_col])[value_col].mean().reset_index()
    pivot = pivot.pivot(index="Row_Factor", columns=group_col, values=value_col)
    pivot.index.name = None
    pivot.columns.name = None

    return pivot
