import logging
import pandas as pd
from models import LEGACY_COLUMN_DEFAULTS

logger = logging.getLogger("NCollector")


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

def generate_header_key(df, group_by=None):
    """Creates the 'Header_Key' column for exporting data."""
    mapping = {"Cell Line": "Cell_Line", "Transfection": "Transfection"}
    exclude_col = mapping.get(group_by)

    # List of Series to combine
    parts = []
    if exclude_col != "Transfection": parts.append(df["Transfection"])
    if exclude_col != "Cell_Line": parts.append(df["Cell_Line"])
    # Only add ligand, if there is more than one
    if df['Ligand'].nunique() > 1: parts.append(df["Ligand"])

    # Check whether this is AUC data
    is_kinetic = df["Time_(min)"].nunique() > 1
    if is_kinetic:
        conc_display = df["Ligand_Conc"].astype(str)
        # Where Is_Vehicle is False, keep Ligand_Conc, otherwise assign Vehicle
        # pandas where - replace if cond is false (~ is logical NOT in pd)
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
