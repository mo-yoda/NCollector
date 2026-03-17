import logging
import pandas as pd

logger = logging.getLogger("NCollector")


def ensure_master_csv_schema(df: pd.DataFrame) -> pd.DataFrame:
    """
    Fills in missing columns and cleans legacy data for Master CSVs from older NCollector versions.
    Add new columns to the defaults dict as the schema evolves.
    """
    # --- Fill missing columns with defaults ---
    defaults = {
        "NCollector_version": "< v2",
        "Path": "undocumented path",
    }
    for col, default_val in defaults.items():
        if col not in df.columns:
            df[col] = default_val
            logger.info(f"'{col}' column missing. Assigned '{default_val}'.")

    # --- Fill missing raw donor and acceptor columns with NaN (added in v2+) ---
    nan_defaults = ["Donor_Raw_kinetic", "Acceptor_Raw_kinetic"]
    for col in nan_defaults:
        if col not in df.columns:
            df[col] = float('nan')

    # --- v1.0.0 legacy cleanup ---
    if "Empty/NoID" in df.get('Transfection', pd.Series()).values:
        logger.info("Removing legacy 'Empty/NoID' data...")
        df = df[df['Transfection'] != "Empty/NoID"].copy()

    # --- Reconstruct Bl_LP if missing (pre-v2 beta CSVs only had Bl_Corrected_BRET) ---
    if 'Bl_LP' not in df.columns and 'Bl_Corrected_BRET' in df.columns:
        logger.info("'Bl_LP' missing in CSV. Reconstructing from 'Bl_Corrected_BRET'...")
        df_sorted = df.sort_values(by=['File_Name', 'Well_ID', 'Time_(min)'])
        last_3 = df_sorted.groupby(['File_Name', 'Well_ID']).tail(3)
        bl_lp_means = last_3.groupby(['File_Name', 'Well_ID'])['Bl_Corrected_BRET'].mean().reset_index()
        bl_lp_means.rename(columns={'Bl_Corrected_BRET': 'Bl_LP'}, inplace=True)
        df = df.merge(bl_lp_means, on=['File_Name', 'Well_ID'], how='left')

    # --- Reconstruct Is_Excluded for pre-v2 CSVs ---
    if 'Is_Excluded' not in df.columns:
        has_exclusion_rules = (
            'Applied_Exclusions' in df.columns
            and (df['Applied_Exclusions'].astype(str) != "None").any()
        )
        if has_exclusion_rules and 'Raw_BRET_kinetic' in df.columns:
            # Exclusion rules were applied: wells with NaN in raw BRET data were excluded
            df['Is_Excluded'] = df['Raw_BRET_kinetic'].isna()
            n_excluded = df['Is_Excluded'].sum()
            logger.info(f"Reconstructed 'Is_Excluded' from NaN in Raw_BRET_kinetic ({n_excluded} excluded rows).")
        else:
            df['Is_Excluded'] = False
            logger.info("'Is_Excluded' column missing and no exclusions were applied. Defaulting to False.")

    return df


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

def build_row_info_str(df):
    """
    Creates the standardized list of row info strings:  'Row A: -9.0 log(M) Ligand'
    Used for both GUI population and default export logic.
    """
    if df is None or df.empty:
        return []

    row_series = (
            "Row " + df['Plate_Row'].astype(str) + ": " +
            df['Ligand_Conc'].astype(str) + " log(M) " +
            df['Ligand'].astype(str)
    )

    # Return unique, sorted values
    return sorted(row_series.unique().tolist())

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
        parts.append(df["Ligand_Conc"].astype(str))

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
