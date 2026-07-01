"""
Turns the live object pipeline into the flat master DataFrame.

compile_master_dataframe melts each processed PrResult (kinetic / CRC / AUC frames)
into the long "one row per well per timepoint" master schema (MASTER_COLUMNS).
process_folder_to_master is the Merge-tab helper that scans + processes a folder into
a master-shaped frame without touching application state. Both are GUI-agnostic; the
app layer supplies application state (experiment/directory/rule history), a log
callback, and the interactive main-plasmids selection.
"""
import os
import logging
import pandas as pd

from models import APP_VERSION, MASTER_COLUMNS, MAIN_ONLY_CONDITION
from parsing import scan_and_load_folders
from processing import process_bret_measurement

logger = logging.getLogger("NCollector")


def compile_master_dataframe(experiment, directory, rule_history_text, log_fn=None):
    """
    Compiles all processing steps into one Master DataFrame.
    Structure: 1 row per well per timepoint.
    Means are repeated for respective technical replicates as AUCs for all timepoints.
    Empty wells (unknwon cell line or empty condition) are dropped.

    Pure function of the (already-processed) experiment, its source directory and the
    running rule-history text. The app layer resolves those from application state and
    passes a log callback.
    """
    log = log_fn or (lambda *a, **k: None)

    if not experiment:
        return None
    log("\n--- Building Master CSV ---")

    all_files_data = []

    for folder in experiment:
        if folder.protocol.main_plasmids:
            main_plasmids = " + ".join(folder.protocol.main_plasmids)
        else:
            main_plasmids = "Unknown"

        for res in folder.results:
            if res.is_excluded: continue

            # --- PREPARE KINETIC DATA ---
            # Use pandas melt function to prepare each df from wide to long format
            def melt_df(df, val_name, time_vec):
                if df is None or df.empty: return pd.DataFrame()

                df_work = df.copy()
                # Check lengths
                if len(df_work) != len(time_vec):
                    logger.warning(
                        f"Length mismatch in {res.file_name}: Data {len(df_work)} vs Time {len(time_vec)}")
                    # Use generic index
                    df_work.index.name = "Time_Idx"
                    id_var = "Time_Idx"
                else:
                    # Set the Time Vector as the Index
                    df_work.index = time_vec
                    df_work.index.name = "Time_(min)"
                    id_var = "Time_(min)"

                return df_work.reset_index().melt(
                    id_vars=id_var,
                    var_name="Well_ID",
                    value_name=val_name)

            # Get the Time Vector for this file
            t_vec = res.time_vector
            if not t_vec:
                # Fallback if time vector calculation failed
                t_vec = range(len(res.raw_bret_ratio_cleaned))

            # Ignore time col in raw bret df
            raw_clean = res.raw_bret_ratio_cleaned.drop(columns=["Time (min)"], errors='ignore')

            # df_og = melt_df(res.raw_bret_ratio_df.drop(columns=["Time (min)"], errors='ignore'),
            #                "OG_BRET_ratio", t_vec)
            df_donor = melt_df(res.donor_df.drop(columns=["Time (min)"], errors='ignore'),
                            "Donor_Raw_kinetic", t_vec)
            df_acceptor = melt_df(res.acceptor_df.drop(columns=["Time (min)"], errors='ignore'),
                            "Acceptor_Raw_kinetic", t_vec)
            df_raw = melt_df(raw_clean, "Raw_BRET_kinetic", t_vec)
            df_lab = melt_df(res.labeling_corr_kinetic, "Lab_BRET_kinetic", t_vec)
            df_bl = melt_df(res.bl_corr_kinetic, "Bl_Corrected_BRET", t_vec)
            df_norm = melt_df(res.kinetic_df, "Veh_Norm_Kinetic", t_vec)

            # Merge on [Time_(min), Well_ID]
            merge_on = [df_donor.columns[0], "Well_ID"]

            merged_df = df_donor.merge(df_acceptor, on=merge_on, how="left") \
                .merge(df_raw, on=merge_on, how="left") \
                .merge(df_lab, on=merge_on, how="left") \
                .merge(df_bl, on=merge_on, how="left") \
                .merge(df_norm, on=merge_on, how="left")

            # --- PRISTINE, EXCLUSION-FREE RAW (uniform provenance for ALL files) ---
            # Raw_BRET_unexcluded = Acceptor / Donor,
            # Donor 0/NaN -> NaN
            # Left UNROUNDED. This column is never NaN-d by exclusion, so it is the
            # pristine source the engine restores from.
            _donor = pd.to_numeric(merged_df["Donor_Raw_kinetic"], errors="coerce")
            _acceptor = pd.to_numeric(merged_df["Acceptor_Raw_kinetic"], errors="coerce")
            merged_df["Raw_BRET_unexcluded"] = (
                _acceptor / _donor.where((_donor != 0) & _donor.notna())
            )

            # --- MAP RAW PLATE READER TIME (file-specific, from res.raw_time) ---
            if res.raw_time is not None and len(res.raw_time) == len(t_vec):
                raw_time_map = dict(zip(t_vec, res.raw_time))
                merged_df["PR_Time(min)"] = merged_df["Time_(min)"].map(raw_time_map)
            else:
                merged_df["PR_Time(min)"] = float('nan')

            # --- MAP DATA FOR CRC ---
            # Last 3x points/AUC is 1 value per well, map data to Well_ID
            # Last 3 time points (lp)
            raw_bret_map = res.raw_bret_points_df.iloc[0].to_dict() if res.raw_bret_points_df is not None else {}
            lp_lab_map = res.labeling_corr_lp_df.iloc[0].to_dict() if res.labeling_corr_lp_df is not None else {}
            lp_bl_map = res.bl_corr_lp_df.iloc[0].to_dict() if res.bl_corr_lp_df is not None else {}
            lp_norm_map = res.lp_df.iloc[0].to_dict() if res.lp_df is not None else {}

            # AUC
            auc_lab_map = res.labeling_corr_auc_df.iloc[0].to_dict() if res.labeling_corr_auc_df is not None else {}
            auc_bl_map = res.bl_corr_auc_df.iloc[0].to_dict() if res.bl_corr_auc_df is not None else {}
            auc_norm_map = res.auc_df.iloc[0].to_dict() if res.auc_df is not None else {}

            merged_df['Raw_BRET_CRC'] = merged_df['Well_ID'].map(raw_bret_map)
            merged_df['Lab_LP'] = merged_df['Well_ID'].map(lp_lab_map)
            merged_df['Bl_LP'] = merged_df['Well_ID'].map(lp_bl_map)
            merged_df['Veh_Norm_LP'] = merged_df['Well_ID'].map(lp_norm_map)

            merged_df['Lab_AUC'] = merged_df['Well_ID'].map(auc_lab_map)
            merged_df['Bl_AUC'] = merged_df['Well_ID'].map(auc_bl_map)
            merged_df['Veh_Norm_AUC'] = merged_df['Well_ID'].map(auc_norm_map)

            # --- PREPARE MEAN KINETIC AND AUC MAPPING ---
            well_to_mean_map = {}
            well_to_lp_mean_map = {}
            well_to_auc_mean_map = {}
            for col_idx, meta in res.column_metadata.items():
                col_str = str(col_idx)
                for row_char in "ABCDEFGH":
                    well_id = f"{row_char}{col_str}"

                    # Construct Key: "Condition|Cell|Ligand|Row"
                    mean_key = f"{meta.condition_name}|{meta.cell_line}|{meta.ligand_identity}|{row_char}"

                    # Grab Kinetic Mean Series
                    if res.kinetic_mean_df is not None and mean_key in res.kinetic_mean_df.columns:
                        well_to_mean_map[well_id] = res.kinetic_mean_df[mean_key].tolist()

                    # Grab Last points Mean Value
                    if res.lp_mean_df is not None and mean_key in res.lp_mean_df.columns:
                        well_to_lp_mean_map[well_id] = res.lp_mean_df[mean_key].iloc[0]

                    # Grab AUC Mean Value
                    if res.auc_mean_df is not None and mean_key in res.auc_mean_df.columns:
                        well_to_auc_mean_map[well_id] = res.auc_mean_df[mean_key].iloc[0]

            merged_df['LP_Mean'] = merged_df['Well_ID'].map(well_to_lp_mean_map)
            merged_df['AUC_Mean'] = merged_df['Well_ID'].map(well_to_auc_mean_map)

            # --- MAP KINETIC MEANS ---
            # Create  specialized DF to merge accurately by Time
            mean_rows = []
            # Use the same t_vec defined above
            for well_id, mean_series in well_to_mean_map.items():
                # Ensure mean series matches time vector length
                if len(mean_series) == len(t_vec):
                    for t, val in zip(t_vec, mean_series):
                        mean_rows.append({
                            'Well_ID': well_id,
                            'Time_(min)': t,  # Using actual time for merge key
                            'Kinetic_Mean': val
                        })
            if mean_rows:
                df_means = pd.DataFrame(mean_rows)
                # Merge on Time and Well
                merge_keys = ['Well_ID', 'Time_(min)']
                merged_df = merged_df.merge(df_means, on=merge_keys, how='left')
            else:
                merged_df['Kinetic_Mean'] = float('nan')

            all_files_data.append(merged_df)

            # --- ADD METADATA ---
            merged_df["NCollector_version"] = APP_VERSION
            merged_df["Path"] = directory
            merged_df["Info_Sheet"] = str(res.info_sheet) if res.info_sheet else ""
            merged_df["File_Name"] = res.file_name
            merged_df["Date"] = res.measurement_date
            merged_df["Main_Plasmids"] = main_plasmids

            # Get the exclusion text (handle empty case)
            exclusion_text = rule_history_text if rule_history_text else "None"
            # v2.0.5: join rules with " || "
            exclusion_text_clean = exclusion_text.replace("\n", " || ")
            merged_df["Applied_Exclusions"] = exclusion_text_clean

            # Mark excluded wells: True if this well was in the exclusion list
            excluded_set = set(res.excluded_wells)
            merged_df["Is_Excluded"] = merged_df["Well_ID"].isin(excluded_set)

            # Meta Lookups (Optimization: Build dicts once per file)
            meta_lookups = {'Transfection': {},
                            'Cell_Line': {},
                            'Ligand': {},
                            'Ligand_Conc': {},
                            'Plate_Row': {},
                            'Replicate': {}}

            for well_id in merged_df['Well_ID'].unique():
                try:
                    c_idx = int(well_id[1:])
                    row_char = well_id[0]
                    meta = res.column_metadata.get(c_idx)
                    if meta:
                        meta_lookups['Transfection'][well_id] = meta.condition_name
                        meta_lookups['Cell_Line'][well_id] = meta.cell_line
                        meta_lookups['Ligand'][well_id] = meta.ligand_identity
                        meta_lookups['Ligand_Conc'][well_id] = meta.ligand_conc.get(
                            row_char, float('nan'))
                        meta_lookups['Plate_Row'][well_id] = row_char
                        meta_lookups['Replicate'][well_id] = meta.replicate
                except: pass

            merged_df['Transfection'] = merged_df['Well_ID'].map(meta_lookups['Transfection'])

            # label main_plasmid-only wells
            _blank_cond = (merged_df['Transfection'].isna()
                           | merged_df['Transfection'].astype(str).str.strip().isin(["", "nan", "None"]))
            if _blank_cond.any():
                merged_df.loc[_blank_cond, 'Transfection'] = MAIN_ONLY_CONDITION
            merged_df['Cell_Line'] = merged_df['Well_ID'].map(meta_lookups['Cell_Line'])
            merged_df['Ligand'] = merged_df['Well_ID'].map(meta_lookups['Ligand'])
            merged_df['Ligand_Conc'] = merged_df['Well_ID'].map(meta_lookups['Ligand_Conc'])
            merged_df['Plate_Row'] = merged_df['Well_ID'].map(meta_lookups['Plate_Row'])
            merged_df['Replicate'] = merged_df['Well_ID'].map(meta_lookups['Replicate'])
            merged_df['Is_Vehicle'] = (
                (merged_df['Plate_Row'] == 'H') & (merged_df['Ligand_Conc'].isna())
            )

            # Drop empty wells (Empty transfections / Unknown cell lines)
            merged_df.drop(
                merged_df[
                    merged_df['Transfection'].astype(str).str.contains("Empty", na=False) |
                    merged_df['Cell_Line'].astype(str).str.startswith("Unknown", na=False)
                ].index, inplace=True
            )

    if not all_files_data:
        return pd.DataFrame()

    # Combine all files
    master_df = pd.concat(all_files_data, ignore_index=True)
    # Cleanup columns
    final_cols = [c for c in MASTER_COLUMNS if c in master_df.columns]
    return master_df[final_cols]


def process_folder_to_master(folder_paths, config, select_main_plasmids_fn, log_fn=None):
    """
    Process one or more experiment subfolders into a master-shaped DataFrame WITHOUT
    touching application state. Used by the Merge tab to add an experiment folder as a
    merge source.

    Mirrors the object pipeline (scan -> process_bret_measurement -> compile) but
    writes to a LOCAL frame and returns it. It DOES run the interactive main-plasmids
    selection via select_main_plasmids_fn(experiment) -> (experiment | None, str), but
    purely on the local experiment list. Ligand dialogs are reused via the passed
    ProcessingConfig callbacks. A folder source is enriched by construction (Donor/
    Acceptor raw channels are always populated by compile), so it always passes the
    merge channel gate.

    Returns the compiled master-shaped DataFrame, an empty DataFrame if there is nothing
    to process, or ``None`` if the user cancels the main-plasmids selection.
    """
    log = log_fn or (lambda *a, **k: None)

    if not folder_paths:
        return pd.DataFrame()

    experiment = scan_and_load_folders(folder_paths, log_fn=log_fn)
    if not experiment:
        return pd.DataFrame()

    # Prompt user to pick main plasmids, folder spans more than one main_plasmids set
    selected_experiment, _selected_str = select_main_plasmids_fn(experiment)
    if selected_experiment is None:
        log("[MERGE] Folder add cancelled at main-plasmids selection.")
        return None
    experiment = selected_experiment

    for folder in experiment:
        if not folder.protocol:
            continue
        for result in folder.results:
            process_bret_measurement(result, folder.protocol, config)

    # Provenance directory for the Path column (common root of the scanned subfolders).
    try:
        directory = (folder_paths[0] if len(folder_paths) == 1
                     else os.path.commonpath(folder_paths))
    except ValueError:
        directory = folder_paths[0]

    return compile_master_dataframe(experiment, directory, "", log_fn=log_fn)
