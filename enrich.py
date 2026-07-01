"""
GUI-agnostic backend for enriching a legacy Master CSV from its source xlsx files.

Reloads the original analysis workbooks to populate Donor_Raw_kinetic,
Acceptor_Raw_kinetic and PR_Time(min) onto an imported master, validating by matching
experimental conditions and Raw_BRET_kinetic data. Pure function of (old master,
source directory) — the folder picker and all GUI refresh live in the app layer.
"""
import os
import logging
import numpy as np
import pandas as pd

from models import (ProcessingConfig, build_plate_layout, ENRICHABLE_COLS,
                    LEGACY_COLUMN_DEFAULTS, APP_VERSION)
from parsing import scan_and_load_folders
from mapping import infer_ligand_info_from_master
from processing import map_plate_metadata, calculate_relative_time

logger = logging.getLogger("NCollector")


def enrich_master(df_old, source_dir, log_fn=None,
                  ask_ligand_choice_fn=None, ask_ligand_layout_fn=None):
    """
    Enrich a legacy Master CSV by reloading source xlsx files to populate
    Donor_Raw_kinetic, Acceptor_Raw_kinetic, and PR_Time(min).
    Validates by matching experimental conditions and Raw_BRET_kinetic data.
    Raw channel data is stored without exclusions for full traceability.

    Returns the enriched master DataFrame on success, or None if enrichment is
    aborted (every abort reason is logged via log_fn). Does not touch any app state.
    ask_ligand_choice_fn / ask_ligand_layout_fn are the GUI dialog fallbacks for
    protocols not resolvable from the old master.
    """
    log = log_fn or (lambda *a, **k: None)

    # --- 2. Load protocols + measurements from path ---
    logger.debug("--- Master Enrichment ---")
    logger.debug("scanning for files...")

    # Determine plate layout from old master (check for labeling control)
    is_labeling = "labeling control" in df_old.get("Replicate", pd.Series()).astype(str).values

    # Extract baseline_end_index from old master: row index where Time_(min) == 0
    baseline_end_idx = None
    if 'Time_(min)' in df_old.columns:
        zero_times = df_old.loc[df_old['Time_(min)'] == 0.0]
        if not zero_times.empty:
            # Get the position within any file (count rows before time==0 for one well)
            sample_file = df_old['File_Name'].iloc[0]
            sample_well = df_old['Well_ID'].iloc[0]
            file_well_mask = (df_old['File_Name'] == sample_file) & (df_old['Well_ID'] == sample_well)
            file_well_times = df_old.loc[file_well_mask, 'Time_(min)'].sort_values()
            zero_idx = (file_well_times == 0.0).values.argmax()
            baseline_end_idx = int(zero_idx)
            logger.debug(f"extracted baseline_end_index={baseline_end_idx} from old master.")

    enrich_config = ProcessingConfig(
        plate_layout=build_plate_layout(is_labeling),
        labeling_correction=is_labeling,
        baseline_end_index=baseline_end_idx,
        ligand_choice_fn=ask_ligand_choice_fn,
        ligand_layout_fn=ask_ligand_layout_fn
    )

    # Get folder paths containing xlsx/xlsm files
    folder_paths = []
    for root, dirs, files in os.walk(source_dir):
        if any(f.endswith(('.xlsx', '.xlsm')) for f in files):
            folder_paths.append(root)

    # Scan and load (debug-level logging — no per-file output to user)
    all_folders = scan_and_load_folders(folder_paths)

    # Only keep folders with both protocol and results
    loaded_folders = [f for f in all_folders if f.protocol and f.results]

    if not loaded_folders:
        log("[ENRICH] ERROR: No valid protocol + measurement pairs found in selected folder.")
        return None

    # Filter to folders matching the Main_Plasmids from the old master
    old_main_plasmids = df_old['Main_Plasmids'].iloc[0] if 'Main_Plasmids' in df_old.columns else None
    if old_main_plasmids and old_main_plasmids != "Unknown":
        matching_folders = []
        for f in loaded_folders:
            folder_mp = " + ".join(f.protocol.main_plasmids) if f.protocol.main_plasmids else "Unknown"
            if folder_mp == old_main_plasmids:
                matching_folders.append(f)
        if matching_folders:
            logger.debug(f"filtered {len(loaded_folders)} folders to "
                         f"{len(matching_folders)} matching Main_Plasmids='{old_main_plasmids}'.")
            loaded_folders = matching_folders
        else:
            log(f"[ENRICH] ERROR: No folders match Main_Plasmids '{old_main_plasmids}'. "
                f"Source files do not match this experiment.")
            return None

    total_files = sum(len(f.results) for f in loaded_folders)
    log(f"[ENRICH] Loaded {total_files} measurement files from "
        f"{len(loaded_folders)} folders. Extracting raw data...")

    # For cases master was created with user-input specified ligand layout:
    # --- 2b. Infer ligand layout and per-plate ligand choices from old master ---
    all_ligand_choices = {}  # file_name -> "L1" or "L2"

    for folder in loaded_folders:
        protocol = folder.protocol
        if not protocol or not protocol.ligand_2:
            continue

        layout, choices = infer_ligand_info_from_master(
            df_old, protocol, len(enrich_config.plate_layout))

        if layout:
            protocol.ligand_layout = layout
            log(f"   [ENRICH] Inferred ligand layout '{layout}' for "
                f"protocol '{protocol.file_name}' from existing master")

        if choices:
            all_ligand_choices.update(choices)
            for fname, choice in choices.items():
                lig_name = (str(protocol.ligand_2) if choice == "L2"
                            else str(protocol.ligand))
                log(f"   [ENRICH] Inferred ligand '{lig_name}' for plate '{fname}' "
                    f"from existing master")

    # Build a ligand_choice_fn that looks up from inferred data,
    # falling back to the dialog for files not found in the old master
    def _enrich_ligand_choice(l1_name, l2_name, plate_info):
        if plate_info in all_ligand_choices:
            return all_ligand_choices[plate_info]
        log(f"   [ENRICH] Plate '{plate_info}' not found in master — asking user")
        return ask_ligand_choice_fn(l1_name, l2_name, plate_info)

    enrich_config.ligand_choice_fn = _enrich_ligand_choice
    # ligand_layout_fn remains as dialog fallback for protocols not in old master

    # --- 3. Map metadata and extract raw data per file ---
    new_file_data = []
    time_col = "Time (min)"
    n_ok = 0
    n_err = 0

    for folder in loaded_folders:
        protocol = folder.protocol
        main_plasmids = " + ".join(protocol.main_plasmids) if protocol.main_plasmids else "Unknown"

        for result in folder.results:
            try:
                # Map conditions onto plate columns
                map_plate_metadata(result, protocol, enrich_config)

                # Calculate relative time vector for this file
                raw_df = result.raw_bret_ratio_df.copy()
                t_vec = calculate_relative_time(
                    raw_df[time_col], enrich_config.baseline_end_index)
                if t_vec is None:
                    t_vec = list(range(len(raw_df)))

                # Persist baseline index for subsequent files
                if enrich_config.baseline_end_index is None and 0 in t_vec:
                    enrich_config.baseline_end_index = t_vec.index(0)

                # PR Time: raw plate reader time mapped by relative time
                raw_time_list = (list(result.raw_time)
                                 if result.raw_time is not None and len(result.raw_time) > 0 else [])
                pr_time_map = (dict(zip(t_vec, raw_time_list))
                               if len(raw_time_list) == len(t_vec) else {})

                # Melt helper — uses calculated Time_(min)
                def melt_df(df_in, val_name):
                    work = df_in.drop(columns=[time_col], errors='ignore').copy()
                    if len(work) == len(t_vec):
                        work.index = t_vec
                    work.index.name = "Time_(min)"
                    return work.reset_index().melt(
                        id_vars="Time_(min)", var_name="Well_ID", value_name=val_name)

                # Raw BRET (for post-merge validation — no exclusions)
                df_raw = melt_df(raw_df, "Raw_BRET_new")
                # Donor and Acceptor (raw, no exclusions — for traceability)
                df_donor = melt_df(result.donor_df, "Donor_Raw_kinetic")
                df_acceptor = melt_df(result.acceptor_df, "Acceptor_Raw_kinetic")

                merge_on = ["Time_(min)", "Well_ID"]
                df_file = df_raw.merge(df_donor, on=merge_on, how="left") \
                                .merge(df_acceptor, on=merge_on, how="left")

                # PR Time
                df_file["PR_Time(min)"] = (df_file["Time_(min)"].map(pr_time_map)
                                           if pr_time_map else float('nan'))

                # Add condition metadata for merge
                df_file["Date"] = result.measurement_date
                df_file["Main_Plasmids"] = main_plasmids
                df_file["Info_Sheet"] = str(result.info_sheet) if result.info_sheet else ""

                meta_maps = {k: {} for k in
                             ['Transfection', 'Cell_Line', 'Ligand']}
                for well_id in df_file['Well_ID'].unique():
                    try:
                        c_idx = int(well_id[1:])
                        meta = result.column_metadata.get(c_idx)
                        if meta:
                            meta_maps['Transfection'][well_id] = meta.condition_name
                            meta_maps['Cell_Line'][well_id] = meta.cell_line
                            meta_maps['Ligand'][well_id] = meta.ligand_identity
                    except:
                        pass

                for col_name, mapping in meta_maps.items():
                    df_file[col_name] = df_file['Well_ID'].map(mapping)

                new_file_data.append(df_file)
                n_ok += 1
                logger.debug(f"Extracted {result.file_name}")

            except Exception as e:
                n_err += 1
                log(f"   [ERROR] Failed to extract {result.file_name}: {e}")

    if not new_file_data:
        log("[ENRICH] ERROR: No data could be extracted. Enrichment aborted.")
        return None

    if n_err > 0:
        log(f"[ENRICH] Extracted {n_ok} files ({n_err} failed).")
    else:
        logger.debug(f"All {n_ok} files extracted successfully.")

    df_new = pd.concat(new_file_data, ignore_index=True)

    # --- 4. Validate condition combinations ---
    logger.debug("Validating conditions...")
    condition_cols = ['Main_Plasmids', 'Transfection', 'Cell_Line', 'Ligand']

    old_combos = set(
        df_old[condition_cols].drop_duplicates().itertuples(index=False, name=None))
    new_combos = set(
        df_new[condition_cols].drop_duplicates().itertuples(index=False, name=None))

    missing_combos = old_combos - new_combos
    if missing_combos:
        log(f"[ENRICH] ERROR: {len(missing_combos)} condition(s) from master CSV "
            f"not found in source files:")
        for combo in list(missing_combos)[:5]:
            log(f"   Missing: {dict(zip(condition_cols, combo))}")
        log("[ENRICH] Enrichment aborted — source files do not match this experiment.")
        return None

    logger.debug(f"All {len(old_combos)} condition combinations found.")

    # --- 5. Validate row count ---
    if len(df_old) != len(df_new):
        log(f"[ENRICH] ERROR: Row count mismatch — master CSV: {len(df_old)}, "
            f"source files: {len(df_new)}. Enrichment aborted.")
        return None

    logger.debug(f"Row count matches ({len(df_old)}).")

    # --- 6. Merge on condition columns + Well_ID + Time_(min) ---
    merge_cols = ['Date', 'Main_Plasmids', 'Transfection', 'Cell_Line', 'Ligand',
                  'Well_ID', 'Time_(min)']

    df_old_work = df_old.copy()

    # Normalize Date to consistent YYYY-MM-DD format
    df_old_work['Date'] = pd.to_datetime(df_old_work['Date']).dt.strftime('%Y-%m-%d')
    df_new['Date'] = pd.to_datetime(df_new['Date']).dt.strftime('%Y-%m-%d')

    # Normalize Time_(min) to float for both sides
    df_old_work['Time_(min)'] = df_old_work['Time_(min)'].astype(float)
    df_new['Time_(min)'] = df_new['Time_(min)'].astype(float)

    # Ensure enrichable columns exist in new data (safety fallback)
    for col in ENRICHABLE_COLS:
        if col not in df_new.columns:
            df_new[col] = LEGACY_COLUMN_DEFAULTS[col]["default"]

    # Select merge keys + enrichable columns + validation column from new data
    validation_col = 'Raw_BRET_new'
    enrich_cols = merge_cols + [validation_col] + ENRICHABLE_COLS
    df_enrich = df_new[enrich_cols].copy()

    # Drop old empty columns before merge
    df_old_work.drop(columns=ENRICHABLE_COLS, inplace=True, errors='ignore')

    df_merged = df_old_work.merge(df_enrich, on=merge_cols, how='left')

    # Check if merge created duplicate rows (many-to-many)
    if len(df_merged) != len(df_old_work):
        log(f"[ENRICH] ERROR: Merge changed row count from {len(df_old_work)} to {len(df_merged)}. "
            f"Likely duplicate merge keys in source data. Enrichment aborted.")
        # Find which keys are duplicated
        dup_keys = df_enrich[df_enrich.duplicated(subset=merge_cols, keep=False)]
        if not dup_keys.empty:
            sample = dup_keys[merge_cols].head(3).to_dict('records')
            logger.debug(f"example duplicate keys: {sample}")
        return None

    # --- 7. Post-merge validation: Raw BRET data should match ---
    # Validation logic: BRET ratio in old master (Raw_BRET_kinetic) is compared to extracted ratio from new path:
    # Difference of the two should be 0 (np.isclose defines decimal tolerance). NAs (excluded values) are ignored.
    if 'Raw_BRET_kinetic' in df_merged.columns and validation_col in df_merged.columns:
        # Compare only non-excluded rows (old master has NaN for excluded wells)
        compare_mask = df_merged['Raw_BRET_kinetic'].notna() & df_merged[validation_col].notna()
        if compare_mask.any():
            old_vals = df_merged.loc[compare_mask, 'Raw_BRET_kinetic'].astype(float).values
            new_vals = df_merged.loc[compare_mask, validation_col].astype(float).values
            diff = ~np.isclose(old_vals, new_vals, rtol=1e-8, atol=1e-8)
            n_mismatched = diff.sum()
            if n_mismatched > 0:
                log(f"[ENRICH] WARNING: {n_mismatched} rows have mismatched Raw BRET values. "
                    f"Source files may not be the originals.")
                # Debug: show mismatched rows
                mismatch_positions = compare_mask[compare_mask].index[diff]
                for idx in mismatch_positions[:5]:
                    row = df_merged.loc[idx]
                    logger.debug(
                        f"Mismatch row {idx}: "
                        f"Well={row.get('Well_ID')} Time={row.get('Time_(min)')} "
                        f"Date={row.get('Date')} Transf={row.get('Transfection')} "
                        f"old={row['Raw_BRET_kinetic']!r} new={row[validation_col]!r} "
                        f"delta={abs(float(row['Raw_BRET_kinetic']) - float(row[validation_col])):.15e}"
                    )
            else:
                logger.debug("Raw BRET validation passed.")
        else:
            log("[ENRICH] WARNING: No overlapping non-NaN Raw BRET data to validate.")

    # Drop the temporary validation column
    df_merged.drop(columns=[validation_col], inplace=True, errors='ignore')

    # --- 8. Final result ---
    n_populated = df_merged['Donor_Raw_kinetic'].notna().sum()
    n_total = len(df_merged)
    logger.info(f"Complete. Populated {n_populated}/{n_total} rows with raw channel data.")

    if n_populated == 0:
        log("[ENRICH] ERROR: No data was populated after merge. Enrichment aborted.")
        return None

    # Update NCollector version to current
    df_merged['NCollector_version'] = APP_VERSION

    return df_merged
