import logging
import pandas as pd

from models import PlateColMetadata, PrResult, ProtocolData, ProcessingConfig
from mapping import get_cell_line_map, get_transfection_map, get_ligand_map, built_conc_dic

logger = logging.getLogger("NCollector")

# --- Calculation Helpers --- #

def calculate_relative_time(raw_time_col: pd.Series, baseline_end_idx: None):
    """
    Calculates a relative time vector. First measurement after baseline is set to 0. Negative time for baseline reads.

    If no index for the number of baseline reads is provided, intervals are used to separate baseline and
    kinetic readings (divergent interval is expected for manual ligand addition).
    If all intervals are the same (e.g. ligand addition by injector) returns None.
    """
    # Clean and convert to numeric
    times = pd.Series(pd.to_numeric(raw_time_col, errors='coerce'))

    # Calculate interval
    # -> Rounding is important since otherwise 1.00 and 1.0 are not considered the same
    intervals = times.diff().dropna().round(2)

    if baseline_end_idx is None:
        if len(intervals.unique()) < 2:
            logger.debug(f"Could not detect baseline reads for {raw_time_col}")
            return None

        # Extract the one interval that is different from others (manual ligand addition)
        unique_interval = intervals.drop_duplicates(keep=False)
        if len(unique_interval) != 1:
            logger.error(f"Multiple intervals found for baseline readings: {unique_interval}")
            return None

        # Extract the index (baseline_end_idx)
        baseline_end_idx = unique_interval.index[0]

    if len(times) < baseline_end_idx:
        logger.warning(f"Data has fewer than {baseline_end_idx} rows.")
        return None

    # Define the plate reader read interval
    read_interval = intervals.mode()[0]

    # Generate time vector for baseline and kinetic reading
    n_rows = len(raw_time_col)
    time_vector = []

    for i in range(n_rows):
        # (current_index - zero_index) * interval
        t = (i - baseline_end_idx) * read_interval
        time_vector.append(t)

    return time_vector

def get_block_start_for_col(col_index: int, plate_blocks: list[range]):
    """Finds the start column of the block that contains col_index."""
    for block in plate_blocks:
        if col_index in block:
            return block[0]  # Return the first column of that block (e.g., 1, 4, 7...)
    return None

def calculate_vehicle_means(bl_corrected_df: pd.DataFrame,
                            plate_blocks: list[range],
                            excluded_wells: list[str]):
    """
    Calculates the mean of the vehicle wells (row H) for each block.
    """
    vehicle_means = {}
    is_kinetic = bl_corrected_df.shape[0] > 1

    for block in plate_blocks:
        start_col = block[0]
        # Vehicle as row H
        wells = [f"H{c}" for c in block]

        # Filter for valid wells present in data and not excluded
        valid_vehicles = [w for w in wells if w in bl_corrected_df.columns and w not in excluded_wells]

        if valid_vehicles:
            # Get vehicle values and make sure that data is numeric
            vehicle_data = bl_corrected_df[valid_vehicles].apply(pd.to_numeric, errors='coerce')
            if is_kinetic:
                # Row-wise mean for kinetic traces (result: series of length = timepoints)
                vehicle_means[start_col] = vehicle_data.mean(axis=1)
            else:
                # Scalar mean for AUC (result: single float)
                vehicle_means[start_col] = vehicle_data.mean(axis=1).iloc[0]
        else:
            vehicle_means[start_col] = None

    return vehicle_means

def calculate_means_on_meta(processed_df: pd.DataFrame,
                            plate_blocks: list[range],
                            col_metadata: dict,  # Dic created from PlateColMetadata
                            grouping_mode: str = "row"):
    """
    Calculates the mean of technical replicates. As in processed_df each col is one well,
    the mean is performed of three cols within one block.
    Header format of returned df: "Condition_Name|Cell_Line|Ligand_Name"
    Grouping mode is either
        "row"   Returns mean per Row (A-H) per Block (Condition)
        "col"   Returns mean of the per Replicate (1-x)
    or block.   Returns mean of the ENTIRE Block (all Rows A-H)
    """
    mean_data = {}
    row_labels = list("ABCDEFGH")

    for block_idx, block_cols in enumerate(plate_blocks):
        # Identify the Metadata for this block (use first col as it is identical to others)
        first_col_in_block = block_cols[0]
        meta = col_metadata.get(first_col_in_block)

        # Create a base name for the condition
        cond_name = meta.condition_name if meta else f"Block_{block_idx + 1}"
        cell_line = meta.cell_line if meta else "Unknown"
        ligand_name = meta.ligand_identity

        # Used to calculate means of replicates
        if grouping_mode == "row":
            # Iterate through plate rows (A-H)
            for row in row_labels:
                # Construct well IDs for this specific condition (e.g., A1, A2, A3)
                replicate_wells = [f"{row}{c}" for c in block_cols]

                # Filter for wells that actually exist in the processed dataframe
                valid_wells = [w for w in replicate_wells if w in processed_df.columns]

                if valid_wells:
                    # Select the data for these wells
                    # axis=1 calculates the mean across columns (replicates) per time point
                    # skipna=True is default, handling excluded wells automatically
                    mean_series = processed_df[valid_wells].apply(pd.to_numeric, errors='coerce').mean(axis=1)

                    # Construct a unique column header
                    header_key = f"{cond_name}|{cell_line}|{ligand_name}|{row}"
                    mean_data[header_key] = mean_series
        # Used for luminescence check
        elif grouping_mode == "column":
            for i, col_idx in enumerate(block_cols):
                repl_num = i + 1
                wells = [f"{r}{col_idx}" for r in row_labels]

                valid = [w for w in wells if w in processed_df.columns]

                if valid:
                    val = processed_df[valid].apply(pd.to_numeric, errors='coerce').mean(axis=1)
                    key = f"{cond_name}|{cell_line}|{ligand_name}|{repl_num}"
                    mean_data[key] = val

        elif grouping_mode == "block":
            all_wells_in_block = []
            for row in row_labels:
                all_wells_in_block.extend([f"{row}{c}" for c in block_cols])

            valid_wells = [w for w in all_wells_in_block if w in processed_df.columns]

            if valid_wells:
                mean_series = processed_df[valid_wells].apply(pd.to_numeric, errors='coerce').mean(axis=1)
                header_key = f"{cond_name}|{cell_line}|{ligand_name}"
                mean_data[header_key] = mean_series

    return pd.DataFrame(mean_data)

# --- Warning Helpers --- #

def format_warning_str(warn_type: str,
                       exp_date: str,
                       cond_name: str,
                       cell_line: str,
                       value: float,
                       ligand: str = "",
                       replicate: str = "",
                       well_id: str = ""):
    """
    Creates the warnings str to be displayed in dialogue for lum and vehicle warnings.
    warn_type: 'Lum' or 'Veh'
    """
    if warn_type == "Lum":
        value_display = round(float(value), 1)
        prefix = "[LOW LUM]"
        suffix = f"| Replicate {replicate} - value: {value_display}"
    elif warn_type == "Veh":
        value_display = round(float(value), 3)
        prefix = "[VEHICLE WARN]"
        suffix = f"| {well_id} - value: {value_display}"
    else:
        prefix, suffix = "", ""

    # Build warning string
    warning_str = f"{prefix}   {ligand} | {cell_line} | {cond_name} | {exp_date} {suffix}"

    return warning_str

def create_warning_record(warn_type: str,
                          exp_date: str,
                          cond_name: str,
                          cell_line: str,
                          value: float,
                          ligand: str = "",
                          replicate: str = "",
                          row: str = "",
                          well_id: str = ""):
    """
    Generates the full warning dictionary from metadata.
    Calls format_warning_str internally to generate the 'Display' key.
    """
    display_text = format_warning_str(
        warn_type=warn_type,
        exp_date=exp_date,
        cond_name=cond_name,
        cell_line=cell_line,
        value=value,
        ligand=ligand,
        replicate= replicate,
        well_id=well_id
    )

    # Return the standardized dictionary structure
    return {
        "Ligand": ligand,
        "Date": exp_date,
        "Cell_Line": cell_line,
        "Condition": cond_name,
        "Replicate": replicate,
        "Row": row,
        "Display": display_text
    }

# --- Processing Sub-Steps --- #

def map_plate_metadata(result: PrResult, protocol: ProtocolData, config: ProcessingConfig):
    """
    Maps cell line, transfection, and ligand identity onto each plate column (1-12).
    Populates result.column_metadata with PlateColMetadata for each column.
    """
    plate_blocks = config.plate_layout
    is_labeling = config.labeling_correction

    # Get cell line map
    cl_map = get_cell_line_map(protocol, result.cell_line, len(plate_blocks))

    # Get transfection map via mapping ID3 info
    raw_ids = [x.strip() for x in str(result.transfection_id).split(',')] if result.transfection_id else []
    mapped_t_ids = get_transfection_map(
        cell_layout_type=protocol.line_layout,
        ligand_layout_type=protocol.ligand_layout,
        t_ids=raw_ids,
        block_count=len(plate_blocks)
    )
    logger.debug(f"Mapped Block Sequence: {mapped_t_ids}")

    # Ligand identity and conc map
    ligand_col_map = get_ligand_map(protocol, len(plate_blocks))
    conc_map_1 = built_conc_dic(protocol.ligand_conc)
    conc_map_2 = built_conc_dic(protocol.ligand_2_conc) if protocol.ligand_2 else {}

    # Apply metadata on cols
    for i, block_cols in enumerate(plate_blocks):
        t_id = mapped_t_ids[i]

        # Resolve ID to Name (using Protocol)
        current_plasmids = []
        if t_id in protocol.transfection_conditions:
            current_plasmids = protocol.transfection_conditions[t_id]
            current_cond_name = " + ".join(sorted(current_plasmids))
        elif t_id == "N/A":
            current_cond_name = "Empty"
        else:
            current_cond_name = f"ID {t_id} (Missing)"

        # Assign to all columns in this block
        for rep_idx, col in enumerate(block_cols):
            current_rep_id = str(rep_idx + 1)
            c_line = cl_map.get(col, "Unknown")

            which_lig = ligand_col_map.get(col, 'L1')
            if which_lig == 'L2' and protocol.ligand_2:
                current_ligand_name = str(protocol.ligand_2)
                current_conc_map = conc_map_2
            else:
                current_ligand_name = str(protocol.ligand)
                current_conc_map = conc_map_1

            # In case of labeling correction, use last col in block as mock labeling control
            if is_labeling and rep_idx == (len(block_cols) - 1):
                current_rep_id = "labeling control"

            result.column_metadata[col] = PlateColMetadata(
                cell_line=c_line,
                transfection_id=t_id,
                condition_name=current_cond_name,
                plasmids=current_plasmids,
                ligand_identity=current_ligand_name,
                ligand_conc=current_conc_map,
                replicate=current_rep_id
            )


def check_luminescence(lum_df: pd.DataFrame,
                       plate_blocks: list[range],
                       column_metadata: dict,
                       lum_threshold: int,
                       date_str: str) -> list[dict]:
    """
    Checks mean luminescence of last 5 timepoints per replicate against threshold.
    Returns list of warning dicts for conditions below threshold.
    """
    warnings = []
    lum_calc_df = lum_df.drop(columns=["Time (min)"], errors='ignore').apply(pd.to_numeric, errors='coerce')
    lum_mean_df = calculate_means_on_meta(lum_calc_df, plate_blocks, column_metadata, grouping_mode="column")

    if lum_mean_df.empty:
        return warnings

    lum_end = lum_mean_df.iloc[-5:] if len(lum_mean_df) >= 5 else lum_mean_df
    mean_lum_end = lum_end.mean()
    low_lum_cond = mean_lum_end[mean_lum_end < lum_threshold]

    for key, val in low_lum_cond.items():
        try:
            # "Cond_Name|Cell|Ligand"
            parts = key.split("|")
            if len(parts) < 4: continue
            cond_name, cell_line, lig_name, repl_num = parts[0], parts[1], parts[2], parts[3]
            if cond_name == "Empty": continue

            warning_dict = create_warning_record(
                warn_type="Lum",
                exp_date=date_str,
                cond_name=cond_name,
                cell_line=cell_line,
                ligand=lig_name,
                value=float(val),
                replicate=str(repl_num),
                row = "All"
            )
            if warning_dict not in warnings:
                warnings.append(warning_dict)
        except Exception as e:
            logger.warning(f"Error parsing lum key {key}: {e}")

    return warnings


def apply_labeling_correction(data_df: pd.DataFrame,
                              plate_blocks: list[range],
                              column_metadata: dict,
                              excluded_wells: list[str]) -> tuple[pd.DataFrame, list[str]]:
    """
    Subtracts mock-labeling control background from each block.
    Returns (labeling_corrected_df, wells_to_drop) where wells_to_drop are the
    control wells that should be excluded from subsequent processing steps.
    """
    logger.info("Applying Labeling Correction (Background Subtraction)")
    wells_to_drop = []

    for block in plate_blocks:
        control_col = [
            c for c in block
            if column_metadata.get(c) and column_metadata[c].replicate == "labeling control"
        ]

        control_wells = [f"{r}{control_col[0]}" for r in "ABCDEFGH"]
        valid_controls = [w for w in control_wells if w in data_df.columns and w not in excluded_wells]

        if not valid_controls:
            logger.warning(f"No valid control wells found for labeling control in block {block}. Skipping.")
            continue

        mock_labeling = data_df[valid_controls].mean(axis=1)

        block_wells = []
        for col_idx in block:
            block_wells.extend([f"{r}{col_idx}" for r in "ABCDEFGH"])
        valid_targets = [w for w in block_wells if w in data_df.columns]

        # .sub(bg_vector, axis=0) ensures alignment on the time index
        data_df[valid_targets] = data_df[valid_targets].sub(mock_labeling, axis=0)
        wells_to_drop.extend(control_wells)

    labeling_corr_df = data_df.copy()
    return labeling_corr_df, wells_to_drop


def apply_baseline_correction(data_df: pd.DataFrame, baseline_end_idx: int) -> pd.DataFrame:
    """
    Divides each well by the mean of its baseline phase (first N reads).
    Returns baseline-corrected DataFrame with infinities replaced by NaN.
    """
    baseline_means = data_df.iloc[0:baseline_end_idx].mean()
    bl_corrected_df = data_df / baseline_means
    bl_corrected_df = bl_corrected_df.replace([float('inf'), -float('inf')], float('nan')).infer_objects()
    return bl_corrected_df


def apply_vehicle_normalization(bl_corrected_df: pd.DataFrame,
                                bl_corr_auc_df: pd.DataFrame,
                                data_df: pd.DataFrame,
                                plate_blocks: list[range],
                                excluded_wells: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Normalises kinetic and AUC data to vehicle (row H) means per block.
    Returns (kinetic_df, auc_df).
    """
    veh_means_kinetic = calculate_vehicle_means(bl_corrected_df, plate_blocks, excluded_wells)
    veh_means_auc = calculate_vehicle_means(bl_corr_auc_df, plate_blocks, excluded_wells)

    kinetic_norm_dict = {}
    auc_norm_dict = {}

    for col in data_df.columns:
        col_num = int(col[1:])
        block_start = get_block_start_for_col(col_num, plate_blocks)

        # Normalize Kinetic
        if veh_means_kinetic.get(block_start) is not None:
            kinetic_norm_dict[col] = bl_corrected_df[col] / veh_means_kinetic[block_start]
        else:
            kinetic_norm_dict[col] = float('nan')

        # Normalize AUC
        v_auc = veh_means_auc.get(block_start)
        if v_auc is not None and v_auc != 0:
            auc_norm_dict[col] = bl_corr_auc_df[col] / v_auc
        else:
            auc_norm_dict[col] = float('nan')

    kinetic_df = pd.DataFrame(kinetic_norm_dict, index=data_df.index)
    auc_df = pd.DataFrame(auc_norm_dict, index=[0])
    return kinetic_df, auc_df


def check_vehicle_wells(kinetic_df: pd.DataFrame,
                        plate_blocks: list[range],
                        column_metadata: dict,
                        excluded_wells: list[str],
                        acc_vehicle_range: float,
                        date_str: str) -> list[dict]:
    """
    Checks vehicle wells (row H) for excessive deviation from 1 after normalisation.
    Returns list of warning dicts for outlier vehicle wells.
    """
    logger.debug("Starting vehicle check")
    warnings = []

    for block in plate_blocks:
        for c_idx in block:
            well_id = f"H{c_idx}"
            if well_id in excluded_wells or well_id not in kinetic_df.columns:
                continue

            val_to_check = kinetic_df[well_id].mean()
            deviation = abs(val_to_check - 1)

            if deviation > acc_vehicle_range:
                meta = column_metadata.get(c_idx)
                if meta.condition_name == "Empty": continue

                warning_dict = create_warning_record(
                    warn_type="Veh",
                    exp_date=date_str,
                    cond_name=meta.condition_name,
                    cell_line=meta.cell_line,
                    ligand=meta.ligand_identity,
                    value=float(val_to_check),
                    replicate=str(meta.replicate),
                    row="H",
                    well_id=well_id
                )
                warnings.append(warning_dict)

    return warnings


# --- Main Processing Pipeline --- #

def process_bret_measurement(result: PrResult, protocol: ProtocolData, config: ProcessingConfig):
    """
    Full processing pipeline for a single plate reader measurement file.
    Maps metadata, then runs: exclusions → lum check → labeling correction →
    baseline correction → AUC → vehicle normalisation → vehicle check → replicate means.
    """
    if result.is_excluded:
        result.kinetic_df = None
        result.kinetic_mean_df = None
        result.auc_df = None
        result.auc_mean_df = None
        result.vehicle_warnings = []
        result.low_lum_warnings = []
        return result

    # Reset for re-run
    result.vehicle_warnings = []
    result.low_lum_warnings = []
    result.low_lum_cond = {}
    result.column_metadata = {}

    # Unpack config
    plate_blocks = config.plate_layout
    baseline_end_idx = config.baseline_end_index
    date_str = result.measurement_date.strftime('%d.%m.%y')

    logger.debug(f"=== Processing File: {result.file_name} ===")

    # 1. MAP METADATA onto plate columns
    map_plate_metadata(result, protocol, config)

    # 2. APPLY EXCLUSIONS to raw BRET ratio (donor/acceptor are kept original for traceability)
    raw_df = result.raw_bret_ratio_df.copy()
    if result.excluded_wells:
        for well in result.excluded_wells:
            if well in raw_df.columns: raw_df[well] = float('nan')

    # 3. CHECK LUMINESCENCE (only when donor channel is 475 nm)
    if result.donor_wavelength == 475:
        # Apply exclusions to a temporary copy for the lum check only
        lum_check_df = result.donor_df.copy()
        for well in result.excluded_wells:
            if well in lum_check_df.columns: lum_check_df[well] = float('nan')
        result.low_lum_warnings = check_luminescence(
            lum_check_df, plate_blocks, result.column_metadata, config.lum_threshold, date_str
        )
    else:
        logger.debug(f"Lum check skipped: donor wavelength is {result.donor_wavelength}, not 475.")

    # 4. BUILD TIME VECTOR
    time_col = "Time (min)"

    time_vec = calculate_relative_time(raw_df[time_col], baseline_end_idx)

    # If auto-detection failed and no manual index was set, ask user
    if time_vec is None and baseline_end_idx is None and config.user_input_fn:
        user_baseline = config.user_input_fn(
            title="Baseline Detection Failed",
            message=(f"Could not automatically detect baseline reads for "
                     f"'{result.file_name}'.\n\n"
                     f"Enter the number of baseline reads "
                     f"(measurement cycles before ligand addition):"),
            input_type="int"
        )
        if user_baseline is not None:
            baseline_end_idx = user_baseline
            time_vec = calculate_relative_time(raw_df[time_col], baseline_end_idx)

    if time_vec is None:
        logger.warning(f"Time vector could not be built for {result.file_name}. Using raw index as fallback.")
        time_vec = range(len(raw_df))

    # Update baseline_end_idx from built time_vec
    baseline_end_idx = time_vec.index(0)

    # Persist in config so subsequent files and re-runs reuse the same index
    if config.baseline_end_index is None:
        config.baseline_end_index = baseline_end_idx

    data_df = raw_df.drop(columns=[time_col]).copy()

    # 5. LABELING CORRECTION (optional)
    if config.labeling_correction:
        labeling_corr_df, wells_to_drop = apply_labeling_correction(
            data_df, plate_blocks, result.column_metadata, result.excluded_wells
        )
        if wells_to_drop:
            data_df.drop(columns=wells_to_drop, inplace=True, errors='ignore')
    else:
        labeling_corr_df = pd.DataFrame(float('nan'), index=data_df.index, columns=data_df.columns)

    # 6. BASELINE CORRECTION
    bl_corrected_df = apply_baseline_correction(data_df, baseline_end_idx)

    # 7. AUC CALCULATION
    # min_count=1 for at least 1 value needed, otherwise return nan -> ensures excluded wells with all-NaN stay NaN
    if labeling_corr_df.isna().all().all():
        labeling_corr_auc_df = pd.DataFrame(float('nan'), index=data_df.index, columns=data_df.columns)
    else:
        labeling_corr_auc_df = labeling_corr_df.iloc[baseline_end_idx:].sum(min_count=1).to_frame().T
    bl_corr_auc_df = bl_corrected_df.iloc[baseline_end_idx:].sum(min_count=1).to_frame().T

    # 8. VEHICLE NORMALISATION
    kinetic_df, auc_df = apply_vehicle_normalization(
        bl_corrected_df, bl_corr_auc_df, data_df, plate_blocks, result.excluded_wells
    )

    # 9. VEHICLE CHECK
    result.vehicle_warnings = check_vehicle_wells(
        kinetic_df, plate_blocks, result.column_metadata,
        result.excluded_wells, config.vehicle_warning_threshold, date_str
    )

    # 10. MEAN OF REPLICATES
    kinetic_mean_df = calculate_means_on_meta(kinetic_df, plate_blocks, result.column_metadata)
    auc_mean_df = calculate_means_on_meta(auc_df, plate_blocks, result.column_metadata)

    # 11. SAVE RESULTS
    result.raw_bret_ratio_cleaned = raw_df
    result.time_vector = time_vec
    result.labeling_corr_kinetic = labeling_corr_df
    result.bl_corr_kinetic = bl_corrected_df
    result.kinetic_df = kinetic_df
    result.kinetic_mean_df = kinetic_mean_df

    result.raw_bret_points_df = raw_df.iloc[-3:].mean().to_frame().T
    result.labeling_corr_lp_df = labeling_corr_df.iloc[-3:].mean().to_frame().T
    result.bl_corr_lp_df = bl_corrected_df.iloc[-3:].mean().to_frame().T
    result.lp_df = kinetic_df.iloc[-3:].mean().to_frame().T
    result.lp_mean_df = kinetic_mean_df.iloc[-3:].mean().to_frame().T

    result.labeling_corr_auc_df = labeling_corr_auc_df
    result.bl_corr_auc_df = bl_corr_auc_df
    result.auc_df = auc_df
    result.auc_mean_df = auc_mean_df

    return result