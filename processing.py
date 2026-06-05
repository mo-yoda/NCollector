import logging
import pandas as pd

from models import (PlateColMetadata, PrResult, ProtocolData, ProcessingConfig,
                    BretMetricsBundle, MASTER_COLUMNS, build_plate_layout)
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
                            excluded_wells: list[str],
                            column_metadata: dict | None = None):
    """
    Calculates the mean of the vehicle wells (row H) for each block.

    Labeling-control columns are always disregarded (use of column_metadata).
    In processing pipeline from new imported xlsx, control columns are already dropped before this runs (so this filter is
    a no-op there). When column_metadata is not supplied the behaviour is unchanged.
    """
    vehicle_means = {}
    is_kinetic = bl_corrected_df.shape[0] > 1

    # Columns flagged as labeling controls — their row-H wells are not vehicles.
    control_cols = set()
    if column_metadata:
        control_cols = {
            c for c, meta in column_metadata.items()
            if meta and getattr(meta, "replicate", None) == "labeling control"
        }

    for block in plate_blocks:
        start_col = block[0]
        # Vehicle as row H, never including a labeling-control column.
        wells = [f"H{c}" for c in block if c not in control_cols]

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
        suffix = f"| {well_id} - value: {value_display}"
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

    # Ligand identity and conc map (before transfection map, as it may set protocol.ligand_layout)
    ligand_col_map = get_ligand_map(protocol, len(plate_blocks), config=config, plate_info=result.file_name)
    conc_map_1 = built_conc_dic(protocol.ligand_conc)
    conc_map_2 = built_conc_dic(protocol.ligand_2_conc) if protocol.ligand_2 else {}

    # Get transfection map via mapping ID3 info
    raw_ids = [x.strip() for x in str(result.transfection_id).split(',')] if result.transfection_id else []
    mapped_t_ids = get_transfection_map(
        cell_layout_type=protocol.line_layout,
        ligand_layout_type=protocol.ligand_layout,
        t_ids=raw_ids,
        block_count=len(plate_blocks)
    )
    logger.debug(f"Mapped Block Sequence: {mapped_t_ids}")

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
        else: # Also if t_id does not match protocol
            current_cond_name = "Empty"

        # Assign to all columns in this block
        for rep_idx, col in enumerate(block_cols):
            current_rep_id = str(rep_idx + 1)
            c_line = cl_map.get(col, "Unknown")

            # Columns with unknown cell lines are treated as empty wells
            col_cond_name = "Empty" if str(c_line).startswith("Unknown") else current_cond_name

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
                condition_name=col_cond_name,
                plasmids=current_plasmids,
                ligand_identity=current_ligand_name,
                ligand_conc=current_conc_map,
                replicate=current_rep_id
            )


def check_luminescence(lum_df: pd.DataFrame,
                       column_metadata: dict,
                       lum_threshold: int,
                       date_str: str) -> list[dict]:
    """
    Checks mean luminescence of last 5 timepoints per well against threshold.
    Returns list of warning dicts for wells below threshold.
    """
    warnings = []
    lum_calc_df = lum_df.drop(columns=["Time (min)"], errors='ignore').apply(pd.to_numeric, errors='coerce')

    if lum_calc_df.empty:
        return warnings

    # Per-well mean over the last 5 timepoints (or all available if fewer)
    lum_end = lum_calc_df.iloc[-5:] if len(lum_calc_df) >= 5 else lum_calc_df
    mean_lum_end = lum_end.mean()
    low_lum_wells = mean_lum_end[mean_lum_end < lum_threshold]

    for well_id, val in low_lum_wells.items():
        try:
            row = well_id[0]
            col_idx = int(well_id[1:])
            meta = column_metadata.get(col_idx)
            if meta is None or meta.condition_name == "Empty":
                continue

            warning_dict = create_warning_record(
                warn_type="Lum",
                exp_date=date_str,
                cond_name=meta.condition_name,
                cell_line=meta.cell_line,
                ligand=meta.ligand_identity,
                value=float(val),
                replicate=str(meta.replicate),
                row=row,
                well_id=well_id
            )
            if warning_dict not in warnings:
                warnings.append(warning_dict)
        except Exception as e:
            logger.warning(f"Error processing lum well {well_id}: {e}")

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
                                excluded_wells: list[str],
                                column_metadata: dict | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Normalises kinetic and AUC data to vehicle (row H) means per block.
    Returns (kinetic_df, auc_df).
    """
    veh_means_kinetic = calculate_vehicle_means(bl_corrected_df, plate_blocks, excluded_wells, column_metadata)
    veh_means_auc = calculate_vehicle_means(bl_corr_auc_df, plate_blocks, excluded_wells, column_metadata)

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


# --- Shared Compute Core (pipeline steps 5-10) --- #

def compute_bret_metrics(raw_bret_wide: pd.DataFrame,
                         time_vec,
                         baseline_end_idx: int,
                         plate_blocks: list[range],
                         column_metadata: dict,
                         excluded_wells: list[str],
                         labeling_correction: bool) -> BretMetricsBundle:
    """
    Shared compute core for the data-math portion of the BRET pipeline
    (steps 5-10 of process_bret_measurement, plus the last-3-timepoint / AUC reductions
    that step 11 derives).

    This is a function of the per-well raw BRET ratio and the set of excluded wells:
    it needs no original xlsx and no live PrResult. The exact same math runs whether the
    caller is the object pipeline (process_bret_measurement) or the master-native engine
    (recompute_master_after_exclusion).

    Args:
        raw_bret_wide:  Wide raw BRET ratio. Rows = timepoints (row ORDER must be
                        ascending time; index labels are irrelevant — all math is
                        positional). Columns = Well_IDs. No "Time (min)" column.
        time_vec:       Relative time vector (same length as raw_bret_wide rows). Not used
                        by the math itself (baseline/AUC/last-3 are positional) — accepted
                        for signature symmetry and a length sanity check.
        baseline_end_idx: Positional index of the first post-baseline read (time == 0).
        plate_blocks:   Column blocks (e.g. TRIPLICATE_LAYOUT / QUADRUPLICATE_LAYOUT).
        column_metadata: {col_idx -> PlateColMetadata} for replicate-mean grouping and
                         labeling-control detection.
        excluded_wells:  Well_IDs to NaN out before computing.
        labeling_correction: Whether to run mock-labeling background subtraction.

    Returns:
        BretMetricsBundle with every processed DataFrame the pipeline derives.
    """
    if time_vec is not None and len(time_vec) != len(raw_bret_wide):
        logger.debug(f"compute_bret_metrics: time_vec len {len(time_vec)} != data rows "
                     f"{len(raw_bret_wide)} (proceeding positionally).")

    # NaN excluded wells in raw (idempotent: a no-op if the caller already applied them).
    data_df = raw_bret_wide.copy()
    for well in excluded_wells:
        if well in data_df.columns:
            data_df[well] = float('nan')

    # Cleaned raw is captured BEFORE labeling correction mutates / drops columns.
    cleaned_raw = data_df.copy()

    # 5. LABELING CORRECTION (optional)
    if labeling_correction:
        labeling_corr_df, wells_to_drop = apply_labeling_correction(
            data_df, plate_blocks, column_metadata, excluded_wells
        )
        if wells_to_drop:
            data_df.drop(columns=wells_to_drop, inplace=True, errors='ignore')
    else:
        labeling_corr_df = pd.DataFrame(float('nan'), index=data_df.index, columns=data_df.columns)

    # 6. BASELINE CORRECTION
    bl_corrected_df = apply_baseline_correction(data_df, baseline_end_idx)

    # 7. AUC CALCULATION
    # min_count=1 ensures all-NaN (excluded) wells stay NaN instead of summing to 0.
    if labeling_corr_df.isna().all().all():
        labeling_corr_auc_df = pd.DataFrame(float('nan'), index=data_df.index, columns=data_df.columns)
    else:
        labeling_corr_auc_df = labeling_corr_df.iloc[baseline_end_idx:].sum(min_count=1).to_frame().T
    bl_corr_auc_df = bl_corrected_df.iloc[baseline_end_idx:].sum(min_count=1).to_frame().T

    # 8. VEHICLE NORMALISATION
    kinetic_df, auc_df = apply_vehicle_normalization(
        bl_corrected_df, bl_corr_auc_df, data_df, plate_blocks, excluded_wells, column_metadata
    )

    # 10. MEAN OF REPLICATES
    kinetic_mean_df = calculate_means_on_meta(kinetic_df, plate_blocks, column_metadata)
    auc_mean_df = calculate_means_on_meta(auc_df, plate_blocks, column_metadata)

    # Last-3-timepoint (CRC) reductions
    raw_bret_points_df = cleaned_raw.iloc[-3:].mean().to_frame().T
    labeling_corr_lp_df = labeling_corr_df.iloc[-3:].mean().to_frame().T
    bl_corr_lp_df = bl_corrected_df.iloc[-3:].mean().to_frame().T
    lp_df = kinetic_df.iloc[-3:].mean().to_frame().T
    lp_mean_df = kinetic_mean_df.iloc[-3:].mean().to_frame().T

    return BretMetricsBundle(
        raw_bret_ratio_cleaned=cleaned_raw,
        labeling_corr_kinetic=labeling_corr_df,
        bl_corr_kinetic=bl_corrected_df,
        kinetic_df=kinetic_df,
        kinetic_mean_df=kinetic_mean_df,
        raw_bret_points_df=raw_bret_points_df,
        labeling_corr_lp_df=labeling_corr_lp_df,
        bl_corr_lp_df=bl_corr_lp_df,
        lp_df=lp_df,
        lp_mean_df=lp_mean_df,
        labeling_corr_auc_df=labeling_corr_auc_df,
        bl_corr_auc_df=bl_corr_auc_df,
        auc_df=auc_df,
        auc_mean_df=auc_mean_df,
    )


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

    # 3. CHECK LUMINESCENCE (always performed — donor/acceptor wavelengths are
    #    informational only; the user knows their channels). The 475 nm gate was
    #    removed so the lum check runs for every dataset (incl. imported masters).
    logger.info(f"Performing luminescence check for {result.file_name} "
                f"(donor={result.donor_wavelength} nm, "
                f"acceptor={result.acceptor_wavelength} nm).")
    # Apply exclusions to a temporary copy for the lum check only
    lum_check_df = result.donor_df.copy()
    for well in result.excluded_wells:
        if well in lum_check_df.columns: lum_check_df[well] = float('nan')
    result.low_lum_warnings = check_luminescence(
        lum_check_df, result.column_metadata, config.lum_threshold, date_str
    )

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

    # 5.-10. DATA MATH — delegated to the shared compute core so there is exactly one
    # copy of the labeling -> baseline -> AUC -> vehicle-norm -> replicate-mean pipeline.
    bundle = compute_bret_metrics(
        raw_bret_wide=data_df,
        time_vec=time_vec,
        baseline_end_idx=baseline_end_idx,
        plate_blocks=plate_blocks,
        column_metadata=result.column_metadata,
        excluded_wells=result.excluded_wells,
        labeling_correction=config.labeling_correction,
    )

    # 9. VEHICLE CHECK (kept OUTSIDE the core: it produces warnings, not data)
    result.vehicle_warnings = check_vehicle_wells(
        bundle.kinetic_df, plate_blocks, result.column_metadata,
        result.excluded_wells, config.vehicle_warning_threshold, date_str
    )

    # 11. SAVE RESULTS
    # raw_bret_ratio_cleaned keeps the time-bearing raw_df (unchanged from before);
    # the bundle's cleaned raw is the same data without the "Time (min)" column.
    result.raw_bret_ratio_cleaned = raw_df
    result.time_vector = time_vec
    result.labeling_corr_kinetic = bundle.labeling_corr_kinetic
    result.bl_corr_kinetic = bundle.bl_corr_kinetic
    result.kinetic_df = bundle.kinetic_df
    result.kinetic_mean_df = bundle.kinetic_mean_df

    result.raw_bret_points_df = raw_df.iloc[-3:].mean().to_frame().T
    result.labeling_corr_lp_df = bundle.labeling_corr_lp_df
    result.bl_corr_lp_df = bundle.bl_corr_lp_df
    result.lp_df = bundle.lp_df
    result.lp_mean_df = bundle.lp_mean_df

    result.labeling_corr_auc_df = bundle.labeling_corr_auc_df
    result.bl_corr_auc_df = bundle.bl_corr_auc_df
    result.auc_df = bundle.auc_df
    result.auc_mean_df = bundle.auc_mean_df

    return result

# --- Master-Native Recompute Engine --- #

def _reconstruct_column_metadata(file_rows: pd.DataFrame) -> dict:
    """
    Rebuilds a {col_idx -> PlateColMetadata} dict from a single master file's rows.

    Only the fields the compute core / replicate-mean grouping actually read are
    reconstructed faithfully: condition_name (Transfection), cell_line (Cell_Line),
    ligand_identity (Ligand), replicate (Replicate) and ligand_conc ({Plate_Row:
    Ligand_Conc}). transfection_id / plasmids are not needed by the math.
    """
    column_metadata = {}
    # Column index is the trailing digits of Well_ID (e.g. "H12" -> 12).
    col_idx_series = file_rows['Well_ID'].str[1:].astype(int)

    for col_idx in sorted(col_idx_series.unique()):
        col_mask = col_idx_series == col_idx
        col_block = file_rows[col_mask]

        first = col_block.iloc[0]

        # Per-row ligand concentration map {Plate_Row -> Ligand_Conc}
        conc_map = {}
        if 'Plate_Row' in col_block.columns and 'Ligand_Conc' in col_block.columns:
            for _, r in col_block.iterrows():
                row_char = r['Plate_Row']
                if isinstance(row_char, str) and row_char:
                    conc_map[row_char] = r['Ligand_Conc']

        column_metadata[int(col_idx)] = PlateColMetadata(
            cell_line=first.get('Cell_Line', 'Unknown'),
            condition_name=first.get('Transfection', 'Empty'),
            ligand_identity=first.get('Ligand', 'N/A'),
            replicate=str(first.get('Replicate', '')),
            ligand_conc=conc_map,
        )
    return column_metadata


def coerce_bool(v) -> bool:
    """
    Robustly coerce a value (bool, int/float, or string like 'True'/'False' that a
    re-read CSV may carry) to a Python bool. Used so the master-native code treats
    Is_Excluded consistently regardless of whether the master came from the live
    pipeline (real bool) or from pd.read_csv (possibly object/str).
    """
    if isinstance(v, bool):
        return v
    if isinstance(v, (int, float)):
        try:
            return bool(v) and not pd.isna(v)
        except Exception:
            return bool(v)
    return str(v).strip().lower() in ("true", "1", "1.0", "yes")


def _pivot_wide(file_rows: pd.DataFrame, value_col: str):
    """Pivot a single file's long rows to a wide (Time_(min) x Well_ID) frame.
    Returns None if the column is absent or the pivot is empty."""
    if value_col not in file_rows.columns:
        return None
    wide = (file_rows.pivot_table(index='Time_(min)', columns='Well_ID',
                                  values=value_col, aggfunc='first')
            .sort_index())
    return wide if not wide.empty else None


def reconstruct_file_inputs(master_df: pd.DataFrame, file_name: str):
    """
    Reconstruct, purely from the flat master DataFrame, the per-file inputs needed to
    re-run the luminescence and vehicle quality checks — no original xlsx and no live
    PrResult required. This is the shared reconstruction helper used by both the
    auto-rerun (after an exclusion) and the on-demand Tab-1 re-run buttons.

    Returns a dict (or None if the file has no rows):
        column_metadata : {col_idx -> PlateColMetadata}
        plate_blocks    : triplicate/quadruplicate layout (labeling iff any
                          Replicate == "labeling control")
        excluded_wells  : Well_IDs currently flagged Is_Excluded for this file
        date_str        : measurement date as '%d.%m.%y' (matches object pipeline)
        veh_norm_wide   : wide Veh_Norm_Kinetic (Time_(min) x Well_ID) or None
        donor_wide      : wide Donor_Raw_kinetic (raw, NOT exclusion-applied) or None
        has_donor       : whether any donor data exists for this file
    """
    if master_df is None or master_df.empty or 'File_Name' not in master_df.columns:
        return None

    file_rows = master_df[master_df['File_Name'] == file_name]
    if file_rows.empty:
        return None

    fr = file_rows.copy()
    fr['Time_(min)'] = pd.to_numeric(fr['Time_(min)'], errors='coerce')

    column_metadata = _reconstruct_column_metadata(fr)

    is_labeling = (fr['Replicate'].astype(str) == "labeling control").any()
    plate_blocks = build_plate_layout(is_labeling)

    if 'Is_Excluded' in fr.columns:
        excl_mask = fr['Is_Excluded'].map(coerce_bool)
        excluded_wells = fr.loc[excl_mask, 'Well_ID'].unique().tolist()
    else:
        excluded_wells = []

    try:
        date_str = pd.to_datetime(fr['Date'].iloc[0]).strftime('%d.%m.%y')
    except Exception:
        date_str = str(fr['Date'].iloc[0]) if 'Date' in fr.columns else ""

    veh_norm_wide = _pivot_wide(fr, 'Veh_Norm_Kinetic')

    has_donor = ('Donor_Raw_kinetic' in fr.columns
                 and fr['Donor_Raw_kinetic'].notna().any())
    donor_wide = _pivot_wide(fr, 'Donor_Raw_kinetic') if has_donor else None

    return {
        'column_metadata': column_metadata,
        'plate_blocks': plate_blocks,
        'excluded_wells': excluded_wells,
        'date_str': date_str,
        'veh_norm_wide': veh_norm_wide,
        'donor_wide': donor_wide,
        'has_donor': has_donor,
    }


def _melt_wide(df_wide: pd.DataFrame, value_name: str) -> pd.DataFrame:
    """Melt a time-indexed wide df (cols = Well_IDs) to long (Time_(min), Well_ID, value)."""
    if df_wide is None or df_wide.empty:
        return pd.DataFrame(columns=['Time_(min)', 'Well_ID', value_name])
    work = df_wide.copy()
    work.index = work.index.astype(float)
    work.index.name = 'Time_(min)'
    return work.reset_index().melt(id_vars='Time_(min)', var_name='Well_ID', value_name=value_name)


def recompute_master_after_exclusion(master_df: pd.DataFrame,
                                     affected_files: list[str],
                                     config: ProcessingConfig) -> pd.DataFrame:
    """
    Recomputes the derived BRET metrics for a set of files DIRECTLY on the flat master
    DataFrame — no live PrResult objects and no original xlsx required.

    Exclusions in this app are add-only (monotonic), so recomputing from the master's
    current Raw_BRET_kinetic (NaN-ing newly excluded wells) reproduces exactly what the
    object pipeline would have produced for the retained wells. The recompute starts from
    Raw_BRET_kinetic (the PRE-correction per-well data), NOT from a corrected column, so
    excluding a labeling-control well correctly re-triggers labeling correction.

    For each File_Name in affected_files this:
      1. Pivots that file's Raw_BRET_kinetic to wide (index = sorted unique Time_(min)).
      2. Reconstructs column_metadata from the file's master rows.
      3. Derives plate_blocks (labeling iff any Replicate == "labeling control") and
         baseline_end_idx (positional index of Time_(min) == 0).
      4. Reads excluded_wells (Is_Excluded == True) for the file.
      5. Calls compute_bret_metrics, then overwrites the derived master columns for that
         file's rows. Donor/Acceptor/PR_Time are left untouched.

    Returns a new master DataFrame (input is not mutated), column-ordered per MASTER_COLUMNS.
    """
    if master_df is None or master_df.empty or not affected_files:
        return master_df

    required = {'File_Name', 'Well_ID', 'Time_(min)', 'Raw_BRET_kinetic',
                'Is_Excluded', 'Transfection', 'Cell_Line', 'Ligand',
                'Plate_Row', 'Replicate'}
    missing = required - set(master_df.columns)
    if missing:
        logger.error(f"recompute_master_after_exclusion: master missing columns {missing}. Aborting.")
        return master_df

    df = master_df.copy()
    df['Time_(min)'] = df['Time_(min)'].astype(float)

    # Kinetic (per Well_ID + Time) columns and their bundle source attribute.
    kinetic_cols = {
        'Raw_BRET_kinetic':   'raw_bret_ratio_cleaned',
        'Lab_BRET_kinetic':   'labeling_corr_kinetic',
        'Bl_Corrected_BRET':  'bl_corr_kinetic',
        'Veh_Norm_Kinetic':   'kinetic_df',
    }
    # Single-value (per Well_ID) columns and their bundle source attribute.
    single_cols = {
        'Raw_BRET_CRC':  'raw_bret_points_df',
        'Lab_LP':        'labeling_corr_lp_df',
        'Bl_LP':         'bl_corr_lp_df',
        'Veh_Norm_LP':   'lp_df',
        'Lab_AUC':       'labeling_corr_auc_df',
        'Bl_AUC':        'bl_corr_auc_df',
        'Veh_Norm_AUC':  'auc_df',
    }

    for fname in affected_files:
        file_mask = df['File_Name'] == fname
        file_rows = df.loc[file_mask]
        if file_rows.empty:
            logger.warning(f"recompute: no rows for file '{fname}'. Skipping.")
            continue

        # 1. Pivot clean Raw_BRET_kinetic to wide. Row order ascending time == object-path order.
        raw_wide = (file_rows.pivot_table(index='Time_(min)', columns='Well_ID',
                                          values='Raw_BRET_kinetic', aggfunc='first')
                    .sort_index())
        if raw_wide.empty:
            logger.warning(f"recompute: empty raw pivot for '{fname}'. Skipping.")
            continue

        times_sorted = [float(t) for t in raw_wide.index.tolist()]

        # 2. Reconstruct column metadata.
        column_metadata = _reconstruct_column_metadata(file_rows)

        # 3. Plate layout + baseline index.
        is_labeling = (file_rows['Replicate'].astype(str) == "labeling control").any()
        plate_blocks = build_plate_layout(is_labeling)

        if 0.0 in times_sorted:
            baseline_end_idx = times_sorted.index(0.0)
        elif config.baseline_end_index is not None:
            baseline_end_idx = config.baseline_end_index
            logger.warning(f"recompute: no Time_(min)==0 for '{fname}'. Using config "
                           f"baseline_end_index={baseline_end_idx}.")
        else:
            baseline_end_idx = 0
            logger.warning(f"recompute: no Time_(min)==0 and no config baseline index for "
                           f"'{fname}'. Defaulting baseline_end_idx=0.")

        # 4. Excluded wells for this file.
        excluded_wells = file_rows.loc[file_rows['Is_Excluded'] == True, 'Well_ID'].unique().tolist()

        # 5. Compute.
        bundle = compute_bret_metrics(
            raw_bret_wide=raw_wide,
            time_vec=times_sorted,
            baseline_end_idx=baseline_end_idx,
            plate_blocks=plate_blocks,
            column_metadata=column_metadata,
            excluded_wells=excluded_wells,
            labeling_correction=is_labeling,
        )

        # --- Write back: kinetic columns mapped on (Well_ID, Time_(min)) ---
        fr_idx = file_rows.index
        well_arr = file_rows['Well_ID'].to_numpy()
        time_arr = file_rows['Time_(min)'].to_numpy()
        keys2 = list(zip(well_arr, (float(t) for t in time_arr)))

        for col, attr in kinetic_cols.items():
            wide = getattr(bundle, attr)
            lookup = {}
            if wide is not None and not wide.empty:
                w = wide.copy()
                w.index = [float(t) for t in w.index]
                for t, srow in w.iterrows():
                    for well, val in srow.items():
                        lookup[(well, float(t))] = val
            df.loc[fr_idx, col] = [lookup.get(k, float('nan')) for k in keys2]

        # --- Write back: single-value-per-well columns mapped on Well_ID ---
        for col, attr in single_cols.items():
            src = getattr(bundle, attr)
            val_map = src.iloc[0].to_dict() if (src is not None and not src.empty) else {}
            df.loc[fr_idx, col] = [val_map.get(w, float('nan')) for w in well_arr]

        # --- Write back: replicate-mean columns (same mean-key grouping as compile) ---
        # mean key = "Transfection|Cell_Line|Ligand|Plate_Row"
        well_to_kin_mean = {}
        well_to_lp_mean = {}
        well_to_auc_mean = {}
        for col_idx, meta in column_metadata.items():
            for row_char in "ABCDEFGH":
                well_id = f"{row_char}{col_idx}"
                mean_key = f"{meta.condition_name}|{meta.cell_line}|{meta.ligand_identity}|{row_char}"
                if bundle.kinetic_mean_df is not None and mean_key in bundle.kinetic_mean_df.columns:
                    well_to_kin_mean[well_id] = bundle.kinetic_mean_df[mean_key].tolist()
                if bundle.lp_mean_df is not None and mean_key in bundle.lp_mean_df.columns:
                    well_to_lp_mean[well_id] = bundle.lp_mean_df[mean_key].iloc[0]
                if bundle.auc_mean_df is not None and mean_key in bundle.auc_mean_df.columns:
                    well_to_auc_mean[well_id] = bundle.auc_mean_df[mean_key].iloc[0]

        df.loc[fr_idx, 'LP_Mean'] = [well_to_lp_mean.get(w, float('nan')) for w in well_arr]
        df.loc[fr_idx, 'AUC_Mean'] = [well_to_auc_mean.get(w, float('nan')) for w in well_arr]

        # Kinetic_Mean is per (well, time): align each well's mean series to times_sorted.
        kin_mean_lookup = {}
        for well, series in well_to_kin_mean.items():
            if len(series) == len(times_sorted):
                for t, val in zip(times_sorted, series):
                    kin_mean_lookup[(well, float(t))] = val
        df.loc[fr_idx, 'Kinetic_Mean'] = [kin_mean_lookup.get(k, float('nan')) for k in keys2]

    # Final column ordering (only columns actually present).
    final_cols = [c for c in MASTER_COLUMNS if c in df.columns]
    # Preserve any extra columns that aren't part of the canonical schema, after the schema cols.
    extra = [c for c in df.columns if c not in final_cols]
    return df[final_cols + extra]