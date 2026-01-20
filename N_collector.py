import os
import tkinter as tk
from tkinter import filedialog, ttk
import pandas as pd
from datetime import datetime, date
from dataclasses import dataclass, field

# --- Dataclass Definition --- #

@dataclass
class PlateColMetadata:
    """Identity of a specific column"""
    cell_line: str = "Unknown"
    transfection_id: str  = "N/A"
    condition_name: str = "Empty"
    plasmids: list[str] = field(default_factory=list)

@dataclass
class PrResult:
    """ Information from a single _analysis file """
    # Information from analysis xlsx itself
    file_name: str
    measurement_date: date
    cell_line: str # ID2
    transfection_id: str # ID3
    raw_bret_ratio_df: pd.DataFrame

    # Connection to protocol file
    # Key = Column Index (1-12), Value = WellMetadata object
    column_metadata: dict[int, PlateColMetadata] = field(default_factory=dict)
    # Stores baseline- and vehicle-normalised BRET ratios
    processed_df: pd.DataFrame | None = None
    kinetic_mean_df: pd.DataFrame | None = None
    auc_mean_df:pd.DataFrame | None = None

    # User interaction for optional exclusion
    is_excluded: bool = False
    excluded_wells: list[str] = field(default_factory=list)

    # Internal check and warnings for helping outlier identification
    vehicle_outliers: dict[str, float] = field(default_factory=dict) # well, value
    warnings : list[str] = field(default_factory=list)

@dataclass
class ProtocolData:
    """Information from a protocol file"""
    file_name: str
    exp_date: date
    n: int
    cell_lines: list[str]
    line_layout: str
    transfection_scheme: pd.DataFrame
    main_plasmids: list[str] # Plasmids transfected in all conditions
    # keys are col strings, values are transfected plasmids
    transfection_conditions: dict[str, list[str]]
    ligand: str
    ligand_conc: pd.DataFrame
    # Second ligand is optional; by | None = None
    ligand_2: str | None = None
    ligand_conc_2: pd.DataFrame | None = None

@dataclass
class MeasurementFolder:
    """A subfolder containing one protocol and multiple result files"""
    folder_name: str
    folder_path: str
    measurement_date: date
    protocol: ProtocolData | None = None # MeasurementFolder is initiated before protocol data is loaded
    results: list[PrResult] = field(default_factory=list) # The default_factory=list initiates this with an empty list
    skipped_files: list[str] = field(default_factory=list)

# --- Tool Functions --- #

def get_location(df: pd.DataFrame, marker: str, match_index: int = 0):
    """
    Finds exact coordinates (row_idx, col_idx) of a text marker in the DataFrame.
    match_index: 0 for 1st occurrence, 1 for 2nd, etc.
    """
    # Create a boolean mask of the whole sheet to search the entire sheet
    # .stack() turns the 2D grid into a 1D Series with double index (row, col)
    mask = df.astype(str).apply(lambda x: x.str.contains(marker, na=False, regex=False))

    # Extract coordinates of all True values, again using .stack()
    all_matches = mask.stack()[mask.stack()].index.tolist()

    if not all_matches:
        print(f"   [ERROR] Marker '{marker}' not found in sheet.")
        return None

    # Check if the requested occurrence exists
    if match_index >= len(all_matches):
        print(
            f"   [ERROR] Requested occurrence #{match_index + 1} of '{marker}' not found. Only {len(all_matches)} found.")
        return None

    # Return the specific (row, col) tuple
    return all_matches[match_index]

def extract_value(
        df: pd.DataFrame,
        marker: str,
        col_offset: int = 1,
        row_offset: int = 0,
        match_index: int = 0
):
    """
    Finds a marker and returns the value of the cell at a relative position.
    match_index: 0 for 1st occurrence, 1 for 2nd, etc.
    Default: Returns the value in the cell immediately to the right (col_offset=1).
    """
    coord = get_location(df, marker, match_index)
    if coord is None:
        return None

    row_idx, col_idx = coord

    # Calculate target coordinates
    target_row = row_idx + row_offset
    target_col = col_idx + col_offset

    # Safety check for bounds of sheet
    if target_row >= df.shape[0] or target_col >= df.shape[1]:
        print(f"   [ERROR] Target for '{marker}' is outside the sheet boundaries.")
        return None

    # .iloc uses [row, col]
    value = df.iloc[target_row, target_col]

    return value


def slice_table(
        df: pd.DataFrame,
        row_marker: str,
        match_index: int = 0
):
    """
    Extracts a table starting at the location of row_marker. Returns pandas dataframe with header.
    """
    # Get marker coordinates
    start_coords = get_location(df, row_marker, match_index)

    if start_coords is None:
        return None

    header_row_idx, col_start_idx = start_coords

    # Look for col cutoff (table width) in header row
    header_row_content = df.iloc[header_row_idx, col_start_idx:].astype(str)

    # Find the first index where the cell is 'nan' or empty; width of table is defined by header content
    col_width = next((i for i, val in enumerate(header_row_content)
                       if val.lower() == "nan" or not val.strip()), len(header_row_content))

    # Calculate the absolute end column
    col_end_idx = col_start_idx + col_width

    # Extract the table from large df
    df_extract = df.iloc[header_row_idx + 1:, col_start_idx:col_end_idx].reset_index(drop=True)
    df_extract.columns = df.iloc[header_row_idx, col_start_idx:col_end_idx].values  # Set header

    # Define row end of transfection scheme (first row with NA in first col)
    is_col1_na = df_extract.iloc[:,0].isna()

    if is_col1_na.any():
        # Find the positional index of the first NA value
        first_na_position = is_col1_na.values.argmax()
        # Slice the DataFrame using .iloc up to the row immediately before the NA row (exclusive)
        df_extract = df_extract.iloc[:first_na_position]
    else:
        # If no NA is found
        print("   [WARNING] Expected table row end delimiter (NaN in first column) not found; using full table.")
        return None

    return df_extract

def process_transfection_scheme(df: pd.DataFrame):
    """
    Gets the BRET pair (returns as list) and the experimental conditions (returns as dic) from the transfection table.
    """
    if df is None or df.empty:
        return [], {}

    # Separate metadata cols from transfection cols
    metadata_cols = {'DNA', 'DB#', 'Conc (ng/uL)', 'vol per transfection', 'vol master'}
    transfection_cols = [c for c in df.columns if c not in metadata_cols]

    if not transfection_cols:
        return [], {}

    # Remove row with empty DNA
    ignored_dna = {"pcdna3.1"}
    # Convert 'DNA' col to string -> strip whitespace -> lowercase -> check if matches ignored list
    # The tilde (~) means "NOT in"
    keep_mask = ~df['DNA'].astype(str).str.strip().str.lower().isin(ignored_dna)
    # Apply filter
    df = df[keep_mask].copy()

    transfection_dic = {}
    for col in transfection_cols:
        # Force column to numeric (coercing errors to NaN but fill them as 0)
        vals = pd.Series(pd.to_numeric(df[col], errors='coerce')).fillna(0)
        active_rows = df[vals > 0]

        # Get the DNA names for these rows
        dna_set = set(active_rows['DNA'].dropna().astype(str).tolist())
        transfection_dic[str(int(col))] = dna_set

    # Find DNA present in all transfection sets (BRET pair)
    if transfection_dic:
        common_dna = set.intersection(*transfection_dic.values())
    else:
        common_dna = set()
    main_plasmids = sorted(list(common_dna))

    # Find variable DNA (conditions)
    variable_dic = {}
    for col, dna_set in transfection_dic.items():
        # Subtract the common BRET pair from the specific set
        unique_dna = dna_set - common_dna
        variable_dic[col] = sorted(list(unique_dna))

    # print(f"identified conditions {variable_dic}")
    return main_plasmids, variable_dic

def extract_protocol_info(xls_obj: pd.ExcelFile):
    """
    Uses already opened pd.ExcelFiles (faster and more flexible than reading from path).
    Reads the 'Protocol' sheet of the protocol file and extracts all needed information.
    Stores and returns ProtocolData class with all info.
    """
    protocol_worksheet = "Protocol"

    try:
        protocol_sheet = pd.read_excel(xls_obj,
                                       sheet_name=protocol_worksheet,
                                       header=None)
    except ValueError:
        # Error if sheet is missing
        print(f"[ERROR] Worksheet '{protocol_worksheet}' not found in file.")
        return None

    file_name = os.path.basename(xls_obj.io)

    exp_date = None
    date_str = extract_value(protocol_sheet, "date of measurement")
    # Transform date str in real
    if date_str is not None:
        try:
            exp_date = datetime.strptime(date_str.strip(), '%d.%m.%y').date()
        except ValueError:
            print(f"   [WARNING] Protocol date '{date_str}' not in DD.MM.YY format.")

    exp_n = extract_value(protocol_sheet, "n =")

    selected_cell_lines = slice_table(protocol_sheet,"Cell line",1)
    if selected_cell_lines is not None:
        # Transform the selected cell lines into a list
        used_cell_lines = selected_cell_lines.iloc[:, 0].dropna().tolist()
    else:
        used_cell_lines = []

    cell_line_layout = extract_value(protocol_sheet,"Cell line layout", col_offset= 0, row_offset= 1)

    df_transfection = slice_table(protocol_sheet, "DNA")
    main_dna = []
    transfection_conditions = {}

    if df_transfection is not None:
        df_transfection = df_transfection.drop(columns=["vol per transfection", "vol master"], errors = 'ignore')

        main_dna, transfection_conditions = process_transfection_scheme(df_transfection)

    ligand_1 = extract_value(protocol_sheet, "Ligand dilution", col_offset=0, row_offset=1)
    ligand_1_conc = slice_table(protocol_sheet, "final concentration in well (log(M))")

    # Check whether second ligand was selected
    check_ligand_2 = extract_value(protocol_sheet, "Ligand dilution", col_offset=0, row_offset=1, match_index=1)
    if check_ligand_2 is not None and str(check_ligand_2).strip().lower() not in ["nan", ""]:
        ligand_2 = check_ligand_2
        ligand_conc_2 = slice_table(protocol_sheet, "final concentration in well (log(M))", match_index=1)
    else:
        ligand_2 = None
        ligand_conc_2 = None

    protocol_info = ProtocolData(file_name=file_name,
                                 exp_date = exp_date,
                                 n = exp_n,
                                 cell_lines = used_cell_lines,
                                 line_layout = cell_line_layout,
                                 transfection_scheme=df_transfection,
                                 main_plasmids= main_dna,
                                 transfection_conditions=transfection_conditions,
                                 ligand = ligand_1,
                                 ligand_conc=ligand_1_conc,
                                 ligand_2 = ligand_2,
                                 ligand_conc_2= ligand_conc_2)
    return protocol_info

def extract_metadata(xls_obj):
    """
    Uses already opened pd.ExcelFiles (faster and more flexible than reading from path).
    Reads the 'Table All Cycles' sheet from an analysis file and extracts metadata from the first column
    (measurement date, ID2: cell line, ID3: transfections #). Returns dictionary of extracted metadata.
    """
    metadata_worksheet = "Table All Cycles"

    try:
        # Read the first column of this sheet
        df_meta = pd.read_excel(xls_obj,
                                sheet_name=metadata_worksheet,
                                header=None,
                                usecols=[0],
                                nrows=30  # Limit rows to read
                                )
        col = df_meta[0].astype(str)  # Transform everything to str

        # Mapping: { "Excel Label": "Desired Key" }
        meta_keys = {"Date:": "measurement_date", "ID2:": "cell_line", "ID3:": "transfections"}
        metadata = {}
        for label, key in meta_keys.items():
            # n=1 split at first ":"; str[-1] select last arg; str-strip() remove spaces; .tolist() convert from pd series
            matches = col[col.str.contains(label, na=False)].str.split(":", n=1).str[-1].str.strip().tolist()
            if matches:
                metadata[key] = matches[0]

        # Handling of different spellings of cell lines
        if 'cell_line' in metadata:
            raw_cell_line = metadata['cell_line']
            # Split by comma to handle potentially multiple lines
            lines = [p.strip() for p in raw_cell_line.split(',')]
            standardised_lines = []

            for cl in lines:
                p_lower = cl.lower()
                if "dq" in p_lower:
                    standardised_lines.append("dQ")
                elif "ar" in p_lower:  # Covers "bArrKO"
                    standardised_lines.append("bArrKO")
                elif "con" in p_lower:  # Covers "Control", "Con", "con"
                    standardised_lines.append("Control")
                else:
                    standardised_lines.append(cl)  # Keep original if no rule matches

            # Re-join unique parts (e.g. "Control, dQ"); set() removed duplicates
            metadata['cell_line'] = ", ".join(list(set(standardised_lines)))

        # Transform date str to actual date
        if 'measurement_date' in metadata:
            date_str = metadata["measurement_date"]
            try:
                metadata["measurement_date"] = datetime.strptime(date_str.strip(), '%d/%m/%Y').date()
            except ValueError:
                print(f"   [WARNING] Analysis date '{date_str}' not in DD/MM/YYYY format.")
                return None  # Fail extraction if date is invalid

        # Ensure all required metadata fields were found
        if 'measurement_date' not in metadata or 'cell_line' not in metadata or 'transfections' not in metadata:
            print("   [WARNING] Missing Date, ID2, or ID3 from metadata sheet.")
            return None

        return metadata

    except ValueError:
        # Error if sheet is missing
        print(f"[ERROR] Worksheet '{metadata_worksheet}' not found in analysis file.")
        return None

def extract_bret_data(xls_obj):
    analysis_worksheet = "Analysis"

    try:
        bret_sheet = pd.read_excel(xls_obj,
                                sheet_name=analysis_worksheet,
                                header=None
                                )
        df_bret = slice_table(bret_sheet, "Time (min)")

        exclude_rows = ["Baseline", "late averg"]
        df_bret = df_bret[~df_bret["Time (min)"].astype(str).isin(exclude_rows)]
        df_bret = df_bret.reset_index(drop=True) # Make sure index is clean

        return df_bret

    except ValueError:
        # Error if sheet is missing
        print(f"[ERROR] Worksheet '{analysis_worksheet}' not found in analysis file.")
        return None

def extract_measurement_data(xls_obj):
    """
    Gets metadata and BRET ratio from analysis files.
    Stores and returns PrResult class with all data.
    """

    file_name = os.path.basename(xls_obj.io)
    metadata_dic = extract_metadata(xls_obj)
    bret_ratio_df = extract_bret_data(xls_obj)

    result_obj = PrResult(
        file_name=file_name,
        measurement_date=metadata_dic['measurement_date'],
        cell_line=metadata_dic['cell_line'],
        transfection_id=metadata_dic['transfections'],
        raw_bret_ratio_df=bret_ratio_df
    )
    return result_obj

def get_cell_line_map(layout_type: str, cell_lines: str):
    """
    Defines the plate layout for cell lines based on the dropdown selection protocol (.line_layout)
    and cell_lines in ID2 of the plate reader metadata (PrResult.cell_line)
    """
    print(f"\n[DEBUG] --- Mapping Cell Lines ---")
    print(f"[DEBUG] Layout Type: '{layout_type}' | Raw ID2: '{cell_lines}'")
    # Split ID2 string to get potentially multiple cell lines
    lines = [x.strip() for x in cell_lines.split(',')]
    print(f"[DEBUG] Parsed Cell Lines: {lines}")

    mapping = {}

    # Handle selection made in dropdown for line layout
    if "one line" in layout_type.lower():
        # Use the first (and likely only) cell line for all columns
        c_name = lines[0] if lines else "Unknown"
        for col in range(1, 13):
            mapping[col] = c_name

    elif "half" in layout_type.lower():
        line_1 = lines[0] if len(lines) > 0 else "Unknown_1"
        line_2 = lines[1] if len(lines) > 1 else "Unknown_2"

        for col in range(1, 7): mapping[col] = line_1
        for col in range(7, 13): mapping[col] = line_2

    elif "alternating" in layout_type.lower():
        line_1 = lines[0] if len(lines) > 0 else "Unknown_1"
        line_2 = lines[1] if len(lines) > 1 else "Unknown_2"

        # Block 1 (1-3) & Block 3 (7-9) -> Line 1
        for col in list(range(1, 4)) + list(range(7, 10)):
            mapping[col] = line_1

        # Block 2 (4-6) & Block 4 (10-12) -> Line 2
        for col in list(range(4, 7)) + list(range(10, 13)):
            mapping[col] = line_2

    else:
        print(f"   [WARNING] Unknown layout type: '{layout_type}'. Defaulting to global.")
        for col in range(1, 13): mapping[col] = cell_lines

    return mapping

def get_transfection_map(layout_type: str, t_ids: list[str], block_count: int = 4):
    """
    Defines the plate layout for blocks of transfection based on the cell line layout
    dropdown selection protocol (.line_layout) and # transfection in ID3 of
    the plate reader metadata (PrResult.transfection_id)
    """
    print(f"[DEBUG] --- Mapping Transfections ---")
    print(f"[DEBUG] Raw ID3 List: {t_ids}")
    # If no IDs, return empty
    if not t_ids:
        return ["N/A"] * block_count

    lower_layout = str(layout_type).lower()
    # IDs 1:1 to blocks
    if "one line" in lower_layout:
        # Extend list if shorter than blocks (fill with last or N/A)
        # Slicing [:block_count] ensures we don't overflow if ID3 has too many
        mapped_ids = (t_ids + ["N/A"] * block_count)[:block_count]
        return mapped_ids
    # IDs alternating in layout to cover all cell line x transfection combinations
    elif "half" in  lower_layout:
        if len(t_ids) >= 2:
            return [t_ids[0], t_ids[1], t_ids[0], t_ids[1]]
        elif len(t_ids) == 1:
            return [t_ids[0]] * 4
        else:
            return ["N/A"] * 4
    # IDs half/half of blocks to cover all cell line x transfection combinations
    elif "alternating" in lower_layout:
        if len(t_ids) >= 2:
            return [t_ids[0], t_ids[0], t_ids[1], t_ids[1]]
        elif len(t_ids) == 1:
            return [t_ids[0]] * 4
        else:
            return ["N/A"] * 4

    # Default fallback
    return (t_ids + ["N/A"] * block_count)[:block_count]

def get_block_start_for_col(col_index: int, plate_blocks: list[range]):
    """Finds the start column of the block that contains col_index."""
    for block in plate_blocks:
        if col_index in block:
            return block[0]  # Return the first column of that block (e.g., 1, 4, 7...)
    return None

def calculate_replicate_means(processed_df: pd.DataFrame,
                              plate_blocks: list[range],
                              col_metadata: dict): # Dic created from PlateColMetadata
    """
    Calculates the mean of technical replicates. As in processed_df each col is one well,
    the mean is performed of three cols within one block.
    Rows (A-H) are treated as distinct conditions (ligand concentration).
    Header format of returned df: "Condition_Name|Cell_Line|Row"
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
                header_key = f"{cond_name}|{cell_line}|{row}"
                mean_data[header_key] = mean_series

    return pd.DataFrame(mean_data)

def process_bret_measurement(result: PrResult, protocol: ProtocolData):
    """
    Maps cell line x transfection plate layout using protocol info.
    Performs Baseline Correction and Vehicle Normalization.
    Checks vehicle for outliers.
    Calculates AUC.(PENDING)
    """
    # TODO: handle ligand identitfy and concentrations
    if result.is_excluded:
        result.processed_df = None
        result.kinetic_mean_df = None
        result.auc_df = None
        return result

    # Reset for re-run
    result.warnings = []
    result.vehicle_outliers = {}

    print(f"\n[DEBUG] === Processing File: {result.file_name} ===")
    # --- CONFIG LAYOUT ---
    result.column_metadata = {}

    # Define triplicates (4 blocks); opt. edit for adding labeling layout
    plate_blocks = [
        range(1, 4),  # Block 1: Cols 1-3
        range(4, 7),  # Block 2: Cols 4-6
        range(7, 10),  # Block 3: Cols 7-9
        range(10, 13)  # Block 4: Cols 10-12
    ]

    # Get cell line map
    cl_map = get_cell_line_map(protocol.line_layout, result.cell_line)
    # Get transfection map via mapping ID3 info
    raw_ids = [x.strip() for x in str(result.transfection_id).split(',')] if result.transfection_id else []
    mapped_t_ids = get_transfection_map(protocol.line_layout, raw_ids, len(plate_blocks))
    print(f"[DEBUG] Mapped Block Sequence: {mapped_t_ids}")

    # Apply metadata on cols
    for i, block_cols in enumerate(plate_blocks):
        # Get the ID assigned to this block
        t_id = mapped_t_ids[i]

        # Resolve ID to Name (using Protocol)
        current_plasmids = []
        if t_id in protocol.transfection_conditions:
            current_plasmids = protocol.transfection_conditions[t_id]
            current_cond_name = " + ".join(sorted(current_plasmids))
        elif t_id == "N/A":
            current_cond_name = "Empty/NoID"
        else:
            current_cond_name = f"ID {t_id} (Missing)"

        # Assign to all columns in this block
        for col in block_cols:
            c_line = cl_map.get(col, "Unknown")

            meta = PlateColMetadata(
                cell_line=c_line,
                transfection_id=t_id,
                condition_name=current_cond_name,
                plasmids=current_plasmids
            )
            result.column_metadata[col] = meta

    # --- CONFIG PROCESSING ---
    # Define accepted vehicle range
    acc_vehicle_range = 0.2
    # Define Baseline: First 5 rows (Index 0-4)
    baseline_end_idx = 5

    # Get raw BRET ratio table
    raw_df = result.raw_bret_ratio_df.copy()
    time_col = "Time (min)"

    # Define Kinetic: Remaining rows (Index 5 onwards)
    if len(raw_df) < baseline_end_idx:
        print(f"   [WARNING] Data has less than {baseline_end_idx} rows.")
        return None, None

    # Prepare the full dataframe for normalization (removing the time col)
    data_df = raw_df.drop(columns=[time_col]).copy()

    # --- BASELINE CORRECTION ---
    # Isolate baseline rows for calculation
    df_baseline_calc = raw_df.iloc[0:baseline_end_idx].copy()
    bl_corrected_df = pd.DataFrame()

    for col in data_df.columns:
        if col in result.excluded_wells:
            bl_corrected_df[col] = None
            continue
        # Mean of baseline rows for this well
        # Force numeric conversion for baseline values to handle potential strings/decimals
        base_mean = pd.to_numeric(df_baseline_calc[col], errors='coerce').mean()

        if base_mean != 0:
            col_vals = pd.to_numeric(data_df[col], errors='coerce')
            # Divide ALL values by the baseline mean
            bl_corrected_df[col] = col_vals / base_mean
        else:
            bl_corrected_df[col] = None

    # --- VEHICLE CORRECTION ---
    # Calculate mean vehicle for each condition
    vehicle_mean = {}

    for block in plate_blocks:
        start_col = block[0]  # e.g., 1, 4, 7, 10 for triplicates

        wells = [f"H{c}" for c in block]
        # Filter for wells that actually exist in the dataframe
        # Logic needed for optionally excluding wells
        valid_vehicles = [w for w in wells if w in bl_corrected_df.columns and w not in result.excluded_wells]

        if valid_vehicles:
            # Get vehicle values and make sure that data is numeric
            vehicle_data = bl_corrected_df[valid_vehicles].apply(pd.to_numeric, errors='coerce')

            # Check bounds of vehicle
            for well in valid_vehicles:
                # Check kinetic mean
                veh_kinetic_mean = vehicle_data[well].mean(axis=0)

                if abs(veh_kinetic_mean - 1) > acc_vehicle_range:
                    # Store well name and value in dic
                    result.vehicle_outliers[well] = float(veh_kinetic_mean)

            # Calculate mean across valid vehicle wells per time point
            vehicle_mean[start_col] = vehicle_data.mean(axis=1)
        else:
            vehicle_mean[start_col] = None

    # Normalise to mean(vehicle)
    vehicle_corr_df = bl_corrected_df.copy()

    for col in data_df.columns:
        # Get well number
        col_num = int(col[1:])

        # Find block start dynamically using the list
        block_start = get_block_start_for_col(col_num, plate_blocks)

        if block_start in vehicle_mean and vehicle_mean[block_start] is not None:
            vehicle_corr_df[col] = bl_corrected_df[col] / vehicle_mean[block_start]
        else:
            vehicle_corr_df[col] = None

    # Store processed df before processing further
    result.processed_df = vehicle_corr_df

    # --- MEAN OF REPLICATES (KINETIC) ---
    result.kinetic_mean_df = calculate_replicate_means(
        processed_df=result.processed_df,
        plate_blocks=plate_blocks,
        col_metadata=result.column_metadata
    )
    print("-"*40)
    print(result.kinetic_mean_df)
    print("-" * 40)

    # TODO: add AUC calculation + subsequent mean of replicates
    return result


# TODO: add function to rearrange cols of processed bret df (flexible for user interaction)

# --- Main Application --- #

class NCollectorApp:
    def __init__(self, main_window):
        self.main_gi = main_window
        main_window.title("N Collector")

        # Separate log window
        self.log_window = tk.Toplevel(main_window)
        self.log_window.title("Processing Log")
        self.log_window.geometry("700x500")

        # Log window to display print statements
        self.log_text = tk.Text(self.log_window)
        self.log_text.pack(expand=True, fill='both')

        # --- GUI State ---
        # Path to folder variable
        self.folder_path = tk.StringVar(value="No folder selected.")
        self.subfolder_paths_with_files = []
        self.experiment: list[MeasurementFolder] = []
        self.master_df = pd.DataFrame()
        self.rule_history_text = ""

        # --- TABS SETUP ---
        self.notebook = ttk.Notebook(main_window)
        self.notebook.pack(expand=True, fill='both')

        # Tab 1: Import Data
        self.tab_import = tk.Frame(self.notebook)
        self.notebook.add(self.tab_import, text="Import & Export Data")

        # Tab 2: Data Selection
        self.tab_select = tk.Frame(self.notebook)
        self.notebook.add(self.tab_select, text="Optional Selection")

        # --- TAB 1 CONTENT ---
        # Select Folder button
        # No self. needed as this does not have to be stored for later changes
        tk.Button(self.tab_import, text="Select folder containing results of experiment",
                  command=self.select_folder).pack(pady=10, padx=10)

        # Display label for path
        self.path_label = tk.Label(self.tab_import,
                                   textvariable=self.folder_path,
                                   wraplength=1000,
                                   justify="left",
                                   font=('Arial', 10))
        self.path_label.pack(pady=10, padx=10)

        # Load button
        self.load_files_button = tk.Button(self.tab_import,
                                           text="Load Files",
                                           state="disabled",
                                           command=self.collect_files
                                           )
        self.load_files_button.pack(pady=15)

        # Display of N summary table
        summary_frame = tk.Frame(self.tab_import)
        summary_frame.pack(pady=10, fill="both", expand=True ,padx=20)
        # Scrollbar for table
        tree_scroll = tk.Scrollbar(summary_frame)
        tree_scroll.pack(side="right", fill="y")

        self.summary_tree = ttk.Treeview(summary_frame,
                                         columns=("Cell", "Cond", "N", "Dates"),
                                         show="headings",
                                         yscrollcommand=tree_scroll.set,
                                         height=6)
        tree_scroll.config(command=self.summary_tree.yview)

        # Define Columns
        self.summary_tree.heading("Cell", text="Cell Line")
        self.summary_tree.heading("Cond", text="Condition")
        self.summary_tree.heading("N", text="N")
        self.summary_tree.heading("Dates", text="Dates")

        self.summary_tree.column("Cell", width=100)
        self.summary_tree.column("Cond", width=250)
        self.summary_tree.column("N", width=30, anchor="center")
        self.summary_tree.column("Dates", width=150)

        self.summary_tree.pack(fill="both", expand=True)

        # Text of applies exclusion rules
        rules_frame = tk.LabelFrame(self.tab_import, text="Applied Exclusion Rules")
        rules_frame.pack(fill="x", padx=20, pady=5)

        self.lbl_rules_summary = tk.Label(rules_frame,
                                          text="No exclusion rules applied",
                                          justify="left",
                                          anchor="w"
                                          )
        self.lbl_rules_summary.pack(fill="x", padx=5, pady=5)

        # Export button
        self.export_button = tk.Button(self.tab_import,
                                       text="Export data",
                                       state="disabled",
                                       command=self.export_data)
        self.export_button.pack(pady=20, padx=10)

        # --- TAB 2 CONTENT ---
        self.setup_exclusion_tab()

    # Helper function to write to that text box
    def log(self, message):
        """Logs to the separate window"""
        try:
            self.log_text.insert(tk.END, message + "\n")
            self.log_text.see(tk.END)
        except tk.TclError:
            # Handle case where user manually closed log window but app is running
            print(message)

    def update_summary_table(self):
        """Fills the summary table with N counts and dates per condition"""
        # Clear existing data
        for i in self.summary_tree.get_children():
            self.summary_tree.delete(i)

        if self.master_df.empty:
            return

        # Group by Cell Line and Condition
        grouped = self.master_df.groupby(['Cell_Line', 'Condition'])

        for (cell, cond), group in grouped:
            # Count unique filenames for N (experiments)
            n_count = group['File_Name'].nunique()

            # Get sorted unique dates
            unique_dates = sorted(group['Date'].unique())
            date_str = ", ".join(unique_dates)

            # Insert into tree
            self.summary_tree.insert("", "end", values=(cell, cond, n_count, date_str))

    def setup_exclusion_tab(self):
        """Builds GUI for Tab 2 data selection"""

        # Frame for dropdowns
        filter_frame = tk.LabelFrame(self.tab_select, text="Exclude Data")
        filter_frame.pack(fill = "x", pady=10, padx=5)

        # Variables
        self.var_date = tk.StringVar(value="All")
        self.var_cell = tk.StringVar(value="All")
        self.var_cond = tk.StringVar(value="All")
        self.var_repl = tk.StringVar(value="")
        self.var_row = tk.StringVar(value="")
        # List to store rules
        self.pending_exclusions = []

        # 1. Date Dropdown
        tk.Label(filter_frame, text="Date:").grid(row=0, column=0, padx=5, pady=5)
        self.cb_date = ttk.Combobox(filter_frame, textvariable=self.var_date, state="readonly")
        self.cb_date.grid(row=0, column=1, padx=5, pady=5)
        self.cb_date.bind("<<ComboboxSelected>>", self.update_cell_options)

        # 2. Cell Line Dropdown
        tk.Label(filter_frame, text="Cell Line:").grid(row=0, column=2, padx=5, pady=5)
        self.cb_cell = ttk.Combobox(filter_frame, textvariable=self.var_cell, state="readonly")
        self.cb_cell.grid(row=0, column=3, padx=5, pady=5)
        self.cb_cell.bind("<<ComboboxSelected>>", self.update_cond_options)

        # 3. Condition Dropdown
        tk.Label(filter_frame, text="Condition:").grid(row=0, column=4, padx=5, pady=5)
        self.cb_cond = ttk.Combobox(filter_frame, textvariable=self.var_cond, state="readonly")
        self.cb_cond.grid(row=0, column=5, padx=5, pady=5)
        self.cb_cond.bind("<<ComboboxSelected>>", self.update_repl_options)

        # 4. Granular Filters Col (replicate) and Row (ligand conc)
        granular_frame = tk.Frame(filter_frame)
        granular_frame.grid(row=1, column=0, columnspan=6, pady=5, sticky="w")

        tk.Label(granular_frame, text="Technical Replicate:").pack(side="left", padx=5)
        self.cb_rep = ttk.Combobox(granular_frame, textvariable=self.var_repl, state="readonly", width=5)
        self.cb_rep.pack(side="left", padx=5)
        self.cb_rep['values'] = ["", "1", "2", "3"]
        self.cb_rep.bind("<<ComboboxSelected>>", self.toggle_row_dropdown)

        # Row dropdown (Initially disabled/hidden until Replicate is picked)
        self.lbl_row = tk.Label(granular_frame, text="Specific Row:")
        self.lbl_row.pack(side="left", padx=5)

        self.cb_row = ttk.Combobox(granular_frame, textvariable=self.var_row, state="disabled", width=5)
        self.cb_row.pack(side="left", padx=5)
        self.cb_row['values'] = ["", "A", "B", "C", "D", "E", "F", "G", "H"]

        # Buttons
        btn_frame = tk.Frame(filter_frame)
        btn_frame.grid(row=2, column=0, columnspan=6, pady=10)

        tk.Button(btn_frame, text="Add Rule to List", command=self.add_exclusion_rule).pack(side="left", padx=5)
        tk.Button(btn_frame, text="Clear List", command=self.clear_exclusion_list).pack(side="left", padx=5)

        # Listbox for Pending Exclusions
        list_frame = tk.LabelFrame(self.tab_select, text="Pending Exclusions (Will be removed upon Apply)")
        list_frame.pack(fill="both", expand=True, padx=10, pady=5)

        self.lb_exclusions = tk.Listbox(list_frame, height=8)
        self.lb_exclusions.pack(fill="both", expand=True, padx=5, pady=5)

        # Apply Button (Bottom)
        tk.Button(self.tab_select, text="Apply exclusions and re-calculate",
                  command=self.apply_exclusions).pack(pady=10, ipadx=10)

    def update_repl_options(self, event=None):
        """Reset replicate when conditions changes"""
        self.var_repl.set("")
        self.toggle_row_dropdown()

    def toggle_row_dropdown(self, event=None):
        """Enable row dropdown only if a specific replicate is selected"""
        if self.var_repl.get() != "":
            self.cb_row.config(state="readonly")
        else:
            self.var_repl.set("")
            self.cb_row.config(state="disabled")

    def refresh_filter_options(self):
        """Called during built master index. Updates dropdown options of date, cell line and condition."""
        if self.master_df.empty: return

        # Get unique dates and add "All"
        dates = sorted(self.master_df['Date'].unique().tolist())
        self.cb_date['values'] = ["All"] + dates
        self.var_date.set("All")

        # Reset others
        self.update_cell_options()

    def update_cell_options(self, event=None):
        """Updates Cell Line options based on selected Date."""
        selected_date = self.var_date.get()

        if self.master_df.empty: return

        if selected_date == "All":
            # Show all cell lines available in the whole dataset
            cells = sorted(self.master_df['Cell_Line'].unique().tolist())
        else:
            # Filter DF by date
            subset = self.master_df[self.master_df['Date'] == selected_date]
            cells = sorted(subset['Cell_Line'].unique().tolist())

        self.cb_cell['values'] = ["All"] + cells
        self.var_cell.set("All")
        self.update_cond_options()

    def update_cond_options(self, event=None):
        """Updates Condition options based on selected Date AND Cell Line."""
        selected_date = self.var_date.get()
        selected_cell = self.var_cell.get()

        if self.master_df.empty: return

        # Start with full DF
        subset = self.master_df.copy()

        # Apply Date Filter
        if selected_date != "All":
            subset = subset[subset['Date'] == selected_date]

        # Apply Cell Filter
        if selected_cell != "All":
            subset = subset[subset['Cell_Line'] == selected_cell]

        conds = sorted(subset['Condition'].unique().tolist())
        self.cb_cond['values'] = ["All"] + conds
        self.var_cond.set("All")

    def add_exclusion_rule(self):
        """
        Adds the current dropdown state to the pending list.
        Handles empty strings for Replicate/Col and Row.
        """
        rep_val = self.var_repl.get()
        row_val = self.var_row.get()

        rule = {
            "Date": self.var_date.get(),
            "Cell_Line": self.var_cell.get(),
            "Condition": self.var_cond.get(),
            "Replicate": rep_val,
            "Row": row_val
        }

        # Create display string
        rep_str = rep_val if rep_val else "All (1-3)"
        row_str = row_val if row_val else "All (A-H)"

        # Check for duplicates or empty
        rule_str = f"Date: {rule['Date']} | Cell: {rule['Cell_Line']} | Cond: {rule['Condition']} | Rep:{rep_str} | Row:{row_str}"

        self.pending_exclusions.append(rule)
        self.lb_exclusions.insert(tk.END, rule_str)

    def clear_exclusion_list(self):
        self.pending_exclusions = []
        self.lb_exclusions.delete(0, tk.END)

    def apply_exclusions(self):
        """
        Iterates through pending rules, finds matching rows in Master DF,
        and updates the 'excluded_wells' list in the Ref_Result objects.
        """
        if not self.pending_exclusions:
            self.log("No exclusion rules defined.")
            return

        self.log(f"\n--- Applying {len(self.pending_exclusions)} Exclusion Rules ---")
        # --- Save applied rules as text for displaying
        pending_lines = self.lb_exclusions.get(0, tk.END)
        new_text_block = "\n".join(pending_lines)

        if self.rule_history_text:
            self.rule_history_text += "\n" + new_text_block
        else:
            self.rule_history_text = new_text_block
        self.lbl_rules_summary.config(text=self.rule_history_text) # Update GUI

        # --- Applying the rules
        count_wells = 0
        count_files = 0

        for rule in self.pending_exclusions:
            # If entire date is excluded
            if (rule['Date'] != "All" and
                    rule['Cell_Line'] == "All" and
                    rule['Condition'] == "All" and
                    rule['Replicate'] == "" and
                    rule['Row'] == ""):

                # Find matching files and exclude them entirely
                for folder in self.experiment:
                    # Date formatting match
                    if folder.measurement_date.strftime('%d.%m.%y') == rule['Date']:
                        for res in folder.results:
                            if not res.is_excluded:
                                res.is_excluded = True
                                count_files += 1
                continue

            # Start with full dataframe
            df = self.master_df.copy()

            # Apply high level filters
            if rule['Date'] != "All":
                df = df[df['Date'] == rule['Date']]
            if rule['Cell_Line'] != "All":
                df = df[df['Cell_Line'] == rule['Cell_Line']]
            if rule['Condition'] != "All":
                df = df[df['Condition'] == rule['Condition']]

            if df.empty:
                self.log(f"   [WARNING] Rule {rule} matched 0 records.")
                continue

            # Logic to handle replicates and wells
            # Grouping to handle different plates separately
            grouped_by_file = df.groupby('File_Name')

            for file_name, group_df in grouped_by_file:
                # group_df should contain one row per plate column of specified block

                if rule['Replicate'] == "":
                    # If no replicate specified, take all columns of this block
                    rows_to_process = group_df
                else:
                    try:
                        repl_index = int(rule['Replicate']) - 1  # Convert "1" -> 0
                        rows_to_process = group_df.iloc[[repl_index]]
                    except IndexError:
                        self.log(f"   [SKIP] File {file_name} does not have replicate {rule['Replicate']}")
                        continue

                # Iterate through the specific rows to store the well ids in PrResult (result_obj)
                for index, row_data in rows_to_process.iterrows():
                    # Assigning to PrResult to store well ids in (PrResult.excluded_wells)
                    # -> by assigning it to result_obj defined previously it globally changes this obj (python logic!)
                    result_obj = row_data['Ref_Result']
                    col_idx = int(row_data['Column_Index'])

                    # Determine Rows (A-H)
                    target_rows = "ABCDEFGH"
                    # If row in plate is specified (not empty), use the specified row to built ID
                    if rule['Row'] != "": target_rows = rule['Row']

                    # Generate Well IDs and Append to Object
                    for r in target_rows:
                        well_id = f"{r}{col_idx}"

                        # Modify the object directly (Objects are mutable, so this updates the global state)
                        if well_id not in result_obj.excluded_wells:
                            result_obj.excluded_wells.append(well_id)
                            count_wells += 1

                            print(f"Excluded {well_id} in {result_obj.file_name}")

        if count_files > 0:
            self.log(f"   [DONE] Excluded {count_files} entire files.")
        if count_wells > 0:
            self.log(f"   [DONE] Excluded {count_wells} specific wells.")

        # Clear list after applying
        self.clear_exclusion_list()

        # Rerun processing to update graphs/stats
        self.run_processing_pipeline()

    def select_folder(self):
            """Opens dialog to select folder to search for xlsx files in"""
            directory = filedialog.askdirectory(title="Select a folder...")
            if not directory: return

            # --- RESET STATE: Clear old data when a new folder is selected
            self.subfolder_paths_with_files = []
            self.experiment = []
            self.master_df = pd.DataFrame()
            # Reset summary table
            for i in self.summary_tree.get_children():
                self.summary_tree.delete(i)
            self.rule_history_text = ""
            self.lbl_rules_summary.config(text="")
            self.export_button.config(state="disabled")


            # Update GUI immediately
            self.folder_path.set(f"Selected Path: {directory}\n\nScanning for files...")
            self.load_files_button.config(state="disabled")
            self.main_gi.update()

            # --- Scan new directory for xlsx or xlsm---
            for root, dirs, files in os.walk(directory):
                if any(f.endswith((".xlsx", ".xlsm")) for f in files):
                    self.subfolder_paths_with_files.append(root) # Save paths of files

            # Update GUI
            if self.subfolder_paths_with_files:
                count = len(self.subfolder_paths_with_files)
                folder_names = [os.path.basename(path) for path in self.subfolder_paths_with_files]
                folder_names_string = "\n ".join(folder_names)
                self.folder_path.set(
                    f"Selected Path: {directory}\n\n Found following subfolders with xlsx/xlsm files:\n {folder_names_string}")
                self.load_files_button.config(state="normal")
                print(f"Found {count} folders: \n {folder_names_string}")
            else:
                self.folder_path.set(f"Error: No .xlsx or .xlsm files found in {directory} or any subfolder.")
                self.load_files_button.config(state="disabled")
                print(f"No .xlsx or .xlsm files found starting from: {directory}")

    def handle_main_plasmids_selection(self):
        """
        Checks main_plasmids consistency. If multiple sets found, user selects one.
        Filters self.experiment to keep only the selected group.
        Returns the string representation of the selected set.
        """
        # Group experiments by their main_plasmids (convert list to tuple for dictionary key)
        main_plasmids_groups = {}
        for folder in self.experiment:
            if folder.protocol and folder.protocol.main_plasmids:
                key = tuple(folder.protocol.main_plasmids)
                if key not in main_plasmids_groups:
                    main_plasmids_groups[key] = []
                main_plasmids_groups[key].append(folder)

        if not main_plasmids_groups:
            return "No Common Plasmids Detected"

        # If only one set exists, return it immediately
        if len(main_plasmids_groups) == 1:
            return " + ".join(list(main_plasmids_groups.keys())[0])

        # --- Multiple Sets Detected: Ask User ---
        # Create a modal dialog window
        dialog = tk.Toplevel(self.main_gi)
        dialog.title("Select Experiment")

        tk.Label(dialog, text="Different experiment set ups detected across folders.\nSelect one to process:",
                 font=("Arial", 11, "bold")).pack(pady=10)

        selected_var = tk.StringVar()
        first_key = list(main_plasmids_groups.keys())[0]
        selected_var.set(str(first_key))  # Set default

        # Helper to map string back to tuple key
        str_to_key_map = {}

        for key in main_plasmids_groups:
            pair_str = " + ".join(key)
            key_val_str = str(key)
            str_to_key_map[key_val_str] = key

            text_label = f"{pair_str} ({len(main_plasmids_groups[key])} folders)"
            tk.Radiobutton(dialog, text=text_label, variable=selected_var, value=key_val_str).pack(anchor="w", padx=20)

        def on_confirm():
            dialog.destroy()
        tk.Button(dialog, text="Confirm", command=on_confirm).pack(pady=20)
        self.main_gi.wait_window(dialog) # Wait until the window is closed

        # Retrieve selection
        selected_key_str = selected_var.get()
        selected_key = str_to_key_map.get(selected_key_str, first_key)

        # Filter the experiment list
        self.experiment = main_plasmids_groups[selected_key]
        print(f"   [FILTER] Keeping {len(self.experiment)} folders matching main plasmids: {selected_key}")

        return " + ".join(selected_key)

    def built_master_index(self):
        """
        Creates a pd dataframe containing all metadata for all wells.
        Enables flexible filtering needed for exclusion of data.
        """

        print("----- building master index ------")

        # collect list of records for each col
        records = []
        rows_str = ["ABCDEFGH"]

        for folder in self.experiment:
            for result in folder.results:
                # Get col metadata first
                if not result.column_metadata: continue
                # Skip if file is excluded
                if result.is_excluded: continue

                # Iterate through the mapped cols (1-12)
                for col_idx, meta in result.column_metadata.items():
                    # Filter out empty cols
                    if not meta.condition_name or "Empty" in meta.condition_name: continue
                    # Check if this column is completely excluded
                    all_wells_excluded = True
                    for r in rows_str:
                        well_id = f"{r}{col_idx}"
                        if well_id not in result.excluded_wells:
                            all_wells_excluded = False
                            break
                    # If all wells are excluded, do not add to summary table (N count decreases)
                    if all_wells_excluded: continue

                    # Create a record for this col
                    record = {
                        "File_Name": result.file_name,
                        "Date": result.measurement_date.strftime('%d.%m.%y'),  # String for dropdowns
                        "Cell_Line": meta.cell_line,
                        "Condition": meta.condition_name,
                        "Transfection_ID": meta.transfection_id,
                        "Column_Index": col_idx,
                        "Ref_Result": result  # Store the actual object to manipulate later
                    }
                    records.append(record)
        # Built df from records
        if records:
            self.master_df = pd.DataFrame(records)

            # --- Summary for verification ---
            summary = self.master_df.groupby(['Cell_Line', 'Condition'])['File_Name'].nunique()
            print("\n[DEBUG] Data Summary:\n", summary)
            self.refresh_filter_options()
            self.update_summary_table()
            return self.master_df
        else:
            self.master_df = pd.DataFrame()
            self.update_summary_table()
            print("No valid data found")
            return self.master_df

    def collect_files(self):
        """
        1. LOAD FILES
        Reads sheet names of all xlsx and xlsm files to identify and separate protocol and result analysis files.
        Validation of correct protocol to analysis files is done via date of measurement in the folder name.
        """
        if not self.subfolder_paths_with_files:
            print("No folders to analyze.")
            return

        self.log("\n--- Starting Data Collection ---")

        # Collect all subfolder info (experiment repeats) in this list
        self.experiment = []

        # Iterate over subfolders
        for folder_path in self.subfolder_paths_with_files:
            folder_name = os.path.basename(folder_path)
            self.log(f"\n--- Processing Folder: {folder_name} ---")

            # Use date in folder name for validation
            try:
                # Transform date str to date obj (YYMMDD)
                folder_date_obj = datetime.strptime(folder_name.split("_")[0], '%y%m%d').date()
            except ValueError:
                self.log(f"   [ERROR] Folder '{folder_name}' invalid date format. Expected YYMMDD. Skipping.")
                continue

            # Create MeasurementFolder object to collect protocol and results
            folder_data = MeasurementFolder(folder_name=folder_name,
                                            folder_path=folder_path,
                                            measurement_date=folder_date_obj)

            files_in_folder = [f for f in os.listdir(folder_path) if f.endswith(('.xlsx', '.xlsm'))]
            for file_name in files_in_folder:
                file_path = os.path.join(folder_path, file_name)
                # Logical variable to keep track of skipped files
                is_imported = False

                try:
                    # Use ExcelFile to check sheet names
                    xls = pd.ExcelFile(file_path) # fast reading of multiple sheet xls files
                    sheet_names = xls.sheet_names

                    # Identify Protocol File
                    if "Protocol" in sheet_names:
                        protocol_info = extract_protocol_info(xls)

                        if protocol_info and protocol_info.exp_date == folder_date_obj:
                            folder_data.protocol = protocol_info
                            self.log(f"   [PROTOCOL] loaded: {file_name}")
                            is_imported = True
                        elif protocol_info:
                            self.log(f"   [MISMATCH] Protocol {protocol_info.exp_date} != Folder {folder_date_obj}")

                    # Identify Analysis File
                    elif "Analysis" in sheet_names:
                        meas_data = extract_measurement_data(xls)

                        if meas_data and meas_data.measurement_date == folder_date_obj:
                            folder_data.results.append(meas_data)
                            self.log(f"   [MEASUREMENT] loaded {file_name}")
                            is_imported = True
                        else:
                            self.log(f"   [MISMATCH] Analysis {meas_data.measurement_date} != Folder {folder_date_obj}")

                    # Files not matching criteria are skipped
                    if not is_imported:
                        folder_data.skipped_files.append(file_name)
                except Exception as e:
                    self.log(f"   [ERROR] Could not read {file_name}: {e}")

            # Store the collected data for this experiment
            self.experiment.append(folder_data)
            print(f"   [SKIPPED]: {folder_data.skipped_files}")

        self.log(f"--- Loading Complete. Loaded {len(self.experiment)} folders. ---")

        # Verification of readings
        print("\n" + "=" * 30)
        print("COLLECTION SUMMARY")
        print("=" * 30)

        for rep in self.experiment:
            print(f"\nFolder: {rep.folder_name}")
            if rep.protocol:
                print(f"  [✓] Protocol: {rep.protocol.file_name}")
            else:
                print(f"  [ ] Protocol: MISSING")

            res_count = len(rep.results)
            print(f"  [i] Results: {res_count} file(s) loaded")
            if rep.skipped_files:
                print(f"   [SKIPPED]: {rep.skipped_files}")

        print("\n" + "=" * 30)

        # Call processing
        self.run_processing_pipeline()

        # Enable export
        self.export_button.config(state="normal")

    def run_processing_pipeline(self):
        """
        2. Processing of raw BRET data and indexing with protocol info
        Calls processing functions and is rerun if data was excluded.
        """
        self.log("\n--- Starting Processing Pipeline ---")

        # Check for plasmids transfected in all conditions (main plasmids) and filter if needed
        selected_exp_name = self.handle_main_plasmids_selection()

        # Iterate through data and perform mapping+calculations
        for folder in self.experiment:
            if not folder.protocol: continue # Protocol is needed for processing
            for result in folder.results:
                # Process each result file within one folder (belonging to one protocol)
                # Also assigns conditions to data
                result = process_bret_measurement(result, folder.protocol)

            # Handle outliers stored in dic
            if result.vehicle_outliers:
                vehicle_out = ", ".join([f"{well} = {val:.2f}" for well, val in result.vehicle_outliers.items()])
                self.log(f"   [VEHICLE WARNING] {result.file_name}: {vehicle_out}")

        # Built master indexing table (needed for flexible data exclusion)
        self.built_master_index()

        self.log("\n--- Processing Complete ---")

    def export_data(self):
        """
        Exports kinetic mean data as multi-sheet xlsx (for now) in user-selected path.
        - one sheet per condition + row, with all cell lines
        - time col is created based on baseline offset
        """
        # TODO: improve time col definition -> handle in process_bret_measurement
        # TODO: improve layout of exporting all data
        # TODO: add plotting helper tab (loading all data exports)

        if not self.experiment:
            self.log("No data to export.")
            return

        export_file_path = filedialog.asksaveasfilename(
            defaultextension=".xlsx",
            filetypes=[("Excel files", "*.xlsx")],
            title="Export Kinetic Means"
        )
        if not export_file_path: return

        self.log(f"--- Exporting to {os.path.basename(export_file_path)} ---")

        # Structure: export_tree[SheetName][CellLine] = [Series_N1, Series_N2...]
        export_tree = {}

        # Setting for baseline measurement count
        baseline_count = 5

        # Collect data from all experiment folders
        for folder in self.experiment:
            for result in folder.results:
                # Skip is whole file is excluded
                if result.is_excluded or result.kinetic_mean_df is None:
                    continue

                # Generate time index
                n_points = len(result.kinetic_mean_df)
                time_index = range(-baseline_count, n_points-baseline_count)

                # Iterate through cols created in kinetic_mean_df
                for col_key in result.kinetic_mean_df.columns:
                    try:
                        cond, cell, row = col_key.split('|')
                    except ValueError:
                        continue # Skip malformed cols

                    # Create sheet name: Condition + Row
                    # (Limit to 31 chars for Excel compatibility!!!)
                    sheet_name = f"{cond}_{row}"
                    bad_chars = "[]:*?/\\"
                    for c in bad_chars: sheet_name = sheet_name.replace(c, "_")
                    sheet_name = sheet_name[:31]

                    # Init storage
                    if sheet_name not in export_tree:
                        export_tree[sheet_name] = {}
                    if cell not in export_tree[sheet_name]:
                        export_tree[sheet_name][cell] = []

                    # Get data and set time index
                    series = result.kinetic_mean_df[col_key].copy()
                    series.index = time_index

                    export_tree[sheet_name][cell].append(series)

        # Write to excel
        try:
            with pd.ExcelWriter(export_file_path) as writer:
                # Iterate through sheets (conditions + row)
                for sheet_name, cell_data in export_tree.items():
                    # Sort cell lines alphabetically
                    sorted_cells = sorted(cell_data.keys())

                    # Create a list of DataFrames to concatenate
                    dfs_to_concat = []
                    header_list = ["Time (min)"]

                    # Built dataframe for this sheet
                    for cell in sorted_cells:
                        replicates = cell_data[cell]

                        # add to header list
                        header_list.append(cell)
                        header_list.extend([""] * (len(replicates) - 1))

                        # 2. Collect Data
                        df_cell = pd.concat(replicates, axis=1)
                        dfs_to_concat.append(df_cell)

                    if not dfs_to_concat: continue

                    # Combine all data side-by-side
                    full_sheet_df = pd.concat(dfs_to_concat, axis=1)

                    # Sort by Index (Time)
                    full_sheet_df.sort_index(inplace=True)

                    # Reset index so "Time" becomes the first column (col 0)
                    full_sheet_df.reset_index(inplace=True)

                    # Force the columns to match our manual header list
                    # (Pandas allows duplicate column names like "" here)
                    full_sheet_df.columns = header_list

                    # Write to Excel without the default index or header processing
                    full_sheet_df.to_excel(writer, sheet_name=sheet_name, index=False)

            self.log("   [SUCCESS] Export complete.")

        except Exception as e:
            self.log(f"   [ERROR] Export failed: {e}")
            print(e)

# TODO: implement showing also errors from tool functions in log window

# --- Main Execution Block ---
if __name__ == "__main__":
    # Create the main window
    root = tk.Tk()

    # Create an instance of the application
    app = NCollectorApp(root)

    # Start the Tkinter event loop
    root.mainloop()