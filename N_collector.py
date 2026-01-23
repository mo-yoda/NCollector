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
    ligand_identity: str = "N/A"
    # Dic mapping row A-H to concentration (float)
    ligand_conc: dict[str, float] = field(default_factory=dict)

@dataclass
class PrResult:
    """ Raw and processed information from a single _analysis file """
    # --- Information from analysis xlsx itself ---
    file_name: str
    measurement_date: date
    cell_line: str # ID2
    transfection_id: str # ID3
    raw_bret_ratio_df: pd.DataFrame

    # --- Connection to protocol file ---
    # Key = Column Index (1-12), Value = WellMetadata object
    column_metadata: dict[int, PlateColMetadata] = field(default_factory=dict)

    # --- Processed BRET data ---
    # Time column
    time_vector: list[float] = field(default_factory=list)
    # Baseline corrected kinetic data
    bl_corr_kinetic: pd.DataFrame | None = None
    # Baseline- and vehicle-normalised kinetic data (technical replicates)
    kinetic_df: pd.DataFrame | None = None
    # Mean of baseline- and vehicle-normalised kinetic data
    kinetic_mean_df: pd.DataFrame | None = None
    # Pre-vehicle norm AUC
    raw_auc_df: pd.DataFrame | None = None
    # Baseline- and vehicle-normalised AUC data (technical replicates)
    auc_df: pd.DataFrame | None = None
    # Mean of baseline- and vehicle-normalised AUC data
    auc_mean_df:pd.DataFrame | None = None

    # --- Export tidy CRC data ---
    raw_auc_tidy_df: pd.DataFrame | None = None
    auc_tidy_df: pd.DataFrame | None = None
    auc_mean_tidy_df: pd.DataFrame | None = None

    # --- Optional exclusion by user interaction  ---
    is_excluded: bool = False
    excluded_wells: list[str] = field(default_factory=list)

    # --- Internal check and warnings for helping vehicle outlier identification ---
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
    ligand_2_conc: pd.DataFrame | None = None

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
        ligand_2_conc = slice_table(protocol_sheet, "final concentration in well (log(M))", match_index=1)
    else:
        ligand_2 = None
        ligand_2_conc = None

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
                                 ligand_2_conc= ligand_2_conc)
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

def built_conc_dic(df_conc: pd.DataFrame):
    """
    Parses the ligand concentration DataFrame (from ProtocolData) into a dict.
    Assumes standard 8-row layout corresponding to A-H
    """

    if df_conc is None or df_conc.empty:
        return {}

    conc_dic = {}
    rows = "ABCDEFGH"

    try:
        # Transform first col in df_conc to list; errors='coerce' turns non-numbers to NaN
        vals = pd.to_numeric(df_conc.iloc[:, 0], errors='coerce').tolist()

        for i, row_char in enumerate(rows):
            if i < len(vals):
                # Store float if valid, else 0.0 (or None if preferred)
                conc_dic[row_char] = float(vals[i]) if not pd.isna(vals[i]) else 0.0
            else:
                conc_dic[row_char] = 0.0 # For vehicle row
    except Exception as e:
        print(f"   [WARNING] Error parsing concentration table: {e}")

    return conc_dic

def get_ligand_map(protocol: ProtocolData):
    """
    Determines which columns contain Ligand 1 and which contain Ligand 2.
    Returns dic of col idx and ligands.
    """
    ligand_map = {}

    # If no second ligand exists, everything is Ligand 1
    if not protocol.ligand_2:
        for c in range(1, 13): ligand_map[c] = 'L1'
        return ligand_map

    # Determine Ligand Layout Strategy based on Cell Layout
    cell_layout = str(protocol.line_layout).lower()

    if "one line" in cell_layout:
        # STRICT RULE: One line layout cannot support 2 ligands in this logic
        print(f"   [ERROR] Protocol '{protocol.file_name}' lists 2 ligands but uses 'One Line' cell layout.")
        print(f"           This configuration is not supported. Defaulting all columns to Ligand 1 {protocol.ligand}.")
        for c in range(1, 13): ligand_map[c] = 'L1'
        return ligand_map

    elif "half" in cell_layout:
        # Cell Layout: Half (1-6 / 7-12) -> Ligand Layout: Alternating Blocks
        # L1: 1-3, 7-9 | L2: 4-6, 10-12
        l1_cols = list(range(1, 4)) + list(range(7, 10))
        l2_cols = list(range(4, 7)) + list(range(10, 13))
        for c in l1_cols: ligand_map[c] = 'L1'
        for c in l2_cols: ligand_map[c] = 'L2'

    elif "alternating" in cell_layout:
        # Cell Layout: Alternating Blocks -> Ligand Layout: Half/Half
        # L1: 1-6 | L2: 7-12
        for c in range(1, 7): ligand_map[c] = 'L1'
        for c in range(7, 13): ligand_map[c] = 'L2'

    else:
        print(f"   [WARNING] Unknown cell layout '{cell_layout}'. Defaulting all to Ligand 1 {protocol.ligand}.")
        for c in range(1, 13): ligand_map[c] = 'L1'

    return ligand_map

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

def convert_to_plate_layout(data_input) -> pd.DataFrame:
    """
    Converts a 1-row df OR a dictionary of {Well_ID: Value} as in AUC data
    into a pandas df representing a 96-well plate (Rows A-H, Cols 1-12).
    """
    # Handle input types:
    # If df (like auc_raw_df), convert 1st row to dict
    if isinstance(data_input, pd.DataFrame):
        if data_input.empty:
            return pd.DataFrame()
        data_dict = data_input.iloc[0].to_dict()
    # If a series, convert to dict
    elif isinstance(data_input, pd.Series):
        data_dict = data_input.to_dict()
    else:
        data_dict = data_input

    if not data_dict: return pd.DataFrame()

    # Check format by checking the first key (well ids as row headers or also conditions)
    first_key = str(list(data_dict.keys())[0])

    rows = list("ABCDEFGH")

    # AUC mean data (header cond|cell|row)
    if "|" in first_key:
        # create {row_char: {col_Header: value}}
        reshaped_data = {r: {} for r in rows}
        for key, value in data_dict.items():
            # "cond|cell|A" -> ["cond|cell", "A"]
            parts = str(key).rsplit('|', 1)

            if len(parts) == 2:
                col_header = parts[0]  # The name without the row letter
                row_char = parts[1]  # The row letter (A, B, etc.)

                if row_char in rows:
                    reshaped_data[row_char][col_header] = value

        # Create df (index=A-H, cols=conditions)
        df_mean = pd.DataFrame.from_dict(reshaped_data, orient='index')
        return df_mean

    # AUC data in other processing steps (header A1, B2...) ---
    else:
        cols = list(range(1, 13))

        # Initialize empty DataFrame with NaN
        plate_df = pd.DataFrame(None, index=rows, columns=cols)

        for well_id, value in data_dict.items():
            # Skip if header is not a string (safety)
            if not isinstance(well_id, str) or len(well_id) < 2:
                continue

            r = well_id[0].upper()
            # Try-except block handles headers that aren't well IDs (like "Time")
            try:
                c = int(well_id[1:])
                # Assign value if coordinates are valid
                if r in rows and c in cols:
                    plate_df.at[r, c] = value
            except ValueError:
                continue

        return plate_df

def calculate_relative_time(raw_time_col: pd.Series, baseline_end_idx: int):
    """
    Calculates a relative time vector. Uses the measuring interval of the kinetic reading to
    set first measurement after baseline (baseline_end_idx) to 0. Negative time for baseline reads.
    """
    # Clean and convert to numeric
    times = pd.Series(pd.to_numeric(raw_time_col, errors='coerce'))
    kinetic_times = times.iloc[baseline_end_idx:].dropna()

    if len(kinetic_times) < 2: return None  # Not enough data points

    # Take the difference to filter out potential jitter
    # -> Rounding is important since otherwise 1.00 and 1.0 are not considered the same
    interval = kinetic_times.diff().dropna().round(2).unique()[0]

    # Generate time vector for baseline and kinetic reading
    n_rows = len(raw_time_col)
    time_vector = []

    for i in range(n_rows):
        # (current_index - zero_index) * interval
        t = (i - baseline_end_idx) * interval
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
                            acc_range: float,
                            outlier_dict: dict,
                            kinetic_reads_count: int = 1):
    """
    Calculates the mean of the vehicle wells (row H) for each block.
    Works for both Kinetic DataFrames (returns Series mean) and AUC DataFrames (returns Float mean).
    """
    vehicle_means = {}

    # Check if this is AUC (1 row) or Kinetic (>1 row)
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

            # --- VEHICLE CHECK ---
            for well in valid_vehicles:
                if is_kinetic:
                    # Kinetic: Mean over time should be close to 1.0
                    val_to_check = vehicle_data[well].mean()
                    target_value = 1
                else:
                    # AUC: Value should be close to (1.0 * number_of_reads)
                    # Use .iloc[0] to get the float from the series
                    val_to_check = vehicle_data[well].iloc[0]
                    target_value = float(kinetic_reads_count)

                threshold = acc_range * target_value
                # Check deviation (only if value is not NaN)
                if pd.notna(val_to_check) and abs(val_to_check - target_value) > threshold:
                    outlier_dict[well] = float(val_to_check)

            # Calculate mean across the valid wells
            if is_kinetic:
                # Row-wise mean for kinetic traces (result: series of length = timepoints)
                vehicle_means[start_col] = vehicle_data.mean(axis=1)
            else:
                # Scalar mean for AUC (result: single float)
                vehicle_means[start_col] = vehicle_data.mean(axis=1).iloc[0]
        else:
            vehicle_means[start_col] = None

    return vehicle_means

def calculate_replicate_means(processed_df: pd.DataFrame,
                              plate_blocks: list[range],
                              col_metadata: dict): # Dic created from PlateColMetadata
    """
    Calculates the mean of technical replicates. As in processed_df each col is one well,
    the mean is performed of three cols within one block.
    Rows (A-H) are treated as distinct conditions (ligand concentration).
    Header format of returned df: "Condition_Name|Cell_Line|Ligand_Name|Row"
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

    return pd.DataFrame(mean_data)

def process_bret_measurement(result: PrResult, protocol: ProtocolData):
    """
    Maps cell line x transfection x ligand plate layout using protocol info.
    Performs baseline correction. Vehicle normalisation with kinetic data and AUC in parallel.
    Checks vehicle for outliers.
    """
    if result.is_excluded:
        result.kinetic_df = None
        result.kinetic_mean_df = None
        result.auc_df = None
        result.auc_mean_df = None
        return result

    # Reset for re-run
    result.warnings = []
    result.vehicle_outliers = {}
    result.column_metadata = {}

    print(f"\n[DEBUG] === Processing File: {result.file_name} ===")
    # --- CONFIG LAYOUT ---
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

    # Ligand identity and conc map
    ligand_col_map = get_ligand_map(protocol)
    conc_map_1 = built_conc_dic(protocol.ligand_conc)
    # Only built if second ligand is defined
    conc_map_2 = built_conc_dic(protocol.ligand_2_conc) if protocol.ligand_2 else {}

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

            # Check the placeholder token ('L1' or 'L2')
            which_lig = ligand_col_map.get(col, 'L1')

            if which_lig == 'L2' and protocol.ligand_2:
                current_ligand_name = str(protocol.ligand_2)
                current_conc_map = conc_map_2
            else:
                current_ligand_name = str(protocol.ligand)
                current_conc_map = conc_map_1

            result.column_metadata[col] = PlateColMetadata(
                cell_line=c_line,
                transfection_id=t_id,
                condition_name=current_cond_name,
                plasmids=current_plasmids,
                ligand_identity=current_ligand_name,
                ligand_conc=current_conc_map
            )

    # --- CONFIG PROCESSING ---
    # Define accepted vehicle range
    acc_vehicle_range = 0.2
    # Define Baseline: First 5 rows (Index 0-4)
    baseline_end_idx = 5

    # Get raw BRET ratio table
    raw_df = result.raw_bret_ratio_df.copy()
    time_col = "Time (min)"
    if len(raw_df) < baseline_end_idx:
        print(f"   [WARNING] Data has less than {baseline_end_idx} rows.")
        return result

    # Built time vector
    time_vec = calculate_relative_time(raw_df["Time (min)"], baseline_end_idx)
    if time_vec is None:
        time_vec = range(len(raw_df))  # Fallback index

    # Prepare the full dataframe for normalization (removing the time col)
    data_df = raw_df.drop(columns=[time_col]).copy()

    # --- BASELINE CORRECTION ---
    # Isolate baseline rows for calculation
    df_baseline_calc = raw_df.iloc[0:baseline_end_idx].copy()
    bl_corrected_df = pd.DataFrame()

    for col in data_df.columns:
        if col in result.excluded_wells:
            bl_corrected_df[col] = float('nan')
            continue
        # Mean of baseline rows for this well
        # Force numeric conversion for baseline values to handle potential strings/decimals
        base_mean = pd.to_numeric(df_baseline_calc[col], errors='coerce').mean()

        if base_mean != 0:
            col_vals = pd.to_numeric(data_df[col], errors='coerce')
            # Divide ALL values by the baseline mean
            bl_corrected_df[col] = col_vals / base_mean
        else:
            bl_corrected_df[col] = float('nan')

    # --- AUC CALCULATION ---
    # Use slicing to sum only the kinetic phase (after baseline)
    auc_raw_series = bl_corrected_df.iloc[baseline_end_idx:].sum(axis=0, skipna=False)
    auc_raw_df = pd.DataFrame([auc_raw_series]) # Convert series to df

    # Determine number of reads for dynamic AUC target
    num_kinetic_points = len(bl_corrected_df.iloc[baseline_end_idx:])

    # --- VEHICLE CORRECTION ---
    # Kinetics (df -> returns series of means over time)
    veh_means_kinetic = calculate_vehicle_means(
        bl_corrected_df, plate_blocks, result.excluded_wells, acc_vehicle_range, result.vehicle_outliers,
        kinetic_reads_count=1
    )
    # For AUC (df with one row -> returns dictionary of scalars)
    veh_means_auc = calculate_vehicle_means(
        auc_raw_df, plate_blocks, result.excluded_wells, acc_vehicle_range, result.vehicle_outliers,
        kinetic_reads_count=num_kinetic_points
    )

    # Apply Normalization (using dictionaries to avoid fragmentation/warnings for pd.Df)
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
            auc_norm_dict[col] = auc_raw_df[col] / v_auc
        else:
            auc_norm_dict[col] = float('nan')

    # Create df from dict
    kinetic_df = pd.DataFrame(kinetic_norm_dict)
    auc_df = pd.DataFrame(auc_norm_dict)

    # --- MEAN OF REPLICATES (KINETIC) ---
    kinetic_mean_df = calculate_replicate_means(
        kinetic_df, plate_blocks, result.column_metadata
    )
    auc_mean_df = calculate_replicate_means(
        auc_df, plate_blocks, result.column_metadata
    )

    # --- SAVE RESULTS ---
    # --- Kinetic data
    result.time_vector = time_vec
    result.bl_corr_kinetic = bl_corrected_df
    result.kinetic_df = pd.DataFrame(kinetic_norm_dict)
    result.kinetic_mean_df = kinetic_mean_df
    # --- AUC data
    result.raw_auc_df = auc_raw_df
    result.raw_auc_tidy_df = convert_to_plate_layout(auc_raw_df)
    result.auc_df = pd.DataFrame(auc_norm_dict)
    result.auc_tidy_df = convert_to_plate_layout(result.auc_df)
    result.auc_mean_df = auc_mean_df
    result.auc_mean_tidy_df = convert_to_plate_layout(result.auc_mean_df)

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

        # --- INTERNAL STORAGE ---
        # Path to folder variable
        self.folder_path = tk.StringVar(value="No folder selected.")
        self.subfolder_paths_with_files = []
        self.experiment: list[MeasurementFolder] = []
        self.master_index = pd.DataFrame() # Index for populating tab 2
        self.rule_history_text = ""
        self.master_df = pd.DataFrame # Used for master csv file storage (by generation or import)

        # --- TABS SETUP ---
        self.notebook = ttk.Notebook(main_window)
        self.notebook.pack(expand=True, fill='both')

        # Tab 1: Import Data
        self.tab_import = tk.Frame(self.notebook)
        self.notebook.add(self.tab_import, text="Import & Export Data")

        # Tab 2: Data Selection
        self.tab_select = tk.Frame(self.notebook)
        self.notebook.add(self.tab_select, text="Optional Selection")

        # Tab 3: Plot Helper
        self.tab_plot_helper = tk.Frame(self.notebook)
        self.notebook.add(self.tab_plot_helper, text="Plot Helper")
        # --- TAB 3 CONTENT ---
        # Setup immediately for optional import of master csv
        self.setup_plot_helper_tab()

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
        # Label to display Main Plasmids
        self.main_plasmids_label = tk.Label(self.tab_import,
                                            text="",
                                            justify="left",
                                            font=("Arial", 10, "bold"))
        self.main_plasmids_label.pack(pady=(0, 5))

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

        # --- EXPORT SECTION ---
        export_frame = tk.LabelFrame(self.tab_import, text="Export Options")
        export_frame.pack(fill="x", padx=20, pady=10)
        # Master CSV button
        self.btn_export_master = tk.Button(export_frame,
                                           text="Export Master CSV",
                                           state="disabled",
                                           command=self.export_master_csv)
        self.btn_export_master.pack(side="left", fill="x", expand=True, padx=5, pady=10)

        # Preview xlsx button
        self.btn_export_excel = tk.Button(export_frame,
                                          text="Export Excel Report (Default)",
                                          state="disabled",
                                          command=self.export_excel_report)
        self.btn_export_excel.pack(side="left", fill="x", expand=True, padx=5, pady=10)

        # --- TAB 2 CONTENT ---
        self.setup_exclusion_tab()

    # Helper function to write to that text box
    def log(self, message):
        """Logs to the separate window"""
        try:
            self.log_text.insert(tk.END, message + "\n")
            self.log_text.see(tk.END)

            # FORCE GUI UPDATE: important for filling log during processing
            self.main_gi.update_idletasks()
        except tk.TclError:
            # Handle case where user manually closed log window but app is running
            print(message)

    def update_summary_table(self):
        """Fills the summary table with N counts and dates per condition"""
        # Clear existing data
        for i in self.summary_tree.get_children():
            self.summary_tree.delete(i)

        if self.master_index.empty:
            return

        # Group by Cell Line and Condition
        grouped = self.master_index.groupby(['Cell_Line', 'Condition'])

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
        if self.master_index.empty: return

        # Get unique dates and add "All"
        dates = sorted(self.master_index['Date'].unique().tolist())
        self.cb_date['values'] = ["All"] + dates
        self.var_date.set("All")

        # Reset others
        self.update_cell_options()

    def update_cell_options(self, event=None):
        """Updates Cell Line options based on selected Date."""
        selected_date = self.var_date.get()

        if self.master_index.empty: return

        if selected_date == "All":
            # Show all cell lines available in the whole dataset
            cells = sorted(self.master_index['Cell_Line'].unique().tolist())
        else:
            # Filter DF by date
            subset = self.master_index[self.master_index['Date'] == selected_date]
            cells = sorted(subset['Cell_Line'].unique().tolist())

        self.cb_cell['values'] = ["All"] + cells
        self.var_cell.set("All")
        self.update_cond_options()

    def update_cond_options(self, event=None):
        """Updates Condition options based on selected Date AND Cell Line."""
        selected_date = self.var_date.get()
        selected_cell = self.var_cell.get()

        if self.master_index.empty: return

        # Start with full DF
        subset = self.master_index.copy()

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
            df = self.master_index.copy()

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

    def setup_plot_helper_tab(self):
        """Builds the GUI for tab 3 plot helper"""

        # --- Import Master CSV ---
        src_frame = tk.Frame(self.tab_plot_helper)
        src_frame.pack(fill="x", padx=10, pady=10)
        tk.Label(src_frame, text="Data Source:", font=("Arial", 9, "bold")).pack(side="left")

        self.lbl_data_source = tk.Label(src_frame, text="No Data Loaded")
        self.lbl_data_source.pack(side="left", padx=10)

        tk.Button(src_frame, text="Import Master CSV",
                  command=self.import_master_csv).pack(side="right")

        # --- Filter Selection ---
        sel_frame = tk.LabelFrame(self.tab_plot_helper, text="Select Data to Include")
        sel_frame.pack(fill="both", expand=True, padx=10, pady=5)

        # Ligand Listbox
        tk.Label(sel_frame, text="Ligands:").grid(row=0, column=0, padx=5, sticky="w")
        self.lb_ligands = tk.Listbox(sel_frame, selectmode="multiple", height=6, exportselection=False)
        self.lb_ligands.grid(row=1, column=0, padx=3, pady=5, sticky="nsew")

        # Cell Lines Listbox
        tk.Label(sel_frame, text="Cell Lines:").grid(row=0, column=1, padx=5, sticky="w")
        self.lb_exp_cells = tk.Listbox(sel_frame, selectmode="multiple", height=6, exportselection=False)
        self.lb_exp_cells.grid(row=1, column=1, padx=5, pady=5, sticky="nsew")

        # Transfections Listbox
        tk.Label(sel_frame, text="Conditions (Transfections):").grid(row=0, column=2, padx=5, sticky="w")
        self.lb_exp_trans = tk.Listbox(sel_frame, selectmode="multiple", height=6, exportselection=False)
        self.lb_exp_trans.grid(row=1, column=2, padx=5, pady=5, sticky="nsew")
        # Add scrollbar
        trans_scroll = tk.Scrollbar(sel_frame, orient="vertical", command=self.lb_exp_trans.yview)
        trans_scroll.grid(row=1, column=3, sticky="ns", pady=5)
        self.lb_exp_trans.config(yscrollcommand=trans_scroll.set) # Update scrollbar

        # Select All Buttons
        tk.Button(sel_frame, text="Select All Ligands", command=lambda: self.select_all_listbox(self.lb_ligands)).grid(
            row=2, column=0)
        tk.Button(sel_frame, text="Select All Cells", command=lambda: self.select_all_listbox(self.lb_exp_cells)).grid(
            row=2, column=1)
        tk.Button(sel_frame, text="Select All Transf.", command=lambda: self.select_all_listbox(self.lb_exp_trans)).grid(
            row=2, column=2)

        sel_frame.columnconfigure(0, weight=1)
        sel_frame.columnconfigure(1, weight=2)
        sel_frame.columnconfigure(2, weight=3)

        # --- Export Data Options ---
        opt_frame = tk.LabelFrame(self.tab_plot_helper, text="Export Data Options")
        opt_frame.pack(fill="x", padx=10, pady=5)

        # Kinetic Options
        tk.Label(opt_frame, text="Kinetic Data:").grid(row=0, column=0, sticky="w", padx=10)
        self.var_exp_kin = tk.StringVar(value="Row A (Max)")
        combo_kin = ttk.Combobox(opt_frame, textvariable=self.var_exp_kin, state="readonly")
        combo_kin['values'] = ["None", "Row A (Max)", "All Rows"]
        combo_kin.grid(row=0, column=1, padx=5, pady=5)

        # AUC Options
        tk.Label(opt_frame, text="AUC Data:").grid(row=1, column=0, sticky="w", padx=10)
        self.var_exp_auc = tk.StringVar(value="Conc Response")
        combo_auc = ttk.Combobox(opt_frame, textvariable=self.var_exp_auc, state="readonly")
        combo_auc['values'] = ["None", "Conc Response"]
        combo_auc.grid(row=1, column=1, padx=5, pady=5)

        # --- 3. Action Button ---
        btn_frame = tk.Frame(self.tab_plot_helper)
        btn_frame.pack(fill="x", padx=10, pady=20)

        self.btn_run_plot_helper = tk.Button(btn_frame, text="Generate Custom Export",
                                            state="disabled", command=self.run_plot_helper)
        self.btn_run_plot_helper.pack(fill="x", ipady=5)

    def select_all_listbox(self, lb):
        lb.select_set(0, tk.END)

    def refresh_plot_helper_options(self):
        """Populates the list boxes in the plot helper from master df (opt. imported csv file)."""
        if self.master_df is None or self.master_df.empty:
            self.btn_run_plot_helper.config(state="disabled")
            self.lbl_data_source.config(text="No Data")
            return

        if not self.experiment:
            # CSV mode -> status already set
            pass
        else:
            #  Fresh Analysis mode
            experiment = self.master_df['Main_Plasmids'].unique()[0]
            self.lbl_data_source.config(text=f"Internal: {experiment}")

        self.btn_run_plot_helper.config(state="normal")

        # Clear
        self.lb_ligands.delete(0, tk.END)
        self.lb_exp_cells.delete(0, tk.END)
        self.lb_exp_trans.delete(0, tk.END)

        df = self.master_df

        # Populate Ligand
        ligands = sorted(set(df['Ligand'].astype(str)))
        for l in ligands: self.lb_ligands.insert(tk.END, l)
        self.select_all_listbox(self.lb_ligands) # Default to all

        # Populate Cells
        cells = sorted(set(df['Cell_Line'].astype(str)))
        for c in cells: self.lb_exp_cells.insert(tk.END, c)
        self.select_all_listbox(self.lb_exp_cells)  # Default to all

        # Populate Transfections
        trans = sorted(set(df['Transfection'].astype(str)))
        for t in trans: self.lb_exp_trans.insert(tk.END, t)
        self.select_all_listbox(self.lb_exp_trans)  # Default to all

    def run_plot_helper(self):
        """Collects GUI selections and calls the core export engine."""
        if self.master_df is None or self.master_df.empty: return

        # Get Selections
        cells = [self.lb_exp_cells.get(i) for i in self.lb_exp_cells.curselection()]
        transfections = [self.lb_exp_trans.get(i) for i in self.lb_exp_trans.curselection()]
        ligands = [self.lb_ligands.get(i) for i in self.lb_ligands.curselection()]

        config = {
            'cells': cells,
            'transfections': transfections,
            'ligands': ligands,
            'kinetic_mode': self.var_exp_kin.get(),
            'auc_mode': self.var_exp_auc.get()
        }

        file_path = filedialog.asksaveasfilename(
            defaultextension=".xlsx",
            filetypes=[("Excel", "*.xlsx")],
            title="Save Custom Export"
        )
        if not file_path: return

        # Call the core engine with the unified DF
        self.write_excel_export(file_path, self.master_df, config)

    def import_master_csv(self):
        """Loads a master csv file directly into the memory for the plot helper."""
        file_path = filedialog.askopenfilename(
            title="Select Master CSV File",
            filetypes=[("CSV Files", "*.csv"), ("All Files", "*.*")]
        )
        if not file_path: return

        try:
            df = pd.read_csv(file_path)

            # Validation
            required = ["Transfection", "Cell_Line", "Ligand", "Kinetic_Mean"]
            if not all(col in df.columns for col in required):
                self.log("[ERROR] Invalid CSV format. Columns missing.")
                return

            # Store in the unified variable
            self.master_df = df

            # Reset raw data references so we know we are in "CSV Mode"
            self.experiment = []
            self.master_index = pd.DataFrame()  # Clear exclusion index

            # Update GUI
            self.lbl_data_source.config(text=f"CSV: {os.path.basename(file_path)}")
            self.refresh_plot_helper_options()
            self.log(f"Loaded Master CSV: {os.path.basename(file_path)}")

        except Exception as e:
            self.log(f"[ERROR] CSV Load Failed: {e}")

    def select_folder(self):
            """Opens dialog to select folder to search for xlsx files in"""
            directory = filedialog.askdirectory(title="Select a folder...")
            if not directory: return

            # --- RESET STATE: Clear old data when a new folder is selected
            self.subfolder_paths_with_files = []
            self.experiment = []
            self.master_index = pd.DataFrame()
            # Reset summary table
            for i in self.summary_tree.get_children():
                self.summary_tree.delete(i)
            self.rule_history_text = ""
            self.main_plasmids_label.config(text="")
            self.lbl_rules_summary.config(text="")
            self.clear_exclusion_list()
            self.btn_export_master.config(state="disabled")
            self.btn_export_excel.config(state="disabled")


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

        # Collect list of records for each col
        records = []
        rows_str = "ABCDEFGH"

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
            self.master_index = pd.DataFrame(records)

            # --- Summary for verification ---
            summary = self.master_index.groupby(['Cell_Line', 'Condition'])['File_Name'].nunique()
            print("\n[DEBUG] Data Summary:\n", summary)
            self.refresh_filter_options()
            self.update_summary_table()
            self.refresh_plot_helper_options()
            return self.master_index
        else:
            self.master_index = pd.DataFrame()
            self.update_summary_table()
            print("No valid data found")
            return self.master_index

    def collect_files(self):
        """
        1. LOAD FILES
        Reads sheet names of all xlsx and xlsm files to identify and separate protocol and result analysis files.
        Validation of correct protocol to analysis files is done via date of measurement in the folder name.
        """
        if not self.subfolder_paths_with_files:
            print("No folders to analyze.")
            return

        # Reset exclusion state
        self.rule_history_text = ""
        self.lbl_rules_summary.config(text="")
        self.clear_exclusion_list()

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
        self.btn_export_master.config(state="normal")
        self.btn_export_excel.config(state="normal")
        self.btn_run_plot_helper.config(state="normal")

    def run_processing_pipeline(self):
        """
        2. Processing of raw BRET data and indexing with protocol info
        Calls processing functions and is rerun if data was excluded.
        """
        self.log("\n--- Starting Processing Pipeline ---")

        # Check for plasmids transfected in all conditions (main plasmids) and filter if needed
        selected_exp_name = self.handle_main_plasmids_selection()
        self.main_plasmids_label.config(text=f"{selected_exp_name}")

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

        self.log("--- Compiling Master Dataframe... ---")
        self.master_df = self.compile_master_dataframe()
        self.refresh_plot_helper_options()

        self.log("\n--- Processing Complete & Plot Helper Ready ---")

    def compile_master_dataframe(self):
        """
        Compiles technical means (Kinetic & AUC) into a tidy Master DataFrame.
        Structure: Long format (1 row per timepoint).
        The scalar AUC value is REPEATED for every timepoint of the same condition.
        """
        if not self.experiment:
            return None
        self.log("\n--- Building Master CSV ---")
        master_rows = []

        for folder in self.experiment:
            if folder.protocol.main_plasmids:
                main_plasmids = main_plasmids_str = " + ".join(folder.protocol.main_plasmids)
            else:
                main_plasmids = "Unknown"

            for res in folder.results:
                if res.is_excluded or res.kinetic_mean_df is None: continue
                # If time_vector is missing, create a generic index
                time_points = res.time_vector if res.time_vector else range(len(res.kinetic_mean_df))

                # Get meta information for each col via PlateColMetadata stored in res.column_metadata
                meta_lookup = {}
                for meta in res.column_metadata.values():
                    # Store key as unique combo
                    key = (meta.condition_name, meta.ligand_identity)
                    if key not in meta_lookup:
                        meta_lookup[key] = meta

                # Iterate through the KEYS of the mean DataFrame
                # Key format from 'calculate_replicate_means': "Condition|Cell_Line|Ligand_Name|Row"
                for col_key in res.kinetic_mean_df.columns:
                    try:
                        parts = col_key.split('|')
                        if len(parts) != 4:
                            print(f"[WARNING] Skipping {col_key}: Format expected 4 parts, got {len(parts)}")
                            continue
                        transfection, cell_line, lig_name, row_char = parts
                    except ValueError:
                        continue

                    # --- RETRIEVE DATA ---
                    # Kinetic Mean Series (Vector)
                    # This list has length = number of timepoints
                    kin_mean_values = res.kinetic_mean_df[col_key].tolist()

                    # AUC Mean Value (Scalar)
                    # Check if this key exists in the AUC Mean DF
                    if res.auc_mean_df is not None and col_key in res.auc_mean_df.columns:
                        auc_mean_val = res.auc_mean_df[col_key].iloc[0]
                    else:
                        print(f"[WARNING] master df compilation: {col_key} "
                              f"no respective AUC mean found")
                        auc_mean_val = float('nan')

                    # Get ligand concentrations
                    target_meta = meta_lookup.get((transfection, lig_name))
                    if target_meta:
                        # Since 'ligand_conc' is already specific to this column (as you confirmed),
                        # we just grab the value for this row.
                        conc_val = target_meta.ligand_conc.get(row_char, 0.0)

                    else:
                        print(f"[WARNING] master df compilation: {col_key} "
                              f"no respective ligand concentration found")
                        conc_val = float('nan')

                    # --- BUILD ROWS (TIDY FORMAT) ---
                    # Zip timepoints with kinetic values
                    for t_val, kin_val in zip(time_points, kin_mean_values):
                        row = {
                            # --- Identifiers ---
                            "File_Name": res.file_name,
                            "Date": res.measurement_date,
                            "Main_Plasmids": main_plasmids,
                            "Cond_Key": col_key,  # Unique ID for this curve
                            "Time_(min)": t_val,

                            # --- Metadata ---
                            "Transfection": transfection,
                            "Cell_Line": cell_line,
                            "Plate_Row": row_char,
                            "Ligand": lig_name,
                            "Ligand_Conc": conc_val,

                            # --- The Data ---
                            "Kinetic_Mean": kin_val,
                            "AUC_Mean": auc_mean_val
                        }
                        master_rows.append(row)

        # Create DataFrame
        df_master = pd.DataFrame(master_rows)
        return df_master

    def export_master_csv(self):
        """Saves compiled master df to csv"""
        if self.master_df is None or self.master_df.empty:
            self.log("No data to export.")
            return

        # Native Dialog handles overwrite warning automatically
        file_path = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV File", "*.csv")],
            title="Save Master CSV"
        )
        if not file_path: return

        try:
            self.master_df.to_csv(file_path, index=False)
            self.log(f"   [SUCCESS] Saved Master CSV: {os.path.basename(file_path)}")
        except Exception as e:
            self.log(f"   [ERROR] Failed to save CSV: {e}")

    def write_excel_export(self, file_path, master_df, config):
        """
        Writes the Excel file based on the config dictionary provided by either tab 1 (default )or tab 3 (user).

        Config Keys:
          - 'cells': list of cell lines to include (or 'All')
          - 'transfections': list of transfections to include (or 'All')
          - 'kinetic_mode': 'None', 'Row A (Max)', or 'All Rows'
          - 'auc_mode': 'None' or 'Conc Response'
        TODO: reformat AUC
        TODO: add which processing step to include -> as mean of techn. replicates then
        TODO: these processing steps would then also have to be included in master csv!
        TODO: add arranging config (group by x) -> keeping in mind opt second ligand
        """
        try:
            df_subset = master_df.copy()

            if config.get('ligands') != 'All':
                df_subset = df_subset[df_subset['Ligand'].isin(config['ligands'])]

            if config.get('cells') != 'All':
                df_subset = df_subset[df_subset['Cell_Line'].isin(config['cells'])]

            if config.get('transfections') != 'All':
                df_subset = df_subset[df_subset['Transfection'].isin(config['transfections'])]

            if df_subset.empty:
                print("[ERROR] Export failed: Filter resulted in no data.")
                return

            with pd.ExcelWriter(file_path) as writer:

                # --- 1. METADATA SHEET ---
                # Include always
                file_names = df_subset["File_Name"].unique().tolist()
                mp_str = "Unknown"
                if 'Main_Plasmids' in df_subset.columns:
                    mp_vals = df_subset['Main_Plasmids'].unique()
                    if len(mp_vals) > 0: mp_str = mp_vals[0]

                meta_dict = {
                    "Export Date": [datetime.now().strftime("%d.%m.%Y - %H:%M:%S")],
                    "Main Plasmids": [mp_str],
                    "Source Files Count": [len(file_names)],
                    "Source Files List": [", ".join(file_names)],
                    "Filter: Ligands": [", ".join(config.get('ligands'))],
                    "Filter: Cells": [", ".join(config.get('cells'))],
                    "Filter: Conditions": [", ".join(config.get('transfections'))]
                }
                pd.DataFrame(meta_dict).transpose().to_excel(writer, sheet_name="Metadata", header=False)

                # --- 2. KINETIC DATA ---
                k_mode = config.get('kinetic_mode', 'None')
                if k_mode != 'None':
                    # Filter Rows based on mode
                    if k_mode == 'Row A (Max)':
                        df_kin = df_subset[df_subset["Plate_Row"] == "A"].copy()
                        sheet_prefix = "Kinetic_Max"
                    else:
                        df_kin = df_subset.copy()
                        sheet_prefix = "Kinetic_All"

                    if not df_kin.empty:
                        # Create Header Key
                        df_kin["Header_Key"] = df_kin["Transfection"] + " | " + df_kin["Cell_Line"] + " | " + df_kin[
                            "Ligand"]
                        if k_mode == 'All Rows':
                            df_kin["Header_Key"] += " | " + df_kin["Plate_Row"]

                        # Merge cells of equal header
                        kin_pivot = df_kin.pivot_table(
                            index="Time_(min)",
                            columns=["Header_Key", "File_Name"],
                            values="Kinetic_Mean"
                        )

                        # Format Headers creating the empty headers
                        new_headers = []
                        last_key = None
                        for key, file_name in kin_pivot.columns:
                            if key != last_key:
                                new_headers.append(key)
                                last_key = key
                            else:
                                new_headers.append("")

                        kin_pivot.columns = new_headers
                        kin_pivot.reset_index(inplace=True)
                        kin_pivot.rename(columns={"Time_(min)": "Time (min)"}, inplace=True)

                        # Save
                        kin_pivot.to_excel(writer, sheet_name=sheet_prefix, index=False)

                # --- 3. AUC DATA ---
                a_mode = config.get('auc_mode', 'None')
                print("___AUC HERE_____")

                if a_mode == 'Conc Response':
                    # Drop duplicates for Scalar AUC
                    df_auc = df_subset.drop_duplicates(subset=["File_Name", "Transfection", "Cell_Line"]).copy()

                    if not df_auc.empty:
                        df_auc["Header_Key"] = df_auc["Transfection"] + " | " + df_auc["Cell_Line"] + " | " + df_auc[
                            "Ligand"]

                        auc_pivot = df_auc.pivot_table(
                            index="Ligand_Conc",
                            columns=["Header_Key", "File_Name"],
                            values="AUC_Mean"
                        )

                        # Format Headers
                        new_headers = []
                        last_key = None
                        for key, file_name in auc_pivot.columns:
                            if key != last_key:
                                new_headers.append(key)
                                last_key = key
                            else:
                                new_headers.append("")

                        auc_pivot.columns = new_headers
                        # auc_pivot.sort_index(inplace=True)
                        # auc_pivot.reset_index(inplace=True)
                        auc_pivot.rename(columns={"Ligand_Conc": "Concentration (logM)"}, inplace=True)

                        auc_pivot.to_excel(writer, sheet_name="AUC_Conc_Response", index=False)

            self.log(f"   [SUCCESS] Exported: {os.path.basename(file_path)}")

        except Exception as e:
            self.log(f"   [ERROR] Export failed: {e}")
            print(e)

    def export_excel_report(self):
        """
        Default Export: All Data, Highest stim kinetics and AUC crc.
        TODO: Grouped for transfection as default
        """
        if self.master_df is None or self.master_df.empty:
            self.log("No data found to export.")
            return

        file_path = filedialog.asksaveasfilename(
            defaultextension=".xlsx",
            filetypes=[("Excel File", "*.xlsx")],
            title="Save Standard Report"
        )
        if not file_path: return

        # Define Standard Config
        default_config = {
            'cells': 'All',
            'transfections': 'All',
            'ligands': 'All',
            'kinetic_mode': 'Row A (Max)',
            'auc_mode': 'Conc Response'
        }

        self.write_excel_export(file_path, self.master_df, default_config)

# TODO: implement showing also errors from tool functions in log window

# --- Main Execution Block ---
if __name__ == "__main__":
    # Create the main window
    root = tk.Tk()

    # Create an instance of the application
    app = NCollectorApp(root)

    # Start the Tkinter event loop
    root.mainloop()