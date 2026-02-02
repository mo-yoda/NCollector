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
    replicate: int = 0

@dataclass
class PrResult:
    """ Raw and processed information from a single _analysis file """
    # --- Information from analysis xlsx itself ---
    file_name: str
    measurement_date: date
    cell_line: str # ID2
    transfection_id: str # ID3
    raw_bret_ratio_df: pd.DataFrame
    lum_df: pd.DataFrame

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
    # Raw BRET from last 3x datapoints
    raw_bret_points_df: pd.DataFrame | None = None
    # Pre-vehicle norm AUC
    bl_corr_auc_df: pd.DataFrame | None = None
    # Baseline- and vehicle-normalised AUC data (technical replicates)
    auc_df: pd.DataFrame | None = None
    # Mean of baseline- and vehicle-normalised AUC data
    auc_mean_df:pd.DataFrame | None = None

    # --- Export tidy CRC data ---
    raw_bret_points_tidy_df: pd.DataFrame | None = None
    bl_corr_auc_tidy_df: pd.DataFrame | None = None
    auc_tidy_df: pd.DataFrame | None = None
    auc_mean_tidy_df: pd.DataFrame | None = None

    # --- Optional exclusion by user interaction  ---
    is_excluded: bool = False
    excluded_wells: list[str] = field(default_factory=list)

    # --- Internal check and warnings for outlier identification ---
    vehicle_outliers: dict[str, float] = field(default_factory=dict) # well, value
    vehicle_warnings: list[dict] = field(default_factory=list)
    low_lum_warnings: list[dict] = field(default_factory=list) # List of dict carrying all needed metadata

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

@dataclass
class ProcessingConfig:
    """All user-defined processing settings"""
    lum_threshold: int=100
    vehicle_warning_threshold: float = 0.2
    baseline_end_index: int=5
    plate_layout: list[range] = field(default_factory=lambda: [
        range(1, 4),  # Block 1: Cols 1-3
        range(4, 7),  # Block 2: Cols 4-6
        range(7, 10),  # Block 3: Cols 7-9
        range(10, 13)  # Block 4: Cols 10-12
    ])

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
        match_index: int = 0,
        row_end_threshold: int = 1 # Needed for extracting bret data
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
                       if str(val).lower() == "nan" or not str(val).strip()), len(header_row_content))

    # Calculate the absolute end column
    col_end_idx = col_start_idx + col_width

    # Extract the table from large df
    df_extract = df.iloc[header_row_idx + 1:, col_start_idx:col_end_idx].reset_index(drop=True)
    df_extract.columns = df.iloc[header_row_idx, col_start_idx:col_end_idx].values  # Set header

    # Check for NaNs in the first column
    is_col1_na = df_extract.iloc[:,0].isna()

    if row_end_threshold > 1:
        # Looks for consecutive rows with NA
        # shift(-1) looks at the next row. fill_value=False ensures end of df doesn't trigger false positives.

        stop_mask = is_col1_na.copy()

        for i in range(1, row_end_threshold):
            stop_mask = stop_mask & is_col1_na.shift(-i, fill_value=False)

        if stop_mask.any():
            # The first True in stop_mask is the start of the gap
            stop_index = stop_mask.values.argmax()
            # Slice the DataFrame using .iloc up to the row immediately before the NA row (exclusive)
            df_extract = df_extract.iloc[:stop_index]

    else:
        if is_col1_na.any():
            # Find the positional index of the first NA value
            first_na_position = is_col1_na.values.argmax()
            # Slice the DataFrame using .iloc up to the row immediately before the NA row (exclusive)
            df_extract = df_extract.iloc[:first_na_position]

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

def extract_protocol_info(xls_obj: pd.ExcelFile, file_name: str):
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

def extract_metadata(pr_export_df):
    """
    Uses df from the 'Table All Cycles' sheet of the PR export and extracts metadata from the first column
    (measurement date, ID2: cell line, ID3: transfections #). Returns dictionary of extracted metadata.
    """
    col = pr_export_df.iloc[:, 0].astype(str)  # Get first col and transform everything to str

    # Mapping: { "Excel Label": "Desired Key" }
    meta_keys = {r"Date\s*:": "measurement_date",
                 r"ID2\s*:": "cell_line",
                 r"ID3\s*:": "transfections"}
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


def extract_bret_data(pr_export_df):
    """
    Uses df from the 'Table All Cycles' sheet of the PR export and extracts raw bret ratio table.
    Reformats df to header Time, wells.
    """
    df = slice_table(pr_export_df, "Well", row_end_threshold=2)
    # Clean up well names from "A01" to "A1"
    # "([A-Za-z])0(\d)" matching what to replace, () groups parts -> r"\1\2" replace with group1 and 2
    df.iloc[:, 0] = df.iloc[:, 0].str.replace(r"([A-Za-z])0(\d)", r"\1\2", regex=True)
    df.iat[0, 0] = "Time (min)"

    # --- Split in Lum count and BRET ratio ---
    cols_to_keep_lum = [c for c in df.columns if str(c).strip().startswith("Raw Data (475")]
    cols_to_keep_lum.append("Well")
    df_lum = df.loc[:, cols_to_keep_lum].copy()

    cols_to_drop_bret = [c for c in df.columns if str(c).strip().startswith("Raw")]
    cols_to_drop_bret.append("Content")
    df_bret = df.drop(columns=cols_to_drop_bret).copy()

    # Transpose
    df_bret = df_bret.set_index("Well").T
    df_bret = df_bret.reset_index(drop=True)  # Make sure index is clean

    df_lum = df_lum.set_index("Well").T
    df_lum = df_lum.reset_index(drop=True)

    return df_bret, df_lum


def extract_measurement_data(xls_obj, file_name: str):
    """
    Gets metadata and BRET ratio from PR export.
    Stores and returns PrResult class with all data.
    """
    worksheet = "Table All Cycles"
    try:
        # Read the first column of this sheet
        pr_export_df = pd.read_excel(xls_obj,
                                     sheet_name=worksheet,
                                     header=None
                                     )
    except ValueError:
        # Error if sheet is missing
        print(f"[ERROR] Worksheet '{worksheet}' not found in file.")
        return None

    metadata_dic = extract_metadata(pr_export_df)
    bret_ratio_df , lum_df= extract_bret_data(pr_export_df)

    result_obj = PrResult(
        file_name=file_name,
        measurement_date=metadata_dic['measurement_date'],
        cell_line=metadata_dic['cell_line'],
        transfection_id=metadata_dic['transfections'],
        raw_bret_ratio_df=bret_ratio_df,
        lum_df=lum_df
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

def process_bret_measurement(result: PrResult, protocol: ProtocolData, config: ProcessingConfig):
    """
    Maps cell line x transfection x ligand plate layout using protocol info.
    Performs baseline correction. Vehicle normalisation with kinetic data and AUC in parallel.
    Checks vehicle for outliers and luminescence count, saves warnings.
    """
    if result.is_excluded:
        result.kinetic_df = None
        result.kinetic_mean_df = None
        result.auc_df = None
        result.auc_mean_df = None
        result.vehicle_warnings = []
        result.low_lum_warnings = []
        result.vehicle_outliers = {}
        return result

    # Reset for re-run
    result.vehicle_warnings = []
    result.vehicle_outliers = {}
    result.low_lum_warnings = []
    result.low_lum_cond = {}
    result.column_metadata = {}

    # --- CONFIG PROCESSING ---
    acc_vehicle_range = config.vehicle_warning_threshold
    baseline_end_idx = config.baseline_end_index
    lum_threshold = config.lum_threshold
    date_str = result.measurement_date.strftime('%d.%m.%y')

    print(f"\n[DEBUG] === Processing File: {result.file_name} ===")

    # --- CONFIG LAYOUT ---
    # TODO: implement layout as part of ProcessingConfig to use for labeling experiments
    plate_blocks = config.plate_layout

    # --- METADATA MAPPING ---
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
        for rep_idx, col in enumerate(block_cols):
            current_rep_id = rep_idx + 1
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
                ligand_conc=current_conc_map,
                replicate=current_rep_id
            )

    # --- ASSIGNING EXCLUDED WELLS ---
    # Get raw BRET ratio table and lum table
    raw_df = result.raw_bret_ratio_df.copy()
    lum_df = result.lum_df.copy()

    if result.excluded_wells:
        # Set entire columns to NaN
        for well in result.excluded_wells:
            if well in raw_df.columns: raw_df[well] = float('nan')
            if well in lum_df.columns: lum_df[well] = float('nan')

    # --- LUM COUNT CHECK ---
    # Drop time col
    lum_calc_df = lum_df.drop(columns=["Time (min)"], errors='ignore').apply(pd.to_numeric, errors='coerce')
    # Use calculate_replicate_means to assign condition keys (and calc mean per row)
    lum_mean_df = calculate_means_on_meta(lum_calc_df, plate_blocks, result.column_metadata, grouping_mode="column")

    if not lum_mean_df.empty:
        # Check last 5 rows (time points)
        lum_end = lum_mean_df.iloc[-5:] if len(lum_mean_df) >= 5 else lum_mean_df
        mean_lum_end = lum_end.mean() # Mean value per condition key

        # Compare with threshold
        low_lum_cond = mean_lum_end[mean_lum_end < lum_threshold]

        for key, val in low_lum_cond.items():
            # Extract cond stats from key
            try:
                # "Cond_Name|Cell|Ligand"
                parts = key.split("|")
                if len(parts) < 4: continue # Safety check
                cond_name, cell_line, lig_name, repl_num = parts[0], parts[1], parts[2], parts[3]
                warning_dict = create_warning_record(
                    warn_type="Lum",
                    exp_date=date_str,
                    cond_name=cond_name,
                    cell_line=cell_line,
                    ligand=lig_name,
                    value=float(val),
                    replicate=str(repl_num)
                )
                if warning_dict not in result.low_lum_warnings:
                    result.low_lum_warnings.append(warning_dict)

            except Exception as e:
                print(f"[WARNING] Error parsing lum key {key}: {e}")

    # --- BUILT TIME VECTOR ---
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

    # --- RAW BRET LAST MEASUREMENT POINTS ---
    lp_raw_bret_df = raw_df.iloc[-3:].mean().to_frame().T

    # --- BASELINE CORRECTION ---
    baseline_means = data_df.iloc[0:baseline_end_idx].mean()
    bl_corrected_df = data_df / baseline_means

    # Handle columns where baseline_mean was 0 (to avoid infinity)
    bl_corrected_df = bl_corrected_df.replace([float('inf'), -float('inf')], float('nan'))

    # --- AUC CALCULATION ---
    # Use slicing to sum only the kinetic phase (after baseline)
    bl_corr_auc_df = bl_corrected_df.iloc[baseline_end_idx:].sum().to_frame().T

    # --- VEHICLE CORRECTION ---
    # Kinetics (df -> returns series of means over time)
    veh_means_kinetic = calculate_vehicle_means(
        bl_corrected_df, plate_blocks, result.excluded_wells
    )
    # For AUC (df with one row -> returns dictionary of scalars)
    veh_means_auc = calculate_vehicle_means(
        bl_corr_auc_df, plate_blocks, result.excluded_wells
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
            auc_norm_dict[col] = bl_corr_auc_df[col] / v_auc
        else:
            auc_norm_dict[col] = float('nan')

    # Create df from dict; index setting is required to handle excluded (nan) data
    kinetic_df = pd.DataFrame(kinetic_norm_dict, index=data_df.index)
    auc_df = pd.DataFrame(auc_norm_dict, index=[0])

    # --- VEHICLE CHECK ---
    print("[DEBUG] starting vehicle check")
    for block in plate_blocks:
        for c_idx in block:
            well_id = f"H{c_idx}"

            # Skip if well is already excluded or doesn't exist
            if well_id in result.excluded_wells or well_id not in kinetic_df.columns:
                continue

            # Check mean of entire kinetic against 1
            val_to_check = kinetic_df[well_id].mean()
            deviation = abs(val_to_check - 1)

            if deviation > acc_vehicle_range:
                meta = result.column_metadata.get(c_idx)

                warning_dict = create_warning_record(
                    warn_type="Veh",
                    exp_date=date_str,
                    cond_name=meta.condition_name,
                    cell_line=meta.cell_line,
                    ligand=meta.ligand_identity,
                    value=float(val_to_check),
                    replicate=str(meta.replicate),
                    row="H",  # Specific row
                    well_id=well_id
                )

                result.vehicle_warnings.append(warning_dict)

    # --- MEAN OF REPLICATES ---
    kinetic_mean_df = calculate_means_on_meta(
        kinetic_df, plate_blocks, result.column_metadata
    )
    auc_mean_df = calculate_means_on_meta(
        auc_df, plate_blocks, result.column_metadata
    )

    # --- SAVE RESULTS ---
    result.raw_bret_points_df = lp_raw_bret_df
    result.raw_bret_points_tidy_df = convert_to_plate_layout(lp_raw_bret_df)

    # --- Kinetic data
    result.time_vector = time_vec
    result.bl_corr_kinetic = bl_corrected_df
    result.kinetic_df = kinetic_df
    result.kinetic_mean_df = kinetic_mean_df

    # --- AUC data (tidy for CRC)
    result.bl_corr_auc_df = bl_corr_auc_df
    result.bl_corr_auc_tidy_df = convert_to_plate_layout(bl_corr_auc_df)
    result.auc_df = auc_df
    result.auc_tidy_df = convert_to_plate_layout(result.auc_df)
    result.auc_mean_df = auc_mean_df
    result.auc_mean_tidy_df = convert_to_plate_layout(result.auc_mean_df)

    return result

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

def create_clean_pivot(df, index_col, value_col, disregard_well_id):
    """
    Pivots the table. If Mean Data: One column per File.
    If Raw Data: One column per Well (Technical Replicates side-by-side)
    """
    is_kinetic = df["Time_(min)"].nunique() > 1

    # Assign a replicate number (1, 2, 3...) per Header_Key
    df = df.copy()
    if disregard_well_id:
        df['Rep_Num'] = df.groupby('Header_Key')['File_Name'].rank(method='dense').astype(int)
    else:
        # Consider File_Name and Well ID or Plate Col for technical replicates
        if is_kinetic:
            df['sort_key'] = df['File_Name'].astype(str) + "_" + df['Well_ID'].astype(str)
        else:
            # Only get Col number instead of complete well id
            df['sort_key'] = df['File_Name'].astype(str) + "_" + df['Well_ID'].astype(str).str[1:].str.zfill(2)
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

# --- Main Application --- #

class NCollectorApp:
    def __init__(self, main_window):
        self.main_gi = main_window
        main_window.title("N Collector")
        main_window.geometry("800x700")

        # --- Data Storage ---
        self.subfolder_paths_with_files = []
        self.experiment: list[MeasurementFolder] = []
        self.master_index = pd.DataFrame()  # Index for populating tab 2
        self.rule_history_text = ""
        self.pending_exclusions = []
        self.master_df = pd.DataFrame  # Used for master csv file storage (by generation or import)
        self.ignored_warnings = set()

        # --- GUI Variables ---
        self.folder_path = tk.StringVar(value="No folder selected.")
        self.var_lum_threshold = tk.IntVar(value=100)
        self.var_lig = tk.StringVar(value="")
        self.var_date = tk.StringVar(value="All")
        self.var_cell = tk.StringVar(value="All")
        self.var_cond = tk.StringVar(value="All")
        self.var_repl = tk.StringVar(value="All")
        self.var_row = tk.StringVar(value="All")
        self.var_data_type = tk.StringVar(value="")
        self.var_group_by = tk.StringVar(value="Transfection")

        # --- GUI Widgets (Initialised to None) ---
        self.log_window = None
        self.log_text = None
        self.notebook = None
        # Tab 1
        self.tab_import = None
        self.path_label = None
        self.load_files_button = None
        self.main_plasmids_label = None
        self.summary_tree = None
        self.lbl_rules_summary = None
        self.btn_export_master = None
        self.btn_export_excel = None
        # Tab 2
        self.tab_select = None
        self.cb_lig = None
        self.cb_date = None
        self.cb_cell = None
        self.cb_cond = None
        self.cb_rep = None
        self.lbl_row = None
        self.cb_row = None
        self.lb_exclusions = None
        # Tab 3
        self.tab_plot_helper = None
        self.lbl_data_source = None
        self.lb_ligands = None
        self.lb_exp_cells = None
        self.lb_exp_trans = None
        self.lb_kin_layout = None
        self.btn_run_plot_helper = None

        # --- Constants ---
        self.data_type_map = {
            "kinetic: raw BRET ratio (kinetic)": "Raw_BRET_kinetic",
            "kinetic: baseline-corrected BRET ratio": "Bl_Corrected_BRET",
            "kinetic: vehicle-normalised BRET ratio, techn. replicates": "Veh_Norm_Kinetic",
            "kinetic: vehicle-normalised BRET ratio, mean of techn. replicates": "Kinetic_Mean",

            "CRC: raw BRET (from last 3x time points)": "Raw_BRET_CRC",
            "CRC: baseline-corrected BRET ratio": "Bl_AUC",
            "CRC: vehicle-normalised BRET ratio, techn. replicates": "Veh_Norm_AUC",
            "CRC: vehicle-normalised BRET ratio, mean of techn. replicates": "AUC_Mean"
        }

        # --- Setup GUI ---
        self.setup_logging()
        self.setup_tabs()

    def setup_logging(self):
        """Setup of log window"""
        self.log_window = tk.Toplevel(self.main_gi)
        self.log_window.title("Processing Log")
        self.log_window.geometry("700x500")
        self.log_text = tk.Text(self.log_window)
        self.log_text.pack(expand=True, fill='both')
        # Opt. save of log file
        btn_frame = tk.Frame(self.log_window)
        btn_frame.pack(fill="x", padx=5, pady=5)
        tk.Button(btn_frame, text="Save Log to File", command=self.save_log_to_file).pack(side="right")
        tk.Button(btn_frame, text="Clear Log", command=lambda: self.log_text.delete('1.0', tk.END)).pack(side="left")

    def setup_tabs(self):
        """Setup tabs"""
        self.notebook = ttk.Notebook(self.main_gi)
        self.notebook.pack(expand=True, fill='both')

        self.tab_import = tk.Frame(self.notebook)
        self.notebook.add(self.tab_import, text="Import & Export Data")

        self.tab_select = tk.Frame(self.notebook)
        self.notebook.add(self.tab_select, text="Exclude Data")

        self.tab_plot_helper = tk.Frame(self.notebook)
        self.notebook.add(self.tab_plot_helper, text="Plot Helper")

        self.setup_import_tab()
        self.setup_exclusion_tab()
        self.setup_plot_helper_tab()

    def setup_import_tab(self):
        # Select Folder button
        tk.Button(self.tab_import, text="Select folder containing results of experiment",
                  command=self.select_folder).pack(pady=10, padx=10)

        # Display label for path
        self.path_label = tk.Label(self.tab_import, textvariable=self.folder_path, wraplength=700, justify="left",
                                   font=('Arial', 10))
        self.path_label.pack(pady=10, padx=10)

        # Load frame
        load_frame = tk.Frame(self.tab_import)
        load_frame.pack(pady=15, fill="x", padx=20)

        # Load button
        self.load_files_button = tk.Button(load_frame, text="Load Files", state="disabled",
                                           command=self.collect_files)
        self.load_files_button.pack(padx=10)

        # Threshold Input
        lum_thresh_entry = tk.Entry(load_frame, textvariable=self.var_lum_threshold, width=10)
        lum_thresh_entry.pack(side="right")
        tk.Label(load_frame, text="Lum. Threshold:").pack(side="right", padx=(10, 5))


        # Label to display Main Plasmids
        self.main_plasmids_label = tk.Label(self.tab_import, text="", justify="left", font=("Arial", 10, "bold"))
        self.main_plasmids_label.pack(pady=(0, 5))

        # Display of N summary table
        summary_frame = tk.Frame(self.tab_import)
        summary_frame.pack(pady=10, fill="both", expand=True, padx=20)
        tree_scroll = tk.Scrollbar(summary_frame) # Scrollbar for table
        tree_scroll.pack(side="right", fill="y")
        self.summary_tree = ttk.Treeview(summary_frame,
                                         columns=("Ligand", "Cell", "Cond", "N", "Dates"),
                                         show="headings",
                                         yscrollcommand=tree_scroll.set, height=6)
        tree_scroll.config(command=self.summary_tree.yview)
        # Define Columns
        self.summary_tree.heading("Ligand", text="Ligand")
        self.summary_tree.heading("Cell", text="Cell Line")
        self.summary_tree.heading("Cond", text="Condition")
        self.summary_tree.heading("N", text="N")
        self.summary_tree.heading("Dates", text="Dates")
        self.summary_tree.column("Ligand", width=80)
        self.summary_tree.column("Cell", width=100)
        self.summary_tree.column("Cond", width=250)
        self.summary_tree.column("N", width=30, anchor="center")
        self.summary_tree.column("Dates", width=150)
        self.summary_tree.pack(fill="both", expand=True)

        # Text of applies exclusion rules
        rules_frame = tk.LabelFrame(self.tab_import, text="Applied Exclusion Rules")
        rules_frame.pack(fill="x", padx=20, pady=5)
        self.lbl_rules_summary = tk.Label(rules_frame, text="No exclusion rules applied", justify="left", anchor="w")
        self.lbl_rules_summary.pack(fill="x", padx=5, pady=5)

        # --- EXPORT SECTION ---
        export_frame = tk.LabelFrame(self.tab_import, text="Export Options")
        export_frame.pack(fill="x", padx=20, pady=10)
        self.btn_export_master = tk.Button(export_frame, text="Export Master CSV", state="disabled",
                                           command=self.export_master_csv)
        self.btn_export_master.pack(side="left", fill="x", expand=True, padx=5, pady=10)
        self.btn_export_excel = tk.Button(export_frame, text="Export Excel Report (Default)", state="disabled",
                                          command=self.export_excel_report)
        self.btn_export_excel.pack(side="left", fill="x", expand=True, padx=5, pady=10)

    def update_summary_table(self):
        """Fills the summary table with N counts and dates per condition, per ligand"""
        # Clear existing data
        for i in self.summary_tree.get_children():
            self.summary_tree.delete(i)

        if self.master_index.empty:
            return

        # Group by Ligand, Cell Line and Condition
        grouped = self.master_index.groupby(['Ligand', 'Cell_Line', 'Condition'])

        for (lig, cell, cond), group in grouped:
            # Count unique filenames for N (experiments)
            n_count = group['File_Name'].nunique()

            # Get sorted unique dates
            unique_dates = sorted(group['Date'].unique())
            date_str = ", ".join(unique_dates)

            # Insert into tree
            self.summary_tree.insert("", "end", values=(lig, cell, cond, n_count, date_str))

    def setup_exclusion_tab(self):
        """Builds GUI for Tab 2 data selection"""

        # Frame for dropdowns
        filter_frame = tk.LabelFrame(self.tab_select, text="Exclude Data")
        filter_frame.pack(fill = "x", pady=10, padx=5)

        # Variables
        self.var_lig = tk.StringVar(value="")
        self.var_date = tk.StringVar(value="All")
        self.var_cell = tk.StringVar(value="All")
        self.var_cond = tk.StringVar(value="All")
        self.var_repl = tk.StringVar(value="")
        self.var_row = tk.StringVar(value="")
        # List to store rules
        self.pending_exclusions = []

        # 1. Ligand Dropdown
        tk.Label(filter_frame, text="Ligand:").grid(row=0, column=0, padx=5, pady=5)
        self.cb_lig = ttk.Combobox(filter_frame, textvariable=self.var_lig, state="readonly", width=12)
        self.cb_lig.grid(row=0, column=1, padx=5, pady=5)
        self.cb_lig.bind("<<ComboboxSelected>>", lambda e: self.update_dropdown_options("Ligand"))

        # 2. Date Dropdown
        tk.Label(filter_frame, text="Date:").grid(row=0, column=2, padx=5, pady=5)
        self.cb_date = ttk.Combobox(filter_frame, textvariable=self.var_date, state="readonly", width=8)
        self.cb_date.grid(row=0, column=3, padx=5, pady=5)
        self.cb_date.bind("<<ComboboxSelected>>", lambda e: self.update_dropdown_options("Date"))

        # 3. Cell Line Dropdown
        tk.Label(filter_frame, text="Cell Line:").grid(row=0, column=4, padx=5, pady=5)
        self.cb_cell = ttk.Combobox(filter_frame, textvariable=self.var_cell, state="readonly", width=8)
        self.cb_cell.grid(row=0, column=5, padx=5, pady=5)
        self.cb_cell.bind("<<ComboboxSelected>>", lambda e: self.update_dropdown_options("Cell_Line"))

        # 4. Condition Dropdown
        tk.Label(filter_frame, text="Condition:").grid(row=0, column=6, padx=5, pady=5)
        self.cb_cond = ttk.Combobox(filter_frame, textvariable=self.var_cond, state="readonly")
        self.cb_cond.grid(row=0, column=7, padx=5, pady=5)
        self.cb_cond.bind("<<ComboboxSelected>>", lambda e: self.update_dropdown_options("Condition"))

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

    def log(self, message):
        """Logs to the separate window"""
        try:
            self.log_text.insert(tk.END, message + "\n")
            self.log_text.see(tk.END)
            self.main_gi.update_idletasks()
        except tk.TclError:
            print(message)

    def show_warning_review(self, all_warnings):
        """Pop-up window allowing users to select warnings to exclude."""
        if not all_warnings:
            return

        dialog = tk.Toplevel(self.main_gi)
        dialog.title("Data Quality Warnings")
        dialog.geometry("700x500")

        tk.Label(dialog, text="Warning - Consider Excluding Items",
                 font=("Arial", 10, "bold")).pack(pady=10)
        tk.Label(dialog, text="Select to exclude:",
                 font=("Arial", 10)).pack(padx=15, pady=5, anchor="w")

        # Select/Deselect button
        vars_to_apply = []

        def select_all():
            for var, _ in vars_to_apply:
                var.set(True)

        def deselect_all():
            for var, _ in vars_to_apply:
                var.set(False)

        toggle_frame = tk.Frame(dialog)
        toggle_frame.pack(fill="x", padx=20, pady=2)

        tk.Button(toggle_frame, text="Select All",
                  command=select_all,
                  font=("Arial", 9)).pack(side="left", padx=(0, 10))

        tk.Button(toggle_frame, text="Deselect All",
                  command=deselect_all,
                  font=("Arial", 9)).pack(side="left")


        # Scrollable Frame Setup
        container = tk.Frame(dialog)
        container.pack(fill="both", expand=True, padx=10)
        canvas = tk.Canvas(container)
        scrollbar = ttk.Scrollbar(container, orient="vertical", command=canvas.yview)
        scrollable_frame = tk.Frame(canvas)

        scrollable_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        # List to track tuples of (BooleanVar, WarningData)
        vars_to_apply = []

        for warn in all_warnings:
            var = tk.BooleanVar(value=True)
            # Checkboxes
            cb = tk.Checkbutton(
                scrollable_frame,
                text=warn["Display"],
                variable=var,
                anchor="w"
            )
            cb.pack(fill="x", padx=10, pady=2, anchor="w")
            vars_to_apply.append((var, warn))

        # Action Buttons
        btn_frame = tk.Frame(dialog)
        btn_frame.pack(fill="x", pady=15)

        tk.Button(
            btn_frame,
            text="Exclude Selected & Re-Calculate",
            command=lambda: self.confirm_warning_selection(vars_to_apply, dialog),
            bg="#ffcccc",
            padx=10,
            relief="raised"
        ).pack(side="left", fill="x", expand=True, padx=5, pady=10)

        tk.Button(
            btn_frame,
            text="Cancel",
            command=dialog.destroy,
        ).pack(side="left", fill="x", expand=True, padx=5, pady=10)

    def confirm_warning_selection(self, vars_to_apply, dialog_window):
        """
        Processes selected check buttons and adds them to pending exclusions.
        Saves unchecked buttons as warnings to be ignored.
        """
        applied_any = False

        for var, data in vars_to_apply:
            if var.get():
                # Add to the pending exclusions list using existing format
                self.pending_exclusions.append({
                    "Ligand": data["Ligand"],
                    "Date": data["Date"],
                    "Cell_Line": data["Cell_Line"],
                    "Condition": data["Condition"],
                    "Row": data["Row"],
                    "Replicate": data["Replicate"]
                })

                # Update the GUI listbox in Tab 2
                display_str = f"AUTO: {data['Display']}"
                self.lb_exclusions.insert(tk.END, display_str)
                applied_any = True
            else:
                self.ignored_warnings.add(data['Display'])

        # Close the pop-up
        dialog_window.destroy()

        if applied_any:
            self.log("Auto-applying selected warnings...")
            # Directly trigger exclusion
            self.apply_exclusions()

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

        # Reset Variables
        self.var_lig.set("All")
        self.var_date.set("All")
        self.var_cell.set("All")
        self.var_cond.set("All")
        self.var_repl.set("")
        self.cb_row.config(state="disabled")

        # If only one ligand exists, default to it and disable the box.
        unique_ligands = sorted(self.master_index['Ligand'].dropna().unique().tolist())
        if len(unique_ligands) == 1:
            single_ligand = unique_ligands[0]
            self.var_lig.set(single_ligand)
            self.cb_lig.config(state="disabled")  # Lock user into this ligand
        else:
            self.var_lig.set("All")
            self.cb_lig.config(state="readonly")  # Allow selection

        # Trigger the update logic (trigger_source=None means full reset)
        trigger = "Ligand" if len(unique_ligands) == 1 else None
        self.update_dropdown_options(trigger_source=trigger)

    def update_dropdown_options(self, trigger_source=None):
        """
        Dynamically updates the values of all dropdowns based on the current selection of others.
        trigger_source: The name of the field that triggered the update.
        """
        if self.master_index.empty: return

        # Map Columns to their UI Components
        field_map = {
            "Ligand": (self.var_lig, self.cb_lig),
            "Date": (self.var_date, self.cb_date),
            "Cell_Line": (self.var_cell, self.cb_cell),
            "Condition": (self.var_cond, self.cb_cond)
        }

        # Get current selection
        current_selections = {col: var.get() for col, (var, _) in field_map.items()}

        # Iterate through each field
        for param, (target_var, target_widget) in field_map.items():
            # Built an all true mask
            mask = pd.Series(True, index=self.master_index.index)

            # Apply filters from ALL OTHER fields
            for curr_param, val in current_selections.items():
                # Skip the changed dropdown (curr_param) so it doesn't filter itself
                if curr_param != param and val != "All" and val != "":
                    mask &= (self.master_index[curr_param] == val)

            # Extract unique values using the mask directly
            valid_options = sorted(self.master_index.loc[mask, param].dropna().unique().tolist())

            # Update Widget
            target_widget['values'] = ["All"] + valid_options

            # If the current selection is no longer valid, reset it to "All"
            # Unless it was the user who just changed it (trigger_source)
            current_val = current_selections[param]
            if current_val != "All" and current_val not in valid_options:
                if param != trigger_source:
                    target_var.set("All")

        # Granular reset of Replicate specific to condition
        if trigger_source == "Condition":
            self.var_repl.set("")
            self.toggle_row_dropdown()

    def add_exclusion_rule(self):
        """
        Adds the current dropdown state to the pending list.
        Handles empty strings for Replicate/Col and Row.
        """
        rep_val = self.var_repl.get()
        row_val = self.var_row.get()

        rule = {
            "Ligand": self.var_lig.get(),
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
        rule_str = f"Ligand: {rule['Ligand']} | Date: {rule['Date']} | "\
                   f"Cell: {rule['Cell_Line']} | Cond: {rule['Condition']} | "\
                   f"Rep:{rep_str} | Row:{row_str}"

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
                    rule['Ligand'] == "All" and
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
            if rule.get('Ligand', 'All') != "All":
                df = df[df['Ligand'] == rule['Ligand']]
            if rule['Date'] != "All":
                df = df[df['Date'] == rule['Date']]
            if rule['Cell_Line'] != "All":
                df = df[df['Cell_Line'] == rule['Cell_Line']]
            if rule['Condition'] != "All":
                df = df[df['Condition'] == rule['Condition']]
            if rule['Replicate'] != "":
                target_rep = int(rule['Replicate'])
                df = df[df['Replicate'] == target_rep]

            if df.empty:
                self.log(f"   [WARNING] Rule {rule} matched 0 records.")
                continue

            # Handle specific Replicates and Rows
            for index, row_data in df.iterrows():
                result_obj = row_data['Ref_Result']
                col_idx = int(row_data['Column_Index'])

                target_rows = "ABCDEFGH"
                if rule['Row'] != "": target_rows = rule['Row']

                for r in target_rows:
                    well_id = f"{r}{col_idx}"
                    if well_id not in result_obj.excluded_wells:
                        result_obj.excluded_wells.append(well_id)
                        count_wells += 1
                        print(f"Excluded {well_id} in {result_obj.file_name}")

        if count_files > 0: self.log(f"   [DONE] Excluded {count_files} entire files.")
        if count_wells > 0: self.log(f"   [DONE] Excluded {count_wells} specific wells.")

        # Clear list after applying
        self.clear_exclusion_list()

        # Rerun processing to update graphs/stats
        self.run_processing_pipeline()
        self.refresh_filter_options()

    def setup_plot_helper_tab(self):
        """Builds the GUI for tab 3 plot helper"""

        # --- Import Master CSV ---
        src_frame = tk.Frame(self.tab_plot_helper)
        src_frame.pack(fill="x", padx=10, pady=10)
        tk.Label(src_frame, text="Data Source", font=("Arial", 9, "bold")).pack(side="top")

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

        # All list boxes expand vertically
        sel_frame.rowconfigure(1, weight=1)

        # Select All Buttons
        tk.Button(sel_frame, text="Select All Ligands", command=lambda: self.lb_ligands.select_set(0, tk.END)).grid(
            row=2, column=0, pady=8)
        tk.Button(sel_frame, text="Select All Cells", command=lambda: self.lb_exp_cells.select_set(0, tk.END)).grid(
            row=2, column=1, pady=8)
        tk.Button(sel_frame, text="Select All Transf.", command=lambda: self.lb_exp_trans.select_set(0, tk.END)).grid(
            row=2, column=2, pady=8)

        sel_frame.columnconfigure(0, weight=1)
        sel_frame.columnconfigure(1, weight=1)
        sel_frame.columnconfigure(2, weight=3)

        # --- Data Type Selection (Single Choice) ---
        type_frame = tk.LabelFrame(self.tab_plot_helper, text="Export Settings")
        type_frame.pack(fill="both", expand=True, padx=10, pady=5)

        # Data Type Dropdown
        tk.Label(type_frame, text="Data Type:").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        default_key = list(self.data_type_map.keys())[6]
        self.var_data_type = tk.StringVar(value=default_key)
        self.var_data_type.trace_add("write", self.toggle_kinetic_options) # Trace kinetic selection

        combo_type = ttk.Combobox(type_frame, textvariable=self.var_data_type, state="readonly", width=57)
        combo_type['values'] = list(self.data_type_map.keys())
        combo_type.grid(row=0, column=1, padx=5, pady=5)

        # --- Layout Selection ---
        # Groupy py
        tk.Label(type_frame, text="Group By:").grid(row=1, column=0, padx=5, pady=5, sticky="w")
        self.var_group_by = tk.StringVar(value="None")
        combo_group = ttk.Combobox(type_frame, textvariable=self.var_group_by, state="readonly", width=15)
        combo_group['values'] = ["None", "Cell Line", "Transfection"]
        combo_group.grid(row=1, column=1, padx=5, pady=5, sticky="w")

        # Kinetic Layout (Only applies if a Kinetic type is chosen)
        tk.Label(type_frame, text="Kinetic Layout:").grid(row=0, column=2, padx=5, pady=5, sticky="e")
        self.lb_kin_layout = tk.Listbox(type_frame, selectmode="multiple", height=6, exportselection=False)
        self.lb_kin_layout.grid(row=0, column=3, rowspan = 2, padx=3, pady=5, sticky="nsew")
        tk.Button(type_frame, text="Select All Rows", command=lambda: self.lb_kin_layout.select_set(0, tk.END)).grid(
            row=1, column=2, pady=5, sticky="nsew")

        type_frame.columnconfigure(1, weight=1)
        type_frame.columnconfigure(3, weight=3)
        type_frame.rowconfigure(0, weight=1)

        # Initialize state based on default value
        self.toggle_kinetic_options()

        # --- Action Button ---
        btn_frame = tk.Frame(self.tab_plot_helper)
        btn_frame.pack(fill="x", padx=10, pady=20)

        self.btn_run_plot_helper = tk.Button(btn_frame, text="Generate Custom Export",
                                            state="disabled", command=self.run_plot_helper)
        self.btn_run_plot_helper.pack(fill="x", ipady=5)

    def toggle_kinetic_options(self, *args):
        """Enables/Disables Kinetic Layout dropdown based on Data Type selection."""
        selection = self.var_data_type.get()

        # Check if "kinetic" is in the selected string
        if "kinetic" in selection.lower():
            self.lb_kin_layout.config(state="normal")
        else:
            self.lb_kin_layout.config(state="disabled")

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
        # Enable kinetic layout to populate list box, afterwards disable
        self.lb_kin_layout.config(state="normal")

        # Clear
        self.lb_ligands.delete(0, tk.END)
        self.lb_exp_cells.delete(0, tk.END)
        self.lb_exp_trans.delete(0, tk.END)
        self.lb_kin_layout.delete(0, tk.END)

        df = self.master_df

        # Populate Ligand
        ligands = sorted(set(df['Ligand'].astype(str)))
        for l in ligands: self.lb_ligands.insert(tk.END, l)
        self.lb_ligands.select_set(0, tk.END) # Default to all

        # Populate Cells
        cells = sorted(set(df['Cell_Line'].astype(str)))
        for c in cells: self.lb_exp_cells.insert(tk.END, c)
        self.lb_exp_cells.select_set(0, tk.END)  # Default to all

        # Populate Transfections
        trans = sorted(set(df['Transfection'].astype(str)))
        for t in trans: self.lb_exp_trans.insert(tk.END, t)
        self.lb_exp_trans.select_set(0, tk.END)  # Default to all

        # Populate Kinetic Layout
        rows = build_row_info_str(df)
        for r in rows: self.lb_kin_layout.insert(tk.END, r)
        self.lb_kin_layout.select_set(0, tk.END)  # Default to all
        # Disable listbox if not kinetic is selected
        self.toggle_kinetic_options()

    def run_plot_helper(self):
        """Collects GUI selections and calls the core export engine."""
        if self.master_df is None or self.master_df.empty: return

        # Get Selections
        cells = [self.lb_exp_cells.get(i) for i in self.lb_exp_cells.curselection()]
        transfections = [self.lb_exp_trans.get(i) for i in self.lb_exp_trans.curselection()]
        ligands = [self.lb_ligands.get(i) for i in self.lb_ligands.curselection()]
        rows = [self.lb_kin_layout.get(i) for i in self.lb_kin_layout.curselection()]

        display_name = self.var_data_type.get()
        if not display_name: return

        # TRANSLATE: Display Name -> Internal Column Name
        internal_name = self.data_type_map.get(display_name)

        if not internal_name:
            print(f"[ERROR] Unknown data type selected: {display_name}")
            return

        config = {
            'cells': cells,
            'transfections': transfections,
            'ligands': ligands,
            'data_types': [internal_name], # As list for engine compatibility with default export
            'group_by': self.var_group_by.get(),
            'kinetic_mode': rows
        }

        file_path = filedialog.asksaveasfilename(
            defaultextension=".xlsx",
            filetypes=[("Excel", "*.xlsx")],
            title="Save Custom Export"
        )
        if file_path:
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
                        "Ligand": meta.ligand_identity,
                        "Transfection_ID": meta.transfection_id,
                        "Column_Index": col_idx,
                        "Replicate": meta.replicate,
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
        self.ignored_warnings.clear()

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
                        protocol_info = extract_protocol_info(xls, file_name)

                        if protocol_info and protocol_info.exp_date == folder_date_obj:
                            folder_data.protocol = protocol_info
                            self.log(f"   [PROTOCOL] loaded: {file_name}")
                            is_imported = True
                        elif protocol_info:
                            self.log(f"   [MISMATCH] Protocol {protocol_info.exp_date} != Folder {folder_date_obj}")

                    # Identify PR export
                    # Must have "Table All Cycles", and the only allowed other sheet is "Protocol Information"
                    elif "Table All Cycles" in sheet_names and set(sheet_names).issubset({"Table All Cycles", "Protocol Information"}):
                        meas_data = extract_measurement_data(xls, file_name)

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

        # Get config from GUI state
        try:
            lum_threshold = self.var_lum_threshold.get()
        except tk.TclError:
            lum_threshold = 100

        current_config = ProcessingConfig(
            lum_threshold=lum_threshold
        )

        # Check for plasmids transfected in all conditions (main plasmids) and filter if needed
        selected_exp_name = self.handle_main_plasmids_selection()
        self.main_plasmids_label.config(text=f"{selected_exp_name}")

        # List to collect all warnings
        all_detected_warnings = []

        # Iterate through data and perform mapping+calculations
        for folder in self.experiment:
            if not folder.protocol: continue # Protocol is needed for processing
            for result in folder.results:
                # Process each result file within one folder (belonging to one protocol)
                # Also assigns conditions to data
                result = process_bret_measurement(result, folder.protocol, current_config)

                # Gather warnings
                all_detected_warnings.extend(result.low_lum_warnings)
                all_detected_warnings.extend(result.vehicle_warnings)

                # Handle outliers stored in dic
                if result.vehicle_outliers:
                     vehicle_out = ", ".join([f"{well} = {val:.2f}" for well, val in result.vehicle_outliers.items()])
                     self.log(f"   [VEHICLE WARNING] {result.file_name}: {vehicle_out}")

        # Built master indexing table (needed for flexible data exclusion)
        self.built_master_index()

        # Only show warnings that were not ignored previously
        new_warnings = [
            w for w in all_detected_warnings
            if w['Display'] not in self.ignored_warnings
        ]

        if new_warnings:
            self.show_warning_review(new_warnings)
        elif all_detected_warnings:
            self.log(f"[INFO] {len(all_detected_warnings)} warnings detected but previously ignored.")

        self.log("\n--- Compiling Master Dataframe... ---")
        self.master_df = self.compile_master_dataframe()
        self.refresh_plot_helper_options()

        self.log("\n--- Processing Complete & Plot Helper Ready ---")

        # Enable Exports
        self.btn_export_master.config(state="normal")
        self.btn_export_excel.config(state="normal")

    def compile_master_dataframe(self):
        """
        Compiles all processing steps into one Master DataFrame.
        Structure: 1 row per well per timepoint.
        Means are repeated for respective technical replicates as AUCs for all timepoints.
        """
        if not self.experiment:
            return None
        self.log("\n--- Building Master CSV ---")

        all_files_data = []

        for folder in self.experiment:
            if folder.protocol.main_plasmids:
                main_plasmids = " + ".join(folder.protocol.main_plasmids)
            else:
                main_plasmids = "Unknown"

            for res in folder.results:
                if res.is_excluded: continue
                if res.kinetic_df is None or res.raw_bret_ratio_df is None:
                    continue # ----------------------add all that is needed!

                # --- PREPARE KINETIC DATA ---
                # Use pandas melt function to prepare each df from wide to long format
                def melt_df(df, val_name, time_vec):
                    if df is None or df.empty: return pd.DataFrame()

                    df_work = df.copy()
                    # Check lengths
                    if len(df_work) != len(time_vec):
                        print(
                            f"[WARNING] Length mismatch in {res.file_name}: Data {len(df_work)} vs Time {len(time_vec)}")
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
                    t_vec = range(len(res.raw_bret_ratio_df))

                # Ignore time col in raw bret df
                raw_clean = res.raw_bret_ratio_df.drop(columns=["Time (min)"], errors='ignore')

                df_raw = melt_df(raw_clean, "Raw_BRET_kinetic", t_vec)
                df_bl = melt_df(res.bl_corr_kinetic, "Bl_Corrected_BRET", t_vec)
                df_norm = melt_df(res.kinetic_df, "Veh_Norm_Kinetic", t_vec)

                # Merge on [Time_(min), Well_ID]
                merge_on = [df_raw.columns[0], "Well_ID"]

                merged_df = df_raw.merge(df_bl, on=merge_on, how="left") \
                    .merge(df_norm, on=merge_on, how="left")

                # --- MAP AUC DATA and RAW BRET POINTS---
                # AUC is 1 value per well. We map it to Well_ID.
                raw_bret_map = res.raw_bret_points_df.iloc[0].to_dict() if res.raw_bret_points_df is not None else {}
                auc_bl_map = res.bl_corr_auc_df.iloc[0].to_dict() if res.bl_corr_auc_df is not None else {}
                auc_norm_map = res.auc_df.iloc[0].to_dict() if res.auc_df is not None else {}

                merged_df['Raw_BRET_CRC'] = merged_df['Well_ID'].map(raw_bret_map)
                merged_df['Bl_AUC'] = merged_df['Well_ID'].map(auc_bl_map)
                merged_df['Veh_Norm_AUC'] = merged_df['Well_ID'].map(auc_norm_map)

                # --- PREPARE MEAN KINETIC AND AUC MAPPING ---
                well_to_mean_map = {}
                well_to_auc_mean_map = {}
                for col_idx, meta in res.column_metadata.items():
                    col_str = str(col_idx)
                    for row_char in "ABCDEFGH":
                        well_id = f"{row_char}{col_idx}"

                        # Construct Key: "Condition|Cell|Ligand|Row"
                        mean_key = f"{meta.condition_name}|{meta.cell_line}|{meta.ligand_identity}|{row_char}"

                        # Grab Kinetic Mean Series
                        if res.kinetic_mean_df is not None and mean_key in res.kinetic_mean_df.columns:
                            well_to_mean_map[well_id] = res.kinetic_mean_df[mean_key].tolist()

                        # Grab AUC Mean Value
                        if res.auc_mean_df is not None and mean_key in res.auc_mean_df.columns:
                            well_to_auc_mean_map[well_id] = res.auc_mean_df[mean_key].iloc[0]

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
                merged_df["File_Name"] = res.file_name
                merged_df["Date"] = res.measurement_date
                merged_df["Main_Plasmids"] = main_plasmids

                # Get the exclusion text (handle empty case)
                exclusion_text = self.rule_history_text if self.rule_history_text else "None"
                # Clean newlines for CSV compatibility
                exclusion_text_clean = exclusion_text.replace("\n", " | ")
                merged_df["Applied_Exclusions"] = exclusion_text_clean

                # Meta Lookups (Optimization: Build dicts once per file)
                meta_lookups = {'Transfection': {}, 'Cell_Line': {}, 'Ligand': {}, 'Ligand_Conc': {}, 'Plate_Row': {}}

                for well_id in merged_df['Well_ID'].unique():
                    try:
                        c_idx = int(well_id[1:])
                        r_char = well_id[0]
                        meta = res.column_metadata.get(c_idx)
                        if meta:
                            meta_lookups['Transfection'][well_id] = meta.condition_name
                            meta_lookups['Cell_Line'][well_id] = meta.cell_line
                            meta_lookups['Ligand'][well_id] = meta.ligand_identity
                            meta_lookups['Ligand_Conc'][well_id] = meta.ligand_conc.get(r_char, 0.0)
                            meta_lookups['Plate_Row'][well_id] = r_char
                    except: pass

                merged_df['Transfection'] = merged_df['Well_ID'].map(meta_lookups['Transfection'])
                merged_df['Cell_Line'] = merged_df['Well_ID'].map(meta_lookups['Cell_Line'])
                merged_df['Ligand'] = merged_df['Well_ID'].map(meta_lookups['Ligand'])
                merged_df['Ligand_Conc'] = merged_df['Well_ID'].map(meta_lookups['Ligand_Conc'])
                merged_df['Plate_Row'] = merged_df['Well_ID'].map(meta_lookups['Plate_Row'])

        if not all_files_data:
            return pd.DataFrame()

        # Combine all files
        master_df = pd.concat(all_files_data, ignore_index=True)
        # Cleanup columns
        cols_order = [
            "File_Name", "Date", "Main_Plasmids", "Applied_Exclusions",
            "Transfection", "Cell_Line", "Ligand",
            "Ligand_Conc", "Plate_Row", "Well_ID", "Time_(min)",
            "Raw_BRET_kinetic", "Bl_Corrected_BRET", "Veh_Norm_Kinetic", "Kinetic_Mean",
            "Raw_BRET_CRC", "Bl_AUC", "Veh_Norm_AUC", "AUC_Mean"
        ]
        final_cols = [c for c in cols_order if c in master_df.columns]
        return master_df[final_cols]

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
        Writes the Excel file based on the config dictionary provided by either tab 1 (default) or tab 3 (user).
        """
        try:
            df_subset = apply_export_filters(self.master_df, config)

            if df_subset.empty:
                print("[ERROR] Export failed: Filter resulted in no data.")
                return

            # Define groups
            group_by = config.get('group_by', 'None')
            data_groups = []  # List of tuples: (Group_Name, DataFrame)

            if group_by == 'Cell Line':
                for name, group in df_subset.groupby('Cell_Line'):
                    data_groups.append((str(name), group))
            elif group_by == 'Transfection':
                for name, group in df_subset.groupby('Transfection'):
                    data_groups.append((str(name), group))
            else:
                data_groups.append(("", df_subset))  # No grouping

            # Definition which is kinetic and what is CRC
            kinetic_types = [val for key, val in self.data_type_map.items() if "kinetic:" in key]
            crc = [val for key, val in self.data_type_map.items() if "CRC:" in key]

            with pd.ExcelWriter(file_path) as writer:
                sheets_written = False

                # --- 1. METADATA SHEET ---
                file_names = df_subset["File_Name"].unique().tolist()
                mp_str = "Unknown"
                if 'Main_Plasmids' in df_subset.columns:
                    mp_vals = df_subset['Main_Plasmids'].unique()
                    if len(mp_vals) > 0: mp_str = mp_vals[0]
                if 'Ligand' in df_subset.columns:
                    ligand = df_subset['Ligand'].unique()
                print(ligand)

                meta_dict = {
                    "Export Date": [datetime.now().strftime("%d.%m.%Y - %H:%M:%S")],
                    "Main Plasmids": [mp_str],
                    "Ligand": [", ".join(ligand)],
                    "Source Files Count": [len(file_names)],
                    "Source Files List": [", ".join(file_names)],
                    "Data Type": [", ".join(config.get('data_types', []))],
                    "Kinetic Layout": [config.get('kinetic_mode')],
                    "Filter: Ligands": [", ".join(config.get('ligands'))],
                    "Filter: Cells": [", ".join(config.get('cells'))],
                    "Filter: Conditions": [", ".join(config.get('transfections'))],
                    "Group By": [group_by]
                }
                pd.DataFrame(meta_dict).transpose().to_excel(writer, sheet_name="Metadata", header=False)
                sheets_written = True

                # Iterate Groups + data types
                for group_name, df_group in data_groups:

                    # --- Loop through selected data types
                    # This handles both Single Selection (Plot Helper) and Default Report (List of 2)
                    selected_types = config.get('data_types', [])

                    # --- KINETIC DATA ---
                    for dtype in selected_types:
                        if dtype in kinetic_types:
                            k_layout = config.get('kinetic_mode', [])
                            if not k_layout:
                                continue

                            filtered_rows = []
                            for row_info in k_layout:
                                row_letter = row_info.split(":")[0].replace("Row ", "").strip()
                                ligand_conc = row_info.split(":")[1].split("log(M)")[0].strip()
                                ligand_name = row_info.split("log(M)")[1].strip()

                                # Filter row
                                row_df = df_group[
                                    (df_group["Plate_Row"] == row_letter) &
                                    (df_group["Ligand_Conc"].astype(str) == ligand_conc) &
                                    (df_group["Ligand"] == ligand_name)
                                    ].copy()

                                filtered_rows.append(row_df)

                            if filtered_rows:
                                df_kin = pd.concat(filtered_rows).drop_duplicates()
                            else:
                                continue

                            if df_kin.empty:
                                continue

                            df_kin = generate_header_key(df_kin, group_by)
                            kin_pivot = create_clean_pivot(df_kin, "Time_(min)", dtype, "Mean" in dtype)
                            kin_pivot.rename(columns={"Time_(min)": "Time (min)"}, inplace=True)

                            # Sheet Name with group_prefix (Max 31 chars)
                            base = f"{group_name}_{dtype}" if group_name else f"{dtype}"
                            sheet_name = base[:31]
                            # Save
                            kin_pivot.to_excel(writer, sheet_name=sheet_name, index=True)
                            sheets_written = True

                    # --- CRC DATA ---
                        elif dtype in crc:
                            df_crc = df_group.drop_duplicates(
                                subset=["File_Name", "Transfection", "Cell_Line", "Well_ID", "Ligand"]).copy()
                            if df_crc.empty: continue
                            df_crc = generate_header_key(df_crc, group_by)

                            if dtype not in df_crc.columns: continue

                            # Always pivot on plate row
                            crc_pivot = create_clean_pivot(df_crc, "Plate_Row", dtype, "Mean" in dtype)

                            # Display ligand conc instead of plate row if there is one ligand
                            if df_crc['Ligand'].nunique() == 1:
                                row_map = df_crc.drop_duplicates("Plate_Row").set_index("Plate_Row")["Ligand_Conc"]
                                crc_pivot.index = crc_pivot.index.map(row_map)
                                crc_pivot.rename(columns={"Plate_Row": "Concentration (logM)"}, inplace=True)

                            else:
                                crc_pivot.index.name = "Plate Row"

                            base = f"{group_name}_AUC" if group_name else f"AUC_{dtype}"
                            sheet_name = base[:31]
                            crc_pivot.to_excel(writer, sheet_name=sheet_name, index=True)
                            sheets_written = True
            if not sheets_written:
                pd.DataFrame({"Info": ["No data"]}).to_excel(writer, sheet_name="Empty")

            self.log(f"   [SUCCESS] Exported: {os.path.basename(file_path)}")

        except Exception as e:
            self.log(f"   [ERROR] Export failed: {e}")

    def export_excel_report(self):
        """
        Default Export: All Data, Highest stim kinetics and AUC crc.
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

        # Built kinetic_mode selection
        df_row_a = self.master_df[self.master_df['Plate_Row'] == 'A']
        kinetic_rows = build_row_info_str(df_row_a)

        # Define Standard Config
        default_config = {
            'cells': 'All',
            'transfections': 'All',
            'ligands': 'All',
            'data_types': ['Kinetic_Mean', 'AUC_Mean'],
            'group_by': 'Transfection',
            'kinetic_mode': kinetic_rows
        }
        self.write_excel_export(file_path, self.master_df, default_config)

    def save_log_to_file(self):
        """Exports the current log to a text file."""
        # Get content from line 1, char 0 to End
        log_content = self.log_text.get("1.0", tk.END)

        if not log_content.strip():
            print("Log is empty, nothing to save.")
            return

        file_path = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text File", "*.txt"), ("All Files", "*.*")],
            title="Save Log File"
        )

        if file_path:
            try:
                with open(file_path, "w", encoding="utf-8") as f:
                    f.write(log_content)
                print(f"Log saved to: {file_path}")
            except Exception as e:
                print(f"Error saving log: {e}")

# --- Main Execution Block ---
if __name__ == "__main__":
    # Create the main window
    root = tk.Tk()

    # Create an instance of the application
    app = NCollectorApp(root)

    # Start the Tkinter event loop
    root.mainloop()