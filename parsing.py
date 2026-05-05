import os
import logging
import pandas as pd
import re
from datetime import datetime

from models import ProtocolData, PrResult, MeasurementFolder

logger = logging.getLogger("NCollector")

# --- Spreadsheet Navigation Helpers --- #

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
        logger.error(f"Marker '{marker}' not found in sheet.")
        return None

    # Check if the requested occurrence exists
    if match_index >= len(all_matches):
        logger.error(
            f"Requested occurrence #{match_index + 1} of '{marker}' not found. Only {len(all_matches)} found.")
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
        logger.error(f"Target for '{marker}' is outside the sheet boundaries.")
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

# --- Protocol & Measurement Extraction --- #

def process_transfection_scheme(df: pd.DataFrame):
    """
    Gets the BRET pair (returns as list) and the experimental conditions (returns as dic) from the transfection table.
    """
    if df is None or df.empty:
        return [], {}

    # Separate metadata cols from transfection cols
    metadata_cols = {'dna', 'db#', 'db#/ flash', 'conc (ng/ul)', 'vol per transfection', 'vol master'}
    transfection_cols = [
        c for c in df.columns
        if str(c).lower() not in metadata_cols
    ]

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
        logger.error(f"Worksheet '{protocol_worksheet}' not found in file.")
        return None

    exp_date = None
    date_str = extract_value(protocol_sheet, "date of measurement")
    # Transform date str in real
    if date_str is not None:
        try:
            exp_date = datetime.strptime(date_str.strip(), '%d.%m.%y').date()
        except ValueError:
            logger.warning(f"Protocol date '{date_str}' not in DD.MM.YY format.")

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

    ligand_1 = extract_value(protocol_sheet, "Ligand dilution", col_offset=0, row_offset=1).strip()
    ligand_1_conc = slice_table(protocol_sheet, "final concentration in well (log(M))")

    # Check whether second ligand was selected
    check_ligand_2 = extract_value(protocol_sheet, "Ligand dilution", col_offset=0, row_offset=1, match_index=1)
    if check_ligand_2 is not None and str(check_ligand_2).strip().lower() not in ["nan", ""]:
        ligand_2 = check_ligand_2.strip()
        ligand_2_conc = slice_table(protocol_sheet, "final concentration in well (log(M))", match_index=1)
        ligand_layout = extract_value(protocol_sheet,"Ligand layout", col_offset= 0, row_offset= 1)
    else:
        ligand_2 = None
        ligand_2_conc = None
        ligand_layout = None

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
                                 ligand_2_conc= ligand_2_conc,
                                 ligand_layout = ligand_layout)
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

        # Re-join unique parts (e.g. "Control, dQ"); dict.fromkeys removed duplicates and keeps order
        metadata['cell_line'] = ", ".join(dict.fromkeys(standardised_lines))

    # Transform date str to actual date
    if 'measurement_date' in metadata:
        date_str = metadata["measurement_date"]
        try:
            metadata["measurement_date"] = datetime.strptime(date_str.strip(), '%d/%m/%Y').date()
        except ValueError:
            logger.warning(f"Analysis date '{date_str}' not in DD/MM/YYYY format.")
            return None  # Fail extraction if date is invalid

    # Ensure all required metadata fields were found
    if 'measurement_date' not in metadata or 'cell_line' not in metadata or 'transfections' not in metadata:
        logger.warning("Missing Date, ID2, or ID3 from metadata sheet.")
        return None

    return metadata


def extract_bret_data(pr_export_df):
    """
    Extracts BRET ratio and both raw data channels from the 'Table All Cycles' sheet.
    Auto-detects wavelengths from column headers (e.g. 'Raw Data (475-30 B)').
    Assigns lower wavelength as donor, higher as acceptor.
    Returns dict with keys: bret_ratio, donor, acceptor, donor_wavelength, acceptor_wavelength.
    """

    # Extract entire table from plate reader export
    df = slice_table(pr_export_df, "Well", row_end_threshold=2)
    # Clean up well names from "A01" to "A1"
    df.iloc[:, 0] = df.iloc[:, 0].str.replace(r"([A-Za-z])0(\d)", r"\1\2", regex=True)
    df.iat[0, 0] = "Time (min)"

    time_values = df.iloc[0, 2:].unique().astype(float)

    # --- Detect Raw Data column groups and extract wavelengths ---
    # Use positional indices to handle duplicate column names correctly
    wl_col_positions = {}  # {wavelength_int: [positional_indices]}
    ratio_positions = []
    well_pos = None # "Well" columns as positional anchor
    content_pos = None

    for i, col_name in enumerate(df.columns):
        c_str = str(col_name).strip()
        if c_str == "Well":
            well_pos = i
        elif c_str == "Content":
            content_pos = i
        elif c_str.startswith("Raw Data ("):
            match = re.search(r"Raw Data \((\d+)", c_str)
            if match:
                wl = int(match.group(1))
                wl_col_positions.setdefault(wl, []).append(i)
        elif c_str.startswith("Ratio"):
            ratio_positions.append(i)

    # Sort wavelengths: lower = donor, higher = acceptor
    sorted_wls = sorted(wl_col_positions.keys())
    if len(sorted_wls) < 2:
        logger.warning(f"Expected 2 raw data channels, found {len(sorted_wls)}: {sorted_wls}")

    donor_wl = sorted_wls[0] if len(sorted_wls) >= 1 else 0
    acceptor_wl = sorted_wls[1] if len(sorted_wls) >= 2 else 0

    # --- Extract each channel by positional index, transpose to Time x Wells ---
    def extract_channel_by_pos(positions):
        """Selects columns by position, transposes to Time x Wells format."""
        col_indices = [well_pos] + positions
        df_ch = df.iloc[:, col_indices].copy()
        df_ch = df_ch.set_index(df_ch.columns[0]).T.reset_index(drop=True)
        return df_ch

    df_donor = extract_channel_by_pos(wl_col_positions[donor_wl]) if donor_wl else pd.DataFrame()
    df_acceptor = extract_channel_by_pos(wl_col_positions[acceptor_wl]) if acceptor_wl else pd.DataFrame()

    # --- Extract BRET ratio by positional index ---
    ratio_col_indices = [well_pos] + ratio_positions
    df_bret = df.iloc[:, ratio_col_indices].copy()
    df_bret = df_bret.set_index(df_bret.columns[0]).T.reset_index(drop=True)

    return {
        "time_min": time_values,
        "bret_ratio": df_bret,
        "donor": df_donor,
        "acceptor": df_acceptor,
        "donor_wavelength": donor_wl,
        "acceptor_wavelength": acceptor_wl,
    }


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
        logger.error(f"Worksheet '{worksheet}' not found in file.")
        return None

    metadata_dic = extract_metadata(pr_export_df)
    bret_data = extract_bret_data(pr_export_df)

    result_obj = PrResult(
        file_name=file_name,
        measurement_date=metadata_dic['measurement_date'],
        cell_line=metadata_dic['cell_line'],
        transfection_id=metadata_dic['transfections'],
        raw_time=bret_data["time_min"],
        raw_bret_ratio_df=bret_data["bret_ratio"],
        donor_df=bret_data["donor"],
        acceptor_df=bret_data["acceptor"],
        donor_wavelength=bret_data["donor_wavelength"],
        acceptor_wavelength=bret_data["acceptor_wavelength"],
    )
    return result_obj


def extract_info_sheet_data(xls_obj) -> dict:
    """
    Extracts key-value pairs from the 'Protocol Information' sheet in a plate reader xlsx.
    Returns a dict of extracted info, or empty dict if the sheet is missing or unreadable.
    """
    worksheet = "Protocol Information"
    try:
        # Read the first two columns of this sheet
        info_sheet_df = pd.read_excel(xls_obj,
                                     sheet_name=worksheet,
                                     header=None)
    except ValueError:
        return {}

    # Extract any info as dict, value is either in col a or b
    info_dict = {}
    for index, row in info_sheet_df.iterrows():
        col_a = row[0]
        col_b = row[1] if len(row) > 1 else None

        if pd.isna(col_a):
            continue # Skip rows where the first column is empty

        # Convert to str only for splitting to ensure NAs are captured correctly
        col_a_str = str(col_a)

        # Check if the string contains our key-value separator
        if ":" in col_a_str:
            # Split by the FIRST colon only (important for Path key)
            key, val_in_col1 = col_a_str.split(":", 1)
            key = key.strip()
            val_in_col1 = val_in_col1.strip()

            if val_in_col1:
                info_dict[key] = val_in_col1

            elif pd.notna(col_b):
                info_dict[key] = col_b.strip() if isinstance(col_b, str) else col_b

    return info_dict


# --- Folder Scanning --- #

def scan_and_load_folders(folder_paths: list[str], log_fn=None) -> list[MeasurementFolder]:
    """
    Scans folder paths for xlsx/xlsm files, classifies them as protocol or measurement,
    and returns a list of MeasurementFolder objects with matched protocol + results.

    Args:
        folder_paths: List of directory paths to scan.
        log_fn: Optional callback for user-facing log messages.
                If None, only logger.debug is used.
    """
    def _log(msg):
        if log_fn:
            log_fn(msg)
        else:
            logger.debug(msg)

    loaded_folders = []

    for folder_path in folder_paths:
        folder_name = os.path.basename(folder_path)

        # Parse date from folder name (YYMMDD_...)
        try:
            folder_date = datetime.strptime(folder_name.split("_")[0], '%y%m%d').date()
        except ValueError:
            _log(f"   [SKIP] Folder '{folder_name}': invalid date format. Expected YYMMDD.")
            continue

        folder_data = MeasurementFolder(
            folder_name=folder_name,
            folder_path=folder_path,
            measurement_date=folder_date
        )

        files = [f for f in os.listdir(folder_path) if f.endswith(('.xlsx', '.xlsm'))]
        for file_name in files:
            file_path = os.path.join(folder_path, file_name)
            is_imported = False

            try:
                xls = pd.ExcelFile(file_path)
                sheet_names = xls.sheet_names

                # Protocol file
                if "Protocol" in sheet_names:
                    protocol = extract_protocol_info(xls, file_name)
                    if protocol and protocol.exp_date == folder_date:
                        folder_data.protocol = protocol
                        _log(f"   [PROTOCOL] loaded: {file_name}")
                        is_imported = True
                    elif protocol:
                        _log(f"   [MISMATCH] Protocol {protocol.exp_date} != Folder {folder_date}")

                # Measurement file
                elif "Table All Cycles" in sheet_names and set(sheet_names).issubset(
                        {"Table All Cycles", "Protocol Information"}):
                    result = extract_measurement_data(xls, file_name)
                    if result and result.measurement_date == folder_date:
                        # Extract optional Protocol Information sheet
                        if "Protocol Information" in sheet_names:
                            result.info_sheet = extract_info_sheet_data(xls)

                        folder_data.results.append(result)
                        _log(f"   [MEASUREMENT] loaded: {file_name}")
                        is_imported = True
                    elif result:
                        _log(f"   [MISMATCH] Analysis {result.measurement_date} != Folder {folder_date}")

                if not is_imported:
                    folder_data.skipped_files.append(file_name)

            except Exception as e:
                _log(f"   [ERROR] Could not read {file_name}: {e}")
                folder_data.skipped_files.append(file_name)

        loaded_folders.append(folder_data)
        logger.debug(f"Folder '{folder_name}': {len(folder_data.results)} results, "
                     f"protocol={'yes' if folder_data.protocol else 'no'}, "
                     f"skipped={len(folder_data.skipped_files)}")

    return loaded_folders