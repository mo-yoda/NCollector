import os
import tkinter as tk
from tkinter import filedialog
import pandas as pd
from datetime import datetime, date
from dataclasses import dataclass, field
from typing import List, Optional

# --- Dataclass Definition --- #

@dataclass
class PRresult:
    """ Information from a single _analysis file """
    file_name: str
    measurement_date: date
    cell_line: str # ID2
    transfection: str # ID3
    raw_bret_ratio_df: pd.DataFrame

@dataclass
class ProtocolData:
    """Information from a protocol file"""
    file_name: str
    exp_date: date
    n: int
    cell_lines: List[str]
    line_layout: str
    transfection_scheme: pd.DataFrame
    ligand: str
    ligand_conc: pd.DataFrame
    # Second ligand is optional
    ligand_2: Optional[str] = None
    ligand_conc_2: Optional[pd.DataFrame] = None

@dataclass
class MeasurementFolder:
    """A subfolder containing one protocol and multiple result files"""
    folder_name: str
    folder_path: str
    measurement_date: date
    protocol: Optional[ProtocolData] = None # MeasurementFolder is initiated before protocol data is loaded
    results: List[PRresult] = field(default_factory=list) # The default_factory=list initiates this with an empty list
    skipped_files: List[str] = field(default_factory=list)

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

def extract_protocol_info(xls_obj: pd.ExcelFile):
    """
    Uses already opened pd.ExcelFiles (faster and more flexible than reading from path).
    Reads the 'Protocol' sheet of the protocol file and extracts all needed information.
    Stores and returns ProtocolData class with all info.

    calls
        extract_transfection_scheme
        extract_ligand_table
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
    if df_transfection is not None:
        df_transfection = df_transfection.drop(columns=["vol per transfection", "vol master"])

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
                                 ligand = ligand_1,
                                 ligand_conc=ligand_1_conc,
                                 ligand_2 = ligand_2,
                                 ligand_conc_2= ligand_conc_2)

    print(protocol_info)

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
        meta_keys = {
            "Date:": "measurement_date",
            "ID2:": "cell_line",
            "ID3:": "transfections"
        }

        metadata = {}
        for label, key in meta_keys.items():
            # n=1 split at first ":"; str[-1] select last arg; str-strip() remove spaces; .tolist() convert from pd series
            matches = col[col.str.contains(label, na=False)].str.split(":", n=1).str[-1].str.strip().tolist()
            if matches:
                metadata[key] = matches[0]

        # Transform date str to actual date
        if 'measurement_date' in metadata:
            try:
                date_str = metadata["measurement_date"]
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
        row_marker = "Time (min)"
        df_bret = slice_table(bret_sheet, row_marker)

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
    Stores and returns PRresult class with all data.

    calls
        extract_metadata
        extract_bret_data
    """

    file_name = os.path.basename(xls_obj.io)
    metadata_dic = extract_metadata(xls_obj)
    bret_ratio_df = extract_bret_data(xls_obj)

    result_obj = PRresult(
        file_name=file_name,
        measurement_date=metadata_dic['measurement_date'],
        cell_line=metadata_dic['cell_line'],
        transfection=metadata_dic['transfections'],
        raw_bret_ratio_df=bret_ratio_df
    )
    return result_obj

# --- Main Application --- #

class NCollectorApp:
    def __init__(self, main_window):
        self.master = main_window
        main_window.title("N Collector")

        # Path to folder variable
        self.folder_path = tk.StringVar()
        self.folder_path.set("No folder selected.")
        self.subfolder_paths_with_files = []
        self.experiment: List[MeasurementFolder] = []

        # Display label for path
        self.path_label = tk.Label(main_window,
                                   textvariable=self.folder_path,
                                   wraplength=1000,
                                   justify="left",
                                   font=('Arial', 10))
        self.path_label.pack(pady=10, padx=10) # placing the text via .pack

        # Select Folder button
        self.select_button = tk.Button(main_window,
                                       text="Select folder containing results of experiment",
                                       command=self.select_folder)
        self.select_button.pack(pady=10, padx=10)

        # Analyse button
        self.collect_button = tk.Button(main_window,
                                        text="Load Files",
                                        state="disabled",
                                        command=self.collect_files
        )
        self.collect_button.pack(pady=15)

    def select_folder(self):
            """Opens dialog to select folder to search for xlsx files in"""
            self.subfolder_paths_with_files = []  # Clear previous results

            directory = filedialog.askdirectory(title="Select a folder...")
            if not directory: # User closed dialog without selecting a folder
                self.folder_path.set(f"No folder selected.")
                self.collect_button.config(state="disabled")
                return

            # Search for xlsx or xlsm in directory tree
            for root, dirs, files in os.walk(directory):
                has_xlsx = any(f.endswith((".xlsx", ".xlsm")) for f in files)

                if has_xlsx: # save path if xlsx/xlsm files are found
                    self.subfolder_paths_with_files.append(root)

            # Update GUI
            if self.subfolder_paths_with_files:
                count = len(self.subfolder_paths_with_files)
                folder_names = [os.path.basename(path) for path in self.subfolder_paths_with_files]
                folder_names_string = "\n ".join(folder_names)
                self.folder_path.set(
                    f"Selected Path: {directory}\n\n Found following subfolders with xlsx/xlsm files:\n {folder_names_string}")
                self.collect_button.config(state="normal")
                print(f"Found {count} folders: {folder_names_string}")
            else:
                self.folder_path.set(f"Error: No .xlsx or .xlsm files found in {directory} or any subfolder.")
                self.collect_button.config(state="disabled")
                print(f"No .xlsx or .xlsm files found starting from: {directory}")

    def collect_files(self):
        """EDIT: Reads sheet names of all xlsx and xlsm files to identify and separate protocol and
        result analysis files. Only these identified files are read. Validation of correct protocol to 
        analysis files is done via date of measurement."""

        if not self.subfolder_paths_with_files:
            print("No folders to analyze.")
            return

        print("\n--- Starting Data Collection ---")

        # Collect all subfolder info (experiment repeats) in this list
        self.experiment = []

        # Iterate over subfolders
        for folder_path in self.subfolder_paths_with_files:
            folder_name = os.path.basename(folder_path)
            print(f"\n--- Processing Folder: {folder_name} ---")

            # Use date in folder name for validation
            try:
                # Transform date str to date obj (YYMMDD)
                folder_date_obj = datetime.strptime(folder_name.split("_")[0], '%y%m%d').date()
            except ValueError:
                print(f"   [ERROR] Folder '{folder_name}' invalid date format. Expected YYMMDD. Skipping.")
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
                            print(f"   [PROTOCOL] Imported: {file_name}")
                            is_imported = True
                        elif protocol_info:
                            print(f"   [MISMATCH] Protocol {protocol_info.exp_date} != Folder {folder_date_obj}")

                    # Identify Analysis File
                    elif "Analysis" in sheet_names:
                        meas_data = extract_measurement_data(xls)

                        if meas_data and meas_data.measurement_date == folder_date_obj:
                            folder_data.results.append(meas_data)
                            print(f"   [RESULT] Imported: {file_name} (ID2: {meas_data.cell_line}, ID3: {meas_data.transfection})")
                            is_imported = True
                        else:
                            print(f"   [MISMATCH] Analysis {meas_data.measurement_date} != Folder {folder_date_obj}")

                    # Files not matching criteria are skipped
                    if not is_imported:
                        folder_data.skipped_files.append(file_name)
                except Exception as e:
                    print(f"   [ERROR] Could not read {file_name}: {e}")

            # Store the collected data for this experiment
            self.experiment.append(folder_data)
            print(f"   [SKIPPED]: {folder_data.skipped_files}")

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


# --- Main Execution Block ---
if __name__ == "__main__":
    # Create the main window
    root = tk.Tk()

    # Create an instance of the application
    app = NCollectorApp(root)

    # Start the Tkinter event loop
    root.mainloop()