import os
import tkinter as tk
from tkinter import filedialog
import pandas as pd
from datetime import datetime
from dataclasses import dataclass, field
from typing import List, Optional

# --- Dataclass Definition --- #

@dataclass
class PRresult:
    """ Information from a single _analysis file """
    file_name: str
    measurement_date: str
    cell_line: str # ID2
    transfection: str # ID3
    raw_bret_ratio_df: pd.DataFrame

@dataclass
class ProtocolData:
    """Information from a protocol file"""
    file_name: str
    transfection_scheme: pd.DataFrame
    ligand_conc: pd.DataFrame

    # add loads HERE-----------

@dataclass
class MeasurementFolder:
    """A subfolder containing one protocol and multiple result files"""
    folder_name: str
    folder_path: str
    measurement_date: str
    protocol: Optional[ProtocolData] = None # MeasurementFolder is initiated before protocol data is loaded
    results: List[PRresult] = field(default_factory=list) # The default_factory=list initiates this with an empty list
    skipped_files: List[str] = field(default_factory=list)

# --- Tool Functions --- #

def slice_table(
        df: pd.DataFrame,
        col_to_search: int,
        row_marker: str,
        second_table: bool = False
):
    """
    Extracts a table within a df without headers. Looks for row_marker in a specified column (col_to_search) to find
    the header line of the table. End of table is defined by first empty cell in header line and the first empty row
    in the first col of this table. Returns pandas dataframe with header.
    second_table: If True, looks for the 2nd occurrence of row_marker.
                  If False (default), uses the 1st occurrence.
    """
    # Find the row index of the marker (can be multiple);
    # returns True, False col; regex for handling row_markers containing ()
    contains_marker = df.iloc[:, col_to_search].astype(str).str.contains(row_marker, na=False, regex=False)
    marker_indices = contains_marker.index[contains_marker].tolist()

    if not marker_indices:
        # The row_marker was not found
        print(f"   [ERROR] {row_marker} was not found in protocol.")
        return None

    # Specify to look for first or second occurrence (needed for ligand tables)
    if second_table:
        if len(marker_indices) <2:
            print(f"   [ERROR] Second occurrence of {row_marker} requested, but only one found.")
            return None
        header_row_idx = marker_indices[1]  # Second occurrence
    else:
        header_row_idx = marker_indices[0]  # First occurrence

    # Look for col cutoff in header row
    header_row_content = df.iloc[header_row_idx, col_to_search:].astype(str)

    # Find the first index where the cell is 'nan' or empty; width of table is defined by header content
    col_width = next((i for i, val in enumerate(header_row_content)
                       if val.lower() == "nan" or not val.strip()), len(header_row_content))

    # Calculate the absolute end column
    col_end = col_to_search + col_width

    # Extract the table from large df
    df_extract = df.iloc[header_row_idx + 1:, col_to_search:col_end].reset_index(drop=True)
    df_extract.columns = df.iloc[header_row_idx, col_to_search:col_end].values  # Set header

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

def extract_transfection_scheme(df):
    """
    Uses the pd imported 'Protocol' sheet and finds the transfection table by identifying
    'DNA' in first column. Returns a pandas dataframe of the transfection scheme.
    """
    # Marker to find transfection table
    row_marker = "DNA"
    df_transfection = slice_table(df, 0, row_marker)
    # Remove not needed cols
    df_transfection_cleaned = df_transfection.drop(columns=["vol per transfection", "vol master"])

    # Validate expected columns are present (DNA, DB#, Conc)
    required_cols = ['DNA', 'DB#', 'Conc (ng/uL)']
    if not all(col in df_transfection_cleaned.columns for col in required_cols):
        print(f"   [ERROR] Transfection table missing expected columns: {required_cols}.")
        return None

    return df_transfection_cleaned

def extract_ligand_table(df):
    """
    Uses the pd imported 'Protocol' sheet and finds the ligand dilution table by identifying
    'dilution (1:)' in third column. Returns a pandas dataframe of the ligand dilution table.
    """
    # Marker to find ligand table
    row_marker = "final concentration in well (log(M))"
    df_ligand = slice_table(df, 10, row_marker)

    return df_ligand

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
    df_transfection = extract_transfection_scheme(protocol_sheet)
    ligand_conc = extract_ligand_table(protocol_sheet)

    protocol_info = ProtocolData(file_name=file_name,
                                 transfection_scheme=df_transfection,
                                 ligand_conc=ligand_conc)

    # TODO: add aspects to collect (date, title, cell lines ...)

    # Date of measurement

    # N
    # Experiment title
    # Cell line layout

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
        df_bret = slice_table(bret_sheet, 1, row_marker)

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
            folder_date = folder_name.split("_")[0]

            # Create MeasurementFolder object to collect protocol and results
            folder_data = MeasurementFolder(folder_name=folder_name, folder_path=folder_path, measurement_date=folder_date)

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
                        # --- TODO: date validation with folder date - only if this matches, store info

                        protocol_info = extract_protocol_info(xls)

                        if protocol_info is not None:
                            folder_data.protocol = protocol_info
                            print(f"   [PROTOCOL] Imported: {file_name}")
                        is_imported = True

                    # Identify Analysis File
                    elif "Analysis" in sheet_names:
                        meas_data = extract_measurement_data(xls)

                        # Date validation with folder_date
                        # Convert folder date format to analysis sheet format
                        folder_date_formated = datetime.strptime(folder_date, '%y%m%d').strftime('%d/%m/%Y')
                        if  meas_data.measurement_date == folder_date_formated:
                            folder_data.results.append(meas_data)
                            is_imported = True
                            print(f"   [RESULT] Imported: {file_name} (ID2: {meas_data.cell_line}, ID3: {meas_data.transfection})")
                        else:
                            print(f"   [DATE MISMATCH ERROR] in file {file_name}: Sheet Date {meas_data.measurement_date} != Folder Date {folder_date_formated}")

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