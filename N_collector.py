import os
import tkinter as tk
from tkinter import filedialog
import pandas as pd
from datetime import datetime
from dataclasses import dataclass, field
from typing import List, Optional

# --- Tool Functions --- #
def extract_transfection_scheme(df):
    """
    Uses the pd imported 'Protocol' sheet and finds the transfection table by identifying
    'Transfection scheme:'. Returns a pandas dataframe of the transfection scheme.
    """
    # Marker to find transfection table
    row_marker = "Transfection scheme:"
    col_end_marker = "vol per transfection"

    # Find the row index of the marker; idxmax() gets first occurrence of the largest value (True as 1 and False as 0)
    contains_marker = df.iloc[:, 0].astype(str).str.contains(row_marker, na=False)

    if not contains_marker.any():
        # The row_marker was not found
        print(f"   [ERROR] Transfection table not found in the Protocol sheet by looking for {row_marker}.")
        return None
    row_match = contains_marker.idxmax()

    # Table starts 3 rows after the marker;
    header_row_idx = row_match + 3
    # Look for col cutoff in header row
    header_row_content = df.iloc[header_row_idx].astype(str).str.lower()
    # Find the first column index that contains col end marker
    # Note: list comprehension is safer than .idxmax() on a row
    col_end_match = [i for i, val in enumerate(header_row_content) if col_end_marker in val]

    if not col_end_match:
        print(f"   [ERROR] Column '{col_end_marker}' not found. No end of transfection table found.")
        return None
    else:
        col_cutoff = col_end_match[0]

    # Get the transfection table from entire sheet
    df_transfection = df.iloc[header_row_idx + 1:, :col_cutoff].reset_index(drop=True)
    df_transfection.columns = df.iloc[header_row_idx, :col_cutoff].values # Set header

    # Define end of transfection scheme (first row with NA in first col)
    is_dna_na = df_transfection['DNA'].isna()
    if is_dna_na.any():
        # Find the positional index of the first NA value
        first_na_position = is_dna_na.values.argmax()

        # Slice the DataFrame using .iloc up to the row immediately before the NA row (exclusive)
        df_transfection = df_transfection.iloc[:first_na_position]
    else:
        # If no NA is found, raise an error as the detection of the table end has failed
        print("   [ERROR] Expected table end delimiter (NaN in 'DNA' column) not found.")
        return None

    # Validate expected columns are present (DNA, DB#, Conc)
    required_cols = ['DNA', 'DB#', 'Conc (ng/uL)']
    if not all(col in df_transfection.columns for col in required_cols):
        print(f"   [ERROR] Transfection table missing expected columns: {required_cols}.")
        return None

    return df_transfection

# TODO: move extract_metadata() and extract_protocol_info() outside of app!

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

    def extract_metadata(self, file_path):
        """ Reads the 'Table All Cycles' sheet from an analysis file and extracts metadata from the first column
        (measurement date, ID2: cell line, ID3: transfections #). Returns dictionary of extracted metadata.
        """
        metadata_worksheet = "Table All Cycles"

        try:
            # Read the first column of this sheet
            df_meta = pd.read_excel(file_path,
                                    sheet_name = metadata_worksheet,
                                    header = None,
                                    usecols = [0],
                                    nrows = 30 # Limit rows to read
            )
            col = df_meta[0].astype(str) # Transform everything to str

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
            print(f"[ERROR] Worksheet '{metadata_worksheet}' not found in file.")
            return None

    def extract_protocol_info(self, file_path):
        """
        Reads the 'Protocol' sheet of the protocol file and extracts all needed information.
        Stores and returns ProtocolData class with all info.

        calls
            extract_transfection_scheme()
        """
        protocol_worksheet = "Protocol"

        try:
            protocol_sheet = pd.read_excel(file_path,
                                    sheet_name=protocol_worksheet,
                                    header=None)
        except ValueError:
            # Error if sheet is missing
            print(f"[ERROR] Worksheet '{protocol_worksheet}' not found in file.")
            return None

        # Transfection df
        df_transfection = extract_transfection_scheme(protocol_sheet)

        # TODO: add apects to collect (date, title, cell lines ...)

        # Load transfection scheme table
        try:
            df_transfection = pd.read_excel(file_path,
                                            sheet_name=protocol_worksheet,
                                            header=header_row)

            # Define end if transfection scheme table cols (first col header starting with "Unnamed")
            cols_to_keep = [col for col in df_transfection.columns if not str(col).startswith('Unnamed')]
            df_transfection = df_transfection[cols_to_keep]

            # Define end of transfection scheme (first row with NA in first col)
            is_dna_na = df_transfection['DNA'].isna()

            if is_dna_na.any():
                # Find the positional index of the first NA value
                first_na_position = is_dna_na.values.argmax()

                # Slice the DataFrame using .iloc up to the row immediately before the NA row (exclusive)
                df_transfection = df_transfection.iloc[:first_na_position]
            else:
                # If no NA is found, raise an error as the detection of the table end has failed
                print("   [ERROR] Expected table end delimiter (NaN in 'DNA' column) not found.")
                return None

            # Validate expected columns are present (DNA, DB#, Conc)
            required_cols = ['DNA', 'DB#', 'Conc (ng/uL)']
            if not all(col in df_transfection.columns for col in required_cols):
                print(f"   [ERROR] Transfection table missing expected columns: {required_cols}.")
                return None

        return df_transfection

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
                # Logical variable to keep track on skipped files
                is_imported = False

                try:
                    # Use ExcelFile to check sheet names
                    xl = pd.ExcelFile(file_path)
                    sheet_names = xl.sheet_names

                    # Identify Protocol File
                    if "Protocol" in sheet_names:
                        # --- TODO: date validation with folder date - only if this matches, store info

                        # Get protocol_info (------ for now only transfection_df in function)
                        transfection_df = self.extract_protocol_info(file_path)

                        if transfection_df is not None:
                            folder_data.protocol = ProtocolData(file_name=file_name,
                                                                transfection_scheme=transfection_df)
                        print(f"   [PROTOCOL] Imported: {file_name}")
                        is_imported = True

                    # Identify Analysis File
                    elif "Analysis" in sheet_names:
                        metadata = self.extract_metadata(file_path)

                        # Date validation with folder_date, import only then
                        # Convert folder date format to analysis sheet format
                        folder_date_formated = datetime.strptime(folder_date, '%y%m%d').strftime('%d/%m/%Y')
                        if metadata['measurement_date'] == folder_date_formated:
                            # --- TODO: Specify here that really just BRET ratio table is read!
                            bret_ratio_df = pd.read_excel(file_path, sheet_name="Analysis")

                            result_obj = PRresult(
                                file_name=file_name,
                                measurement_date=metadata['measurement_date'],
                                cell_line=metadata['cell_line'],
                                transfection=metadata['transfections'],
                                raw_bret_ratio_df=bret_ratio_df
                            )
                            folder_data.results.append(result_obj)

                            is_imported = True
                            print(f"   [RESULT] Imported: {file_name} (ID2: {metadata['cell_line']}, ID3: {metadata['transfections']})")

                        else:
                            print(f"   [DATE MISMATCH ERROR] in file {file_name}: Sheet Date {metadata['date_sheet']} != Folder Date {folder_date_formated}")

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