import os
import tkinter as tk
from tkinter import filedialog
import pandas as pd
from datetime import datetime

class NCollectorApp:
    def __init__(self, main_window):
        self.master = main_window
        main_window.title("N Collector")

        # Path to folder variable
        self.folder_path = tk.StringVar()
        self.folder_path.set("No folder selected.")
        self.subfolder_paths_with_files = []

        # Display label for path
        self.path_label = tk.Label(main_window,
                                   textvariable=self.folder_path,
                                   wraplength=1000,
                                   justify=tk.LEFT,
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
        metadata = {}
        metadata_worksheet = "Table All Cycles"

        try:
            # Read the first column of this sheet
            df_meta = pd.read_excel(file_path,
                                    sheet_name = metadata_worksheet,
                                    header = None,
                                    usecols = [0]
            )
        except ValueError:
            # Error if sheet is missing
            print(f"[ERROR] Worksheet '{metadata_worksheet}' not found in file.")
            return None

        # Iterate through rows of the first column to find metadata
        for index, row in df_meta.iterrows():
            line = str(row[0]).strip()  # Get the string value from the cell

            # Extract Date
            if line.startswith("Date:"):
                # Expects format like "Date: 21/11/2025"
                metadata['measurement_date'] = line.split(":", 1)[-1].strip()

            # Extract ID2 (Condition 1)
            elif line.startswith("ID2:"):
                # Expects format like "ID2: Con"
                metadata['cell_line'] = line.split(":", 1)[-1].strip()

            # Extract ID3 (Condition 2)
            elif line.startswith("ID3:"):
                # Expects format like "ID3: 1,2,3,4"
                metadata['transfections'] = line.split(":", 1)[-1].strip()

            # Ensure all required metadata fields were found
        if 'measurement_date' not in metadata or 'cell_line' not in metadata or 'transfections' not in metadata:
            print("   [WARNING] Missing Date, ID2, or ID3 from metadata sheet.")
            return None

        return metadata

    def extract_transfection_scheme(self, file_path):
        """
        Reads the 'Protocol' sheet of the protocol file and finds the transfection table by identifying
        'Transfection scheme:'. Returns a pandas dataframe of the transfection scheme.
        """
        protocol_worksheet = "Protocol"
        row_marker = "Transfection scheme:"

        try:
            df_col1 = pd.read_excel(file_path,
                                    sheet_name=protocol_worksheet,
                                    header=None,
                                    usecols=[0])
        except ValueError:
            # Error if sheet is missing
            print(f"[ERROR] Worksheet '{protocol_worksheet}' not found in file.")
            return None

        # Find the row index of the marker; idxmax() gets first occurrence of the largest value (True as 1 and False as 0)
        match = df_col1[0].astype(str).str.contains(row_marker, na=False).idxmax()
        if df_col1.loc[match, 0] != row_marker:
            # Check if idxmax found the marker or just the first row
            print(f"   [ERROR] Transfection table not found in the Protocol sheet by looking for {row_marker}.")
            return None

        # Table starts 3 rows after the marker;
        header_row = match + 3
        print(f"   [PROTOCOL] Transfection scheme marker found at row {match}. Loading table from row {header_row}.")

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
                # Find the positional index (0, 1, 2, ...) of the first NA value
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

        except Exception as e:
            print(f"   [ERROR] Failed to load transfection scheme table: {e}")
            return None


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

        # Create a dictionary to store all imported data for later processing
        # Dic structure: {'folder_name': {'protocol': df, 'results':[df1, df2, ...]
        all_data = {}

        # Iterate over subfolders
        for folder_path in self.subfolder_paths_with_files:
            folder_name = os.path.basename(folder_path)
            print(f"\n--- Processing Folder: {folder_name} ---")

            # Use date in folder name for validation
            folder_date = folder_name.split("_")[0]

            # Dic storage for this folder (one protocol, list of results)
            folder_data = {'protocol': None, 'results': []}
            # Keep track of skipped files
            skipped_files = []

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
                        ### Add extract_protocol_info() here
                        # Date validation with folder date - only if this matches, store all info in dic

                        # Get transfection scheme
                        transfection_df = self.extract_transfection_scheme(file_path)

                        if transfection_df is not None:
                            folder_data['protocol'] = dict(file_name=file_name, transfection_scheme=transfection_df)
                        print(f"   [PROTOCOL] Imported: {file_name}")
                        is_imported = True

                    # Identify Analysis File
                    elif "Analysis" in sheet_names:
                        metadata = self.extract_metadata(file_path)

                        # Date validation with folder_date, import only then
                        # Convert folder date format to analysis sheet format
                        folder_date_formated = datetime.strptime(folder_date, '%y%m%d').strftime('%d/%m/%Y')
                        if metadata['measurement_date'] == folder_date_formated:
                            result_analysis = pd.read_excel(file_path, sheet_name="Analysis")
                            folder_data['results'].append({
                                'file_name': file_name,
                                'Analysis': result_analysis,
                                'meta': metadata
                            })
                            is_imported = True
                            print(f"   [RESULT] Imported: {file_name} (ID2: {metadata['cell_line']}, ID3: {metadata['transfections']})")

                        else:
                            print(f"   [DATE MISMATCH ERROR] in file {file_name}: Sheet Date {metadata['date_sheet']} != Folder Date {folder_date_formated}")

                    # Files not matching criteria are skipped
                    if not is_imported:
                        skipped_files.append(file_name)
                except Exception as e:
                    print(f"   [ERROR] Could not read {file_name}: {e}")

            # Store the collected data for this experiment
            all_data[folder_name] = folder_data
            print(f"   [SKIPPED]: {skipped_files}")


# --- Main Execution Block ---
if __name__ == "__main__":
    # Create the main window
    root = tk.Tk()

    # Create an instance of the application
    app = NCollectorApp(root)

    # Start the Tkinter event loop
    root.mainloop()