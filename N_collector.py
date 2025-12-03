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

    def dissect_folder_name(self, folder_name):
        """
        Seperates a folder name into components (expected format: Date_Experiment_Ncount).
        Returns (Date, Experiment, Ncount)
        """
        parts = folder_name.split("_")
        # Check expected structure
        if len(parts) >= 3 and parts[-1].startswith("n") and parts[0].isdigit():
            date = parts[0] # yymmdd format
            n_count = int(parts[-1][-1])
            # Experiment is everything in between date and N count
            experiment = "_".join(parts[1:-1])
            return date, experiment, n_count
        return None, folder_name, None

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

            # Search for xlsx in directory tree
            for root, dirs, files in os.walk(directory):
                has_xlsx = any(f.endswith(".xlsx") for f in files)

                if has_xlsx: # save path if xlsx files are found
                    self.subfolder_paths_with_files.append(root)

            # Update GUI
            if self.subfolder_paths_with_files:
                count = len(self.subfolder_paths_with_files)
                folder_names = [os.path.basename(path) for path in self.subfolder_paths_with_files]
                folder_names_string = "\n ".join(folder_names)
                self.folder_path.set(
                    f"Selected Path: {directory}\n\n Found following subfolders with xlsx files:\n {folder_names_string}")
                self.collect_button.config(state="normal")
                print(f"Found {count} folders: {folder_names_string}")
            else:
                self.folder_path.set(f"Error: No .xlsx files found in {directory} or any subfolder.")
                self.collect_button.config(state="disabled")
                print(f"No .xlsx files found starting from: {directory}")

    def collect_files(self):
        """Imports xlsx files found in the subfolders, separating protocol and result
        analysis files based on name matching rules."""
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

            # Get folder details and initialize storage for this folder
            folder_date, experiment, n_count = self.dissect_folder_name(folder_name)

            print(f"\n--- Processing Folder: {folder_name} ---")
            print(f"   Details: Date={folder_date}, Experiment='{experiment}', N={n_count}")

            protocol_file_name = f"{folder_name}.xlsx"
            # Initialise dic per folder
            folder_data = {'protocol': {'transfection_scheme': None}, 'results': []}
            skipped_files = []

            # Iterate over files in subfolder
            files_in_folder = os.listdir(folder_path)
            for file_name in files_in_folder:
                # Only look at xlsx files
                if not file_name.endswith('.xlsx'):
                    continue
                file_path = os.path.join(folder_path, file_name)
                is_imported = False

                try:
                    # Check for the protocol file (same name as folder)
                    if file_name == protocol_file_name:
                        # Validate that date of folder and protocol match
                        if file_name.startswith(folder_date):
                            df = pd.read_excel(file_path)
                            folder_data['protocol'] = df
                            # Get transfection scheme
                            transfection_df = self.extract_transfection_scheme(file_path)

                            if transfection_df is not None:
                                folder_data['protocol'] = {'transfection_scheme': transfection_df}
                            print(f"   [PROTOCOL] Imported: {file_name}")
                            is_imported = True
                        else:
                            print(f"   [DATE MISMATCH ERROR] Protocol file name {file_name} does not start with folder date {folder_date}.")

                    # Check for results file (_analysis in filename)
                    elif "_analysis" in file_name:
                        df = pd.read_excel(file_path)
                        # Validate that date of folder and measuring data match
                        metadata = self.extract_metadata(file_path)

                        # Convert folder date format to analysis sheet format
                        folder_date_formated = datetime.strptime(folder_date, '%y%m%d').strftime('%d/%m/%Y')

                        if metadata['measurement_date'] == folder_date_formated:
                             print(f"   [RESULT] Imported: {file_name} (ID2: {metadata['cell_line']}, ID3: {metadata['transfections']})")
                             # Store the data and its metadata only then
                             folder_data['results'].append({'df': df, 'meta': metadata})
                             is_imported = True
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