import os
import tkinter as tk
from tkinter import filedialog

import pandas as pd

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
            date = parts[0]
            n_count = int(parts[-1][-1])
            # Experiment is everything in between date and N count
            experiment = "_".join(parts[1:-1])
            return date, experiment, n_count
        return None, folder_name, None


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
                    print(print(f"Found {count} folders: {folder_names}"))
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
            date, experiment, n_count = self.dissect_folder_name(folder_name)

            print(f"\n--- Processing Folder: {folder_name} ---")
            print(f"   Details: Date={date}, Experiment='{experiment}', N={n_count}")

            protocol_file_name = f"{folder_name}.xlsx"
            # Initialise dic per folder
            folder_data = {'protocol': None, 'results': []}
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
                        df = pd.read_excel(file_path)
                        folder_data['protocol'] = df
                        print(f"   [PROTOCOL] Imported: {file_name}")
                        is_imported = True

                    # Check for results file (_analysis in filename)
                    elif "_analysis" in file_name:
                        df = pd.read_excel(file_path)
                        folder_data['results'].append(df)
                        print(f"   [RESULT] Imported: {file_name}")
                        is_imported = True
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