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


    def select_folder(self):
            """Opens dialog to select folder to search for xlsx files in"""
            self.subfolder_paths_with_files = []  # Clear previous results

            directory = filedialog.askdirectory(title="Select a folder...")
            if not directory: # User closed dialog without selecting a folder
                self.folder_path.set(f"No folder selected.")
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
                    print(print(f"Found {count} folders: {folder_names}"))
                else:
                    self.folder_path.set(f"Error: No .xlsx files found in {directory} or any subfolder.")
                    print(f"No .xlsx files found starting from: {directory}")




# --- Main Execution Block ---
if __name__ == "__main__":
    # Create the main window
    root = tk.Tk()

    # Create an instance of the application
    app = NCollectorApp(root)

    # Start the Tkinter event loop
    root.mainloop()