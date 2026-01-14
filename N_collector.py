import os
import tkinter as tk
from tkinter import filedialog, ttk
import pandas as pd
from datetime import datetime, date
from dataclasses import dataclass, field

# --- Dataclass Definition --- #

@dataclass
class PrResult:
    """ Information from a single _analysis file """
    # Information from analysis xlsx itself
    file_name: str
    measurement_date: date
    cell_line: str # ID2
    transfection: str # ID3
    raw_bret_ratio_df: pd.DataFrame

    # Connection to protocol file
    exp_conditions: list[str] = field(default_factory=list)
    # Stores baseline- and vehicle-normalised BRET ratios
    processed_df: pd.DataFrame | None = None

    # User interaction (optionally excluding a plate)
    is_excluded: bool = False

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
    transfection_conditions: dict[str, list[str]]
    ligand: str
    ligand_conc: pd.DataFrame
    # Second ligand is optional; by | None = None
    ligand_2: str | None = None
    ligand_conc_2: pd.DataFrame | None = None

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
                                 main_plasmids= main_dna,
                                 transfection_conditions=transfection_conditions,
                                 ligand = ligand_1,
                                 ligand_conc=ligand_1_conc,
                                 ligand_2 = ligand_2,
                                 ligand_conc_2= ligand_conc_2)
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

            # Re-join unique sorted parts (e.g. "Control, dQ")
            metadata['cell_line'] = ", ".join(sorted(list(set(standardised_lines))))

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
        transfection=metadata_dic['transfections'],
        raw_bret_ratio_df=bret_ratio_df
    )
    return result_obj

def process_bret_measurement(result: PrResult, protocol: ProtocolData):
    """
    Performs Baseline Correction and Vehicle Normalization.
    Calculates AUC.
    Returns (processed_df, stats_df)
    """
    # Get raw BRET ratio table
    raw_df = result.raw_bret_ratio_df.copy()
    time_col = "Time (min)"

    # Define Baseline: First 5 rows (Index 0-4)
    baseline_end_idx = 5
    # Define Kinetic: Remaining rows (Index 5 onwards)
    if len(raw_df) < (baseline_end_idx):
        print(f"   [WARNING] Data has less than {(baseline_end_idx)} rows.")
        return None, None

    # Prepare the full dataframe for normalization (removing the time col)
    data_df = raw_df.drop(columns=[time_col]).copy()

    # --- BASELINE CORRECTION ---
    # Isolate baseline rows for calculation
    df_baseline_calc = raw_df.iloc[0:baseline_end_idx].copy()
    bl_corrected_df = data_df.copy()

    for col in data_df.columns:
        # Mean of baseline rows for this well
        # Force numeric conversion for baseline values to handle potential strings/decimals
        base_vals = pd.to_numeric(df_baseline_calc[col], errors='coerce')
        base_mean = base_vals.mean()

        if base_mean != 0:
            col_vals = pd.to_numeric(data_df[col], errors='coerce')
            # Divide ALL values by the baseline mean
            bl_corrected_df[col] = col_vals / base_mean
        else:
            bl_corrected_df[col] = None

    # --- SPECIFY LAYOUT ---
    # Define starting inx of replicate cols
    triplicate_block_starts = [1,4,7,10]
    # Layout map; which col # belong to which starting block
    plate_layout_map = {}
    for start in triplicate_block_starts:
        for offset in range(3):  # 0, 1, 2 (Triplicates)
            col_idx = start + offset
            plate_layout_map[col_idx] = start

    # --- VEHICLE CORRECTION ---
    # Calculate mean vehicle for each condition
    vehicle_mean = {}

    for start in triplicate_block_starts:
        wells = [f"H{start + k}" for k in range(3)]
        # Filter for wells that actually exist in the dataframe
        # Logic needed for optionally excluding wells
        valid_wells = [w for w in wells if w in bl_corrected_df.columns]

        if valid_wells:
            # Get vehicle values and make sure that data is numeric
            vehicle_data = bl_corrected_df[valid_wells].apply(pd.to_numeric, errors='coerce')
            # TODO: check if vehicle_date is out of bounds -> recommend exclusion of specific wells; based on mean or for each time point?
            # Average across valid vehicle wells per time point
            vehicle_mean[start] = vehicle_data.mean(axis=1)
            # print(f"vehicle values {vehicle_mean[start]}")

    # Normalise to mean(vehicle)
    vehicle_corr_df = bl_corrected_df.copy()

    for col in data_df.columns:
        # Get well number
        col_num_str = col[1:]

        try:
            col_num = int(col_num_str)

            # Check map
            if col_num in plate_layout_map:
                block_start = plate_layout_map[col_num]

                # Normalize to mean(vehicle) for each block
                vehicle_corr_df[col] = bl_corrected_df[col] / vehicle_mean[block_start]

        except ValueError:
            continue
    # TODO: add AUC calculation
    return vehicle_corr_df

# TODO: add function to rearrange cols of processed bret df (flexible for user interaction)

# --- Main Application --- #

class NCollectorApp:
    def __init__(self, main_window):
        self.master = main_window
        main_window.title("N Collector")

        # --- GUI State ---
        # Path to folder variable
        self.folder_path = tk.StringVar(value="No folder selected.")
        self.subfolder_paths_with_files = []
        self.experiment: list[MeasurementFolder] = []

        # --- TABS SETUP ---
        self.notebook = ttk.Notebook(main_window)
        self.notebook.pack(expand=True, fill='both')

        # Tab 1: Import Data
        self.tab_import = tk.Frame(self.notebook)
        self.notebook.add(self.tab_import, text="1. Import Data")

        # Tab 2: Data Selection
        self.tab_select = tk.Frame(self.notebook)
        self.notebook.add(self.tab_select, text="2. Data Selection")

        # --- TAB 1 CONTENT ---
        # Display label for path
        self.path_label = tk.Label(self.tab_import,
                                   textvariable=self.folder_path,
                                   wraplength=1000,
                                   justify="left",
                                   font=('Arial', 10))
        self.path_label.pack(pady=10, padx=10) # placing the text via .pack

        # Select Folder button
        # No self. needed as this does not have to be stored for later changes
        tk.Button(self.tab_import, text="Select folder containing results of experiment",
                  command=self.select_folder).pack(pady=10, padx=10)

        # Analyse button
        self.collect_button = tk.Button(self.tab_import,
                                        text="Load Files",
                                        state="disabled",
                                        command=self.collect_files
        )
        self.collect_button.pack(pady=15)

    def select_folder(self):
            """Opens dialog to select folder to search for xlsx files in"""
            directory = filedialog.askdirectory(title="Select a folder...")
            if not directory: return

            # --- RESET STATE: Clear old data when a new folder is selected
            self.subfolder_paths_with_files = []
            self.experiment = []

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
                self.collect_button.config(state="normal")
                print(f"Found {count} folders: \n {folder_names_string}")
            else:
                self.folder_path.set(f"Error: No .xlsx or .xlsm files found in {directory} or any subfolder.")
                self.collect_button.config(state="disabled")
                print(f"No .xlsx or .xlsm files found starting from: {directory}")

    def map_conditions_to_results(self):
        """
        Iterates through all loaded experiments and resolves the numerical ID3
        into actual conditions using the Protocol information.
        """
        print("\n--- Resolving Experimental Conditions ---")
        for folder in self.experiment:
            if not folder.protocol:
                continue

            # Dic of # transfection to condition {'1': ['plasmid_A', 'plasmid_B'], '2': ...}
            mapping = folder.protocol.transfection_conditions

            for result in folder.results:
                if result.transfection:
                    # Split ID3 by comma and strip whitespace
                    ids = [x.strip() for x in str(result.transfection).split(',')]
                else:
                    ids = []

                condition_found = []
                for i in ids:
                    # Look up ID in the protocol mapping, [] list as fallback
                    plasmids = mapping.get(i, [f"Unknown_ID3_part_{i}"])
                    cond_name = " + ".join(sorted(plasmids)) # Handling co-transfection
                    condition_found.append(cond_name)

                # Sort to ensure "Rab5 + b2AR" is treated same as "b2AR + Rab5" if order implies same condition
                result.exp_conditions = sorted(condition_found)

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
        # TODO: test this function as user!

        # Create a modal dialog window
        dialog = tk.Toplevel(self.master)
        dialog.title("Select Experiment")
        dialog.geometry("400x300")

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
        self.master.wait_window(dialog) # Wait until the window is closed

        # Retrieve selection
        selected_key_str = selected_var.get()
        selected_key = str_to_key_map.get(selected_key_str, first_key)

        # Filter the experiment list
        self.experiment = main_plasmids_groups[selected_key]
        print(f"   [FILTER] Keeping {len(self.experiment)} folders matching main plasmids: {selected_key}")

        return " + ".join(selected_key)

    def aggregate_experiments(self, main_plasmids_name):
        """
        Groups results by (Cell Line, Condition). Returns a dictionary of groups as preparation
        for optional exclusion by User.
        """
        print(f"\n--- Collecting Ns for measurements with {main_plasmids_name} ---")

        # Nested dic as planned treeview GUI expects this
        grouped_data = {}

        # sth here takes ages
        for folder in self.experiment:
            for result in folder.results:
                # TODO: for CKs layout, there are several cell lines possible -> edit later to handle this
                # Get ID2: cell_lines
                cell_line = result.cell_line

                if not result.exp_conditions:
                    # Handle case where no conditions were mapped or ID3 was empty
                    cond_list = ["Undefined Condition"]
                else:
                    cond_list = result.exp_conditions

                if cell_line not in grouped_data:
                    grouped_data[cell_line] = {}

                # Iterate through each separate condition found in the result file
                for cond_name in cond_list:
                    if cond_name not in grouped_data[cell_line]:
                        grouped_data[cell_line][cond_name] = []

                    # TODO: not the entire result has to be appended, only the specific condition! -> BRET processing has to be done first
                    grouped_data[cell_line][cond_name].append(result)
        return grouped_data

    def collect_files(self):
        """
        Reads sheet names of all xlsx and xlsm files to identify and separate protocol and result analysis files.
        Validation of correct protocol to analysis files is done via date of measurement in the folder name.
        """

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
                            meas_data.processed_df = process_bret_measurement(meas_data, folder_data.protocol)
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

        # --- Connecting Protocol and Analysis files ---
        # Get the transfected plasmids to assign conditions
        self.map_conditions_to_results()
        # Check for plasmids transfected in all conditions (main plasmids) and filter if needed
        selected_exp_name = self.handle_main_plasmids_selection()
        # Collect Ns
        grouped_results = self.aggregate_experiments(selected_exp_name)

        print("\n" + "=" * 40)
        print("AGGREGATED DATA SUMMARY (N COUNTS)")
        print("=" * 40)

        if not grouped_results:
            print("No data aggregated.")

        for cell_line, conditions in grouped_results.items():
            print(f"\nCell Line: {cell_line}")
            for cond_name, results_list in conditions.items():
                n_count = len(results_list)
                print(f"  • Condition: {cond_name}")
                print(f"      -> N = {n_count}")
                # Optional: Show which days contributed
                days = sorted([r.measurement_date.strftime('%y%m%d') for r in results_list])
                print(f"      -> Days: {', '.join(days)}")

        print("\n" + "=" * 40)

# TODO: implement window to show any ERROR messages + add optional export of log file

# --- Main Execution Block ---
if __name__ == "__main__":
    # Create the main window
    root = tk.Tk()

    # Create an instance of the application
    app = NCollectorApp(root)

    # Start the Tkinter event loop
    root.mainloop()