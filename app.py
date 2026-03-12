import os
import logging
import tkinter as tk
from tkinter import filedialog, ttk
import pandas as pd
from datetime import datetime, date

from models import MeasurementFolder, ProcessingConfig
from parsing import extract_protocol_info, extract_measurement_data
from processing import process_bret_measurement
from export import apply_export_filters, build_row_info_str, generate_header_key, create_clean_pivot

logger = logging.getLogger("NCollector")

# --- Main Application --- #

class NCollectorApp:
    def __init__(self, main_window):
        self.main_gi = main_window
        self.version = "N Collector v2.0 Beta"
        main_window.title(self.version)
        main_window.geometry("800x700")

        # --- Data Storage ---
        self.directory = ""
        self.subfolder_paths_with_files = []
        self.experiment: list[MeasurementFolder] = []
        self.master_index = pd.DataFrame()  # Index for populating tab 2
        self.rule_history_text = ""
        self.pending_exclusions = []
        self.master_df = pd.DataFrame()  # Used for master csv file storage (by generation or import)
        self.ignored_warnings = set()
        self.current_config = None

        # --- GUI Variables ---
        self.folder_path = tk.StringVar(value="No folder selected.")
        self.var_labeling_is_checked = tk.BooleanVar(value=False)
        self.var_lum_threshold = tk.IntVar(value=100)
        self.var_lig = tk.StringVar(value="")
        self.var_date = tk.StringVar(value="All")
        self.var_cell = tk.StringVar(value="All")
        self.var_cond = tk.StringVar(value="All")
        self.var_repl = tk.StringVar(value="All")
        self.var_row = tk.StringVar(value="All")
        self.var_category = tk.StringVar(value="")
        self.var_specific_type = tk.StringVar(value="")
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
        self.combo_category = None
        self.combo_specific = None
        self.lb_kin_layout = None
        self.btn_run_plot_helper = None

        # --- Constants ---
        self.data_type_map = {
            "kinetic": {
                "raw BRET ratio": "Raw_BRET_kinetic",
                "labeling-corrected BRET ratio": "Lab_BRET_kinetic",
                "baseline-corrected BRET ratio": "Bl_Corrected_BRET",
                "vehicle-normalised BRET ratio, techn. replicates": "Veh_Norm_Kinetic",
                "vehicle-normalised BRET ratio, mean of techn. replicates": "Kinetic_Mean"
            },
            "CRC": {
                "last 3x timepoints: raw BRET": "Raw_BRET_CRC",
                "last 3x timepoints: labeling-corrected BRET ratio": "Lab_LP",
                "last 3x timepoints: baseline-corrected BRET ratio": "Bl_LP",
                "last 3x tp: vehicle-normalised BRET ratio, techn. replicates": "Veh_Norm_LP",
                "last 3x tp: vehicle-normalised BRET ratio, mean of techn. replicates": "LP_Mean",
                "AUC: vehicle-normalised BRET ratio, techn. replicates": "Veh_Norm_AUC",
                "AUC: vehicle-normalised BRET ratio, mean of techn. replicates": "AUC_Mean"
            }
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
        self.path_label.pack(pady=0, padx=15, anchor="nw")

        # Load frame
        load_frame = tk.Frame(self.tab_import)
        load_frame.pack(pady=5, fill="x", padx=10)
        # Weight setting to place load button in the middle
        load_frame.columnconfigure(0, weight=4)
        load_frame.columnconfigure(1, weight=1)

        # Load button
        self.load_files_button = tk.Button(load_frame, text="Load Files", state="disabled",
                                           command=self.collect_files)
        # self.load_files_button.pack(padx=10)
        self.load_files_button.grid(row=0, column=0, padx=5, pady=0, sticky="sew")

        # Loading specs
        load_settings_frame = tk.LabelFrame(load_frame, text="Loading Specs")
        load_settings_frame.grid(row=0, column=1, sticky="ew")

        # Labeling correction checkbox
        chk_container = tk.Frame(load_settings_frame)
        chk_container.pack(anchor="ne")
        # Container needed as text is only supported right of chbx

        lbl_correction = tk.Label(chk_container, text="Labeling correction")
        lbl_correction.pack(side="left", padx=(0, 5))
        labeling_chk = tk.Checkbutton(chk_container,
                                      variable=self.var_labeling_is_checked,
                                      command=self.on_checkbox_toggle)
        labeling_chk.pack(side="right")

        # Threshold Input
        lum_thresh_entry = tk.Entry(load_settings_frame, textvariable=self.var_lum_threshold, width=10)
        lum_thresh_entry.pack(side="right", padx=(10, 5), pady=(0,8))
        tk.Label(load_settings_frame, text="Lum. Threshold:").pack(side="right", padx=(10, 5))

        # Collected Ns frame
        loaded_data_frame = tk.LabelFrame(self.tab_import, text="Loaded Data")
        loaded_data_frame.pack(fill="both", expand=True,  padx=10, pady=5)

        # Label to display Main Plasmids
        self.main_plasmids_label = tk.Label(loaded_data_frame, text="", justify="left", font=("Arial", 10, "bold"))
        self.main_plasmids_label.pack(pady=(0, 5))

        # Display of N summary table
        summary_frame = tk.Frame(loaded_data_frame)
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
        rules_frame.pack(fill="x", padx=10, pady=5)
        self.lbl_rules_summary = tk.Label(rules_frame, text="No exclusion rules applied", justify="left", anchor="w")
        self.lbl_rules_summary.pack(fill="x", padx=5, pady=5)

        # --- EXPORT SECTION ---
        export_frame = tk.LabelFrame(self.tab_import, text="Export Options")
        export_frame.pack(fill="x", padx=10, pady=10)
        self.btn_export_master = tk.Button(export_frame, text="Export Master CSV", state="disabled",
                                           command=self.export_master_csv)
        self.btn_export_master.pack(side="left", fill="x", expand=True, padx=5, pady=10)
        self.btn_export_excel = tk.Button(export_frame, text="Export Excel Report (Default)", state="disabled",
                                          command=self.export_excel_report)
        self.btn_export_excel.pack(side="left", fill="x", expand=True, padx=5, pady=10)

    def on_checkbox_toggle(self):
        if self.var_labeling_is_checked.get():
            self.log("[LABELING CORRECTION ENABLED]   Press load files to apply.")
        else:
            self.log("[LABELING CORRECTION DISABLED]")

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
        self.var_repl = tk.StringVar(value="All")
        self.var_row = tk.StringVar(value="All")
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
        self.cb_rep = ttk.Combobox(granular_frame, textvariable=self.var_repl, state="readonly", width=12)
        self.cb_rep.pack(side="left", padx=5)
        self.cb_rep.bind("<<ComboboxSelected>>", self.toggle_row_dropdown)

        # Row dropdown (Initially disabled/hidden until Replicate is picked)
        self.lbl_row = tk.Label(granular_frame, text="Specific Row:")
        self.lbl_row.pack(side="left", padx=5)

        self.cb_row = ttk.Combobox(granular_frame, textvariable=self.var_row, state="disabled", width=5)
        self.cb_row.pack(side="left", padx=5)
        self.cb_row['values'] = ["All", "A", "B", "C", "D", "E", "F", "G", "H"]

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
            logger.info(message)

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
        self.update_dropdown_options("Replicate")

    def refresh_filter_options(self):
        """Called during built master index. Updates dropdown options of date, cell line and condition."""
        if self.master_index.empty: return

        # Reset Variables
        self.var_lig.set("All")
        self.var_date.set("All")
        self.var_cell.set("All")
        self.var_cond.set("All")
        self.var_repl.set("All")
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
            "Condition": (self.var_cond, self.cb_cond),
            "Replicate": (self.var_repl, self.cb_rep)
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
                    mask &= (self.master_index[curr_param].astype(str) == str(val))

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
            self.var_repl.set("All")
            self.toggle_row_dropdown()

    def add_exclusion_rule(self):
        """
        Adds the current dropdown state to the pending list.
        Handles empty strings for Replicate/Col and Row.
        """
        rule = {
            "Ligand": self.var_lig.get(),
            "Date": self.var_date.get(),
            "Cell_Line": self.var_cell.get(),
            "Condition": self.var_cond.get(),
            "Replicate": self.var_repl.get(),
            "Row": self.var_row.get()
        }

        # Check for duplicates or empty
        rule_str = f"Ligand: {rule['Ligand']} | Date: {rule['Date']} | "\
                   f"Cell: {rule['Cell_Line']} | Cond: {rule['Condition']} | "\
                   f"Rep:{rule['Replicate']} | Row:{rule['Row']}"

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
            # Check whether all ligands or the only one possible is selected
            ligand_is_all = (rule['Ligand'] == "All" or str(self.cb_lig['state']) == 'disabled')

            # If entire date is excluded
            if (rule['Date'] != "All" and
                    ligand_is_all and
                    rule['Cell_Line'] == "All" and
                    rule['Condition'] == "All" and
                    rule['Replicate'] == "All" and
                    rule['Row'] == "All"):

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
            if rule['Replicate'] != "All":
                target_rep = str(rule['Replicate'])
                df = df[df['Replicate'] == target_rep]

            if df.empty:
                self.log(f"   [WARNING] Rule {rule} matched 0 records.")
                continue

            # Handle specific Replicates and Rows
            for index, row_data in df.iterrows():
                result_obj = row_data['Ref_Result']
                col_idx = int(row_data['Column_Index'])

                target_rows = "ABCDEFGH"
                if rule['Row'] != "All": target_rows = rule['Row']

                for r in target_rows:
                    well_id = f"{r}{col_idx}"
                    if well_id not in result_obj.excluded_wells:
                        result_obj.excluded_wells.append(well_id)
                        count_wells += 1
                        logger.debug(f"Excluded {well_id} in {result_obj.file_name}")

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

        # Data Category Dropdown
        tk.Label(type_frame, text="Category:").grid(row=0, column=0, padx=5, pady=2, sticky="w")
        self.combo_category = ttk.Combobox(type_frame, textvariable=self.var_category, state="readonly", width=10)
        self.combo_category.grid(row=0, column=1, padx=5, pady=5, sticky="w")
        # Bind change event to update the second dropdown
        self.combo_category.bind("<<ComboboxSelected>>", self.update_subtype_options)

        # Specific Data Type
        tk.Label(type_frame, text="Subtype:").grid(row=1, column=0, padx=5, pady=2, sticky="w")
        self.combo_specific = ttk.Combobox(type_frame, textvariable=self.var_specific_type, state="readonly", width=60)
        self.combo_specific.grid(row=1, column=1, padx=5, pady=5, sticky="w")

        # --- Layout Selection ---
        # Groupy by
        tk.Label(type_frame, text="Group By:").grid(row=3, column=0, padx=5, pady=5, sticky="w")
        self.var_group_by = tk.StringVar(value="None")
        combo_group = ttk.Combobox(type_frame, textvariable=self.var_group_by, state="readonly", width=15)
        combo_group['values'] = ["None", "Cell Line", "Transfection"]
        combo_group.grid(row=3, column=1, padx=5, pady=5, sticky="w")

        # Kinetic Layout (Only applies if a Kinetic type is chosen)
        tk.Label(type_frame, text="Kinetic Layout:").grid(row=0, column=2, padx=5, pady=5, sticky="e")
        self.lb_kin_layout = tk.Listbox(type_frame, selectmode="multiple", height=6, exportselection=False)
        self.lb_kin_layout.grid(row=0, column=3, rowspan = 2, padx=3, pady=5, sticky="nsew")
        tk.Button(type_frame, text="Select All Rows", command=lambda: self.lb_kin_layout.select_set(0, tk.END)).grid(
            row=2, column=3, pady=5, sticky="nsew")
        tk.Button(type_frame, text="Deselect All Rows", command=lambda: self.lb_kin_layout.selection_clear(0, tk.END)).grid(
            row=3, column=3, pady=5, sticky="nsew")

        type_frame.columnconfigure(1, weight=1)
        type_frame.columnconfigure(3, weight=3)
        type_frame.rowconfigure(0, weight=2)
        type_frame.rowconfigure(1, weight=2)

        # --- Action Button ---
        btn_frame = tk.Frame(self.tab_plot_helper)
        btn_frame.pack(fill="x", padx=10, pady=20)

        self.btn_run_plot_helper = tk.Button(btn_frame, text="Generate Custom Export",
                                            state="disabled", command=self.run_plot_helper)
        self.btn_run_plot_helper.pack(fill="x", ipady=5)

    def update_subtype_options(self, event=None):
        """
        Triggered when Data type Category changes.
        Populates specific type dropdown based on category and availability in master_df.
        Toggles Kinetic Layout listbox visibility.
        """
        if self.master_df is None or self.master_df.empty: return

        selected_cat = self.var_category.get()  # "kinetic" or "CRC"
        subtype_map = self.data_type_map.get(selected_cat, {})

        # Filter to find valid specific options in master_df
        valid_specifics = []

        # Iterate over the sub-dictionary
        for display_name, internal_col in subtype_map.items():
            # Check if data exists for this specific type
            if internal_col in self.master_df.columns and self.master_df[internal_col].notna().any():
                valid_specifics.append(display_name)

        self.combo_specific['values'] = valid_specifics

        if valid_specifics:
            self.combo_specific.current(0)  # Select first
        else:
            self.var_specific_type.set("")

        # Toggle Kinetic Layout Listbox
        if selected_cat == "kinetic":
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
            experiment = self.master_df['Main_Plasmids'].unique()[0] if 'Main_Plasmids' in self.master_df.columns else "Experiment"
            self.lbl_data_source.config(text=f"Internal: {experiment}")

        self.btn_run_plot_helper.config(state="normal")

        available_categories = set()

        # Iterate over high-level keys ("kinetic", "CRC")
        for cat, sub_map in self.data_type_map.items():
            # Check if ANY column in this category exists in the dataframe
            for internal_col in sub_map.values():
                # And at least one value is not NA
                if internal_col in self.master_df.columns and self.master_df[internal_col].notna().any():
                    available_categories.add(cat)
                    break  # Found one valid column, so this category is valid

        sorted_cats = sorted(list(available_categories), reverse=True)
        self.combo_category['values'] = sorted_cats

        # Trigger update of specific types if categories exist
        if sorted_cats:
            self.update_subtype_options()

        # Enable kinetic layout to populate list box
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

    def run_plot_helper(self):
        """Collects GUI selections and calls the core export engine."""
        if self.master_df is None or self.master_df.empty: return

        # Get Selections
        cells = [self.lb_exp_cells.get(i) for i in self.lb_exp_cells.curselection()]
        transfections = [self.lb_exp_trans.get(i) for i in self.lb_exp_trans.curselection()]
        ligands = [self.lb_ligands.get(i) for i in self.lb_ligands.curselection()]
        rows = [self.lb_kin_layout.get(i) for i in self.lb_kin_layout.curselection()]

        category = self.var_category.get()
        specific_type = self.var_specific_type.get()

        if not category or not specific_type: return
        internal_name = self.data_type_map.get(category, {}).get(specific_type)

        if not internal_name:
            logger.error(f"Unknown data type selected: {category} - {specific_type}")
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

            # For Backward compatibility: Assign default version and path for older CSVs (< v2)
            if "NCollector_version" not in df.columns:
                df["NCollector_version"] = "< v2"
                logger.info("'NCollector_version' column missing. Assigned '< v2'.")
            if "Path" not in df.columns:
                df["Path"] = "undocumented path"
                logger.info("'Path' column missing. Assigned 'undocumented path'.")

            # For Backward compatibility with NCollector v1.0.0
            if "Empty/NoID" in df['Transfection'].values:
                logger.info("Removing legacy 'Empty/NoID' data...")
                df = df[df['Transfection'] != "Empty/NoID"].copy()

            if 'Bl_LP' not in df.columns and 'Bl_Corrected_BRET' in df.columns:
                logger.info("'Bl_LP' missing in CSV. Reconstructing from 'Bl_Corrected_BRET'...")
                # Sort to guarantee chronological order
                df_sorted = df.sort_values(by=['File_Name', 'Well_ID', 'Time_(min)'])
                # Grab the last 3 timepoints per file and well
                last_3 = df_sorted.groupby(['File_Name', 'Well_ID']).tail(3)

                # Calculate the mean of those last 3 points
                bl_lp_means = last_3.groupby(['File_Name', 'Well_ID'])['Bl_Corrected_BRET'].mean().reset_index()
                bl_lp_means.rename(columns={'Bl_Corrected_BRET': 'Bl_LP'}, inplace=True)

                # Merge the calculated values back into the main DataFrame
                df = df.merge(bl_lp_means, on=['File_Name', 'Well_ID'], how='left')

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
            self.directory = filedialog.askdirectory(title="Select a folder...")
            if not self.directory: return

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
            self.folder_path.set(f"Selected Path: {self.directory}\n\nScanning for files...")
            self.load_files_button.config(state="disabled")
            self.main_gi.update()

            # --- Scan new directory for xlsx or xlsm---
            for root, dirs, files in os.walk(self.directory):
                if any(f.endswith((".xlsx", ".xlsm")) for f in files):
                    self.subfolder_paths_with_files.append(root) # Save paths of files

            # Update GUI
            if self.subfolder_paths_with_files:
                count = len(self.subfolder_paths_with_files)
                folder_names = [os.path.basename(path) for path in self.subfolder_paths_with_files]
                folder_names_string = "\n ".join(folder_names)
                self.folder_path.set(
                    f"Selected Path: {self.directory}\n\n Found following subfolders with xlsx/xlsm files:\n {folder_names_string}")
                self.load_files_button.config(state="normal")
                logger.info(f"Found {count} folders: \n {folder_names_string}")
            else:
                self.folder_path.set(f"Error: No .xlsx or .xlsm files found in {self.directory} or any subfolder.")
                self.load_files_button.config(state="disabled")
                logger.warning(f"No .xlsx or .xlsm files found starting from: {self.directory}")

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
        logger.info(f"Keeping {len(self.experiment)} folders matching main plasmids: {selected_key}")

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
            logger.debug(f"Data Summary:\n{summary}")
            self.refresh_filter_options()
            self.update_summary_table()
            self.refresh_plot_helper_options()
            return self.master_index
        else:
            self.master_index = pd.DataFrame()
            self.update_summary_table()
            logger.warning("No valid data found")
            return self.master_index

    def collect_files(self):
        """
        1. LOAD FILES
        Reads sheet names of all xlsx and xlsm files to identify and separate protocol and result analysis files.
        Validation of correct protocol to analysis files is done via date of measurement in the folder name.
        """
        if not self.subfolder_paths_with_files:
            logger.warning("No folders to analyze.")
            return

        # Create config
        try:
            is_labeling = self.var_labeling_is_checked.get()
            lum_threshold = self.var_lum_threshold.get()
        except tk.TclError:
            is_labeling = False
            lum_threshold = 100

        if is_labeling:
            # Quadruplicates
            selected_layout = [
                range(1, 5),  # Block 1: Cols 1-4
                range(5, 9),  # Block 2: Cols 5-8
                range(9, 13)   # Block 3: Cols 9-12
            ]
        else:
            # Standard triplicate Layout
            selected_layout = [
                range(1, 4),  # Block 1: Cols 1-3
                range(4, 7),  # Block 2: Cols 4-6
                range(7, 10),  # Block 3: Cols 7-9
                range(10, 13)  # Block 4: Cols 10-12
            ]

        self.current_config = ProcessingConfig(
            lum_threshold=lum_threshold,
            labeling_correction=is_labeling,
            plate_layout=selected_layout
        )

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
            logger.debug(f"Skipped files: {folder_data.skipped_files}")

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

        # Get current configuration
        current_config = self.current_config
        if current_config.labeling_correction:
            self.log("[PROCESSING CONFIG]   Labeling correction is applied")

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

                empty_cols = [str(idx) for idx, meta in result.column_metadata.items() if
                              meta.condition_name == "Empty"]
                if empty_cols:
                    self.log(f"   [EMPTY WELLS] Empty columns {', '.join(empty_cols)} detected in {result.file_name}")

                # Gather warnings
                all_detected_warnings.extend(result.low_lum_warnings)
                all_detected_warnings.extend(result.vehicle_warnings)

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

                # --- PREPARE KINETIC DATA ---
                # Use pandas melt function to prepare each df from wide to long format
                def melt_df(df, val_name, time_vec):
                    if df is None or df.empty: return pd.DataFrame()

                    df_work = df.copy()
                    # Check lengths
                    if len(df_work) != len(time_vec):
                        logger.warning(
                            f"Length mismatch in {res.file_name}: Data {len(df_work)} vs Time {len(time_vec)}")
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
                    t_vec = range(len(res.raw_bret_ratio_cleaned))

                # Ignore time col in raw bret df
                raw_clean = res.raw_bret_ratio_cleaned.drop(columns=["Time (min)"], errors='ignore')

                df_og = melt_df(res.raw_bret_ratio_df.drop(columns=["Time (min)"], errors='ignore'),
                                "OG_BRET_ratio", t_vec)
                df_raw = melt_df(raw_clean, "Raw_BRET_kinetic", t_vec)
                df_lab = melt_df(res.labeling_corr_kinetic, "Lab_BRET_kinetic", t_vec)
                df_bl = melt_df(res.bl_corr_kinetic, "Bl_Corrected_BRET", t_vec)
                df_norm = melt_df(res.kinetic_df, "Veh_Norm_Kinetic", t_vec)

                # Merge on [Time_(min), Well_ID]
                merge_on = [df_og.columns[0], "Well_ID"]

                merged_df = df_og.merge(df_raw, on=merge_on, how="left") \
                    .merge(df_lab, on=merge_on, how="left") \
                    .merge(df_bl, on=merge_on, how="left") \
                    .merge(df_norm, on=merge_on, how="left")

                # --- MAP DATA FOR CRC ---
                # Last 3x points/AUC is 1 value per well, map data to Well_ID
                # Last 3 time points (lp)
                raw_bret_map = res.raw_bret_points_df.iloc[0].to_dict() if res.raw_bret_points_df is not None else {}
                lp_lab_map = res.labeling_corr_lp_df.iloc[0].to_dict() if res.labeling_corr_lp_df is not None else {}
                lp_bl_map = res.bl_corr_lp_df.iloc[0].to_dict() if res.bl_corr_lp_df is not None else {}
                lp_norm_map = res.lp_df.iloc[0].to_dict() if res.lp_df is not None else {}

                # AUC
                auc_lab_map = res.labeling_corr_auc_df.iloc[0].to_dict() if res.labeling_corr_auc_df is not None else {}
                auc_bl_map = res.bl_corr_auc_df.iloc[0].to_dict() if res.bl_corr_auc_df is not None else {}
                auc_norm_map = res.auc_df.iloc[0].to_dict() if res.auc_df is not None else {}

                merged_df['Raw_BRET_CRC'] = merged_df['Well_ID'].map(raw_bret_map)
                merged_df['Lab_LP'] = merged_df['Well_ID'].map(lp_lab_map)
                merged_df['Bl_LP'] = merged_df['Well_ID'].map(lp_bl_map)
                merged_df['Veh_Norm_LP'] = merged_df['Well_ID'].map(lp_norm_map)

                merged_df['Lab_AUC'] = merged_df['Well_ID'].map(auc_lab_map)
                merged_df['Bl_AUC'] = merged_df['Well_ID'].map(auc_bl_map)
                merged_df['Veh_Norm_AUC'] = merged_df['Well_ID'].map(auc_norm_map)

                # --- PREPARE MEAN KINETIC AND AUC MAPPING ---
                well_to_mean_map = {}
                well_to_lp_mean_map = {}
                well_to_auc_mean_map = {}
                for col_idx, meta in res.column_metadata.items():
                    col_str = str(col_idx)
                    for row_char in "ABCDEFGH":
                        well_id = f"{row_char}{col_str}"

                        # Construct Key: "Condition|Cell|Ligand|Row"
                        mean_key = f"{meta.condition_name}|{meta.cell_line}|{meta.ligand_identity}|{row_char}"

                        # Grab Kinetic Mean Series
                        if res.kinetic_mean_df is not None and mean_key in res.kinetic_mean_df.columns:
                            well_to_mean_map[well_id] = res.kinetic_mean_df[mean_key].tolist()

                        # Grab Last points Mean Value
                        if res.lp_mean_df is not None and mean_key in res.lp_mean_df.columns:
                            well_to_lp_mean_map[well_id] = res.lp_mean_df[mean_key].iloc[0]

                        # Grab AUC Mean Value
                        if res.auc_mean_df is not None and mean_key in res.auc_mean_df.columns:
                            well_to_auc_mean_map[well_id] = res.auc_mean_df[mean_key].iloc[0]

                merged_df['LP_Mean'] = merged_df['Well_ID'].map(well_to_lp_mean_map)
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
                merged_df["NCollector_version"] = self.version
                merged_df["Path"] = self.directory
                merged_df["File_Name"] = res.file_name
                merged_df["Date"] = res.measurement_date
                merged_df["Main_Plasmids"] = main_plasmids

                # Get the exclusion text (handle empty case)
                exclusion_text = self.rule_history_text if self.rule_history_text else "None"
                # Clean newlines for CSV compatibility
                exclusion_text_clean = exclusion_text.replace("\n", " | ")
                merged_df["Applied_Exclusions"] = exclusion_text_clean

                # Meta Lookups (Optimization: Build dicts once per file)
                meta_lookups = {'Transfection': {},
                                'Cell_Line': {},
                                'Ligand': {},
                                'Ligand_Conc': {},
                                'Plate_Row': {},
                                'Replicate': {}}

                for well_id in merged_df['Well_ID'].unique():
                    try:
                        c_idx = int(well_id[1:])
                        row_char = well_id[0]
                        meta = res.column_metadata.get(c_idx)
                        if meta:
                            meta_lookups['Transfection'][well_id] = meta.condition_name
                            meta_lookups['Cell_Line'][well_id] = meta.cell_line
                            meta_lookups['Ligand'][well_id] = meta.ligand_identity
                            meta_lookups['Ligand_Conc'][well_id] = meta.ligand_conc.get(row_char, 0.0)
                            meta_lookups['Plate_Row'][well_id] = row_char
                            meta_lookups['Replicate'][well_id] = meta.replicate
                    except: pass

                merged_df['Transfection'] = merged_df['Well_ID'].map(meta_lookups['Transfection'])
                merged_df['Cell_Line'] = merged_df['Well_ID'].map(meta_lookups['Cell_Line'])
                merged_df['Ligand'] = merged_df['Well_ID'].map(meta_lookups['Ligand'])
                merged_df['Ligand_Conc'] = merged_df['Well_ID'].map(meta_lookups['Ligand_Conc'])
                merged_df['Plate_Row'] = merged_df['Well_ID'].map(meta_lookups['Plate_Row'])
                merged_df['Replicate'] = merged_df['Well_ID'].map(meta_lookups['Replicate'])

        if not all_files_data:
            return pd.DataFrame()

        # Combine all files
        master_df = pd.concat(all_files_data, ignore_index=True)
        # Cleanup columns
        cols_order = [
            "NCollector_version", "Path",
            "File_Name", "Date", "Main_Plasmids", "Applied_Exclusions",
            "Transfection", "Cell_Line", "Ligand",
            "Ligand_Conc", "Plate_Row", "Replicate", "Well_ID", "Time_(min)",
            "Raw_BRET_kinetic", "Lab_BRET_kinetic", "Bl_Corrected_BRET", "Veh_Norm_Kinetic", "Kinetic_Mean",
            "Raw_BRET_CRC", "Lab_LP", "Bl_LP", "Veh_Norm_LP", "LP_Mean",
            "Lab_AUC", "Bl_AUC", "Veh_Norm_AUC", "AUC_Mean"
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
                logger.error("Export failed: Filter resulted in no data.")
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

            # Check whether labeling control was applied
            is_labeling = True if "labeling control" in df_subset["Replicate"].values else False

            # Definition which is kinetic and what is CRC
            kinetic_types = list(self.data_type_map["kinetic"].values())
            crc = list(self.data_type_map["CRC"].values())

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
                        # Only keep the labeling control column if raw BRET ratio is exported
                        if dtype in ["Raw_BRET_kinetic", "Raw_BRET_CRC"] and is_labeling:
                            drop_labeling_col = False
                        else:
                            drop_labeling_col = True

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
                            kin_pivot = create_clean_pivot(df_kin, "Time_(min)",
                                                           dtype, "Mean" in dtype,
                                                           drop_labeling_col)
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
                            crc_pivot = create_clean_pivot(df_crc, "Plate_Row",
                                                           dtype, "Mean" in dtype,
                                                           drop_labeling_col)

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
            logger.info("Log is empty, nothing to save.")
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
                logger.info(f"Log saved to: {file_path}")
            except Exception as e:
                logger.error(f"Error saving log: {e}")

# --- Main Execution Block ---
if __name__ == "__main__":
    # Configure logging: DEBUG to console, adjust level as needed
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S"
    )

    # Create the main window
    root = tk.Tk()

    # Create an instance of the application
    app = NCollectorApp(root)

    # Start the Tkinter event loop
    root.mainloop()