import os
import logging
import tkinter as tk
from tkinter import filedialog, ttk, messagebox
import numpy as np
import pandas as pd
from datetime import datetime, date

from models import (MeasurementFolder, ProcessingConfig, APP_VERSION, MASTER_COLUMNS,
                    DATA_TYPE_MAP, build_plate_layout, ENRICHABLE_COLS, LEGACY_COLUMN_DEFAULTS,
                    SINGLE_CONC_CATEGORIES, REQUIRES_GROUP_BY)
from parsing import scan_and_load_folders
from processing import (process_bret_measurement, calculate_relative_time, map_plate_metadata,
                        recompute_master_after_exclusion, reconstruct_file_inputs,
                        check_luminescence, check_vehicle_wells, coerce_bool, parse_date_series)
from mapping import infer_ligand_info_from_master
from export import (apply_export_filters, build_row_info, generate_header_key,
                    create_clean_pivot, create_bargraph_table, create_heatmap_table, filter_by_conc,
                    ensure_master_csv_schema)
from restore import (list_active_exclusions, restore_rule, restore_wells,
                     build_resolve_ctx, well_token)
from dialogs import (ask_user_parameter, ask_ligand_choice, ask_ligand_layout,
                     ask_filename_collision, warn_and_abort, ask_main_plasmids_selection)
from plasmid_selection import (group_by_main_plasmids, resolve_selection)
import merge
import crc_window

logger = logging.getLogger("NCollector")

# --- Lightweight tooltip helper --- #

class _Tooltip:
    """Minimal hover tooltip for a Tk widget (used to explain disabled buttons)."""
    def __init__(self, widget):
        self.widget = widget
        self.tip = None
        self.text = ""
        widget.bind("<Enter>", self._show)
        widget.bind("<Leave>", self._hide)

    def set_text(self, text):
        self.text = text or ""

    def _show(self, _e=None):
        if self.tip or not self.text:
            return
        try:
            x = self.widget.winfo_rootx() + 20
            y = self.widget.winfo_rooty() + self.widget.winfo_height() + 5
            self.tip = tk.Toplevel(self.widget)
            self.tip.wm_overrideredirect(True)
            self.tip.wm_geometry(f"+{x}+{y}")
            tk.Label(self.tip, text=self.text, background="#ffffe0",
                     relief="solid", borderwidth=1, font=("Arial", 9),
                     justify="left").pack(ipadx=3, ipady=2)
        except tk.TclError:
            self.tip = None

    def _hide(self, _e=None):
        if self.tip:
            self.tip.destroy()
            self.tip = None


# --- Main Application --- #

class NCollectorApp:
    def __init__(self, main_window):
        self.main_gi = main_window
        main_window.title(APP_VERSION)
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
        self.conc_row_lookup = {}  # Maps display strings to filter criteria for conc layout

        # --- GUI Variables ---
        self.var_labeling_is_checked = tk.BooleanVar(value=False)
        self.var_lum_threshold = tk.IntVar(value=100)
        self.var_vehicle_threshold = tk.DoubleVar(value=0.2)

        # --- GUI Widgets (Initialised to None) ---
        self.log_window = None
        self.log_text = None
        self.notebook = None
        # Tab 1
        self.tab_import = None
        self.subfolders_label = None
        self.load_files_button = None
        self.main_plasmids_label = None
        self.summary_tree = None
        self.btn_export_master = None
        self.btn_export_excel = None
        self.btn_rerun_lum = None
        self.btn_rerun_vehicle = None
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
        # Tab 2 — Exclude/Revert
        self.var_excl_mode = None          # StringVar in {"exclude", "revert"} (default "exclude")
        self.rb_exclude = None
        self.rb_revert = None
        self.btn_apply_excl = None         # "Apply exclusions and re-calculate" (exclude mode only)
        self.btn_revert_excl = None        # "Revert exclusions and re-calculate" (revert mode only)
        self.lb_active_excl_tab2 = None    # Tab 2 multi-select active-exclusions box (revert mode only)
        # Shared "Active Exclusions" state (both tabs are fed by refresh_active_exclusions)
        self.lb_active_excl_tab1 = None    # Tab 1 display-only active-exclusions box (replaces lbl_rules_summary)
        self._active_excl_entries = []     # Parallel list: visible row index -> entry dict (both boxes share order)
        # Tab 3
        self.tab_plot_helper = None
        self.lbl_data_source = None
        self.lbl_data_source_tab1 = None  # Tab 1 mirror of the data-source label
        self.lbl_csv_source = None         # Tab 1 CSV-only label (imported CSV name)
        self.lb_ligands = None
        self.lb_exp_cells = None
        self.lb_exp_trans = None
        self.combo_category = None
        self.combo_specific = None
        self.lb_conc_layout = None
        self.btn_run_plot_helper = None
        # Tab 4 — Merge
        self.tab_merge = None
        # List of source entries:
        #   {"label": str, "kind": "master"|"folder", "path": str,
        #    "df": pd.DataFrame, "n_files": int, "n_rows": int}
        self.merge_sources: list[dict] = []
        self.merge_tree = None
        self.merge_summary_tree = None
        self.merge_main_plasmids_label = None
        self.btn_merge_run = None
        self.merge_main_plasmids_choice = None
        self.crc_window = None # exclude/restore refreshes this window live

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

        self.tab_merge = tk.Frame(self.notebook)
        self.notebook.add(self.tab_merge, text="Merge")

        self.setup_import_tab()
        self.setup_exclusion_tab()
        self.setup_plot_helper_tab()
        self.setup_merge_tab()

    def setup_import_tab(self):
        # --- Three side-by-side boxes: Loading Specs | Folder Import | Master Import ---
        io_frame = tk.Frame(self.tab_import)
        io_frame.pack(fill="x", padx=10, pady=10)
        io_frame.columnconfigure(0, weight=0)  # loading specs (natural width, far left)
        io_frame.columnconfigure(1, weight=1)  # folder box (expands)
        io_frame.columnconfigure(2, weight=0)  # master box (natural width)

        # === BOX 1 (far left): Loading Specs ===
        load_settings_frame = tk.LabelFrame(io_frame, text="Loading Specs")
        load_settings_frame.grid(row=0, column=0, sticky="nsew", padx=(0, 5))
        load_settings_frame.columnconfigure(2, weight=1)  # let buttons stretch

        # Row 0 — Labeling correction (no button on this row)
        tk.Label(load_settings_frame, text="Labeling correction").grid(
            row=0, column=0, sticky="w", padx=(8, 5), pady=(6, 4))
        labeling_chk = tk.Checkbutton(load_settings_frame,
                                      variable=self.var_labeling_is_checked,
                                      command=self.on_checkbox_toggle)
        labeling_chk.grid(row=0, column=1, sticky="w", pady=(6, 4))

        # Row 1 — Luminescence threshold + its re-run button (same line)
        tk.Label(load_settings_frame, text="Lum. Threshold:").grid(
            row=1, column=0, sticky="w", padx=(8, 5), pady=4)
        lum_thresh_entry = tk.Entry(load_settings_frame,
                                    textvariable=self.var_lum_threshold, width=10)
        lum_thresh_entry.grid(row=1, column=1, sticky="w", pady=4)
        self.btn_rerun_lum = tk.Button(load_settings_frame, text="Re-run Lum Check",
                                       state="disabled", command=self.rerun_lum_check)
        self.btn_rerun_lum.grid(row=1, column=2, sticky="ew", padx=(10, 8), pady=4)

        # Row 2 — Vehicle warning threshold + its re-run button (same line).
        tk.Label(load_settings_frame, text="Veh. Threshold:").grid(
            row=2, column=0, sticky="w", padx=(8, 5), pady=(4, 6))
        veh_thresh_entry = tk.Entry(load_settings_frame,
                                    textvariable=self.var_vehicle_threshold, width=10)
        veh_thresh_entry.grid(row=2, column=1, sticky="w", pady=(4, 6))
        self.btn_rerun_vehicle = tk.Button(load_settings_frame, text="Re-run Vehicle Check",
                                           state="disabled", command=self.rerun_vehicle_check)
        self.btn_rerun_vehicle.grid(row=2, column=2, sticky="ew", padx=(10, 8), pady=(4, 6))

        # === BOX 2 (center): Import from experiment folder ===
        box_folder = tk.LabelFrame(io_frame, text="Import from Experiment Folder")
        box_folder.grid(row=0, column=1, sticky="nsew", padx=5)

        # Select Folder button
        tk.Button(box_folder, text="Select folder containing results of experiment",
                  command=self.select_folder).pack(fill="x", anchor="w", padx=8, pady=(8, 4))

        # Found subfolders display (path itself is logged, not shown here)
        self.subfolders_label = tk.Label(box_folder, text="No folder selected.", wraplength=320,
                                         justify="left", font=('Arial', 10))
        self.subfolders_label.pack(anchor="nw", padx=8, pady=0)

        # Load Files button
        self.load_files_button = tk.Button(box_folder, text="Load Files", state="disabled",
                                           command=self.collect_files)
        self.load_files_button.pack(fill="x", padx=8, pady=5)

        # === BOX 3 (right): Import Master CSV ===
        box_master = tk.LabelFrame(io_frame, text="Import Master CSV")
        box_master.grid(row=0, column=2, sticky="nsew", padx=(5, 0))

        # Import Master CSV button — same handler as Tab 3
        tk.Button(box_master, text="Import Master CSV",
                  command=self.import_master_csv).pack(anchor="w", padx=8, pady=(8, 4))

        # CSV-specific label: shows the imported CSV file name (CSV imports only,
        # not experiment-folder loads). "No CSV loaded" when no master CSV imported.
        self.lbl_csv_source = tk.Label(box_master, text="No CSV loaded",
                                       justify="left", wraplength=180)
        self.lbl_csv_source.pack(anchor="w", padx=8, pady=(0, 8))

        # Collected Ns frame
        loaded_data_frame = tk.LabelFrame(self.tab_import, text="Loaded Data")
        loaded_data_frame.pack(fill="both", expand=True,  padx=10, pady=5)

        # Data Source (mirrors the Tab 3 data-source label), above Main Plasmids
        ds_frame = tk.Frame(loaded_data_frame)
        ds_frame.pack(anchor="w", padx=8, pady=(6, 0))
        tk.Label(ds_frame, text="Data Source", font=("Arial", 9, "bold")).pack(side="left")
        self.lbl_data_source_tab1 = tk.Label(ds_frame, text="No Data Loaded", justify="left")
        self.lbl_data_source_tab1.pack(side="left", padx=(6, 0))

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

        # Active exclusions (display-only)
        # the SAME shared "Active Exclusions" view as Tab 2, fed by refresh_active_exclusions()
        # off list_active_exclusions(master_df)
        rules_frame = tk.LabelFrame(self.tab_import, text="Active Exclusions")
        rules_frame.pack(fill="both", expand=True, padx=10, pady=5)
        active_box1 = tk.Frame(rules_frame)
        active_box1.pack(fill="both", expand=True, padx=5, pady=5)
        active_scroll1 = ttk.Scrollbar(active_box1, orient="vertical")
        active_scroll1.pack(side="right", fill="y")
        # Display-only: selectmode "none"
        self.lb_active_excl_tab1 = tk.Listbox(active_box1, height=4, activestyle="none",
                                              selectmode="none", exportselection=False,
                                              yscrollcommand=active_scroll1.set)
        self.lb_active_excl_tab1.pack(side="left", fill="both", expand=True)
        active_scroll1.config(command=self.lb_active_excl_tab1.yview)
        # click/drag is specifically silenced (display only)
        for _seq in ("<Button-1>", "<B1-Motion>", "<Double-Button-1>"):
            self.lb_active_excl_tab1.bind(_seq, lambda e: "break")

        # --- EXPORT SECTION ---
        export_frame = tk.LabelFrame(self.tab_import, text="Export Options")
        export_frame.pack(fill="x", padx=10, pady=10)
        self.btn_export_master = tk.Button(export_frame, text="Export Master CSV", state="disabled",
                                           command=self.export_master_csv)
        self.btn_export_master.pack(side="left", fill="x", expand=True, padx=5, pady=10)
        self.btn_export_excel = tk.Button(export_frame, text="Export Excel Report (Default)", state="disabled",
                                          command=self.export_excel_report)
        self.btn_export_excel.pack(side="left", fill="x", expand=True, padx=5, pady=10)

    def _set_data_source(self, text):
        """Update both the Tab 3 and Tab 1 data-source labels in sync."""
        if self.lbl_data_source is not None:
            self.lbl_data_source.config(text=text)
        if self.lbl_data_source_tab1 is not None:
            self.lbl_data_source_tab1.config(text=text)

    def _set_csv_label(self, name):
        """Tab 1 CSV-only label: show the imported CSV file name, or 'No CSV loaded'.
        Reflects master-CSV imports only — never experiment-folder loads."""
        if self.lbl_csv_source is not None:
            self.lbl_csv_source.config(text=name if name else "No CSV loaded")

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

        # Live concentration-response plot in separate window
        plot_frame = tk.Frame(self.tab_select)
        plot_frame.pack(fill="x", pady=(10, 0))
        tk.Button(plot_frame, text="Plot concentration-response curves",
                  command=self.open_crc_window).pack(fill="x")

        # Mode toggle (Exclude vs Revert)
        # Mutually-exclusive, exactly one selected, default "exclude"
        self.var_excl_mode = tk.StringVar(value="exclude")
        mode_frame = tk.LabelFrame(self.tab_select, text="Mode")
        mode_frame.pack(fill="x", pady=(5, 0), padx=5)
        self.rb_exclude = ttk.Radiobutton(mode_frame, text="Exclude Data",
                                           variable=self.var_excl_mode, value="exclude",
                                           command=self._on_excl_mode_change)
        self.rb_exclude.pack(side="left", padx=10, pady=5)
        self.rb_revert = ttk.Radiobutton(mode_frame, text="Revert Exclusions",
                                          variable=self.var_excl_mode, value="revert",
                                          command=self._on_excl_mode_change, state="disabled")
        self.rb_revert.pack(side="left", padx=10, pady=5)

        # Frame for dropdowns
        filter_frame = tk.LabelFrame(self.tab_select, text="Exclude & Revert Data")
        filter_frame.pack(fill = "x", pady=10, padx=5)

        # Variables
        self.var_lig = tk.StringVar(value="")
        self.var_date = tk.StringVar(value="All")
        self.var_cell = tk.StringVar(value="All")
        self.var_cond = tk.StringVar(value="All")
        self.var_repl = tk.StringVar(value="All")
        self.var_row = tk.StringVar(value="All")
        # Scope (Main_Plasmids + File_Name) selectors shown only on post-merge masters (see_refresh_scope_selectors)
        self.var_main = tk.StringVar(value="All")
        self.var_file = tk.StringVar(value="All")

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

        # Scope filters (Main_Plasmids + File_Name) shown only on post-merge masters
        # (populated by _refresh_scope_selectors)
        scope_frame = tk.Frame(filter_frame)
        scope_frame.grid(row=2, column=0, columnspan=8, pady=5, sticky="w")

        self.lbl_main = tk.Label(scope_frame, text="Main Plasmids:")
        self.lbl_main.grid(row=0, column=0, padx=5)
        self.cb_main = ttk.Combobox(scope_frame, textvariable=self.var_main,
                                    state="readonly", width=18)
        self.cb_main.grid(row=0, column=1, padx=5)
        self.cb_main['values'] = ["All"]
        self.cb_main.bind("<<ComboboxSelected>>",
                          lambda e: self.update_dropdown_options("Main_Plasmids"))

        self.lbl_file = tk.Label(scope_frame, text="Source File:")
        self.lbl_file.grid(row=0, column=2, padx=5)
        self.cb_file = ttk.Combobox(scope_frame, textvariable=self.var_file,
                                    state="readonly", width=22)
        self.cb_file.grid(row=0, column=3, padx=5)
        self.cb_file['values'] = ["All"]
        self.cb_file.bind("<<ComboboxSelected>>",
                          lambda e: self.update_dropdown_options("File_Name"))

        # Hidden until a master proves them meaningful (grid_remove preserves placement).
        self.lbl_main.grid_remove(); self.cb_main.grid_remove()
        self.lbl_file.grid_remove(); self.cb_file.grid_remove()

        # Buttons
        btn_frame = tk.Frame(filter_frame)
        btn_frame.grid(row=3, column=0, columnspan=6, pady=10)

        tk.Button(btn_frame, text="Add Rule to List", command=self.add_exclusion_rule).pack(side="left", padx=5)
        tk.Button(btn_frame, text="Clear List", command=self.clear_exclusion_list).pack(side="left", padx=5)

        # Listbox for Pending Exclusions
        list_frame = tk.LabelFrame(self.tab_select, text="Pending Exclusions (Will be removed upon Apply)")
        list_frame.pack(fill="both", expand=True, padx=10, pady=5)

        self.lb_exclusions = tk.Listbox(list_frame, height=8)
        self.lb_exclusions.pack(fill="both", expand=True, padx=5, pady=5)

        # Apply Exclusions is enabled only in EXCLUDE mode; Revert only in REVERT mode.
        action_frame = tk.Frame(self.tab_select)
        action_frame.pack(pady=10)
        self.btn_apply_excl = tk.Button(action_frame, text="Apply exclusions and re-calculate",
                                        command=self.apply_exclusions)
        self.btn_apply_excl.pack(side="left", padx=5, ipadx=10)
        self.btn_revert_excl = tk.Button(action_frame, text="Revert exclusions and re-calculate",
                                         command=self.revert_exclusions, state="disabled")
        self.btn_revert_excl.pack(side="left", padx=5, ipadx=10)

        # Active Exclusions box (Tab 2 copy)
        # Multi-select; selectable ONLY in revert mode
        active_frame2 = tk.LabelFrame(self.tab_select, text="Active Exclusions")
        active_frame2.pack(fill="both", expand=True, padx=10, pady=5)
        active_box2 = tk.Frame(active_frame2)
        active_box2.pack(fill="both", expand=True, padx=5, pady=5)
        active_scroll2 = ttk.Scrollbar(active_box2, orient="vertical")
        active_scroll2.pack(side="right", fill="y")
        self.lb_active_excl_tab2 = tk.Listbox(active_box2, height=5, selectmode="extended",
                                              exportselection=False,
                                              yscrollcommand=active_scroll2.set)
        self.lb_active_excl_tab2.pack(side="left", fill="both", expand=True)
        active_scroll2.config(command=self.lb_active_excl_tab2.yview)
        # Selection guard: in exclude mode, swallow click/drag so nothing can be selected.
        self.lb_active_excl_tab2.bind("<Button-1>", self._tab2_active_select_guard)
        self.lb_active_excl_tab2.bind("<B1-Motion>", self._tab2_active_select_guard)

        # Initialise button/selectability state and prime both active-exclusion boxes.
        self._update_excl_action_buttons()
        self._apply_tab2_selectability()

    # --- Exclude/Revert mode ----------------------------------- #

    def _on_excl_mode_change(self):
        """
        Mode change handler (TASK 1). Bound to both radiobuttons' command and also invoked
        when the Revert option is force-disabled. On every change:
          (a) repopulate the comboboxes from the correct subset (refresh_filter_options
              -> mode-aware _option_source_df),
          (b) CLEAR the pending list (a pending exclude rule is meaningless in revert mode
              and vice versa),
          (c) toggle the two action buttons,
          (d) toggle the Tab-2 Active-Exclusions box selectability.
        """
        self.refresh_filter_options()          # (a)
        self.clear_exclusion_list()             # (b)
        self._update_excl_action_buttons()      # (c)
        self._apply_tab2_selectability()        # (d)

    def _tab2_active_select_guard(self, _event):
        """Swallow clicks/drags on the Tab-2 Active-Exclusions box unless in revert mode."""
        if self.var_excl_mode is None or self.var_excl_mode.get() != "revert":
            return "break"
        return None

    def _apply_tab2_selectability(self):
        """Clear any lingering Tab-2 selection when not in revert mode (selection is only
        meaningful while reverting; the click guard prevents new selections)."""
        if self.lb_active_excl_tab2 is None:
            return
        if self.var_excl_mode is None or self.var_excl_mode.get() != "revert":
            self.lb_active_excl_tab2.selection_clear(0, tk.END)

    def _update_excl_action_buttons(self):
        """Apply enabled only in EXCLUDE mode, Revert only in REVERT mode (TASK 4).
        Both additionally require data to be present."""
        mode = self.var_excl_mode.get() if self.var_excl_mode is not None else "exclude"
        has_data = self.master_df is not None and not self.master_df.empty
        if self.btn_apply_excl is not None:
            self.btn_apply_excl.config(
                state="normal" if (mode == "exclude" and has_data) else "disabled")
        if self.btn_revert_excl is not None:
            self.btn_revert_excl.config(
                state="normal" if (mode == "revert" and has_data) else "disabled")

    def _update_revert_mode_state(self, has_entries):
        """
        Enable/disable the "Revert Exclusions" radiobutton (TASK 1). Greyed whenever there
        are no active exclusions. If revert was the active mode and it just became empty,
        snap back to exclude mode (and run the mode-change side effects).
        """
        if self.rb_revert is None:
            return
        if has_entries:
            self.rb_revert.config(state="normal")
        else:
            self.rb_revert.config(state="disabled")
            if self.var_excl_mode is not None and self.var_excl_mode.get() == "revert":
                # Programmatic var.set() does NOT fire the radiobutton command, so run the
                # mode-change side effects explicitly.
                self.var_excl_mode.set("exclude")
                self._on_excl_mode_change()

    @staticmethod
    def _display_label(label):
        """
        Human-friendly text for an active-exclusion entry. For restore-split tokens
        (emitted by restore.py when a rule is partially reverted) the trailing
        "- value: n/a" carries no information, so it is dropped for display. The original
        label is left untouched in self._active_excl_entries for the revert dispatch.
        """
        s = str(label)
        if "[RESTORE-SPLIT]" in s and " - value:" in s:
            s = s.split(" - value:")[0].rstrip()
        return s

    def refresh_active_exclusions(self):
        """
        Shared refresh for the Active-Exclusions view (TASK 3). Calls
        list_active_exclusions(master_df) ONCE and repopulates BOTH listboxes (Tab 1
        display-only + Tab 2 multi-select), keeping self._active_excl_entries as the
        parallel visible-row-index -> entry-dict map (identical order in both boxes). Also
        re-evaluates the Revert-mode enable state and re-toggles the action buttons +
        Tab-2 selectability.

        This is the single method BOTH tabs' events call after anything that mutates
        exclusions or swaps the data source.
        """
        entries = []
        if self.master_df is not None and not self.master_df.empty:
            try:
                entries = list_active_exclusions(self.master_df)
            except Exception as e:
                logger.warning(f"Could not list active exclusions: {e}")
                entries = []
        self._active_excl_entries = entries

        # Compact row text: "<label>  (N wells)" or "(N wells; M not revertable)".
        # The displayed label is cleaned for restore-split tokens (drops the meaningless
        # "- value: n/a" tail); the REAL label is kept in self._active_excl_entries for the
        # dispatch in revert_exclusions, so this only affects what the user sees.
        rows = []
        for e in entries:
            n = len(e.get("wells", []))
            m = len(e.get("non_restorable", []))
            suffix = f"({n} wells; {m} not revertable)" if m else f"({n} wells)"
            rows.append(f"{self._display_label(e['label'])}  {suffix}")

        # Tab 1 (display-only)
        if self.lb_active_excl_tab1 is not None:
            self.lb_active_excl_tab1.delete(0, tk.END)
            if rows:
                for r in rows:
                    self.lb_active_excl_tab1.insert(tk.END, r)
            else:
                self.lb_active_excl_tab1.insert(tk.END, "No exclusion rules applied")

        # Tab 2 (multi-select; selectable only in revert mode)
        if self.lb_active_excl_tab2 is not None:
            self.lb_active_excl_tab2.delete(0, tk.END)
            for r in rows:
                self.lb_active_excl_tab2.insert(tk.END, r)

        # Revert-mode availability + button/selectability state under the current mode.
        self._update_revert_mode_state(bool(entries))
        self._update_excl_action_buttons()
        self._apply_tab2_selectability()

    def _collect_revert_inputs(self):
        """
        Gather the two revert input (exclusion rules or dropwdown selected) sources separately,
        because they mean different things and must dispatch differently:

          (A) WHOLE-RULE selection — the rows multi-selected in the Tab-2 Active-Exclusions
              box. Returned as a list of entry LABELS so each can be reverted as a unit via
              restore_rule, which already leaves wells still claimed by ANOTHER active rule
              excluded (and that other rule intact).

          (B) PER-WELL criteria — the pending combobox rules. Resolved via the SHARED resolver
              (_resolve_rule_to_wells) and intersected with the currently-excluded set, then
              released from whichever rule(s) own them (releasing a shared well does require
              releasing it from every claiming rule — that is the intended well-level semantic).

        Returns (selected_labels: list[str], pending_target: set[(file, well)]).
        """
        selected_labels = []
        if self.lb_active_excl_tab2 is not None:
            for i in self.lb_active_excl_tab2.curselection():
                if 0 <= i < len(self._active_excl_entries):
                    selected_labels.append(self._active_excl_entries[i]["label"])

        pending_target = set()
        if self.pending_exclusions:
            df = self.master_df
            norm_dates = parse_date_series(df['Date'], context="_collect_revert_inputs").dt.strftime('%d.%m.%y')
            excl = (df['Is_Excluded'].map(coerce_bool)
                    if 'Is_Excluded' in df.columns else pd.Series(False, index=df.index))
            excluded_set = set(zip(df.loc[excl, 'File_Name'].astype(str),
                                   df.loc[excl, 'Well_ID'].astype(str)))
            for rule in self.pending_exclusions:
                wells = self._resolve_rule_to_wells(rule, norm_dates)
                pending_target |= (wells & excluded_set)

        return selected_labels, pending_target

    def revert_exclusions(self):
        """
        Revert dispatch, via restore.py's label-based API. Two input sources are
        handled with DIFFERENT semantics (see _collect_revert_inputs):

          (A) Whole-rule selections -> restore_rule per selected entry. restore_rule already
              leaves wells still claimed by ANOTHER active rule excluded and leaves that
              other rule intact, so reverting one rule never silently drops an overlapping
              rule. These wells are NOT cross-dispatched into other rules.

          (B) Pending combobox criteria -> a per-well target released from whichever rule(s)
              own each well (full coverage -> restore_rule, partial -> restore_wells),
              evaluated AFTER the whole-rule reverts so it sees the updated blob.

        ALL restore decisions live in restore.py. Non-restorable / blocked wells are logged.

        Efficiency: each entry's restorable/non_restorable split is ALREADY computed by
        list_active_exclusions (it drives the "(N wells; M not revertable)" suffix) and held
        in self._active_excl_entries. We consult that here to SKIP rules that can revert
        nothing.
        """
        if self.master_df is None or self.master_df.empty:
            self.log("No data loaded to revert exclusions from.")
            return

        # Use the entries from the last refresh (already carry restorable/non_restorable and
        # match the listbox the user selected from) instead of recomputing the whole list.
        entries_now = self._active_excl_entries
        if not entries_now:
            self.log("No active exclusions to revert.")
            self.refresh_active_exclusions()
            return

        selected_labels, pending_target = self._collect_revert_inputs()
        if not selected_labels and not pending_target:
            self.log("[REVERT] Nothing selected and no pending rule matched a currently-excluded well.")
            return

        self.log(f"\n--- Reverting {len(selected_labels)} selected rule(s)"
                 f" + {len(pending_target)} criteria-targeted well(s) ---")

        config = self._build_recompute_config()
        before_excluded = self._excluded_well_set()   # to detect which wells get reverted

        # One shared resolution context for the whole batch revert. Precomputes the
        # invariant resolution work (normalised dates, stripped criteria columns, per-rule
        # criteria masks, restorability) ONCE so reverting N rules is not O(N x N) full-master
        # re-scans. Built on self.master_df (the object the restore calls mutate in place).
        resolve_ctx = build_resolve_ctx(self.master_df)

        total_restored = 0
        total_blocked = 0
        all_nonrestorable = []
        affected_files = set()
        any_decomposed = False

        # Collect for restoring summary
        def _tally(report):
            nonlocal total_restored, total_blocked, any_decomposed
            total_restored += report.get("restored_count", 0)
            total_blocked += report.get("blocked_by_other_rule_count", 0)
            all_nonrestorable.extend(report.get("nonrestorable", []))
            affected_files.update(report.get("affected_files", []))
            any_decomposed = any_decomposed or report.get("rule_decomposed", False)

        # (A) Whole-rule reverts
        # Labels captured before mutation; each restore_rule
        # only rewrites its OWN token, so the remaining selected labels stay resolvable
        entry_by_label = {e["label"]: e for e in entries_now}
        for label in selected_labels:
            e = entry_by_label.get(label)
            if e is not None and not e.get("restorable"):
                nrk = e.get("non_restorable", [])
                all_nonrestorable.extend(nrk)
                self.log(f"   [SKIP] '{self._display_label(label)}' has no revertable wells "
                         f"({len(nrk)} non-revertable) — not attempted.")
                continue
            _tally(restore_rule(self.master_df, label, config, recompute=False, ctx=resolve_ctx))

        # (B) Per-well criteria reverts
        if pending_target:
            for entry in list_active_exclusions(self.master_df, ctx=resolve_ctx):
                entry_wells = set(entry["wells"])
                inter = pending_target & entry_wells
                if not inter:
                    continue
                restorable_inter = inter & {(f, w) for (f, w) in entry.get("restorable", [])}
                if not restorable_inter:
                    all_nonrestorable.extend(
                        [(f, w, r) for (f, w, r) in entry.get("non_restorable", []) if (f, w) in inter])
                    continue
                if inter == entry_wells:
                    _tally(restore_rule(self.master_df, entry["label"], config,
                                        recompute=False, ctx=resolve_ctx))
                else:
                    _tally(restore_wells(self.master_df, entry["label"], inter, config,
                                         recompute=False, ctx=resolve_ctx))

        self.log(f"   [DONE] Restored {total_restored} well(s) across "
                 f"{len(affected_files)} file(s).")
        if any_decomposed:
            self.log("   [INFO] One or more rules were decomposed to per-well tokens "
                     "(partial restore).")

        # --- Single batched recompute over the UNION of affected files ---
        # The per-rule restore calls above ran with recompute=False: they flipped
        # Is_Excluded and rewrote the Applied_Exclusions blob in place, but deferred the
        # expensive derived-column recompute to here. Recomputing the union ONCE means each
        # affected file is recomputed a single time instead of once per rule (the previous
        # behaviour, which made a many-rule revert run the whole BRET pipeline N times).
        # Must run BEFORE built_master_index / refresh_* below, which read derived columns.
        if affected_files:
            self.master_df = recompute_master_after_exclusion(
                self.master_df, sorted(affected_files), config)

        # --- Refresh everything the revert touched ---
        # IMPORTANT: restore.py rewrote master_df['Applied_Exclusions'] in place but does not
        # touch rule_history_text. Re-sync it from the authoritative blob now so a later
        # apply_exclusions (which rebuilds the blob from rule_history_text) does not resurrect
        # the rule that was just reverted.
        self._sync_rule_history_from_master()
        self.clear_exclusion_list()
        self.built_master_index(source="master")
        self.refresh_plot_helper_options()
        self._update_quality_buttons_state()
        self.refresh_active_exclusions()

        # --- Non-restorable / blocked feedback: logged summary, not a hard block ---
        if total_blocked:
            self.log(f"   [INFO] {total_blocked} well(s) left excluded — still claimed by "
                     f"another active exclusion rule.")
        if all_nonrestorable:
            # De-duplicate and log each well that could not be reverted (with its reason).
            uniq = sorted(set(all_nonrestorable))
            self.log(f"   [WARNING] {len(uniq)} well(s) could not be reverted:")
            for (f, w, reason) in uniq:
                self.log(f"       - {f} / {w}: {reason}")

        # Re-run vehicle check
        # if reverted wells include vehicle wells + only for the recalculated conditions
        affected_keys = before_excluded ^ self._excluded_well_set()
        veh_detected = self._vehicle_warnings_for_affected(affected_keys)
        new_veh = [w for w in veh_detected if w['Display'] not in self.ignored_warnings]
        if new_veh:
            self.show_warning_review(new_veh)

        # Update the crc plot (if open) against the recomputed data
        self._refresh_crc_window()

        self.log("--- Exclusions reverted & data recomputed ---")

    # --- Concentration-response plot window (Exclude tab) ---
    def open_crc_window(self):
        """Open the concentration-response plot window for the current master."""
        if self.master_df is None or self.master_df.empty:
            self.log("[CRC] No data loaded to plot.")
            return
        if self.crc_window is not None:
            try:
                self.crc_window.win.deiconify()
                self.crc_window.win.lift()
                self.crc_window.refresh()
                return
            except tk.TclError:
                self.crc_window = None  # window was closed/destroyed; recreate below
        self.crc_window = crc_window.CRCWindow(
            self.main_gi,
            get_master_df=lambda: self.master_df,
            on_close=lambda: setattr(self, "crc_window", None),
            log_fn=self.log)

    def _refresh_crc_window(self):
        """Redraw the open CRC plot (if any) against the current master — called after every
        exclude/restore and whenever a new master is loaded/merged."""
        if self.crc_window is not None:
            try:
                self.crc_window.refresh()
            except tk.TclError:
                self.crc_window = None

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

    def _option_source_df(self):
        """
        Build the Tab-2 dropdown OPTION SOURCE, filtered by the current mode:
          exclude mode -> rows still in play   (Is_Excluded == False)
          revert  mode -> already-excluded rows (Is_Excluded == True)

        Derived directly from master_df (using processing.coerce_bool for the Is_Excluded
        test) so that REVERT mode can surface even fully-excluded columns, which
        built_master_index drops from self.master_index. Returns a frame carrying exactly
        the vocabulary columns the cascading logic expects: Ligand, Date ('%d.%m.%y'),
        Cell_Line, Condition (<- Transfection), Replicate, plus the scope columns
        Main_Plasmids and File_Name (used by the cascade only when their selectors are
        shown; otherwise they sit at "All" and act as no-ops).

        Bootstrap fallback: at the very first object-mode load, built_master_index ->
        refresh_filter_options runs BEFORE master_df is compiled. In that window master_df
        is empty, so we fall back to self.master_index (which at load time already holds the
        clean, exclusion-free vocabulary). Only valid for exclude mode (no exclusions yet).
        """
        cols = ["Ligand", "Date", "Cell_Line", "Condition", "Replicate",
                "Main_Plasmids", "File_Name"]
        want_excluded = (self.var_excl_mode is not None
                         and self.var_excl_mode.get() == "revert")

        df = self.master_df
        if df is not None and not df.empty:
            excl = (df['Is_Excluded'].map(coerce_bool)
                    if 'Is_Excluded' in df.columns else pd.Series(False, index=df.index))
            work = df[excl] if want_excluded else df[~excl]
            # Drop empty/unknown columns (consistent with _build_index_from_master).
            if not work.empty:
                work = work[~work['Transfection'].astype(str).str.contains("Empty", na=False)]
                work = work[~work['Cell_Line'].astype(str).str.startswith("Unknown", na=False)]
            if work.empty:
                return pd.DataFrame(columns=cols)
            return pd.DataFrame({
                "Ligand": work['Ligand'].astype(str),
                "Date": parse_date_series(work['Date'], context="_option_source_df").dt.strftime('%d.%m.%y'),
                "Cell_Line": work['Cell_Line'].astype(str),
                "Condition": work['Transfection'].astype(str),
                "Replicate": work['Replicate'].astype(str),
                "Main_Plasmids": work['Main_Plasmids'].astype(str),
                "File_Name": work['File_Name'].astype(str),
            })

        # --- Bootstrap fallback (object-mode load, master_df not yet compiled) ---
        if want_excluded or self.master_index is None or self.master_index.empty:
            return pd.DataFrame(columns=cols)
        return self.master_index[cols].copy()

    def refresh_filter_options(self):
        """Called during built master index. Updates dropdown options of date, cell line and condition."""
        idx = self._option_source_df()
        if idx.empty: return

        # Reset Variables
        self.var_lig.set("All")
        self.var_date.set("All")
        self.var_cell.set("All")
        self.var_cond.set("All")
        self.var_repl.set("All")
        self.var_row.set("All")
        self.var_main.set("All")
        self.var_file.set("All")
        self.cb_row.config(state="disabled")

        # Show + populate the Main_Plasmids + File_Name scope selectors iff meaningful
        self._refresh_scope_selectors()

        # If only one ligand exists, default to it and disable the box.
        unique_ligands = sorted(idx['Ligand'].dropna().unique().tolist())
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

    def _refresh_scope_selectors(self):
        """
        Show/populate the Main_Plasmids and File_Name scope selectors only when meaningful
        (case only possible for merged masters), and otherwise hide them and reset their vars to "All".

        Visibility rules (each independent):
          * Main_Plasmids — shown iff master_df['Main_Plasmids'].nunique() > 1, i.e. the
            master spans more than one experiment (a merged blob).
          * File_Name — shown iff some File_Name carries the merge auto-rename marker '#'
            (a FILENAME_DATA_COLLISION resolved during a blob merge).

        Option lists are seeded from the master index (Main_Plasmids, File_Name) with an
        "All" sentinel first; update_dropdown_options then refines them as part of the
        cascade.
        """
        df = self.master_df
        idx = self.master_index

        # --- Main_Plasmids ---
        show_main = (df is not None and not df.empty and 'Main_Plasmids' in df.columns
                     and df['Main_Plasmids'].nunique(dropna=True) > 1)
        main_vals = (sorted(idx['Main_Plasmids'].dropna().astype(str).unique().tolist())
                     if idx is not None and not idx.empty and 'Main_Plasmids' in idx.columns
                     else [])
        self.cb_main['values'] = ["All"] + main_vals
        if show_main:
            self.lbl_main.grid(); self.cb_main.grid()
        else:
            self.var_main.set("All")
            self.lbl_main.grid_remove(); self.cb_main.grid_remove()

        # --- File_Name ---
        show_file = (df is not None and not df.empty and 'File_Name' in df.columns
                     and df['File_Name'].astype(str).str.contains('#', regex=False, na=False).any())
        file_vals = (sorted(idx['File_Name'].dropna().astype(str).unique().tolist())
                     if idx is not None and not idx.empty and 'File_Name' in idx.columns
                     else [])
        self.cb_file['values'] = ["All"] + file_vals
        if show_file:
            self.lbl_file.grid(); self.cb_file.grid()
        else:
            self.var_file.set("All")
            self.lbl_file.grid_remove(); self.cb_file.grid_remove()

    def update_dropdown_options(self, trigger_source=None):
        """
        Dynamically updates the values of all dropdowns based on the current selection of others.
        trigger_source: The name of the field that triggered the update.

        The OPTION SOURCE is mode-filtered via _option_source_df (exclude -> not-excluded
        rows, revert -> excluded rows); everything else (cascading narrowing, Row
        enable/disable) is unchanged.
        """
        idx = self._option_source_df()
        if idx.empty: return

        # Map Columns to their UI Components
        field_map = {
            "Ligand": (self.var_lig, self.cb_lig),
            "Date": (self.var_date, self.cb_date),
            "Cell_Line": (self.var_cell, self.cb_cell),
            "Condition": (self.var_cond, self.cb_cond),
            "Replicate": (self.var_repl, self.cb_rep),
            # Scope selectors, only shown if relevant (otherwise hidden with "All")
            "Main_Plasmids": (self.var_main, self.cb_main),
            "File_Name": (self.var_file, self.cb_file),
        }

        # Get current selection
        current_selections = {col: var.get() for col, (var, _) in field_map.items()}

        # Iterate through each field
        for param, (target_var, target_widget) in field_map.items():
            # Built an all true mask
            mask = pd.Series(True, index=idx.index)

            # Apply filters from ALL OTHER fields
            for curr_param, val in current_selections.items():
                # Skip the changed dropdown (curr_param) so it doesn't filter itself
                if curr_param != param and val != "All" and val != "":
                    mask &= (idx[curr_param].astype(str) == str(val))

            # Extract unique values using the mask directly
            valid_options = sorted(idx.loc[mask, param].dropna().unique().tolist())

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
            "Row": self.var_row.get(),
            "Main_Plasmids": self.var_main.get(),
            "File_Name": self.var_file.get(),
        }

        # Check for duplicates or empty
        rule_str = f"Ligand: {rule['Ligand']} | Date: {rule['Date']} | "\
                   f"Cell: {rule['Cell_Line']} | Cond: {rule['Condition']} | "\
                   f"Rep:{rule['Replicate']} | Row:{rule['Row']}"
        # Scope segments are appended ONLY when meaningful (not "All"/empty)
        if rule['Main_Plasmids'] not in ("All", ""):
            rule_str += f" | Main: {rule['Main_Plasmids']}"
        if rule['File_Name'] not in ("All", ""):
            rule_str += f" | File: {rule['File_Name']}"

        self.pending_exclusions.append(rule)
        self.lb_exclusions.insert(tk.END, rule_str)

    def clear_exclusion_list(self):
        self.pending_exclusions = []
        self.lb_exclusions.delete(0, tk.END)

    # Columns that together identify which master rows an exclusion target refers
    # to. File_Name + Well_ID alone could collide if two loaded files happen to
    # share a name (e.g. same filename in two folders), so the match is widened to
    # the well's biological identity. Both _resolve_rule_to_targets (which builds
    # the target tuples) and apply_exclusions (which builds the match mask) read
    # this list, so the two keys are guaranteed to line up by construction.
    _EXCLUSION_KEY_COLS = ('File_Name', 'Well_ID', 'Ligand', 'Transfection',
                           'Cell_Line', 'Main_Plasmids')

    def _exclusion_key_cols(self, df):
        """The exclusion-key columns actually present in df, in fixed order.
        File_Name + Well_ID are always present; the identity columns are added
        when available (they are part of MASTER_COLUMNS, so normally all six)."""
        return [c for c in self._EXCLUSION_KEY_COLS if c in df.columns]

    def _rule_mask(self, rule, norm_dates):
        """
        Boolean mask over master_df for ONE pending criteria rule — the single source of
        the pending-rule filter semantics, shared by both the exclude path
        (_resolve_rule_to_targets) and the revert path (_resolve_rule_to_wells) so the two
        can never drift. "All ligands" also covers the single-locked-ligand case
        (cb_lig disabled). norm_dates is master_df['Date'] pre-normalised to '%d.%m.%y'.
        """
        df = self.master_df
        mask = pd.Series(True, index=df.index)

        ligand_is_all = (rule.get('Ligand', 'All') in ("All", "")
                         or str(self.cb_lig['state']) == 'disabled')
        if not ligand_is_all:
            mask &= (df['Ligand'].astype(str) == str(rule['Ligand']))
        if rule.get('Date', 'All') != "All":
            mask &= (norm_dates == rule['Date'])
        if rule.get('Cell_Line', 'All') != "All":
            mask &= (df['Cell_Line'].astype(str) == str(rule['Cell_Line']))
        if rule.get('Condition', 'All') != "All":
            # Index vocabulary "Condition" maps to master_df "Transfection".
            mask &= (df['Transfection'].astype(str) == str(rule['Condition']))
        rep = rule.get('Replicate', 'All')
        if rep not in ("All", ""):
            mask &= (df['Replicate'].astype(str) == str(rep))
        row = rule.get('Row', 'All')
        if row not in ("All", ""):
            mask &= (df['Plate_Row'].astype(str) == str(row))
        # Scope narrowing (post-merge masters). Skipped when "All"/empty
        main = rule.get('Main_Plasmids', 'All')
        if main not in ("All", "") and 'Main_Plasmids' in df.columns:
            mask &= (df['Main_Plasmids'].astype(str) == str(main))
        fname = rule.get('File_Name', 'All')
        if fname not in ("All", "") and 'File_Name' in df.columns:
            mask &= (df['File_Name'].astype(str) == str(fname))
        return mask

    def _resolve_rule_to_targets(self, rule, norm_dates):
        """
        Resolves one pending exclusion rule to a set of identity tuples directly on
        master_df (the single source of truth). Same rule semantics as before,
        including the "whole date" branch — which now resolves to every matching
        well on that date rather than flagging PrResult.is_excluded.

        Each returned tuple is keyed on self._EXCLUSION_KEY_COLS — i.e. not just
        (File_Name, Well_ID) but also Ligand / Transfection / Cell_Line /
        Main_Plasmids — so that a duplicate file name cannot cause the wrong rows
        to be flagged. apply_exclusions builds its match mask from the same column
        list, keeping the two sides consistent.

        norm_dates is master_df['Date'] pre-normalised to the '%d.%m.%y' dropdown form.
        """
        df = self.master_df
        mask = self._rule_mask(rule, norm_dates)

        sub = df.loc[mask]
        if sub.empty:
            self.log(f"   [WARNING] Rule {rule} matched 0 records.")
            return set()
        key_cols = self._exclusion_key_cols(df)
        return set(zip(*[sub[c].astype(str) for c in key_cols]))

    def _resolve_rule_to_wells(self, rule, norm_dates):
        """
        Resolve one pending criteria rule to a set of (File_Name, Well_ID) tuples, using
        the SAME mask as _resolve_rule_to_targets (via _rule_mask). Used by the revert
        path, which needs plain (file, well) pairs for the restore.py API rather than the
        full exclusion-identity tuples the add-only flag mask uses.
        """
        df = self.master_df
        mask = self._rule_mask(rule, norm_dates)
        sub = df.loc[mask]
        if sub.empty:
            return set()
        return set(zip(sub['File_Name'].astype(str), sub['Well_ID'].astype(str)))

    def _build_recompute_config(self):
        """
        Returns a ProcessingConfig for the master-native recompute engine. In fresh
        mode this is the load-time config (carries the persisted baseline index). In
        CSV/import mode it is derived from master_df (the engine derives per-file
        baseline/labeling itself; config only supplies fallbacks + thresholds).
        """
        if self.current_config is not None and self.experiment:
            return self.current_config

        df = self.master_df
        is_labeling = False
        if df is not None and not df.empty and 'Replicate' in df.columns:
            is_labeling = (df['Replicate'].astype(str) == "labeling control").any()

        baseline_idx = None
        if df is not None and not df.empty and {'File_Name', 'Well_ID', 'Time_(min)'}.issubset(df.columns):
            sample_file = df['File_Name'].iloc[0]
            sample_well = df['Well_ID'].iloc[0]
            m = (df['File_Name'] == sample_file) & (df['Well_ID'] == sample_well)
            times = pd.to_numeric(df.loc[m, 'Time_(min)'], errors='coerce').sort_values().tolist()
            if 0.0 in times:
                baseline_idx = times.index(0.0)

        lum_thr, veh_thr = self._get_thresholds()
        return ProcessingConfig(
            labeling_correction=is_labeling,
            plate_layout=build_plate_layout(is_labeling),
            baseline_end_index=baseline_idx,
            lum_threshold=lum_thr,
            vehicle_warning_threshold=veh_thr,
        )

    def _sync_rule_history_from_master(self):
        """
        Re-derive self.rule_history_text from the AUTHORITATIVE Applied_Exclusions blob on
        master_df.

        rule_history_text is the app's running provenance string. apply_exclusions rebuilds
        the whole blob from it, and the Excel export writes it out. restore.py, however,
        rewrites master_df['Applied_Exclusions'] in place on a revert WITHOUT touching
        rule_history_text — Calling this after any revert (and after a data-
        source swap that brings its own blob) keeps rule_history_text equal to the blob.
        """
        blob = ""
        if (self.master_df is not None and not self.master_df.empty
                and 'Applied_Exclusions' in self.master_df.columns):
            s = self.master_df['Applied_Exclusions'].dropna()
            if not s.empty:
                blob = str(s.iloc[0]).strip()
        if not blob or blob.lower() == "none":
            self.rule_history_text = ""
        else:
            self.rule_history_text = blob.replace(" || ", "\n")

    def apply_exclusions(self):
        """
        Unified exclusion entry point for BOTH fresh-load and CSV/import modes.

        master_df is the single source of truth for exclusion state. Pending rules are
        resolved to (File_Name, Well_ID) targets, those rows are flagged
        Is_Excluded=True (add-only), and the affected files are recomputed via
        recompute_master_after_exclusion — the one and only re-application path. The
        object pipeline (process_bret_measurement) is NOT re-run here (it runs only at
        initial load). Only the vehicle check is auto-re-run on the recomputed data (gated to
        when vehicle / labeling-control wells were affected, scoped to the recalculated
        conditions) and filtered against previously-ignored warnings. The lum check is NOT
        re-run on exclusion: use the manual "Re-run Lum Check" button for an on-demand full scan).
        """
        if not self.pending_exclusions:
            self.log("No exclusion rules defined.")
            return
        if self.master_df is None or self.master_df.empty:
            self.log("No data loaded to apply exclusions to.")
            return

        self.log(f"\n--- Applying {len(self.pending_exclusions)} Exclusion Rules ---")

        # --- Save applied rules as text for display ---
        pending_lines = self.lb_exclusions.get(0, tk.END)
        new_text_block = "\n".join(pending_lines)
        if self.rule_history_text:
            self.rule_history_text += "\n" + new_text_block
        else:
            self.rule_history_text = new_text_block

        # --- Resolve every pending rule to (File_Name, Well_ID) targets on master_df ---
        norm_dates = parse_date_series(self.master_df['Date'], context="apply_exclusions").dt.strftime('%d.%m.%y')
        targets = set()
        for rule in self.pending_exclusions:
            targets |= self._resolve_rule_to_targets(rule, norm_dates)

        if not targets:
            self.log("   [WARNING] Pending rules matched 0 records.")
            self.clear_exclusion_list()
            return

        # --- Flag Is_Excluded=True on the resolved rows (add-only / monotonic) ---
        # The match key is the full exclusion-identity tuple (File_Name + Well_ID +
        # Ligand/Transfection/Cell_Line/Main_Plasmids), built from the exact same
        # column list _resolve_rule_to_targets used, so a shared file name cannot
        # cause the wrong rows to be flagged.
        key_cols = self._exclusion_key_cols(self.master_df)
        key_index = pd.MultiIndex.from_arrays(
            [self.master_df[c].astype(str) for c in key_cols])
        target_mask = pd.Series(key_index.isin(list(targets)), index=self.master_df.index)

        before_excluded = self._excluded_well_set()   # to detect which wells actually flip
        self.master_df.loc[target_mask, 'Is_Excluded'] = True
        affected_files = sorted(self.master_df.loc[target_mask, 'File_Name'].astype(str).unique().tolist())
        self.log(f"   [DONE] Flagged {len(targets)} well(s) across "
                 f"{len(affected_files)} file(s) as excluded.")

        # --- Append the full applied-exclusion provenance to all rows ---
        # v2.0.5: rules are joined with " || "
        self.master_df['Applied_Exclusions'] = self.rule_history_text.replace("\n", " || ")

        # --- Recompute affected files on the flat master (single engine path) ---
        config = self._build_recompute_config()
        self.master_df = recompute_master_after_exclusion(self.master_df, affected_files, config)

        # --- Rebuild the exclusion index + summary from the updated master ---
        self.built_master_index(source="master")

        # Clear pending list now that the rules are applied
        self.clear_exclusion_list()

        # --- AUTO-RERUN vehicle check on the recomputed data ---
        # only if excluded wells include vehicle wells + only for the recalculated conditions
        affected_keys = before_excluded ^ self._excluded_well_set()
        detected = self._vehicle_warnings_for_affected(affected_keys)

        new_warnings = [w for w in detected if w['Display'] not in self.ignored_warnings]
        if new_warnings:
            self.show_warning_review(new_warnings)
        elif detected:
            self.log(f"[INFO] {len(detected)} warnings detected but previously ignored.")

        # --- Refresh the plot helper + re-evaluate quality-button availability ---
        self.refresh_plot_helper_options()
        self._update_quality_buttons_state()

        # --- Refresh the shared Active-Exclusions view (both tabs) + Tab-2 comboboxes ---
        # repopulates both active boxes and re-evaluates the Revert-mode enable state now that the excluded set changed.
        self.refresh_active_exclusions() #

        # Update the crc plot (if open) against the recomputed data
        self._refresh_crc_window()

        self.log("--- Exclusions applied & data recomputed ---")

    # --- Quality-check helpers (shared by auto-rerun and manual buttons) --- #

    def _get_thresholds(self):
        """Reads the lum + vehicle thresholds from the Tab-1 entries (single source)."""
        try:
            lum = self.var_lum_threshold.get()
        except tk.TclError:
            lum = 100
        try:
            veh = self.var_vehicle_threshold.get()
        except tk.TclError:
            veh = 0.2
        return lum, veh

    def _collect_quality_warnings(self, files, run_lum, run_vehicle,
                                  lum_threshold, vehicle_threshold, log_lum_skips=False):
        """
        Re-runs the requested checks across `files` by reconstructing each file's
        inputs from master_df (shared reconstruction helper) and calling
        check_luminescence / check_vehicle_wells VERBATIM. Returns the full detected
        warning list (no ignored-filtering — callers decide whether to filter).

        Lum check: skips files lacking donor data (graceful, with a log line). Donor
        counts are reconstructed raw and the file's currently-excluded wells are NaN'd
        in a copy before the check (mirrors process_bret_measurement step 3).
        Vehicle check: always runs (needs only Veh_Norm); it self-skips excluded wells.
        """
        warnings = []
        missing_donor = []

        for fname in files:
            inputs = reconstruct_file_inputs(self.master_df, fname)
            if inputs is None:
                continue

            if run_vehicle and inputs['veh_norm_wide'] is not None:
                warnings.extend(check_vehicle_wells(
                    inputs['veh_norm_wide'], inputs['plate_blocks'],
                    inputs['column_metadata'], inputs['excluded_wells'],
                    vehicle_threshold, inputs['date_str']))

            if run_lum:
                if not inputs['has_donor']:
                    missing_donor.append(fname)
                    continue
                donor = inputs['donor_wide'].copy()
                for w in inputs['excluded_wells']:
                    if w in donor.columns:
                        donor[w] = float('nan')
                warnings.extend(check_luminescence(
                    donor, inputs['column_metadata'], lum_threshold, inputs['date_str']))

        if run_lum and log_lum_skips and missing_donor:
            self.log(f"   [LUM] Skipped {len(missing_donor)} file(s) lacking donor data: "
                     f"{', '.join(missing_donor)}")
        return warnings

    def _excluded_well_set(self):
        """Current set of (File_Name, Well_ID) flagged Is_Excluded in master_df."""
        df = self.master_df
        if df is None or df.empty or 'Is_Excluded' not in df.columns:
            return set()
        mask = df['Is_Excluded'].map(coerce_bool)
        return set(zip(df.loc[mask, 'File_Name'].astype(str),
                       df.loc[mask, 'Well_ID'].astype(str)))

    def _vehicle_warnings_for_affected(self, affected_keys):
        """Vehicle check restricted to what an exclude/revert actually changed.

        `affected_keys` is the set of (File_Name, Well_ID) whose Is_Excluded state flipped.
        A block's vehicle-check result can move only when one of ITS OWN wells that feeds the
        vehicle normalisation changed:
          * a vehicle well (Is_Vehicle == True) — it sets the block's vehicle mean; OR
          * a labeling-control well (Replicate == "labeling control") — under labeling
            correction it sets the block's background subtraction, which shifts every well's
            baseline-corrected value (incl. the vehicle wells), hence the vehicle mean.
        So the check runs ONLY if a changed well is one of those, and only its warnings for
        the RECALCULATED conditions (the Ligand / Transfection / Cell_Line of those changed
        wells) are kept. Returns warning dicts (NOT filtered against ignored_warnings;
        callers decide).
        """
        df = self.master_df
        if df is None or df.empty or not affected_keys:
            return []
        key = pd.MultiIndex.from_arrays(
            [df['File_Name'].astype(str), df['Well_ID'].astype(str)])
        changed = df[pd.Series(key.isin(list(affected_keys)), index=df.index)]
        if changed.empty:
            return []
        relevant_mask = pd.Series(False, index=changed.index)
        if 'Is_Vehicle' in changed.columns:
            relevant_mask |= changed['Is_Vehicle'].map(coerce_bool)
        if 'Replicate' in changed.columns:
            relevant_mask |= changed['Replicate'].astype(str) == "labeling control"
        relevant = changed[relevant_mask]
        if relevant.empty:
            # no vehicle / labeling-control wells changed -> skip the check
            self.log("   [VEHICLE] Vehicle check skipped — no vehicle or labeling-control "
                     "wells affected.")
            return []
        affected_conditions = set(zip(relevant['Ligand'].astype(str),
                                      relevant['Transfection'].astype(str),
                                      relevant['Cell_Line'].astype(str)))
        affected_files = sorted(relevant['File_Name'].astype(str).unique())
        _, veh_thr = self._get_thresholds()
        self.log(f"   [VEHICLE] Re-running vehicle check on {len(affected_files)} file(s), "
                 f"{len(affected_conditions)} condition(s) (threshold={veh_thr}).")
        warnings = self._collect_quality_warnings(
            affected_files, run_lum=False, run_vehicle=True,
            lum_threshold=0, vehicle_threshold=veh_thr)
        return [w for w in warnings
                if (str(w['Ligand']), str(w['Condition']), str(w['Cell_Line']))
                in affected_conditions]

    def rerun_lum_check(self):
        """
        On-demand luminescence re-run (Tab 1). Reads the lum threshold from Tab 1,
        checks every file with donor data (skipping & logging those without), and
        shows ALL detected warnings — deliberately bypassing self.ignored_warnings
        so previously-dismissed warnings resurface. Excluded wells are still
        suppressed (NaN'd before the check). Reuses show_warning_review unchanged.
        """
        if self.master_df is None or self.master_df.empty:
            self.log("[LUM] No data loaded.")
            return
        lum_thr, _ = self._get_thresholds()
        self.log(f"\n--- Re-running Luminescence Check (threshold={lum_thr}) ---")
        files = self.master_df['File_Name'].astype(str).unique().tolist()
        detected = self._collect_quality_warnings(
            files, run_lum=True, run_vehicle=False,
            lum_threshold=lum_thr, vehicle_threshold=0.0, log_lum_skips=True)
        if detected:
            self.show_warning_review(detected)
        else:
            self.log("[LUM] No low-luminescence warnings detected.")

    def rerun_vehicle_check(self):
        """
        On-demand vehicle re-run (Tab 1). Reads the vehicle threshold from Tab 1 and
        shows ALL detected vehicle warnings — bypassing self.ignored_warnings (same
        deliberate divergence from the auto-rerun). Reuses show_warning_review.
        """
        if self.master_df is None or self.master_df.empty:
            self.log("[VEHICLE] No data loaded.")
            return
        _, veh_thr = self._get_thresholds()
        self.log(f"\n--- Re-running Vehicle Check (threshold={veh_thr}) ---")
        files = self.master_df['File_Name'].astype(str).unique().tolist()
        detected = self._collect_quality_warnings(
            files, run_lum=False, run_vehicle=True,
            lum_threshold=0, vehicle_threshold=veh_thr)
        if detected:
            self.show_warning_review(detected)
        else:
            self.log("[VEHICLE] No vehicle warnings detected.")

    def _set_tooltip(self, widget, text):
        """Attach (or update) a hover tooltip on a widget."""
        if widget is None:
            return
        tip = getattr(widget, "_tooltip_obj", None)
        if tip is None:
            tip = _Tooltip(widget)
            widget._tooltip_obj = tip
        tip.set_text(text)

    def _update_quality_buttons_state(self):
        """
        Re-evaluates the enabled/disabled state of the two re-run buttons. Called
        whenever the data source changes (load, import, enrich, exclusion). The
        vehicle button is enabled whenever data is present (Veh_Norm always exists);
        the lum button is enabled only if the master carries any donor data.
        """
        has_data = self.master_df is not None and not self.master_df.empty

        if self.btn_rerun_vehicle is not None:
            self.btn_rerun_vehicle.config(state="normal" if has_data else "disabled")
            self._set_tooltip(self.btn_rerun_vehicle,
                              "Re-run the vehicle-deviation check using the Tab-1 threshold."
                              if has_data else "Load or import data to enable.")

        has_donor = (has_data and 'Donor_Raw_kinetic' in self.master_df.columns
                     and self.master_df['Donor_Raw_kinetic'].notna().any())
        if self.btn_rerun_lum is not None:
            if has_donor:
                self.btn_rerun_lum.config(state="normal")
                self._set_tooltip(self.btn_rerun_lum,
                                  "Re-run the low-luminescence check on donor counts.")
            else:
                self.btn_rerun_lum.config(state="disabled")
                if has_data:
                    msg = ("Donor (luminescence) data is missing from this dataset "
                           "— the lum check cannot run.")
                    logger.info("Lum re-run disabled: no Donor_Raw_kinetic data in master.")
                else:
                    msg = "Load or import data to enable."
                self._set_tooltip(self.btn_rerun_lum, msg)

    def setup_plot_helper_tab(self):
        """Builds the GUI for tab 3 plot helper"""

        # --- Import Master CSV ---
        src_frame = tk.Frame(self.tab_plot_helper)
        src_frame.pack(fill="x", padx=10, pady=10)

        # Import button on the right; data-source title + value left-aligned above it
        tk.Button(src_frame, text="Import Master CSV",
                  command=self.import_master_csv).pack(side="right")

        src_label_frame = tk.Frame(src_frame)
        src_label_frame.pack(side="left", anchor="w")
        tk.Label(src_label_frame, text="Data Source",
                 font=("Arial", 9, "bold")).pack(anchor="w")
        self.lbl_data_source = tk.Label(src_label_frame, text="No Data Loaded",
                                        justify="left")
        self.lbl_data_source.pack(anchor="w")

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
        self.var_category = tk.StringVar(value="")
        tk.Label(type_frame, text="Category:").grid(row=0, column=0, padx=5, pady=2, sticky="w")
        self.combo_category = ttk.Combobox(type_frame, textvariable=self.var_category, state="readonly", width=10)
        self.combo_category.grid(row=0, column=1, padx=5, pady=5, sticky="w")
        # Bind change event to update the second dropdown
        self.combo_category.bind("<<ComboboxSelected>>", self.update_subtype_options)

        # Specific Data Type
        self.var_specific_type = tk.StringVar(value="")
        tk.Label(type_frame, text="Subtype:").grid(row=1, column=0, padx=5, pady=2, sticky="w")
        self.combo_specific = ttk.Combobox(type_frame, textvariable=self.var_specific_type, state="readonly", width=60)
        self.combo_specific.grid(row=1, column=1, padx=5, pady=5, sticky="w")

        # --- Layout Selection ---
        # Groupy by
        tk.Label(type_frame, text="Group By:").grid(row=3, column=0, padx=5, pady=5, sticky="w")
        self.var_group_by = tk.StringVar(value="None")
        self.combo_group = ttk.Combobox(type_frame, textvariable=self.var_group_by, state="readonly", width=15)
        self.combo_group['values'] = ["None", "Cell Line", "Transfection"]
        self.combo_group.grid(row=3, column=1, padx=5, pady=5, sticky="w")

        # Conc layout (Only applies if a Kinetic or Bargraph type is chosen)
        tk.Label(type_frame, text="Conc. Layout:").grid(row=0, column=2, padx=5, pady=5, sticky="e")
        self.lb_conc_layout = tk.Listbox(type_frame, selectmode="multiple", height=6, exportselection=False)
        self.lb_conc_layout.grid(row=0, column=3, rowspan = 2, padx=3, pady=5, sticky="nsew")
        tk.Button(type_frame, text="Select All Rows", command=lambda: self.lb_conc_layout.select_set(0, tk.END)).grid(
            row=2, column=3, pady=5, sticky="nsew")
        tk.Button(type_frame, text="Deselect All Rows", command=lambda: self.lb_conc_layout.selection_clear(0, tk.END)).grid(
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
        Toggles Conc. Layout listbox mode and Group By options based on category.
        """
        if self.master_df is None or self.master_df.empty: return

        selected_cat = self.var_category.get()  # e.g. "CRC"
        subtype_map = DATA_TYPE_MAP.get(selected_cat, {})

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

        # --- Toggle Conc. Layout Listbox ---
        if selected_cat in ("kinetic", "bargraph", "heatmap"):
            self.lb_conc_layout.config(state="normal")
            # Single select for bargraph and heatmap
            if selected_cat in SINGLE_CONC_CATEGORIES:
                self.lb_conc_layout.config(selectmode="browse")  # Single selection
                # Auto-select first if nothing selected
                if not self.lb_conc_layout.curselection() and self.lb_conc_layout.size() > 0:
                    self.lb_conc_layout.select_set(0)
            else:
                self.lb_conc_layout.config(selectmode="multiple")
        else:
            self.lb_conc_layout.config(state="disabled")

        # --- Toggle Group By options ---
        if selected_cat in REQUIRES_GROUP_BY:
            self.combo_group['values'] = ["Cell Line", "Transfection"]
            if self.var_group_by.get() == "None":
                self.var_group_by.set("Cell Line")
        else:
            self.combo_group['values'] = ["None", "Cell Line", "Transfection"]

    def refresh_plot_helper_options(self):
        """Populates the list boxes in the plot helper from master df (opt. imported csv file)."""
        if self.master_df is None or self.master_df.empty:
            self.btn_run_plot_helper.config(state="disabled")
            self._set_data_source("No Data")
            return

        if not self.experiment:
            # CSV mode -> status already set
            pass
        else:
            #  Fresh Analysis mode
            experiment = self.master_df['Main_Plasmids'].unique()[0] if 'Main_Plasmids' in self.master_df.columns else "Experiment"
            self._set_data_source(f"Internal: {experiment}")

        self.btn_run_plot_helper.config(state="normal")

        available_categories = set()

        # Iterate over high-level keys ("kinetic", "CRC")
        for cat, sub_map in DATA_TYPE_MAP.items():
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

        # Temporarily enable conc layout to populate list box
        self.lb_conc_layout.config(state="normal")

        # Clear
        self.lb_ligands.delete(0, tk.END)
        self.lb_exp_cells.delete(0, tk.END)
        self.lb_exp_trans.delete(0, tk.END)
        self.lb_conc_layout.delete(0, tk.END)

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

        # Populate Conc. Layout
        self.conc_row_lookup = build_row_info(df)
        for r in self.conc_row_lookup: self.lb_conc_layout.insert(tk.END, r)
        self.lb_conc_layout.select_set(0, tk.END)  # Default to all

        # Disable Conc. layout if CRC is currently selected
        if self.var_category.get() != "kinetic":
            self.lb_conc_layout.config(state="disabled")

    def run_plot_helper(self):
        """Collects GUI selections and calls the core export engine."""
        if self.master_df is None or self.master_df.empty: return

        # Get Selections
        cells = [self.lb_exp_cells.get(i) for i in self.lb_exp_cells.curselection()]
        transfections = [self.lb_exp_trans.get(i) for i in self.lb_exp_trans.curselection()]
        ligands = [self.lb_ligands.get(i) for i in self.lb_ligands.curselection()]
        rows = [self.lb_conc_layout.get(i) for i in self.lb_conc_layout.curselection()]
        # Resolve display strings to structured filter criteria
        conc_filters = [self.conc_row_lookup[r] for r in rows if r in self.conc_row_lookup]

        category = self.var_category.get()
        specific_type = self.var_specific_type.get()
        group_by = self.var_group_by.get()

        if not category or not specific_type: return

        # Validate: heatmap requires group_by
        if category in REQUIRES_GROUP_BY and group_by == "None":
            self.log("[ERROR] Heatmap export requires a 'Group By' selection.")
            return

        # Validate: bargraph/heatmap require exactly one conc selection
        if category in SINGLE_CONC_CATEGORIES and len(conc_filters) != 1:
            self.log(f"[ERROR] {category.capitalize()} export requires exactly one concentration selection.")
            return

        internal_name = DATA_TYPE_MAP.get(category, {}).get(specific_type)

        if not internal_name:
            logger.error(f"Unknown data type selected: {category} - {specific_type}")
            return

        config = {
            'category': category,
            'cells': cells,
            'transfections': transfections,
            'ligands': ligands,
            'data_types': [internal_name], # As list for engine compatibility with default export
            'group_by': group_by,
            'conc_mode': conc_filters
        }

        file_path = filedialog.asksaveasfilename(
            defaultextension=".xlsx",
            filetypes=[("Excel", "*.xlsx")],
            title="Save Custom Export"
        )
        if file_path:
            # Call the core engine with the unified DF
            self.write_excel_export(file_path, self.master_df, config)

    #  Shared loaders (used by Import & by the Merge tab)
    def _load_master_df(self, path):
        """
        Read a master CSV from disk, validate its core columns, migrate it to the current
        schema, and coerce Is_Excluded to real booleans. Shared by import_master_csv and
        the Merge tab's Add-Master-CSV path.

        Returns (df, was_modified, was_fixed).
        Raises ValueError if the file lacks required master columns.
        """
        df = pd.read_csv(path, low_memory=False)

        required = ["Transfection", "Cell_Line", "Ligand", "Kinetic_Mean"]
        if not all(col in df.columns for col in required):
            raise ValueError("Invalid CSV format. Required master columns are missing.")

        # Backward compatibility: fill missing columns and clean legacy data.
        df, was_modified, was_fixed = ensure_master_csv_schema(df, log_fn=self.log)

        # Normalize Is_Excluded to real booleans (CSV may carry "True"/"False" strings,
        # 1/0, etc.) so the unified engine / merge backend can trust it directly.
        if 'Is_Excluded' in df.columns:
            df['Is_Excluded'] = df['Is_Excluded'].map(coerce_bool)

        return df, was_modified, was_fixed

    def _process_folder_to_master(self, folder_paths, config):
        """
        Process one or more experiment subfolders into a master-shaped DataFrame WITHOUT
        touching application state. Used by the Merge tab to add an experiment folder as a
        merge source.

        Mirrors the object pipeline (scan -> process_bret_measurement -> compile) but
        writes to a LOCAL frame and returns it: it does NOT assign self.master_df, does NOT
        flip the warning/export buttons.
        It DOES run the interactive main-plasmids selection — but purely on the local experiment
        list (it never reads or mutates self.experiment). Ligand dialogs are reused via the
        passed ProcessingConfig callbacks, exactly as collect_files builds them. A folder
        source is enriched by construction (Donor/Acceptor raw channels are always populated
        by compile), so it always passes the merge channel gate.

        Returns the compiled master-shaped DataFrame, an empty DataFrame if there is nothing
        to process, or ``None`` if the user cancels the main-plasmids selection.
        """
        if not folder_paths:
            return pd.DataFrame()

        experiment = scan_and_load_folders(folder_paths, log_fn=self.log)
        if not experiment:
            return pd.DataFrame()

        # Prompt user to pick main plasmids, folder spans more than one main_plasmids set
        selected_experiment, _selected_str = self._select_main_plasmids(experiment)
        if selected_experiment is None:
            self.log("[MERGE] Folder add cancelled at main-plasmids selection.")
            return None
        experiment = selected_experiment

        for folder in experiment:
            if not folder.protocol:
                continue
            for result in folder.results:
                process_bret_measurement(result, folder.protocol, config)

        # Provenance directory for the Path column (common root of the scanned subfolders).
        try:
            directory = (folder_paths[0] if len(folder_paths) == 1
                         else os.path.commonpath(folder_paths))
        except ValueError:
            directory = folder_paths[0]

        return self.compile_master_dataframe(experiment=experiment, directory=directory,
                                             rule_history_text="")

    def setup_merge_tab(self):
        """Builds the Merge tab: a Sources list, a Summary preview (mirrors Tab 1, with an
        extra Source column), and the Merge action. Adding sources never mutates
        self.master_df — only Merge does. Compatibility is checked automatically at merge
        time and reported to the app log."""
        # --- Sources ---
        src_frame = tk.LabelFrame(self.tab_merge, text="Sources")
        src_frame.pack(fill="both", expand=True, padx=10, pady=(10, 5))

        tree_box = tk.Frame(src_frame)
        tree_box.pack(fill="both", expand=True, padx=5, pady=5)
        tree_scroll = ttk.Scrollbar(tree_box, orient="vertical")
        tree_scroll.pack(side="right", fill="y")
        self.merge_tree = ttk.Treeview(
            tree_box, columns=("Label", "Kind", "Files", "Rows"),
            show="headings", height=6, yscrollcommand=tree_scroll.set)
        tree_scroll.config(command=self.merge_tree.yview)
        for col, txt, w in (("Label", "Label", 240), ("Kind", "Kind", 80),
                            ("Files", "Files", 60), ("Rows", "Rows", 80)):
            self.merge_tree.heading(col, text=txt)
            self.merge_tree.column(col, width=w,
                                   anchor=("center" if col in ("Kind", "Files", "Rows") else "w"))
        self.merge_tree.pack(side="left", fill="both", expand=True)

        src_btns = tk.Frame(src_frame)
        src_btns.pack(fill="x", padx=5, pady=(0, 6))
        tk.Button(src_btns, text="Add Master CSV…",
                  command=self.merge_add_master_csv).pack(side="left", padx=4)
        tk.Button(src_btns, text="Add Experiment Folder…",
                  command=self.merge_add_folder).pack(side="left", padx=4)
        tk.Button(src_btns, text="Remove Selected",
                  command=self.merge_remove_selected).pack(side="left", padx=4)
        tk.Button(src_btns, text="Clear All",
                  command=self.merge_clear_all).pack(side="left", padx=4)

        # --- Summary preview (mirrors Tab 1's Loaded Data table + Main Plasmids) ---
        summary_frame = tk.LabelFrame(self.tab_merge, text="Sources Summary")
        summary_frame.pack(fill="both", expand=True, padx=10, pady=5)

        tk.Button(summary_frame, text="Display Summary",
                  command=self.merge_display_summary).pack(anchor="w", padx=6, pady=(6, 4))

        # Main Plasmids label (same style as Tab 1)
        self.merge_main_plasmids_label = tk.Label(summary_frame, text="", justify="left",
                                                  font=("Arial", 10, "bold"))
        self.merge_main_plasmids_label.pack(pady=(0, 5))

        sum_box = tk.Frame(summary_frame)
        sum_box.pack(pady=(0, 8), fill="both", expand=True, padx=10)
        sum_scroll = ttk.Scrollbar(sum_box, orient="vertical")
        sum_scroll.pack(side="right", fill="y")
        sum_scroll_x = ttk.Scrollbar(sum_box, orient="horizontal")
        sum_scroll_x.pack(side="bottom", fill="x")
        # Base columns; merge_display_summary appends one check column per source at runtime.
        self.merge_summary_tree = ttk.Treeview(
            sum_box, columns=("Ligand", "Cell", "Cond", "N", "Dates"),
            show="headings", yscrollcommand=sum_scroll.set,
            xscrollcommand=sum_scroll_x.set, height=6)
        sum_scroll.config(command=self.merge_summary_tree.yview)
        sum_scroll_x.config(command=self.merge_summary_tree.xview)
        for col, txt, w, anchor in (
                ("Ligand", "Ligand", 70, "w"), ("Cell", "Cell Line", 90, "w"),
                ("Cond", "Condition", 190, "w"), ("N", "N", 30, "center"),
                ("Dates", "Dates", 130, "w")):
            self.merge_summary_tree.heading(col, text=txt)
            self.merge_summary_tree.column(col, width=w, anchor=anchor)
        self.merge_summary_tree.pack(fill="both", expand=True)

        # --- Merge action ---
        merge_frame = tk.LabelFrame(self.tab_merge, text="Merge")
        merge_frame.pack(fill="x", padx=10, pady=(5, 10))
        self.btn_merge_run = tk.Button(merge_frame, text="Merge into Working Master",
                                       state="disabled", command=self.merge_run)
        self.btn_merge_run.pack(side="left", padx=6, pady=8)
        tk.Label(merge_frame,
                 text="Replaces the current working master with the merged result. "
                      "Compatibility is checked automatically (see the log).",
                 fg="#555555", justify="left").pack(side="left", padx=8)

    # --- Source-list helpers ---
    def _merge_unique_label(self, base):
        """Return `base`, or `base (2)`, `base (3)`, … so labels are unique across sources.
        merge.py keys collision suffixes and dup_keep:: resolutions off label, so duplicate
        labels would silently corrupt resolution routing."""
        existing = {s["label"] for s in self.merge_sources}
        if base not in existing:
            return base
        i = 2
        while f"{base} ({i})" in existing:
            i += 1
        return f"{base} ({i})"

    def _refresh_merge_tree(self):
        """Repopulate the Sources Treeview from self.merge_sources."""
        if self.merge_tree is None:
            return
        for item in self.merge_tree.get_children():
            self.merge_tree.delete(item)
        for s in self.merge_sources:
            self.merge_tree.insert("", "end",
                                   values=(s["label"], s["kind"], s["n_files"], s["n_rows"]))

    def _on_merge_sources_changed(self):
        """Any source-list mutation: refresh the source tree, clear the now-stale summary
        preview, and re-evaluate whether Merge is allowed (>= 2 sources). Compatibility is
        no longer pre-checked — it runs automatically inside merge_run."""
        self._refresh_merge_tree()
        if self.merge_summary_tree is not None:
            for item in self.merge_summary_tree.get_children():
                self.merge_summary_tree.delete(item)
        if self.merge_main_plasmids_label is not None:
            self.merge_main_plasmids_label.config(text="")
        self.merge_main_plasmids_choice = None
        self._update_merge_button_state()

    # --- Adding sources (never touches self.master_df) ---
    def merge_add_master_csv(self):
        """Add a master CSV as a merge source. Gates enrichment strictly: a master with any
        NaN in Donor_Raw_kinetic / Acceptor_Raw_kinetic is refused (enrich it first) — the
        same predicate merge.py uses for RAW_CHANNELS_REQUIRED, so passing here means the
        source will not trip that forbidden issue later."""
        path = filedialog.askopenfilename(
            title="Select Master CSV to add as a merge source",
            filetypes=[("CSV Files", "*.csv"), ("All Files", "*.*")])
        if not path:
            return
        try:
            df, _was_modified, _was_fixed = self._load_master_df(path)
        except Exception as e:
            self.log(f"[MERGE] Could not add '{os.path.basename(path)}': {e}")
            return

        # Channel gate (fail fast): refuse unenriched masters.
        donor = pd.to_numeric(df["Donor_Raw_kinetic"], errors="coerce")
        acceptor = pd.to_numeric(df["Acceptor_Raw_kinetic"], errors="coerce")
        if donor.isna().any() or acceptor.isna().any():
            self.log(f"[MERGE] Refused '{os.path.basename(path)}': missing raw channel "
                     f"values (Donor/Acceptor). Enrich this master first via the import "
                     f"+ enrich path, then add it.")
            messagebox.showwarning(
                "Master not enriched",
                f"'{os.path.basename(path)}' has missing raw channel values "
                f"(Donor_Raw_kinetic / Acceptor_Raw_kinetic).\n\n"
                f"Only fully enriched masters can be merged. Import it on the "
                f"Import tab and enrich it from source files first, then add it here.")
            return

        label = self._merge_unique_label(os.path.basename(path))
        self.merge_sources.append({
            "label": label, "kind": "master", "path": path, "df": df,
            "n_files": int(df['File_Name'].nunique()), "n_rows": int(len(df))})
        self.log(f"[MERGE] Added master source '{label}' "
                 f"({df['File_Name'].nunique()} files, {len(df)} rows).")
        self._on_merge_sources_changed()

    def merge_add_folder(self):
        """Add an experiment folder as a merge source. Processes it into a local
        master-shaped frame via _process_folder_to_master (no application state touched).
        Folder sources are enriched by construction."""
        directory = filedialog.askdirectory(
            title="Select experiment folder to add as a merge source")
        if not directory:
            return

        # Scan for subfolders containing xlsx/xlsm (mirror select_folder's discovery).
        subfolder_paths = []
        for root, _dirs, files in os.walk(directory):
            if any(f.endswith((".xlsx", ".xlsm")) for f in files):
                subfolder_paths.append(root)
        if not subfolder_paths:
            self.log(f"[MERGE] No .xlsx/.xlsm files found under '{directory}'.")
            messagebox.showwarning("No data found",
                                   "No .xlsx or .xlsm files found in that folder.")
            return

        # Build the processing config exactly as collect_files does (same dialogs).
        try:
            is_labeling = self.var_labeling_is_checked.get()
        except tk.TclError:
            is_labeling = False
        lum_threshold, vehicle_threshold = self._get_thresholds()
        config = ProcessingConfig(
            lum_threshold=lum_threshold,
            vehicle_warning_threshold=vehicle_threshold,
            labeling_correction=is_labeling,
            plate_layout=build_plate_layout(is_labeling),
            user_input_fn=lambda **kwargs: ask_user_parameter(self.main_gi, **kwargs),
            ligand_choice_fn=self._ask_ligand_choice_logged,
            ligand_layout_fn=self._ask_ligand_layout_logged)

        self.log(f"[MERGE] Processing folder '{os.path.basename(directory)}' "
                 f"({len(subfolder_paths)} subfolder(s))…")
        try:
            df = self._process_folder_to_master(subfolder_paths, config)
        except Exception as e:
            self.log(f"[MERGE] Could not process '{os.path.basename(directory)}': {e}")
            return

        if df is None:
            # User cancelled at the main-plasmids selection; abort without adding a source.
            return

        if df.empty:
            self.log(f"[MERGE] '{os.path.basename(directory)}' produced no master rows; "
                     f"not added.")
            return

        label = self._merge_unique_label(os.path.basename(directory) or directory)
        self.merge_sources.append({
            "label": label, "kind": "folder", "path": directory, "df": df,
            "n_files": int(df['File_Name'].nunique()), "n_rows": int(len(df))})
        self.log(f"[MERGE] Added folder source '{label}' "
                 f"({df['File_Name'].nunique()} files, {len(df)} rows).")
        self._on_merge_sources_changed()

    def merge_remove_selected(self):
        """Remove the selected source(s) from the list."""
        if self.merge_tree is None:
            return
        sel = self.merge_tree.selection()
        if not sel:
            return
        labels = {self.merge_tree.item(i, "values")[0] for i in sel}
        self.merge_sources = [s for s in self.merge_sources if s["label"] not in labels]
        self._on_merge_sources_changed()

    def merge_clear_all(self):
        """Empty the source list."""
        if not self.merge_sources:
            return
        self.merge_sources = []
        self._on_merge_sources_changed()

    # --- Sources summary preview (mirrors Tab 1, with a Source column) ---
    def _build_merge_sources(self):
        return [merge.MergeSource(label=s["label"], df=s["df"], kind=s["kind"], path=s["path"])
                for s in self.merge_sources]

    def _master_df_to_index_records(self, df, source_str=None):
        """Tab-1-style per-column index records from a master-shaped df: one record per
        (File_Name, column) that is not Empty/Unknown and not fully excluded, carrying
        Ligand/Cell_Line/Condition(=Transfection)/Date/Main_Plasmids. If source_str is given,
        each record is tagged Source=source_str. Read-only — never mutates df."""
        records = []
        if df is None or df.empty:
            return records
        work = df.copy()
        try:
            work['_DateStr'] = parse_date_series(
                work['Date'], context="merge index").dt.strftime('%d.%m.%y')
        except Exception:
            work['_DateStr'] = work['Date'].astype(str)
        work['_ColIdx'] = work['Well_ID'].astype(str).str[1:]
        work['_Excl'] = (work['Is_Excluded'].map(coerce_bool)
                         if 'Is_Excluded' in work.columns else False)
        work = work[~work['Transfection'].astype(str).str.contains("Empty", na=False)]
        work = work[~work['Cell_Line'].astype(str).str.startswith("Unknown", na=False)]
        for (fname, _col), g in work.groupby(['File_Name', '_ColIdx'], sort=False):
            if g['_Excl'].all():
                continue
            first = g.iloc[0]
            rec = {
                "File_Name": fname,
                "Date": first['_DateStr'],
                "Cell_Line": first['Cell_Line'],
                "Condition": first['Transfection'],
                "Ligand": first['Ligand'],
                "Main_Plasmids": (str(first['Main_Plasmids'])
                                  if 'Main_Plasmids' in g.columns else "Unknown"),
            }
            if source_str is not None:
                rec["Source"] = source_str
            records.append(rec)
        return records

    def _compute_merge_preview(self, sources, report):
        """Produce the merged frame for the summary preview WITHOUT touching self.master_df.

        Mirrors merge_run's conflict handling so the summary reflects the ACTUAL merged
        result: hard-forbidden (no-override) issues -> None (cannot preview); filename
        collisions are auto-resolved by 'rename' for the preview (so both copies show);
        conflicting exclusions are unioned (Is_Excluded only — enough for accurate N, no
        recompute needed for counts)."""
        hard = [i for i in report.forbidden
                if i["code"] not in (merge.FILENAME_DATA_COLLISION,
                                     merge.MULTIPLE_MAIN_PLASMIDS)]
        if hard:
            codes = ", ".join(sorted({i["code"] for i in hard}))
            self.log(f"[MERGE] Cannot preview merged summary — unresolved forbidden "
                     f"issue(s): {codes}. Resolve these before merging.")
            return None

        resolutions = {}
        if any(i["code"] == merge.FILENAME_DATA_COLLISION for i in report.forbidden):
            resolutions["filename_collision"] = "rename"
            self.log("[MERGE] Summary preview: same-named files with different data are "
                     "shown as separate (renamed) copies.")

        # MULTIPLE_MAIN_PLASMIDS: use selected main plasmids for display summary
        canon = getattr(report, "canon_plan", None)
        if canon is not None and canon.requires_selection:
            selected_mp = self.merge_main_plasmids_choice or list(canon.preselected)
            if not selected_mp:
                self.log("[MERGE] Cannot preview merged summary — no Main Plasmids chosen.")
                return None
            resolutions["main_plasmids"] = list(selected_mp)

        try:
            merged_df, _r = merge.merge_masters(sources, resolutions, log_fn=None)
        except merge.MergeForbidden as e:
            codes = ", ".join(sorted({i["code"] for i in e.issues}))
            self.log(f"[MERGE] Cannot preview merged summary: {codes}.")
            return None

        # Union conflicting exclusions (Is_Excluded only) so fully-excluded columns drop
        # from N exactly as they will after a real merge.
        file_col = merged_df['File_Name'].astype(str)
        well_col = merged_df['Well_ID'].astype(str)
        for iss in report.needs_input:
            if iss["code"] != merge.FILENAME_EXCLUSION_CONFLICT:
                continue
            ctx = iss["context"]
            fname = ctx.get("file_name")
            wells = {str(w) for lst in ctx.get("excluded_wells_by_source", {}).values()
                     for w in lst}
            for well in wells:
                m = (file_col == str(fname)) & (well_col == str(well))
                if m.any():
                    merged_df.loc[m, 'Is_Excluded'] = True
        return merged_df

    def merge_display_summary(self):
        """Run the compatibility check (logged) and show the summary of the MERGED result —
        one row per (Ligand, Cell Line, Condition) with the merged N + dates, plus one check
        column per source marking which sources contribute to that row. Conflicts are
        resolved exactly as a real merge would (dedup, collision-rename, union exclusions).
        Read-only — does not merge into or mutate self.master_df."""
        tree = self.merge_summary_tree
        if tree is not None:
            for item in tree.get_children():
                tree.delete(item)
        if self.merge_main_plasmids_label is not None:
            self.merge_main_plasmids_label.config(text="")

        if not self.merge_sources:
            self.log("[MERGE] No sources to summarize.")
            return

        sources = self._build_merge_sources()

        # Run compatibility (to the log), so the summary reflects a conflict-resolved merge.
        self.log("[MERGE] Checking compatibility ...")
        try:
            report = merge.classify_sources(sources)
        except Exception as e:
            self.log(f"[MERGE] Compatibility check failed: {e}")
            return
        self._log_merge_report(report)

        # MULTIPLE_MAIN_PLASMIDS: if Main_Plasmids selection is needed, ask here (on Display Summary)
        canon = getattr(report, "canon_plan", None)
        if canon is not None and canon.requires_selection:
            default = self.merge_main_plasmids_choice or canon.preselected
            selection = ask_main_plasmids_selection(
                self.main_gi, canon.candidate_tokens, default,
                context=("Sources use different Main Plasmids. Choose the shared plasmids "
                         "for the merged data."))
            if not selection:
                self.log("[MERGE] Summary cancelled: no Main Plasmids chosen.")
                return
            self.merge_main_plasmids_choice = selection # selection is stored and reused by merge
            self.log(f"[MERGE] Main Plasmids selected: {', '.join(selection)}.")

        merged_df = self._compute_merge_preview(sources, report)
        if merged_df is None:
            return  # hard-forbidden; reason already logged

        merged_records = self._master_df_to_index_records(merged_df)
        if not merged_records:
            self.log("[MERGE] Merged result produced no displayable rows.")
            return
        merged_idx = pd.DataFrame(merged_records)

        # Per-source contribution keys: which (Ligand, Cell_Line, Condition) each source
        # provides (from its own non-excluded data). Group-level attribution is robust to
        # dedup/rename (a deduped-away copy still counts as that source contributing).
        src_labels = [s["label"] for s in self.merge_sources]
        # Canonicalize each source the SAME way as the merged frame before building the keys
        canon_plan = (merge.CanonPlan(selected_main=list(self.merge_main_plasmids_choice))
                      if self.merge_main_plasmids_choice else None)
        source_keys = {}
        for s in self.merge_sources:
            src_df = s["df"]
            if canon_plan is not None:
                src_df = merge.apply_canonicalization(src_df, canon_plan)
            recs = self._master_df_to_index_records(src_df)
            source_keys[s["label"]] = {(r["Ligand"], r["Cell_Line"], r["Condition"])
                                       for r in recs}

        # Main Plasmids of the merged result, same presentation as Tab 1.
        mp_vals = [str(v) for v in pd.unique(merged_idx['Main_Plasmids'].dropna())
                   if str(v) not in ("", "Unknown")]
        if self.merge_main_plasmids_label is not None:
            self.merge_main_plasmids_label.config(text=", ".join(mp_vals))

        # Configure columns dynamically: fixed summary columns + one check column per source.
        fixed = [("Ligand", "Ligand", 70, "w"), ("Cell", "Cell Line", 90, "w"),
                 ("Cond", "Condition", 190, "w"), ("N", "N", 30, "center"),
                 ("Dates", "Dates", 130, "w")]
        src_col_ids = [f"src::{lbl}" for lbl in src_labels]
        tree["columns"] = [c[0] for c in fixed] + src_col_ids
        for cid, txt, w, anchor in fixed:
            tree.heading(cid, text=txt)
            tree.column(cid, width=w, anchor=anchor, stretch=False)
        for lbl, cid in zip(src_labels, src_col_ids):
            tree.heading(cid, text=lbl)
            tree.column(cid, width=110, anchor="center", stretch=False)

        # One row per merged (Ligand, Cell, Condition); N = unique merged files; per-source ✓.
        grouped = merged_idx.groupby(['Ligand', 'Cell_Line', 'Condition'])
        for (lig, cell, cond), group in grouped:
            n_count = group['File_Name'].nunique()
            date_str = ", ".join(sorted(group['Date'].unique()))
            checks = ["✓" if (lig, cell, cond) in source_keys[lbl] else ""
                      for lbl in src_labels]
            tree.insert("", "end", values=[lig, cell, cond, n_count, date_str] + checks)

        self.log(f"[MERGE] Summary displayed: merged view across "
                 f"{len(self.merge_sources)} source(s).")

    # --- Compatibility logging + Merge gate ---
    def _log_merge_report(self, report):
        """Send the compatibility analysis to the app log, by severity (never by code), with
        one deliberate exception: FILENAME_EXCLUSION_CONFLICT is no longer a user-input
        decision (it is auto-unioned, see _apply_merge_exclusion_union), so its stale
        'pick which source...' message is suppressed here — the union handler logs its own
        clear, per-file message instead."""
        for iss in report.forbidden:
            self.log(f"   [COMPAT][FORBIDDEN] {iss['message']}")
        for iss in report.needs_input:
            if iss["code"] == merge.FILENAME_EXCLUSION_CONFLICT:
                continue  # auto-unioned; see merge_run / _apply_merge_exclusion_union
            self.log(f"   [COMPAT][NEEDS INPUT] {iss['message']}")
        for iss in report.auto:
            self.log(f"   [COMPAT][INFO] {iss['message']}")
        if not report.all_issues():
            self.log("   [COMPAT] No compatibility issues detected.")
        s = report.summary or {}
        if s:
            self.log(f"   [COMPAT] Projected: {s.get('total_sources', 0)} sources, "
                     f"{s.get('total_files', 0)} files, {s.get('total_rows', 0)} rows.")

    def _update_merge_button_state(self):
        """Enable Merge iff there are >= 2 sources. Compatibility (and any forbidden block)
        is evaluated at merge time, inside merge_run, not here."""
        enable = len(self.merge_sources) >= 2
        if self.btn_merge_run is not None:
            self.btn_merge_run.config(state="normal" if enable else "disabled")

    def _clear_folder_load_state(self):
        """Clear any leftover experiment-folder display on Tab 1 (subfolders discovered by a
        previous folder selection). Called after a master import or merge swaps in a CSV/merged
        master so Tab 1 does not show stale folder candidates or a live Load Files button."""
        self.subfolder_paths_with_files = []
        self.directory = ""
        if self.subfolders_label is not None:
            self.subfolders_label.config(text="No folder selected.")
        if self.load_files_button is not None:
            self.load_files_button.config(state="disabled")

    def _apply_merge_exclusion_union(self, merged_df, conflict_union):
        """For each duplicate file with conflicting exclusions, force Is_Excluded=True on the
        UNION of wells, add matching per-well File-pinned tokens to the Applied_Exclusions blob
        so the added wells are shown AND revertable in the Exclude tab, then recompute the affected
        files via the canonical engine path. """
        file_col = merged_df['File_Name'].astype(str)
        well_col = merged_df['Well_ID'].astype(str)

        blob = "None"
        if ('Applied_Exclusions' in merged_df.columns
                and merged_df['Applied_Exclusions'].notna().any()):
            blob = str(merged_df['Applied_Exclusions'].dropna().iloc[0])
        existing_tokens = ([t.strip() for t in blob.split("||")]
                           if blob and blob != "None" else [])
        token_set = set(existing_tokens)

        new_tokens = []
        affected_files = set()
        added_wells = 0
        for fname, wells in conflict_union.items():
            for well in sorted(wells):
                row_mask = (file_col == str(fname)) & (well_col == str(well))
                if not row_mask.any():
                    continue
                already = bool(merged_df.loc[row_mask, 'Is_Excluded'].map(coerce_bool).all())
                if already:
                    continue  # keeper already excludes this well; its rule is in the blob
                merged_df.loc[row_mask, 'Is_Excluded'] = True
                affected_files.add(str(fname))
                added_wells += 1
                tok = well_token(merged_df, str(fname), str(well))
                if tok not in token_set:
                    token_set.add(tok)
                    new_tokens.append(tok)

        if new_tokens:
            merged_df['Applied_Exclusions'] = " || ".join(existing_tokens + new_tokens)

        if affected_files:
            config = self._build_recompute_config()
            merged_df = recompute_master_after_exclusion(
                merged_df, sorted(affected_files), config)
            self.log(f"   [MERGE] Applied union exclusions: {added_wells} added well(s) across "
                     f"{len(affected_files)} file(s); recomputed.")
        return merged_df

    # --- Merge action ---
    def merge_run(self):
        """Auto-check compatibility (to the log), resolve any filename collision via dialog,
        union conflicting exclusions (no prompt), call merge.merge_masters once
        (defensively), and on success swap in the merged master and run the refresh trio."""
        if len(self.merge_sources) < 2:
            self.log("[MERGE] Need at least two sources to merge.")
            return

        sources = self._build_merge_sources()

        self.log("[MERGE] Checking compatibility ...")
        try:
            report = merge.classify_sources(sources)
        except Exception as e:
            self.log(f"[MERGE] Compatibility check failed: {e}")
            return
        self._log_merge_report(report)

        resolutions = {}

        # MULTIPLE_MAIN_PLASMIDS: if user merged without viewing summary, prompt here
        canon = getattr(report, "canon_plan", None)
        if canon is not None and canon.requires_selection:
            selection = self.merge_main_plasmids_choice
            if not selection:
                selection = ask_main_plasmids_selection(
                    self.main_gi, canon.candidate_tokens, canon.preselected,
                    context=("Sources use different Main Plasmids. Choose the shared plasmids "
                             "for the merged data."))
                if not selection:
                    warn_and_abort(
                        self.main_gi, merge.MULTIPLE_MAIN_PLASMIDS,
                        "Merge cancelled: a Main Plasmids must be chosen to "
                        "reconcile the sources.")
                    self.log("[MERGE] Aborted: Main Plasmids selection cancelled.")
                    return
                self.merge_main_plasmids_choice = selection
            resolutions["main_plasmids"] = selection
            self.log(f"[MERGE] Main Plasmids: {', '.join(selection)}.")

        # FILENAME_DATA_COLLISION
        collision_issues = [i for i in report.forbidden
                            if i["code"] == merge.FILENAME_DATA_COLLISION]
        if collision_issues:
            colliding = sorted({i["context"].get("file_name")
                                for i in collision_issues if i["context"].get("file_name")})
            choice = ask_filename_collision(self.main_gi, colliding)
            if choice == "rename":
                resolutions["filename_collision"] = "rename"
            else:
                self.log("[MERGE] Aborted: filename collision left unresolved.")
                return

        # Any OTHER forbidden code (no override) -> generic warn_and_abort sink + abort.
        other_forbidden = [i for i in report.forbidden
                          if i["code"] not in (merge.FILENAME_DATA_COLLISION,
                                               merge.MULTIPLE_MAIN_PLASMIDS)]
        if other_forbidden:
            iss = other_forbidden[0]
            warn_and_abort(self.main_gi, iss["code"], iss["message"])
            self.log(f"[MERGE] Aborted: {iss['code']}.")
            return

        # FILENAME_EXCLUSION_CONFLICT
        conflict_union = {}  # file_name -> set(well_id)
        for iss in report.needs_input:
            if iss["code"] != merge.FILENAME_EXCLUSION_CONFLICT:
                continue
            ctx = iss["context"]
            fname = ctx.get("file_name")
            by_source = ctx.get("excluded_wells_by_source", {})
            wells = {str(w) for lst in by_source.values() for w in lst}
            if fname is not None and wells:
                conflict_union[fname] = wells
                self.log(f"   [MERGE] Exclusion conflict on '{fname}': applying the UNION of "
                         f"all sources' exclusions ({len(wells)} well(s): "
                         f"{', '.join(sorted(wells))}). These remain revertable in the "
                         f"Exclude tab.")

        # Call the backend exactly once, defensively.
        try:
            merged_df, report = merge.merge_masters(sources, resolutions, log_fn=self.log)
        except merge.MergeForbidden as e:
            for iss in e.issues:
                self.log(f"[MERGE BLOCKED] {iss['code']}: {iss['message']}")
            return

        # Apply the union of exclusions for conflicting duplicates (rows + blob), recompute.
        if conflict_union:
            merged_df = self._apply_merge_exclusion_union(merged_df, conflict_union)

        # Swap in the merged master + run the documented refresh trio (mirror import_master_csv)
        self.master_df = merged_df
        self.experiment = []                       # CSV/merged mode, no live objects
        self.master_index = pd.DataFrame()
        self._sync_rule_history_from_master()
        self._set_data_source(f"Merged: {len(self.merge_sources)} sources")
        self._set_csv_label(None)                  # merged result is not a CSV import
        self._clear_folder_load_state()            # clear any leftover Tab-1 folder display
        if 'Main_Plasmids' in merged_df.columns:
            mp_vals = [str(v) for v in pd.unique(merged_df['Main_Plasmids'].dropna())]
            self.main_plasmids_label.config(text=", ".join(mp_vals))
        else:
            self.main_plasmids_label.config(text="")
        self.built_master_index(source="master")
        self._update_quality_buttons_state()
        self.refresh_active_exclusions()
        self.refresh_plot_helper_options()
        self._refresh_crc_window()
        self.btn_export_master.config(state="normal")
        self.btn_export_excel.config(state="normal")

        # Concise human summary; offer to export.
        s = report.summary
        self.log(f"[MERGE] Merged {s['total_sources']} sources -> {s['total_files']} files, "
                 f"{s['total_rows']} rows "
                 f"({s['collisions_handled']} renamed, {s['duplicates_dropped']} deduped).")
        if messagebox.askyesno("Merge complete",
                               "Merged into the working master.\n\n"
                               "Export the merged master to CSV now?"):
            self.export_master_csv()

    def import_master_csv(self):
        """Loads a master csv file directly into the memory for the plot helper."""
        file_path = filedialog.askopenfilename(
            title="Select Master CSV File",
            filetypes=[("CSV Files", "*.csv"), ("All Files", "*.*")]
        )
        if not file_path: return

        try:
            self.log(f"Loading Master CSV: {os.path.basename(file_path)}")
            # Shared read + validate + schema-migrate + Is_Excluded coercion.
            df, was_modified, was_fixed = self._load_master_df(file_path)

            # Store in the unified variable
            self.master_df = df

            # Reset raw data references so we know we are in "CSV Mode"
            self.experiment = []
            self.master_index = pd.DataFrame()  # Clear exclusion index
            # The imported master carries its own Applied_Exclusions blob -> align
            # rule_history_text with it
            self._sync_rule_history_from_master()

            # Update GUI
            self._clear_folder_load_state()
            self._set_data_source(f"CSV: {os.path.basename(file_path)}")
            self._set_csv_label(os.path.basename(file_path))
            # Populate the Loaded Data "Main Plasmids" label from the imported master
            if 'Main_Plasmids' in df.columns:
                mp_vals = [str(v) for v in pd.unique(df['Main_Plasmids'].dropna())]
                self.main_plasmids_label.config(text=", ".join(mp_vals))
            else:
                self.main_plasmids_label.config(text="")
            self.refresh_plot_helper_options()

            # Build the dropdown/summary index directly from master_df (CSV mode)
            # and re-evaluate the on-demand quality-check buttons for this source.
            self.built_master_index(source="master")
            self._update_quality_buttons_state()
            self.refresh_active_exclusions()
            self._refresh_crc_window()

            # Enable exports for the imported master (both run purely off master_df)
            self.btn_export_master.config(state="normal")
            self.btn_export_excel.config(state="normal")

            # If schema was updated, offer to save and optionally enrich
            if was_modified:
                self.log(f"[MASTER UPDATED] Master was updated to current version.")
                self.show_csv_updated_dialog(was_fixed)

        except Exception as e:
            self.log(f"[ERROR] CSV Load Failed: {e}")

    def show_csv_updated_dialog(self, was_fixed = False):
        """
        Shows a dialog after importing CSV that was updated to current app version state,
        offering to save and/or enrich from source files.
        """
        dialog = tk.Toplevel(self.main_gi)
        dialog.title("Master CSV Updated")
        dialog.geometry("520x220")
        dialog.transient(self.main_gi)
        dialog.grab_set()

        txt_frame = tk.Frame(dialog)
        txt_frame.pack(pady=10, fill="both", padx=5, expand=True)

        tk.Label(txt_frame,
                 text="Imported Master CSV was Updated",
                 wraplength=480, justify="left", font=("Arial", 10, "bold")).pack(side="top")

        # Check if raw data columns are missing (all NaN)
        has_raw_data = (
                'Donor_Raw_kinetic' in self.master_df.columns
                and self.master_df['Donor_Raw_kinetic'].notna().any()
        )
        if not has_raw_data:
            tk.Label(txt_frame,
                     text="Raw data columns (Donor, Acceptor, PR Time) are empty.\n"
                          "You can enrich this CSV by pointing to the original source files.",
                     wraplength=480, justify="left", font=("Arial", 9)).pack(pady=(0, 10), padx=5)

        if was_fixed:
            tk.Label(txt_frame,
                     text="Bugs were fixed. Please save updated Master.",
                     wraplength=480, justify="left", font=("Arial", 9, "bold"),
                     fg="red").pack(side="bottom")

        btn_frame = tk.Frame(dialog)
        btn_frame.pack(pady=10, fill="x", padx=5)

        def on_save():
            dialog.destroy()
            self.export_master_csv(is_updated=True)

        def on_enrich():
            dialog.destroy()
            self.enrich_master_from_source_files()

        def on_skip():
            dialog.destroy()

        tk.Button(btn_frame, text="Save Updated CSV",
                  command=on_save).pack(side="left", fill="x", expand=True, padx=3)
        if not has_raw_data:
            tk.Button(btn_frame, text="Enrich from Source Files",
                      command=on_enrich).pack(side="left", fill="x", expand=True, padx=3)
        tk.Button(btn_frame, text="Skip",
                  command=on_skip).pack(side="left", fill="x", expand=True, padx=3)

        self.main_gi.wait_window(dialog)

    def enrich_master_from_source_files(self):
        """
        Enriches a legacy Master CSV by reloading source xlsx files to populate
        Donor_Raw_kinetic, Acceptor_Raw_kinetic, and PR_Time(min).
        Validates by matching experimental conditions and Raw_BRET_kinetic data.
        Raw channel data is stored without exclusions for full traceability.
        """
        df_old = self.master_df
        if df_old is None or df_old.empty:
            return

        # --- 1. Resolve source path ---
        stored_path = df_old['Path'].iloc[0] if 'Path' in df_old.columns else None
        source_dir = None

        if stored_path and stored_path != "undocumented path" and os.path.isdir(stored_path):
            source_dir = stored_path
            self.log(f"[ENRICH] Using stored path: {source_dir}")
        else:
            self.log("[ENRICH] Original path unavailable. Please select the experiment folder.")
            source_dir = filedialog.askdirectory(title="Select folder containing source xlsx files")

        if not source_dir:
            self.log("[ENRICH] Cancelled — no folder selected.")
            return

        # --- 2. Load protocols + measurements from path ---
        logger.debug("--- Master Enrichment ---")
        logger.debug("scanning for files...")

        # Determine plate layout from old master (check for labeling control)
        is_labeling = "labeling control" in df_old.get("Replicate", pd.Series()).astype(str).values

        # Extract baseline_end_index from old master: row index where Time_(min) == 0
        baseline_end_idx = None
        if 'Time_(min)' in df_old.columns:
            zero_times = df_old.loc[df_old['Time_(min)'] == 0.0]
            if not zero_times.empty:
                # Get the position within any file (count rows before time==0 for one well)
                sample_file = df_old['File_Name'].iloc[0]
                sample_well = df_old['Well_ID'].iloc[0]
                file_well_mask = (df_old['File_Name'] == sample_file) & (df_old['Well_ID'] == sample_well)
                file_well_times = df_old.loc[file_well_mask, 'Time_(min)'].sort_values()
                zero_idx = (file_well_times == 0.0).values.argmax()
                baseline_end_idx = int(zero_idx)
                logger.debug(f"extracted baseline_end_index={baseline_end_idx} from old master.")

        enrich_config = ProcessingConfig(
            plate_layout=build_plate_layout(is_labeling),
            labeling_correction=is_labeling,
            baseline_end_index=baseline_end_idx,
            ligand_choice_fn=self._ask_ligand_choice_logged,
            ligand_layout_fn=self._ask_ligand_layout_logged
        )

        # Get folder paths containing xlsx/xlsm files
        folder_paths = []
        for root, dirs, files in os.walk(source_dir):
            if any(f.endswith(('.xlsx', '.xlsm')) for f in files):
                folder_paths.append(root)

        # Scan and load (debug-level logging — no per-file output to user)
        all_folders = scan_and_load_folders(folder_paths)

        # Only keep folders with both protocol and results
        loaded_folders = [f for f in all_folders if f.protocol and f.results]

        if not loaded_folders:
            self.log("[ENRICH] ERROR: No valid protocol + measurement pairs found in selected folder.")
            return

        # Filter to folders matching the Main_Plasmids from the old master
        old_main_plasmids = df_old['Main_Plasmids'].iloc[0] if 'Main_Plasmids' in df_old.columns else None
        if old_main_plasmids and old_main_plasmids != "Unknown":
            matching_folders = []
            for f in loaded_folders:
                folder_mp = " + ".join(f.protocol.main_plasmids) if f.protocol.main_plasmids else "Unknown"
                if folder_mp == old_main_plasmids:
                    matching_folders.append(f)
            if matching_folders:
                logger.debug(f"filtered {len(loaded_folders)} folders to "
                             f"{len(matching_folders)} matching Main_Plasmids='{old_main_plasmids}'.")
                loaded_folders = matching_folders
            else:
                self.log(f"[ENRICH] ERROR: No folders match Main_Plasmids '{old_main_plasmids}'. "
                         f"Source files do not match this experiment.")
                return

        total_files = sum(len(f.results) for f in loaded_folders)
        self.log(f"[ENRICH] Loaded {total_files} measurement files from "
                 f"{len(loaded_folders)} folders. Extracting raw data...")

        # For cases master was created with user-input specified ligand layout:
        # --- 2b. Infer ligand layout and per-plate ligand choices from old master ---
        all_ligand_choices = {}  # file_name -> "L1" or "L2"

        for folder in loaded_folders:
            protocol = folder.protocol
            if not protocol or not protocol.ligand_2:
                continue

            layout, choices = infer_ligand_info_from_master(
                df_old, protocol, len(enrich_config.plate_layout))

            if layout:
                protocol.ligand_layout = layout
                self.log(f"   [ENRICH] Inferred ligand layout '{layout}' for "
                         f"protocol '{protocol.file_name}' from existing master")

            if choices:
                all_ligand_choices.update(choices)
                for fname, choice in choices.items():
                    lig_name = (str(protocol.ligand_2) if choice == "L2"
                                else str(protocol.ligand))
                    self.log(f"   [ENRICH] Inferred ligand '{lig_name}' for plate '{fname}' "
                             f"from existing master")

        # Build a ligand_choice_fn that looks up from inferred data,
        # falling back to the dialog for files not found in the old master
        def _enrich_ligand_choice(l1_name, l2_name, plate_info):
            if plate_info in all_ligand_choices:
                return all_ligand_choices[plate_info]
            self.log(f"   [ENRICH] Plate '{plate_info}' not found in master — asking user")
            return self._ask_ligand_choice_logged(l1_name, l2_name, plate_info)

        enrich_config.ligand_choice_fn = _enrich_ligand_choice
        # ligand_layout_fn remains as dialog fallback for protocols not in old master

        # --- 3. Map metadata and extract raw data per file ---
        new_file_data = []
        time_col = "Time (min)"
        n_ok = 0
        n_err = 0

        for folder in loaded_folders:
            protocol = folder.protocol
            main_plasmids = " + ".join(protocol.main_plasmids) if protocol.main_plasmids else "Unknown"

            for result in folder.results:
                try:
                    # Map conditions onto plate columns
                    map_plate_metadata(result, protocol, enrich_config)

                    # Calculate relative time vector for this file
                    raw_df = result.raw_bret_ratio_df.copy()
                    t_vec = calculate_relative_time(
                        raw_df[time_col], enrich_config.baseline_end_index)
                    if t_vec is None:
                        t_vec = list(range(len(raw_df)))

                    # Persist baseline index for subsequent files
                    if enrich_config.baseline_end_index is None and 0 in t_vec:
                        enrich_config.baseline_end_index = t_vec.index(0)

                    # PR Time: raw plate reader time mapped by relative time
                    raw_time_list = (list(result.raw_time)
                                     if result.raw_time is not None and len(result.raw_time) > 0 else [])
                    pr_time_map = (dict(zip(t_vec, raw_time_list))
                                   if len(raw_time_list) == len(t_vec) else {})

                    # Melt helper — uses calculated Time_(min)
                    def melt_df(df_in, val_name):
                        work = df_in.drop(columns=[time_col], errors='ignore').copy()
                        if len(work) == len(t_vec):
                            work.index = t_vec
                        work.index.name = "Time_(min)"
                        return work.reset_index().melt(
                            id_vars="Time_(min)", var_name="Well_ID", value_name=val_name)

                    # Raw BRET (for post-merge validation — no exclusions)
                    df_raw = melt_df(raw_df, "Raw_BRET_new")
                    # Donor and Acceptor (raw, no exclusions — for traceability)
                    df_donor = melt_df(result.donor_df, "Donor_Raw_kinetic")
                    df_acceptor = melt_df(result.acceptor_df, "Acceptor_Raw_kinetic")

                    merge_on = ["Time_(min)", "Well_ID"]
                    df_file = df_raw.merge(df_donor, on=merge_on, how="left") \
                                    .merge(df_acceptor, on=merge_on, how="left")

                    # PR Time
                    df_file["PR_Time(min)"] = (df_file["Time_(min)"].map(pr_time_map)
                                               if pr_time_map else float('nan'))

                    # Add condition metadata for merge
                    df_file["Date"] = result.measurement_date
                    df_file["Main_Plasmids"] = main_plasmids
                    df_file["Info_Sheet"] = str(result.info_sheet) if result.info_sheet else ""

                    meta_maps = {k: {} for k in
                                 ['Transfection', 'Cell_Line', 'Ligand']}
                    for well_id in df_file['Well_ID'].unique():
                        try:
                            c_idx = int(well_id[1:])
                            meta = result.column_metadata.get(c_idx)
                            if meta:
                                meta_maps['Transfection'][well_id] = meta.condition_name
                                meta_maps['Cell_Line'][well_id] = meta.cell_line
                                meta_maps['Ligand'][well_id] = meta.ligand_identity
                        except:
                            pass

                    for col_name, mapping in meta_maps.items():
                        df_file[col_name] = df_file['Well_ID'].map(mapping)

                    new_file_data.append(df_file)
                    n_ok += 1
                    logger.debug(f"Extracted {result.file_name}")

                except Exception as e:
                    n_err += 1
                    self.log(f"   [ERROR] Failed to extract {result.file_name}: {e}")

        if not new_file_data:
            self.log("[ENRICH] ERROR: No data could be extracted. Enrichment aborted.")
            return

        if n_err > 0:
            self.log(f"[ENRICH] Extracted {n_ok} files ({n_err} failed).")
        else:
            logger.debug(f"All {n_ok} files extracted successfully.")

        df_new = pd.concat(new_file_data, ignore_index=True)

        # --- 4. Validate condition combinations ---
        logger.debug("Validating conditions...")
        condition_cols = ['Main_Plasmids', 'Transfection', 'Cell_Line', 'Ligand']

        old_combos = set(
            df_old[condition_cols].drop_duplicates().itertuples(index=False, name=None))
        new_combos = set(
            df_new[condition_cols].drop_duplicates().itertuples(index=False, name=None))

        missing_combos = old_combos - new_combos
        if missing_combos:
            self.log(f"[ENRICH] ERROR: {len(missing_combos)} condition(s) from master CSV "
                     f"not found in source files:")
            for combo in list(missing_combos)[:5]:
                self.log(f"   Missing: {dict(zip(condition_cols, combo))}")
            self.log("[ENRICH] Enrichment aborted — source files do not match this experiment.")
            return

        logger.debug(f"All {len(old_combos)} condition combinations found.")

        # --- 5. Validate row count ---
        if len(df_old) != len(df_new):
            self.log(f"[ENRICH] ERROR: Row count mismatch — master CSV: {len(df_old)}, "
                     f"source files: {len(df_new)}. Enrichment aborted.")
            return

        logger.debug(f"Row count matches ({len(df_old)}).")

        # --- 6. Merge on condition columns + Well_ID + Time_(min) ---
        merge_cols = ['Date', 'Main_Plasmids', 'Transfection', 'Cell_Line', 'Ligand',
                      'Well_ID', 'Time_(min)']

        df_old_work = df_old.copy()

        # Normalize Date to consistent YYYY-MM-DD format
        df_old_work['Date'] = pd.to_datetime(df_old_work['Date']).dt.strftime('%Y-%m-%d')
        df_new['Date'] = pd.to_datetime(df_new['Date']).dt.strftime('%Y-%m-%d')

        # Normalize Time_(min) to float for both sides
        df_old_work['Time_(min)'] = df_old_work['Time_(min)'].astype(float)
        df_new['Time_(min)'] = df_new['Time_(min)'].astype(float)

        # Ensure enrichable columns exist in new data (safety fallback)
        for col in ENRICHABLE_COLS:
            if col not in df_new.columns:
                df_new[col] = LEGACY_COLUMN_DEFAULTS[col]["default"]

        # Select merge keys + enrichable columns + validation column from new data
        validation_col = 'Raw_BRET_new'
        enrich_cols = merge_cols + [validation_col] + ENRICHABLE_COLS
        df_enrich = df_new[enrich_cols].copy()

        # Drop old empty columns before merge
        df_old_work.drop(columns=ENRICHABLE_COLS, inplace=True, errors='ignore')

        df_merged = df_old_work.merge(df_enrich, on=merge_cols, how='left')

        # Check if merge created duplicate rows (many-to-many)
        if len(df_merged) != len(df_old_work):
            self.log(f"[ENRICH] ERROR: Merge changed row count from {len(df_old_work)} to {len(df_merged)}. "
                     f"Likely duplicate merge keys in source data. Enrichment aborted.")
            # Find which keys are duplicated
            dup_keys = df_enrich[df_enrich.duplicated(subset=merge_cols, keep=False)]
            if not dup_keys.empty:
                sample = dup_keys[merge_cols].head(3).to_dict('records')
                logger.debug(f"example duplicate keys: {sample}")
            return

        # --- 7. Post-merge validation: Raw BRET data should match ---
        # Validation logic: BRET ratio in old master (Raw_BRET_kinetic) is compared to extracted ratio from new path:
        # Difference of the two should be 0 (np.isclose defines decimal tolerance). NAs (excluded values) are ignored.
        if 'Raw_BRET_kinetic' in df_merged.columns and validation_col in df_merged.columns:
            # Compare only non-excluded rows (old master has NaN for excluded wells)
            compare_mask = df_merged['Raw_BRET_kinetic'].notna() & df_merged[validation_col].notna()
            if compare_mask.any():
                old_vals = df_merged.loc[compare_mask, 'Raw_BRET_kinetic'].astype(float).values
                new_vals = df_merged.loc[compare_mask, validation_col].astype(float).values
                diff = ~np.isclose(old_vals, new_vals, rtol=1e-8, atol=1e-8)
                n_mismatched = diff.sum()
                if n_mismatched > 0:
                    self.log(f"[ENRICH] WARNING: {n_mismatched} rows have mismatched Raw BRET values. "
                             f"Source files may not be the originals.")
                    # Debug: show mismatched rows
                    mismatch_positions = compare_mask[compare_mask].index[diff]
                    for idx in mismatch_positions[:5]:
                        row = df_merged.loc[idx]
                        logger.debug(
                            f"Mismatch row {idx}: "
                            f"Well={row.get('Well_ID')} Time={row.get('Time_(min)')} "
                            f"Date={row.get('Date')} Transf={row.get('Transfection')} "
                            f"old={row['Raw_BRET_kinetic']!r} new={row[validation_col]!r} "
                            f"delta={abs(float(row['Raw_BRET_kinetic']) - float(row[validation_col])):.15e}"
                        )
                else:
                    logger.debug("Raw BRET validation passed.")
            else:
                self.log("[ENRICH] WARNING: No overlapping non-NaN Raw BRET data to validate.")

        # Drop the temporary validation column
        df_merged.drop(columns=[validation_col], inplace=True, errors='ignore')

        # --- 8. Final result ---
        n_populated = df_merged['Donor_Raw_kinetic'].notna().sum()
        n_total = len(df_merged)
        logger.info(f"Complete. Populated {n_populated}/{n_total} rows with raw channel data.")

        if n_populated == 0:
            self.log("[ENRICH] ERROR: No data was populated after merge. Enrichment aborted.")
            return

        # Update master
        self.master_df = df_merged

        # Update NCollector version to current
        self.master_df['NCollector_version'] = APP_VERSION

        # Update path if it was missing
        if stored_path in (None, "undocumented path"):
            self.master_df['Path'] = source_dir
            logger.info(f"Updated Path column to: {source_dir}")

        # Refresh GUI and offer to save
        self.refresh_plot_helper_options()

        # Rebuild the index from the now-enriched master and re-evaluate the
        # quality-check buttons (donor data may now be present after enrichment).
        self.built_master_index(source="master")
        self._update_quality_buttons_state()
        # Enrichment may make previously non-restorable wells restorable (Donor/Acceptor
        # channels now present), so refresh the shared Active-Exclusions view.
        self.refresh_active_exclusions()

        self.export_master_csv(is_updated=True)

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
            self._set_csv_label(None)  # selecting a folder is not a CSV import
            self.refresh_active_exclusions()
            self.clear_exclusion_list()
            self.btn_export_master.config(state="disabled")
            self.btn_export_excel.config(state="disabled")


            # Log selected path
            self.log(f"Selected path: {self.directory}\n   Scanning for files...")
            self.subfolders_label.config(text="Scanning for files...")
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
                folder_names_string = "\n   ".join(folder_names)
                self.log(f"Selected path: {self.directory}\n   "
                         f"Found {count} subfolder(s) with xlsx/xlsm files:\n   {folder_names_string}")
                self.subfolders_label.config(
                    text=f"Found {count} subfolder(s) with xlsx/xlsm files:\n   " + "\n   ".join(folder_names))
                self.load_files_button.config(state="normal")
                logger.info(f"Found {count} folders: \n {folder_names_string}")
            else:
                self.log(f"[ERROR] No .xlsx or .xlsm files found in {self.directory} or any subfolder.")
                self.subfolders_label.config(text="No .xlsx or .xlsm files found in this folder.")
                self.load_files_button.config(state="disabled")
                logger.warning(f"No .xlsx or .xlsm files found starting from: {self.directory}")

    def _select_main_plasmids(self, experiment):
        """Group ``experiment`` by ``tuple(protocol.main_plasmids)`` and, when more than
        one set is present, prompt the user (same modal as before) to pick one.

        This is the reusable core shared by the Import tab and the Merge tab. It operates
        purely on the passed ``experiment`` and NEVER reads or mutates ``self.experiment``;
        only ``self.main_gi`` is used, as the dialog's parent.

        Returns:
            ``(filtered_experiment, selected_str)`` on success:
              * no set declared    -> ``(experiment, "No Common Plasmids Detected")``
              * exactly one set     -> ``(experiment, "<a + b>")`` (no prompt)
              * multiple, confirmed -> ``(folders_for_chosen_set, "<a + b>")``
            ``(None, None)`` if the user cancels/closes the selection dialog.
        """
        return resolve_selection(experiment, self._prompt_main_plasmids_choice)

    def _prompt_main_plasmids_choice(self, groups):
        """Open the modal main-plasmids picker for the given ``{tuple: [folders]}`` map.

        Returns the chosen tuple key, or ``None`` if the user cancels or closes the
        dialog. Called by ``resolve_selection`` only when more than one set exists.
        """
        dialog = tk.Toplevel(self.main_gi)
        dialog.title("Select Experiment")

        tk.Label(dialog,
                 text="Different experiment set ups detected across folders.\n"
                      "Select one to process:",
                 font=("Arial", 11, "bold")).pack(pady=10)

        selected_var = tk.StringVar()
        first_key = next(iter(groups))
        selected_var.set(str(first_key))  # Set default

        # Map the StringVar's string form back to the real tuple key.
        str_to_key_map = {}
        for key in groups:
            key_val_str = str(key)
            str_to_key_map[key_val_str] = key
            text_label = f"{' + '.join(key)} ({len(groups[key])} folders)"
            tk.Radiobutton(dialog, text=text_label, variable=selected_var,
                           value=key_val_str).pack(anchor="w", padx=20)

        # Cancel or closing the window aborts, so the merge path can bail cleanly.
        outcome = {"confirmed": False}

        def on_confirm():
            outcome["confirmed"] = True
            dialog.destroy()

        def on_cancel():
            outcome["confirmed"] = False
            dialog.destroy()

        btn_frame = tk.Frame(dialog)
        btn_frame.pack(pady=20)
        tk.Button(btn_frame, text="Confirm", command=on_confirm).pack(side="left", padx=10)
        tk.Button(btn_frame, text="Cancel", command=on_cancel).pack(side="left", padx=10)

        self.main_gi.wait_window(dialog)  # Wait until the window is closed

        if not outcome["confirmed"]:
            return None
        return str_to_key_map.get(selected_var.get(), first_key)

    def handle_main_plasmids_selection(self):
        """
        Checks main_plasmids consistency across ``self.experiment``. If multiple sets are
        found, the user selects one and ``self.experiment`` is filtered to that set.
        Returns the string representation of the selected set.

        Thin Import-tab wrapper over ``_select_main_plasmids``. The Import tab has no abort
        path: if the user cancels/closes the dialog, we preserve the historical behaviour
        and default to the first detected set (matching the old "closed window -> first key"
        fallback).
        """
        passed = self.experiment
        filtered, selected_str = self._select_main_plasmids(passed)

        if filtered is None:
            # User cancelled. Import has no abort path: default to the first set.
            groups = group_by_main_plasmids(passed)
            first_key = next(iter(groups))
            self.experiment = groups[first_key]
            logger.info(f"Keeping {len(self.experiment)} folders matching main plasmids: "
                        f"{first_key}")
            return " + ".join(first_key)

        self.experiment = filtered
        # Only an actual multi-set filter returns a new list object
        if filtered is not passed:
            logger.info(f"Keeping {len(self.experiment)} folders matching main plasmids: "
                        f"{selected_str}")
        return selected_str

    def built_master_index(self, source="auto"):
        """
        Builds the per-column metadata index that drives the Tab-2 exclusion dropdowns
        and the Tab-1 summary table.

        source:
          "object" — build from live PrResult.column_metadata. Used ONLY at initial
                     load (clean, no exclusions yet).
          "master" — build directly from the flat master_df rows. Used in CSV/import
                     mode AND after every exclusion (master_df is the single source of
                     truth for exclusion state; no Ref_Result object needed).
          "auto"   — master if experiment empty & master_df populated, else object.

        Both branches emit the same schema (File_Name, Date, Cell_Line, Condition,
        Ligand, Replicate) so refresh_filter_options / update_dropdown_options /
        update_summary_table work unchanged in both modes. The old Ref_Result /
        Column_Index / Transfection_ID columns are dropped — re-application no longer
        walks live objects.
        """
        if source == "auto":
            source = ("master" if (not self.experiment and self.master_df is not None
                                   and not self.master_df.empty) else "object")
        if source == "master":
            return self._build_index_from_master()
        return self._build_index_from_objects()

    def _finalize_index(self, records):
        """Common tail: store records, refresh dropdowns + summary table."""
        if records:
            self.master_index = pd.DataFrame(records)
            summary = self.master_index.groupby(['Cell_Line', 'Condition'])['File_Name'].nunique()
            logger.debug(f"Data Summary:\n{summary}")
            self.refresh_filter_options()
            self.update_summary_table()
        else:
            self.master_index = pd.DataFrame()
            self.update_summary_table()
            logger.warning("No valid data found")
        return self.master_index

    def _build_index_from_objects(self):
        """Object-path index build (initial load). Reads live column_metadata."""
        records = []
        rows_str = "ABCDEFGH"

        for folder in self.experiment:
            # Match the Main_Plasmids string exactly as the compile path writes it onto master_df
            proto = getattr(folder, "protocol", None)
            main_plasmids = (" + ".join(proto.main_plasmids)
                             if proto and proto.main_plasmids else "Unknown")
            for result in folder.results:
                if not result.column_metadata: continue
                if result.is_excluded: continue

                for col_idx, meta in result.column_metadata.items():
                    # Filter out empty cols
                    if meta.condition_name is None or "Empty" in meta.condition_name: continue
                    # Skip columns whose every well is excluded (N count drops). At
                    # initial load excluded_wells is empty, so nothing is skipped here.
                    all_wells_excluded = all(
                        f"{r}{col_idx}" in result.excluded_wells for r in rows_str)
                    if all_wells_excluded: continue

                    records.append({
                        "File_Name": result.file_name,
                        "Date": result.measurement_date.strftime('%d.%m.%y'),  # dropdown string
                        "Cell_Line": meta.cell_line,
                        "Condition": meta.condition_name,
                        "Ligand": meta.ligand_identity,
                        "Replicate": meta.replicate,
                        "Main_Plasmids": main_plasmids,
                    })

        result_df = self._finalize_index(records)
        # Object path also primes the plot helper (master_df may not exist yet).
        self.refresh_plot_helper_options()
        return result_df

    def _build_index_from_master(self):
        """
        Master-path index build (CSV/import mode + after every exclusion). Derives the
        same dropdown vocabulary directly from master_df. Note the column mapping:
        index "Condition" <- master_df "Transfection".
        """
        df = self.master_df
        if df is None or df.empty:
            return self._finalize_index([])

        work = df.copy()
        # Dropdown date string (consistent with the object branch '%d.%m.%y').
        work['_DateStr'] = parse_date_series(work['Date'], context="_build_index_from_master").dt.strftime('%d.%m.%y')
        work['_ColIdx'] = work['Well_ID'].astype(str).str[1:]
        work['_Excl'] = (work['Is_Excluded'].map(coerce_bool)
                         if 'Is_Excluded' in work.columns else False)

        # Drop empty/unknown columns (defensive — compile already removes them).
        work = work[~work['Transfection'].astype(str).str.contains("Empty", na=False)]
        work = work[~work['Cell_Line'].astype(str).str.startswith("Unknown", na=False)]

        records = []
        for (fname, col_idx), g in work.groupby(['File_Name', '_ColIdx'], sort=False):
            # Skip a column whose every well is currently excluded (N count drops).
            if g['_Excl'].all():
                continue
            first = g.iloc[0]
            records.append({
                "File_Name": fname,
                "Date": first['_DateStr'],
                "Cell_Line": first['Cell_Line'],
                "Condition": first['Transfection'],
                "Ligand": first['Ligand'],
                "Replicate": str(first['Replicate']),
                "Main_Plasmids": str(first['Main_Plasmids']) if 'Main_Plasmids' in g.columns else "Unknown",
            })

        return self._finalize_index(records)

    def _ask_ligand_choice_logged(self, ligand_1_name, ligand_2_name, plate_info):
        """Wrapper around ask_ligand_choice that logs the user's selection."""
        choice = ask_ligand_choice(self.main_gi, ligand_1_name, ligand_2_name, plate_info)
        if choice is not None:
            chosen_name = ligand_2_name if choice == "L2" else ligand_1_name
            self.log(f"   [LIGAND] Plate '{plate_info}': manually assigned to '{chosen_name}'")
        return choice

    def _ask_ligand_layout_logged(self, ligand_1_name, ligand_2_name, protocol_name):
        """Wrapper around ask_ligand_layout that logs the user's selection."""
        layout = ask_ligand_layout(self.main_gi, ligand_1_name, ligand_2_name, protocol_name)
        if layout is not None:
            self.log(f"   [LIGAND LAYOUT] User selected '{layout}' for '{ligand_1_name}' / "
                     f"'{ligand_2_name}' — applied to all plates of protocol '{protocol_name}'")
        return layout

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
        except tk.TclError:
            is_labeling = False
        lum_threshold, vehicle_threshold = self._get_thresholds()

        self.current_config = ProcessingConfig(
            lum_threshold=lum_threshold,
            vehicle_warning_threshold=vehicle_threshold,
            labeling_correction=is_labeling,
            plate_layout=build_plate_layout(is_labeling),
            # Method for optionally needed dialogs are stored
            user_input_fn=lambda **kwargs: ask_user_parameter(self.main_gi, **kwargs),
            ligand_choice_fn=self._ask_ligand_choice_logged,
            ligand_layout_fn=self._ask_ligand_layout_logged
        )

        # Reset exclusion state
        self.rule_history_text = ""
        self.refresh_active_exclusions()
        self.clear_exclusion_list()
        self.ignored_warnings.clear()

        self.log("\n--- Starting Data Collection ---")

        # Scan and load all folders
        self.experiment = scan_and_load_folders(
            self.subfolder_paths_with_files, log_fn=self.log
        )

        self.log(f"--- Loading Complete. Loaded {len(self.experiment)} folders. ---")

        # Log detected wavelengths once (from the first measurement file found).
        # The lum check is always performed regardless of wavelength (the user knows
        # their channels); wavelengths are logged for reference only.
        for folder in self.experiment:
            for res in folder.results:
                self.log(f"   [CHANNELS] Donor: {res.donor_wavelength} nm | "
                         f"Acceptor: {res.acceptor_wavelength} nm")
                self.log("   [CHANNELS] Luminescence check will be performed for all files.")
                break
            else:
                continue
            break

        # Call processing
        self.run_processing_pipeline()

        # Enable export
        self.btn_export_master.config(state="normal")
        self.btn_export_excel.config(state="normal")
        self.btn_run_plot_helper.config(state="normal")
        self._update_quality_buttons_state()

    def run_processing_pipeline(self):
        """
        2. Processing of raw BRET data and indexing with protocol info.

        This runs the OBJECT pipeline (process_bret_measurement) and is invoked ONLY
        at initial load. Re-application of exclusions no longer comes through here — it
        is handled entirely on master_df by apply_exclusions via the master-native
        engine. Because no exclusions exist at load, the resulting master is clean
        (Raw_BRET_kinetic complete), satisfying the engine's add-only assumption.
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

        # Built master indexing table from the live objects (clean, initial load).
        self.built_master_index(source="object")

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
        self._update_quality_buttons_state()
        # possible applied auto-warning exclusions
        self.refresh_active_exclusions()

        self.log("\n--- Processing Complete & Plot Helper Ready ---")

        # Enable Exports
        self.btn_export_master.config(state="normal")
        self.btn_export_excel.config(state="normal")

    def compile_master_dataframe(self, experiment=None, directory=None, rule_history_text=None):
        """
        Compiles all processing steps into one Master DataFrame.
        Structure: 1 row per well per timepoint.
        Means are repeated for respective technical replicates as AUCs for all timepoints.
        Empty wells (unknwon cell line or empty condition) are dropped.

        Defaults pull from application state (self.experiment / self.directory /
        self.rule_history_text) for the normal object pipeline. The Merge tab passes a
        LOCAL experiment, its source directory, and an empty rule history so a folder can
        be compiled into a master-shaped frame without mutating any application state.
        """
        experiment = self.experiment if experiment is None else experiment
        directory = self.directory if directory is None else directory
        rule_history_text = self.rule_history_text if rule_history_text is None else rule_history_text

        if not experiment:
            return None
        self.log("\n--- Building Master CSV ---")

        all_files_data = []

        for folder in experiment:
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

                # df_og = melt_df(res.raw_bret_ratio_df.drop(columns=["Time (min)"], errors='ignore'),
                #                "OG_BRET_ratio", t_vec)
                df_donor = melt_df(res.donor_df.drop(columns=["Time (min)"], errors='ignore'),
                                "Donor_Raw_kinetic", t_vec)
                df_acceptor = melt_df(res.acceptor_df.drop(columns=["Time (min)"], errors='ignore'),
                                "Acceptor_Raw_kinetic", t_vec)
                df_raw = melt_df(raw_clean, "Raw_BRET_kinetic", t_vec)
                df_lab = melt_df(res.labeling_corr_kinetic, "Lab_BRET_kinetic", t_vec)
                df_bl = melt_df(res.bl_corr_kinetic, "Bl_Corrected_BRET", t_vec)
                df_norm = melt_df(res.kinetic_df, "Veh_Norm_Kinetic", t_vec)

                # Merge on [Time_(min), Well_ID]
                merge_on = [df_donor.columns[0], "Well_ID"]

                merged_df = df_donor.merge(df_acceptor, on=merge_on, how="left") \
                    .merge(df_raw, on=merge_on, how="left") \
                    .merge(df_lab, on=merge_on, how="left") \
                    .merge(df_bl, on=merge_on, how="left") \
                    .merge(df_norm, on=merge_on, how="left")

                # --- PRISTINE, EXCLUSION-FREE RAW (uniform provenance for ALL files) ---
                # Raw_BRET_unexcluded = Acceptor / Donor,
                # Donor 0/NaN -> NaN
                # Left UNROUNDED. This column is never NaN-d by exclusion, so it is the
                # pristine source the engine restores from.
                _donor = pd.to_numeric(merged_df["Donor_Raw_kinetic"], errors="coerce")
                _acceptor = pd.to_numeric(merged_df["Acceptor_Raw_kinetic"], errors="coerce")
                merged_df["Raw_BRET_unexcluded"] = (
                    _acceptor / _donor.where((_donor != 0) & _donor.notna())
                )

                # --- MAP RAW PLATE READER TIME (file-specific, from res.raw_time) ---
                if res.raw_time is not None and len(res.raw_time) == len(t_vec):
                    raw_time_map = dict(zip(t_vec, res.raw_time))
                    merged_df["PR_Time(min)"] = merged_df["Time_(min)"].map(raw_time_map)
                else:
                    merged_df["PR_Time(min)"] = float('nan')

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
                merged_df["NCollector_version"] = APP_VERSION
                merged_df["Path"] = directory
                merged_df["Info_Sheet"] = str(res.info_sheet) if res.info_sheet else ""
                merged_df["File_Name"] = res.file_name
                merged_df["Date"] = res.measurement_date
                merged_df["Main_Plasmids"] = main_plasmids

                # Get the exclusion text (handle empty case)
                exclusion_text = rule_history_text if rule_history_text else "None"
                # v2.0.5: join rules with " || "
                exclusion_text_clean = exclusion_text.replace("\n", " || ")
                merged_df["Applied_Exclusions"] = exclusion_text_clean

                # Mark excluded wells: True if this well was in the exclusion list
                excluded_set = set(res.excluded_wells)
                merged_df["Is_Excluded"] = merged_df["Well_ID"].isin(excluded_set)

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
                            meta_lookups['Ligand_Conc'][well_id] = meta.ligand_conc.get(
                                row_char, float('nan'))
                            meta_lookups['Plate_Row'][well_id] = row_char
                            meta_lookups['Replicate'][well_id] = meta.replicate
                    except: pass

                merged_df['Transfection'] = merged_df['Well_ID'].map(meta_lookups['Transfection'])
                merged_df['Cell_Line'] = merged_df['Well_ID'].map(meta_lookups['Cell_Line'])
                merged_df['Ligand'] = merged_df['Well_ID'].map(meta_lookups['Ligand'])
                merged_df['Ligand_Conc'] = merged_df['Well_ID'].map(meta_lookups['Ligand_Conc'])
                merged_df['Plate_Row'] = merged_df['Well_ID'].map(meta_lookups['Plate_Row'])
                merged_df['Replicate'] = merged_df['Well_ID'].map(meta_lookups['Replicate'])
                merged_df['Is_Vehicle'] = (
                    (merged_df['Plate_Row'] == 'H') & (merged_df['Ligand_Conc'].isna())
                )

                # Drop empty wells (Empty transfections / Unknown cell lines)
                merged_df.drop(
                    merged_df[
                        merged_df['Transfection'].astype(str).str.contains("Empty", na=False) |
                        merged_df['Cell_Line'].astype(str).str.startswith("Unknown", na=False)
                    ].index, inplace=True
                )

        if not all_files_data:
            return pd.DataFrame()

        # Combine all files
        master_df = pd.concat(all_files_data, ignore_index=True)
        # Cleanup columns
        final_cols = [c for c in MASTER_COLUMNS if c in master_df.columns]
        return master_df[final_cols]

    def export_master_csv(self, is_updated = False):
        """Saves compiled master df to csv"""
        if self.master_df is None or self.master_df.empty:
            self.log("No data to export.")
            return

        # Native Dialog handles overwrite warning automatically
        file_path = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV File", "*.csv")],
            title=f"Save {'updated ' if is_updated else ''}Master CSV"
        )
        if not file_path: return

        try:
            self.master_df.to_csv(file_path, index=False)
            csv_name = os.path.basename(file_path)
            self.log(f"   [SUCCESS] Saved {'updated ' if is_updated else ''}Master CSV: {csv_name}")
            self._set_data_source(f"CSV: {csv_name}")
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

            # Category from config (set by plot helper), or detect from column membership for default export
            export_category = config.get('category', None)
            kinetic_types = list(DATA_TYPE_MAP["kinetic"].values())
            crc_types = list(DATA_TYPE_MAP["CRC"].values())

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
                    "Conc. Selection": [", ".join(
                        f"Row {c['row']}: Vehicle {c['ligand']}" if c.get('is_vehicle')
                        else f"Row {c['row']}: {c['conc']} log(M) {c['ligand']}"
                        for c in config.get('conc_mode', [])
                    )],
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

                    for dtype in selected_types:
                        # Determine labeling column handling
                        if dtype in ["Raw_BRET_kinetic", "Raw_BRET_CRC"] and is_labeling:
                            drop_labeling_col = False
                        else:
                            drop_labeling_col = True

                        # Determine which category this dtype belongs to
                        if export_category:
                            cat = export_category
                        elif dtype in kinetic_types:
                            cat = "kinetic"
                        elif dtype in crc_types:
                            cat = "CRC"
                        else:
                            cat = "unknown"

                        # --- KINETIC ---
                        if cat == "kinetic":
                            k_layout = config.get('conc_mode', [])
                            if not k_layout:
                                continue

                            df_kin = filter_by_conc(df_group, k_layout)
                            if df_kin.empty:
                                continue

                            df_kin = generate_header_key(df_kin, group_by, include_conc=True)
                            kin_pivot = create_clean_pivot(df_kin, "Time_(min)",
                                                           dtype, "Mean" in dtype,
                                                           drop_labeling_col)
                            kin_pivot.rename(columns={"Time_(min)": "Time (min)"}, inplace=True)

                            # Sheet Name with group_prefix (Max 31 chars)
                            base = f"{group_name}_{dtype}" if group_name else f"{dtype}"
                            kin_pivot.to_excel(writer, sheet_name=base[:31], index=True)
                            sheets_written = True

                        # --- CRC ---
                        elif cat == "CRC":
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
                                deduplicated = df_crc.drop_duplicates("Plate_Row").set_index("Plate_Row")
                                row_map = deduplicated["Ligand_Conc"].copy().astype(object)
                                # astype(object) to handle float and str (vehicle)
                                vehicle_rows = deduplicated["Is_Vehicle"].astype(bool)
                                row_map[vehicle_rows] = "Vehicle"
                                crc_pivot.index = crc_pivot.index.map(row_map)
                                crc_pivot.index.name = f"{df_crc['Ligand'].iloc[0]} (logM)"

                            else:
                                crc_pivot.index.name = "Plate Row"

                            base = f"{group_name}_AUC" if group_name else f"AUC_{dtype}"
                            crc_pivot.to_excel(writer, sheet_name=base[:31], index=True)
                            sheets_written = True

                        # --- BARGRAPH ---
                        elif cat == "bargraph":
                            conc_criteria = config.get('conc_mode', [])
                            if not conc_criteria:
                                continue

                            df_bar = filter_by_conc(df_group, conc_criteria)
                            if df_bar.empty:
                                continue

                            df_bar = generate_header_key(df_bar, group_by)

                            bar_table = create_bargraph_table(df_bar, dtype, group_by,
                                                              drop_labeling_control=drop_labeling_col)
                            if bar_table.empty:
                                continue

                            base = f"{group_name}_bargraph" if group_name else "Bargraph"
                            bar_table.to_excel(writer, sheet_name=base[:31], index=False)
                            sheets_written = True

                        # --- HEATMAP ---
                        elif cat == "heatmap":
                            conc_criteria = config.get('conc_mode', [])
                            if not conc_criteria:
                                continue

                            # Use df_subset (full dataset), NOT df_group (already split by group_by)
                            df_hm = filter_by_conc(df_subset, conc_criteria)
                            if df_hm.empty:
                                continue

                            hm_table = create_heatmap_table(df_hm, dtype, group_by,
                                                            drop_labeling_control=drop_labeling_col)
                            if hm_table.empty:
                                continue

                            base = "Heatmap"
                            hm_table.to_excel(writer, sheet_name=base[:31], index=True)
                            sheets_written = True
                            break  # Heatmap handles grouping internally, skip other groups

                    # If heatmap was written, break out of data_groups loop too
                    if export_category == "heatmap" and sheets_written:
                        break

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

        # Build conc_mode selection for default export (highest concentration: row A)
        df_row_a = self.master_df[self.master_df['Plate_Row'] == 'A']
        conc_lookup = build_row_info(df_row_a)
        conc_filters = list(conc_lookup.values())

        # Define Standard Config
        default_config = {
            'cells': 'All',
            'transfections': 'All',
            'ligands': 'All',
            'data_types': ['Kinetic_Mean', 'AUC_Mean'],
            'group_by': 'Transfection',
            'conc_mode': conc_filters
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