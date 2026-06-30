"""
crc_window.py — Tkinter + matplotlib window for the concentration-response plot.

Opened from the Exclude tab. Reads the CURRENT master via a get_master_df callback (the
app reassigns master_df on every recompute, so we must fetch it live, never cache it), so
calling refresh() after an exclude/restore redraws and refits against the new data.

UI: ligand / cell-line / transfection dropdowns on top (kept in sync), a split middle with
the embedded figure on the left and a data-preview table on the right (the Plot Helper CRC
for the selected condition, with a subtype selector AUC_Mean / Veh_Norm_AUC), and < / > arrows
that step through every condition present in the data, plus Save buttons.
"""

import logging
import tkinter as tk
from tkinter import ttk, filedialog

import pandas as pd
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import crc_plot
import export
from models import DATA_TYPE_MAP

logger = logging.getLogger("NCollector")
logging.getLogger("matplotlib").setLevel(logging.WARNING) # do not show warnings from matplotlib

def _colors_for(n):
    """n visually-distinct colors: the categorical tab10 for <=10 dates, else sampled viridis."""
    if n <= 0:
        return []
    if n <= 10:
        cmap = matplotlib.colormaps['tab10']
        return [cmap(i) for i in range(n)]
    cmap = matplotlib.colormaps['viridis']
    return [cmap(i / (n - 1)) for i in range(n)]


# Only allow AUC subtype selection
_AUC_SUBTYPES = {k: v for k, v in DATA_TYPE_MAP["CRC"].items() if "AUC" in k}

class CRCWindow:
    def __init__(self, parent, get_master_df, on_close=None, log_fn=None, export_fn=None):
        self.get_master_df = get_master_df
        self.on_close = on_close
        self.log = log_fn or (lambda m: None)
        # export_fn(file_path, master_df, config) assigned with write_excel_export by app
        # this way, shared engine with plot helper export
        self.export_fn = export_fn
        self.conditions = []
        self.index = 0

        self.win = tk.Toplevel(parent)
        self.win.title("Concentration-Response")
        self.win.geometry("1180x680")
        # Intercept the window's "X" to run cleanup (on_close) before destroying the window
        self.win.protocol("WM_DELETE_WINDOW", self._close)

        # --- Top: selectors ---
        top = tk.Frame(self.win)
        top.pack(side="top", fill="x", padx=8, pady=(8, 4))
        self.var_lig = tk.StringVar()
        self.var_cell = tk.StringVar()
        self.var_trans = tk.StringVar()
        self.cb_lig = self._add_selector(top, "Ligand:", self.var_lig, "ligand", width=16)
        self.cb_cell = self._add_selector(top, "Cell line:", self.var_cell, "cell_line", width=14)
        self.cb_trans = self._add_selector(top, "Transfection:", self.var_trans,
                                           "transfection", width=22)

        # --- Bottom bar (packed before the figure so it's always visible):
        #     Save plot (left) · condition arrows (center) · Export Table (right) ---
        bottom = tk.Frame(self.win)
        bottom.pack(side="bottom", fill="x", padx=8, pady=(4, 8))
        tk.Button(bottom, text="Save plot", command=self._save).pack(side="left")
        tk.Button(bottom, text="Export Table", command=self._export_table).pack(side="right")
        center = tk.Frame(bottom)
        center.pack(side="left", expand=True)        # expands into the middle
        tk.Button(center, text="◀", width=4,
                  command=lambda: self._step(-1)).pack(side="left")
        self.lbl_pos = tk.Label(center, text="")
        self.lbl_pos.pack(side="left", padx=10)
        tk.Button(center, text="▶", width=4,
                  command=lambda: self._step(1)).pack(side="left")

        # --- Middle: figure (left) | data preview table (right), resizable split ---
        mid = ttk.PanedWindow(self.win, orient="horizontal")
        mid.pack(side="top", fill="both", expand=True, padx=8, pady=4)

        left = tk.Frame(mid)
        self.fig = Figure(figsize=(7.0, 5.2), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=left)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        mid.add(left, weight=3)

        right = tk.Frame(mid)
        sel = tk.Frame(right)
        sel.pack(side="top", fill="x", pady=(0, 4))
        tk.Label(sel, text="Plotted data:").pack(side="left")
        default_key = next(k for k, v in _AUC_SUBTYPES.items() if v == "AUC_Mean")
        self.var_subtype = tk.StringVar(value = default_key)
        cb_sub = ttk.Combobox(sel, textvariable=self.var_subtype, state="readonly",
                              values=list(_AUC_SUBTYPES), width=60)
        cb_sub.pack(side="left", padx=(4, 0))
        cb_sub.bind("<<ComboboxSelected>>", lambda e: self._draw_preview_table())

        tvf = tk.Frame(right)
        tvf.pack(side="top", fill="both", expand=True)
        self.table = ttk.Treeview(tvf, show="headings")
        vsb = ttk.Scrollbar(tvf, orient="vertical", command=self.table.yview)
        hsb = ttk.Scrollbar(tvf, orient="horizontal", command=self.table.xview)
        self.table.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        self.table.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")
        tvf.grid_rowconfigure(0, weight=1)
        tvf.grid_columnconfigure(0, weight=1)
        mid.add(right, weight=2)

        self.refresh()

    def _add_selector(self, parent, label, var, field, width):
        tk.Label(parent, text=label).pack(side="left")
        cb = ttk.Combobox(parent, textvariable=var, state="readonly", width=width)
        cb.pack(side="left", padx=(2, 12))
        cb.bind("<<ComboboxSelected>>", lambda e, f=field: self._on_field_change(f, var.get()))
        return cb

    # --- Data / state ---
    def refresh(self):
        """Re-read the master, rebuild the condition list (preserving the current selection
        if it still exists), and redraw. Safe to call after every exclude/restore."""
        try:
            df = self.get_master_df()
        except Exception:
            df = None

        prev = (self.conditions[self.index]
                if self.conditions and 0 <= self.index < len(self.conditions) else None)
        self.conditions = crc_plot.list_crc_conditions(df)

        # Populate dropdown vocabularies from the available conditions.
        self.cb_lig['values'] = sorted({c['ligand'] for c in self.conditions})
        self.cb_cell['values'] = sorted({c['cell_line'] for c in self.conditions})
        self.cb_trans['values'] = sorted({c['transfection'] for c in self.conditions})

        if not self.conditions:
            self.index = 0
            self._draw_message("No concentration-response data in the loaded master.")
            return

        if prev is not None:
            key = (prev['ligand'], prev['cell_line'], prev['transfection'])
            found = next((i for i, c in enumerate(self.conditions)
                          if (c['ligand'], c['cell_line'], c['transfection']) == key), None)
            self.index = found if found is not None else min(self.index, len(self.conditions) - 1)
        else:
            self.index = 0
        self._draw()

    def _step(self, delta):
        if not self.conditions:
            return
        self.index = (self.index + delta) % len(self.conditions)
        self._draw()

    def _on_field_change(self, field, value):
        """Jump to the condition matching the changed dropdown, keeping the other two fields
        when such a combination exists; otherwise the first condition with the new value."""
        if not self.conditions:
            return
        cur = self.conditions[self.index]
        desired = {'ligand': cur['ligand'], 'cell_line': cur['cell_line'],
                   'transfection': cur['transfection']}
        desired[field] = value
        key = (desired['ligand'], desired['cell_line'], desired['transfection'])
        idx = next((i for i, c in enumerate(self.conditions)
                    if (c['ligand'], c['cell_line'], c['transfection']) == key), None)
        if idx is None:
            idx = next((i for i, c in enumerate(self.conditions) if c[field] == value), None)
        if idx is not None:
            self.index = idx
            self._draw()

    # --- Drawing ---
    def _draw_message(self, msg):
        self.ax.clear()
        self.ax.text(0.5, 0.5, msg, ha="center", va="center",
                     transform=self.ax.transAxes, fontsize=11, color="gray")
        self.ax.set_axis_off()
        self.lbl_pos.config(text="")
        self.canvas.draw()
        if getattr(self, "table", None) is not None:
            self._populate_table(None)

    def _draw(self):
        c = self.conditions[self.index]
        # Keep dropdowns in sync with the current condition.
        self.var_lig.set(c['ligand'])
        self.var_cell.set(c['cell_line'])
        self.var_trans.set(c['transfection'])
        self.lbl_pos.config(text=f"Condition {self.index + 1} / {len(self.conditions)}")

        data = crc_plot.extract_condition_data(
            self.get_master_df(), c['ligand'], c['cell_line'], c['transfection'])

        self.ax.clear()
        self.ax.set_axis_on()
        if data is None or not data.replicates:
            self._draw_message("No (non-excluded) data for this condition.")
            return

        # Individual biological replicates, one color per measurement date.
        colors = _colors_for(len(data.replicates))
        for rep, color in zip(data.replicates, colors):
            self.ax.scatter(rep.concs, rep.values, color=color, s=36,
                            label=rep.date_label, zorder=3)

        # Across-date mean.
        self.ax.scatter(data.concs, data.mean_values, color="black", s=28,
                        marker="v", label=f"Mean (N={data.n})", zorder=4)

        # 4PL fit to the mean (fall back to connecting the means if it can't be fit).
        fit = crc_plot.fit_four_pl(data.concs, data.mean_values)
        if fit is not None:
            _popt, xs, ys = fit
            self.ax.plot(xs, ys, color="black", linewidth=1.5, label="4PL fit", zorder=5)
        else:
            self.ax.text(0.02, 0.98, "4PL fit not available", transform=self.ax.transAxes,
                         va="top", ha="left", fontsize=8, color="black")

        main = c['main_plasmids']
        title = f"{main} + {c['transfection']}" if main else c['transfection']
        title += f"\n{c['cell_line']}"
        self.ax.set_title(title, fontsize=11)
        self.ax.set_xlabel(f"{c['ligand']}  [log₁₀(M)]")
        self.ax.set_ylabel("Vehicle-normalised AUC (mean of techn. replicates)")
        self.ax.legend(fontsize=8, loc="best")
        self.ax.margins(y=0.08)            # Y autoscales to the data with a little headroom
        self.fig.tight_layout()
        self.canvas.draw()

        self._draw_preview_table()

    # --- Preview table (Plot Helper CRC "second sheet") ---
    def _draw_preview_table(self):
        """Rebuild the right-hand preview table for the current condition + selected subtype.
        Reuses the export builder (build_crc_preview_table), reading the live master so it
        reflects the current exclusion state."""
        if not self.conditions:
            self._populate_table(None)
            return
        c = self.conditions[self.index]
        value_col = _AUC_SUBTYPES.get(self.var_subtype.get(), "AUC_Mean")
        config = {"ligands": [c['ligand']], "cells": [c['cell_line']],
                  "transfections": [c['transfection']]}
        try:
            df_sub = export.apply_export_filters(self.get_master_df(), config)
            pivot = export.build_crc_preview_table(df_sub, value_col)
        except Exception as e:
            logger.debug(f"CRC preview table failed: {e}")
            pivot = None
        self._populate_table(pivot)

    def _populate_table(self, pivot):
        tv = self.table
        tv.delete(*tv.get_children())
        if pivot is None or pivot.empty or len(pivot.columns) == 0:
            tv["columns"] = ("_msg",)
            tv.heading("_msg", text="")
            tv.column("_msg", width=160, anchor="w")
            tv.insert("", "end", values=("No data",))
            return
        col_ids = ["_idx"] + [f"c{i}" for i in range(len(pivot.columns))]
        tv["columns"] = col_ids
        tv.heading("_idx", text=str(pivot.index.name or "Row"))
        tv.column("_idx", width=120, anchor="w", stretch=False)
        for i, col in enumerate(pivot.columns):
            cid = f"c{i}"
            tv.heading(cid, text=str(col))
            tv.column(cid, width=95, anchor="center", stretch=False)
        for idx, row in pivot.iterrows():
            cells = ["" if (isinstance(idx, float) and pd.isna(idx)) else idx]
            cells += ["" if pd.isna(v) else f"{v:.3f}" for v in row]
            tv.insert("", "end", values=cells)

    # --- Actions ---
    def _save(self):
        if not self.conditions:
            return
        path = filedialog.asksaveasfilename(
            parent=self.win, defaultextension=".png",
            filetypes=[("PNG image", "*.png"), ("PDF", "*.pdf"), ("SVG", "*.svg")],
            title="Save concentration-response plot")
        if not path:
            return
        try:
            self.fig.savefig(path, dpi=200, bbox_inches="tight")
            self.log(f"[CRC] Saved plot: {path}")
        except Exception as e:
            self.log(f"[CRC] Could not save plot: {e}")

    def _export_table(self):
        """Run the Plot-Helper export for the current condition, category CRC, and the
        selected subtype — i.e. the same workbook the Plot Helper tab would produce."""
        if not self.conditions or self.export_fn is None:
            return
        c = self.conditions[self.index]
        value_col = _AUC_SUBTYPES.get(self.var_subtype.get(), "AUC_Mean")
        config = {
            "category": "CRC",
            "cells": [c['cell_line']],
            "transfections": [c['transfection']],
            "ligands": [c['ligand']],
            "data_types": [value_col],
            "group_by": "None",
            "conc_mode": [],
        }
        path = filedialog.asksaveasfilename(
            parent=self.win, defaultextension=".xlsx",
            filetypes=[("Excel", "*.xlsx")], title="Export CRC table")
        if not path:
            return
        # write_excel_export logs its own success/failure.
        self.export_fn(path, self.get_master_df(), config)

    def _close(self):
        if self.on_close:
            self.on_close()
        self.win.destroy()
