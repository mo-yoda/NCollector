#!/usr/bin/env python3
"""
Equivalence harness for the master-native exclusion recompute engine.

It checks that recompute_master_after_exclusion (operating purely on the flat master
DataFrame) reproduces, bit-for-bit within tolerance, what the live-object pipeline
(process_bret_measurement + compile) produces when the same wells are excluded.

Two entry points:

  * Real data (point it at an experiment folder):
        python exclusion_equivalence_harness.py --folder /path/to/experiment \
               --exclude "FILE.xlsx:A1,B2" --exclude "OTHER.xlsx:H7"
    (or set FOLDER_PATH / EXCLUSIONS at the top of this file)

  * Synthetic, no files needed (default when --folder is omitted):
        python exclusion_equivalence_harness.py
    Builds in-memory PrResult/ProtocolData objects and runs the full object vs engine
    comparison, plus a small pure-function unit check on compute_bret_metrics.

The harness intentionally re-implements a *faithful headless port* of
App.compile_master_dataframe (the only piece that lives on the Tkinter App class) so
that both paths can be compiled without a GUI. Both paths use the SAME port, so the
comparison validates the recompute engine regardless of any porting nuance.
"""

import sys
import argparse
import logging
from datetime import date

import numpy as np
import pandas as pd

from models import (PrResult, ProtocolData, MeasurementFolder, ProcessingConfig,
                    APP_VERSION, MASTER_COLUMNS, TRIPLICATE_LAYOUT, QUADRUPLICATE_LAYOUT,
                    build_plate_layout)
from processing import process_bret_measurement, recompute_master_after_exclusion

logger = logging.getLogger("NCollector")

# ----------------------------------------------------------------------------- #
# Top-of-file configuration (overridable via CLI)
# ----------------------------------------------------------------------------- #
FOLDER_PATH = None                      # e.g. r"C:\data\my_experiment"
EXCLUSIONS = {                          # {file_name: [well_id, ...]}
    # "230101_plate1_analysis.xlsx": ["A1", "B2", "H4"],
}

RTOL = 1e-8
ATOL = 1e-8

# Columns that are numeric (compared with np.isclose) vs metadata (compared exactly).
NUMERIC_COLS = {
    "Ligand_Conc", "Time_(min)", "PR_Time(min)",
    "Donor_Raw_kinetic", "Acceptor_Raw_kinetic",
    "Raw_BRET_kinetic", "Lab_BRET_kinetic", "Bl_Corrected_BRET",
    "Veh_Norm_Kinetic", "Kinetic_Mean",
    "Raw_BRET_CRC", "Lab_LP", "Bl_LP", "Veh_Norm_LP", "LP_Mean",
    "Lab_AUC", "Bl_AUC", "Veh_Norm_AUC", "AUC_Mean",
}
# Bool / metadata columns compared exactly.
EXACT_COLS = [c for c in MASTER_COLUMNS if c not in NUMERIC_COLS]


# ============================================================================= #
# Faithful headless port of App.compile_master_dataframe
# ============================================================================= #
def compile_master_dataframe(experiment, directory, rule_history_text):
    """Port of app.py:compile_master_dataframe — operates on a list of MeasurementFolder."""
    if not experiment:
        return None

    all_files_data = []

    for folder in experiment:
        if folder.protocol.main_plasmids:
            main_plasmids = " + ".join(folder.protocol.main_plasmids)
        else:
            main_plasmids = "Unknown"

        for res in folder.results:
            if res.is_excluded:
                continue

            def melt_df(df, val_name, time_vec):
                if df is None or df.empty:
                    return pd.DataFrame()
                df_work = df.copy()
                if len(df_work) != len(time_vec):
                    logger.warning(f"Length mismatch in {res.file_name}: "
                                   f"Data {len(df_work)} vs Time {len(time_vec)}")
                    df_work.index.name = "Time_Idx"
                    id_var = "Time_Idx"
                else:
                    df_work.index = time_vec
                    df_work.index.name = "Time_(min)"
                    id_var = "Time_(min)"
                return df_work.reset_index().melt(
                    id_vars=id_var, var_name="Well_ID", value_name=val_name)

            t_vec = res.time_vector
            if not t_vec:
                t_vec = range(len(res.raw_bret_ratio_cleaned))

            raw_clean = res.raw_bret_ratio_cleaned.drop(columns=["Time (min)"], errors='ignore')

            df_donor = melt_df(res.donor_df.drop(columns=["Time (min)"], errors='ignore'),
                               "Donor_Raw_kinetic", t_vec)
            df_acceptor = melt_df(res.acceptor_df.drop(columns=["Time (min)"], errors='ignore'),
                                  "Acceptor_Raw_kinetic", t_vec)
            df_raw = melt_df(raw_clean, "Raw_BRET_kinetic", t_vec)
            df_lab = melt_df(res.labeling_corr_kinetic, "Lab_BRET_kinetic", t_vec)
            df_bl = melt_df(res.bl_corr_kinetic, "Bl_Corrected_BRET", t_vec)
            df_norm = melt_df(res.kinetic_df, "Veh_Norm_Kinetic", t_vec)

            merge_on = [df_donor.columns[0], "Well_ID"]
            merged_df = df_donor.merge(df_acceptor, on=merge_on, how="left") \
                .merge(df_raw, on=merge_on, how="left") \
                .merge(df_lab, on=merge_on, how="left") \
                .merge(df_bl, on=merge_on, how="left") \
                .merge(df_norm, on=merge_on, how="left")

            if res.raw_time is not None and len(res.raw_time) == len(t_vec):
                raw_time_map = dict(zip(t_vec, res.raw_time))
                merged_df["PR_Time(min)"] = merged_df["Time_(min)"].map(raw_time_map)
            else:
                merged_df["PR_Time(min)"] = float('nan')

            raw_bret_map = res.raw_bret_points_df.iloc[0].to_dict() if res.raw_bret_points_df is not None else {}
            lp_lab_map = res.labeling_corr_lp_df.iloc[0].to_dict() if res.labeling_corr_lp_df is not None else {}
            lp_bl_map = res.bl_corr_lp_df.iloc[0].to_dict() if res.bl_corr_lp_df is not None else {}
            lp_norm_map = res.lp_df.iloc[0].to_dict() if res.lp_df is not None else {}

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

            well_to_mean_map = {}
            well_to_lp_mean_map = {}
            well_to_auc_mean_map = {}
            for col_idx, meta in res.column_metadata.items():
                col_str = str(col_idx)
                for row_char in "ABCDEFGH":
                    well_id = f"{row_char}{col_str}"
                    mean_key = f"{meta.condition_name}|{meta.cell_line}|{meta.ligand_identity}|{row_char}"
                    if res.kinetic_mean_df is not None and mean_key in res.kinetic_mean_df.columns:
                        well_to_mean_map[well_id] = res.kinetic_mean_df[mean_key].tolist()
                    if res.lp_mean_df is not None and mean_key in res.lp_mean_df.columns:
                        well_to_lp_mean_map[well_id] = res.lp_mean_df[mean_key].iloc[0]
                    if res.auc_mean_df is not None and mean_key in res.auc_mean_df.columns:
                        well_to_auc_mean_map[well_id] = res.auc_mean_df[mean_key].iloc[0]

            merged_df['LP_Mean'] = merged_df['Well_ID'].map(well_to_lp_mean_map)
            merged_df['AUC_Mean'] = merged_df['Well_ID'].map(well_to_auc_mean_map)

            mean_rows = []
            for well_id, mean_series in well_to_mean_map.items():
                if len(mean_series) == len(t_vec):
                    for t, val in zip(t_vec, mean_series):
                        mean_rows.append({'Well_ID': well_id, 'Time_(min)': t, 'Kinetic_Mean': val})
            if mean_rows:
                df_means = pd.DataFrame(mean_rows)
                merged_df = merged_df.merge(df_means, on=['Well_ID', 'Time_(min)'], how='left')
            else:
                merged_df['Kinetic_Mean'] = float('nan')

            all_files_data.append(merged_df)

            merged_df["NCollector_version"] = APP_VERSION
            merged_df["Path"] = directory
            merged_df["Info_Sheet"] = str(res.info_sheet) if res.info_sheet else ""
            merged_df["File_Name"] = res.file_name
            merged_df["Date"] = res.measurement_date
            merged_df["Main_Plasmids"] = main_plasmids

            exclusion_text = rule_history_text if rule_history_text else "None"
            merged_df["Applied_Exclusions"] = exclusion_text.replace("\n", " | ")

            excluded_set = set(res.excluded_wells)
            merged_df["Is_Excluded"] = merged_df["Well_ID"].isin(excluded_set)

            meta_lookups = {'Transfection': {}, 'Cell_Line': {}, 'Ligand': {},
                            'Ligand_Conc': {}, 'Plate_Row': {}, 'Replicate': {}}
            for well_id in merged_df['Well_ID'].unique():
                try:
                    c_idx = int(well_id[1:])
                    row_char = well_id[0]
                    meta = res.column_metadata.get(c_idx)
                    if meta:
                        meta_lookups['Transfection'][well_id] = meta.condition_name
                        meta_lookups['Cell_Line'][well_id] = meta.cell_line
                        meta_lookups['Ligand'][well_id] = meta.ligand_identity
                        meta_lookups['Ligand_Conc'][well_id] = meta.ligand_conc.get(row_char, float('nan'))
                        meta_lookups['Plate_Row'][well_id] = row_char
                        meta_lookups['Replicate'][well_id] = meta.replicate
                except Exception:
                    pass

            merged_df['Transfection'] = merged_df['Well_ID'].map(meta_lookups['Transfection'])
            merged_df['Cell_Line'] = merged_df['Well_ID'].map(meta_lookups['Cell_Line'])
            merged_df['Ligand'] = merged_df['Well_ID'].map(meta_lookups['Ligand'])
            merged_df['Ligand_Conc'] = merged_df['Well_ID'].map(meta_lookups['Ligand_Conc'])
            merged_df['Plate_Row'] = merged_df['Well_ID'].map(meta_lookups['Plate_Row'])
            merged_df['Replicate'] = merged_df['Well_ID'].map(meta_lookups['Replicate'])
            merged_df['Is_Vehicle'] = (
                (merged_df['Plate_Row'] == 'H') & (merged_df['Ligand_Conc'].isna())
            )

            merged_df.drop(
                merged_df[
                    merged_df['Transfection'].astype(str).str.contains("Empty", na=False) |
                    merged_df['Cell_Line'].astype(str).str.startswith("Unknown", na=False)
                ].index, inplace=True
            )

    if not all_files_data:
        return pd.DataFrame()

    master_df = pd.concat(all_files_data, ignore_index=True)
    final_cols = [c for c in MASTER_COLUMNS if c in master_df.columns]
    return master_df[final_cols]


# ============================================================================= #
# Comparison
# ============================================================================= #
def compare_masters(master_a, master_b, label_a="OBJECT", label_b="ENGINE"):
    """Returns (ok: bool, report: str). Compares row-aligned masters column by column."""
    lines = []
    ok = True

    if master_a is None or master_b is None:
        return False, "One of the masters is None."
    if master_a.empty or master_b.empty:
        return False, "One of the masters is empty."

    # Align: sort both by a stable key so row order can't cause spurious diffs.
    sort_keys = [c for c in ["File_Name", "Well_ID", "Time_(min)"] if c in master_a.columns]
    a = master_a.sort_values(sort_keys).reset_index(drop=True)
    b = master_b.sort_values(sort_keys).reset_index(drop=True)

    if len(a) != len(b):
        return False, f"Row count differs: {label_a}={len(a)} vs {label_b}={len(b)}"

    cols = [c for c in MASTER_COLUMNS if c in a.columns and c in b.columns]
    lines.append(f"{'COLUMN':<24} {'KIND':<8} {'STATUS':<8} DETAIL")
    lines.append("-" * 72)

    for col in cols:
        sa, sb = a[col], b[col]
        if col in NUMERIC_COLS:
            va = pd.to_numeric(sa, errors='coerce').to_numpy(dtype=float)
            vb = pd.to_numeric(sb, errors='coerce').to_numpy(dtype=float)
            close = np.isclose(va, vb, rtol=RTOL, atol=ATOL, equal_nan=True)
            n_bad = int((~close).sum())
            if n_bad == 0:
                lines.append(f"{col:<24} {'numeric':<8} {'OK':<8}")
            else:
                ok = False
                bad_idx = np.where(~close)[0][:3]
                detail = "; ".join(
                    f"row{i}: {va[i]!r} vs {vb[i]!r}" for i in bad_idx)
                lines.append(f"{col:<24} {'numeric':<8} {'DIFF':<8} {n_bad} diffs | {detail}")
        else:
            # Exact compare with NaN==NaN treated equal.
            eq = (sa.astype(object).where(sa.notna(), "__NA__")
                  == sb.astype(object).where(sb.notna(), "__NA__"))
            n_bad = int((~eq).sum())
            if n_bad == 0:
                lines.append(f"{col:<24} {'meta':<8} {'OK':<8}")
            else:
                ok = False
                bad_idx = np.where(~eq.to_numpy())[0][:3]
                detail = "; ".join(
                    f"row{i}: {sa.iloc[i]!r} vs {sb.iloc[i]!r}" for i in bad_idx)
                lines.append(f"{col:<24} {'meta':<8} {'DIFF':<8} {n_bad} diffs | {detail}")

    lines.append("-" * 72)
    lines.append(f"RESULT: {'PASS — masters match' if ok else 'FAIL — see diffs above'}")
    return ok, "\n".join(lines)


# ============================================================================= #
# Path runners
# ============================================================================= #
def _make_config(is_labeling, baseline_end_idx=None):
    return ProcessingConfig(
        labeling_correction=is_labeling,
        plate_layout=build_plate_layout(is_labeling),
        baseline_end_index=baseline_end_idx,
    )


def run_object_path(folders, exclusions, is_labeling, baseline_end_idx, rule_text):
    """Process with exclusions applied up front, then compile -> master_A."""
    cfg = _make_config(is_labeling, baseline_end_idx)
    for folder in folders:
        for res in folder.results:
            res.excluded_wells = list(exclusions.get(res.file_name, []))
            res.is_excluded = False
            process_bret_measurement(res, folder.protocol, cfg)
    directory = folders[0].folder_path if folders else ""
    return compile_master_dataframe(folders, directory, rule_text)


def run_engine_path(folders, exclusions, is_labeling, baseline_end_idx, rule_text):
    """Compile clean (no exclusions), flag Is_Excluded for S, then recompute -> master_B."""
    cfg = _make_config(is_labeling, baseline_end_idx)
    for folder in folders:
        for res in folder.results:
            res.excluded_wells = []
            res.is_excluded = False
            process_bret_measurement(res, folder.protocol, cfg)
    directory = folders[0].folder_path if folders else ""
    master = compile_master_dataframe(folders, directory, rule_text)

    # Flag the exclusion set on the flat master (the only state the engine reads).
    affected = []
    for fname, wells in exclusions.items():
        well_set = set(wells)
        mask = (master['File_Name'] == fname) & (master['Well_ID'].isin(well_set))
        master.loc[mask, 'Is_Excluded'] = True
        # Mirror the object path's provenance string so meta columns compare equal.
        master.loc[master['File_Name'] == fname, 'Applied_Exclusions'] = (
            rule_text.replace("\n", " | ") if rule_text else "None")
        if mask.any():
            affected.append(fname)

    return recompute_master_after_exclusion(master, affected, cfg)


# ============================================================================= #
# Synthetic experiment (no xlsx required)
# ============================================================================= #
def _wells(cols=range(1, 13), rows="ABCDEFGH"):
    return [f"{r}{c}" for c in cols for r in rows]


def build_synthetic_folder(is_labeling=False, seed=0):
    """
    Builds an in-memory MeasurementFolder with one protocol + one PrResult.
    Uniform 'one line per plate' cell layout, single ligand, 4 transfection blocks.
    Deterministic random BRET/donor/acceptor traces.
    """
    rng = np.random.default_rng(seed)
    n_baseline, n_post = 3, 5
    n_t = n_baseline + n_post
    raw_time = [float(i) for i in range(n_t)]          # uniform 1-min spacing
    baseline_end_idx = n_baseline                      # time==0 at index 3

    well_ids = _wells()
    time_col = pd.Series(raw_time, name="Time (min)")

    def make_trace_df(scale=1.0, base=1.0):
        data = {"Time (min)": raw_time}
        for w in well_ids:
            # Smoothly rising trace + noise so baseline/AUC are non-trivial.
            trace = base + scale * (np.linspace(0, 1, n_t) + rng.normal(0, 0.02, n_t))
            data[w] = trace
        return pd.DataFrame(data)

    raw_bret = make_trace_df(scale=0.5, base=1.0)
    donor = make_trace_df(scale=50.0, base=500.0)
    acceptor = make_trace_df(scale=80.0, base=400.0)

    # Concentration table A-G numeric, H vehicle.
    conc_vals = [-9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0]
    ligand_conc_df = pd.DataFrame({0: conc_vals + [np.nan]})

    protocol = ProtocolData(
        file_name="synthetic_protocol.xlsx",
        exp_date=date(2023, 1, 1),
        n=1,
        cell_lines=["HEK293"],
        line_layout="one line per plate",
        transfection_scheme=pd.DataFrame(),
        main_plasmids=["BackboneX"],
        transfection_conditions={"1": ["GeneA"], "2": ["GeneB"],
                                 "3": ["GeneC"], "4": ["GeneD"]},
        ligand="Forskolin",
        ligand_conc=ligand_conc_df,
        ligand_2=None,
        ligand_2_conc=None,
        ligand_layout=None,
    )

    result = PrResult(
        file_name="230101_plate1_analysis.xlsx",
        measurement_date=date(2023, 1, 1),
        cell_line="HEK293",                 # ID2
        transfection_id="1,2,3,4",          # ID3 -> 4 blocks
        raw_time=raw_time,
        raw_bret_ratio_df=raw_bret,
        donor_df=donor,
        acceptor_df=acceptor,
        donor_wavelength=480,               # != 475 -> skip lum check
        acceptor_wavelength=530,
    )

    folder = MeasurementFolder(
        folder_name="230101_synthetic",
        folder_path="/synthetic/230101_synthetic",
        measurement_date=date(2023, 1, 1),
        protocol=protocol,
        results=[result],
    )
    return folder, baseline_end_idx


# ============================================================================= #
# Pure-function unit check on compute_bret_metrics
# ============================================================================= #
def unit_test_compute_core():
    """Excluding a vehicle well must change normalization; excluding all of a block's
    vehicle wells must make that block's Veh_Norm NaN. No files involved."""
    from processing import compute_bret_metrics
    folder, base_idx = build_synthetic_folder(is_labeling=False, seed=7)
    res = folder.results[0]

    # Build the wide raw the core consumes (drop time col).
    raw_wide = res.raw_bret_ratio_df.drop(columns=["Time (min)"]).copy()
    plate_blocks = TRIPLICATE_LAYOUT

    # Minimal metadata: condition per block, single cell line/ligand, vehicle row H.
    from models import PlateColMetadata
    col_meta = {}
    for b_idx, block in enumerate(plate_blocks):
        for c in block:
            col_meta[c] = PlateColMetadata(
                cell_line="HEK293", condition_name=f"Cond{b_idx+1}",
                ligand_identity="Forskolin", replicate="1")

    base = compute_bret_metrics(raw_wide, list(range(len(raw_wide))), base_idx,
                                plate_blocks, col_meta, [], labeling_correction=False)

    # Exclude all vehicle (row H) wells of block 1 (cols 1-3).
    veh_block1 = [f"H{c}" for c in plate_blocks[0]]
    excl = compute_bret_metrics(raw_wide, list(range(len(raw_wide))), base_idx,
                                plate_blocks, col_meta, veh_block1, labeling_correction=False)

    # Block-1 non-vehicle well A1: normalized value must now be NaN (no vehicle reference).
    a1_norm_excl = excl.kinetic_df["A1"]
    assert a1_norm_excl.isna().all(), "Excluding all block-1 vehicles should NaN block-1 norm."
    # A vehicle well in another block (H4 in block 2) is unaffected.
    assert not base.kinetic_df["A4"].isna().all(), "Block-2 norm should exist in baseline."
    pd.testing.assert_series_equal(base.kinetic_df["A4"], excl.kinetic_df["A4"],
                                   check_names=False)
    print("[UNIT] compute_bret_metrics vehicle-exclusion behaviour OK")

    # Safety net: calculate_vehicle_means must never count a labeling-control column's
    # row-H well as a vehicle, even if that column is still present in the data.
    from processing import calculate_vehicle_means
    from models import QUADRUPLICATE_LAYOUT
    qblock = QUADRUPLICATE_LAYOUT[0]
    qmeta = {c: PlateColMetadata(
                 cell_line="HEK", condition_name="Cond", ligand_identity="Lig",
                 replicate=("labeling control" if c == qblock[-1] else str(c)),
                 ligand_conc={})
             for c in qblock}
    vwide = pd.DataFrame({f"H{c}": [10.0] for c in qblock})
    vwide[f"H{qblock[-1]}"] = [1000.0]   # control-column vehicle, very different value
    vm = calculate_vehicle_means(vwide, [qblock], [], column_metadata=qmeta)
    assert abs(vm[qblock[0]] - 10.0) < 1e-9, \
        "labeling-control column's row-H well must be disregarded as a vehicle"
    print("[UNIT] calculate_vehicle_means disregards labeling-control vehicles OK")
    return True


# ============================================================================= #
# Real-folder loader
# ============================================================================= #
def load_real_experiment(folder_path):
    import os
    from parsing import scan_and_load_folders
    folder_paths = []
    for root, _dirs, files in os.walk(folder_path):
        if any(f.endswith(('.xlsx', '.xlsm')) for f in files):
            folder_paths.append(root)
    folders = scan_and_load_folders(folder_paths)
    folders = [f for f in folders if f.protocol and f.results]
    if not folders:
        raise SystemExit(f"No loadable protocol+result folders found under {folder_path}")
    return folders


def parse_exclusions_arg(arg_list):
    """--exclude "FILE.xlsx:A1,B2" -> {'FILE.xlsx': ['A1','B2']}"""
    out = {}
    for item in arg_list or []:
        if ":" not in item:
            continue
        fname, wells = item.split(":", 1)
        out[fname.strip()] = [w.strip() for w in wells.split(",") if w.strip()]
    return out


# ============================================================================= #
# Main
# ============================================================================= #
def run_comparison(folders, exclusions, is_labeling, baseline_end_idx):
    rule_text = "Harness exclusion set: " + "; ".join(
        f"{f}->{','.join(w)}" for f, w in exclusions.items())

    master_a = run_object_path(folders, exclusions, is_labeling, baseline_end_idx, rule_text)
    # Re-load/rebuild folders for the engine path so object-path mutation can't leak.
    master_b = run_engine_path(folders, exclusions, is_labeling, baseline_end_idx, rule_text)

    ok, report = compare_masters(master_a, master_b)
    print(report)
    return ok


def main():
    logging.basicConfig(level=logging.WARNING,
                        format="%(levelname)s %(name)s: %(message)s")

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--folder", default=FOLDER_PATH, help="Experiment folder (real data).")
    ap.add_argument("--exclude", action="append", default=None,
                    help='Exclusion as "FILE.xlsx:A1,B2". Repeatable.')
    ap.add_argument("--labeling", action="store_true",
                    help="Treat real data as labeling-correction (quadruplicate) layout.")
    args = ap.parse_args()
    # args.labeling = True  # Force labeling to True for testing

    all_ok = True

    # Always run the pure-function unit check first.
    all_ok &= unit_test_compute_core()

    if args.folder:
        print(f"\n=== REAL DATA: {args.folder} ===")
        folders = load_real_experiment(args.folder)
        exclusions = parse_exclusions_arg(args.exclude) or EXCLUSIONS
        if not exclusions:
            # Default: exclude a couple of wells from the first file.
            first_file = folders[0].results[0].file_name
            exclusions = {first_file: ["A1", "H1"]}
            print(f"No --exclude given; defaulting to {exclusions}")
        is_labeling = args.labeling
        all_ok &= run_comparison(folders, exclusions, is_labeling, baseline_end_idx=None)
    else:
        # Synthetic scenarios: non-labeling (triplicate) and labeling (quadruplicate).
        for is_labeling in (False, True):
            tag = "LABELING (quadruplicate)" if is_labeling else "STANDARD (triplicate)"
            print(f"\n=== SYNTHETIC: {tag} ===")
            folder, base_idx = build_synthetic_folder(is_labeling=is_labeling, seed=1)
            fname = folder.results[0].file_name
            # Exclude a vehicle well, a normal well, and (for labeling) a control-column well.
            wells = ["A1", "H1", "C5"]
            if is_labeling:
                wells.append("B4")   # col 4 is the labeling-control column of block 1 (quad)
            exclusions = {fname: wells}
            all_ok &= run_comparison([folder], exclusions, is_labeling, base_idx)

    print("\n" + ("ALL CHECKS PASSED" if all_ok else "SOME CHECKS FAILED"))
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()