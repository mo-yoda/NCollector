"""
Headless CSV-mode simulation. Exercises the REAL refactored methods (bound to a
lightweight fake `self`) plus the shared reconstruction helper, with NO Tkinter GUI
and NO source xlsx — proving the unified master-native exclusion path works in
CSV/import mode.
"""
import sys, types
sys.path.insert(0, "_tkstub")

import numpy as np
import pandas as pd
from types import SimpleNamespace

import app
from app import NCollectorApp
from processing import reconstruct_file_inputs, recompute_master_after_exclusion, coerce_bool
from models import ProcessingConfig, TRIPLICATE_LAYOUT, build_plate_layout
from exclusions import exclusion_key_cols

ROWS = "ABCDEFGH"
TIMES = [0.0, 1.0, 2.0, 3.0, 4.0]  # one baseline read at t=0


def build_master_df():
    """One file, triplicate layout. Cols 1-3 = block 1 (CondA/HEK), cols 4-6 =
    block 2 (CondB/HEK). Row H = vehicle. Donor present. Engineered so:
      - H4 vehicle deviates strongly (vehicle warning on CondB)
      - A1 donor is very low (lum warning on CondA)
    """
    rng = np.random.default_rng(0)
    records = []
    cond_for_col = {}
    for col in range(1, 13):
        block = (col - 1) // 3  # 0..3
        cond_for_col[col] = ["CondA", "CondB", "CondC", "CondD"][block]
    for col in list(range(1, 13)):
        block_cond = cond_for_col[col]
        for r in ROWS:
            well = f"{r}{col}"
            is_veh = (r == "H")
            for t in TIMES:
                # Raw BRET ~ 0.5 baseline, rising after t=0
                raw = 0.5 + (0.0 if t == 0 else 0.05 * t)
                donor = 800.0 + rng.normal(0, 5)
                if well == "A1":
                    donor = 5.0  # below default lum threshold 100 -> lum warning
                veh_norm = 1.0 + (0.0 if t == 0 else 0.01 * t)
                if well == "H4":
                    veh_norm = 1.0 + (0.0 if t == 0 else 0.9)  # big deviation -> veh warning
                records.append({
                    "NCollector_version": "N Collector v2.0.4",
                    "Path": "undocumented path",
                    "Info_Sheet": "",
                    "File_Name": "exp1.xlsx",
                    "Date": "2024-02-01",
                    "Main_Plasmids": "P1",
                    "Applied_Exclusions": "",
                    "Is_Excluded": False,
                    "Is_Vehicle": is_veh,
                    "Transfection": block_cond,
                    "Cell_Line": "HEK",
                    "Ligand": "Iso",
                    "Ligand_Conc": float("nan") if is_veh else 1.0,
                    "Plate_Row": r,
                    "Replicate": "1",
                    "Well_ID": well,
                    "Time_(min)": t,
                    "PR_Time(min)": t,
                    "Donor_Raw_kinetic": donor,
                    "Acceptor_Raw_kinetic": donor * 0.8,
                    "Raw_BRET_kinetic": raw,
                    "Lab_BRET_kinetic": raw,
                    "Bl_Corrected_BRET": np.nan,
                    "Veh_Norm_Kinetic": veh_norm,
                    "Kinetic_Mean": np.nan,
                    "Raw_BRET_CRC": np.nan, "Lab_LP": np.nan, "Bl_LP": np.nan,
                    "Veh_Norm_LP": np.nan, "LP_Mean": np.nan,
                    "Lab_AUC": np.nan, "Bl_AUC": np.nan,
                    "Veh_Norm_AUC": np.nan, "AUC_Mean": np.nan,
                })
    return pd.DataFrame.from_records(records)


def make_fake_self(df):
    """A SimpleNamespace carrying just the attributes the pure methods touch,
    with GUI calls stubbed."""
    captured = {"index_records": None, "warnings_shown": None, "logs": []}

    s = SimpleNamespace()
    s.master_df = df
    s.experiment = []
    s.current_config = None
    s.ignored_warnings = set()
    s.var_lum_threshold = SimpleNamespace(get=lambda: 100)
    s.var_vehicle_threshold = SimpleNamespace(get=lambda: 0.2)
    # cb_lig['state'] -> 'normal' (ligand filter active, not single-locked)
    s.cb_lig = {"state": "normal"}
    s.btn_rerun_lum = None
    s.btn_rerun_vehicle = None
    s.log = lambda m: captured["logs"].append(str(m))


    # Stub GUI-touching methods used by the methods under test
    s.refresh_filter_options = lambda: None
    s.update_summary_table = lambda: None
    s.refresh_plot_helper_options = lambda: None
    s._update_quality_buttons_state = lambda: None
    s.clear_exclusion_list = lambda: None

    def _show(warns):
        captured["warnings_shown"] = warns
    s.show_warning_review = _show

    s._captured = captured
    return s


def bind(method_name, s):
    """Bind a real NCollectorApp method to the fake self."""
    return getattr(NCollectorApp, method_name).__get__(s, NCollectorApp)


def test_reconstruct():
    df = build_master_df()
    inp = reconstruct_file_inputs(df, "exp1.xlsx")
    assert inp is not None
    assert bool(inp["has_donor"]) is True
    assert inp["donor_wide"] is not None and inp["donor_wide"].shape[0] == len(TIMES)
    assert inp["veh_norm_wide"] is not None
    assert inp["date_str"] == "01.02.24", inp["date_str"]
    assert inp["excluded_wells"] == []
    # metadata maps column -> condition
    assert inp["column_metadata"][1].condition_name == "CondA"
    assert inp["column_metadata"][4].condition_name == "CondB"
    print("[TEST] reconstruct_file_inputs OK")


def test_collect_quality_warnings():
    df = build_master_df()
    s = make_fake_self(df)
    collect = bind("_collect_quality_warnings", s)
    files = ["exp1.xlsx"]
    warns = collect(files, run_lum=True, run_vehicle=True,
                    lum_threshold=100, vehicle_threshold=0.2, log_lum_skips=True)
    displays = [w["Display"] for w in warns]
    assert any("[LOW LUM]" in d for d in displays), f"expected Lum warning, got {displays}"
    assert any("[VEHICLE WARN]" in d for d in displays), f"expected Veh warning, got {displays}"
    # A1 low-lum should be present
    assert any("[LOW LUM]" in d and "A1" in d for d in displays)
    # H4 vehicle should be present
    assert any("[VEHICLE WARN]" in d and "H4" in d for d in displays)
    print(f"[TEST] _collect_quality_warnings OK (detected {len(warns)} warnings)")


def test_resolve_rule_to_targets():
    df = build_master_df()
    s = make_fake_self(df)
    resolve = bind("_resolve_rule_to_targets", s)
    norm_dates = pd.to_datetime(df["Date"], errors="coerce").dt.strftime("%d.%m.%y")

    # Rule: exclude Condition CondA (maps to Transfection) -> cols 1-3, all rows
    rule = {"Ligand": "All", "Date": "All", "Cell_Line": "All",
            "Condition": "CondA", "Replicate": "All", "Row": "All"}
    targets = resolve(rule, norm_dates)
    # Target tuples are keyed on EXCLUSION_KEY_COLS:
    # (File_Name, Well_ID, Ligand, Transfection, Cell_Line, Main_Plasmids)
    cols_hit = {int(t[1][1:]) for t in targets}
    assert cols_hit == {1, 2, 3}, cols_hit
    assert all(t[0] == "exp1.xlsx" for t in targets)
    # Identity columns are carried in the key
    assert all(t[3] == "CondA" for t in targets)  # Transfection

    # Whole-date rule -> every well of that file
    rule_date = {"Ligand": "All", "Date": "01.02.24", "Cell_Line": "All",
                 "Condition": "All", "Replicate": "All", "Row": "All"}
    targets_date = resolve(rule_date, norm_dates)
    assert len(targets_date) == 12 * 8, len(targets_date)  # 12 cols x 8 rows
    print("[TEST] _resolve_rule_to_targets OK (Condition->Transfection + whole-date)")


def test_build_index_from_master():
    df = build_master_df()
    # Exclude all of column 1 -> column should drop from the index (N count drops)
    df.loc[df["Well_ID"].str.match(r"^[A-H]1$"), "Is_Excluded"] = True
    s = make_fake_self(df)
    # bind real _finalize_index too (its groupby is part of what we test)
    s._finalize_index = bind("_finalize_index", s)
    build = bind("_build_index_from_master", s)
    idx = build()
    assert not idx.empty
    assert "Condition" in idx.columns
    assert set(idx["Condition"].unique()) <= {"CondA", "CondB", "CondC", "CondD"}
    cols_present = {int(str(w)) for w in idx["_ColIdx"]} if "_ColIdx" in idx.columns else None
    # Column 1 fully excluded -> not represented; cols 2-6 remain
    present_cols = set()
    for fname in idx["File_Name"]:
        pass
    # index rows are per (file,col); reconstruct from grouping count
    assert len(idx) == 11, f"expected 11 surviving columns (2-12), got {len(idx)}"
    print("[TEST] _build_index_from_master OK (fully-excluded column dropped)")


def test_full_apply_exclusions_then_recompute():
    """End-to-end-ish: flag CondA via the real apply_exclusions, ensure master is
    recomputed (Raw_BRET_kinetic NaN'd for excluded wells) and a warning review is
    triggered, all headless."""
    df = build_master_df()
    s = make_fake_self(df)

    # apply_exclusions reads pending_exclusions + lb_exclusions + rule_history_text
    s.pending_exclusions = [{"Ligand": "All", "Date": "All", "Cell_Line": "All",
                             "Condition": "CondA", "Replicate": "All", "Row": "All"}]
    s.rule_history_text = ""
    s.lbl_rules_summary = SimpleNamespace(config=lambda **k: None)
    s.lb_exclusions = SimpleNamespace(get=lambda a, b=None: ["Exclude CondA"])

    # bind the real helpers apply_exclusions depends on
    s._resolve_rule_to_targets = bind("_resolve_rule_to_targets", s)
    s._build_recompute_config = bind("_build_recompute_config", s)
    s._get_thresholds = bind("_get_thresholds", s)
    s._collect_quality_warnings = bind("_collect_quality_warnings", s)
    s._finalize_index = bind("_finalize_index", s)
    s._build_index_from_master = bind("_build_index_from_master", s)
    s.built_master_index = bind("built_master_index", s)
    # apply_exclusions now snapshots the excluded set and re-checks vehicles / refreshes
    # views afterwards. Bind the data helper; stub the GUI/warning side-effects the
    # assertions below don't cover.
    s._excluded_well_set = bind("_excluded_well_set", s)
    s._vehicle_warnings_for_affected = lambda affected_keys: []
    s.refresh_active_exclusions = lambda: None
    s.crc_window = None
    s._refresh_crc_window = bind("_refresh_crc_window", s)

    apply_exc = bind("apply_exclusions", s)
    apply_exc()

    # CondA wells (cols 1-3) should now be Is_Excluded and Raw_BRET NaN
    excl = s.master_df[s.master_df["Well_ID"].str.match(r"^[A-H][123]$")]
    assert excl["Is_Excluded"].all(), "CondA wells not flagged excluded"
    assert excl["Raw_BRET_kinetic"].isna().all(), "excluded Raw_BRET not NaN'd by recompute"
    # CondB (cols 4-6) untouched
    keep = s.master_df[s.master_df["Well_ID"].str.match(r"^[A-H][456]$")]
    assert not keep["Is_Excluded"].any()
    # Applied_Exclusions provenance written
    assert (s.master_df["Applied_Exclusions"].astype(str) != "").all()
    print("[TEST] apply_exclusions (CSV-mode, master-native) OK")


def test_duplicate_filename_disambiguation():
    """
    The concern: two loaded files share a File_Name but are different experiments
    (different Ligand / Main_Plasmids). With a (File_Name, Well_ID)-only key, a rule
    targeting one of them would wrongly flag the other's identically-named wells.
    The widened key (File_Name + Well_ID + Ligand + Transfection + Cell_Line +
    Main_Plasmids) must flag ONLY the intended rows.
    """
    rng = np.random.default_rng(1)
    recs = []
    # Two groups, SAME File_Name, distinguished by Ligand + Main_Plasmids.
    for ligand, plasmid in [("Iso", "P1"), ("Adr", "P2")]:
        for col in (1, 2, 3):
            for r in ROWS:
                for t in TIMES:
                    recs.append({
                        "File_Name": "dup.xlsx", "Date": "2024-02-01",
                        "Main_Plasmids": plasmid, "Is_Excluded": False,
                        "Is_Vehicle": (r == "H"), "Transfection": "CondX",
                        "Cell_Line": "HEK", "Ligand": ligand,
                        "Plate_Row": r, "Replicate": "1", "Well_ID": f"{r}{col}",
                        "Time_(min)": t, "Raw_BRET_kinetic": 0.5,
                        "Donor_Raw_kinetic": 800.0, "Veh_Norm_Kinetic": 1.0,
                        "Applied_Exclusions": "",
                    })
    df = pd.DataFrame.from_records(recs)
    s = make_fake_self(df)
    resolve = bind("_resolve_rule_to_targets", s)
    norm_dates = pd.to_datetime(df["Date"], errors="coerce").dt.strftime("%d.%m.%y")

    # Rule targets ONLY the "Iso" experiment.
    rule = {"Ligand": "Iso", "Date": "All", "Cell_Line": "All",
            "Condition": "All", "Replicate": "All", "Row": "All"}
    targets = resolve(rule, norm_dates)

    # Build the mask exactly as apply_exclusions does, via the real helper.
    key_cols = exclusion_key_cols(df)
    key_index = pd.MultiIndex.from_arrays([df[c].astype(str) for c in key_cols])
    mask = pd.Series(key_index.isin(list(targets)), index=df.index)

    flagged = df.loc[mask]
    assert (flagged["Ligand"] == "Iso").all(), "non-Iso rows were wrongly flagged"
    assert (flagged["Main_Plasmids"] == "P1").all()
    # Every Iso row is flagged; no Adr row is.
    assert mask.sum() == (df["Ligand"] == "Iso").sum()
    assert not mask[df["Ligand"] == "Adr"].any()
    print("[TEST] duplicate File_Name disambiguation OK "
          f"(flagged {mask.sum()} Iso rows, 0 Adr rows)")


if __name__ == "__main__":
    test_reconstruct()
    test_collect_quality_warnings()
    test_resolve_rule_to_targets()
    test_build_index_from_master()
    test_full_apply_exclusions_then_recompute()
    test_duplicate_filename_disambiguation()
    print("\nALL CSV-MODE TESTS PASSED")