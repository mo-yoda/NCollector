#!/usr/bin/env python3
"""
restore_harness.py — standalone equivalence / round-trip harness for the reversible
exclusion backend (Prompt 4). Extends the Prompt-1 equivalence harness.

Run headless:  python restore_harness.py

Checks
------
1. CHANNEL EQUIVALENCE: reconstructed-from-channels ratio (Acceptor/Donor) ≈ the
   imported instrument Raw_BRET_kinetic within rtol=1e-4, atol=1e-6 (the same gate the
   merge uses; the user verified max abs diff ~5e-6 in R).
2. ROUND-TRIP INVERSE: snapshot master_0 (no exclusions); apply an exclusion set S that
   includes a vehicle (row H) well AND a labeling-control well -> master_1; restore all
   of S -> master_2; assert master_2 ≈ master_0 (numeric isclose; exact for string/meta).
3. PARTIAL / OVERLAP: exclude via two overlapping rules; restore one; assert wells still
   matched by the OTHER rule remain excluded, and wells matched only by the restored rule
   come back.
4. LEGACY PATH: import a master with a " | " blob and NCollector_version < 2.0.5; confirm
   ensure_master_csv_schema rewrites the blob to " || " and back-fills Raw_BRET_unexcluded
   from channels; confirm the active-rules parser recovers the same wells the original
   blob implied.
"""

import sys
import numpy as np
import pandas as pd

from models import ProcessingConfig, TRIPLICATE_LAYOUT
from processing import recompute_master_after_exclusion
import exclusions
from exclusions import (parse_exclusion_blob, list_active_exclusions,
                     restore_rule, is_restorable, _resolve_rule_wells,
                     _resolve_criteria, _excluded_mask)
from export import ensure_master_csv_schema

RTOL, ATOL = 1e-4, 1e-6

# Columns that should match exactly (not via isclose).
_STR_COLS = ["NCollector_version", "File_Name", "Date", "Cell_Line", "Transfection",
             "Ligand", "Plate_Row", "Replicate", "Well_ID", "Is_Vehicle"]

_PASS, _FAIL = "  [PASS]", "  [FAIL]"
_failures = []


def check(cond, msg):
    print((_PASS if cond else _FAIL), msg)
    if not cond:
        _failures.append(msg)
    return cond


# --------------------------------------------------------------------------- #
# Synthetic master builders
# --------------------------------------------------------------------------- #

def make_file(fname, date, cell, cond_map, rep_map, ligand, times, seed):
    """One file's long-format rows. Raw_BRET_unexcluded = exact Acceptor/Donor;
    Raw_BRET_kinetic = that ratio + tiny instrument-like noise (the 'imported' value)."""
    rng = np.random.default_rng(seed)
    rows = []
    for col in sorted(cond_map.keys()):
        for ri, row_char in enumerate("ABCDEFGH"):
            well = f"{row_char}{col}"
            conc = float("nan") if row_char == "H" else float((ri + 1) * 10)
            for t in times:
                donor = 1000.0 + rng.uniform(0, 200)
                acceptor = donor * (0.3 + 0.05 * ri) + rng.uniform(0, 5)
                ratio = acceptor / donor
                stored = ratio * (1.0 + rng.normal(0, 1e-6))   # instrument rounding
                rows.append({
                    "NCollector_version": "N Collector v2.0.5",
                    "File_Name": fname, "Date": date, "Cell_Line": cell,
                    "Transfection": cond_map[col], "Ligand": ligand,
                    "Ligand_Conc": conc, "Plate_Row": row_char,
                    "Replicate": rep_map[col], "Well_ID": well, "Time_(min)": float(t),
                    "PR_Time(min)": float(t),
                    "Donor_Raw_kinetic": donor, "Acceptor_Raw_kinetic": acceptor,
                    "Raw_BRET_unexcluded": ratio, "Raw_BRET_kinetic": stored,
                    "Is_Excluded": False, "Is_Vehicle": (row_char == "H"),
                    "Applied_Exclusions": "None", "Main_Plasmids": "P1 + P2",
                    "Path": "synthetic", "Info_Sheet": "",
                })
    return pd.DataFrame(rows)


def build_master():
    times = [-1.0, 0.0, 1.0, 2.0, 3.0]
    # File A: triplicate, two conditions (cols 1-3 CondX, cols 4-6 CondY).
    condA = {1: "CondX", 2: "CondX", 3: "CondX", 4: "CondY", 5: "CondY", 6: "CondY"}
    repA = {1: "1", 2: "2", 3: "3", 4: "1", 5: "2", 6: "3"}
    fa = make_file("fileA.xlsx", "01.01.24", "CellA", condA, repA, "LigA", times, 1)
    # File B: a labeling block (cols 1-4): three replicates + one labeling control.
    condB = {1: "CondZ", 2: "CondZ", 3: "CondZ", 4: "CondZ"}
    repB = {1: "1", 2: "2", 3: "3", 4: "labeling control"}
    fb = make_file("fileB.xlsx", "02.01.24", "CellB", condB, repB, "LigA", times, 2)
    return pd.concat([fa, fb], ignore_index=True)


def cfg():
    return ProcessingConfig(labeling_correction=False, plate_layout=TRIPLICATE_LAYOUT,
                            baseline_end_index=None, lum_threshold=100,
                            vehicle_warning_threshold=0.2)


def apply_rule_tokens(master, tokens, config):
    """Set Is_Excluded for the union of the given rule tokens, write the blob, recompute."""
    m = master.copy()
    union = set()
    for entry in parse_exclusion_blob(" || ".join(tokens)):
        union |= _resolve_rule_wells(m, entry, only_excluded=False)
    key = pd.MultiIndex.from_arrays([m["File_Name"].astype(str), m["Well_ID"].astype(str)])
    m.loc[pd.Series(key.isin(list(union)), index=m.index), "Is_Excluded"] = True
    m["Applied_Exclusions"] = " || ".join(tokens)
    affected = sorted({f for f, w in union})
    m = recompute_master_after_exclusion(m, affected, config)
    return m, union


def frames_equiv(a, b, label):
    """master_2 ≈ master_0: align on (File_Name, Well_ID, Time), numeric isclose, exact str."""
    keys = ["File_Name", "Well_ID", "Time_(min)"]
    a2 = a.sort_values(keys).reset_index(drop=True)
    b2 = b.sort_values(keys).reset_index(drop=True)
    if not check(list(a2.columns) == list(b2.columns), f"{label}: same columns"):
        return
    if not check(len(a2) == len(b2), f"{label}: same row count"):
        return
    ok = True
    for col in a2.columns:
        if col in _STR_COLS or a2[col].dtype == object:
            same = (a2[col].astype(str).fillna("nan") == b2[col].astype(str).fillna("nan")).all()
        else:
            x = pd.to_numeric(a2[col], errors="coerce").to_numpy()
            y = pd.to_numeric(b2[col], errors="coerce").to_numpy()
            same = np.allclose(x, y, rtol=RTOL, atol=ATOL, equal_nan=True)
        if not same:
            ok = False
            print(f"        diff in column: {col}")
    check(ok, f"{label}: all columns match within tolerance")


# --------------------------------------------------------------------------- #
# Tests
# --------------------------------------------------------------------------- #

def test_channel_equivalence(raw_master):
    print("\n[1] Channel reconstruction ≈ imported instrument ratio")
    donor = pd.to_numeric(raw_master["Donor_Raw_kinetic"], errors="coerce")
    acceptor = pd.to_numeric(raw_master["Acceptor_Raw_kinetic"], errors="coerce")
    recon = acceptor / donor.where((donor != 0) & donor.notna())
    imported = pd.to_numeric(raw_master["Raw_BRET_kinetic"], errors="coerce")
    diff = (recon - imported).abs()
    check(np.allclose(recon, imported, rtol=RTOL, atol=ATOL, equal_nan=True),
          f"channel ratio within tolerance (max abs diff {diff.max():.2e})")


def test_round_trip(raw_master, config):
    print("\n[2] Round-trip exclude -> restore is the identity")
    master_0 = recompute_master_after_exclusion(raw_master.copy(),
                                                sorted(raw_master["File_Name"].unique()),
                                                config)
    # S includes a vehicle (row H) well AND a labeling-control well (fileB col4).
    veh_token = exclusions._well_token(master_0, "fileA.xlsx", "H1")
    lab_token = exclusions._well_token(master_0, "fileB.xlsx", "A4")
    master_1, S = apply_rule_tokens(master_0, [veh_token, lab_token], config)
    check(master_1.loc[(master_1["File_Name"] == "fileA.xlsx") &
                       (master_1["Well_ID"] == "H1"), "Raw_BRET_kinetic"].isna().all()
          and master_1.loc[(master_1["File_Name"] == "fileB.xlsx") &
                           (master_1["Well_ID"] == "A4"), "Raw_BRET_kinetic"].isna().all(),
          "S wells are NaN-d in master_1 (vehicle + labeling-control)")

    master_2 = master_1.copy()
    for entry in list_active_exclusions(master_2):
        restore_rule(master_2, entry["label"], config)
    check(not _excluded_mask(master_2).any(), "all of S restored (no Is_Excluded remain)")
    frames_equiv(master_0, master_2, "round-trip")


def test_partial_overlap(raw_master, config):
    print("\n[3] Partial restore with two overlapping rules")
    master_0 = recompute_master_after_exclusion(raw_master.copy(),
                                                sorted(raw_master["File_Name"].unique()),
                                                config)
    # Rule1: all of CondX (fileA cols 1-3, rows A-H). Rule2: all of row A (every file/col).
    rule1 = "Ligand: All | Date: All | Cell: CellA | Cond: CondX | Rep:All | Row:All"
    rule2 = "Ligand: All | Date: All | Cell: All | Cond: All | Rep:All | Row:A"
    master_1, _ = apply_rule_tokens(master_0, [rule1, rule2], config)

    overlap = _resolve_rule_wells(master_1,
                                  parse_exclusion_blob(rule1)[0], only_excluded=True) & \
              _resolve_rule_wells(master_1,
                                  parse_exclusion_blob(rule2)[0], only_excluded=True)
    check(len(overlap) > 0, f"rules overlap on {sorted(overlap)}")

    only_rule1 = _resolve_rule_wells(master_1,
                                     parse_exclusion_blob(rule1)[0], only_excluded=True) - overlap
    rep = restore_rule(master_1, rule1, config)
    print(f"        report: restored={rep['restored_count']} "
          f"blocked={rep['blocked_by_other_rule_count']} decomposed={rep['rule_decomposed']}")

    excl_now = exclusions._all_excluded_wells(master_1)
    check(overlap.issubset(excl_now),
          "overlap wells (matched by rule2) STAY excluded after restoring rule1")
    check(only_rule1.isdisjoint(excl_now),
          "wells matched ONLY by rule1 are restored")
    # rule2 must still re-resolve to its wells (incl. the overlap), invariant preserved.
    rule2_now = _resolve_rule_wells(master_1, parse_exclusion_blob(rule2)[0], only_excluded=True)
    check(overlap.issubset(rule2_now), "rule2 still re-resolves to the overlap wells")


def test_legacy_path(config):
    print("\n[4] Legacy ' | ' blob + missing Raw_BRET_unexcluded migration")
    times = [-1.0, 0.0, 1.0, 2.0, 3.0]
    condA = {1: "CondX", 2: "CondX", 3: "CondX", 4: "CondY", 5: "CondY", 6: "CondY"}
    repA = {1: "1", 2: "2", 3: "3", 4: "1", 5: "2", 6: "3"}
    legacy = make_file("legacy.xlsx", "01.01.24", "CellA", condA, repA, "LigA", times, 7)
    legacy["NCollector_version"] = "N Collector v2.0.4"
    legacy = legacy.drop(columns=["Raw_BRET_unexcluded"])     # pre-v2.0.5: column absent

    # Two manual rules, written in the LEGACY collapsed form (inter-rule sep == " | ").
    r1 = "Ligand: All | Date: 01.01.24 | Cell: CellA | Cond: CondX | Rep:All | Row:All"
    r2 = "Ligand: All | Date: 01.01.24 | Cell: CellA | Cond: CondY | Rep:All | Row:All"
    legacy_blob = r1 + " | " + r2
    legacy["Applied_Exclusions"] = legacy_blob
    # The rules imply CondX + CondY wells excluded; NaN their Raw_BRET_kinetic like a real
    # legacy master would (channels stay intact -> still reconstructable).
    implied = set()
    for entry in parse_exclusion_blob(legacy_blob):
        implied |= _resolve_rule_wells(legacy, entry, only_excluded=False)
    key = pd.MultiIndex.from_arrays([legacy["File_Name"].astype(str),
                                     legacy["Well_ID"].astype(str)])
    m = pd.Series(key.isin(list(implied)), index=legacy.index)
    legacy.loc[m, "Is_Excluded"] = True
    legacy.loc[m, "Raw_BRET_kinetic"] = float("nan")

    migrated, was_mod, _ = ensure_master_csv_schema(legacy.copy())
    check(was_mod, "ensure_master_csv_schema reported modifications")
    new_blob = migrated["Applied_Exclusions"].dropna().iloc[0]
    check(" || " in new_blob, "blob migrated to ' || ' separator")
    check(len(parse_exclusion_blob(new_blob)) == 2, "parser recovers 2 rules")

    # Raw_BRET_unexcluded back-filled from channels even for currently-excluded wells.
    check("Raw_BRET_unexcluded" in migrated.columns
          and migrated["Raw_BRET_unexcluded"].notna().all(),
          "Raw_BRET_unexcluded back-filled from channels (incl. excluded wells)")

    # Active-rules parser recovers the SAME wells the original blob implied.
    recovered = set()
    for entry in list_active_exclusions(migrated):
        recovered |= set(entry["wells"])
    check(recovered == implied,
          f"parser recovers original wells ({len(recovered)} == {len(implied)})")

    # And those legacy-excluded wells are now restorable (channels present).
    f, w = sorted(implied)[0]
    ok, reason = is_restorable(migrated, f, w)
    check(ok, f"legacy-excluded well {f}/{w} is restorable via channels")


def main():
    raw = build_master()
    config = cfg()
    test_channel_equivalence(raw)
    test_round_trip(raw, config)
    test_partial_overlap(raw, config)
    test_legacy_path(config)

    print("\n" + "=" * 60)
    if _failures:
        print(f"FAILED ({len(_failures)}):")
        for f in _failures:
            print("   -", f)
        sys.exit(1)
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
