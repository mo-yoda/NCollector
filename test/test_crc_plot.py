"""
Tests for crc_plot.py — the GUI-agnostic data extraction + 4PL fit behind the
concentration-response plot. No Tkinter / matplotlib needed here.
"""
import numpy as np
import pandas as pd

import crc_plot

ROWS_CONC = {"A": -9.0, "B": -8.0, "C": -7.0, "D": -6.0, "E": -5.0, "F": -4.0, "G": -3.0}


def _resp(x, scale=1.0):
    # Clean 4PL: bottom 1, top 5, logEC50 -6, hill 1.
    return scale * (1.0 + 4.0 / (1.0 + 10.0 ** ((-6.0 - x) * 1.0)))


def _crc_master(conditions=(("ATP", "HEK", "WT", "pA"),),
                dates=(("F1", "2026-06-05", 1.0), ("F2", "2026-06-20", 1.05))):
    """Small master: for each condition, each date contributes one dose row per concentration
    (A-G) plus a vehicle row H (no concentration)."""
    rows = []
    for (ligand, cell, transf, main) in conditions:
        for fbase, date, scale in dates:
            fname = f"{fbase}_{ligand}_{cell}_{transf}"
            for r, x in ROWS_CONC.items():
                rows.append({"File_Name": fname, "Date": date, "Main_Plasmids": main,
                             "Transfection": transf, "Cell_Line": cell, "Ligand": ligand,
                             "Ligand_Conc": x, "Plate_Row": r, "Well_ID": f"{r}1",
                             "Is_Vehicle": False, "AUC_Mean": _resp(x, scale)})
            rows.append({"File_Name": fname, "Date": date, "Main_Plasmids": main,
                         "Transfection": transf, "Cell_Line": cell, "Ligand": ligand,
                         "Ligand_Conc": float("nan"), "Plate_Row": "H", "Well_ID": "H1",
                         "Is_Vehicle": True, "AUC_Mean": 1.0})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# list_crc_conditions
# --------------------------------------------------------------------------- #
def test_list_conditions_orders_and_dedups():
    df = _crc_master(conditions=(("ATP", "HEK", "WT", "pA"),
                                 ("DA", "COS", "MUT", "pB")))
    conds = crc_plot.list_crc_conditions(df)
    assert len(conds) == 2
    assert conds[0] == {"ligand": "ATP", "cell_line": "HEK",
                        "transfection": "WT", "main_plasmids": "pA"}
    # Sorted by (ligand, cell, transfection): ATP before DA.
    assert [c["ligand"] for c in conds] == ["ATP", "DA"]


def test_list_conditions_skips_empty_and_unknown():
    df = _crc_master(conditions=(("ATP", "HEK", "WT", "pA"),
                                 ("ATP", "Unknown_1", "WT", "pA"),
                                 ("ATP", "HEK", "Empty", "pA")))
    conds = crc_plot.list_crc_conditions(df)
    assert conds == [{"ligand": "ATP", "cell_line": "HEK",
                      "transfection": "WT", "main_plasmids": "pA"}]


def test_list_conditions_empty_master():
    assert crc_plot.list_crc_conditions(pd.DataFrame()) == []
    assert crc_plot.list_crc_conditions(None) == []


# --------------------------------------------------------------------------- #
# extract_condition_data
# --------------------------------------------------------------------------- #
def test_extract_basic_two_replicates_vehicle_dropped():
    df = _crc_master()
    data = crc_plot.extract_condition_data(df, "ATP", "HEK", "WT")
    assert data is not None
    assert data.n == 2                                   # two dates = two replicates
    assert [r.date_label for r in data.replicates] == ["05.06.26", "20.06.26"]  # date-ordered
    # 7 dose concentrations, vehicle (NaN conc) excluded.
    assert data.concs == sorted(ROWS_CONC.values())
    assert len(data.replicates[0].concs) == 7
    assert all(np.isfinite(v) for v in data.replicates[0].values)
    # Across-date mean equals the mean of the two dates at the lowest concentration.
    lo = min(ROWS_CONC.values())
    expected = (_resp(lo, 1.0) + _resp(lo, 1.05)) / 2
    assert abs(data.mean_values[0] - expected) < 1e-9


def test_extract_reflects_exclusion_via_nan():
    df = _crc_master()
    # Emulate an exclusion: AUC_Mean NaN at conc -6 (row D) for BOTH dates -> that point drops.
    df.loc[df["Plate_Row"] == "D", "AUC_Mean"] = float("nan")
    data = crc_plot.extract_condition_data(df, "ATP", "HEK", "WT")
    assert -6.0 not in data.concs                        # dropped from the mean
    for rep in data.replicates:
        assert -6.0 not in rep.concs                     # and from each replicate
        assert len(rep.concs) == 6


def test_extract_partial_exclusion_one_date_only():
    df = _crc_master()
    # Exclude conc -6 for only ONE date -> mean at -6 still exists (from the other date).
    one = df["Plate_Row"].eq("D") & df["File_Name"].str.startswith("F1")
    df.loc[one, "AUC_Mean"] = float("nan")
    data = crc_plot.extract_condition_data(df, "ATP", "HEK", "WT")
    assert -6.0 in data.concs
    by_label = {r.date_label: r for r in data.replicates}
    assert -6.0 not in by_label["05.06.26"].concs        # F1 lost it
    assert -6.0 in by_label["20.06.26"].concs            # F2 kept it


def test_extract_ragged_concentrations_kept_at_true_x():
    # Two dates with DIFFERENT concentration sets: concentrations must never be averaged —
    # each date's points stay at their true x, and the mean curve is the union of x positions
    # (a concentration unique to one date appears at its true x with that date's value).
    rows = []
    base = dict(Main_Plasmids="pA", Transfection="WT", Cell_Line="HEK", Ligand="ATP",
                Is_Vehicle=False)
    # Date 1 (F1): concentrations -9, -8, -7
    for x in (-9.0, -8.0, -7.0):
        rows.append({**base, "File_Name": "F1", "Date": "2026-06-05",
                     "Ligand_Conc": x, "Plate_Row": "?", "Well_ID": "A1",
                     "AUC_Mean": _resp(x, 1.0)})
    # Date 2 (F2): concentrations -9, -8, -6  (has -6 instead of -7)
    for x in (-9.0, -8.0, -6.0):
        rows.append({**base, "File_Name": "F2", "Date": "2026-06-20",
                     "Ligand_Conc": x, "Plate_Row": "?", "Well_ID": "A1",
                     "AUC_Mean": _resp(x, 1.2)})
    df = pd.DataFrame(rows)

    data = crc_plot.extract_condition_data(df, "ATP", "HEK", "WT")
    by_label = {r.date_label: r for r in data.replicates}
    # Each date keeps its OWN concentrations at their true x.
    assert by_label["05.06.26"].concs == [-9.0, -8.0, -7.0]
    assert by_label["20.06.26"].concs == [-9.0, -8.0, -6.0]
    # Mean x = union of all concentrations (sorted), nothing averaged away.
    assert data.concs == [-9.0, -8.0, -7.0, -6.0]
    mean_by_x = dict(zip(data.concs, data.mean_values))
    # Shared concentrations: mean of both dates. Unique ones: that single date's value.
    assert abs(mean_by_x[-9.0] - (_resp(-9.0, 1.0) + _resp(-9.0, 1.2)) / 2) < 1e-9
    assert abs(mean_by_x[-7.0] - _resp(-7.0, 1.0)) < 1e-9   # only F1
    assert abs(mean_by_x[-6.0] - _resp(-6.0, 1.2)) < 1e-9   # only F2


def test_extract_none_for_missing_condition():
    df = _crc_master()
    assert crc_plot.extract_condition_data(df, "NOPE", "HEK", "WT") is None


def test_extract_averages_files_within_one_date():
    # Two files on the SAME date = one replicate; their per-conc values are averaged.
    df = _crc_master(dates=(("Fa", "2026-06-05", 1.0), ("Fb", "2026-06-05", 1.2)))
    data = crc_plot.extract_condition_data(df, "ATP", "HEK", "WT")
    assert data.n == 1                                   # one date -> one replicate
    lo = min(ROWS_CONC.values())
    expected = (_resp(lo, 1.0) + _resp(lo, 1.2)) / 2     # averaged across the two files
    assert abs(data.replicates[0].values[0] - expected) < 1e-9


# --------------------------------------------------------------------------- #
# fit_four_pl
# --------------------------------------------------------------------------- #
def test_fit_four_pl_recovers_sigmoid():
    x = sorted(ROWS_CONC.values())
    y = [_resp(v) for v in x]
    fit = crc_plot.fit_four_pl(x, y)
    assert fit is not None
    popt, xs, ys = fit
    assert len(popt) == 4
    bottom, top, logec50, hill = popt
    assert abs(bottom - 1.0) < 0.1
    assert abs(top - 5.0) < 0.1
    assert abs(logec50 - (-6.0)) < 0.1
    assert len(xs) == len(ys) == 200


def test_fit_four_pl_too_few_points():
    assert crc_plot.fit_four_pl([-9.0, -7.0, -5.0], [1.0, 3.0, 5.0]) is None


def test_fit_four_pl_handles_nan():
    x = [-9.0, -8.0, -7.0, -6.0, -5.0, float("nan")]
    y = [1.0, 1.5, 3.0, 4.0, 5.0, 9.0]
    fit = crc_plot.fit_four_pl(x, y)          # NaN dropped, 5 valid points remain
    assert fit is not None
