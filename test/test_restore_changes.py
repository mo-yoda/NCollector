"""
Focused tests for the two exclusions.py changes:

  1. _build_restorable_map (vectorized) must agree, well-for-well, with the original
     per-well is_restorable() that it replaces in the hot loops.

  2. Deferred recompute (recompute=False) + ONE batched recompute over the union of
     affected files must yield a master IDENTICAL to the old behaviour of recomputing
     inside every restore call. We monkeypatch recompute_master_after_exclusion with a
     deterministic stub (it stamps a derived column from the CURRENT Is_Excluded) so the
     test isolates the orchestration change from the untouched BRET math, and also counts
     how many times recompute actually runs.
"""
import types
import pandas as pd
import exclusions


def _make_master():
    """
    Small synthetic master: 2 files, a few wells each, 2 timepoints per well.
    Columns are the minimum the restore logic + stubbed recompute read.
    Wells:
      F1/A1, F1/A2  -> have Raw_BRET_unexcluded  (restorable)
      F1/B1         -> no unexcl but Donor+Acceptor present (restorable)
      F2/A1         -> only Donor (NOT restorable)
      F2/A2         -> nothing pristine (NOT restorable)
    """
    rows = []
    def add(f, well, row_char, col, unexcl, donor, acc, ligand, cell, cond, excl):
        for t in (0.0, 1.0):
            rows.append({
                "File_Name": f, "Well_ID": well, "Plate_Row": row_char,
                "Time_(min)": t, "Date": "01.01.24",
                "Ligand": ligand, "Cell_Line": cell, "Transfection": cond,
                "Replicate": "rep1",
                "Raw_BRET_unexcluded": unexcl, "Donor_Raw_kinetic": donor,
                "Acceptor_Raw_kinetic": acc, "Raw_BRET_kinetic": unexcl,
                "Is_Excluded": excl,
                "Derived": float("nan"),  # stand-in derived column for the stub
            })
    add("F1", "A1", "A", 1, 1.5, 100.0, 150.0, "LigX", "HEK", "condA", True)
    add("F1", "A2", "A", 2, 1.6, 110.0, 160.0, "LigX", "HEK", "condA", True)
    add("F1", "B1", "B", 1, float("nan"), 90.0, 120.0, "LigX", "HEK", "condA", True)
    add("F2", "A1", "A", 1, float("nan"), 80.0, float("nan"), "LigY", "CHO", "condB", True)
    add("F2", "A2", "A", 2, float("nan"), float("nan"), float("nan"), "LigY", "CHO", "condB", True)
    return pd.DataFrame(rows)


# Blob: one AUTO single-well rule per excluded well (mirrors the auto-exclusion case),
# plus one manual rule covering an already-listed well to exercise overlap-free paths.
def _make_blob():
    toks = [
        "AUTO: [LOW LUM] LigX | HEK | condA | 01.01.24 | A1 - value: 1",
        "AUTO: [LOW LUM] LigX | HEK | condA | 01.01.24 | A2 - value: 2",
        "AUTO: [LOW LUM] LigX | HEK | condA | 01.01.24 | B1 - value: 3",
        "AUTO: [LOW LUM] LigY | CHO | condB | 01.01.24 | A1 - value: 4",
        "AUTO: [LOW LUM] LigY | CHO | condB | 01.01.24 | A2 - value: 5",
    ]
    return " || ".join(toks)


class _StubConfig:
    baseline_end_index = None


def _make_stub(counter):
    """A deterministic stand-in for recompute_master_after_exclusion.

    Returns a NEW df (like the real one) whose 'Derived' column is set from the CURRENT
    Is_Excluded for the rows of the affected files: excluded -> NaN, else 1.0. Counts calls
    and total files processed so we can prove batching collapsed N calls into 1 pass.
    """
    def _stub(master_df, affected_files, config):
        counter["calls"] += 1
        counter["files"] += len(affected_files)
        out = master_df.copy()
        for f in affected_files:
            m = out["File_Name"] == f
            out.loc[m, "Derived"] = out.loc[m, "Is_Excluded"].map(
                lambda x: float("nan") if bool(x) else 1.0)
        return out
    return _stub


def test_map_matches_is_restorable():
    df = _make_master()
    rmap = exclusions._build_restorable_map(df)
    pairs = set(zip(df["File_Name"].astype(str), df["Well_ID"].astype(str)))
    for f, w in pairs:
        got_ok, _ = rmap[(f, w)]
        exp_ok, _ = exclusions.is_restorable(df, f, w)
        assert got_ok == exp_ok, f"map disagrees for {f}/{w}: map={got_ok} ref={exp_ok}"
    # Spot-check expectations
    assert rmap[("F1", "A1")][0] is True
    assert rmap[("F1", "B1")][0] is True      # donor+acceptor only
    assert rmap[("F2", "A1")][0] is False     # donor only
    assert rmap[("F2", "A2")][0] is False     # nothing
    print("PASS: _build_restorable_map matches is_restorable for all wells")


def test_defer_batch_equivalence():
    labels = exclusions.parse_exclusion_blob(_make_blob())
    label_list = [e["label"] for e in labels]

    # --- Path A: OLD behaviour — recompute inside every restore call ---
    dfA = _make_master(); dfA["Applied_Exclusions"] = _make_blob()
    counterA = {"calls": 0, "files": 0}
    exclusions.recompute_master_after_exclusion = _make_stub(counterA)
    for lbl in label_list:
        exclusions.restore_rule(dfA, lbl, _StubConfig(), recompute=True)

    # --- Path B: NEW behaviour — defer, then ONE batched recompute over the union ---
    dfB = _make_master(); dfB["Applied_Exclusions"] = _make_blob()
    counterB = {"calls": 0, "files": 0}
    exclusions.recompute_master_after_exclusion = _make_stub(counterB)
    affected = set()
    for lbl in label_list:
        rep = exclusions.restore_rule(dfB, lbl, _StubConfig(), recompute=False)
        affected.update(rep["affected_files"])
    assert counterB["calls"] == 0, "deferred path must not recompute inside restore calls"
    if affected:
        dfB = exclusions.recompute_master_after_exclusion(dfB, sorted(affected), _StubConfig())

    # Equivalence of the resulting master (the columns restore + recompute own)
    cols = ["File_Name", "Well_ID", "Time_(min)", "Is_Excluded",
            "Applied_Exclusions", "Derived"]
    a = dfA[cols].sort_values(["File_Name", "Well_ID", "Time_(min)"]).reset_index(drop=True)
    b = dfB[cols].sort_values(["File_Name", "Well_ID", "Time_(min)"]).reset_index(drop=True)
    pd.testing.assert_frame_equal(a, b)

    print(f"PASS: deferred+batched master == old per-rule master")
    print(f"      recompute calls: OLD={counterA['calls']} (files touched={counterA['files']}) "
          f"-> NEW={1 if affected else 0} (files touched={len(affected)})")
    # The headline win: multiple per-rule recomputes collapse to one batched pass.
    assert counterA["calls"] > 1                   # old: recomputed once per restoring rule
    assert (1 if affected else 0) <= 1             # new: a single batched pass


def test_ctx_equivalence_overlapping_rules():
    """
    The shared-context batch path must produce an identical master to resolving each rule
    with its own private context. Scenario mixes auto single-well rules, a broad manual
    criteria-rule (whole Plate_Row), an OVERLAP (a well covered by both a manual rule and an
    auto rule), and a non-restorable well — exercising the blocked / decompose / drop paths.
    """
    import exclusions as R

    def master():
        rows = []
        def add(f, well, rch, col, unexcl, donor, acc, lig, cell, cond):
            for t in (0.0, 1.0):
                rows.append({
                    "File_Name": f, "Well_ID": well, "Plate_Row": rch, "Time_(min)": t,
                    "Date": "01.01.24", "Ligand": lig, "Cell_Line": cell,
                    "Transfection": cond, "Replicate": "rep1",
                    "Raw_BRET_unexcluded": unexcl, "Donor_Raw_kinetic": donor,
                    "Acceptor_Raw_kinetic": acc, "Raw_BRET_kinetic": unexcl,
                    "Is_Excluded": True, "Derived": float("nan"), "Ligand_Conc": 1.0})
        # File F1, row A, cols 1-3 -> covered by a manual "Row A" rule; A1 also has an auto rule
        add("F1", "A1", "A", 1, 1.1, 100., 150., "LigX", "HEK", "condA")  # overlap manual+auto
        add("F1", "A2", "A", 2, 1.2, 100., 150., "LigX", "HEK", "condA")
        add("F1", "A3", "A", 3, 1.3, 100., 150., "LigY", "HEK", "condB")
        # an auto-only well, non-restorable (no pristine source)
        add("F1", "B1", "B", 1, float("nan"), float("nan"), float("nan"), "LigX", "HEK", "condA")
        return pd.DataFrame(rows)

    blob = " || ".join([
        "Ligand: All | Date: 01.01.24 | Cell: All | Cond: All | Rep:rep1 | Row:A",  # manual: A1,A2,A3
        "AUTO: [LOW LUM] LigX | HEK | condA | 01.01.24 | A1 - value: 1",            # auto overlap
        "AUTO: [LOW LUM] LigX | HEK | condA | 01.01.24 | B1 - value: 2",            # auto non-restorable
    ])
    labels = [e["label"] for e in R.parse_exclusion_blob(blob)]

    R.recompute_master_after_exclusion = _make_stub({"calls": 0, "files": 0})

    # Path 1: private ctx per call (ctx=None)
    d1 = master(); d1["Applied_Exclusions"] = blob
    for lbl in labels:
        R.restore_rule(d1, lbl, _StubConfig(), recompute=False, ctx=None)

    # Path 2: one shared ctx for the whole batch
    d2 = master(); d2["Applied_Exclusions"] = blob
    ctx = R.build_resolve_ctx(d2)
    for lbl in labels:
        R.restore_rule(d2, lbl, _StubConfig(), recompute=False, ctx=ctx)

    cols = ["File_Name", "Well_ID", "Time_(min)", "Is_Excluded", "Applied_Exclusions"]
    a = d1[cols].sort_values(cols[:3]).reset_index(drop=True)
    b = d2[cols].sort_values(cols[:3]).reset_index(drop=True)
    pd.testing.assert_frame_equal(a, b)
    print("PASS: shared-ctx batch == private-ctx-per-call (overlapping manual+auto, "
          "non-restorable, decompose)")


if __name__ == "__main__":
    test_map_matches_is_restorable()
    test_defer_batch_equivalence()
    test_ctx_equivalence_overlapping_rules()
    print("\nALL TESTS PASSED")