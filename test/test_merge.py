"""
Test suite for merge.py (the GUI-agnostic master-merge engine).

Covers EVERY Issue.code path in the safety taxonomy plus the seven acceptance
criteria from the spec:

  Issue codes
    LABELING_MISMATCH            (forbidden, no override)
    TIME_VECTOR_DIVERGENCE       (forbidden, no override)
    FILENAME_DATA_COLLISION      (forbidden + rename / cancel resolutions)
    RAW_CHANNELS_REQUIRED        (forbidden, no override)
    FILENAME_DUPLICATE           (identical copies -> info dedupe)
    FILENAME_EXCLUSION_CONFLICT  (same file, sources disagree on excluded wells ->
                                  needs_input + dup_keep resolution)
    CONDITION_PARTITION_OVERLAP  (deferred; hook is a no-op, emits nothing)
    SCHEMA_MIGRATED              (info)
    NEW_CONDITIONS               (info)
    EXTRA_COLUMNS_DROPPED        (info)
    FILENAME_COLLISION_RENAMED   (info, on a resolved rename)

  Acceptance criteria
    1. self-merge identity (dedupe to original; one FILENAME_DUPLICATE->info per file)
    2. labeling master + non-labeling master -> MergeForbidden (no override)
    3. disjoint files + conditions -> len(A)+len(B) rows; NEW_CONDITIONS == B's count
    4. collision rename keeps both files distinct AND recompute runs cleanly on the
       merged frame for a sanity excluded-well flip; cancel/unset raises
    5. exclusion-rule scoping: a manual rule from source A never matches source B's
       rows; a manual row-rule stays ONE rule (verified via list_active_exclusions +
       a targeted restore_rule); per-well decomposition only for the pathological
       shared-everything overlap
    6. one NaN raw-channel value -> MergeForbidden (no override); folder source passes
    7. (this file) a test per Issue.code
"""

import numpy as np
import pandas as pd
import pytest

import models
from models import MASTER_COLUMNS, ProcessingConfig
import merge
import restore
from merge import MergeSource, MergeForbidden

ROWS = "ABCDEFGH"


# --------------------------------------------------------------------------- #
# Fixture builder
# --------------------------------------------------------------------------- #
def make_master(file_name="F1", date="2026-05-21", cell="HEK", ligand="ATP",
                main="pA", conds=None, cols=None, times=(0.0, 1.0, 2.0, 3.0),
                excluded_wells=(), labeling=False, blob="None", seed=0):
    """Build a small valid master: columns 1..N, rows A-H, given timepoints.

    Mirrors a freshly compiled master: Donor/Acceptor channels on every row,
    Raw_BRET_unexcluded = Acceptor/Donor, Time 0.0 present (baseline), row H
    flagged Is_Vehicle. Deterministic per `seed`.
    """
    if cols is None:
        cols = [1, 2, 3]                       # one triplicate block
    if conds is None:
        conds = {c: ("WT" if not labeling
                     else ("labeling control" if c == cols[-1] else "WT"))
                 for c in cols}
    rng = np.random.default_rng(seed)
    rows_out = []
    for c in cols:
        cond = conds[c]
        replicate = "labeling control" if (labeling and cond == "labeling control") else "1"
        for rc in ROWS:
            well = f"{rc}{c}"
            is_veh = (rc == "H")
            base = rng.uniform(0.5, 1.5)
            for t in times:
                donor = 1000.0 + 10 * t + rng.uniform(0, 1)
                acceptor = donor * (base + 0.01 * t)
                raw = acceptor / donor
                rows_out.append({
                    "NCollector_version": models.APP_VERSION,
                    "Path": f"/data/{file_name}",
                    "Info_Sheet": "",
                    "File_Name": file_name,
                    "Date": date,
                    "Main_Plasmids": main,
                    "Applied_Exclusions": blob,
                    "Is_Excluded": well in excluded_wells,
                    "Is_Vehicle": is_veh,
                    "Transfection": cond,
                    "Cell_Line": cell,
                    "Ligand": ligand,
                    "Ligand_Conc": float('nan') if is_veh else float(-9 + ROWS.index(rc)),
                    "Plate_Row": rc,
                    "Replicate": replicate,
                    "Well_ID": well,
                    "Time_(min)": t,
                    "PR_Time(min)": t,
                    "Donor_Raw_kinetic": donor,
                    "Acceptor_Raw_kinetic": acceptor,
                    "Raw_BRET_unexcluded": raw,
                    "Raw_BRET_kinetic": float('nan') if well in excluded_wells else raw,
                })
    df = pd.DataFrame(rows_out)
    for col in MASTER_COLUMNS:
        if col not in df.columns:
            df[col] = float('nan')
    return df[MASTER_COLUMNS].copy()


def codes(issues):
    return [i["code"] for i in issues]


def find(issues, code):
    return [i for i in issues if i["code"] == code]


# --------------------------------------------------------------------------- #
# LABELING_MISMATCH
# --------------------------------------------------------------------------- #
def test_labeling_mismatch_forbidden():
    A = make_master("F1", labeling=False)
    B = make_master("G1", labeling=True, cols=[1, 2, 3, 4])
    rep = merge.classify_sources([MergeSource("A", A), MergeSource("B", B)])
    assert merge.LABELING_MISMATCH in codes(rep.forbidden)
    with pytest.raises(MergeForbidden) as ei:
        merge.merge_masters([MergeSource("A", A), MergeSource("B", B)])
    assert merge.LABELING_MISMATCH in codes(ei.value.issues)


def test_labeling_mismatch_has_no_override():
    A = make_master("F1", labeling=False)
    B = make_master("G1", labeling=True, cols=[1, 2, 3, 4])
    # There is no override: even passing the old key must still raise.
    with pytest.raises(MergeForbidden) as ei:
        merge.merge_masters([MergeSource("A", A), MergeSource("B", B)],
                            {"allow_labeling_mismatch": True})
    assert merge.LABELING_MISMATCH in codes(ei.value.issues)


# --------------------------------------------------------------------------- #
# TIME_VECTOR_DIVERGENCE
# --------------------------------------------------------------------------- #
def test_time_vector_divergence_forbidden_no_override():
    A = make_master("F1", times=(0.0, 1.0, 2.0, 3.0))
    B = make_master("F2", times=(0.0, 1.0, 2.0, 3.0, 4.0))
    rep = merge.classify_sources([MergeSource("A", A), MergeSource("B", B)])
    assert merge.TIME_VECTOR_DIVERGENCE in codes(rep.forbidden)
    iss = find(rep.forbidden, merge.TIME_VECTOR_DIVERGENCE)[0]
    # context records the differing vectors (lengths at least).
    assert "context" in iss
    # No override exists -> still raises even if someone tries to pass one.
    with pytest.raises(MergeForbidden):
        merge.merge_masters([MergeSource("A", A), MergeSource("B", B)],
                            {"allow_time_divergence": True})


# --------------------------------------------------------------------------- #
# FILENAME_DATA_COLLISION
# --------------------------------------------------------------------------- #
def test_filename_data_collision_forbidden():
    A = make_master("F1", seed=1)
    B = make_master("F1", seed=999)          # same name, different raw data
    rep = merge.classify_sources([MergeSource("A", A), MergeSource("B", B)])
    assert merge.FILENAME_DATA_COLLISION in codes(rep.forbidden)


def test_filename_data_collision_cancel_or_unset_raises():
    A = make_master("F1", seed=1)
    B = make_master("F1", seed=999)
    with pytest.raises(MergeForbidden):
        merge.merge_masters([MergeSource("A", A), MergeSource("B", B)])
    with pytest.raises(MergeForbidden):
        merge.merge_masters([MergeSource("A", A), MergeSource("B", B)],
                            {"filename_collision": "cancel"})


def test_filename_data_collision_rename_keeps_both_distinct():
    A = make_master("F1", seed=1)
    B = make_master("F1", seed=999)
    m, rep = merge.merge_masters([MergeSource("A", A), MergeSource("B", B)],
                                 {"filename_collision": "rename"})
    assert sorted(m["File_Name"].unique()) == ["F1", "F1#B"]
    assert len(m) == len(A) + len(B)
    assert merge.FILENAME_COLLISION_RENAMED in codes(rep.auto)


# --------------------------------------------------------------------------- #
# RAW_CHANNELS_REQUIRED
# --------------------------------------------------------------------------- #
def test_raw_channels_required_forbidden_no_override():
    A = make_master("F1")
    bad = make_master("F2")
    bad.loc[bad.index[0], "Donor_Raw_kinetic"] = float('nan')
    rep = merge.classify_sources([MergeSource("A", A), MergeSource("B", bad)])
    iss = find(rep.forbidden, merge.RAW_CHANNELS_REQUIRED)
    assert iss, "expected RAW_CHANNELS_REQUIRED"
    ctx = iss[0]["context"]["files"]
    assert "F2" in ctx and ctx["F2"]["nan_rows"] >= 1 and ctx["F2"]["sample_wells"]
    # No override.
    with pytest.raises(MergeForbidden):
        merge.merge_masters([MergeSource("A", A), MergeSource("B", bad)],
                            {"allow_missing_channels": True})


def test_raw_channels_acceptor_nan_also_forbidden():
    A = make_master("F1")
    bad = make_master("F2")
    bad.loc[bad.index[5], "Acceptor_Raw_kinetic"] = float('nan')
    rep = merge.classify_sources([MergeSource("A", A), MergeSource("B", bad)])
    assert merge.RAW_CHANNELS_REQUIRED in codes(rep.forbidden)


def test_raw_channels_clean_masters_pass():
    A = make_master("F1")
    B = make_master("F2")
    rep = merge.classify_sources([MergeSource("A", A), MergeSource("B", B)])
    assert merge.RAW_CHANNELS_REQUIRED not in codes(rep.forbidden)


# --------------------------------------------------------------------------- #
# FILENAME_DUPLICATE
# --------------------------------------------------------------------------- #
def test_filename_duplicate_agree_downgraded_to_info():
    A = make_master("F1", seed=3)
    rep = merge.classify_sources([MergeSource("A", A), MergeSource("A2", A.copy())])
    assert not rep.forbidden
    dup = find(rep.auto, merge.FILENAME_DUPLICATE)
    assert dup and dup[0]["severity"] == merge.INFO


def test_filename_exclusion_conflict_needs_input():
    A = make_master("F1", seed=3, excluded_wells=())
    B = A.copy()
    # Same content-identity key (raw data identical) but DIFFERENT exclusion state.
    B.loc[B["Well_ID"] == "A1", "Is_Excluded"] = True
    rep = merge.classify_sources([MergeSource("A", A), MergeSource("B", B)])
    conf = find(rep.needs_input, merge.FILENAME_EXCLUSION_CONFLICT)
    assert conf, "expected FILENAME_EXCLUSION_CONFLICT on exclusion disagreement"
    assert conf[0]["context"]["resolution_key"] == "dup_keep::F1"
    # The identical-data case is NOT a forbidden data collision.
    assert merge.FILENAME_DATA_COLLISION not in codes(rep.forbidden)


def test_filename_exclusion_conflict_resolution_picks_source():
    A = make_master("F1", seed=3)
    B = A.copy()
    B.loc[B["Well_ID"] == "A1", "Is_Excluded"] = True
    # Choose B's exclusion state to win.
    m, rep = merge.merge_masters([MergeSource("A", A), MergeSource("B", B)],
                                 {"dup_keep::F1": "B"})
    assert m["File_Name"].nunique() == 1                 # deduped to one file
    a1_excluded = m.loc[m["Well_ID"] == "A1", "Is_Excluded"].all()
    assert a1_excluded, "B's exclusion of A1 should have won"


def test_filename_exclusion_conflict_default_keeps_first():
    A = make_master("F1", seed=3)                         # A1 not excluded
    B = A.copy()
    B.loc[B["Well_ID"] == "A1", "Is_Excluded"] = True
    # No dup_keep resolution -> default keeps first (A), so A1 stays NOT excluded.
    m, rep = merge.merge_masters([MergeSource("A", A), MergeSource("B", B)])
    assert m["File_Name"].nunique() == 1
    assert not m.loc[m["Well_ID"] == "A1", "Is_Excluded"].any()


# --------------------------------------------------------------------------- #
# CONDITION_PARTITION_OVERLAP — deferred (hook is currently a no-op)
# --------------------------------------------------------------------------- #
def test_condition_partition_overlap_is_deferred_noop():
    # Two sources record the same biology under different Main_Plasmids/Transfection
    # splits. Canonicalization is deferred, so the hook must emit NOTHING and the merge
    # must proceed normally (splits kept distinct).
    A = make_master("F1", cell="HEK", ligand="ATP", main="pA+pB",
                    conds={1: "pC", 2: "pC", 3: "pC"})
    B = make_master("F2", cell="HEK", ligand="ATP", main="pA",
                    conds={1: "pB+pC", 2: "pB+pC", 3: "pB+pC"})
    rep = merge.classify_sources([MergeSource("A", A), MergeSource("B", B)])
    assert merge.CONDITION_PARTITION_OVERLAP not in codes(rep.all_issues())
    m, rep2 = merge.merge_masters([MergeSource("A", A), MergeSource("B", B)])
    assert len(m) == len(A) + len(B)
    assert not rep2.forbidden


# --------------------------------------------------------------------------- #
# SCHEMA_MIGRATED
# --------------------------------------------------------------------------- #
def test_schema_migrated_info():
    A = make_master("F1")
    legacy = make_master("F2")
    # Simulate an older master: drop a post-v1 column so ensure_master_csv_schema
    # has to backfill it (was_modified -> True).
    legacy = legacy.drop(columns=["PR_Time(min)"])
    rep = merge.classify_sources([MergeSource("A", A), MergeSource("legacy", legacy)])
    mig = find(rep.auto, merge.SCHEMA_MIGRATED)
    assert mig, "expected SCHEMA_MIGRATED info for the legacy source"
    assert mig[0]["severity"] == merge.INFO


# --------------------------------------------------------------------------- #
# NEW_CONDITIONS
# --------------------------------------------------------------------------- #
def test_new_conditions_count():
    A = make_master("F1", cell="HEK", ligand="ATP", main="pA")
    B = make_master("F2", cell="COS", ligand="DA", main="pB")
    rep = merge.classify_sources([MergeSource("A", A), MergeSource("B", B)])
    nc = find(rep.auto, merge.NEW_CONDITIONS)
    assert nc and nc[0]["context"]["count"] == 1   # B contributes exactly one new condition


# --------------------------------------------------------------------------- #
# EXTRA_COLUMNS_DROPPED
# --------------------------------------------------------------------------- #
def test_extra_columns_dropped_info():
    A = make_master("F1")
    B = make_master("F2")
    B = B.copy()
    B["Bogus_Column"] = "junk"
    rep = merge.classify_sources([MergeSource("A", A), MergeSource("B", B)])
    ex = find(rep.auto, merge.EXTRA_COLUMNS_DROPPED)
    assert ex, "expected EXTRA_COLUMNS_DROPPED info"
    assert "Bogus_Column" in ex[0]["context"]["columns"]
    # And the merged frame is strictly schema-shaped.
    m, _ = merge.merge_masters([MergeSource("A", A), MergeSource("B", B)])
    assert list(m.columns) == MASTER_COLUMNS


# --------------------------------------------------------------------------- #
# Acceptance #1 — self-merge identity
# --------------------------------------------------------------------------- #
def test_self_merge_identity():
    A = make_master("F1", seed=7)
    m, rep = merge.merge_masters([MergeSource("A", A), MergeSource("A2", A.copy())])
    assert not rep.forbidden
    # one FILENAME_DUPLICATE -> info per file
    dup = find(rep.auto, merge.FILENAME_DUPLICATE)
    assert len(dup) == 1 and dup[0]["severity"] == merge.INFO
    # All files dedupe -> same row count and same files as the single source.
    assert len(m) == len(A)
    assert sorted(m["File_Name"].unique()) == sorted(A["File_Name"].unique())
    assert rep.summary["duplicates_dropped"] == 1


# --------------------------------------------------------------------------- #
# Acceptance #3 — disjoint concat
# --------------------------------------------------------------------------- #
def test_disjoint_concat_len_and_conditions():
    A = make_master("F1", cell="HEK", ligand="ATP", main="pA")
    B = make_master("F2", cell="COS", ligand="DA", main="pB")
    m, rep = merge.merge_masters([MergeSource("A", A), MergeSource("B", B)])
    assert len(m) == len(A) + len(B)
    assert not rep.forbidden
    assert rep.summary["conditions_added"] == 1
    assert set(m["File_Name"].unique()) == {"F1", "F2"}


# --------------------------------------------------------------------------- #
# Acceptance #4 — rename + recompute sanity on the merged frame
# --------------------------------------------------------------------------- #
def test_rename_then_recompute_runs_cleanly_on_merged_frame():
    A = make_master("F1", seed=1)
    B = make_master("F1", seed=999)          # collision (different data)
    m, rep = merge.merge_masters([MergeSource("A", A), MergeSource("B", B)],
                                 {"filename_collision": "rename"})
    assert sorted(m["File_Name"].unique()) == ["F1", "F1#B"]

    config = ProcessingConfig()
    # Flip an excluded well on the renamed file and recompute just that file.
    target_file = "F1#B"
    well = "B1"
    m2 = m.copy()
    m2.loc[(m2["File_Name"] == target_file) & (m2["Well_ID"] == well),
           "Is_Excluded"] = True
    out = __import__("processing").recompute_master_after_exclusion(
        m2, [target_file], config)
    # Excluded well's kinetic value is NaN'd; the OTHER file is untouched.
    excl_vals = out.loc[(out["File_Name"] == target_file) & (out["Well_ID"] == well),
                        "Raw_BRET_kinetic"]
    assert excl_vals.isna().all()
    untouched = out.loc[out["File_Name"] == "F1", "Raw_BRET_kinetic"]
    assert untouched.notna().any()       # the non-recomputed file kept its data


# --------------------------------------------------------------------------- #
# Acceptance #5 — exclusion-rule scoping
# --------------------------------------------------------------------------- #
def _manual_row_rule(date="21.05.26", cell="HEK", ligand="ATP", cond="WT", row="A"):
    """Build a manual row-rule blob token matching restore's canonical format."""
    return (f"Ligand: {ligand} | Date: {date} | Cell: {cell} | "
            f"Cond: {cond} | Rep: | Row:{row}")


def test_manual_rule_scopes_by_discriminator_no_bleed():
    # Source A: Main=pA, row A wells excluded by a manual Row:A rule.
    # Source B: Main=pB, same cell/ligand/cond/date/row -> would bleed without a
    # discriminator. The scoper should inject Main: pA (one readable rule).
    blobA = _manual_row_rule(row="A")
    A = make_master("F1", main="pA", excluded_wells=("A1", "A2", "A3"), blob=blobA)
    B = make_master("F2", main="pB")     # same HEK/ATP/WT, row A present, NOT excluded

    m, rep = merge.merge_masters([MergeSource("A", A), MergeSource("B", B)])

    # The merged blob must re-resolve A's rule to EXACTLY A's three origin wells.
    active = restore.list_active_exclusions(m)
    # Exactly one rule entry, covering only F1's row-A wells.
    assert len(active) == 1
    wells = set(active[0]["wells"])
    assert wells == {("F1", "A1"), ("F1", "A2"), ("F1", "A3")}
    # No bleed onto source B.
    assert all(f == "F1" for (f, _) in wells)
    # The rule stayed ONE readable manual rule (not decomposed into per-well tokens):
    # its label is a single 'Ligand: ... | Row:A | Main: pA' token, not three AUTO tokens.
    assert "AUTO:" not in active[0]["label"]
    assert "Main: pA" in active[0]["label"]


def test_scoped_manual_rule_targeted_restore_round_trips():
    blobA = _manual_row_rule(row="A")
    A = make_master("F1", main="pA", excluded_wells=("A1", "A2", "A3"), blob=blobA)
    B = make_master("F2", main="pB")
    m, _ = merge.merge_masters([MergeSource("A", A), MergeSource("B", B)])

    config = ProcessingConfig()
    active = restore.list_active_exclusions(m)
    label = active[0]["label"]
    before_excluded = int(m["Is_Excluded"].sum())
    assert before_excluded == 12          # 3 wells x 4 timepoints (per-row flag)

    report = restore.restore_rule(m, label, config)
    # restore_rule mutates `m` in place and returns a report scoreboard.
    assert report["restored_count"] == 3         # 3 wells (well-granular count)
    # All three of A's wells restored; B never touched.
    assert int(m["Is_Excluded"].sum()) == 0
    assert restore.list_active_exclusions(m) == []


def test_manual_rule_does_not_bleed_when_b_also_has_excluded_elsewhere():
    # A excludes row A (pA); B independently excludes row C (pB). After merge each
    # rule must resolve to its own source only.
    A = make_master("F1", main="pA", excluded_wells=("A1", "A2", "A3"),
                    blob=_manual_row_rule(row="A"))
    B = make_master("F2", main="pB", excluded_wells=("C1", "C2", "C3"),
                    blob=_manual_row_rule(row="C"))
    m, _ = merge.merge_masters([MergeSource("A", A), MergeSource("B", B)])
    active = restore.list_active_exclusions(m)
    by_file = {}
    for entry in active:
        for (f, w) in entry["wells"]:
            by_file.setdefault(f, set()).add(w)
    assert by_file["F1"] == {"A1", "A2", "A3"}
    assert by_file["F2"] == {"C1", "C2", "C3"}


def test_pathological_overlap_forces_per_well_decomposition():
    # Directly exercise restore.scope_blob_for_merge with a provisional merged frame
    # where origin and bleed rows are IDENTICAL in every meaningful field AND share a
    # File_Name, and the bleed well is itself excluded -> File scoping over-resolves,
    # so the scoper must fall back to per-well File-pinned tokens for THIS rule only.
    base = make_master("SHARED", main="pX", cell="HEK", ligand="ATP",
                       excluded_wells=("A1", "A2"),
                       blob=_manual_row_rule(row="A"))
    merged = base.copy().reset_index(drop=True)
    # Tag A1's rows as origin, A2's rows as bleed (same file, same everything else).
    origin_mask = merged["Well_ID"] == "A1"

    tokens = restore.scope_blob_for_merge(merged, _manual_row_rule(row="A"),
                                          origin_mask)
    # Fallback emits per-well tokens (AUTO well-pinned), NOT a broad manual rule.
    assert tokens, "scoper returned no tokens"
    assert all(t.startswith("AUTO:") for t in tokens)
    # And they resolve to exactly A1 within the origin (the bleed well A2 is excluded
    # only by its own source's tokens, which we did not request here).
    ctx = restore.build_resolve_ctx(merged)
    resolved = set()
    for t in tokens:
        for e in restore.parse_exclusion_blob(t):
            resolved |= ctx.rule_wells(e, only_excluded=True)
    assert resolved == {("SHARED", "A1")}


def test_whole_date_rule_with_cell_discriminator_stays_one_rule():
    # Spec example: A's Date rule where B measured the same date but a different cell
    # line -> scopes to Date + Cell, ONE readable rule, no file names.
    blobA = "Ligand: ATP | Date: 21.05.26 | Cell: HEK | Cond: All | Rep: | Row:"
    A = make_master("F1", date="2026-05-21", cell="HEK", main="pA",
                    excluded_wells=tuple(f"{r}{c}" for r in ROWS for c in (1, 2, 3)),
                    blob=blobA)
    B = make_master("F2", date="2026-05-21", cell="COS", main="pB")  # same date, other cell
    m, _ = merge.merge_masters([MergeSource("A", A), MergeSource("B", B)])
    active = restore.list_active_exclusions(m)
    # One rule, all wells in F1, none in F2, and no File: token / no AUTO decomposition.
    assert len(active) == 1
    assert all(f == "F1" for (f, _) in active[0]["wells"])
    assert "AUTO:" not in active[0]["label"]
    assert "File:" not in active[0]["label"]


def test_date_discriminates_date_agnostic_rule_across_sources():
    # A date-AGNOSTIC manual rule (Date: All) over two sources identical in every other
    # field but measured on different days. Date is the only thing that separates them, so
    # the scoper must inject the origin's Date (normalized %d.%m.%y) and keep ONE readable
    # rule — never falling through to File scoping or per-well decomposition.
    blobA = "Ligand: ATP | Date: All | Cell: HEK | Cond: WT | Rep: | Row:A"
    A = make_master("F1", date="2026-05-21", cell="HEK", ligand="ATP", main="pA",
                    excluded_wells=("A1", "A2", "A3"), blob=blobA)
    B = make_master("F2", date="2026-05-22", cell="HEK", ligand="ATP", main="pA")  # other day
    m, _ = merge.merge_masters([MergeSource("A", A), MergeSource("B", B)])
    active = restore.list_active_exclusions(m)
    assert len(active) == 1
    assert set(active[0]["wells"]) == {("F1", "A1"), ("F1", "A2"), ("F1", "A3")}
    # Scoped by Date, as one readable manual rule (not File-pinned, not per-well AUTO).
    assert "AUTO:" not in active[0]["label"]
    assert "File:" not in active[0]["label"]
    assert "Date: 21.05.26" in active[0]["label"]
def test_folder_source_passes_raw_channel_gate():
    # A freshly processed folder source is enriched by construction; emulate it with a
    # clean master and kind="folder". Empty Info_Sheet must not trip the gate.
    A = make_master("F1")
    A["Info_Sheet"] = ""
    rep = merge.classify_sources([MergeSource("folder", A, kind="folder")])
    assert merge.RAW_CHANNELS_REQUIRED not in codes(rep.forbidden)


# --------------------------------------------------------------------------- #
# Misc: report structure + summary integrity
# --------------------------------------------------------------------------- #
def test_report_summary_keys_present():
    A = make_master("F1", cell="HEK", ligand="ATP", main="pA")
    B = make_master("F2", cell="COS", ligand="DA", main="pB")
    m, rep = merge.merge_masters([MergeSource("A", A), MergeSource("B", B)])
    for k in ("total_sources", "total_files", "total_rows",
              "conditions_added", "collisions_handled", "duplicates_dropped"):
        assert k in rep.summary
    assert rep.summary["total_sources"] == 2
    assert rep.summary["total_files"] == 2
    assert rep.summary["total_rows"] == len(A) + len(B)


def test_single_source_merge_is_noop_like():
    A = make_master("F1")
    m, rep = merge.merge_masters([MergeSource("A", A)])
    assert len(m) == len(A)
    assert not rep.forbidden
    assert list(m.columns) == MASTER_COLUMNS


def test_merge_forbidden_carries_report():
    A = make_master("F1", times=(0.0, 1.0, 2.0))
    B = make_master("F2", times=(0.0, 1.0, 2.0, 3.0))
    with pytest.raises(MergeForbidden) as ei:
        merge.merge_masters([MergeSource("A", A), MergeSource("B", B)])
    assert ei.value.report is not None
    assert merge.TIME_VECTOR_DIVERGENCE in codes(ei.value.report.forbidden)