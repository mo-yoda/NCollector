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
import exclusions
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


# Backbone resolution helpers --------------------------------------------------
# Under the current definition, any two OVERLAPPING sources (sharing >=1 plasmid token)
# must confirm a single Main_Plasmids backbone before merging — even when they already
# agree on Main_Plasmids. Tests that aren't about canonicalization pass the (matching)
# backbone so canon is a no-op and the rest of the assertions are unchanged.
def _bb(*tokens, **extra):
    """resolutions dict selecting `tokens` as the Main_Plasmids backbone."""
    return {"main_plasmids": list(tokens), **extra}


def _concat(*dfs):
    return pd.concat(dfs, ignore_index=True)


def _scope_tokens(merged, blob, origin_files):
    """Scope ONE source's blob against a provisional merged frame, via the restore engine
    (the function merge_masters delegates per-source exclusion scoping to)."""
    origin = merged["File_Name"].astype(str).isin([str(f) for f in origin_files])
    return exclusions.scope_blob_for_merge(merged, blob, origin)


def _resolve_wells(merged, tokens):
    """(File, Well) set that `tokens` re-resolve to (criteria ∩ Is_Excluded)."""
    ctx = exclusions.build_resolve_ctx(merged)
    wells = set()
    for t in tokens:
        for e in exclusions.parse_exclusion_blob(t):
            wells |= ctx.rule_wells(e, only_excluded=True)
    return wells


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
    # Identical copies share the same Main_Plasmids -> no backbone prompt; not a collision.
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
# CONDITION_PARTITION_OVERLAP / Main_Plasmids canonicalization
# (Was test_condition_partition_overlap_is_deferred_noop, which asserted the feature
#  was deferred. This implements the feature, so the test now exercises it. The old
#  fixture used "pA+pB" without the " + " spacing the compiler emits, so it incidentally
#  dodged tokenization; with proper spacing the two splits ARE the same biology.)
# --------------------------------------------------------------------------- #
def test_condition_partition_overlap_now_canonicalizes():
    # Two sources record the SAME biology (HEK/ATP, tokens {pA,pB,pC}) under different
    # Main_Plasmids/Transfection splits. Unresolved -> forbidden; with a backbone chosen
    # -> reconciled into one condition.
    A = make_master("F1", cell="HEK", ligand="ATP", main="pA + pB",
                    conds={1: "pC", 2: "pC", 3: "pC"})
    B = make_master("F2", cell="HEK", ligand="ATP", main="pA",
                    conds={1: "pB + pC", 2: "pB + pC", 3: "pB + pC"})
    rep = merge.classify_sources([MergeSource("A", A), MergeSource("B", B)])
    assert rep.canon_plan.requires_selection
    assert merge.MULTIPLE_MAIN_PLASMIDS in codes(rep.forbidden)
    # Unresolved without a backbone:
    with pytest.raises(MergeForbidden):
        merge.merge_masters([MergeSource("A", A), MergeSource("B", B)])
    # Resolved with backbone "pA": both collapse to Main="pA", Transf="pB + pC".
    m, rep2 = merge.merge_masters([MergeSource("A", A), MergeSource("B", B)],
                                  {"main_plasmids": ["pA"]})
    assert not rep2.forbidden
    assert set(m["Main_Plasmids"].unique()) == {"pA"}
    assert set(m["Transfection"].unique()) == {"pB + pC"}
    assert len(m) == len(A) + len(B)


# --------------------------------------------------------------------------- #
# Main_Plasmids canonicalization — spec acceptance tests
# --------------------------------------------------------------------------- #
def _canon_pair(excluded_wells=(), blobA="None"):
    """The spec's b2AR example: A and B encode one biology (HEK/ATP, tokens
    {b2AR, D44KE, CAMYEL}) under different Main/Transfection splits."""
    A = make_master("F1", cell="HEK", ligand="ATP", main="b2AR + D44KE",
                    conds={1: "CAMYEL", 2: "CAMYEL", 3: "CAMYEL"},
                    excluded_wells=excluded_wells, blob=blobA)
    B = make_master("F2", cell="HEK", ligand="ATP", main="b2AR",
                    conds={1: "CAMYEL + D44KE", 2: "CAMYEL + D44KE", 3: "CAMYEL + D44KE"})
    return A, B


def test_canon_select_b2ar_collapses_to_one_condition():
    # (spec test 1)
    A, B = _canon_pair()
    m, rep = merge.merge_masters([MergeSource("A", A), MergeSource("B", B)],
                                 {"main_plasmids": ["b2AR"]})
    assert set(m["Main_Plasmids"].unique()) == {"b2AR"}
    assert set(m["Transfection"].unique()) == {"CAMYEL + D44KE"}
    # One biological condition (HEK, ATP, b2AR, CAMYEL + D44KE), two files.
    bio = m[["Main_Plasmids", "Transfection", "Cell_Line", "Ligand"]].drop_duplicates()
    assert len(bio) == 1
    assert m["File_Name"].nunique() == 2
    assert len(m) == len(A) + len(B)


def test_canon_select_b2ar_d44ke_keeps_d44ke_in_backbone():
    # (spec test 2)
    A, B = _canon_pair()
    m, rep = merge.merge_masters([MergeSource("A", A), MergeSource("B", B)],
                                 {"main_plasmids": ["b2AR", "D44KE"]})
    assert set(m["Main_Plasmids"].unique()) == {"b2AR + D44KE"}
    assert set(m["Transfection"].unique()) == {"CAMYEL"}


def test_canon_plan_requires_selection_and_orders_candidates():
    # (spec test 3) Add a third source so "b2AR" appears in 2 distinct identities while
    # the others appear in 1 -> a meaningful descending-occurrence ordering.
    A, B = _canon_pair()
    C = make_master("F3", cell="HEK", ligand="DA", main="b2AR",
                    conds={1: "GsX", 2: "GsX", 3: "GsX"})
    plan = merge.plan_canonicalization([A, B, C])
    assert plan.requires_selection
    assert plan.selected_main is None
    assert not plan.rewrites                      # no rewrites until a backbone is chosen
    # Candidates ordered by descending MEASUREMENT (file) count: b2AR is in all 3 files,
    # CAMYEL/D44KE in 2 (F1,F2), GsX in 1 (F3) -> b2AR ranks first.
    assert plan.candidate_tokens[0] == ("b2AR", 3)
    others = {tok for tok, _ in plan.candidate_tokens[1:]}
    assert others == {"CAMYEL", "D44KE", "GsX"}
    # Default backbone = the top candidate (most measurements, then alphabetical).
    assert plan.preselected == ["b2AR"]


def test_canon_superset_backbone_requires_selection():
    # Regression for the real-world report: A="b2AR-nLuc", B="b2AR-nLuc + miniG" where B's
    # conditions carry MORE total plasmid content than A's (so there is NO identical-content
    # collision), yet the overlapping backbone still requires a single-backbone selection.
    A = make_master("F1", cell="HEK", ligand="ISO", main="b2AR-nLuc",
                    conds={1: "GRK2", 2: "GRK2", 3: "GRK2"})
    B = make_master("F2", cell="HEK", ligand="ISO", main="b2AR-nLuc + miniG",
                    conds={1: "GRK2", 2: "GRK2", 3: "GRK2"})
    plan = merge.plan_canonicalization([A, B])
    assert plan.requires_selection, "overlapping-backbone mismatch must prompt"
    assert plan.preselected == ["b2AR-nLuc"]            # the shared backbone is the default
    # Unresolved without a backbone choice.
    rep = merge.classify_sources([MergeSource("A", A), MergeSource("B", B)])
    assert merge.MULTIPLE_MAIN_PLASMIDS in codes(rep.forbidden)
    with pytest.raises(MergeForbidden):
        merge.merge_masters([MergeSource("A", A), MergeSource("B", B)])
    # Choosing "b2AR-nLuc" unifies the backbone; miniG moves into B's Transfection.
    m, _ = merge.merge_masters([MergeSource("A", A), MergeSource("B", B)],
                               {"main_plasmids": ["b2AR-nLuc"]})
    assert set(m["Main_Plasmids"].unique()) == {"b2AR-nLuc"}
    assert set(m[m["File_Name"] == "F1"]["Transfection"].unique()) == {"GRK2"}
    assert set(m[m["File_Name"] == "F2"]["Transfection"].unique()) == {"GRK2 + miniG"}


def test_canon_overlap_via_shared_transfection_prompts_with_default():
    # Main_Plasmids differ (X vs Y) and the ONLY shared token is a Transfection (T) present
    # in both sources. Treating Main and Transfection tokens at the same level, this overlaps
    # -> prompt (not forbidden), and the default backbone is the shared token T (the only
    # backbone choice that keeps both sources). Regression for an empty default that broke the
    # merge preview.
    A = make_master("F1", cell="HEK", ligand="ATP", main="X", conds={1: "T", 2: "T", 3: "T"})
    B = make_master("F2", cell="COS", ligand="DA", main="Y", conds={1: "T", 2: "T", 3: "T"})
    plan = merge.plan_canonicalization([A, B])
    assert plan.requires_selection
    assert plan.preselected == ["T"]                    # shared token is the default backbone
    rep = merge.classify_sources([MergeSource("A", A), MergeSource("B", B)])
    assert merge.MULTIPLE_MAIN_PLASMIDS in codes(rep.forbidden)
    assert merge.NO_PLASMID_OVERLAP not in codes(rep.forbidden)
    # Picking T as the backbone keeps both sources (X and Y move into Transfection).
    m, _ = merge.merge_masters([MergeSource("A", A), MergeSource("B", B)],
                               {"main_plasmids": ["T"]})
    assert set(m["Main_Plasmids"].unique()) == {"T"}
    assert set(m["File_Name"]) == {"F1", "F2"}
    assert set(m[m["File_Name"] == "F1"]["Transfection"].unique()) == {"X"}
    assert set(m[m["File_Name"] == "F2"]["Transfection"].unique()) == {"Y"}


def test_canon_cross_role_overlap_default_backbone():
    # "b2AR" is Main_Plasmids in A but a Transfection in B -> overlap via b2AR; the default
    # backbone must be b2AR (not empty), so the merge preview can render.
    A = make_master("F1", main="b2AR", conds={1: "GRK2", 2: "GRK2", 3: "GRK2"})
    B = make_master("F2", main="miniG", conds={1: "b2AR", 2: "b2AR", 3: "b2AR"})
    plan = merge.plan_canonicalization([A, B])
    assert plan.requires_selection
    assert plan.preselected == ["b2AR"]


def test_no_plasmid_overlap_forbidden():
    # Definition rule 1: sources sharing NO plasmid token at all (across Main_Plasmids ∪
    # Transfection) are different experiments -> forbidden, no override.
    A = make_master("F1", cell="HEK", ligand="ATP", main="pA",
                    conds={1: "X", 2: "X", 3: "X"})
    B = make_master("F2", cell="COS", ligand="DA", main="pB",
                    conds={1: "Y", 2: "Y", 3: "Y"})
    rep = merge.classify_sources([MergeSource("A", A), MergeSource("B", B)])
    assert merge.NO_PLASMID_OVERLAP in codes(rep.forbidden)
    assert not rep.canon_plan.requires_selection        # forbidden, not a selectable prompt
    with pytest.raises(MergeForbidden) as ei:
        merge.merge_masters([MergeSource("A", A), MergeSource("B", B)])
    assert merge.NO_PLASMID_OVERLAP in codes(ei.value.issues)
    # No override (even supplying a backbone cannot force it).
    with pytest.raises(MergeForbidden):
        merge.merge_masters([MergeSource("A", A), MergeSource("B", B)], _bb("pA"))


def test_canon_matching_main_no_prompt():
    # Definition rule 1: when the sources already share the SAME Main_Plasmids backbone there
    # is nothing to reconcile -> NO prompt, no forbidden, merge proceeds untouched. Token
    # order/casing in the Main_Plasmids string does not matter (set-based comparison).
    A = make_master("F1", cell="HEK", main="b2AR-NanoLuc + bArr2")
    B = make_master("F2", cell="COS", main="bArr2 + b2AR-NanoLuc")   # same backbone, reordered
    plan = merge.plan_canonicalization([A, B])
    assert not plan.requires_selection                  # mains match -> no prompt
    m, rep = merge.merge_masters([MergeSource("A", A), MergeSource("B", B)])
    assert not rep.forbidden                            # no MULTIPLE_MAIN_PLASMIDS gate
    assert merge.MULTIPLE_MAIN_PLASMIDS not in codes(rep.forbidden)
    assert len(m) == len(A) + len(B)                    # merged as-is, nothing dropped
    # apply is a no-op when no backbone is carried.
    out = merge.apply_canonicalization(A, merge.plan_canonicalization(A))
    assert out.equals(A)


def test_canon_exclusion_survival_on_swapped_condition():
    # (spec test 5) An excluded well on A's swapped condition must re-resolve to the SAME
    # (File, Well) set in the post-merge unified blob.
    blobA = ("Ligand: ATP | Date: 21.05.26 | Cell: HEK | "
             "Cond: CAMYEL | Rep: | Row:A")
    A, B = _canon_pair(excluded_wells=("A1", "A2", "A3"), blobA=blobA)
    m, _ = merge.merge_masters([MergeSource("A", A), MergeSource("B", B)],
                               {"main_plasmids": ["b2AR"]})
    # Is_Excluded preserved (3 wells x 4 timepoints), only on F1.
    assert int(m["Is_Excluded"].sum()) == 12
    excl_files = set(m.loc[m["Is_Excluded"], "File_Name"])
    assert excl_files == {"F1"}
    # The unified blob re-resolves to exactly A's three origin wells (no bleed to F2).
    active = exclusions.list_active_exclusions(m)
    wells = set()
    for e in active:
        wells |= set(e["wells"])
    assert wells == {("F1", "A1"), ("F1", "A2"), ("F1", "A3")}


def test_canon_excluded_wells_revert_after_merge():
    # (exclusion survival, restore round-trip) After canonicalization strands the original
    # manual rule, the backstop keeps the wells excluded AND revertable: restoring every
    # active rule clears Is_Excluded and recomputes cleanly.
    blobA = ("Ligand: ATP | Date: 21.05.26 | Cell: HEK | "
             "Cond: CAMYEL | Rep: | Row:A")
    A, B = _canon_pair(excluded_wells=("A1", "A2", "A3"), blobA=blobA)
    m, _ = merge.merge_masters([MergeSource("A", A), MergeSource("B", B)],
                               {"main_plasmids": ["b2AR"]})
    assert int(m["Is_Excluded"].sum()) == 12

    config = ProcessingConfig()
    for entry in list(exclusions.list_active_exclusions(m)):
        exclusions.restore_rule(m, entry["label"], config)
    # All exclusions reverted; the canonical raw came back.
    assert int(m["Is_Excluded"].sum()) == 0
    assert exclusions.list_active_exclusions(m) == []
    restored = m.loc[(m["File_Name"] == "F1") & (m["Well_ID"] == "A1"), "Raw_BRET_kinetic"]
    assert restored.notna().any()


def test_canon_canonicalizes_labeling_control_rows():
    # A labeling-control column (marker in Replicate, real condition in Transfection) belongs
    # to a condition and must be canonicalized like the rest of it: its Main_Plasmids /
    # Transfection are rewritten, while the Replicate marker is preserved (canon never touches
    # Replicate). Otherwise the control column desyncs from its block — each unique
    # (Main_Plasmids + Transfection) per cell line has its own labeling correction.
    A = make_master("F1", cell="HEK", ligand="ATP", main="b2AR + D44KE",
                    conds={1: "CAMYEL", 2: "CAMYEL", 3: "CAMYEL"})
    B = make_master("F2", cell="HEK", ligand="ATP", main="b2AR",
                    conds={1: "CAMYEL + D44KE", 2: "CAMYEL + D44KE", 3: "CAMYEL + D44KE"})
    # Mark column 3 of each as the labeling-control column (real condition kept in Transfection).
    for df in (A, B):
        ctrl = df["Well_ID"].astype(str).str.endswith("3")
        df.loc[ctrl, "Replicate"] = "labeling control"

    m, _ = merge.merge_masters([MergeSource("A", A), MergeSource("B", B)],
                               {"main_plasmids": ["b2AR"]})
    ctrl_rows = m[m["Replicate"] == "labeling control"]
    assert not ctrl_rows.empty                                  # marker preserved
    # Control rows carry the canonical backbone + transfection (same as their condition)...
    assert set(ctrl_rows["Main_Plasmids"].unique()) == {"b2AR"}
    assert set(ctrl_rows["Transfection"].unique()) == {"CAMYEL + D44KE"}
    # ...and the merged master is internally consistent (controls match the other replicates).
    assert set(m["Main_Plasmids"].unique()) == {"b2AR"}
    assert set(m["Transfection"].unique()) == {"CAMYEL + D44KE"}


def test_canon_unresolved_condition_dropped_as_no_metadata():
    # (spec test 6) An OVERLAPPING source whose condition cannot honor the SELECTED backbone
    # (lacks one of the chosen tokens) is dropped as no-metadata (CONDITION_PARTITION_OVERLAP),
    # while the merge itself proceeds. A and B share b2AR + GRK2 (overlap OK); the chosen
    # backbone "b2AR + miniG" includes miniG, which B does not have -> B is dropped.
    A = make_master("F1", cell="HEK", ligand="ATP", main="b2AR + miniG",
                    conds={1: "GRK2", 2: "GRK2", 3: "GRK2"})
    B = make_master("F2", cell="HEK", ligand="ATP", main="b2AR",
                    conds={1: "GRK2", 2: "GRK2", 3: "GRK2"})
    m, rep = merge.merge_masters([MergeSource("A", A), MergeSource("B", B)],
                                 {"main_plasmids": ["b2AR", "miniG"]})
    # B is gone from the merged master (its conditions can't honor the miniG backbone).
    assert set(m["File_Name"]) == {"F1"}
    assert merge.CONDITION_PARTITION_OVERLAP in codes(rep.all_issues())
    assert set(m["Main_Plasmids"].unique()) == {"b2AR + miniG"}


def test_canon_selection_withheld_is_forbidden():
    # (spec test 7)
    A, B = _canon_pair()
    with pytest.raises(MergeForbidden) as ei:
        merge.merge_masters([MergeSource("A", A), MergeSource("B", B)])
    assert merge.MULTIPLE_MAIN_PLASMIDS in codes(ei.value.issues)


def test_canon_case_folding_identity_and_casing_preserved():
    # (spec test 8) Tokens differing only in case are one token for identity/ordering;
    # emitted strings keep a single canonical casing.
    A = make_master("F1", cell="HEK", ligand="ATP", main="B2ar + d44ke",
                    conds={1: "camyel", 2: "camyel", 3: "camyel"})
    B = make_master("F2", cell="HEK", ligand="ATP", main="b2AR",
                    conds={1: "CAMYEL + D44KE", 2: "CAMYEL + D44KE", 3: "CAMYEL + D44KE"})
    # Case-folded identities match -> collision detected.
    plan = merge.plan_canonicalization([A, B])
    assert plan.requires_selection
    m, _ = merge.merge_masters([MergeSource("A", A), MergeSource("B", B)],
                               {"main_plasmids": ["b2AR"]})
    # Single backbone (user's casing) and a single, consistently-cased Transfection.
    assert set(m["Main_Plasmids"].unique()) == {"b2AR"}
    assert m["Transfection"].nunique() == 1            # case-insensitive collapse worked
    # One biological condition despite the mixed input casing.
    bio = m[["Main_Plasmids", "Transfection", "Cell_Line", "Ligand"]].drop_duplicates()
    assert len(bio) == 1


def test_canon_multiple_ligands_preserved():
    # (spec test 9) Different ligands are never collapsed; each keeps its own condition.
    A, B = _canon_pair()                               # ATP collision
    C = make_master("F3", cell="HEK", ligand="DA", main="b2AR + D44KE",
                    conds={1: "CAMYEL", 2: "CAMYEL", 3: "CAMYEL"})
    m, _ = merge.merge_masters(
        [MergeSource("A", A), MergeSource("B", B), MergeSource("C", C)],
        {"main_plasmids": ["b2AR"]})
    assert set(m["Ligand"].unique()) == {"ATP", "DA"}
    # One transfection per ligand, backbone unified.
    per_ligand = m.groupby("Ligand")["Transfection"].nunique()
    assert (per_ligand == 1).all()
    assert set(m["Main_Plasmids"].unique()) == {"b2AR"}


def test_canon_new_conditions_not_inflated_by_cosmetic_split():
    # (spec test 10) Post-canon NEW_CONDITIONS excludes cosmetic duplicates.
    A, B = _canon_pair()                               # same biology, different split
    rep = merge.classify_sources([MergeSource("A", A), MergeSource("B", B)])
    # Provisional canonical view (preselected backbone) collapses A and B -> no new condition.
    assert not find(rep.auto, merge.NEW_CONDITIONS)
    m, rep2 = merge.merge_masters([MergeSource("A", A), MergeSource("B", B)],
                                  {"main_plasmids": ["b2AR"]})
    assert not find(rep2.auto, merge.NEW_CONDITIONS)
    assert rep2.summary["conditions_added"] == 0


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
# Acceptance #3 — matching-backbone concat (different conditions, no prompt)
# --------------------------------------------------------------------------- #
def test_matching_main_concat_len_and_conditions():
    # Two sources sharing the SAME backbone (pX) but different cell/ligand -> no prompt;
    # nothing is dropped and the rows concatenate.
    A = make_master("F1", cell="HEK", ligand="ATP", main="pX")
    B = make_master("F2", cell="COS", ligand="DA", main="pX")
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


# NOTE: these exercise exclusions.scope_blob_for_merge directly (the function merge_masters
# delegates per-source exclusion scoping to). Since canonicalization now forces a single
# Main_Plasmids backbone on every merge, the "Main as discriminator" scenario only arises
# at the engine level, so it is unit-tested here on a hand-built provisional merged frame
# (mirroring test_pathological_overlap_forces_per_well_decomposition).
def test_manual_rule_scopes_by_discriminator_no_bleed():
    # Source A: Main=pA, row A wells excluded by a manual Row:A rule.
    # Source B: Main=pB, same cell/ligand/cond/date/row -> would bleed without a
    # discriminator. The scoper should inject Main: pA (one readable rule).
    blobA = _manual_row_rule(row="A")
    A = make_master("F1", main="pA", excluded_wells=("A1", "A2", "A3"), blob=blobA)
    B = make_master("F2", main="pB")     # same HEK/ATP/WT, row A present, NOT excluded
    merged = _concat(A, B)

    tokens = _scope_tokens(merged, blobA, ["F1"])
    # ONE readable manual rule (not per-well AUTO), scoped by Main: pA.
    assert len(tokens) == 1
    assert "AUTO:" not in tokens[0]
    assert "Main: pA" in tokens[0]
    # Re-resolves to EXACTLY A's three origin wells; no bleed onto B.
    assert _resolve_wells(merged, tokens) == {("F1", "A1"), ("F1", "A2"), ("F1", "A3")}


def test_scoped_manual_rule_targeted_restore_round_trips():
    blobA = _manual_row_rule(row="A")
    A = make_master("F1", main="pA", excluded_wells=("A1", "A2", "A3"), blob=blobA)
    B = make_master("F2", main="pB")
    merged = _concat(A, B)
    merged["Applied_Exclusions"] = " || ".join(_scope_tokens(merged, blobA, ["F1"]))

    config = ProcessingConfig()
    active = exclusions.list_active_exclusions(merged)
    assert len(active) == 1
    label = active[0]["label"]
    assert int(merged["Is_Excluded"].sum()) == 12   # 3 wells x 4 timepoints (per-row flag)

    report = exclusions.restore_rule(merged, label, config)
    assert report["restored_count"] == 3            # 3 wells (well-granular count)
    # All three of A's wells restored; B never touched.
    assert int(merged["Is_Excluded"].sum()) == 0
    assert exclusions.list_active_exclusions(merged) == []


def test_manual_rule_does_not_bleed_when_b_also_has_excluded_elsewhere():
    # A excludes row A (pA); B independently excludes row C (pB). Each rule must resolve to
    # its own source only.
    A = make_master("F1", main="pA", excluded_wells=("A1", "A2", "A3"),
                    blob=_manual_row_rule(row="A"))
    B = make_master("F2", main="pB", excluded_wells=("C1", "C2", "C3"),
                    blob=_manual_row_rule(row="C"))
    merged = _concat(A, B)
    tokens = (_scope_tokens(merged, A["Applied_Exclusions"].iloc[0], ["F1"])
              + _scope_tokens(merged, B["Applied_Exclusions"].iloc[0], ["F2"]))
    by_file = {}
    for (f, w) in _resolve_wells(merged, tokens):
        by_file.setdefault(f, set()).add(w)
    assert by_file["F1"] == {"A1", "A2", "A3"}
    assert by_file["F2"] == {"C1", "C2", "C3"}


def test_pathological_overlap_forces_per_well_decomposition():
    # Directly exercise exclusions.scope_blob_for_merge with a provisional merged frame
    # where origin and bleed rows are IDENTICAL in every meaningful field AND share a
    # File_Name, and the bleed well is itself excluded -> File scoping over-resolves,
    # so the scoper must fall back to per-well File-pinned tokens for THIS rule only.
    base = make_master("SHARED", main="pX", cell="HEK", ligand="ATP",
                       excluded_wells=("A1", "A2"),
                       blob=_manual_row_rule(row="A"))
    merged = base.copy().reset_index(drop=True)
    # Tag A1's rows as origin, A2's rows as bleed (same file, same everything else).
    origin_mask = merged["Well_ID"] == "A1"

    tokens = exclusions.scope_blob_for_merge(merged, _manual_row_rule(row="A"),
                                          origin_mask)
    # Fallback emits per-well tokens (AUTO well-pinned), NOT a broad manual rule.
    assert tokens, "scoper returned no tokens"
    assert all(t.startswith("AUTO:") for t in tokens)
    # And they resolve to exactly A1 within the origin (the bleed well A2 is excluded
    # only by its own source's tokens, which we did not request here).
    ctx = exclusions.build_resolve_ctx(merged)
    resolved = set()
    for t in tokens:
        for e in exclusions.parse_exclusion_blob(t):
            resolved |= ctx.rule_wells(e, only_excluded=True)
    assert resolved == {("SHARED", "A1")}


def test_whole_date_rule_with_discriminator_stays_one_rule():
    # A's whole-plate rule where B differs by a single discriminating field -> ONE readable
    # rule, all wells in F1, none in F2, no File: token / no AUTO decomposition.
    blobA = "Ligand: ATP | Date: 21.05.26 | Cell: HEK | Cond: All | Rep: | Row:"
    A = make_master("F1", date="2026-05-21", cell="HEK", main="pA",
                    excluded_wells=tuple(f"{r}{c}" for r in ROWS for c in (1, 2, 3)),
                    blob=blobA)
    B = make_master("F2", date="2026-05-21", cell="COS", main="pB")  # same date, other cell
    merged = _concat(A, B)
    tokens = _scope_tokens(merged, blobA, ["F1"])
    assert len(tokens) == 1
    assert "AUTO:" not in tokens[0]
    assert "File:" not in tokens[0]
    assert all(f == "F1" for (f, _) in _resolve_wells(merged, tokens))


def test_date_discriminates_date_agnostic_rule_across_sources():
    # A date-AGNOSTIC manual rule (Date: All) over two sources identical in every other
    # field but measured on different days. Date is the only thing that separates them, so
    # the scoper must inject the origin's Date (normalized %d.%m.%y) and keep ONE readable
    # rule — never falling through to File scoping or per-well decomposition.
    blobA = "Ligand: ATP | Date: All | Cell: HEK | Cond: WT | Rep: | Row:A"
    A = make_master("F1", date="2026-05-21", cell="HEK", ligand="ATP", main="pA",
                    excluded_wells=("A1", "A2", "A3"), blob=blobA)
    B = make_master("F2", date="2026-05-22", cell="HEK", ligand="ATP", main="pA")  # other day
    merged = _concat(A, B)
    tokens = _scope_tokens(merged, blobA, ["F1"])
    assert len(tokens) == 1
    assert _resolve_wells(merged, tokens) == {("F1", "A1"), ("F1", "A2"), ("F1", "A3")}
    # Scoped by Date, as one readable manual rule (not File-pinned, not per-well AUTO).
    assert "AUTO:" not in tokens[0]
    assert "File:" not in tokens[0]
    assert "Date: 21.05.26" in tokens[0]


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
    A = make_master("F1", cell="HEK", ligand="ATP", main="pX")
    B = make_master("F2", cell="COS", ligand="DA", main="pX")
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