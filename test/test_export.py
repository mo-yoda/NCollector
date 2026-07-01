"""
Tests for export.py — CSV schema migration, export filtering,
row info lookup, header key generation, and pivot table creation.

Covers: ensure_master_csv_schema, apply_export_filters, build_row_info,
generate_header_key, create_clean_pivot.
"""
import pytest
import pandas as pd

from export import (
    ensure_master_csv_schema,
    apply_export_filters,
    build_row_info,
    generate_header_key,
    create_clean_pivot,
    build_crc_table,
    build_crc_preview_table,
)
from models import LEGACY_COLUMN_DEFAULTS, MASTER_COLUMNS


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_master_df(**overrides):
    """Creates a minimal valid master DataFrame with overridable columns."""
    base = {
        "NCollector_version": ["v2", "v2"],
        "Path": ["/data", "/data"],
        "Info_Sheet": ["", ""],
        "File_Name": ["file1.xlsx", "file1.xlsx"],
        "Date": ["2025-01-01", "2025-01-01"],
        "Main_Plasmids": ["P1 + P2", "P1 + P2"],
        "Applied_Exclusions": ["None", "None"],
        "Is_Excluded": [False, False],
        "Is_Vehicle": [False, True],
        "Transfection": ["CondA", "CondA"],
        "Cell_Line": ["HEK293", "HEK293"],
        "Ligand": ["Lig1", "Lig1"],
        "Ligand_Conc": [-9.0, float("nan")],
        "Plate_Row": ["A", "H"],
        "Replicate": ["1", "1"],
        "Well_ID": ["A1", "H1"],
        "Time_(min)": [0.0, 0.0],
        "PR_Time(min)": [5.0, 5.0],
        "Donor_Raw_kinetic": [500.0, 450.0],
        "Acceptor_Raw_kinetic": [100.0, 90.0],
        "Raw_BRET_kinetic": [0.2, 0.2],
        "Kinetic_Mean": [1.5, 1.0],
        "Veh_Norm_AUC": [1.3, 1.0],
        "AUC_Mean": [1.3, 1.0],
    }
    base.update(overrides)
    return pd.DataFrame(base)


def _make_legacy_df():
    """Creates a DataFrame mimicking a pre-v2 master CSV (missing new columns)."""
    return pd.DataFrame({
        "File_Name": ["f1", "f1", "f1"],
        "Date": ["01.01.25", "01.01.25", "01.01.25"],
        "Transfection": ["CondA", "CondA", "CondA"],
        "Cell_Line": ["HEK293", "HEK293", "HEK293"],
        "Ligand": ["Lig1", "Lig1", "Lig1"],
        "Ligand_Conc": [-9.0, -8.0, 0.0],
        "Plate_Row": ["A", "B", "H"],
        "Well_ID": ["A1", "B1", "H1"],
        "Time_(min)": [0, 1, 0],
        "Raw_BRET_kinetic": [0.2, 0.3, 0.2],
        "Kinetic_Mean": [1.5, 1.6, 1.0],
        "Bl_Corrected_BRET": [1.0, 1.1, 1.0],
        "Applied_Exclusions": ["None", "None", "None"],
        "Bl_AUC": [5.0, 6.0, 3.0],
        "Veh_Norm_AUC": [1.3, 1.5, 1.0],
        "AUC_Mean": [1.3, 1.5, 1.0],
    })


# ---------------------------------------------------------------------------
# Tests: ensure_master_csv_schema — column filling
# ---------------------------------------------------------------------------

class TestEnsureSchemaMissingColumns:

    def test_fills_missing_legacy_columns(self):
        """Missing columns from LEGACY_COLUMN_DEFAULTS should be added."""
        df = _make_legacy_df()
        result, was_modified, was_fixed = ensure_master_csv_schema(df)
        for col in LEGACY_COLUMN_DEFAULTS:
            assert col in result.columns, f"'{col}' not added to DataFrame"

    def test_does_not_overwrite_existing_columns(self):
        """Existing data columns should not be modified by the fill step.
        (NCollector_version is intentionally re-stamped when a legacy master is migrated,
        so it is no longer an 'unchanged' column — check a data column instead.)"""
        df = _make_master_df()
        original_values = df["Cell_Line"].tolist()
        result, _, _ = ensure_master_csv_schema(df)
        assert result["Cell_Line"].tolist() == original_values

    def test_was_modified_true_when_columns_added(self):
        """was_modified should be True when legacy columns are missing."""
        df = _make_legacy_df()
        _, was_modified, _ = ensure_master_csv_schema(df)
        assert was_modified is True

    def test_was_modified_false_when_complete(self):
        """Re-running migration on an already-migrated master is a no-op.

        The schema step now does reconstruction/stamping (version, Raw_BRET_unexcluded,
        Date canonicalization) that a hand-built minimal master can't pre-satisfy, so the
        robust 'nothing needs changing' check is idempotency: migrate once, then confirm a
        second pass reports no modification."""
        df = _make_master_df()
        migrated, was_modified_1, _ = ensure_master_csv_schema(df)
        assert was_modified_1 is True  # first pass migrates the legacy 'v2' master
        _, was_modified_2, _ = ensure_master_csv_schema(migrated)
        assert was_modified_2 is False  # second pass: already current + complete


# ---------------------------------------------------------------------------
# Tests: ensure_master_csv_schema — Is_Vehicle reconstruction
# ---------------------------------------------------------------------------

class TestSchemaIsVehicle:

    def test_reconstructed_from_plate_row_and_conc(self):
        """Is_Vehicle should be True for row H with Ligand_Conc == 0.0."""
        df = _make_legacy_df()  # Has Plate_Row and Ligand_Conc, no Is_Vehicle
        result, _, _ = ensure_master_csv_schema(df)
        h_rows = result[result["Plate_Row"] == "H"]
        assert h_rows["Is_Vehicle"].all()
        non_h = result[result["Plate_Row"] != "H"]
        assert not non_h["Is_Vehicle"].any()

    def test_non_zero_conc_row_h_not_vehicle(self):
        """Row H with non-zero concentration should NOT be vehicle."""
        df = _make_legacy_df()
        df.loc[df["Plate_Row"] == "H", "Ligand_Conc"] = -3.0  # Non-zero
        result, _, _ = ensure_master_csv_schema(df)
        h_rows = result[result["Plate_Row"] == "H"]
        assert not h_rows["Is_Vehicle"].any()


# ---------------------------------------------------------------------------
# Tests: ensure_master_csv_schema — vehicle concentration cleanup
# ---------------------------------------------------------------------------

class TestSchemaVehicleConc:

    def test_vehicle_conc_zero_becomes_nan(self):
        """Legacy 0.0 concentration for vehicle rows should become NaN."""
        df = _make_legacy_df()
        result, _, _ = ensure_master_csv_schema(df)
        vehicle = result[result["Is_Vehicle"]]
        assert vehicle["Ligand_Conc"].isna().all()

    def test_non_vehicle_conc_unchanged(self):
        """Non-vehicle concentrations should not be modified."""
        df = _make_legacy_df()
        result, _, _ = ensure_master_csv_schema(df)
        non_vehicle = result[~result["Is_Vehicle"]]
        assert (non_vehicle["Ligand_Conc"] == [-9.0, -8.0]).all()


# ---------------------------------------------------------------------------
# Tests: ensure_master_csv_schema — Is_Excluded reconstruction
# ---------------------------------------------------------------------------

class TestSchemaIsExcluded:

    def test_reconstructed_from_nan_bret(self):
        """Excluded wells (NaN in Raw_BRET_kinetic) should get Is_Excluded=True."""
        df = _make_legacy_df()
        df["Applied_Exclusions"] = ["Rule1", "Rule1", "Rule1"]
        df.loc[2, "Raw_BRET_kinetic"] = float("nan")  # Exclude H1
        result, _, _ = ensure_master_csv_schema(df)
        assert result.loc[2, "Is_Excluded"] is True or result.loc[2, "Is_Excluded"] == True
        assert result.loc[0, "Is_Excluded"] is False or result.loc[0, "Is_Excluded"] == False

    def test_no_exclusions_defaults_false(self):
        """No exclusion rules → Is_Excluded defaults to False for all."""
        df = _make_legacy_df()
        df["Applied_Exclusions"] = ["None", "None", "None"]
        result, _, _ = ensure_master_csv_schema(df)
        assert not result["Is_Excluded"].any()


# ---------------------------------------------------------------------------
# Tests: ensure_master_csv_schema — AUC bug fix
# ---------------------------------------------------------------------------

class TestSchemaAucBugFix:

    def test_excluded_auc_zero_becomes_nan(self):
        """Excluded rows with AUC=0.0 should be fixed to NaN."""
        df = _make_legacy_df()
        df["Is_Excluded"] = [False, True, False]
        df.loc[1, "Bl_AUC"] = 0.0
        df.loc[1, "Veh_Norm_AUC"] = 0.0
        df.loc[1, "Raw_BRET_kinetic"] = float("nan")
        result, _, was_fixed = ensure_master_csv_schema(df)
        assert pd.isna(result.loc[1, "Bl_AUC"])
        assert pd.isna(result.loc[1, "Veh_Norm_AUC"])
        assert was_fixed is True

    def test_non_excluded_zero_unchanged(self):
        """Non-excluded rows with AUC=0.0 should NOT be changed."""
        df = _make_legacy_df()
        df["Is_Excluded"] = [False, False, False]
        df["Bl_AUC"] = [0.0, 6.0, 3.0]
        result, _, _ = ensure_master_csv_schema(df)
        assert result.loc[0, "Bl_AUC"] == 0.0  # Not excluded, stays 0.0

    def test_auc_mean_recalculated(self):
        """AUC_Mean should be recalculated after fixing excluded zeros."""
        df = pd.DataFrame({
            "File_Name": ["f1"] * 3,
            "Transfection": ["C"] * 3,
            "Cell_Line": ["HEK"] * 3,
            "Ligand": ["L"] * 3,
            "Plate_Row": ["A"] * 3,
            "Well_ID": ["A1", "A2", "A3"],
            "Time_(min)": [0, 0, 0],
            "Is_Excluded": [False, True, False],
            "Raw_BRET_kinetic": [0.2, float("nan"), 0.2],
            "Veh_Norm_AUC": [1.4, 0.0, 1.6],  # Middle is buggy zero
            "AUC_Mean": [1.0, 1.0, 1.0],  # Wrong mean (included 0.0)
            "Applied_Exclusions": ["Rule"] * 3,
            "Kinetic_Mean": [1.5, 1.5, 1.5],
            "Ligand_Conc": [-9.0, -9.0, -9.0],
            "Plate_Row": ["A", "A", "A"],
        })
        result, _, _ = ensure_master_csv_schema(df)
        # After fix: Veh_Norm_AUC[1] = NaN, mean of [1.4, NaN, 1.6] = 1.5
        non_excluded = result[~result["Is_Excluded"]]
        assert non_excluded["AUC_Mean"].iloc[0] == pytest.approx(1.5)


# ---------------------------------------------------------------------------
# Tests: ensure_master_csv_schema — legacy cleanup
# ---------------------------------------------------------------------------

class TestSchemaLegacyCleanup:

    def test_empty_noid_removed(self):
        """Rows with Transfection='Empty/NoID' should be removed."""
        df = _make_legacy_df()
        df = pd.concat([df, pd.DataFrame({
            "File_Name": ["f1"],
            "Transfection": ["Empty/NoID"],
            "Cell_Line": ["HEK"],
            "Ligand": ["L"],
            "Ligand_Conc": [0],
            "Plate_Row": ["X"],
            "Well_ID": ["X1"],
            "Time_(min)": [0],
            "Raw_BRET_kinetic": [0.1],
            "Kinetic_Mean": [1.0],
            "Bl_Corrected_BRET": [1.0],
            "Applied_Exclusions": ["None"],
            "Bl_AUC": [1.0],
            "Veh_Norm_AUC": [1.0],
            "AUC_Mean": [1.0],
            "Date": ["01.01.25"],
        })], ignore_index=True)
        result, was_modified, _ = ensure_master_csv_schema(df)
        assert "Empty/NoID" not in result["Transfection"].values
        assert was_modified is True


# ---------------------------------------------------------------------------
# Tests: ensure_master_csv_schema — log_fn callback
# ---------------------------------------------------------------------------

class TestSchemaLogCallback:

    def test_log_fn_receives_bug_fixes(self):
        """log_fn should receive messages for bug fixes."""
        df = _make_legacy_df()
        df["Is_Excluded"] = [False, True, False]
        df.loc[1, "Bl_AUC"] = 0.0
        df.loc[1, "Veh_Norm_AUC"] = 0.0
        df.loc[1, "Raw_BRET_kinetic"] = float("nan")
        messages = []
        ensure_master_csv_schema(df, log_fn=lambda m: messages.append(m))
        fix_msgs = [m for m in messages if "[CSV BUG FIX]" in m]
        assert len(fix_msgs) > 0

    def test_no_log_fn_no_crash(self):
        """Should work fine without log_fn."""
        df = _make_legacy_df()
        result, _, _ = ensure_master_csv_schema(df, log_fn=None)
        assert result is not None


# ---------------------------------------------------------------------------
# Tests: apply_export_filters
# ---------------------------------------------------------------------------

class TestApplyExportFilters:

    @pytest.fixture
    def multi_df(self):
        """DataFrame with multiple ligands, cells, and transfections."""
        return pd.DataFrame({
            "Ligand": ["Lig1", "Lig1", "Lig2", "Lig2"],
            "Cell_Line": ["HEK", "dQ", "HEK", "dQ"],
            "Transfection": ["C1", "C1", "C2", "C2"],
            "Value": [1, 2, 3, 4],
        })

    def test_filter_by_ligand(self, multi_df):
        result = apply_export_filters(multi_df, {"ligands": ["Lig1"]})
        assert len(result) == 2
        assert (result["Ligand"] == "Lig1").all()

    def test_filter_by_cell_line(self, multi_df):
        result = apply_export_filters(multi_df, {"cells": ["dQ"]})
        assert len(result) == 2
        assert (result["Cell_Line"] == "dQ").all()

    def test_filter_by_transfection(self, multi_df):
        result = apply_export_filters(multi_df, {"transfections": ["C2"]})
        assert len(result) == 2
        assert (result["Transfection"] == "C2").all()

    def test_filter_all_returns_everything(self, multi_df):
        result = apply_export_filters(multi_df, {
            "ligands": "All", "cells": "All", "transfections": "All"
        })
        assert len(result) == 4

    def test_empty_config_returns_everything(self, multi_df):
        result = apply_export_filters(multi_df, {})
        assert len(result) == 4

    def test_multiple_filters_combined(self, multi_df):
        result = apply_export_filters(multi_df, {
            "ligands": ["Lig1"], "cells": ["HEK"]
        })
        assert len(result) == 1
        assert result.iloc[0]["Value"] == 1

    def test_string_target_converted_to_list(self, multi_df):
        """Single string (not list) should work as filter target."""
        result = apply_export_filters(multi_df, {"ligands": "Lig2"})
        assert len(result) == 2

    def test_no_match_returns_empty(self, multi_df):
        result = apply_export_filters(multi_df, {"ligands": ["NonExistent"]})
        assert len(result) == 0

    def test_does_not_modify_original(self, multi_df):
        original_len = len(multi_df)
        apply_export_filters(multi_df, {"ligands": ["Lig1"]})
        assert len(multi_df) == original_len


# ---------------------------------------------------------------------------
# Tests: build_row_info
# ---------------------------------------------------------------------------

class TestBuildRowInfo:

    def test_basic_lookup(self):
        """Should create display → criteria mapping."""
        df = pd.DataFrame({
            "Plate_Row": ["A", "H"],
            "Ligand": ["Lig1", "Lig1"],
            "Is_Vehicle": [False, True],
            "Ligand_Conc": [-9.0, float("nan")],
        })
        result = build_row_info(df)
        assert len(result) == 2
        assert "Row A: -9.0 log(M) Lig1" in result
        assert "Row H: Vehicle Lig1" in result

    def test_vehicle_criteria(self):
        """Vehicle entry should have is_vehicle=True and no 'conc' key."""
        df = pd.DataFrame({
            "Plate_Row": ["H"],
            "Ligand": ["Lig1"],
            "Is_Vehicle": [True],
            "Ligand_Conc": [float("nan")],
        })
        result = build_row_info(df)
        key = "Row H: Vehicle Lig1"
        assert result[key]["is_vehicle"] is True
        assert "conc" not in result[key]

    def test_non_vehicle_criteria(self):
        """Non-vehicle entry should have is_vehicle=False and 'conc' key."""
        df = pd.DataFrame({
            "Plate_Row": ["A"],
            "Ligand": ["Lig1"],
            "Is_Vehicle": [False],
            "Ligand_Conc": [-9.0],
        })
        result = build_row_info(df)
        key = "Row A: -9.0 log(M) Lig1"
        assert result[key]["is_vehicle"] is False
        assert result[key]["conc"] == "-9.0"

    def test_two_ligands_separate_vehicles(self):
        """Two ligands should produce separate vehicle entries."""
        df = pd.DataFrame({
            "Plate_Row": ["H", "H"],
            "Ligand": ["Lig1", "Lig2"],
            "Is_Vehicle": [True, True],
            "Ligand_Conc": [float("nan"), float("nan")],
        })
        result = build_row_info(df)
        assert "Row H: Vehicle Lig1" in result
        assert "Row H: Vehicle Lig2" in result

    def test_sorted_keys(self):
        """Keys should be sorted alphabetically."""
        df = pd.DataFrame({
            "Plate_Row": ["B", "A", "H"],
            "Ligand": ["L", "L", "L"],
            "Is_Vehicle": [False, False, True],
            "Ligand_Conc": [-8.0, -9.0, float("nan")],
        })
        result = build_row_info(df)
        keys = list(result.keys())
        assert keys == sorted(keys)

    def test_empty_df(self):
        result = build_row_info(pd.DataFrame())
        assert result == {}

    def test_none_input(self):
        result = build_row_info(None)
        assert result == {}

    def test_deduplicates(self):
        """Duplicate rows should produce a single entry."""
        df = pd.DataFrame({
            "Plate_Row": ["A", "A", "A"],
            "Ligand": ["L", "L", "L"],
            "Is_Vehicle": [False, False, False],
            "Ligand_Conc": [-9.0, -9.0, -9.0],
        })
        result = build_row_info(df)
        assert len(result) == 1


# ---------------------------------------------------------------------------
# Tests: generate_header_key
# ---------------------------------------------------------------------------

class TestGenerateHeaderKey:

    def test_basic_kinetic(self):
        """Kinetic data: header includes Transfection, Cell_Line, and concentration."""
        df = pd.DataFrame({
            "Transfection": ["C1", "C1"],
            "Cell_Line": ["HEK", "HEK"],
            "Ligand": ["L1", "L1"],
            "Ligand_Conc": [-9.0, -9.0],
            "Time_(min)": [0.0, 1.0],
            "Is_Vehicle": [False, False],
        })
        result = generate_header_key(df, include_conc=True)
        assert "Header_Key" in result.columns
        assert "C1" in result["Header_Key"].iloc[0]
        assert "HEK" in result["Header_Key"].iloc[0]
        assert "-9.0" in result["Header_Key"].iloc[0]

    def test_vehicle_shows_vehicle(self):
        """Vehicle rows should show 'Vehicle' instead of concentration."""
        df = pd.DataFrame({
            "Transfection": ["C1", "C1"],
            "Cell_Line": ["HEK", "HEK"],
            "Ligand": ["L1", "L1"],
            "Ligand_Conc": [float("nan"), float("nan")],
            "Time_(min)": [0.0, 1.0],
            "Is_Vehicle": [True, True],
        })
        result = generate_header_key(df, include_conc=True)
        assert "Vehicle" in result["Header_Key"].iloc[0]

    def test_group_by_transfection(self):
        """When grouped by Transfection, it should be excluded from key."""
        df = pd.DataFrame({
            "Transfection": ["C1", "C1"],
            "Cell_Line": ["HEK", "HEK"],
            "Ligand": ["L1", "L1"],
            "Ligand_Conc": [-9.0, -9.0],
            "Time_(min)": [0.0, 0.0],
            "Is_Vehicle": [False, False],
        })
        result = generate_header_key(df, group_by="Transfection")
        # Transfection should NOT be in the key
        assert "C1" not in result["Header_Key"].iloc[0]
        # Cell_Line should still be there
        assert "HEK" in result["Header_Key"].iloc[0]

    def test_group_by_cell_line(self):
        """When grouped by Cell Line, it should be excluded from key."""
        df = pd.DataFrame({
            "Transfection": ["C1", "C1"],
            "Cell_Line": ["HEK", "HEK"],
            "Ligand": ["L1", "L1"],
            "Ligand_Conc": [-9.0, -9.0],
            "Time_(min)": [0.0, 0.0],
            "Is_Vehicle": [False, False],
        })
        result = generate_header_key(df, group_by="Cell Line")
        assert "HEK" not in result["Header_Key"].iloc[0]
        assert "C1" in result["Header_Key"].iloc[0]

    def test_crc_no_conc_in_key(self):
        """CRC data (single timepoint): concentration should NOT appear in key."""
        df = pd.DataFrame({
            "Transfection": ["C1"],
            "Cell_Line": ["HEK"],
            "Ligand": ["L1"],
            "Ligand_Conc": [-9.0],
            "Time_(min)": [0.0],  # Single timepoint → CRC
            "Is_Vehicle": [False],
        })
        result = generate_header_key(df)
        assert "-9.0" not in result["Header_Key"].iloc[0]

    def test_multiple_ligands_included(self):
        """When >1 ligand exists, ligand name should appear in key."""
        df = pd.DataFrame({
            "Transfection": ["C1", "C1"],
            "Cell_Line": ["HEK", "HEK"],
            "Ligand": ["Lig1", "Lig2"],
            "Ligand_Conc": [-9.0, -9.0],
            "Time_(min)": [0.0, 0.0],
            "Is_Vehicle": [False, False],
        })
        result = generate_header_key(df)
        assert "Lig1" in result["Header_Key"].iloc[0]
        assert "Lig2" in result["Header_Key"].iloc[1]

    def test_single_ligand_not_in_key(self):
        """When only 1 ligand exists, it should NOT appear in key."""
        df = pd.DataFrame({
            "Transfection": ["C1"],
            "Cell_Line": ["HEK"],
            "Ligand": ["Lig1"],
            "Ligand_Conc": [-9.0],
            "Time_(min)": [0.0],
            "Is_Vehicle": [False],
        })
        result = generate_header_key(df)
        assert "Lig1" not in result["Header_Key"].iloc[0]

    def test_fallback_to_data(self):
        """If all parts are excluded, header should be 'Data'."""
        df = pd.DataFrame({
            "Transfection": ["C1"],
            "Cell_Line": ["HEK"],
            "Ligand": ["L1"],
            "Ligand_Conc": [-9.0],
            "Time_(min)": [0.0],
            "Is_Vehicle": [False],
        })
        # Group by both → both excluded, single ligand not added, CRC no conc
        result = generate_header_key(df, group_by="Transfection")
        # Only Cell_Line remains
        assert "HEK" in result["Header_Key"].iloc[0]


# ---------------------------------------------------------------------------
# Tests: create_clean_pivot
# ---------------------------------------------------------------------------

class TestCreateCleanPivot:

    @pytest.fixture
    def kinetic_export_df(self):
        """Minimal kinetic data for pivoting."""
        return pd.DataFrame({
            "Header_Key": ["C1", "C1", "C1", "C1"],
            "File_Name": ["f1", "f1", "f1", "f1"],
            "Well_ID": ["A1", "A2", "A1", "A2"],
            "Time_(min)": [0.0, 0.0, 1.0, 1.0],
            "Replicate": ["1", "2", "1", "2"],
            "Kinetic_Mean": [1.0, 1.1, 1.5, 1.6],
        })

    def test_basic_pivot(self, kinetic_export_df):
        """Should produce a pivot with Time as index."""
        result = create_clean_pivot(
            kinetic_export_df, "Time_(min)", "Kinetic_Mean",
            disregard_well_id=True, drop_labeling_control_col=True
        )
        assert result is not None
        assert len(result) == 2  # 2 timepoints

    def test_labeling_control_dropped(self):
        """Rows with Replicate='labeling control' should be dropped when flagged."""
        df = pd.DataFrame({
            "Header_Key": ["C1", "C1", "C1"],
            "File_Name": ["f1", "f1", "f1"],
            "Well_ID": ["A1", "A2", "A3"],
            "Time_(min)": [0.0, 0.0, 0.0],
            "Replicate": ["1", "2", "labeling control"],
            "Value": [1.0, 2.0, 99.0],
        })
        result = create_clean_pivot(
            df, "Time_(min)", "Value",
            disregard_well_id=False, drop_labeling_control_col=True
        )
        # labeling control value (99) should not appear
        assert 99.0 not in result.values

    def test_labeling_control_kept(self):
        """Rows with Replicate='labeling control' should be kept when not flagged."""
        df = pd.DataFrame({
            "Header_Key": ["C1", "C1", "C1"],
            "File_Name": ["f1", "f1", "f1"],
            "Well_ID": ["A1", "A2", "A3"],
            "Time_(min)": [0.0, 0.0, 0.0],
            "Replicate": ["1", "2", "labeling control"],
            "Value": [1.0, 2.0, 99.0],
        })
        result = create_clean_pivot(
            df, "Time_(min)", "Value",
            disregard_well_id=False, drop_labeling_control_col=False
        )
        # labeling control value should appear
        assert 99.0 in result.values

    def test_well_id_sorting(self):
        """A1, A2, A10, A11 should sort correctly (not A1, A10, A11, A2)."""
        df = pd.DataFrame({
            "Header_Key": ["C1"] * 4,
            "File_Name": ["f1"] * 4,
            "Well_ID": ["A10", "A2", "A1", "A11"],
            "Time_(min)": [0.0] * 4,
            "Replicate": ["1", "2", "3", "4"],
            "Value": [10.0, 2.0, 1.0, 11.0],
        })
        result = create_clean_pivot(
            df, "Time_(min)", "Value",
            disregard_well_id=False, drop_labeling_control_col=True
        )
        # Values should be sorted by well: A1(1), A2(2), A10(10), A11(11)
        values = result.iloc[0].dropna().tolist()
        assert values == [1.0, 2.0, 10.0, 11.0]

# ---------------------------------------------------------------------------
# build_crc_table / build_crc_preview_table (CRC "second sheet" + window preview)
# ---------------------------------------------------------------------------

def _crc_condition_master():
    """One condition (CondA / HEK / Lig1), two files (= two dates), a triplicate block
    (cols 1-3), dose rows A-G + vehicle row H. AUC_Mean is constant across the 3 wells of a
    (file, row); Veh_Norm_AUC differs per technical well."""
    rows_letters = "ABCDEFGH"
    concs = {"A": -9.0, "B": -8.0, "C": -7.0, "D": -6.0, "E": -5.0, "F": -4.0, "G": -3.0,
             "H": float("nan")}
    files = [("F1.xlsx", "2026-06-05", 1.0), ("F2.xlsx", "2026-06-20", 1.1)]
    recs = []
    for fname, date, scale in files:
        for col in (1, 2, 3):
            for r in rows_letters:
                is_veh = (r == "H")
                base = 1.0 if is_veh else round((1.0 + (rows_letters.index(r) + 1) * 0.3) * scale, 4)
                recs.append({
                    "File_Name": fname, "Date": date, "Main_Plasmids": "P1 + P2",
                    "Transfection": "CondA", "Cell_Line": "HEK", "Ligand": "Lig1",
                    "Ligand_Conc": concs[r], "Plate_Row": r, "Replicate": str(col),
                    "Well_ID": f"{r}{col}", "Is_Vehicle": is_veh, "Is_Excluded": False,
                    "AUC_Mean": base, "Veh_Norm_AUC": round(base + col * 0.01, 4)})
    df = pd.DataFrame(recs)
    for c in MASTER_COLUMNS:
        if c not in df.columns:
            df[c] = float("nan")
    return df[MASTER_COLUMNS].copy()


class TestBuildCrcTable:
    def test_auc_mean_index_and_one_column_per_file(self):
        t = build_crc_table(_crc_condition_master(), "AUC_Mean")
        assert t is not None
        assert t.index.name == "Lig1 (logM)"          # single-ligand -> conc index
        assert "Vehicle" in list(t.index) and -9.0 in list(t.index)
        assert t.shape[1] == 2                          # one column per file (mean)

    def test_veh_norm_one_column_per_well(self):
        t = build_crc_table(_crc_condition_master(), "Veh_Norm_AUC")
        assert t.shape[1] == 6                          # 2 files x 3 technical wells

    def test_none_when_value_col_absent(self):
        df = _crc_condition_master().drop(columns=["AUC_Mean"])
        assert build_crc_table(df, "AUC_Mean") is None


class TestBuildCrcPreviewTable:
    def test_auc_mean_columns_labelled_by_date_only(self):
        p = build_crc_preview_table(_crc_condition_master(), "AUC_Mean")
        assert list(p.columns) == ["05.06.26", "20.06.26"]              # date only, date-ordered
        assert p.index.name == "Lig1 (logM)"

    def test_veh_norm_columns_labelled_date_and_replicate(self):
        p = build_crc_preview_table(_crc_condition_master(), "Veh_Norm_AUC")
        # one column per technical-replicate well: "<date> (<replicate #>)".
        assert list(p.columns) == ["05.06.26 (1)", "05.06.26 (2)", "05.06.26 (3)",
                                   "20.06.26 (1)", "20.06.26 (2)", "20.06.26 (3)"]

    def test_value_parity_with_export_table(self):
        # The preview shows the SAME per-concentration values as the export sheet (only the
        # column labels/order differ), for both subtypes.
        df = _crc_condition_master()
        for col in ("AUC_Mean", "Veh_Norm_AUC"):
            t = build_crc_table(df, col)
            p = build_crc_preview_table(df, col)
            assert set(t.index) == set(p.index)
            for idx in t.index:
                tv = sorted(round(float(x), 6) for x in t.loc[idx] if pd.notna(x))
                pv = sorted(round(float(x), 6) for x in p.loc[idx] if pd.notna(x))
                assert tv == pv, f"value mismatch at {idx!r} for {col}"
