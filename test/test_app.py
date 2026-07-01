"""
Tests for app.py — compile_master_dataframe, write_excel_export,
and enrichment validation logic.

Most processing/parsing/export logic is tested in their own test files.
These tests cover the app-level orchestration and data flow.

Uses a MockApp that provides the same attributes as NCollectorApp
without requiring a tkinter display.
"""
import os
import pytest
import numpy as np
import pandas as pd
from datetime import date
from unittest.mock import MagicMock, patch, PropertyMock

from models import (
    PlateColMetadata, PrResult, ProtocolData, MeasurementFolder,
    ProcessingConfig, MASTER_COLUMNS, APP_VERSION, TRIPLICATE_LAYOUT,
)

# Import just the methods we need — we'll bind them to a mock
from app import NCollectorApp


class MockApp:
    """Lightweight stand-in for NCollectorApp without tkinter dependency."""

    def __init__(self):
        self.directory = ""
        self.experiment = []
        self.master_df = pd.DataFrame()
        self.master_index = pd.DataFrame()
        self.rule_history_text = ""
        self.lbl_data_source = MagicMock()       # Mock tkinter label (Tab 3)
        self.lbl_data_source_tab1 = None         # Tab 1 mirror label (skipped when None)

    def log(self, message):
        pass  # Suppress log output during tests

    # Stubs for GUI methods called by the real methods
    def refresh_filter_options(self): pass
    def update_summary_table(self): pass
    def refresh_plot_helper_options(self): pass

    # Bind the real methods from NCollectorApp
    compile_master_dataframe = NCollectorApp.compile_master_dataframe
    write_excel_export = NCollectorApp.write_excel_export
    built_master_index = NCollectorApp.built_master_index
    _build_index_from_objects = NCollectorApp._build_index_from_objects
    _build_index_from_master = NCollectorApp._build_index_from_master
    _finalize_index = NCollectorApp._finalize_index
    export_master_csv = NCollectorApp.export_master_csv
    _set_data_source = NCollectorApp._set_data_source


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def app():
    """Creates a MockApp with the same interface as NCollectorApp."""
    return MockApp()


@pytest.fixture
def processed_experiment():
    """
    Creates a minimal experiment structure (1 folder, 1 protocol, 1 result)
    with all processing results populated, ready for compile_master_dataframe.
    """
    n_tp = 5
    baseline = 2
    t_vec = [-2.0, -1.0, 0.0, 1.0, 2.0]
    wells = [f"{r}{c}" for r in "ABCDEFGH" for c in range(1, 13)]

    # Synthetic BRET data (timepoints × wells)
    rng = np.random.RandomState(42)
    raw_data = {w: rng.uniform(0.18, 0.22, n_tp) for w in wells}
    raw_df = pd.DataFrame(raw_data)

    # Build kinetic_df, bl_corr, etc. as simple copies with slight modifications
    kinetic_df = pd.DataFrame({w: np.ones(n_tp) + rng.normal(0, 0.01, n_tp) for w in wells})
    bl_corr_df = pd.DataFrame({w: np.ones(n_tp) + rng.normal(0, 0.01, n_tp) for w in wells})
    lab_corr_df = pd.DataFrame({w: [float("nan")] * n_tp for w in wells})

    # Mean DFs: "Condition|Cell|Ligand|Row"
    mean_cols = {}
    for block_idx, block in enumerate(TRIPLICATE_LAYOUT):
        cond = f"Cond{block_idx + 1}"
        for row in "ABCDEFGH":
            key = f"{cond}|HEK293|TestLigand|{row}"
            mean_cols[key] = np.ones(n_tp) + rng.normal(0, 0.005, n_tp)
    kinetic_mean_df = pd.DataFrame(mean_cols)

    # AUC (single row)
    auc_df = pd.DataFrame({w: [rng.uniform(6, 8)] for w in wells})
    auc_mean_df = pd.DataFrame({k: [v.mean()] for k, v in mean_cols.items()})

    # LP (single row)
    lp_df = pd.DataFrame({w: [rng.uniform(0.9, 1.1)] for w in wells})
    lp_mean_df = pd.DataFrame({k: [v.mean()] for k, v in mean_cols.items()})

    # Build metadata
    col_metadata = {}
    conditions = ["Cond1", "Cond2", "Cond3", "Cond4"]
    for i, block in enumerate(TRIPLICATE_LAYOUT):
        for rep_idx, col in enumerate(block):
            col_metadata[col] = PlateColMetadata(
                cell_line="HEK293",
                transfection_id=str(i + 1),
                condition_name=conditions[i],
                plasmids=[f"Plasmid_{conditions[i]}"],
                ligand_identity="TestLigand",
                ligand_conc={r: float(-9 + ord(r) - ord("A")) for r in "ABCDEFG"} | {"H": float("nan")},
                replicate=str(rep_idx + 1),
            )

    result = PrResult(
        file_name="test_analysis.xlsx",
        measurement_date=date(2025, 1, 21),
        cell_line="HEK293",
        transfection_id="1,2,3,4",
        raw_time=[0.0, 1.0, 2.0, 4.0, 5.0],
        raw_bret_ratio_df=raw_df.copy(),
        donor_df=raw_df.copy(),
        acceptor_df=raw_df.copy(),
        donor_wavelength=475,
        acceptor_wavelength=535,
        column_metadata=col_metadata,
        time_vector=t_vec,
        raw_bret_ratio_cleaned=raw_df.copy(),
        labeling_corr_kinetic=lab_corr_df,
        bl_corr_kinetic=bl_corr_df,
        kinetic_df=kinetic_df,
        kinetic_mean_df=kinetic_mean_df,
        raw_bret_points_df=raw_df.iloc[-3:].mean().to_frame().T,
        labeling_corr_lp_df=lab_corr_df.iloc[-3:].mean().to_frame().T,
        bl_corr_lp_df=bl_corr_df.iloc[-3:].mean().to_frame().T,
        lp_df=lp_df,
        lp_mean_df=lp_mean_df,
        labeling_corr_auc_df=pd.DataFrame({w: [float("nan")] for w in wells}),
        bl_corr_auc_df=pd.DataFrame({w: [rng.uniform(5, 8)] for w in wells}),
        auc_df=auc_df,
        auc_mean_df=auc_mean_df,
    )

    protocol = ProtocolData(
        file_name="test_protocol.xlsx",
        exp_date=date(2025, 1, 21),
        n=1,
        cell_lines=["HEK293"],
        line_layout="one line",
        transfection_scheme=pd.DataFrame(),
        main_plasmids=["PlasmidA", "PlasmidB"],
        transfection_conditions={
            "1": ["Cond1"], "2": ["Cond2"],
            "3": ["Cond3"], "4": ["Cond4"],
        },
        ligand="TestLigand",
        ligand_conc=pd.DataFrame({"conc": [-9, -8, -7, -6, -5, -4, -3, 0]}),
    )

    folder = MeasurementFolder(
        folder_name="250121_test",
        folder_path="/data/250121_test",
        measurement_date=date(2025, 1, 21),
        protocol=protocol,
        results=[result],
    )

    return [folder]


# ---------------------------------------------------------------------------
# Tests: compile_master_dataframe
# ---------------------------------------------------------------------------

class TestCompileMasterDataframe:

    def test_produces_dataframe(self, app, processed_experiment):
        """Should return a non-empty DataFrame."""
        app.experiment = processed_experiment
        app.directory = "/data"
        result = app.compile_master_dataframe()
        assert result is not None
        assert not result.empty

    def test_columns_from_master_schema(self, app, processed_experiment):
        """All columns should be from MASTER_COLUMNS (no extras)."""
        app.experiment = processed_experiment
        app.directory = "/data"
        result = app.compile_master_dataframe()
        for col in result.columns:
            assert col in MASTER_COLUMNS, f"Unexpected column '{col}' not in MASTER_COLUMNS"

    def test_essential_columns_present(self, app, processed_experiment):
        """Key columns should always be populated."""
        app.experiment = processed_experiment
        app.directory = "/data"
        result = app.compile_master_dataframe()
        for col in ["File_Name", "Date", "Well_ID", "Time_(min)",
                     "Transfection", "Cell_Line", "Ligand", "Plate_Row"]:
            assert col in result.columns
            assert result[col].notna().any(), f"'{col}' is all NaN"

    def test_row_count(self, app, processed_experiment):
        """Should have timepoints × wells rows per file."""
        app.experiment = processed_experiment
        app.directory = "/data"
        result = app.compile_master_dataframe()
        n_timepoints = 5
        n_wells = 96  # 8 rows × 12 cols
        expected = n_timepoints * n_wells
        assert len(result) == expected

    def test_metadata_populated(self, app, processed_experiment):
        """Metadata columns should be filled from protocol mapping."""
        app.experiment = processed_experiment
        app.directory = "/data"
        result = app.compile_master_dataframe()
        assert result["NCollector_version"].iloc[0] == APP_VERSION
        assert result["Main_Plasmids"].iloc[0] == "PlasmidA + PlasmidB"
        assert "Cond1" in result["Transfection"].values

    def test_vehicle_rows_flagged(self, app, processed_experiment):
        """Row H with NaN concentration should have Is_Vehicle=True."""
        app.experiment = processed_experiment
        app.directory = "/data"
        result = app.compile_master_dataframe()
        h_rows = result[result["Plate_Row"] == "H"]
        assert h_rows["Is_Vehicle"].all(), "Not all H rows are flagged as vehicle"

    def test_non_vehicle_rows_not_flagged(self, app, processed_experiment):
        """Rows A-G should have Is_Vehicle=False."""
        app.experiment = processed_experiment
        app.directory = "/data"
        result = app.compile_master_dataframe()
        non_h = result[result["Plate_Row"] != "H"]
        assert not non_h["Is_Vehicle"].any(), "Non-H rows incorrectly flagged as vehicle"

    def test_excluded_wells_marked(self, app, processed_experiment):
        """Excluded wells should have Is_Excluded=True."""
        processed_experiment[0].results[0].excluded_wells = ["A1", "B1"]
        app.experiment = processed_experiment
        app.directory = "/data"
        result = app.compile_master_dataframe()
        excluded = result[result["Well_ID"].isin(["A1", "B1"])]
        assert excluded["Is_Excluded"].all()

    def test_excluded_file_skipped(self, app, processed_experiment):
        """Files with is_excluded=True should not appear in master."""
        processed_experiment[0].results[0].is_excluded = True
        app.experiment = processed_experiment
        app.directory = "/data"
        result = app.compile_master_dataframe()
        assert result.empty or len(result) == 0

    def test_empty_experiment_returns_empty(self, app):
        """No experiment data should return empty DataFrame."""
        app.experiment = []
        result = app.compile_master_dataframe()
        assert result is None or (hasattr(result, 'empty') and result.empty)

    def test_pr_time_mapped(self, app, processed_experiment):
        """PR_Time(min) should be populated from raw_time."""
        app.experiment = processed_experiment
        app.directory = "/data"
        result = app.compile_master_dataframe()
        assert "PR_Time(min)" in result.columns
        assert result["PR_Time(min)"].notna().any()

    def test_info_sheet_populated(self, app, processed_experiment):
        """Info_Sheet should be populated."""
        processed_experiment[0].results[0].info_sheet = {"Key": "Value"}
        app.experiment = processed_experiment
        app.directory = "/data"
        result = app.compile_master_dataframe()
        assert "Info_Sheet" in result.columns
        assert result["Info_Sheet"].iloc[0] != ""


# ---------------------------------------------------------------------------
# Tests: write_excel_export
# ---------------------------------------------------------------------------

class TestWriteExcelExport:

    @pytest.fixture
    def export_df(self):
        """Minimal master DataFrame for export testing."""
        rows = []
        for t in [0.0, 1.0, 2.0]:
            for well in ["A1", "A2", "A3", "H1", "H2", "H3"]:
                row_char = well[0]
                rows.append({
                    "File_Name": "f1",
                    "Transfection": "CondA",
                    "Cell_Line": "HEK",
                    "Ligand": "Lig1",
                    "Ligand_Conc": -9.0 if row_char == "A" else float("nan"),
                    "Plate_Row": row_char,
                    "Well_ID": well,
                    "Time_(min)": t,
                    "Replicate": "1",
                    "Is_Vehicle": row_char == "H",
                    "Is_Excluded": False,
                    "Kinetic_Mean": 1.5 if row_char == "A" else 1.0,
                    "AUC_Mean": 1.3 if row_char == "A" else 1.0,
                    "Veh_Norm_AUC": 1.3 if row_char == "A" else 1.0,
                    "Raw_BRET_kinetic": 0.2,
                })
        return pd.DataFrame(rows)

    def test_creates_excel_file(self, app, export_df, tmp_path):
        """Should create an xlsx file."""
        app.master_df = export_df
        path = str(tmp_path / "test_export.xlsx")

        config = {
            "cells": "All",
            "transfections": "All",
            "ligands": "All",
            "data_types": ["AUC_Mean"],
            "group_by": "None",
            "conc_select": [{"row": "A", "ligand": "Lig1", "is_vehicle": False, "conc": "-9.0"}],
        }
        app.write_excel_export(path, export_df, config)
        assert os.path.exists(path)

    def test_metadata_sheet_created(self, app, export_df, tmp_path):
        """Exported Excel should contain a Metadata sheet."""
        app.master_df = export_df
        path = str(tmp_path / "test_meta.xlsx")

        config = {
            "cells": "All",
            "transfections": "All",
            "ligands": "All",
            "data_types": ["AUC_Mean"],
            "group_by": "None",
            "conc_select": [{"row": "A", "ligand": "Lig1", "is_vehicle": False, "conc": "-9.0"}],
        }
        app.write_excel_export(path, export_df, config)

        xls = pd.ExcelFile(path)
        assert "Metadata" in xls.sheet_names

    def test_kinetic_export(self, app, export_df, tmp_path):
        """Kinetic data type should produce a sheet."""
        app.master_df = export_df
        path = str(tmp_path / "test_kinetic.xlsx")

        config = {
            "cells": "All",
            "transfections": "All",
            "ligands": "All",
            "data_types": ["Kinetic_Mean"],
            "group_by": "None",
            "conc_mode": [
                {"row": "A", "ligand": "Lig1", "is_vehicle": False, "conc": "-9.0"},
                {"row": "H", "ligand": "Lig1", "is_vehicle": True},
            ],
        }
        app.write_excel_export(path, export_df, config)

        xls = pd.ExcelFile(path)
        # Should have Metadata + at least one data sheet
        assert len(xls.sheet_names) >= 2

    def test_crc_export(self, app, export_df, tmp_path):
        """CRC data type should produce a sheet."""
        app.master_df = export_df
        path = str(tmp_path / "test_crc.xlsx")

        config = {
            "cells": "All",
            "transfections": "All",
            "ligands": "All",
            "data_types": ["AUC_Mean"],
            "group_by": "None",
            "conc_select": [],
        }
        app.write_excel_export(path, export_df, config)

        xls = pd.ExcelFile(path)
        assert len(xls.sheet_names) >= 2

    def test_group_by_produces_multiple_sheets(self, app, tmp_path):
        """Grouping by Transfection should create one data sheet per group."""
        rows = []
        for cond in ["CondA", "CondB"]:
            for well in ["A1", "H1"]:
                row_char = well[0]
                rows.append({
                    "File_Name": "f1",
                    "Transfection": cond,
                    "Cell_Line": "HEK",
                    "Ligand": "Lig1",
                    "Ligand_Conc": -9.0 if row_char == "A" else float("nan"),
                    "Plate_Row": row_char,
                    "Well_ID": well,
                    "Time_(min)": 0.0,
                    "Replicate": "1",
                    "Is_Vehicle": row_char == "H",
                    "Is_Excluded": False,
                    "AUC_Mean": 1.3,
                    "Veh_Norm_AUC": 1.3,
                })
        df = pd.DataFrame(rows)
        app.master_df = df
        path = str(tmp_path / "test_grouped.xlsx")

        config = {
            "cells": "All",
            "transfections": "All",
            "ligands": "All",
            "data_types": ["AUC_Mean"],
            "group_by": "Transfection",
            "conc_select": [],
        }
        app.write_excel_export(path, df, config)

        xls = pd.ExcelFile(path)
        # Metadata + CondA sheet + CondB sheet
        assert len(xls.sheet_names) >= 3

    def test_empty_filter_produces_fallback(self, app, export_df, tmp_path):
        """Filtering to nothing should not crash."""
        app.master_df = export_df
        path = str(tmp_path / "test_empty.xlsx")

        config = {
            "cells": ["NonExistent"],
            "transfections": "All",
            "ligands": "All",
            "data_types": ["AUC_Mean"],
            "group_by": "None",
            "conc_select": [],
        }
        # Should not raise
        app.write_excel_export(path, export_df, config)


# ---------------------------------------------------------------------------
# Tests: built_master_index
# ---------------------------------------------------------------------------

class TestBuiltMasterIndex:

    def test_produces_index(self, app, processed_experiment):
        """Should populate master_index from experiment data."""
        app.experiment = processed_experiment
        result = app.built_master_index()
        assert not result.empty

    def test_index_columns(self, app, processed_experiment):
        """Index should contain expected columns."""
        app.experiment = processed_experiment
        result = app.built_master_index()
        for col in ["File_Name", "Date", "Cell_Line", "Condition", "Ligand", "Replicate"]:
            assert col in result.columns

    def test_excludes_empty_conditions(self, app, processed_experiment):
        """Columns with condition_name='Empty' should not appear."""
        # Set one column to Empty
        processed_experiment[0].results[0].column_metadata[1].condition_name = "Empty"
        app.experiment = processed_experiment
        result = app.built_master_index()
        assert "Empty" not in result["Condition"].values

    def test_excluded_files_skipped(self, app, processed_experiment):
        """Files with is_excluded=True should not appear in index."""
        processed_experiment[0].results[0].is_excluded = True
        app.experiment = processed_experiment
        result = app.built_master_index()
        assert result.empty

    def test_empty_experiment(self, app):
        """Empty experiment should produce empty index."""
        app.experiment = []
        result = app.built_master_index()
        assert result.empty


# ---------------------------------------------------------------------------
# Tests: export_master_csv
# ---------------------------------------------------------------------------

class TestExportMasterCsv:

    def test_saves_csv(self, app, processed_experiment, tmp_path):
        """Should save master_df to CSV."""
        app.experiment = processed_experiment
        app.directory = "/data"
        app.master_df = app.compile_master_dataframe()

        path = str(tmp_path / "master.csv")
        with patch("app.filedialog.asksaveasfilename", return_value=path):
            app.export_master_csv()

        assert os.path.exists(path)
        df = pd.read_csv(path)
        assert not df.empty

    def test_updates_label(self, app, processed_experiment, tmp_path):
        """After save, lbl_data_source should be updated with the new filename."""
        app.experiment = processed_experiment
        app.directory = "/data"
        app.master_df = app.compile_master_dataframe()

        path = str(tmp_path / "my_master.csv")
        with patch("app.filedialog.asksaveasfilename", return_value=path):
            app.export_master_csv()

        # MagicMock: check that config was called with text containing the filename
        app.lbl_data_source.config.assert_called()
        call_kwargs = app.lbl_data_source.config.call_args[1]
        assert "my_master.csv" in call_kwargs["text"]

    def test_cancel_saves_nothing(self, app, processed_experiment):
        """Cancelling the dialog should not crash or save."""
        app.experiment = processed_experiment
        app.directory = "/data"
        app.master_df = app.compile_master_dataframe()

        with patch("app.filedialog.asksaveasfilename", return_value=""):
            app.export_master_csv()  # Should not raise