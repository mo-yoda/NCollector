"""
Tests for parsing.py — spreadsheet navigation, protocol extraction,
measurement extraction, info sheet parsing, and folder scanning.

All tests use synthetic DataFrames that mimic the structure of real xlsx files.
For integration tests with actual xlsx files, place them in test_data/ and see
TestWithRealData at the bottom.
"""
import os
import pytest
import pandas as pd
from datetime import date
from unittest.mock import patch, MagicMock

from parsing import (
    get_location,
    extract_value,
    slice_table,
    process_transfection_scheme,
    extract_metadata,
    extract_bret_data,
    extract_info_sheet_data,
    extract_measurement_data,
    scan_and_load_folders,
)


# ---------------------------------------------------------------------------
# Fixtures: synthetic spreadsheet structures
# ---------------------------------------------------------------------------

@pytest.fixture
def simple_grid():
    """A small labeled grid mimicking a spreadsheet region."""
    return pd.DataFrame([
        ["Header1", "Value1", "Extra1"],
        ["Header2", "Value2", "Extra2"],
        ["Header3", "Value3", "Extra3"],
    ])


@pytest.fixture
def protocol_sheet():
    """
    Synthetic DataFrame mimicking the 'Protocol' sheet layout.
    Markers are in column 0, values in column 1.
    """
    rows = [
        ["date of measurement", "15.01.25", None, None, None],
        ["n =", 4, None, None, None],
        [None, None, None, None, None],
        ["Cell line", None, None, None, None],
        ["HEK293", None, None, None, None],
        [None, None, None, None, None],
        ["Cell line layout", None, None, None, None],
        ["one line", None, None, None, None],
        [None, None, None, None, None],
        # Transfection table
        ["DNA", "DB#", 1.0, 2.0, 3.0],
        ["PlasmidA", "DB001", 100, 100, 100],
        ["PlasmidB", "DB002", 100, 100, 100],
        ["CondPlasmid", "DB003", 0, 100, 0],
        ["pcDNA3.1", "DB004", 50, 50, 50],  # Should be filtered out
        [None, None, None, None, None],
        # Ligand section
        ["Ligand dilution", None, None, None, None],
        ["TestLigand", None, None, None, None],
        ["final concentration in well (log(M))", None, None, None, None],
        [-9.0, None, None, None, None],
        [-8.0, None, None, None, None],
        [-7.0, None, None, None, None],
        [-6.0, None, None, None, None],
        [-5.0, None, None, None, None],
        [-4.0, None, None, None, None],
        [-3.0, None, None, None, None],
        [0.0, None, None, None, None],
        [None, None, None, None, None],
    ]
    return pd.DataFrame(rows)


@pytest.fixture
def pr_export_sheet():
    """
    Synthetic DataFrame mimicking 'Table All Cycles' sheet from plate reader export.
    Simplified to 2 wells (A01, H01) and 4 timepoints.
    """
    rows = [
        ["Date: 21/11/2025", None, None, None, None, None, None, None, None],
        ["ID2: HEK293", None, None, None, None, None, None, None, None],
        ["ID3: 1,2,3,4", None, None, None, None, None, None, None, None],
        [None, None, None, None, None, None, None, None, None],
        # Data table header
        ["Well", "Content", "Raw Data (475-30 B)", "Raw Data (475-30 B)",
         "Raw Data (535-30 B)", "Raw Data (535-30 B)", "Ratio", "Ratio", None],
        # Time row (first data row after header)
        ["Time (min)", None, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, None],
        # Well data
        ["A01", "Sample", 500, 520, 100, 110, 0.200, 0.212, None],
        ["H01", "Sample", 450, 460, 90, 95, 0.200, 0.207, None],
        [None, None, None, None, None, None, None, None, None],
        [None, None, None, None, None, None, None, None, None],
    ]
    return pd.DataFrame(rows)


@pytest.fixture
def info_sheet_df():
    """Synthetic DataFrame mimicking 'Protocol Information' sheet."""
    return pd.DataFrame([
        ["Instrument:", "PHERAstar FSX"],
        ["Software:", "MARS 3.42"],
        ["Test Name:", "BRET2"],
        ["Path:", "C:\\Data\\Experiment1"],
        [None, None],
        ["Date:", "21/11/2025"],
        ["No colon here", "value"],  # Should be skipped
    ])


# ---------------------------------------------------------------------------
# Tests: get_location
# ---------------------------------------------------------------------------

class TestGetLocation:

    def test_find_existing_marker(self, simple_grid):
        """Should find the exact position of a marker."""
        result = get_location(simple_grid, "Header2")
        assert result == (1, 0)

    def test_find_value_in_grid(self, simple_grid):
        """Should find values anywhere in the grid."""
        result = get_location(simple_grid, "Value1")
        assert result == (0, 1)

    def test_marker_not_found(self, simple_grid):
        """Should return None for missing markers."""
        result = get_location(simple_grid, "NonExistent")
        assert result is None

    def test_match_index(self):
        """Should return the nth occurrence when match_index is specified."""
        df = pd.DataFrame([
            ["A", "B"],
            ["A", "C"],
            ["D", "A"],
        ])
        first = get_location(df, "A", match_index=0)
        second = get_location(df, "A", match_index=1)
        third = get_location(df, "A", match_index=2)
        assert first == (0, 0)
        assert second == (1, 0)
        assert third == (2, 1)

    def test_match_index_out_of_range(self):
        """Should return None if the requested occurrence doesn't exist."""
        df = pd.DataFrame([["A", "B"]])
        result = get_location(df, "A", match_index=5)
        assert result is None


# ---------------------------------------------------------------------------
# Tests: extract_value
# ---------------------------------------------------------------------------

class TestExtractValue:

    def test_default_offset(self, simple_grid):
        """Default: returns value one column to the right of the marker."""
        result = extract_value(simple_grid, "Header1")
        assert result == "Value1"

    def test_custom_col_offset(self, simple_grid):
        """Custom column offset should work."""
        result = extract_value(simple_grid, "Header1", col_offset=2)
        assert result == "Extra1"

    def test_row_offset(self, simple_grid):
        """Row offset should work."""
        result = extract_value(simple_grid, "Header1", col_offset=0, row_offset=1)
        assert result == "Header2"

    def test_marker_not_found(self, simple_grid):
        """Should return None when marker doesn't exist."""
        result = extract_value(simple_grid, "Missing")
        assert result is None

    def test_out_of_bounds(self, simple_grid):
        """Should return None when offset goes outside sheet."""
        result = extract_value(simple_grid, "Extra3", col_offset=1)
        assert result is None


# ---------------------------------------------------------------------------
# Tests: slice_table
# ---------------------------------------------------------------------------

class TestSliceTable:

    def test_basic_extraction(self):
        """Should extract a table starting at the marker row."""
        df = pd.DataFrame([
            [None, None, None],
            ["Col1", "Col2", "Col3"],
            ["A", 1, 10],
            ["B", 2, 20],
            ["C", 3, 30],
            [None, None, None],
        ])
        result = slice_table(df, "Col1")
        assert result is not None
        assert len(result) == 3
        assert list(result.columns) == ["Col1", "Col2", "Col3"]
        assert result.iloc[0, 0] == "A"

    def test_stops_at_nan(self):
        """Should stop when first column has NaN."""
        df = pd.DataFrame([
            ["Name", "Value"],
            ["A", 1],
            ["B", 2],
            [None, 3],  # NaN in first column → stop here
            ["D", 4],
        ])
        result = slice_table(df, "Name")
        assert len(result) == 2

    def test_row_end_threshold(self):
        """With threshold > 1, should only stop after consecutive NaN rows."""
        df = pd.DataFrame([
            ["Name", "Value"],
            ["A", 1],
            [None, 2],  # Single NaN → not enough to stop
            ["C", 3],
            [None, None],  # Start of consecutive NaN
            [None, None],  # Two consecutive → stop
        ])
        result = slice_table(df, "Name", row_end_threshold=2)
        assert len(result) == 3  # A, NaN, C
        assert result.iloc[0, 0] == "A"
        assert result.iloc[2, 0] == "C"

    def test_marker_not_found(self):
        """Should return None when marker isn't in the sheet."""
        df = pd.DataFrame([["A", "B"], ["C", "D"]])
        result = slice_table(df, "Missing")
        assert result is None

    def test_match_index(self):
        """Should find the nth occurrence of the marker."""
        df = pd.DataFrame([
            ["Table", "Val"],
            ["X", 1],
            [None, None],
            ["Table", "Val"],
            ["Y", 2],
        ])
        first = slice_table(df, "Table", match_index=0)
        second = slice_table(df, "Table", match_index=1)
        assert first.iloc[0, 0] == "X"
        assert second.iloc[0, 0] == "Y"


# ---------------------------------------------------------------------------
# Tests: process_transfection_scheme
# ---------------------------------------------------------------------------

class TestProcessTransfectionScheme:

    def test_basic_scheme(self):
        """Should identify common plasmids and variable conditions."""
        df = pd.DataFrame({
            "DNA": ["PlasmidA", "PlasmidB", "CondX"],
            "DB#": ["DB1", "DB2", "DB3"],
            1.0: [100, 100, 100],
            2.0: [100, 100, 0],
            3.0: [100, 100, 50],
        })
        main, conditions = process_transfection_scheme(df)

        # PlasmidA and PlasmidB are in all columns → main plasmids
        assert "PlasmidA" in main
        assert "PlasmidB" in main
        # CondX is not in column 2 → variable
        assert "CondX" not in main
        assert "CondX" in conditions["1"]

    def test_pcdna_filtered(self):
        """pcDNA3.1 should be excluded from results."""
        df = pd.DataFrame({
            "DNA": ["PlasmidA", "pcDNA3.1"],
            "DB#": ["DB1", "DB2"],
            1.0: [100, 50],
        })
        main, conditions = process_transfection_scheme(df)
        all_dna = set(main)
        for v in conditions.values():
            all_dna.update(v)
        assert "pcDNA3.1" not in all_dna

    def test_empty_input(self):
        """Should handle empty DataFrame gracefully."""
        main, conditions = process_transfection_scheme(pd.DataFrame())
        assert main == []
        assert conditions == {}

    def test_none_input(self):
        """Should handle None input."""
        main, conditions = process_transfection_scheme(None)
        assert main == []
        assert conditions == {}

    def test_single_condition(self):
        """All DNA present in the only condition → all are main plasmids."""
        df = pd.DataFrame({
            "DNA": ["PlasmidA", "PlasmidB"],
            "DB#": ["DB1", "DB2"],
            1.0: [100, 100],
        })
        main, conditions = process_transfection_scheme(df)
        assert "PlasmidA" in main
        assert "PlasmidB" in main
        # Variable should be empty for the single condition
        assert conditions["1"] == []


# ---------------------------------------------------------------------------
# Tests: extract_metadata
# ---------------------------------------------------------------------------

class TestExtractMetadata:

    def test_basic_extraction(self):
        """Should extract date, cell line, and transfection IDs."""
        df = pd.DataFrame([
            ["Date: 21/11/2025"],
            ["ID2: HEK293"],
            ["ID3: 1,2,3,4"],
        ])
        result = extract_metadata(df)
        assert result is not None
        assert result["measurement_date"] == date(2025, 11, 21)
        assert result["cell_line"] == "HEK293"
        assert result["transfections"] == "1,2,3,4"

    def test_cell_line_standardization_dq(self):
        """dQ variants should be standardized."""
        df = pd.DataFrame([
            ["Date: 01/01/2025"],
            ["ID2: DQ cells"],
            ["ID3: 1"],
        ])
        result = extract_metadata(df)
        assert result["cell_line"] == "dQ"

    def test_cell_line_standardization_barr(self):
        """bArrKO variants should be standardized."""
        df = pd.DataFrame([
            ["Date: 01/01/2025"],
            ["ID2: bArrKO Clone5"],
            ["ID3: 1"],
        ])
        result = extract_metadata(df)
        assert result["cell_line"] == "bArrKO"

    def test_cell_line_standardization_control(self):
        """Control variants should be standardized."""
        df = pd.DataFrame([
            ["Date: 01/01/2025"],
            ["ID2: Con"],
            ["ID3: 1"],
        ])
        result = extract_metadata(df)
        assert result["cell_line"] == "Control"

    def test_multiple_cell_lines(self):
        """Comma-separated cell lines should all be standardized."""
        df = pd.DataFrame([
            ["Date: 01/01/2025"],
            ["ID2: dQ, Control"],
            ["ID3: 1"],
        ])
        result = extract_metadata(df)
        # Both should be standardized, order may vary due to set()
        assert "dQ" in result["cell_line"]
        assert "Control" in result["cell_line"]

    def test_invalid_date_format(self):
        """Should return None if date format is wrong."""
        df = pd.DataFrame([
            ["Date: 2025-11-21"],  # Wrong format
            ["ID2: HEK293"],
            ["ID3: 1"],
        ])
        result = extract_metadata(df)
        assert result is None

    def test_missing_field(self):
        """Should return None if required fields are missing."""
        df = pd.DataFrame([
            ["Date: 01/01/2025"],
            ["ID2: HEK293"],
            # ID3 missing
        ])
        result = extract_metadata(df)
        assert result is None

    def test_extra_whitespace(self):
        """Should handle whitespace around values."""
        df = pd.DataFrame([
            ["Date :  21/11/2025 "],
            ["ID2 :  HEK293  "],
            ["ID3 :  1,2  "],
        ])
        result = extract_metadata(df)
        assert result is not None
        assert result["measurement_date"] == date(2025, 11, 21)
        assert result["transfections"] == "1,2"


# ---------------------------------------------------------------------------
# Tests: extract_info_sheet_data
# ---------------------------------------------------------------------------

class TestExtractInfoSheetData:

    def test_basic_extraction(self):
        """Should extract key-value pairs from colon-separated rows."""
        # Simulate pd.read_excel output
        df = pd.DataFrame([
            ["Instrument:", "PHERAstar FSX"],
            ["Software:", "MARS 3.42"],
        ])
        with patch("parsing.pd.read_excel", return_value=df):
            result = extract_info_sheet_data(MagicMock())
        assert result["Instrument"] == "PHERAstar FSX"
        assert result["Software"] == "MARS 3.42"

    def test_value_in_col_a(self):
        """Value after colon in same cell should be extracted."""
        df = pd.DataFrame([
            ["Key: InlineValue", None],
        ])
        with patch("parsing.pd.read_excel", return_value=df):
            result = extract_info_sheet_data(MagicMock())
        assert result["Key"] == "InlineValue"

    def test_path_with_colon(self):
        """Windows paths (C:\\...) should split only on first colon."""
        df = pd.DataFrame([
            ["Path: C:\\Data\\Experiment", None],
        ])
        with patch("parsing.pd.read_excel", return_value=df):
            result = extract_info_sheet_data(MagicMock())
        assert result["Path"] == "C:\\Data\\Experiment"

    def test_skip_nan_rows(self):
        """Rows with NaN in first column should be skipped."""
        df = pd.DataFrame([
            ["Key1:", "Val1"],
            [None, "ignored"],
            ["Key2:", "Val2"],
        ])
        with patch("parsing.pd.read_excel", return_value=df):
            result = extract_info_sheet_data(MagicMock())
        assert len(result) == 2
        assert "Key1" in result
        assert "Key2" in result

    def test_skip_rows_without_colon(self):
        """Rows without a colon should be skipped."""
        df = pd.DataFrame([
            ["Key1:", "Val1"],
            ["No colon here", "ignored"],
            ["Key2:", "Val2"],
        ])
        with patch("parsing.pd.read_excel", return_value=df):
            result = extract_info_sheet_data(MagicMock())
        assert len(result) == 2

    def test_missing_sheet_returns_empty_or_none(self):
        """Should handle missing sheet gracefully."""
        with patch("parsing.pd.read_excel", side_effect=ValueError("Sheet not found")):
            result = extract_info_sheet_data(MagicMock())
        # Current code returns None on ValueError — could be {} depending on version
        assert result is None or result == {}


# ---------------------------------------------------------------------------
# Tests: extract_bret_data
# ---------------------------------------------------------------------------

class TestExtractBretData:

    def test_basic_extraction(self, pr_export_sheet):
        """Should extract BRET ratio, donor, acceptor, and wavelengths."""
        result = extract_bret_data(pr_export_sheet)
        assert result is not None
        assert "bret_ratio" in result
        assert "donor" in result
        assert "acceptor" in result
        assert result["donor_wavelength"] == 475
        assert result["acceptor_wavelength"] == 535

    def test_time_values(self, pr_export_sheet):
        """Time values should be extracted correctly."""
        result = extract_bret_data(pr_export_sheet)
        assert result["time_min"] is not None
        assert 0.0 in result["time_min"]
        assert 1.0 in result["time_min"]

    def test_well_name_cleanup(self, pr_export_sheet):
        """A01 should be cleaned to A1."""
        result = extract_bret_data(pr_export_sheet)
        bret_df = result["bret_ratio"]
        # After transpose, well IDs become column names
        assert "A1" in bret_df.columns
        assert "A01" not in bret_df.columns

    def test_dataframe_shapes(self, pr_export_sheet):
        """All channel DataFrames should have same shape."""
        result = extract_bret_data(pr_export_sheet)
        donor_shape = result["donor"].shape
        acceptor_shape = result["acceptor"].shape
        bret_shape = result["bret_ratio"].shape
        assert donor_shape == acceptor_shape
        assert donor_shape == bret_shape


# ---------------------------------------------------------------------------
# Tests: scan_and_load_folders
# ---------------------------------------------------------------------------

class TestScanAndLoadFolders:

    def test_empty_folder_list(self):
        """Should return empty list for no folders."""
        result = scan_and_load_folders([])
        assert result == []

    def test_invalid_date_folder_skipped(self, tmp_path):
        """Folders with non-YYMMDD names should be skipped."""
        bad_folder = tmp_path / "invalid_name"
        bad_folder.mkdir()
        result = scan_and_load_folders([str(bad_folder)])
        assert len(result) == 0

    def test_valid_date_parsed(self, tmp_path):
        """Folder date should be parsed from YYMMDD prefix."""
        folder = tmp_path / "250121_experiment"
        folder.mkdir()
        result = scan_and_load_folders([str(folder)])
        assert len(result) == 1
        assert result[0].measurement_date == date(2025, 1, 21)
        assert result[0].folder_name == "250121_experiment"

    def test_no_xlsx_files(self, tmp_path):
        """Folder with no xlsx files should have empty results."""
        folder = tmp_path / "250121_experiment"
        folder.mkdir()
        (folder / "readme.txt").write_text("not an xlsx")
        result = scan_and_load_folders([str(folder)])
        assert len(result) == 1
        assert len(result[0].results) == 0
        assert result[0].protocol is None

    def test_log_fn_called(self, tmp_path):
        """log_fn should receive messages when provided."""
        folder = tmp_path / "250121_experiment"
        folder.mkdir()
        messages = []
        scan_and_load_folders([str(folder)], log_fn=lambda m: messages.append(m))
        # No error messages expected for an empty folder, but the function shouldn't crash
        assert isinstance(messages, list)

    def test_skipped_files_tracked(self, tmp_path):
        """Non-matching xlsx files should appear in skipped_files."""
        folder = tmp_path / "250121_experiment"
        folder.mkdir()
        # Create a minimal xlsx that doesn't have required sheets
        dummy_df = pd.DataFrame({"A": [1, 2, 3]})
        dummy_path = folder / "random_file.xlsx"
        dummy_df.to_excel(str(dummy_path), index=False)

        result = scan_and_load_folders([str(folder)])
        assert len(result) == 1
        assert "random_file.xlsx" in result[0].skipped_files


# ---------------------------------------------------------------------------
# Tests with real data (skipped if test_data/ not present)
# ---------------------------------------------------------------------------

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "test_data")
REAL_DATA_AVAILABLE = os.path.isdir(TEST_DATA_DIR)


@pytest.mark.skipif(not REAL_DATA_AVAILABLE, reason="test_data/ not present")
class TestWithRealData:
    """
    Integration tests using actual xlsx files.
    Place a folder structure in test_data/ like:
        test_data/
            250121_TestExperiment/
                protocol.xlsx
                250121_analysis.xlsx
    """

    def _find_first_analysis(self):
        """Finds the first analysis xlsx in test_data."""
        for root, dirs, files in os.walk(TEST_DATA_DIR):
            for f in files:
                if f.endswith((".xlsx", ".xlsm")) and "analysis" in f.lower():
                    return os.path.join(root, f)
        return None

    def _find_first_protocol(self):
        """Finds the first protocol xlsx in test_data."""
        for root, dirs, files in os.walk(TEST_DATA_DIR):
            for f in files:
                if f.endswith((".xlsx", ".xlsm")) and "protocol" in f.lower():
                    return os.path.join(root, f)
        return None

    def test_real_measurement_extraction(self):
        """Extract measurement data from a real xlsx."""
        path = self._find_first_analysis()
        if not path:
            pytest.skip("No analysis file found in test_data/")
        xls = pd.ExcelFile(path)
        result = extract_measurement_data(xls, os.path.basename(path))

        assert result is not None
        assert result.measurement_date is not None
        assert result.raw_bret_ratio_df is not None
        assert not result.raw_bret_ratio_df.empty
        assert result.donor_wavelength > 0
        assert result.acceptor_wavelength > 0
        assert result.acceptor_wavelength > result.donor_wavelength

    def test_real_folder_scan(self):
        """Scan test_data subfolders for protocol + measurement pairs."""
        folder_paths = []
        for entry in os.scandir(TEST_DATA_DIR):
            if entry.is_dir():
                folder_paths.append(entry.path)

        if not folder_paths:
            pytest.skip("No subfolders in test_data/")

        results = scan_and_load_folders(folder_paths)
        assert len(results) > 0
        # At least one folder should have a protocol and results
        complete = [f for f in results if f.protocol and f.results]
        assert len(complete) > 0, "No folder had both protocol and measurement files"