"""
Tests for the 2.0.6 parsing additions — ID3 range parsing and checks, optional protocol
markers, the plate-layout field, concentration lists and the folder-scan order.

Covers: parse_id_list, check_id3, get_location / extract_value (required=False),
_read_plate_mode, extract_protocol_info (plate_mode, plate_mode_inferred, optional
'Ligand layout'), slice_table (concentration lists), extract_measurement_data
(PlateFileSkip), scan_and_load_folders, and the ID3 mapping in
processing.map_plate_metadata.

Unit tests use synthetic data. The end-to-end tests read the local fixtures in
test/data/layout_2.0.6/ (git-ignored; found by position, not by name; are skipped
when that folder is absent).
"""
import logging
import os
import shutil

import pandas as pd
import pytest

from models import DEFAULT_PLATE_MODE, ProcessingConfig, TRIPLICATE_LAYOUT
from mapping import built_conc_dic
from parsing import (
    check_id3,
    extract_protocol_info,
    extract_value,
    get_location,
    parse_id_list,
    scan_and_load_folders,
    slice_table,
)
from processing import map_plate_metadata

FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "layout_2.0.6")
HAS_FIXTURES = os.path.isdir(FIXTURES)
needs_fixtures = pytest.mark.skipif(not HAS_FIXTURES, reason="layout_2.0.6 fixtures not present")


def _experiment_folder(category: str) -> str | None:
    """The (single) YYMMDD_… experiment folder of a fixture category, found by position."""
    base = os.path.join(FIXTURES, category)
    if not os.path.isdir(base):
        return None
    folders = sorted(d for d in os.listdir(base) if os.path.isdir(os.path.join(base, d)))
    return os.path.join(base, folders[0]) if folders else None


# Experiment folders (real data, local only): old triplicate layout (one ligand,
# transfections 1-8), one-concentration layout, second ligand, labeling control.
OLD_FOLDER = _experiment_folder("old_triplicate")
NEW_FOLDER = _experiment_folder("new_one_concentration")
SECOND_LIGAND_FOLDER = _experiment_folder("second_ligand")
LABELING_FOLDER = _experiment_folder("old_labeling")
BROKEN = os.path.join(FIXTURES, "broken")

PROTOCOL_IDS = [str(i) for i in range(1, 9)]  # old triplicate protocol: transfections 1-8


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _scan(folder, logs=None):
    logs = [] if logs is None else logs
    return scan_and_load_folders([str(folder)], log_fn=logs.append), logs


def _folder_with(tmp_path, name, files):
    """tmp experiment folder `name` containing copies of `files` [(src, target_name)]."""
    folder = tmp_path / name
    folder.mkdir()
    for src, target in files:
        shutil.copy(src, folder / target)
    return folder


def _protocol_name(folder) -> str:
    """File name of the workbook with a 'Protocol' sheet in an experiment folder."""
    for name in sorted(os.listdir(folder)):
        if name.endswith((".xlsx", ".xlsm")):
            with pd.ExcelFile(os.path.join(folder, name)) as xls:
                if "Protocol" in xls.sheet_names:
                    return name
    raise FileNotFoundError(f"no protocol workbook in {folder}")


def _load_protocol(folder):
    name = _protocol_name(folder)
    with pd.ExcelFile(os.path.join(folder, name)) as xls:
        return extract_protocol_info(xls, name)


def _old_folder_with_plate(tmp_path, plate_src, plate_name=None):
    """Copy of the old triplicate protocol + one plate file in a folder with the same date."""
    protocol = _protocol_name(OLD_FOLDER)
    date_prefix = os.path.basename(OLD_FOLDER).split("_")[0]  # YYMMDD must match the protocol
    return _folder_with(tmp_path, f"{date_prefix}_step1", [
        (os.path.join(OLD_FOLDER, protocol), protocol),
        (plate_src, plate_name or os.path.basename(plate_src)),
    ])


def _block_conditions(result):
    """Condition name of the first column of every triplicate block."""
    return [result.column_metadata[b[0]].condition_name for b in TRIPLICATE_LAYOUT]


# ---------------------------------------------------------------------------
# parse_id_list
# ---------------------------------------------------------------------------

class TestParseIdList:

    def test_full_range(self):
        assert parse_id_list("1-12") == [str(i) for i in range(1, 13)]

    def test_range_and_single(self):
        assert parse_id_list("1-4,7") == ["1", "2", "3", "4", "7"]

    def test_range_with_spaces(self):
        assert parse_id_list(" 1 - 3 , 5 ") == ["1", "2", "3", "5"]

    def test_single_value_range(self):
        assert parse_id_list("2-2") == ["2"]

    def test_plain_list_unchanged(self):
        """Identical to the 2.0.5 split for plain lists (no behaviour change)."""
        assert parse_id_list("1,2,3,4") == ["1", "2", "3", "4"]

    def test_descending_range_raises(self):
        with pytest.raises(ValueError, match="descending range '4-1'"):
            parse_id_list("4-1")

    def test_descending_range_inside_list_raises(self):
        with pytest.raises(ValueError):
            parse_id_list("1,2,6-5")

    def test_empty_entries_keep_position(self):
        assert parse_id_list("1,,3,4") == ["1", "", "3", "4"]

    def test_trailing_comma_keeps_position(self):
        assert parse_id_list("1,2,3,") == ["1", "2", "3", ""]

    def test_non_numeric_returned_for_check(self):
        assert parse_id_list("1a,2") == ["1a", "2"]

    @pytest.mark.parametrize("raw", [None, "", "   ", float("nan")])
    def test_empty_input(self, raw):
        assert parse_id_list(raw) == []


# ---------------------------------------------------------------------------
# check_id3
# ---------------------------------------------------------------------------

class TestCheckId3:

    def test_valid_list_no_warning(self):
        ids, warnings = check_id3(["1", "2", "3", "4"], PROTOCOL_IDS, 4)
        assert ids == ["1", "2", "3", "4"]
        assert warnings == []

    def test_non_numeric_dropped_position_kept(self):
        ids, warnings = check_id3(["1a", "2"], PROTOCOL_IDS, 4)
        assert ids == ["", "2"]
        assert len(warnings) == 1
        assert "'1a'" in warnings[0] and "not a number" in warnings[0]

    def test_id_not_in_protocol_dropped_position_kept(self):
        ids, warnings = check_id3(["1", "13", "3"], PROTOCOL_IDS, 4)
        assert ids == ["1", "", "3"]
        assert len(warnings) == 1
        assert "'13'" in warnings[0] and "not a transfection of the protocol" in warnings[0]

    def test_more_ids_than_conditions_extras_ignored(self):
        ids, warnings = check_id3([str(i) for i in range(1, 13)], [str(i) for i in range(1, 13)], 4)
        assert ids == ["1", "2", "3", "4"]
        assert len(warnings) == 1
        assert "12 IDs for 4 conditions" in warnings[0]
        assert "5, 6, 7, 8, 9, 10, 11, 12" in warnings[0]

    def test_extras_are_not_checked_against_protocol(self):
        """Ignored extras produce only the 'extra IDs' warning, no 'not in protocol' warnings."""
        ids, warnings = check_id3([str(i) for i in range(1, 13)], PROTOCOL_IDS, 4)
        assert ids == ["1", "2", "3", "4"]
        assert len(warnings) == 1

    def test_empty_entry_no_warning(self):
        ids, warnings = check_id3(["1", "", "3", "4"], PROTOCOL_IDS, 4)
        assert ids == ["1", "", "3", "4"]
        assert warnings == []

    def test_trailing_empty_extra_no_warning(self):
        """'1,2,3,4,' on a 4-condition plate: the empty 5th entry is ignored silently."""
        ids, warnings = check_id3(["1", "2", "3", "4", ""], PROTOCOL_IDS, 4)
        assert ids == ["1", "2", "3", "4"]
        assert warnings == []

    def test_fewer_ids_than_conditions_no_warning(self):
        ids, warnings = check_id3(["1", "2"], PROTOCOL_IDS, 4)
        assert ids == ["1", "2"]
        assert warnings == []

    def test_no_ids(self):
        assert check_id3([], PROTOCOL_IDS, 4) == ([], [])

    def test_n_conditions_is_a_parameter(self):
        """Labeling layout today = 3 blocks: the 4th ID is an extra."""
        ids, warnings = check_id3(["1", "2", "3", "4"], PROTOCOL_IDS, 3)
        assert ids == ["1", "2", "3"]
        assert len(warnings) == 1


# ---------------------------------------------------------------------------
# Optional-marker lookup
# ---------------------------------------------------------------------------

class TestOptionalMarker:

    @pytest.fixture
    def grid(self):
        return pd.DataFrame([["Header1", "Value1"], ["Header2", "Value2"]])

    def test_required_missing_logs_error(self, grid, caplog):
        with caplog.at_level(logging.DEBUG, logger="NCollector"):
            assert get_location(grid, "Missing") is None
        assert any(r.levelno == logging.ERROR for r in caplog.records)

    def test_optional_missing_no_error(self, grid, caplog):
        with caplog.at_level(logging.DEBUG, logger="NCollector"):
            assert get_location(grid, "Missing", required=False) is None
            assert extract_value(grid, "Missing", required=False) is None
        assert not any(r.levelno >= logging.WARNING for r in caplog.records)

    def test_optional_occurrence_missing_no_error(self, grid, caplog):
        with caplog.at_level(logging.DEBUG, logger="NCollector"):
            assert get_location(grid, "Header1", match_index=3, required=False) is None
        assert not any(r.levelno >= logging.WARNING for r in caplog.records)

    def test_optional_out_of_bounds_no_error(self, grid, caplog):
        with caplog.at_level(logging.DEBUG, logger="NCollector"):
            assert extract_value(grid, "Value1", col_offset=5, required=False) is None
        assert not any(r.levelno >= logging.WARNING for r in caplog.records)

    def test_optional_found_returns_value(self, grid):
        assert extract_value(grid, "Header2", required=False) == "Value2"


# ---------------------------------------------------------------------------
# Concentration list: eight rows by position
# ---------------------------------------------------------------------------

class TestConcentrationList:
    """The list is read as in 2.0.5: it ends at the first empty cell (decided 14.10.26, so
    layouts with a different number of concentrations need no fixed row count)."""

    def test_list_ends_at_first_empty_cell(self):
        sheet = pd.DataFrame([["final concentration in well (log(M))"],
                              [-6], [-7], [None], [-9]])
        df = slice_table(sheet, "final concentration in well (log(M))")
        assert df.iloc[:, 0].tolist() == [-6, -7]

    @needs_fixtures
    def test_old_protocol_list_and_concentration_map(self):
        """Seven values (A-G); per-row concentrations unchanged (H NaN)."""
        protocol = _load_protocol(OLD_FOLDER)
        assert protocol.ligand_conc.iloc[:, 0].tolist() == [-7, -8, -9, -10, -11, -12, -13]
        conc = built_conc_dic(protocol.ligand_conc)
        assert [conc[r] for r in "ABCDEFG"] == [-7, -8, -9, -10, -11, -12, -13]
        assert pd.isna(conc["H"])

    @needs_fixtures
    def test_one_concentration_protocol_list(self):
        protocol = _load_protocol(NEW_FOLDER)
        assert protocol.ligand and protocol.ligand_2 is None
        assert protocol.ligand_conc.iloc[:, 0].tolist() == [-6, -6, -6]

    @needs_fixtures
    def test_second_ligand_list(self):
        protocol = _load_protocol(SECOND_LIGAND_FOLDER)
        assert protocol.ligand_2 and protocol.ligand_2 != protocol.ligand
        assert protocol.ligand_conc.iloc[:, 0].tolist() == [-5, -6, -6.5, -7, -7.5, -8, -9]
        assert protocol.ligand_2_conc.iloc[:, 0].tolist() == [-6, -7, -8, -9, -10, -11, -12]


# ---------------------------------------------------------------------------
# 'Specify plate layout:' (plate_mode)
# ---------------------------------------------------------------------------

@needs_fixtures
class TestPlateMode:

    def test_new_protocol_reads_field(self):
        protocol = _load_protocol(NEW_FOLDER)
        assert protocol.plate_mode == "one concentration over time"
        assert protocol.plate_mode_inferred is False

    def test_old_protocol_inferred(self, caplog):
        with caplog.at_level(logging.DEBUG, logger="NCollector"):
            protocol = _load_protocol(OLD_FOLDER)
        assert protocol.plate_mode == DEFAULT_PLATE_MODE == "concentration-response over time"
        assert protocol.plate_mode_inferred is True
        # optional markers: no error for the missing field
        assert not any(r.levelno >= logging.ERROR for r in caplog.records)

    def test_old_protocol_scan_logs_layout_warning(self):
        result, logs = _scan(OLD_FOLDER)
        assert result and result[0].protocol.plate_mode_inferred
        expected = (f"   [LAYOUT] Protocol '{_protocol_name(OLD_FOLDER)}': no 'Specify plate layout' field "
                    f"– assumed 'concentration-response over time'")
        assert expected in logs

    def test_new_protocol_scan_no_layout_warning(self):
        result, logs = _scan(NEW_FOLDER)
        assert result[0].protocol.plate_mode == "one concentration over time"
        assert not any("[LAYOUT]" in m for m in logs)

    def test_old_protocol_without_ligand_layout_no_error(self, caplog):
        """'Ligand layout' is optional (old protocols with a second ligand)."""
        with caplog.at_level(logging.DEBUG, logger="NCollector"):
            protocol = _load_protocol(SECOND_LIGAND_FOLDER)
        assert protocol.ligand_layout is None
        assert not any(r.levelno >= logging.ERROR for r in caplog.records)

    def test_empty_field_value_is_inferred(self):
        from parsing import _read_plate_mode
        sheet = pd.DataFrame([["Specify plate layout:"], [None], ["Day 1"]])
        assert _read_plate_mode(sheet) == (DEFAULT_PLATE_MODE, True)
        sheet = pd.DataFrame([["Specify plate layout:"], ["  one concentration over time "]])
        assert _read_plate_mode(sheet) == ("one concentration over time", False)


# ---------------------------------------------------------------------------
# Folder scan: protocol checks before any plate file, skip folder
# ---------------------------------------------------------------------------

@needs_fixtures
class TestScanOrder:

    def test_valid_folders_load(self):
        for folder in (OLD_FOLDER, NEW_FOLDER, SECOND_LIGAND_FOLDER, LABELING_FOLDER):
            result, logs = _scan(folder)
            assert len(result) == 1, folder
            assert result[0].protocol is not None
            assert result[0].results
            assert not any("[SKIP]" in m or "[ERROR]" in m for m in logs), logs

    def test_protocol_read_before_plate_files(self):
        """The protocol is read first, whatever order os.listdir returns."""
        _result, logs = _scan(OLD_FOLDER)
        first_protocol = next(i for i, m in enumerate(logs) if "[PROTOCOL]" in m)
        first_plate = next(i for i, m in enumerate(logs) if "[MEASUREMENT]" in m)
        assert first_protocol < first_plate


# ---------------------------------------------------------------------------
# Plate files: skip (descending ID3 range, no ID3 line) and ID3 mapping end-to-end
# ---------------------------------------------------------------------------

@needs_fixtures
class TestPlateFileId3:

    def test_descending_range_skips_plate_file(self, tmp_path):
        name = "plate with descending ID3.xlsx"
        folder = _old_folder_with_plate(tmp_path, os.path.join(BROKEN, name))
        result, logs = _scan(folder)
        assert len(result) == 1 and result[0].results == []
        assert name in result[0].skipped_files
        assert f"   [SKIP] Plate file '{name}': ID3 '4-1' contains a descending range '4-1'" in logs
        assert not any("[ERROR]" in m for m in logs)

    def test_no_id3_line_skips_plate_file(self, tmp_path):
        name = "plate_no_ID3.xlsx"
        folder = _old_folder_with_plate(tmp_path, os.path.join(BROKEN, name))
        result, logs = _scan(folder)
        assert len(result) == 1 and result[0].results == []
        assert name in result[0].skipped_files
        assert f"   [SKIP] Plate file '{name}': no 'ID3:' line" in logs
        assert not any("[ERROR]" in m for m in logs)

    def _process_single(self, tmp_path, plate_name, caplog):
        folder = _old_folder_with_plate(tmp_path, os.path.join(BROKEN, plate_name))
        result, logs = _scan(folder)
        assert len(result[0].results) == 1, logs
        res, protocol = result[0].results[0], result[0].protocol
        config = ProcessingConfig(labeling_correction=False, plate_layout=TRIPLICATE_LAYOUT)
        with caplog.at_level(logging.WARNING, logger="NCollector"):
            map_plate_metadata(res, protocol, config)
        warnings = [r.getMessage() for r in caplog.records if "[WARNING]" in r.getMessage()]
        return res, protocol, warnings

    def _name(self, protocol, t_id):
        return " + ".join(sorted(protocol.transfection_conditions[t_id]))

    def test_non_numeric_entry_dropped_position_kept(self, tmp_path, caplog):
        res, protocol, warnings = self._process_single(tmp_path, "plate_ID3_non_numeric.xlsx", caplog)
        assert _block_conditions(res) == ["Empty", self._name(protocol, "2"), "Empty", "Empty"]
        assert len(warnings) == 1
        assert warnings[0].startswith("[WARNING] 'plate_ID3_non_numeric.xlsx': ID3 entry '1a' (position 1)")

    def test_id_not_in_protocol_dropped(self, tmp_path, caplog):
        res, protocol, warnings = self._process_single(tmp_path, "plate_ID3_13.xlsx", caplog)
        assert _block_conditions(res) == ["Empty"] * 4
        assert len(warnings) == 1
        assert "ID3 entry '13' (position 1) is not a transfection of the protocol" in warnings[0]

    def test_too_many_ids_extras_ignored(self, tmp_path, caplog):
        res, protocol, warnings = self._process_single(tmp_path, "plate too many ID3.xlsx", caplog)
        assert _block_conditions(res) == [self._name(protocol, t) for t in "1234"]
        assert len(warnings) == 1
        assert "ID3 lists 12 IDs for 4 conditions" in warnings[0]

    def test_empty_entry_partial_plate_no_warning(self, tmp_path, caplog):
        res, protocol, warnings = self._process_single(tmp_path, "plate_ID3_empty_entry.xlsx", caplog)
        assert _block_conditions(res) == [self._name(protocol, "1"), "Empty",
                                          self._name(protocol, "3"), self._name(protocol, "4")]
        assert warnings == []

    def test_valid_plates_no_id3_warning(self, caplog):
        result, _logs = _scan(OLD_FOLDER)
        config = ProcessingConfig(labeling_correction=False, plate_layout=TRIPLICATE_LAYOUT)
        with caplog.at_level(logging.WARNING, logger="NCollector"):
            for res in result[0].results:
                map_plate_metadata(res, result[0].protocol, config)
        assert not any("ID3" in r.getMessage() for r in caplog.records)
        for res in result[0].results:
            assert "Empty" not in _block_conditions(res)