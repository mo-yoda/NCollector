"""
Tests for mapping.py — plate layout mapping for cell lines, ligands,
transfections, and concentration dictionaries.

Covers: generate_col_mapping, get_cell_line_map, built_conc_dic,
get_ligand_map, get_transfection_map.
"""
import math
import pytest
import pandas as pd
from datetime import date

from models import ProtocolData
from mapping import (
    generate_col_mapping,
    get_cell_line_map,
    built_conc_dic,
    get_ligand_map,
    get_transfection_map,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_protocol(**overrides) -> ProtocolData:
    """Creates a minimal ProtocolData with sensible defaults, overridable."""
    defaults = dict(
        file_name="test_protocol.xlsx",
        exp_date=date(2025, 1, 1),
        n=3,
        cell_lines=["HEK293"],
        line_layout="one line",
        transfection_scheme=pd.DataFrame(),
        main_plasmids=["PlasmidA"],
        transfection_conditions={"1": ["CondA"]},
        ligand="Ligand1",
        ligand_conc=pd.DataFrame({"conc": [-9, -8, -7, -6, -5, -4, -3, 0]}),
        ligand_2=None,
        ligand_2_conc=None,
        ligand_layout=None,
    )
    defaults.update(overrides)
    return ProtocolData(**defaults)


# ---------------------------------------------------------------------------
# Tests: generate_col_mapping
# ---------------------------------------------------------------------------

class TestGenerateColMapping:

    def test_one_line_all_item1(self):
        """'one line' layout: all 12 columns get item_1."""
        result = generate_col_mapping("one line", "A", "B", block_count=4)
        assert all(result[c] == "A" for c in range(1, 13))

    def test_half_plate(self):
        """'half' layout: cols 1-6 get item_1, cols 7-12 get item_2."""
        result = generate_col_mapping("half plate", "A", "B")
        for c in range(1, 7):
            assert result[c] == "A"
        for c in range(7, 13):
            assert result[c] == "B"

    def test_alternating_4_blocks(self):
        """'alternating' with 4 blocks: A,B,A,B pattern in groups of 3."""
        result = generate_col_mapping("alternating", "A", "B", block_count=4)
        # Block 1 (cols 1-3): A
        assert result[1] == "A"
        assert result[2] == "A"
        assert result[3] == "A"
        # Block 2 (cols 4-6): B
        assert result[4] == "B"
        assert result[5] == "B"
        assert result[6] == "B"
        # Block 3 (cols 7-9): A
        assert result[7] == "A"
        # Block 4 (cols 10-12): B
        assert result[10] == "B"

    def test_alternating_3_blocks(self):
        """'alternating' with 3 blocks: A,B,A pattern in groups of 4."""
        result = generate_col_mapping("alternating", "A", "B", block_count=3)
        # Block 1 (cols 1-4): A
        for c in range(1, 5):
            assert result[c] == "A"
        # Block 2 (cols 5-8): B
        for c in range(5, 9):
            assert result[c] == "B"
        # Block 3 (cols 9-12): A
        for c in range(9, 13):
            assert result[c] == "A"

    def test_covers_all_12_columns(self):
        """Every layout should map exactly columns 1-12."""
        for style in ["one line", "half", "alternating"]:
            result = generate_col_mapping(style, "A", "B", block_count=4)
            assert sorted(result.keys()) == list(range(1, 13)), f"Failed for '{style}'"

    def test_unrecognized_fallback(self):
        """Unknown layout string should default to all item_1."""
        result = generate_col_mapping("some_unknown_layout", "A", "B")
        assert all(result[c] == "A" for c in range(1, 13))

    def test_case_insensitive(self):
        """Layout matching should be case-insensitive."""
        result = generate_col_mapping("HALF Plate", "A", "B")
        assert result[1] == "A"
        assert result[7] == "B"


# ---------------------------------------------------------------------------
# Tests: get_cell_line_map
# ---------------------------------------------------------------------------

class TestGetCellLineMap:

    def test_single_cell_line_one_line(self):
        """Single cell line with 'one line' layout → all columns same."""
        protocol = _make_protocol(line_layout="one line")
        result = get_cell_line_map(protocol, "HEK293", block_count=4)
        assert all(result[c] == "HEK293" for c in range(1, 13))

    def test_two_cell_lines_half(self):
        """Two cell lines with 'half' layout → split at column 7."""
        protocol = _make_protocol(line_layout="half plate")
        result = get_cell_line_map(protocol, "HEK293, dQ", block_count=4)
        for c in range(1, 7):
            assert result[c] == "HEK293"
        for c in range(7, 13):
            assert result[c] == "dQ"

    def test_two_cell_lines_alternating(self):
        """Two cell lines with 'alternating' layout → A,B,A,B blocks."""
        protocol = _make_protocol(line_layout="alternating")
        result = get_cell_line_map(protocol, "HEK293, dQ", block_count=4)
        assert result[1] == "HEK293"
        assert result[4] == "dQ"
        assert result[7] == "HEK293"
        assert result[10] == "dQ"

    def test_single_cell_line_with_half(self):
        """Single cell line with 'half' layout → second line is Unknown_2."""
        protocol = _make_protocol(line_layout="half plate")
        result = get_cell_line_map(protocol, "HEK293", block_count=4)
        assert result[1] == "HEK293"
        assert result[7] == "Unknown_2"

    def test_whitespace_handling(self):
        """Whitespace around cell line names should be stripped."""
        protocol = _make_protocol(line_layout="half plate")
        result = get_cell_line_map(protocol, "  HEK293 ,  dQ  ", block_count=4)
        assert result[1] == "HEK293"
        assert result[7] == "dQ"


# ---------------------------------------------------------------------------
# Tests: built_conc_dic
# ---------------------------------------------------------------------------

class TestBuiltConcDic:

    def test_standard_8_row(self):
        """Standard 8-value concentration table → rows A-G get values, H gets NaN."""
        df = pd.DataFrame({"conc": [-9, -8, -7, -6, -5, -4, -3, 0]})
        result = built_conc_dic(df)
        assert result["A"] == -9.0
        assert result["B"] == -8.0
        assert result["G"] == -3.0
        assert math.isnan(result["H"])  # Vehicle row

    def test_vehicle_row_is_nan(self):
        """Row H should always be NaN regardless of input value."""
        df = pd.DataFrame({"conc": [-9, -8, -7, -6, -5, -4, -3, 99]})
        result = built_conc_dic(df)
        assert math.isnan(result["H"])

    def test_fewer_than_8_values(self):
        """Missing rows should get NaN."""
        df = pd.DataFrame({"conc": [-9, -8, -7]})
        result = built_conc_dic(df)
        assert result["A"] == -9.0
        assert result["C"] == -7.0
        assert math.isnan(result["D"])  # Missing → NaN
        assert math.isnan(result["H"])  # Vehicle

    def test_non_numeric_values(self):
        """Non-numeric values should become NaN."""
        df = pd.DataFrame({"conc": [-9, "invalid", -7, -6, -5, -4, -3, 0]})
        result = built_conc_dic(df)
        assert result["A"] == -9.0
        assert math.isnan(result["B"])  # "invalid" → NaN
        assert result["C"] == -7.0

    def test_empty_dataframe(self):
        """Empty DataFrame should return empty dict."""
        result = built_conc_dic(pd.DataFrame())
        assert result == {}

    def test_none_input(self):
        """None input should return empty dict."""
        result = built_conc_dic(None)
        assert result == {}

    def test_all_rows_present(self):
        """Result should have exactly keys A-H."""
        df = pd.DataFrame({"conc": [-9, -8, -7, -6, -5, -4, -3, 0]})
        result = built_conc_dic(df)
        assert sorted(result.keys()) == list("ABCDEFGH")


# ---------------------------------------------------------------------------
# Tests: get_ligand_map
# ---------------------------------------------------------------------------

class TestGetLigandMap:

    def test_single_ligand(self):
        """No second ligand → all columns get L1."""
        protocol = _make_protocol(ligand="Lig1", ligand_2=None)
        result = get_ligand_map(protocol, block_count=4)
        assert all(result[c] == "L1" for c in range(1, 13))

    def test_two_ligands_with_half_layout(self):
        """Explicit 'half' ligand layout → L1 cols 1-6, L2 cols 7-12."""
        protocol = _make_protocol(
            ligand="Lig1", ligand_2="Lig2", ligand_layout="half plate"
        )
        result = get_ligand_map(protocol, block_count=4)
        for c in range(1, 7):
            assert result[c] == "L1"
        for c in range(7, 13):
            assert result[c] == "L2"

    def test_two_ligands_with_alternating_layout(self):
        """Explicit 'alternating' ligand layout."""
        protocol = _make_protocol(
            ligand="Lig1", ligand_2="Lig2", ligand_layout="alternating"
        )
        result = get_ligand_map(protocol, block_count=4)
        assert result[1] == "L1"
        assert result[4] == "L2"
        assert result[7] == "L1"
        assert result[10] == "L2"

    def test_two_ligands_inferred_from_half_cells(self):
        """No ligand_layout, cells=half → ligands should be alternating."""
        protocol = _make_protocol(
            line_layout="half plate",
            ligand="Lig1", ligand_2="Lig2", ligand_layout=None,
        )
        result = get_ligand_map(protocol, block_count=4)
        # Half cells → alternating ligands
        assert result[1] == "L1"
        assert result[4] == "L2"
        assert result[7] == "L1"
        assert result[10] == "L2"

    def test_two_ligands_inferred_from_alternating_cells(self):
        """No ligand_layout, cells=alternating → ligands should be half."""
        protocol = _make_protocol(
            line_layout="alternating",
            ligand="Lig1", ligand_2="Lig2", ligand_layout=None,
        )
        result = get_ligand_map(protocol, block_count=4)
        # Alternating cells → half ligands
        for c in range(1, 7):
            assert result[c] == "L1"
        for c in range(7, 13):
            assert result[c] == "L2"

    def test_two_ligands_one_line_fallback(self):
        """No ligand_layout, cells=one line → fallback to all L1."""
        protocol = _make_protocol(
            line_layout="one line",
            ligand="Lig1", ligand_2="Lig2", ligand_layout=None,
        )
        result = get_ligand_map(protocol, block_count=4)
        assert all(result[c] == "L1" for c in range(1, 13))


# ---------------------------------------------------------------------------
# Tests: get_transfection_map
# ---------------------------------------------------------------------------

class TestGetTransfectionMap:

    # --- No ligand layout ---

    def test_one_line_no_ligand_unique(self):
        """one line, no ligand layout → unique: [1,2,3,4]."""
        result = get_transfection_map("one line", None, ["1", "2", "3", "4"])
        assert result == ["1", "2", "3", "4"]

    def test_half_no_ligand_seq_repeat(self):
        """half, no ligand layout → seq_repeat: [1,2,1,2]."""
        result = get_transfection_map("half plate", None, ["1", "2"])
        assert result == ["1", "2", "1", "2"]

    def test_alternating_no_ligand_elem_repeat(self):
        """alternating, no ligand layout → elem_repeat: [1,1,2,2]."""
        result = get_transfection_map("alternating", None, ["1", "2"])
        assert result == ["1", "1", "2", "2"]

    def test_empty_ids_returns_na(self):
        """Empty ID list → all N/A."""
        result = get_transfection_map("one line", None, [])
        assert result == ["N/A"] * 4

    def test_single_id_repeated(self):
        """Single ID → repeated for all blocks."""
        result = get_transfection_map("half plate", None, ["1"])
        assert result == ["1", "1", "1", "1"]

    def test_more_ids_than_blocks(self):
        """More IDs than blocks → truncated to block_count."""
        result = get_transfection_map("one line", None, ["1", "2", "3", "4", "5", "6"])
        assert len(result) == 4
        assert result == ["1", "2", "3", "4"]

    # --- With ligand layout ---

    def test_one_line_one_ligand(self):
        """one line + one ligand → unique."""
        result = get_transfection_map("one line", "one ligand", ["1", "2", "3", "4"])
        assert result == ["1", "2", "3", "4"]

    def test_one_line_half_ligand(self):
        """one line + half ligand → seq_repeat."""
        result = get_transfection_map("one line", "half", ["1", "2"])
        assert result == ["1", "2", "1", "2"]

    def test_one_line_alternating_ligand(self):
        """one line + alternating ligand → elem_repeat."""
        result = get_transfection_map("one line", "alternating", ["1", "2"])
        assert result == ["1", "1", "2", "2"]

    def test_half_one_ligand(self):
        """half + one ligand → seq_repeat."""
        result = get_transfection_map("half", "one ligand", ["1", "2"])
        assert result == ["1", "2", "1", "2"]

    def test_half_alternating_ligand(self):
        """half + alternating → single (all same). The single mapping is only applied
        when exactly one transfection ID is provided (safety guard); more IDs fall back
        to unique."""
        result = get_transfection_map("half", "alternating", ["1"])
        assert result == ["1", "1", "1", "1"]

    def test_alternating_one_ligand(self):
        """alternating + one ligand → elem_repeat."""
        result = get_transfection_map("alternating", "one ligand", ["1", "2"])
        assert result == ["1", "1", "2", "2"]

    def test_alternating_half_ligand(self):
        """alternating + half → single (only when exactly one transfection ID is given)."""
        result = get_transfection_map("alternating", "half", ["1"])
        assert result == ["1", "1", "1", "1"]

    # --- Block count variations ---

    def test_3_blocks(self):
        """Should work with 3 blocks (labeling correction layout)."""
        result = get_transfection_map("one line", None, ["1", "2", "3"], block_count=3)
        assert result == ["1", "2", "3"]
        assert len(result) == 3

    def test_3_blocks_seq_repeat(self):
        """seq_repeat with 3 blocks and single ID → all same."""
        result = get_transfection_map("half", None, ["1"], block_count=3)
        assert result == ["1", "1", "1"]

    def test_3_blocks_seq_repeat_2_ids_padded(self):
        """seq_repeat with 3 blocks and 2 IDs → pad with N/A."""
        result = get_transfection_map("half", None, ["1", "2"], block_count=3)
        assert result == ["1", "2", "N/A"]

    # --- Half + half special case ---

    def test_half_half_fewer_ids(self):
        """half + half, fewer IDs than blocks → seq_repeat."""
        result = get_transfection_map("half", "half", ["1", "2"])
        assert result == ["1", "2", "1", "2"]

    def test_half_half_enough_ids(self):
        """half + half, enough IDs → unique."""
        result = get_transfection_map("half", "half", ["1", "2", "3", "4"])
        assert result == ["1", "2", "3", "4"]

    # --- Alternating + alternating special case ---

    def test_alt_alt_fewer_ids(self):
        """alternating + alternating, fewer IDs → elem_repeat."""
        result = get_transfection_map("alternating", "alternating", ["1", "2"])
        assert result == ["1", "1", "2", "2"]

    def test_alt_alt_enough_ids(self):
        """alternating + alternating, enough IDs → unique."""
        result = get_transfection_map("alternating", "alternating", ["1", "2", "3", "4"])
        assert result == ["1", "2", "3", "4"]

    # --- Result length always matches block_count ---

    @pytest.mark.parametrize("cell,ligand,ids,bc", [
        ("one line", None, ["1", "2"], 4),
        ("half", None, ["1"], 4),
        ("alternating", None, ["1", "2", "3"], 4),
        ("one line", "half", ["1", "2"], 3),
        ("half", "alternating", ["1"], 3),
    ])
    def test_result_length_matches_block_count(self, cell, ligand, ids, bc):
        """Output length should always equal block_count."""
        result = get_transfection_map(cell, ligand, ids, block_count=bc)
        assert len(result) == bc