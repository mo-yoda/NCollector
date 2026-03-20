"""
Tests for models.py — constants consistency, build_plate_layout,
dataclass defaults and construction, schema validation.

These tests act as guardrails: if someone adds a new column to MASTER_COLUMNS
but forgets to add it to LEGACY_COLUMN_DEFAULTS (or vice versa), a test fails.
"""
import pytest
import pandas as pd
from datetime import date

from models import (
    APP_VERSION,
    TRIPLICATE_LAYOUT,
    QUADRUPLICATE_LAYOUT,
    build_plate_layout,
    MASTER_COLUMNS,
    LEGACY_COLUMN_DEFAULTS,
    ENRICHABLE_COLS,
    DATA_TYPE_MAP,
    PlateColMetadata,
    PrResult,
    ProtocolData,
    MeasurementFolder,
    ProcessingConfig,
)


# ---------------------------------------------------------------------------
# Tests: Plate Layout Constants
# ---------------------------------------------------------------------------

class TestPlateLayouts:

    def test_triplicate_covers_12_columns(self):
        """Triplicate layout should cover columns 1-12 exactly."""
        all_cols = []
        for block in TRIPLICATE_LAYOUT:
            all_cols.extend(block)
        assert sorted(all_cols) == list(range(1, 13))

    def test_quadruplicate_covers_12_columns(self):
        """Quadruplicate layout should cover columns 1-12 exactly."""
        all_cols = []
        for block in QUADRUPLICATE_LAYOUT:
            all_cols.extend(block)
        assert sorted(all_cols) == list(range(1, 13))

    def test_triplicate_has_4_blocks(self):
        assert len(TRIPLICATE_LAYOUT) == 4

    def test_quadruplicate_has_3_blocks(self):
        assert len(QUADRUPLICATE_LAYOUT) == 3

    def test_triplicate_block_size_is_3(self):
        for block in TRIPLICATE_LAYOUT:
            assert len(block) == 3

    def test_quadruplicate_block_size_is_4(self):
        for block in QUADRUPLICATE_LAYOUT:
            assert len(block) == 4

    def test_no_overlapping_columns_triplicate(self):
        all_cols = []
        for block in TRIPLICATE_LAYOUT:
            all_cols.extend(block)
        assert len(all_cols) == len(set(all_cols))

    def test_no_overlapping_columns_quadruplicate(self):
        all_cols = []
        for block in QUADRUPLICATE_LAYOUT:
            all_cols.extend(block)
        assert len(all_cols) == len(set(all_cols))


# ---------------------------------------------------------------------------
# Tests: build_plate_layout
# ---------------------------------------------------------------------------

class TestBuildPlateLayout:

    def test_labeling_false_returns_triplicate(self):
        result = build_plate_layout(is_labeling=False)
        assert result == TRIPLICATE_LAYOUT

    def test_labeling_true_returns_quadruplicate(self):
        result = build_plate_layout(is_labeling=True)
        assert result == QUADRUPLICATE_LAYOUT


# ---------------------------------------------------------------------------
# Tests: MASTER_COLUMNS
# ---------------------------------------------------------------------------

class TestMasterColumns:

    def test_no_duplicates(self):
        """MASTER_COLUMNS should have no duplicate entries."""
        assert len(MASTER_COLUMNS) == len(set(MASTER_COLUMNS))

    def test_essential_columns_present(self):
        """Core columns that the app depends on should exist."""
        essential = [
            "File_Name", "Date", "Transfection", "Cell_Line", "Ligand",
            "Well_ID", "Time_(min)", "Plate_Row", "Replicate",
            "Raw_BRET_kinetic", "Kinetic_Mean", "AUC_Mean",
            "Is_Excluded", "Is_Vehicle",
        ]
        for col in essential:
            assert col in MASTER_COLUMNS, f"Essential column '{col}' missing from MASTER_COLUMNS"

    def test_all_data_type_map_columns_in_master(self):
        """Every internal column name in DATA_TYPE_MAP should exist in MASTER_COLUMNS."""
        for category, subtypes in DATA_TYPE_MAP.items():
            for display_name, col_name in subtypes.items():
                assert col_name in MASTER_COLUMNS, (
                    f"DATA_TYPE_MAP['{category}']['{display_name}'] = '{col_name}' "
                    f"not in MASTER_COLUMNS"
                )


# ---------------------------------------------------------------------------
# Tests: LEGACY_COLUMN_DEFAULTS
# ---------------------------------------------------------------------------

class TestLegacyColumnDefaults:

    def test_all_keys_are_in_master_columns(self):
        """Every legacy column should exist in MASTER_COLUMNS."""
        for col in LEGACY_COLUMN_DEFAULTS:
            assert col in MASTER_COLUMNS, (
                f"LEGACY_COLUMN_DEFAULTS key '{col}' not in MASTER_COLUMNS"
            )

    def test_all_entries_have_required_keys(self):
        """Each entry should have 'default', 'reconstructable', and 'enrichable'."""
        required = {"default", "reconstructable", "enrichable"}
        for col, config in LEGACY_COLUMN_DEFAULTS.items():
            assert required.issubset(config.keys()), (
                f"LEGACY_COLUMN_DEFAULTS['{col}'] missing keys: {required - config.keys()}"
            )

    def test_reconstructable_and_enrichable_are_bool(self):
        """Flags should be boolean."""
        for col, config in LEGACY_COLUMN_DEFAULTS.items():
            assert isinstance(config["reconstructable"], bool), f"'{col}' reconstructable not bool"
            assert isinstance(config["enrichable"], bool), f"'{col}' enrichable not bool"

    def test_no_column_is_both_reconstructable_and_enrichable(self):
        """A column shouldn't be both — these are different migration paths."""
        for col, config in LEGACY_COLUMN_DEFAULTS.items():
            if config["reconstructable"] and config["enrichable"]:
                pytest.fail(f"'{col}' is both reconstructable and enrichable")


# ---------------------------------------------------------------------------
# Tests: ENRICHABLE_COLS
# ---------------------------------------------------------------------------

class TestEnrichableCols:

    def test_derived_correctly(self):
        """ENRICHABLE_COLS should match the enrichable=True entries."""
        expected = [col for col, cfg in LEGACY_COLUMN_DEFAULTS.items() if cfg["enrichable"]]
        assert ENRICHABLE_COLS == expected

    def test_not_empty(self):
        """Should have at least one enrichable column."""
        assert len(ENRICHABLE_COLS) > 0

    def test_all_are_in_master_columns(self):
        for col in ENRICHABLE_COLS:
            assert col in MASTER_COLUMNS


# ---------------------------------------------------------------------------
# Tests: DATA_TYPE_MAP
# ---------------------------------------------------------------------------

class TestDataTypeMap:

    def test_has_kinetic_and_crc(self):
        assert "kinetic" in DATA_TYPE_MAP
        assert "CRC" in DATA_TYPE_MAP

    def test_no_empty_categories(self):
        for cat, subtypes in DATA_TYPE_MAP.items():
            assert len(subtypes) > 0, f"Category '{cat}' has no subtypes"

    def test_no_duplicate_internal_columns_within_category(self):
        """No two display names within the same category should map to the same column."""
        for cat, subtypes in DATA_TYPE_MAP.items():
            cols = list(subtypes.values())
            assert len(cols) == len(set(cols)), (
                f"Duplicate internal columns within DATA_TYPE_MAP['{cat}']: "
                f"{[c for c in cols if cols.count(c) > 1]}"
            )

    def test_kinetic_mean_exists(self):
        """Kinetic_Mean is used as default export — must be present."""
        assert "Kinetic_Mean" in DATA_TYPE_MAP["kinetic"].values()

    def test_auc_mean_exists(self):
        """AUC_Mean is used as default export — must be present."""
        assert "AUC_Mean" in DATA_TYPE_MAP["CRC"].values()


# ---------------------------------------------------------------------------
# Tests: PlateColMetadata defaults
# ---------------------------------------------------------------------------

class TestPlateColMetadata:

    def test_defaults(self):
        meta = PlateColMetadata()
        assert meta.cell_line == "Unknown"
        assert meta.transfection_id == "N/A"
        assert meta.condition_name == "Empty"
        assert meta.plasmids == []
        assert meta.ligand_identity == "N/A"
        assert meta.ligand_conc == {}
        assert meta.replicate == ""

    def test_mutable_defaults_independent(self):
        """Each instance should get its own list/dict, not shared references."""
        a = PlateColMetadata()
        b = PlateColMetadata()
        a.plasmids.append("X")
        a.ligand_conc["A"] = -9.0
        assert b.plasmids == []
        assert b.ligand_conc == {}


# ---------------------------------------------------------------------------
# Tests: PrResult defaults
# ---------------------------------------------------------------------------

class TestPrResult:

    @pytest.fixture
    def minimal_result(self):
        return PrResult(
            file_name="test.xlsx",
            measurement_date=date(2025, 1, 1),
            cell_line="HEK293",
            transfection_id="1,2",
            raw_time=[0.0, 1.0],
            raw_bret_ratio_df=pd.DataFrame({"A1": [0.2, 0.3]}),
            donor_df=pd.DataFrame({"A1": [500, 520]}),
            acceptor_df=pd.DataFrame({"A1": [100, 110]}),
        )

    def test_required_fields(self, minimal_result):
        assert minimal_result.file_name == "test.xlsx"
        assert minimal_result.measurement_date == date(2025, 1, 1)
        assert minimal_result.cell_line == "HEK293"

    def test_optional_defaults(self, minimal_result):
        assert minimal_result.donor_wavelength == 0
        assert minimal_result.acceptor_wavelength == 0
        assert minimal_result.info_sheet == {}
        assert minimal_result.column_metadata == {}
        assert minimal_result.time_vector == []
        assert minimal_result.is_excluded is False
        assert minimal_result.excluded_wells == []
        assert minimal_result.vehicle_warnings == []
        assert minimal_result.low_lum_warnings == []

    def test_processed_fields_default_none(self, minimal_result):
        assert minimal_result.kinetic_df is None
        assert minimal_result.kinetic_mean_df is None
        assert minimal_result.auc_df is None
        assert minimal_result.auc_mean_df is None
        assert minimal_result.bl_corr_kinetic is None
        assert minimal_result.raw_bret_ratio_cleaned is None

    def test_mutable_defaults_independent(self):
        a = PrResult(
            file_name="a.xlsx", measurement_date=date(2025, 1, 1),
            cell_line="X", transfection_id="1",
            raw_time=[], raw_bret_ratio_df=pd.DataFrame(),
            donor_df=pd.DataFrame(), acceptor_df=pd.DataFrame(),
        )
        b = PrResult(
            file_name="b.xlsx", measurement_date=date(2025, 1, 1),
            cell_line="X", transfection_id="1",
            raw_time=[], raw_bret_ratio_df=pd.DataFrame(),
            donor_df=pd.DataFrame(), acceptor_df=pd.DataFrame(),
        )
        a.excluded_wells.append("A1")
        a.info_sheet["key"] = "val"
        assert b.excluded_wells == []
        assert b.info_sheet == {}


# ---------------------------------------------------------------------------
# Tests: ProcessingConfig defaults
# ---------------------------------------------------------------------------

class TestProcessingConfig:

    def test_defaults(self):
        config = ProcessingConfig()
        assert config.labeling_correction is False
        assert config.lum_threshold == 100
        assert config.vehicle_warning_threshold == 0.2
        assert config.baseline_end_index is None
        assert config.plate_layout == TRIPLICATE_LAYOUT
        assert config.user_input_fn is None

    def test_plate_layout_default_matches_triplicate(self):
        """Default plate_layout should equal TRIPLICATE_LAYOUT."""
        config = ProcessingConfig()
        assert config.plate_layout == TRIPLICATE_LAYOUT

    def test_custom_config(self):
        config = ProcessingConfig(
            labeling_correction=True,
            lum_threshold=50,
            baseline_end_index=5,
            plate_layout=QUADRUPLICATE_LAYOUT,
        )
        assert config.labeling_correction is True
        assert config.lum_threshold == 50
        assert config.baseline_end_index == 5
        assert config.plate_layout == QUADRUPLICATE_LAYOUT


# ---------------------------------------------------------------------------
# Tests: MeasurementFolder defaults
# ---------------------------------------------------------------------------

class TestMeasurementFolder:

    def test_defaults(self):
        folder = MeasurementFolder(
            folder_name="250121_test",
            folder_path="/data/250121_test",
            measurement_date=date(2025, 1, 21),
        )
        assert folder.protocol is None
        assert folder.results == []
        assert folder.skipped_files == []

    def test_mutable_defaults_independent(self):
        a = MeasurementFolder("a", "/a", date(2025, 1, 1))
        b = MeasurementFolder("b", "/b", date(2025, 1, 1))
        a.results.append("x")
        a.skipped_files.append("y")
        assert b.results == []
        assert b.skipped_files == []


# ---------------------------------------------------------------------------
# Tests: ProtocolData construction
# ---------------------------------------------------------------------------

class TestProtocolData:

    def test_minimal_construction(self):
        proto = ProtocolData(
            file_name="proto.xlsx",
            exp_date=date(2025, 1, 1),
            n=3,
            cell_lines=["HEK293"],
            line_layout="one line",
            transfection_scheme=pd.DataFrame(),
            main_plasmids=["PlasmidA"],
            transfection_conditions={"1": ["CondA"]},
            ligand="TestLigand",
            ligand_conc=pd.DataFrame({"conc": [-9, -8]}),
        )
        assert proto.ligand_2 is None
        assert proto.ligand_2_conc is None
        assert proto.ligand_layout is None

    def test_two_ligand_construction(self):
        proto = ProtocolData(
            file_name="proto.xlsx",
            exp_date=date(2025, 1, 1),
            n=3,
            cell_lines=["HEK293"],
            line_layout="half plate",
            transfection_scheme=pd.DataFrame(),
            main_plasmids=["PlasmidA"],
            transfection_conditions={"1": ["CondA"]},
            ligand="Ligand1",
            ligand_conc=pd.DataFrame({"conc": [-9]}),
            ligand_2="Ligand2",
            ligand_2_conc=pd.DataFrame({"conc": [-8]}),
            ligand_layout="alternating",
        )
        assert proto.ligand_2 == "Ligand2"
        assert proto.ligand_layout == "alternating"


# ---------------------------------------------------------------------------
# Tests: APP_VERSION
# ---------------------------------------------------------------------------

class TestAppVersion:

    def test_is_string(self):
        assert isinstance(APP_VERSION, str)

    def test_not_empty(self):
        assert len(APP_VERSION) > 0