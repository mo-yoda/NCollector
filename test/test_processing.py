"""
Tests for the core BRET data processing pipeline.

Covers: calculate_relative_time, get_block_start_for_col, calculate_vehicle_means,
calculate_means_on_meta, apply_baseline_correction, apply_vehicle_normalization,
apply_labeling_correction, check_vehicle_wells, check_luminescence,
and the integration of process_bret_measurement.
"""
import pytest
import numpy as np
import pandas as pd
from datetime import date

from models import PlateColMetadata, PrResult, ProtocolData, ProcessingConfig, TRIPLICATE_LAYOUT
from processing import (
    calculate_relative_time,
    get_block_start_for_col,
    calculate_vehicle_means,
    calculate_means_on_meta,
    apply_baseline_correction,
    apply_vehicle_normalization,
    apply_labeling_correction,
    check_vehicle_wells,
    check_luminescence,
    process_bret_measurement,
)


# ---------------------------------------------------------------------------
# Fixtures: synthetic plate data
# ---------------------------------------------------------------------------

PLATE_BLOCKS = TRIPLICATE_LAYOUT  # [range(1,4), range(4,7), range(7,10), range(10,13)]
ROWS = "ABCDEFGH"
N_TIMEPOINTS = 10
BASELINE_END = 3  # First 3 reads are baseline


def _well_ids():
    """All 96 well IDs: A1..H12."""
    return [f"{r}{c}" for r in ROWS for c in range(1, 13)]


def _make_kinetic_df(n_timepoints=N_TIMEPOINTS, base_value=0.2, noise_seed=42):
    """
    Creates a synthetic BRET-ratio DataFrame (timepoints × wells).
    Mimics raw plate reader output with a baseline phase and kinetic phase.
    Wells in rows A-G get a ligand-dependent increase; row H (vehicle) stays flat.
    """
    rng = np.random.RandomState(noise_seed)
    wells = _well_ids()
    data = {}

    for w in wells:
        row_char = w[0]
        baseline = np.full(BASELINE_END, base_value)
        if row_char == "H":
            # Vehicle: flat
            kinetic = np.full(n_timepoints - BASELINE_END, base_value)
        else:
            # Stimulated: ramp up
            row_idx = ord(row_char) - ord("A")
            stim_factor = 1.0 + (7 - row_idx) * 0.1  # A gets most, G gets least
            kinetic = np.linspace(base_value, base_value * stim_factor,
                                  n_timepoints - BASELINE_END)
        trace = np.concatenate([baseline, kinetic])
        trace += rng.normal(0, 0.001, len(trace))  # Small noise
        data[w] = trace

    return pd.DataFrame(data)


def _make_raw_time(n_timepoints=N_TIMEPOINTS, interval=1.0, gap_at=BASELINE_END):
    """
    Creates a raw time column with uniform intervals and one gap at the baseline end.
    The gap simulates manual ligand addition (interval doubles at that point).
    """
    times = []
    t = 0.0
    for i in range(n_timepoints):
        times.append(t)
        if i == gap_at - 1:
            t += interval * 2  # Gap for ligand addition
        else:
            t += interval
    return pd.Series(times)


def _make_metadata(plate_blocks=PLATE_BLOCKS):
    """Creates PlateColMetadata for each column (1-12)."""
    metadata = {}
    conditions = ["CondA", "CondB", "CondC", "CondD"]
    for i, block in enumerate(plate_blocks):
        cond = conditions[i] if i < len(conditions) else f"Cond{i}"
        for rep_idx, col in enumerate(block):
            metadata[col] = PlateColMetadata(
                cell_line="HEK293",
                transfection_id=str(i + 1),
                condition_name=cond,
                plasmids=[f"Plasmid_{cond}"],
                ligand_identity="TestLigand",
                ligand_conc={r: float(-9 + ord(r) - ord("A")) for r in "ABCDEFG"} | {"H": float("nan")},
                replicate=str(rep_idx + 1),
            )
    return metadata


@pytest.fixture
def kinetic_df():
    """Synthetic kinetic BRET ratio DataFrame."""
    return _make_kinetic_df()


@pytest.fixture
def raw_time():
    """Synthetic raw time series with a gap at baseline end."""
    return _make_raw_time()


@pytest.fixture
def col_metadata():
    """Column metadata for a standard 96-well plate."""
    return _make_metadata()


@pytest.fixture
def processing_config():
    """Standard processing config."""
    return ProcessingConfig(
        plate_layout=PLATE_BLOCKS,
        baseline_end_index=BASELINE_END,
    )


# ---------------------------------------------------------------------------
# Tests: calculate_relative_time
# ---------------------------------------------------------------------------

class TestCalculateRelativeTime:
    """Tests for baseline detection and time vector generation."""

    def test_auto_detect_baseline(self, raw_time):
        """Auto-detects baseline from gap in time intervals."""
        result = calculate_relative_time(raw_time, baseline_end_idx=None)
        assert result is not None
        assert result[BASELINE_END] == 0.0  # First post-baseline read is t=0
        assert all(t < 0 for t in result[:BASELINE_END])  # Baseline is negative
        assert all(t >= 0 for t in result[BASELINE_END:])  # Post-baseline is non-negative

    def test_manual_baseline_index(self, raw_time):
        """Uses provided baseline index directly."""
        result = calculate_relative_time(raw_time, baseline_end_idx=BASELINE_END)
        assert result is not None
        assert result[BASELINE_END] == 0.0

    def test_uniform_intervals_returns_none(self):
        """If all intervals are identical, auto-detection fails."""
        uniform_time = pd.Series([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
        result = calculate_relative_time(uniform_time, baseline_end_idx=None)
        assert result is None

    def test_multiple_irregular_intervals_returns_none(self):
        """When no single unique gap exists (all intervals repeat), returns None."""
        # Intervals: [1.0, 2.5, 1.0, 2.5, 1.0] — both 1.0 and 2.5 repeat
        # drop_duplicates(keep=False) gives empty series → should return None, not crash
        times = pd.Series([0.0, 1.0, 3.5, 4.5, 7.0, 8.0])
        result = calculate_relative_time(times, baseline_end_idx=None)
        assert result is None

    def test_two_distinct_gaps_returns_none(self):
        """Two genuinely unique gaps (neither repeats) should return None."""
        # Intervals: [1.0, 1.0, 3.0, 1.0, 5.0] — two unique non-repeating gaps (3.0 and 5.0)
        times = pd.Series([0.0, 1.0, 2.0, 5.0, 6.0, 11.0])
        result = calculate_relative_time(times, baseline_end_idx=None)
        assert result is None

    def test_baseline_index_too_large(self):
        """Returns None if data has fewer rows than baseline index."""
        short_time = pd.Series([0.0, 1.0, 2.0])
        result = calculate_relative_time(short_time, baseline_end_idx=10)
        assert result is None

    def test_time_vector_length(self, raw_time):
        """Output length matches input length."""
        result = calculate_relative_time(raw_time, baseline_end_idx=BASELINE_END)
        assert len(result) == len(raw_time)

    def test_time_vector_spacing(self, raw_time):
        """Post-baseline spacing should be uniform."""
        result = calculate_relative_time(raw_time, baseline_end_idx=BASELINE_END)
        post_baseline = result[BASELINE_END:]
        diffs = [post_baseline[i + 1] - post_baseline[i] for i in range(len(post_baseline) - 1)]
        assert all(abs(d - diffs[0]) < 1e-10 for d in diffs)

    def test_negative_baseline_index(self):
        """Negative baseline index should still work (treated as valid int)."""
        times = pd.Series([0.0, 1.0, 2.0, 3.0, 4.0])
        # baseline_end_idx=0 means no baseline phase
        result = calculate_relative_time(times, baseline_end_idx=0)
        assert result is not None
        assert result[0] == 0.0  # First read is t=0


# ---------------------------------------------------------------------------
# Tests: get_block_start_for_col
# ---------------------------------------------------------------------------

class TestGetBlockStartForCol:

    def test_first_block(self):
        assert get_block_start_for_col(1, PLATE_BLOCKS) == 1
        assert get_block_start_for_col(2, PLATE_BLOCKS) == 1
        assert get_block_start_for_col(3, PLATE_BLOCKS) == 1

    def test_second_block(self):
        assert get_block_start_for_col(4, PLATE_BLOCKS) == 4
        assert get_block_start_for_col(6, PLATE_BLOCKS) == 4

    def test_last_block(self):
        assert get_block_start_for_col(10, PLATE_BLOCKS) == 10
        assert get_block_start_for_col(12, PLATE_BLOCKS) == 10

    def test_invalid_column(self):
        assert get_block_start_for_col(13, PLATE_BLOCKS) is None
        assert get_block_start_for_col(0, PLATE_BLOCKS) is None


# ---------------------------------------------------------------------------
# Tests: apply_baseline_correction
# ---------------------------------------------------------------------------

class TestBaselineCorrection:

    def test_basic_correction(self):
        """Values should be divided by the mean of baseline rows."""
        data = pd.DataFrame({
            "A1": [2.0, 2.0, 2.0, 4.0, 6.0],  # Baseline mean = 2.0
            "A2": [1.0, 1.0, 1.0, 3.0, 5.0],  # Baseline mean = 1.0
        })
        result = apply_baseline_correction(data, baseline_end_idx=3)

        # A1: 4/2=2.0, 6/2=3.0
        assert result["A1"].iloc[3] == pytest.approx(2.0)
        assert result["A1"].iloc[4] == pytest.approx(3.0)
        # A2: values/1.0 = unchanged
        assert result["A2"].iloc[3] == pytest.approx(3.0)

    def test_baseline_rows_become_one(self):
        """Baseline rows themselves should approximate 1.0 (mean / mean)."""
        data = pd.DataFrame({
            "A1": [2.0, 2.0, 2.0, 4.0],
        })
        result = apply_baseline_correction(data, baseline_end_idx=3)
        for i in range(3):
            assert result["A1"].iloc[i] == pytest.approx(1.0)

    def test_zero_baseline_becomes_nan(self):
        """Division by zero baseline should produce NaN, not infinity."""
        data = pd.DataFrame({
            "A1": [0.0, 0.0, 0.0, 4.0],
        })
        result = apply_baseline_correction(data, baseline_end_idx=3)
        assert pd.isna(result["A1"].iloc[3])

    def test_nan_in_data_propagates(self):
        """NaN values (excluded wells) should remain NaN."""
        data = pd.DataFrame({
            "A1": [2.0, 2.0, 2.0, float("nan"), 6.0],
        })
        result = apply_baseline_correction(data, baseline_end_idx=3)
        assert pd.isna(result["A1"].iloc[3])
        assert result["A1"].iloc[4] == pytest.approx(3.0)


# ---------------------------------------------------------------------------
# Tests: calculate_vehicle_means
# ---------------------------------------------------------------------------

class TestCalculateVehicleMeans:

    def test_kinetic_vehicle_means(self):
        """Kinetic mode: returns Series per block (mean of vehicle wells over time)."""
        data = pd.DataFrame({
            "H1": [1.0, 1.1, 1.2],
            "H2": [1.0, 1.3, 1.0],
            "H3": [1.0, 0.9, 1.2],
            "H4": [2.0, 2.0, 2.0],
            "H5": [2.0, 2.0, 2.0],
            "H6": [2.0, 2.0, 2.0],
        })
        blocks = [range(1, 4), range(4, 7)]
        result = calculate_vehicle_means(data, blocks, excluded_wells=[])

        assert 1 in result  # Block starting at col 1
        assert 4 in result  # Block starting at col 4
        # Block 1: mean of H1, H2, H3 at each timepoint
        assert result[1].iloc[0] == pytest.approx(1.0)
        # Block 2: all 2.0
        assert result[4].iloc[0] == pytest.approx(2.0)

    def test_auc_vehicle_means(self):
        """AUC mode (single row): returns scalar per block."""
        data = pd.DataFrame({
            "H1": [10.0],
            "H2": [12.0],
            "H3": [11.0],
        })
        blocks = [range(1, 4)]
        result = calculate_vehicle_means(data, blocks, excluded_wells=[])
        assert isinstance(result[1], float)
        assert result[1] == pytest.approx(11.0)

    def test_excluded_wells_ignored(self):
        """Excluded vehicle wells should not contribute to the mean."""
        data = pd.DataFrame({
            "H1": [1.0, 1.0],
            "H2": [100.0, 100.0],  # Will be excluded
            "H3": [1.0, 1.0],
        })
        blocks = [range(1, 4)]
        result = calculate_vehicle_means(data, blocks, excluded_wells=["H2"])
        # Mean of H1 and H3 only
        assert result[1].iloc[0] == pytest.approx(1.0)

    def test_all_vehicles_excluded(self):
        """If all vehicle wells are excluded, result is None."""
        data = pd.DataFrame({"H1": [1.0], "H2": [1.0], "H3": [1.0]})
        blocks = [range(1, 4)]
        result = calculate_vehicle_means(data, blocks, excluded_wells=["H1", "H2", "H3"])
        assert result[1] is None


# ---------------------------------------------------------------------------
# Tests: apply_vehicle_normalization
# ---------------------------------------------------------------------------

class TestVehicleNormalization:

    def test_vehicle_row_normalizes_to_one(self):
        """Vehicle wells (row H) should be approximately 1.0 after normalization."""
        # Baseline-corrected data where vehicle is flat at 1.0
        bl_data = pd.DataFrame({
            "A1": [1.0, 2.0, 3.0],
            "H1": [1.0, 1.0, 1.0],
            "A2": [1.0, 2.0, 3.0],
            "H2": [1.0, 1.0, 1.0],
            "A3": [1.0, 2.0, 3.0],
            "H3": [1.0, 1.0, 1.0],
        })
        bl_auc = pd.DataFrame({
            "A1": [6.0], "H1": [3.0],
            "A2": [6.0], "H2": [3.0],
            "A3": [6.0], "H3": [3.0],
        })
        data_df = bl_data.copy()  # Only used for column iteration
        blocks = [range(1, 4)]

        kinetic, auc = apply_vehicle_normalization(
            bl_data, bl_auc, data_df, blocks, excluded_wells=[]
        )
        # Vehicle wells should be 1.0
        assert kinetic["H1"].iloc[0] == pytest.approx(1.0)
        assert kinetic["H2"].iloc[1] == pytest.approx(1.0)
        # Stimulated wells: 2.0/1.0 = 2.0
        assert kinetic["A1"].iloc[1] == pytest.approx(2.0)
        # AUC: 6.0/3.0 = 2.0
        assert auc["A1"].iloc[0] == pytest.approx(2.0)

    def test_excluded_wells_produce_nan_auc(self):
        """Excluded wells with NaN AUC should stay NaN after normalization."""
        bl_data = pd.DataFrame({
            "A1": [float("nan")], "H1": [1.0],
            "A2": [1.0], "H2": [1.0],
            "A3": [1.0], "H3": [1.0],
        })
        bl_auc = pd.DataFrame({
            "A1": [float("nan")], "H1": [1.0],
            "A2": [1.0], "H2": [1.0],
            "A3": [1.0], "H3": [1.0],
        })
        data_df = bl_data.copy()
        blocks = [range(1, 4)]

        kinetic, auc = apply_vehicle_normalization(
            bl_data, bl_auc, data_df, blocks, excluded_wells=[]
        )
        assert pd.isna(auc["A1"].iloc[0])


# ---------------------------------------------------------------------------
# Tests: calculate_means_on_meta
# ---------------------------------------------------------------------------

class TestCalculateMeansOnMeta:

    @pytest.fixture
    def simple_data_and_meta(self):
        """Processed kinetic DF with 3 replicates per condition and metadata."""
        data = pd.DataFrame({
            "A1": [1.0, 2.0], "A2": [3.0, 4.0], "A3": [5.0, 6.0],
            "B1": [10.0, 20.0], "B2": [30.0, 40.0], "B3": [50.0, 60.0],
        })
        meta = {}
        for col in [1, 2, 3]:
            meta[col] = PlateColMetadata(
                condition_name="TestCond",
                cell_line="HEK",
                ligand_identity="Lig1",
                replicate=str(col),
            )
        return data, meta, [range(1, 4)]

    def test_row_mode(self, simple_data_and_meta):
        """Row mode: mean of replicates per row per timepoint."""
        data, meta, blocks = simple_data_and_meta
        result = calculate_means_on_meta(data, blocks, meta, grouping_mode="row")

        key_a = "TestCond|HEK|Lig1|A"
        key_b = "TestCond|HEK|Lig1|B"
        assert key_a in result.columns
        # Mean of A1, A2, A3 at t=0: (1+3+5)/3 = 3.0
        assert result[key_a].iloc[0] == pytest.approx(3.0)
        # Mean of B1, B2, B3 at t=0: (10+30+50)/3 = 30.0
        assert result[key_b].iloc[0] == pytest.approx(30.0)

    def test_block_mode(self, simple_data_and_meta):
        """Block mode: mean of ALL wells in the block."""
        data, meta, blocks = simple_data_and_meta
        result = calculate_means_on_meta(data, blocks, meta, grouping_mode="block")

        key = "TestCond|HEK|Lig1"
        assert key in result.columns
        # Mean of all 6 wells at t=0: (1+3+5+10+30+50)/6 = 16.5
        assert result[key].iloc[0] == pytest.approx(16.5)

    def test_nan_wells_skipped(self):
        """NaN values (excluded wells) should be ignored in mean calculation."""
        data = pd.DataFrame({
            "A1": [float("nan"), float("nan")],
            "A2": [4.0, 6.0],
            "A3": [4.0, 6.0],
        })
        meta = {col: PlateColMetadata(
            condition_name="C", cell_line="X", ligand_identity="L", replicate=str(col)
        ) for col in [1, 2, 3]}

        result = calculate_means_on_meta(data, [range(1, 4)], meta, grouping_mode="row")
        key = "C|X|L|A"
        # Mean of A2, A3 at t=0 (A1 is NaN): (4+4)/2 = 4.0
        assert result[key].iloc[0] == pytest.approx(4.0)


# ---------------------------------------------------------------------------
# Tests: apply_labeling_correction
# ---------------------------------------------------------------------------

class TestLabelingCorrection:

    def test_background_subtracted(self):
        """Labeling control column's mean is subtracted from all wells in the block."""
        # Use a single block of 4 columns (labeling = last col)
        blocks = [range(1, 5)]
        # Include all rows A-H for the block (function iterates ABCDEFGH)
        data = {}
        for r in "ABCDEFGH":
            data[f"{r}1"] = [10.0, 20.0]
            data[f"{r}2"] = [10.0, 20.0]
            data[f"{r}3"] = [10.0, 20.0]
            data[f"{r}4"] = [2.0, 4.0]  # Labeling control
        data = pd.DataFrame(data)

        # Col 4 is the labeling control (last in block of 4)
        meta = {}
        for col in [1, 2, 3]:
            meta[col] = PlateColMetadata(replicate=str(col))
        meta[4] = PlateColMetadata(replicate="labeling control")

        corrected, wells_to_drop = apply_labeling_correction(data, blocks, meta, [])

        # A1: 10 - 2 = 8 at t=0, 20 - 4 = 16 at t=1
        assert corrected["A1"].iloc[0] == pytest.approx(8.0)
        assert corrected["A1"].iloc[1] == pytest.approx(16.0)
        # B1: same values, same correction
        assert corrected["B1"].iloc[0] == pytest.approx(8.0)
        # Control wells are returned for dropping
        assert "A4" in wells_to_drop
        assert "H4" in wells_to_drop


# ---------------------------------------------------------------------------
# Tests: check_vehicle_wells
# ---------------------------------------------------------------------------

class TestCheckVehicleWells:

    def test_no_warnings_for_normal_vehicles(self, col_metadata):
        """Vehicle wells at exactly 1.0 should produce no warnings."""
        data = pd.DataFrame({f"H{c}": [1.0, 1.0, 1.0] for c in range(1, 13)})
        warnings = check_vehicle_wells(data, PLATE_BLOCKS, col_metadata,
                                       excluded_wells=[], acc_vehicle_range=0.2,
                                       date_str="01.01.25")
        assert len(warnings) == 0

    def test_warning_for_deviant_vehicle(self, col_metadata):
        """Vehicle well deviating > threshold should produce a warning."""
        data = pd.DataFrame({f"H{c}": [1.0, 1.0, 1.0] for c in range(1, 13)})
        data["H1"] = [1.5, 1.5, 1.5]  # Deviates by 0.5 > 0.2

        warnings = check_vehicle_wells(data, PLATE_BLOCKS, col_metadata,
                                       excluded_wells=[], acc_vehicle_range=0.2,
                                       date_str="01.01.25")
        assert len(warnings) >= 1
        assert any("H1" in w.get("Display", "") for w in warnings)

    def test_excluded_vehicle_not_checked(self, col_metadata):
        """Excluded vehicle wells should not produce warnings."""
        data = pd.DataFrame({f"H{c}": [1.0, 1.0, 1.0] for c in range(1, 13)})
        data["H1"] = [5.0, 5.0, 5.0]  # Very deviant

        warnings = check_vehicle_wells(data, PLATE_BLOCKS, col_metadata,
                                       excluded_wells=["H1"], acc_vehicle_range=0.2,
                                       date_str="01.01.25")
        # H1 is excluded, so no warning for it
        assert all("H1" not in w.get("Display", "") for w in warnings)


# ---------------------------------------------------------------------------
# Tests: check_luminescence
# ---------------------------------------------------------------------------

class TestCheckLuminescence:

    def test_no_warnings_above_threshold(self, col_metadata):
        """All wells above threshold should produce no warnings."""
        data = pd.DataFrame({f"{r}{c}": [500.0] * 10
                             for r in ROWS for c in range(1, 13)})
        warnings = check_luminescence(data, col_metadata,
                                      lum_threshold=100, date_str="01.01.25")
        assert len(warnings) == 0

    def test_warning_below_threshold(self, col_metadata):
        """Wells below threshold should produce warnings."""
        data = pd.DataFrame({f"{r}{c}": [500.0] * 10
                             for r in ROWS for c in range(1, 13)})
        # Set column 1 very low
        for r in ROWS:
            data[f"{r}1"] = [10.0] * 10

        warnings = check_luminescence(data, col_metadata,
                                      lum_threshold=100, date_str="01.01.25")
        assert len(warnings) > 0


# ---------------------------------------------------------------------------
# Tests: AUC calculation (min_count=1 fix)
# ---------------------------------------------------------------------------

class TestAucExcludedWells:

    def test_excluded_well_auc_is_nan(self):
        """Excluded wells (all NaN) should produce NaN AUC, not 0.0."""
        data = pd.DataFrame({
            "A1": [1.0, 2.0, 3.0, 4.0, 5.0],
            "A2": [float("nan")] * 5,  # Excluded well
        })
        baseline_end = 2
        auc = data.iloc[baseline_end:].sum(min_count=1)
        assert auc["A1"] == pytest.approx(12.0)  # 3+4+5
        assert pd.isna(auc["A2"])  # Not 0.0

    def test_non_excluded_well_auc_correct(self):
        """Non-excluded wells should have correct AUC sum."""
        data = pd.DataFrame({
            "A1": [0.5, 0.5, 1.0, 2.0, 3.0],
        })
        auc = data.iloc[2:].sum(min_count=1)
        assert auc["A1"] == pytest.approx(6.0)


# ---------------------------------------------------------------------------
# Tests: Integration — process_bret_measurement
# ---------------------------------------------------------------------------

class TestProcessBretMeasurement:
    """Integration tests for the full processing pipeline."""

    @pytest.fixture
    def mock_protocol(self):
        """Minimal ProtocolData for the pipeline."""
        conc_df = pd.DataFrame({"conc": [-9, -8, -7, -6, -5, -4, -3, 0]})
        return ProtocolData(
            file_name="test_protocol.xlsx",
            exp_date=date(2025, 1, 1),
            n=1,
            cell_lines=["HEK293"],
            line_layout="one line",
            transfection_scheme=pd.DataFrame(),
            main_plasmids=["PlasmidA", "PlasmidB"],
            transfection_conditions={"1": ["CondA"], "2": ["CondB"],
                                     "3": ["CondC"], "4": ["CondD"]},
            ligand="TestLigand",
            ligand_conc=conc_df,
        )

    @pytest.fixture
    def mock_result(self):
        """Minimal PrResult with synthetic data."""
        kinetic = _make_kinetic_df()
        raw_time = _make_raw_time()

        # Add "Time (min)" column to match expected format
        kinetic_with_time = kinetic.copy()
        kinetic_with_time.insert(0, "Time (min)", raw_time.values)

        donor = kinetic_with_time.copy()
        acceptor = kinetic_with_time.copy()

        return PrResult(
            file_name="test_analysis.xlsx",
            measurement_date=date(2025, 1, 1),
            cell_line="HEK293",
            transfection_id="1,2,3,4",
            raw_time=raw_time.tolist(),
            raw_bret_ratio_df=kinetic_with_time,
            donor_df=donor,
            acceptor_df=acceptor,
            donor_wavelength=475,
            acceptor_wavelength=535,
        )

    def test_pipeline_produces_results(self, mock_result, mock_protocol, processing_config):
        """The full pipeline should populate all result fields."""
        result = process_bret_measurement(mock_result, mock_protocol, processing_config)

        assert result.kinetic_df is not None
        assert result.kinetic_mean_df is not None
        assert result.auc_df is not None
        assert result.auc_mean_df is not None
        assert result.bl_corr_kinetic is not None
        assert result.time_vector is not None
        assert len(result.time_vector) == N_TIMEPOINTS

    def test_pipeline_time_vector(self, mock_result, mock_protocol, processing_config):
        """Time vector should have zero at baseline end."""
        result = process_bret_measurement(mock_result, mock_protocol, processing_config)
        assert result.time_vector[BASELINE_END] == 0.0

    def test_excluded_file_clears_results(self, mock_result, mock_protocol, processing_config):
        """An excluded file should have None for all processed data."""
        mock_result.is_excluded = True
        result = process_bret_measurement(mock_result, mock_protocol, processing_config)

        assert result.kinetic_df is None
        assert result.kinetic_mean_df is None
        assert result.auc_df is None
        assert result.auc_mean_df is None

    def test_excluded_wells_have_nan(self, mock_result, mock_protocol, processing_config):
        """Excluded wells should have NaN in processed data."""
        mock_result.excluded_wells = ["A1", "B1"]
        result = process_bret_measurement(mock_result, mock_protocol, processing_config)

        assert result.raw_bret_ratio_cleaned is not None
        assert pd.isna(result.raw_bret_ratio_cleaned["A1"]).all()
        assert pd.isna(result.raw_bret_ratio_cleaned["B1"]).all()

    def test_baseline_persisted_in_config(self, mock_result, mock_protocol):
        """Baseline index should be persisted in config for subsequent files."""
        config = ProcessingConfig(plate_layout=PLATE_BLOCKS, baseline_end_index=None)
        process_bret_measurement(mock_result, mock_protocol, config)
        assert config.baseline_end_index is not None
        assert config.baseline_end_index == BASELINE_END

    def test_vehicle_warnings_generated(self, mock_result, mock_protocol, processing_config):
        """Pipeline should detect deviant vehicle wells when present."""
        # Make H1 ramp up after baseline (other H wells stay flat ~0.2)
        # After baseline correction: H1 goes to ~2.5 while others stay ~1.0
        # After vehicle norm: H1 deviates significantly from 1.0
        h1_trace = [0.2] * BASELINE_END + [0.5] * (N_TIMEPOINTS - BASELINE_END)
        mock_result.raw_bret_ratio_df["H1"] = h1_trace
        result = process_bret_measurement(mock_result, mock_protocol, processing_config)
        # Should have at least one vehicle warning for H1
        assert len(result.vehicle_warnings) > 0

    def test_lum_check_fires_for_475(self, mock_result, mock_protocol, processing_config):
        """Luminescence check should run when donor wavelength is 475."""
        # Set column 1 donor counts very low
        for r in ROWS:
            mock_result.donor_df[f"{r}1"] = [5.0] * N_TIMEPOINTS
        result = process_bret_measurement(mock_result, mock_protocol, processing_config)
        assert len(result.low_lum_warnings) > 0

    def test_lum_check_runs_for_non_475(self, mock_result, mock_protocol, processing_config):
        """Luminescence check now runs regardless of donor wavelength (the old 475 nm
        gating was removed in v2.0.5; wavelengths are still logged for reference)."""
        mock_result.donor_wavelength = 450
        # Low counts should still be flagged even when the donor is not 475 nm.
        for r in ROWS:
            mock_result.donor_df[f"{r}1"] = [5.0] * N_TIMEPOINTS
        result = process_bret_measurement(mock_result, mock_protocol, processing_config)
        assert len(result.low_lum_warnings) > 0


# ---------------------------------------------------------------------------
# Tests: Edge cases
# ---------------------------------------------------------------------------

class TestEdgeCases:

    def test_single_timepoint(self):
        """Pipeline helpers should handle single-row data (AUC mode)."""
        data = pd.DataFrame({
            "A1": [2.0], "A2": [2.0], "A3": [2.0],
            "H1": [1.0], "H2": [1.0], "H3": [1.0],
        })
        means = calculate_vehicle_means(data, [range(1, 4)], excluded_wells=[])
        assert isinstance(means[1], float)
        assert means[1] == pytest.approx(1.0)

    def test_all_nan_column(self):
        """A fully excluded (all-NaN) column should not crash any function."""
        data = pd.DataFrame({
            "A1": [float("nan")] * 5,
            "A2": [1.0, 2.0, 3.0, 4.0, 5.0],
            "A3": [1.0, 2.0, 3.0, 4.0, 5.0],
            "H1": [float("nan")] * 5,
            "H2": [1.0] * 5,
            "H3": [1.0] * 5,
        })
        bl = apply_baseline_correction(data, baseline_end_idx=2)
        assert pd.isna(bl["A1"]).all()
        assert pd.isna(bl["H1"]).all()