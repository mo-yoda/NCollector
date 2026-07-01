import pandas as pd
from datetime import date
from dataclasses import dataclass, field

# --- Application Constants --- #
APP_VERSION = "N Collector v2.0.5"

# Plate column layouts
TRIPLICATE_LAYOUT = [range(1, 4), # Block 1: Cols 1-3
                     range(4, 7), # Block 2: Cols 4-6
                     range(7, 10), # Block 3: Cols 7-9
                     range(10, 13)] # Block 4: Cols 10-12
QUADRUPLICATE_LAYOUT = [range(1, 5), range(5, 9), range(9, 13)]

def build_plate_layout(is_labeling: bool) -> list[range]:
    """Returns the plate column layout. Quadruplicates for labeling correction, triplicates otherwise."""
    return QUADRUPLICATE_LAYOUT if is_labeling else TRIPLICATE_LAYOUT

# Master DataFrame Column Schema
# Column ordering for the master CSV. Referenced by compile_master_dataframe
# for building the DF and ensure_master_csv_schema for importing older CSVs.
MASTER_COLUMNS = [
    "NCollector_version", "Path", "Info_Sheet",
    "File_Name", "Date", "Main_Plasmids", "Applied_Exclusions", "Is_Excluded", "Is_Vehicle",
    "Transfection", "Cell_Line", "Ligand",
    "Ligand_Conc", "Plate_Row", "Replicate", "Well_ID", "Time_(min)", "PR_Time(min)",
    "Donor_Raw_kinetic", "Acceptor_Raw_kinetic",
    "Raw_BRET_unexcluded",     # Pristine, exclusion-free raw BRET ratio (Acceptor/Donor) - to allow reversability of exclusions
    "Raw_BRET_kinetic", "Lab_BRET_kinetic", "Bl_Corrected_BRET", "Veh_Norm_Kinetic", "Kinetic_Mean",
    "Raw_BRET_CRC", "Lab_LP", "Bl_LP", "Veh_Norm_LP", "LP_Mean",
    "Lab_AUC", "Bl_AUC", "Veh_Norm_AUC", "AUC_Mean"
]

# Condition label for wells that carry only the Main_Plasmids backbone (blank Transfection).
# Must be non-empty (the dropdown cascade treats "" as "nothing selected") and
# not a plasmid token (so Main_Plasmids needs no rewrite). See ensure_master_csv_schema.
MAIN_ONLY_CONDITION = "-"


# Columns added after v1
# with their default values for legacy CSVs and whether they are reconstructable without loading original files again
# enrichable: True if the column can be populated by reloading source xlsx files
LEGACY_COLUMN_DEFAULTS = {
    "NCollector_version": {"default": "< v2", "reconstructable": False, "enrichable": False},
    "Path":               {"default": "undocumented path", "reconstructable": False, "enrichable": False},
    "Info_Sheet":         {"default": "", "reconstructable": False, "enrichable": True},
    "Donor_Raw_kinetic":  {"default": float('nan'), "reconstructable": False, "enrichable": True},
    "Acceptor_Raw_kinetic": {"default": float('nan'), "reconstructable": False, "enrichable": True},
    "Raw_BRET_unexcluded": {"default": float('nan'), "reconstructable": True, "enrichable": False}, # reconstructable only guaranteeable from channels
    "PR_Time(min)":       {"default": float('nan'), "reconstructable": False, "enrichable": True},
    "Is_Vehicle":         {"default": False, "reconstructable": True, "enrichable": False},
    "Is_Excluded":        {"default": False, "reconstructable": True, "enrichable": False},
}

# Columns that can be populated by reloading source xlsx files (derived from LEGACY_COLUMN_DEFAULTS)
ENRICHABLE_COLS = [col for col, cfg in LEGACY_COLUMN_DEFAULTS.items() if cfg["enrichable"]]


# Dictionary defining which dropdown option corresponds to which column in Master df
DATA_TYPE_MAP = {
    "kinetic": {
        "raw BRET ratio": "Raw_BRET_kinetic",
        "labeling-corrected BRET ratio": "Lab_BRET_kinetic",
        "baseline-corrected BRET ratio": "Bl_Corrected_BRET",
        "vehicle-normalised BRET ratio, techn. replicates": "Veh_Norm_Kinetic",
        "vehicle-normalised BRET ratio, mean of techn. replicates": "Kinetic_Mean"
    },
    "CRC": {
        "last 3x timepoints: raw BRET": "Raw_BRET_CRC",
        "last 3x timepoints: labeling-corrected BRET ratio": "Lab_LP",
        "last 3x timepoints: baseline-corrected BRET ratio": "Bl_LP",
        "last 3x tp: vehicle-normalised BRET ratio, techn. replicates": "Veh_Norm_LP",
        "last 3x tp: vehicle-normalised BRET ratio, mean of techn. replicates": "LP_Mean",
        "AUC: vehicle-normalised BRET ratio, techn. replicates": "Veh_Norm_AUC",
        "AUC: vehicle-normalised BRET ratio, mean of techn. replicates": "AUC_Mean"
    },
    "bargraph": {
        "vehicle-normalised AUC, techn. replicates": "Veh_Norm_AUC",
        "vehicle-normalised AUC, mean of techn. replicates": "AUC_Mean"
    },
    "heatmap": {
        "vehicle-normalised AUC, mean of techn. replicates": "AUC_Mean"
    }
}

# Categories that require exactly one concentration selection
SINGLE_CONC_CATEGORIES = {"bargraph", "heatmap"}
# Categories that require a group_by selection (no "None" option)
REQUIRES_GROUP_BY = {"heatmap"}

@dataclass
class PlateColMetadata:
    """Identity of a specific column"""
    cell_line: str = "Unknown"
    transfection_id: str  = "N/A"
    condition_name: str = "Empty"
    plasmids: list[str] = field(default_factory=list)
    ligand_identity: str = "N/A"
    # Dic mapping row A-H to concentration (float)
    ligand_conc: dict[str, float] = field(default_factory=dict)
    replicate: str = ""

@dataclass
class PrResult:
    """ Raw and processed information from a single _analysis file """
    # --- Information from analysis xlsx itself ---
    file_name: str
    measurement_date: date
    cell_line: str # ID2
    transfection_id: str # ID3
    raw_time: list[float] # Time as extracted
    raw_bret_ratio_df: pd.DataFrame # RET ratio as extracted, no exclusions applied
    donor_df: pd.DataFrame # Raw counts from donor channel (lower wavelength)
    acceptor_df: pd.DataFrame # Raw counts from acceptor channel (higher wavelength)
    donor_wavelength: int = 0
    acceptor_wavelength: int = 0

    # --- Optional metadata from "Protocol Information" sheet in analysis xlsx ---
    info_sheet: dict = field(default_factory=dict)

    # --- Connection to protocol file ---
    # Key = Column Index (1-12), Value = WellMetadata object
    column_metadata: dict[int, PlateColMetadata] = field(default_factory=dict)

    # --- Processed BRET data ---
    # Time column
    time_vector: list[float] = field(default_factory=list)
    # Raw BRET ratio with applied exclusions
    raw_bret_ratio_cleaned: pd.DataFrame | None = None
    # Labeling corrected kinetic data
    labeling_corr_kinetic: pd.DataFrame | None = None
    # Baseline corrected kinetic data
    bl_corr_kinetic: pd.DataFrame | None = None
    # Baseline- and vehicle-normalised kinetic data (technical replicates)
    kinetic_df: pd.DataFrame | None = None
    # Mean of baseline- and vehicle-normalised kinetic data
    kinetic_mean_df: pd.DataFrame | None = None

    # Raw BRET from last 3x datapoints
    raw_bret_points_df: pd.DataFrame | None = None
    # From last 3x datapoints: pre-baseline and vehicle norm but after labeling correction
    labeling_corr_lp_df: pd.DataFrame | None = None
    # From last 3x datapoints: pre-vehicle norm
    bl_corr_lp_df: pd.DataFrame | None = None
    # From last 3x datapoints: baseline- and vehicle-normalised data
    lp_df: pd.DataFrame | None = None
    # From last 3x datapoints: mean of baseline- and vehicle-normalised data
    lp_mean_df: pd.DataFrame | None = None

    # Pre-baseline and vehicle norm AUC but after labeling correction
    labeling_corr_auc_df: pd.DataFrame | None = None
    # Pre-vehicle norm AUC
    bl_corr_auc_df: pd.DataFrame | None = None

    # Baseline- and vehicle-normalised AUC data (technical replicates)
    auc_df: pd.DataFrame | None = None
    # Mean of baseline- and vehicle-normalised AUC data
    auc_mean_df:pd.DataFrame | None = None

    # --- Optional exclusion by user interaction  ---
    is_excluded: bool = False
    excluded_wells: list[str] = field(default_factory=list)

    # --- Internal check and warnings for outlier identification ---
    vehicle_warnings: list[dict] = field(default_factory=list)
    low_lum_warnings: list[dict] = field(default_factory=list) # List of dict carrying all needed metadata

@dataclass
class ProtocolData:
    """Information from a protocol file"""
    file_name: str
    exp_date: date
    n: int
    cell_lines: list[str]
    line_layout: str
    transfection_scheme: pd.DataFrame
    main_plasmids: list[str] # Plasmids transfected in all conditions
    # keys are col strings, values are transfected plasmids
    transfection_conditions: dict[str, list[str]]
    ligand: str
    ligand_conc: pd.DataFrame
    # Second ligand is optional; by | None = None
    ligand_2: str | None = None
    ligand_2_conc: pd.DataFrame | None = None
    ligand_layout: str | None = None # Only available in protocols > 1.03

@dataclass
class MeasurementFolder:
    """A subfolder containing one protocol and multiple result files"""
    folder_name: str
    folder_path: str
    measurement_date: date
    protocol: ProtocolData | None = None # MeasurementFolder is initiated before protocol data is loaded
    results: list[PrResult] = field(default_factory=list) # The default_factory=list initiates this with an empty list
    skipped_files: list[str] = field(default_factory=list)

@dataclass
class ProcessingConfig:
    """All user-defined processing settings"""
    labeling_correction: bool = field(default=False)
    lum_threshold: int = 100
    vehicle_warning_threshold: float = 0.2
    baseline_end_index: int | None = None
    # Default to triplicates
    plate_layout: list[range] = field(default_factory=lambda: TRIPLICATE_LAYOUT)

    # Optional callback for requesting user input (set by GUI layer)
    # Signature: fn(title: str, message: str, input_type: str, default) -> value | None
    user_input_fn: object = field(default=None, repr=False)
    # Signature: fn(ligand_1_name: str, ligand_2_name: str, plate_info: str) -> "L1" | "L2" | None
    ligand_choice_fn: object = field(default=None, repr=False)
    # Signature: fn(ligand_1_name: str, ligand_2_name: str, protocol_name: str) -> "half" | "alternating" | "one ligand" | None
    ligand_layout_fn: object = field(default=None, repr=False)


@dataclass
class BretMetricsBundle:
    """
    Return bundle for processing.compute_bret_metrics.

    Holds every processed DataFrame that step 11 of process_bret_measurement assigns
    onto a PrResult, plus the cleaned raw BRET ratio. All of these are pure functions
    of the per-well raw BRET ratio and the set of excluded wells, so the same bundle
    can be produced either from a live PrResult (object pipeline) or reconstructed
    purely from a flat master DataFrame (master-native recompute) — no original xlsx
    files required.

    Note: raw_bret_ratio_cleaned here has NO "Time (min)" column (it is the wide
    per-well matrix the math operates on). process_bret_measurement keeps assigning its
    own time-bearing raw_bret_ratio_cleaned onto the PrResult for backwards compatibility.
    """
    # Cleaned raw BRET ratio (excluded wells NaN'd, columns = Well_IDs, no time column)
    raw_bret_ratio_cleaned: pd.DataFrame

    # --- Kinetic traces (one row per timepoint, one column per well) ---
    labeling_corr_kinetic: pd.DataFrame
    bl_corr_kinetic: pd.DataFrame
    kinetic_df: pd.DataFrame
    kinetic_mean_df: pd.DataFrame

    # --- Last-3-timepoint reductions (CRC: one value per well / per replicate-mean) ---
    raw_bret_points_df: pd.DataFrame
    labeling_corr_lp_df: pd.DataFrame
    bl_corr_lp_df: pd.DataFrame
    lp_df: pd.DataFrame
    lp_mean_df: pd.DataFrame

    # --- AUC (one value per well / per replicate-mean) ---
    labeling_corr_auc_df: pd.DataFrame
    bl_corr_auc_df: pd.DataFrame
    auc_df: pd.DataFrame
    auc_mean_df: pd.DataFrame