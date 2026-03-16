import pandas as pd
from datetime import date
from dataclasses import dataclass, field

# --- Application Constants --- #
APP_VERSION = "N Collector v2.0 Beta"

# Master DataFrame Column Schema
# Column ordering for the master CSV. Referenced by compile_master_dataframe
# for building the DF and ensure_master_csv_schema for importing older CSVs.
MASTER_COLUMNS = [
    "NCollector_version", "Path",
    "File_Name", "Date", "Main_Plasmids", "Applied_Exclusions",
    "Transfection", "Cell_Line", "Ligand",
    "Ligand_Conc", "Plate_Row", "Replicate", "Well_ID", "Time_(min)",
    "Raw_BRET_kinetic", "Lab_BRET_kinetic", "Bl_Corrected_BRET", "Veh_Norm_Kinetic", "Kinetic_Mean",
    "Raw_BRET_CRC", "Lab_LP", "Bl_LP", "Veh_Norm_LP", "LP_Mean",
    "Lab_AUC", "Bl_AUC", "Veh_Norm_AUC", "AUC_Mean"
]

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
    }
}

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
    raw_bret_ratio_df: pd.DataFrame
    lum_df: pd.DataFrame

    # --- Connection to protocol file ---
    # Key = Column Index (1-12), Value = WellMetadata object
    column_metadata: dict[int, PlateColMetadata] = field(default_factory=dict)

    # --- Processed BRET data ---
    # Time column
    time_vector: list[float] = field(default_factory=list)
    # Raw BRET ratio after with applied exclusions
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
    lum_threshold: int=100
    vehicle_warning_threshold: float = 0.2
    baseline_end_index: int=5
    # Default to triplicates
    plate_layout: list[range] = field(default_factory=lambda: [
        range(1, 4),  # Block 1: Cols 1-3
        range(4, 7),  # Block 2: Cols 4-6
        range(7, 10),  # Block 3: Cols 7-9
        range(10, 13)  # Block 4: Cols 10-12
    ])
