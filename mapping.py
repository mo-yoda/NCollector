import logging
import pandas as pd

from models import ProtocolData, ProcessingConfig

logger = logging.getLogger("NCollector")

# --- Plate Layout Mapping --- #

def generate_col_mapping(layout_style: str, item_1: str, item_2: str, block_count: int = 4) -> dict:
    """Shared helper to map two items across 12 columns based on standard layouts."""
    mapping = {}
    style = str(layout_style).lower()

    if "half" in style:
        # First 6 cols Item 1, Last 6 cols Item 2
        for col in range(1, 7): mapping[col] = item_1
        for col in range(7, 13): mapping[col] = item_2

    elif "alternating" in style:
        # Calculate how many columns are in each block
        block_size = 12 // block_count

        # Loop through each block index (0, 1, 2, etc.)
        for i in range(block_count):
            # Even blocks get item_1, odd blocks get item_2
            current_item = item_1 if i % 2 == 0 else item_2

            # Calculate the start and end column for this specific block
            start_col = (i * block_size) + 1
            end_col = start_col + block_size

            # Map the columns
            for col in range(start_col, end_col):
                mapping[col] = current_item

    else:
        # "one line", "one ligand", or unrecognized fallback -> All cols Item 1
        for col in range(1, 13): mapping[col] = item_1

    return mapping

def get_cell_line_map(protocol: ProtocolData, cell_lines: str, block_count: int = 4):
    """
    Defines the plate layout for cell lines based on the dropdown selection protocol (.line_layout)
    and cell_lines in ID2 of the plate reader metadata (PrResult.cell_line)
    """
    logger.debug(f"--- Mapping Cell Lines ---")
    logger.debug(f"Layout Type: '{protocol.line_layout}' | Raw ID2: '{cell_lines}'")
    # Split ID2 string to get potentially multiple cell lines
    lines = [x.strip() for x in cell_lines.split(',')]

    # Safely assign line names or default to Unknown
    line_1 = lines[0] if len(lines) > 0 else "Unknown_1"
    line_2 = lines[1] if len(lines) > 1 else "Unknown_2"

    return generate_col_mapping(protocol.line_layout, line_1, line_2, block_count)

def built_conc_dic(df_conc: pd.DataFrame):
    """
    Parses the ligand concentration DataFrame (from ProtocolData) into a dict.
    Assumes standard 8-row layout corresponding to A-H
    """

    if df_conc is None or df_conc.empty:
        return {}

    conc_dic = {}
    rows = "ABCDEFGH"

    try:
        # Transform first col in df_conc to list; errors='coerce' turns non-numbers to NaN
        vals = pd.to_numeric(df_conc.iloc[:, 0], errors='coerce').tolist()

        for i, row_char in enumerate(rows):
            if row_char == "H":
                # Vehicle row: no ligand, concentration is not meaningful
                conc_dic[row_char] = float('nan')
            elif i < len(vals):
                # Store float if valid, else NaN (value missing in protocol)
                conc_dic[row_char] = float(vals[i]) if not pd.isna(vals[i]) else float('nan')
            else:
                conc_dic[row_char] = float('nan')
    except Exception as e:
        logger.warning(f"Error parsing concentration table: {e}")

    return conc_dic

def get_ligand_map(protocol: ProtocolData, block_count: int = 4,
                   config: ProcessingConfig | None = None, plate_info: str = ""):
    """
    Determines which columns contain Ligand 1 and which contain Ligand 2.
    Returns dic of col idx and ligands.
    """

    # If no second ligand exists, everything is Ligand 1
    if not protocol.ligand_2:
        return generate_col_mapping("one", 'L1', 'L2', block_count)

    # Unselected ligand_layout will be float('nan'); unfound ligand layout (old prtl) returns None
    # Normalize NaN to None
    if protocol.ligand_layout is not None and str(protocol.ligand_layout).lower() == "nan":
        protocol.ligand_layout = None

    # If a specific Ligand Layout is provided (v1.03+)
    if protocol.ligand_layout:
        l_layout = str(protocol.ligand_layout).lower()
        if "one ligand" in l_layout:
            chosen = None
            if config and config.ligand_choice_fn:
                chosen = config.ligand_choice_fn(
                    protocol.ligand or "Ligand 1",
                    protocol.ligand_2 or "Ligand 2",
                    plate_info)
            if chosen == "L2":
                return generate_col_mapping("one", 'L2', 'L1', block_count)
            if chosen is None and config and config.ligand_choice_fn:
                logger.warning(f"Ligand selection skipped for '{plate_info}'. Defaulting to Ligand 1 ({protocol.ligand}).")
            return generate_col_mapping("one", 'L1', 'L2', block_count)
        return generate_col_mapping(l_layout, 'L1', 'L2', block_count)

    # Fallback for older protocols (Inferred from Cell Layout)
    c_layout = str(protocol.line_layout).lower()

    if "one line" in c_layout:
        logger.error(f"Protocol '{protocol.file_name}' lists 2 ligands but uses 'One Line' "
              f"cell layout without specifying ligand layout.")

        # Ask user to choose a ligand layout, stored on protocol so subsequent plates reuse it
        layout_choice = None
        if config and config.ligand_layout_fn:
            layout_choice = config.ligand_layout_fn(
                protocol.ligand or "Ligand 1",
                protocol.ligand_2 or "Ligand 2",
                protocol.file_name)
        if layout_choice is None:
            logger.warning(f"Ligand layout selection skipped for protocol '{protocol.file_name}'. "
                           f"Defaulting to all Ligand 1 ({protocol.ligand}).")
            protocol.ligand_layout = "one ligand"
        else:
            logger.info(f"User selected ligand layout '{layout_choice}' for protocol '{protocol.file_name}'.")
            protocol.ligand_layout = layout_choice

        # Once protocol.ligand_layout is set, re-enter the function
        return get_ligand_map(protocol, block_count, config, plate_info)

    elif "half" in c_layout:
        # If Cells are Half, Ligands are Alternating
        return generate_col_mapping("alternating", 'L1', 'L2', block_count)

    elif "alternating" in c_layout:
        # If Cells are Alternating, Ligands are Half
        return generate_col_mapping("half", 'L1', 'L2', block_count)

    else:
        logger.warning(f"Unknown cell layout '{c_layout}'. Defaulting all to Ligand 1 {protocol.ligand}.")
        return generate_col_mapping("one", 'L1', 'L2', block_count)


def get_transfection_map(cell_layout_type: str, ligand_layout_type: str | None, t_ids: list[str], block_count: int = 4):
    """
    Defines the plate layout for blocks of transfection (taken from ID3 of the plate reader metadata).
    Depending on whether ligand layout is provided, transfection layout is decided on
    cell line layout or cell line and ligand layout.
    """
    logger.debug(f"--- Mapping Transfections ---")
    logger.debug(f"Raw ID3 List: {t_ids}")
    logger.debug(f"Layouts: Cell='{cell_layout_type}' / Ligand='{ligand_layout_type}'")
    # If no IDs, return empty
    if not t_ids:
        return ["N/A"] * block_count

    # Normalize missing ligand_layout to 'one ligand' (equivalent behavior)
    ligand_layout = str(ligand_layout_type).lower() if ligand_layout_type else "one ligand"
    cell_layout = str(cell_layout_type).lower()

    # Transfection ids are padded with N/A to be safe with indexing
    safe_ids = t_ids + ["N/A"] * block_count
    pair_count = block_count // 2  # expected t_ids length for paired patterns

    def _unique():
        """1 per block (e.g. 1, 2, 3, 4)"""
        return safe_ids[:block_count]

    def _single():
        """Single ID repeated (e.g. 1, 1, 1, 1)"""
        return [safe_ids[0]] * block_count

    def _seq_repeat():
        """Sequence repeat (e.g. 1,2,1,2). Caller guarantees len(t_ids) == pair_count."""
        return t_ids * 2

    def _elem_repeat():
        """Element repeat (e.g. 1,1,2,2). Caller guarantees len(t_ids) == pair_count."""
        return [x for x in t_ids for _ in range(2)]

    # Safety guards:  _single is correct only when exactly 1 transfection ID is provided, paired patterns
    # (_seq_repeat and _elem_repeat) are only valid when exactly pair_count IDs are provided
    # -> otherwise fall back to _unique
    def _single_safe():
        if len(t_ids) == 1: return _single()
        return _unique()

    def _seq_safe():
        if len(t_ids) == pair_count: return _seq_repeat()
        if len(t_ids) == 1:          return _single()
        return _unique()

    def _elem_safe():
        if len(t_ids) == pair_count: return _elem_repeat()
        if len(t_ids) == 1:          return _single()
        return _unique()

    # Decision matrix
    if "one line" in cell_layout:
        if "one ligand" in ligand_layout:    return _unique()
        if "half" in ligand_layout:          return _seq_safe()
        if "alternating" in ligand_layout:   return _elem_safe()

    elif "half" in cell_layout:
        if "one ligand" in ligand_layout:    return _seq_safe()
        if "half" in ligand_layout:          return _seq_safe()
        if "alternating" in ligand_layout:   return _single_safe()

    elif "alternating" in cell_layout:
        if "one ligand" in ligand_layout:    return _elem_safe()
        if "half" in ligand_layout:          return _single_safe()
        if "alternating" in ligand_layout:   return _elem_safe()

    # Fallback
    logger.warning(f"Unhandled layout combination: {cell_layout} + {ligand_layout}")
    return _unique()


# --- Ligand Info Inference from Master CSV --- #

def infer_ligand_info_from_master(df_master: pd.DataFrame, protocol: ProtocolData, block_count: int = 4):
    """
    Infers ligand layout and per-plate ligand choices from an existing master CSV.
    Used during enrichment to avoid re-prompting the user for information that was
    already decided when the master was originally created.

    Returns:
        inferred_layout: "half", "alternating", or "one ligand" (or None if not inferrable)
        file_ligand_choices: dict mapping file_name -> "L1" or "L2" (for one-ligand-per-plate cases)
    """
    if not protocol.ligand_2:
        return None, {}

    required_cols = {'File_Name', 'Well_ID', 'Ligand'}
    if not required_cols.issubset(df_master.columns):
        logger.warning("Master CSV missing columns needed for ligand inference. "
                       f"Required: {required_cols}")
        return None, {}

    # Only consider rows that have ligand info
    df = df_master[df_master['Ligand'].notna() & (df_master['Ligand'] != "")].copy()
    if df.empty:
        return None, {}

    # Extract column index from Well_ID
    df = df.copy()
    df['_col_idx'] = df['Well_ID'].str[1:].astype(int)

    file_ligand_choices = {}
    inferred_layout = None

    for fname in df['File_Name'].unique():
        file_rows = df[df['File_Name'] == fname]

        # Build col_idx -> ligand_name mapping (use first occurrence per column)
        col_ligand = {}
        for _, row in file_rows.drop_duplicates(subset=['_col_idx']).iterrows():
            col_ligand[row['_col_idx']] = row['Ligand']

        if not col_ligand:
            continue

        unique_ligands = set(col_ligand.values())

        if len(unique_ligands) == 1:
            # Single ligand on this plate — record which one
            lig_name = unique_ligands.pop() # pop returns arbitrary element of py set
            if lig_name == str(protocol.ligand_2):
                file_ligand_choices[fname] = "L2"
            else:
                file_ligand_choices[fname] = "L1"

            if inferred_layout is None:
                inferred_layout = "one ligand"

        elif len(unique_ligands) == 2 and inferred_layout is None:
            # Two ligands — determine the spatial layout
            try:
                inferred_layout = _detect_two_ligand_layout(col_ligand, protocol, block_count)
            except ValueError as e:
                logger.error(str(e))
                # If layout cannot be inferred, fallback to dialog prompt

    return inferred_layout, file_ligand_choices


def _detect_two_ligand_layout(col_ligand: dict, protocol: ProtocolData, block_count: int) -> str:
    """
    Detects whether columns follow a 'half' or 'alternating' ligand layout
    based on observed col->ligand assignments.
    Returns:
        "half" or "alternating" (or empty str if pattern is unclear).
    """
    # Map ligand names to L1/L2 labels
    col_labels = {}
    for col_idx, lig_name in col_ligand.items():
        if lig_name == str(protocol.ligand_2):
            col_labels[col_idx] = "L2"
        else:
            col_labels[col_idx] = "L1"

    # Check "half" pattern: cols 1-6 = L1, cols 7-12 = L2
    left_labels = {col_labels.get(c) for c in range(1, 7) if c in col_labels}
    right_labels = {col_labels.get(c) for c in range(7, 13) if c in col_labels}

    if left_labels == {"L1"} and right_labels == {"L2"}:
        return "half"

    # Check "alternating" pattern: blocks alternate between L1 and L2
    block_size = 12 // block_count
    block_labels = []
    for i in range(block_count):
        start = (i * block_size) + 1
        end = start + block_size
        labels_in_block = {col_labels.get(c) for c in range(start, end) if c in col_labels}
        if len(labels_in_block) == 1:
            block_labels.append(labels_in_block.pop())
        else:
            # Mixed block — doesn't fit either pattern cleanly
            raise ValueError(
                f"Cannot infer ligand layout from master: mixed ligands within block {i + 1} "
                f"(cols {start}-{end - 1}). Column assignments: {col_labels}")

    # Alternating
    if all(block_labels[i] != block_labels[i + 1] for i in range(len(block_labels) - 1)):
        return "alternating"

    raise ValueError(
        f"Cannot infer ligand layout from master: block pattern {block_labels} "
        f"does not match 'half' or 'alternating'. Column assignments: {col_labels}")