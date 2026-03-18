import logging
import pandas as pd

from models import ProtocolData

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

def get_ligand_map(protocol: ProtocolData, block_count: int = 4):
    """
    Determines which columns contain Ligand 1 and which contain Ligand 2.
    Returns dic of col idx and ligands.
    """

    # If no second ligand exists, everything is Ligand 1
    if not protocol.ligand_2:
        return generate_col_mapping("one", 'L1', 'L2', block_count)

    # If a specific Ligand Layout is provided (v1.03+)
    if protocol.ligand_layout:
        l_layout = str(protocol.ligand_layout).lower()
        if "one ligand" in l_layout:
            pass  # TODO: create pop up to ask which ligand was used for each plate
        return generate_col_mapping(l_layout, 'L1', 'L2', block_count)

    # Fallback for older protocols (Inferred from Cell Layout)
    c_layout = str(protocol.line_layout).lower()

    if "one line" in c_layout:
        logger.error(f"Protocol '{protocol.file_name}' lists 2 ligands but uses 'One Line' "
              f"cell layout without specifying ligand layout.")
        # TODO: create pop up to ask which ligand was used for each plate
        return generate_col_mapping("one", 'L1', 'L2', block_count)

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

    # Transfection ids are padded with N/A to be safe with indexing
    safe_ids = t_ids + ["N/A"] * block_count
    cell_layout = str(cell_layout_type).lower()

    def _unique():
        """1 per block (e.g. 1, 2, 3, 4)"""
        return safe_ids[:block_count]

    def _single():
        """Single ID repeated (e.g. 1, 1, 1, 1)"""
        return [safe_ids[0]] * block_count

    def _seq_repeat():
        """Repeat sequence (e.g. 1, 2, 1, 2)"""
        if len(t_ids) == 1:
            return _single()
        if len(t_ids) >= block_count:
            return _unique()
        if len(t_ids) == block_count / 2:
            return t_ids * 2

        # Fallback for uneven lengths
        return (safe_ids[:block_count // 2] * 2)[:block_count]

    def _elem_repeat():
        """Repeat elements (e.g. 1, 1, 2, 2)"""
        if len(t_ids) == 1:
            return _single()
        if len(t_ids) >= block_count:
            return _unique()
        if len(t_ids) == block_count / 2:
            res = []
            for x in t_ids:
                res.extend([x, x])
            return res

        # Fallback for uneven lengths
        res = []
        for x in safe_ids[:(block_count + 1) // 2]:
            res.extend([x, x])
        return res[:block_count]

    if not ligand_layout_type:
        if "one line" in cell_layout:
            return _unique()
        elif "half" in cell_layout:
            return _seq_repeat()
        elif "alternating" in cell_layout:
            return _elem_repeat()
        return _unique()  # Default fallback

    ligand_layout = str(ligand_layout_type).lower()
    if "one line" in cell_layout:
        if "one ligand" in ligand_layout:    return _unique()
        if "half" in ligand_layout:          return _seq_repeat()
        if "alternating" in ligand_layout:   return _elem_repeat()

    elif "half" in cell_layout:
        if "one ligand" in ligand_layout:    return _seq_repeat()
        if "half" in ligand_layout:
            # "1,2,1,2 or 1,2,3,4" -> depends on t_id count vs block_count
            return _unique() if len(t_ids) >= block_count else _seq_repeat()
        if "alternating" in ligand_layout:   return _single()

    elif "alternating" in cell_layout:
        if "one ligand" in ligand_layout:    return _elem_repeat()
        if "half" in ligand_layout:          return _single()
        if "alternating" in ligand_layout:
            # "1,1,2,2 or 1,2,3,4" -> depends on t_id count vs block_count
            return _unique() if len(t_ids) >= block_count else _elem_repeat()

        # Fallback
    logger.warning(f"Unhandled layout combination: {cell_layout} + {ligand_layout}")
    return _unique()
