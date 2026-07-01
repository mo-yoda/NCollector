"""
Builders for the per-column metadata index that drives the Tab-2 exclusion dropdowns and the Tab-1 summary table.

GUI-agnostic: pure record builders from either live PrResult objects (initial load)
or a flat master DataFrame (CSV/import mode, after every exclusion, and the merge
preview). The GUI tail (storing the DataFrame, refreshing dropdowns) stays in app.py.
"""
import logging

from models import MAIN_ONLY_CONDITION
from processing import parse_date_series, coerce_bool

logger = logging.getLogger("NCollector")


def build_records_from_objects(experiment):
    """Object-path index records (initial load). Reads live column_metadata.

    Emits one record per non-Empty, not-fully-excluded plate column carrying
    File_Name / Date / Cell_Line / Condition / Ligand / Replicate / Main_Plasmids.
    """
    records = []
    rows_str = "ABCDEFGH"

    for folder in experiment:
        # Match the Main_Plasmids string exactly as the compile path writes it onto master_df
        proto = getattr(folder, "protocol", None)
        main_plasmids = (" + ".join(proto.main_plasmids)
                         if proto and proto.main_plasmids else "Unknown")
        for result in folder.results:
            if not result.column_metadata: continue
            if result.is_excluded: continue

            for col_idx, meta in result.column_metadata.items():
                # Filter out empty cols
                if meta.condition_name is None or "Empty" in meta.condition_name: continue
                # Skip columns whose every well is excluded (N count drops). At
                # initial load excluded_wells is empty, so nothing is skipped here.
                all_wells_excluded = all(
                    f"{r}{col_idx}" in result.excluded_wells for r in rows_str)
                if all_wells_excluded: continue

                # Main plasmid-only condition (empty plasmids list -> condition_name "") gets
                # the same sentinel as compile_master_dataframe / the CSV-import path, so the
                # initial-load dropdowns match the master_df.
                cond_label = (meta.condition_name
                              if str(meta.condition_name).strip() not in ("", "nan", "None")
                              else MAIN_ONLY_CONDITION)

                records.append({
                    "File_Name": result.file_name,
                    "Date": result.measurement_date.strftime('%d.%m.%y'),  # dropdown string
                    "Cell_Line": meta.cell_line,
                    "Condition": cond_label,
                    "Ligand": meta.ligand_identity,
                    "Replicate": meta.replicate,
                    "Main_Plasmids": main_plasmids,
                })

    return records


def build_records_from_master(master_df, source_str=None, include_replicate=False,
                              context="index"):
    """Master-path index records from a master-shaped df: one record per
    (File_Name, column) that is not Empty/Unknown and not fully excluded, carrying
    File_Name / Date / Cell_Line / Condition(=Transfection) / Ligand / Main_Plasmids.

    include_replicate adds the Replicate field (Tab-2 dropdowns); source_str, when
    given, tags each record Source=source_str (merge preview). Read-only — never
    mutates the df. `context` is only a label passed to parse_date_series for logging.
    """
    records = []
    if master_df is None or master_df.empty:
        return records
    work = master_df.copy()
    # Dropdown date string (consistent with the object branch '%d.%m.%y').
    try:
        work['_DateStr'] = parse_date_series(
            work['Date'], context=context).dt.strftime('%d.%m.%y')
    except Exception:
        work['_DateStr'] = work['Date'].astype(str)
    work['_ColIdx'] = work['Well_ID'].astype(str).str[1:]
    work['_Excl'] = (work['Is_Excluded'].map(coerce_bool)
                     if 'Is_Excluded' in work.columns else False)

    # Drop empty/unknown columns (defensive — compile already removes them).
    work = work[~work['Transfection'].astype(str).str.contains("Empty", na=False)]
    work = work[~work['Cell_Line'].astype(str).str.startswith("Unknown", na=False)]

    for (fname, _col), g in work.groupby(['File_Name', '_ColIdx'], sort=False):
        # Skip a column whose every well is currently excluded (N count drops).
        if g['_Excl'].all():
            continue
        first = g.iloc[0]
        rec = {
            "File_Name": fname,
            "Date": first['_DateStr'],
            "Cell_Line": first['Cell_Line'],
            "Condition": first['Transfection'],
            "Ligand": first['Ligand'],
        }
        if include_replicate:
            rec["Replicate"] = str(first['Replicate'])
        rec["Main_Plasmids"] = (str(first['Main_Plasmids'])
                                if 'Main_Plasmids' in g.columns else "Unknown")
        if source_str is not None:
            rec["Source"] = source_str
        records.append(rec)
    return records
