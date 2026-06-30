"""
crc_plot.py — GUI-agnostic data extraction + 4-parameter sigmoidal fit for the
concentration-response plot in Exclude tab.

A "condition" is one (Ligand, Cell_Line, Transfection) combination. For each condition
the plot shows, per measurement Date (= one biological replicate, N), the response
(AUC_Mean, the vehicle-normalised AUC averaged over technical replicates) against the
ligand concentration (Ligand_Conc, log10 M), plus the across-date mean and a 4PL fit to
that mean.

Everything here is a pure function of the master DataFrame (no Tkinter, no matplotlib).
Only `fit_four_pl` needs scipy, imported lazily so listing/extraction work without it.
Because AUC_Mean is recomputed on every exclude/restore, reading it here automatically
reflects the current exclusion state (fully-excluded concentrations are NaN -> dropped).
"""

import logging
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from processing import parse_date_series

logger = logging.getLogger("NCollector")

# Columns the plot needs from the master.
_REQUIRED = {'Ligand', 'Cell_Line', 'Transfection', 'Ligand_Conc', 'AUC_Mean',
             'File_Name', 'Date'}


@dataclass
class CRCReplicate:
    """One biological replicate (one measurement Date) of a condition."""
    date_label: str
    concs: list           # log10(M) concentrations, ascending
    values: list          # AUC_Mean at each concentration


@dataclass
class CRCData:
    """All series needed to draw one condition's concentration-response plot."""
    ligand: str
    cell_line: str
    transfection: str
    main_plasmids: str
    replicates: list = field(default_factory=list)   # list[CRCReplicate], date-ordered
    concs: list = field(default_factory=list)    # distinct concentrations x
    mean_values: list = field(default_factory=list)   # across-date mean: y
    n: int = 0                                        # number of biological replicates


def _has_required(df) -> bool:
    return df is not None and not df.empty and _REQUIRED.issubset(df.columns)


def list_crc_conditions(master_df) -> list:
    """Ordered list of plottable conditions in the master.

    Each entry is {ligand, cell_line, transfection, main_plasmids}. Only conditions with at
    least one real dose concentration carrying a value are included; Empty/Unknown rows are
    skipped. Ordered by (ligand, cell_line, transfection) for stable < / > navigation.
    """
    if not _has_required(master_df):
        return []
    df = master_df
    conc = pd.to_numeric(df['Ligand_Conc'], errors='coerce')
    auc = pd.to_numeric(df['AUC_Mean'], errors='coerce')
    mask = (conc.notna() & auc.notna()
            & ~df['Transfection'].astype(str).str.contains("Empty", na=False)
            & ~df['Cell_Line'].astype(str).str.startswith("Unknown", na=False))
    sub = df[mask]
    if sub.empty:
        return []

    out = []
    for (lig, cell, trans), g in sub.groupby(['Ligand', 'Cell_Line', 'Transfection'],
                                             sort=True):
        mp = ""
        if 'Main_Plasmids' in g.columns:
            mps = [str(v) for v in pd.unique(g['Main_Plasmids'].dropna())
                   if str(v) not in ("", "nan")]
            mp = " / ".join(mps)
        out.append({"ligand": str(lig), "cell_line": str(cell),
                    "transfection": str(trans), "main_plasmids": mp})
    return out


def extract_condition_data(master_df, ligand, cell_line, transfection) -> "CRCData | None":
    """Build the per-replicate + mean series for one condition, or None if no rows match.

    One biological replicate = one measurement Date (files sharing a Date are averaged into
    that date's value per concentration). The across-date mean is the mean over dates at each
    concentration. Vehicle rows (no concentration) are dropped; NaN AUC_Mean (fully-excluded
    concentrations) is dropped, so the result reflects the current exclusion state.
    """
    if not _has_required(master_df):
        return None
    df = master_df
    m = ((df['Ligand'].astype(str) == str(ligand))
         & (df['Cell_Line'].astype(str) == str(cell_line))
         & (df['Transfection'].astype(str) == str(transfection)))
    sub = df[m].copy()
    if sub.empty:
        return None

    sub['_conc'] = pd.to_numeric(sub['Ligand_Conc'], errors='coerce')
    sub['_auc'] = pd.to_numeric(sub['AUC_Mean'], errors='coerce')
    sub = sub[sub['_conc'].notna()]                      # drop vehicle / no-dose rows
    if sub.empty:
        return None

    main_plasmids = ""
    if 'Main_Plasmids' in sub.columns:
        mps = [str(v) for v in pd.unique(sub['Main_Plasmids'].dropna())
               if str(v) not in ("", "nan")]
        main_plasmids = " / ".join(mps)

    # One value per (File, concentration): AUC_Mean is constant across a file's wells/times
    # for a given concentration, but average defensively and keep the file's Date.
    per_file = sub.groupby(['File_Name', '_conc'], as_index=False).agg(
        _auc=('_auc', 'mean'), Date=('Date', 'first'))
    per_file = per_file[per_file['_auc'].notna()]
    if per_file.empty:
        return CRCData(str(ligand), str(cell_line), str(transfection), main_plasmids)

    # Date string -> parsed datetime (for ordering + label); group replicates by Date.
    per_file['_datestr'] = per_file['Date'].astype(str)
    parsed = parse_date_series(per_file['Date'], context="crc_plot")
    datestr_to_dt = dict(zip(per_file['_datestr'], parsed))

    # Per (date, concentration): mean over the date's files.
    grp = per_file.groupby(['_datestr', '_conc'], as_index=False)['_auc'].mean()

    def _sort_key(s):
        dt = datestr_to_dt.get(s)
        return (1, str(s)) if dt is None or pd.isna(dt) else (0, dt)

    def _label(s):
        dt = datestr_to_dt.get(s)
        if dt is not None and not pd.isna(dt):
            return pd.Timestamp(dt).strftime('%d.%m.%y')
        return str(s)

    date_strs = sorted(grp['_datestr'].unique(), key=_sort_key)

    replicates = []
    for ds in date_strs:
        g = grp[grp['_datestr'] == ds].sort_values('_conc')
        replicates.append(CRCReplicate(
            date_label=_label(ds),
            concs=[float(x) for x in g['_conc']],
            values=[float(y) for y in g['_auc']]))

    # Across-date mean per concentration,  grouping by '_conc' keeps each distinct
    # concentration as its own x
    mean_grp = grp.groupby('_conc', as_index=False)['_auc'].mean().sort_values('_conc')
    return CRCData(
        ligand=str(ligand), cell_line=str(cell_line), transfection=str(transfection),
        main_plasmids=main_plasmids,
        replicates=replicates,
        concs=[float(x) for x in mean_grp['_conc']],
        mean_values=[float(y) for y in mean_grp['_auc']],
        n=len(date_strs))


def four_pl(x, bottom, top, logec50, hill):
    """4-parameter logistic (log-dose response):
        y = bottom + (top - bottom) / (1 + 10^((logEC50 - x) * hill))
    x is log10(concentration)."""
    return bottom + (top - bottom) / (1.0 + 10.0 ** ((logec50 - x) * hill))


def fit_four_pl(x, y):
    """Fit the 4PL to (x = log10 conc, y = response).

    Returns (popt, x_smooth, y_smooth) for drawing, or None if the fit cannot be computed
    (fewer than 4 distinct concentrations, or no convergence). scipy is imported lazily.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(np.unique(x)) < 4:
        return None
    try:
        from scipy.optimize import curve_fit
    except Exception as e:
        logger.warning(f"scipy unavailable, skipping 4PL fit: {e}")
        return None

    ymin, ymax = float(np.min(y)), float(np.max(y))
    p0 = [ymin, ymax, float(np.median(x)), 1.0]
    try:
        popt, _ = curve_fit(four_pl, x, y, p0=p0, maxfev=10000)
    except Exception as e:
        logger.debug(f"4PL fit did not converge: {e}")
        return None

    xs = np.linspace(float(np.min(x)), float(np.max(x)), 200)
    ys = four_pl(xs, *popt)
    return popt, xs, ys
