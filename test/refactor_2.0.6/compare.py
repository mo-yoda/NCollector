#!/usr/bin/env python3
"""
Step-0 comparer: compares master outputs with the recorded baseline, column by column.

    # re-run every baseline case with the CURRENT code and compare with test/data/baseline_2.0.5/
    uv run python test/refactor_2.0.6/compare.py --run [--no-legacy]

    # compare two files directly (csv or csv.gz)
    uv run python test/refactor_2.0.6/compare.py NEW.csv BASE.csv

Rows are aligned on (File_Name, Date, Well_ID, Time_(min)) when that key is unique,
otherwise by position. Numeric columns are compared with a float tolerance (NaN == NaN),
everything else as text. For each differing column, every differing well is printed
("File_Name:Well_ID") with a few example values. Missing / extra columns and row-count
differences are reported too. The provenance column Path (folder location) is not
compared. Exit code 0 = identical, 1 = differences.

TEMPORARY tooling for the 2.0.6 refactor (test/refactor_2.0.6/): tracked on the development
branch, deleted before the merge into main (plan: release cleanup). The data it reads and
writes (test/data/) is real lab data, local only and git-ignored.
"""
import argparse
import io
import os
import sys

import numpy as np
import pandas as pd

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

RTOL = 1e-9
ATOL = 1e-12
KEY = ["File_Name", "Date", "Well_ID", "Time_(min)"]
# Provenance only: the experiment folder path changes whenever the data folder moves
# (another machine, renamed folder), so it is not compared.
IGNORED_COLS = {"Path"}
MAX_EXAMPLES = 3


def _csv_roundtrip(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise an in-memory frame exactly like the baseline was stored (CSV)."""
    buf = io.StringIO()
    df.to_csv(buf, index=False)
    buf.seek(0)
    return pd.read_csv(buf, low_memory=False)


def _read(path: str) -> pd.DataFrame:
    return pd.read_csv(path, low_memory=False)


def _align(new: pd.DataFrame, base: pd.DataFrame):
    key = [k for k in KEY if k in new.columns and k in base.columns]
    if key and not new.duplicated(key).any() and not base.duplicated(key).any():
        n = new.set_index(key, drop=False)
        b = base.set_index(key, drop=False)
        only_new = n.index.difference(b.index)
        only_base = b.index.difference(n.index)
        common = b.index.intersection(n.index)
        return n.loc[common].reset_index(drop=True), b.loc[common].reset_index(drop=True), \
            len(only_new), len(only_base), "key " + "+".join(key)
    m = min(len(new), len(base))
    return new.iloc[:m].reset_index(drop=True), base.iloc[:m].reset_index(drop=True), \
        max(0, len(new) - m), max(0, len(base) - m), "position"


def _is_numeric(s: pd.Series) -> bool:
    return pd.api.types.is_numeric_dtype(s) and not pd.api.types.is_bool_dtype(s)


def _diff_mask(a: pd.Series, b: pd.Series) -> np.ndarray:
    if _is_numeric(a) and _is_numeric(b):
        av = a.to_numpy(dtype=float)
        bv = b.to_numpy(dtype=float)
        return ~np.isclose(av, bv, rtol=RTOL, atol=ATOL, equal_nan=True)
    an, bn = a.isna().to_numpy(), b.isna().to_numpy()
    eq = (a.astype(str).to_numpy() == b.astype(str).to_numpy())
    return ~((an & bn) | (~an & ~bn & eq))


def compare_frames(new: pd.DataFrame, base: pd.DataFrame, label: str = "", out=print) -> bool:
    """Prints every difference; returns True when the frames are identical."""
    ok = True
    title = f"=== {label} ===" if label else "==="
    out(title)
    missing = [c for c in base.columns if c not in new.columns]
    extra = [c for c in new.columns if c not in base.columns]
    if missing:
        ok = False
        out(f"  columns missing in new: {missing}")
    if extra:
        ok = False
        out(f"  columns only in new:    {extra}")
    if [c for c in new.columns if c in base.columns] != [c for c in base.columns if c in new.columns]:
        # Note only, not a failure: 2.0.5's ensure_master_csv_schema appends missing
        # legacy columns in set order, which varies between Python processes.
        out("  note: column order differs (values compared by name)")
    if len(new) != len(base):
        ok = False
        out(f"  row count: new {len(new)} vs baseline {len(base)}")

    n, b, n_only_new, n_only_base, how = _align(new, base)
    if n_only_new or n_only_base:
        ok = False
        out(f"  rows aligned by {how}: {n_only_new} only in new, {n_only_base} only in baseline")

    well_key = None
    if "File_Name" in b.columns and "Well_ID" in b.columns:
        well_key = b["File_Name"].astype(str) + ":" + b["Well_ID"].astype(str)

    for col in [c for c in base.columns if c in new.columns and c not in IGNORED_COLS]:
        mask = _diff_mask(n[col], b[col])
        if not mask.any():
            continue
        ok = False
        idx = np.flatnonzero(mask)
        wells = list(dict.fromkeys(well_key.iloc[idx])) if well_key is not None else []
        out(f"  DIFF {col}: {mask.sum()} rows, {len(wells)} wells")
        if wells:
            out("       wells: " + ", ".join(wells))
        for i in idx[:MAX_EXAMPLES]:
            where = well_key.iloc[i] if well_key is not None else f"row {i}"
            out(f"       e.g. {where}: new={n[col].iloc[i]!r} baseline={b[col].iloc[i]!r}")
    if ok:
        out("  identical")
    return ok


def run_all(baseline_dir: str, include_legacy: bool) -> bool:
    import make_baseline as mb  # same folder
    if not os.path.isdir(baseline_dir):
        # Baseline and fixtures are real lab data and not in the public repository.
        print(f"No baseline folder {baseline_dir} - run make_baseline.py on the unchanged "
              f"code first (needs the local fixtures in {mb.FIXTURES}).")
        return False
    os.chdir(REPO_ROOT)
    outputs = mb.build_outputs(include_legacy=include_legacy)
    all_ok = True
    for name in sorted(outputs):
        p = os.path.join(baseline_dir, name + ".csv.gz")
        if not os.path.exists(p):
            print(f"=== {name} ===\n  no baseline file {p}")
            all_ok = False
            continue
        all_ok &= compare_frames(_csv_roundtrip(outputs[name]), _read(p), label=name)
    for f in sorted(os.listdir(baseline_dir)):
        name = f[:-len(".csv.gz")] if f.endswith(".csv.gz") else None
        if name and name not in outputs and (include_legacy or not name.startswith("legacy__")):
            print(f"=== {name} ===\n  baseline exists but the case was not produced")
            all_ok = False
    print("\nRESULT: " + ("NO DIFFERENCES" if all_ok else "DIFFERENCES FOUND"))
    return all_ok


def main():
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="*", help="NEW BASE")
    ap.add_argument("--run", action="store_true", help="re-run all cases and compare with --baseline")
    ap.add_argument("--baseline", default=os.path.join(REPO_ROOT, "test", "data", "baseline_2.0.5"))
    ap.add_argument("--no-legacy", action="store_true")
    args = ap.parse_args()
    if args.run:
        ok = run_all(os.path.abspath(args.baseline), include_legacy=not args.no_legacy)
    elif len(args.files) == 2:
        ok = compare_frames(_read(args.files[0]), _read(args.files[1]),
                            label=f"{args.files[0]} vs {args.files[1]}")
    else:
        ap.error("give NEW BASE, or --run")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
