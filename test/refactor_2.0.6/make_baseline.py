#!/usr/bin/env python3
"""
Step 0 of the 2.0.6 implementation plan: record what the app produces, so later steps
can prove they change nothing except the decided differences.

Run from the repository root (paths in the output stay relative and therefore stable):

    uv run python test/refactor_2.0.6/make_baseline.py                 # -> test/data/baseline_2.0.5/
    uv run python test/refactor_2.0.6/make_baseline.py --out some/dir  # any other folder

The baseline in test/data/baseline_2.0.5/ must be produced with UNCHANGED 2.0.5 code.
test/refactor_2.0.6/compare.py re-runs `build_outputs()` with the current code and compares.

TEMPORARY tooling for the 2.0.6 refactor (test/refactor_2.0.6/): tracked on the development
branch, deleted before the merge into main (plan: release cleanup). The data it reads and
writes (test/data/) is real lab data, local only and git-ignored.

Outputs (gzip CSV, one per case):
  old_triplicate__nolab            master of the old_triplicate fixture, labeling off
  old_triplicate__lab              same, labeling correction on (quadruplicate layout)
  old_labeling__lab                the real labeling-control experiment, labeling on
  <each of the above>__excl_B2_H5  after excluding B2 and H5 in every plate file
                                   (recompute_master_after_exclusion)
  legacy__<csv name>               each legacy master after ensure_master_csv_schema

The ONLY place that knows how the current code is configured is `make_config()`, so a
later step that changes ProcessingConfig adapts that one function.
"""
import argparse
import os
import sys

import pandas as pd

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

FIXTURES = os.path.join("test", "data", "layout_2.0.6")
DEFAULT_OUT = os.path.join("test", "data", "baseline_2.0.5")

FOLDER_CASES = {
    # name: (fixture sub-path, labeling on?)
    "old_triplicate__nolab": ("old_triplicate", False),
    "old_triplicate__lab": ("old_triplicate", True),
    "old_labeling__lab": ("old_labeling", True),
}
EXCLUDED_WELLS = ["B2", "H5"]
EXCL_SUFFIX = "__excl_B2_H5"


def make_config(labeling: bool):
    """ProcessingConfig as the app builds it for a fresh load (no dialogs)."""
    from models import ProcessingConfig
    try:  # 2.0.5: explicit block layout on the config
        from models import build_plate_layout
        return ProcessingConfig(labeling_correction=labeling,
                                plate_layout=build_plate_layout(labeling))
    except ImportError:  # later steps: the layout comes from the presets
        return ProcessingConfig(labeling_correction=labeling)


def _experiment_folders(sub: str) -> list[str]:
    """The YYMMDD_* experiment folders below one fixture directory (relative paths)."""
    base = os.path.join(FIXTURES, sub)
    return sorted(os.path.join(base, d) for d in os.listdir(base)
                  if os.path.isdir(os.path.join(base, d)))


def build_folder_master(sub: str, labeling: bool, log_fn=None) -> pd.DataFrame:
    from master_builder import process_folder_to_master
    cfg = make_config(labeling)
    df = process_folder_to_master(_experiment_folders(sub), cfg,
                                  select_main_plasmids_fn=lambda exp: (exp, ""),
                                  log_fn=log_fn)
    return df.reset_index(drop=True)


def apply_fixed_exclusion(master: pd.DataFrame, labeling: bool) -> pd.DataFrame:
    from processing import recompute_master_after_exclusion
    df = master.copy()
    df.loc[df["Well_ID"].isin(EXCLUDED_WELLS), "Is_Excluded"] = True
    files = list(dict.fromkeys(df["File_Name"].tolist()))
    return recompute_master_after_exclusion(df, files, make_config(labeling)).reset_index(drop=True)


def migrate_legacy(path: str) -> pd.DataFrame:
    """Same steps as app._load_master_df (read -> ensure_master_csv_schema -> bool)."""
    from export import ensure_master_csv_schema
    from processing import coerce_bool
    df = pd.read_csv(path, low_memory=False)
    df, _modified, _fixed = ensure_master_csv_schema(df, log_fn=lambda *a, **k: None)
    if "Is_Excluded" in df.columns:
        df["Is_Excluded"] = df["Is_Excluded"].map(coerce_bool)
    return df.reset_index(drop=True)


def legacy_paths() -> dict[str, str]:
    base = os.path.join(FIXTURES, "legacy_masters")
    if not os.path.isdir(base):
        return {}
    out = {}
    for f in sorted(os.listdir(base)):
        if f.lower().endswith(".csv"):
            out["legacy__" + os.path.splitext(f)[0].replace(" ", "_")] = os.path.join(base, f)
    return out


def build_outputs(log_fn=None, include_legacy: bool = True) -> dict[str, pd.DataFrame]:
    """Every baseline case, computed with the code currently on sys.path."""
    out = {}
    for name, (sub, labeling) in FOLDER_CASES.items():
        if not os.path.isdir(os.path.join(FIXTURES, sub)):
            print(f"[make_baseline] fixture '{sub}' missing - case '{name}' skipped")
            continue
        master = build_folder_master(sub, labeling, log_fn=log_fn)
        out[name] = master
        out[name + EXCL_SUFFIX] = apply_fixed_exclusion(master, labeling)
    if include_legacy:
        for name, path in legacy_paths().items():
            out[name] = migrate_legacy(path)
    return out


def write_outputs(outputs: dict[str, pd.DataFrame], out_dir: str) -> None:
    os.makedirs(out_dir, exist_ok=True)
    for name, df in outputs.items():
        p = os.path.join(out_dir, name + ".csv.gz")
        df.to_csv(p, index=False)
        print(f"[make_baseline] {p}: {len(df)} rows x {df.shape[1]} cols")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=DEFAULT_OUT)
    ap.add_argument("--no-legacy", action="store_true", help="skip the (large) legacy masters")
    args = ap.parse_args()
    os.chdir(REPO_ROOT)
    write_outputs(build_outputs(include_legacy=not args.no_legacy), args.out)


if __name__ == "__main__":
    main()
