# N Collector

**N Collector** is a desktop application for documenting, processing, and
analysing plate reader resonance-energy-transfer (RET) experiments.

It uses the raw spreadsheets exported from the plate reader, together
with the experiment's protocol file, and turns them into:

1. a single **master CSV** that captures every well, timepoint, and piece of experimental
   metadata in one tidy, machine-readable table &rarr; suitable for downstream
   analysis via programmatic pipeline
2. **GraphPad Prism-ready export tables** (concentration-response, bar-graph, and heatmap
   layouts) &rarr; so users can directly create tidied plots

The app runs a BRET processing pipeline ((optional) labeling correction &rarr;
baseline correction &rarr; vehicle normalisation)
and provides an interactive GUI for excluding/restoring data, merging datasets, and
plotting concentration-response curves.

> **Download latest version: [v2.0.5](https://github.com/mo-yoda/NCollector/releases/tag/v2.0.5)**

---

This work is funded by the [Wellcome Trust](https://wellcome.org/research-funding/funding-portfolio/funded-grants/spatiotemporal-assessment-b-arrestin-centred)

---

## Infos for users
Currently, the setup was done for PHERAstar exports of concentration-response measurements. 
User-driven configuration options are planned updates (see below).

### Input

- **Experiment folders** named `YYMMDD_...` (the date prefix is required and it links the
  protocol to its measurements).
- Each folder contains:
  - **one protocol workbook** &rarr; an `.xlsx`/`.xlsm` file with a **`Protocol`** sheet
    (specifying cell lines, transfections, plasmids, ligand identities and concentrations, plate
    layout); and
  - **one or more measurement workbooks** &rarr; the plate reader exports, each with a
    **`Table All Cycles`** sheet (and an optional **`Protocol Information`** sheet).

Protocol and measurement files are matched automatically by the date parsed from the
folder name.

### Output

- **Master CSV**: one flat, long table (one row per well per timepoint) carrying all raw
  channels, every processed metric, and full metadata. This is the canonical,
  machine-readable record of the experiment, designed to feed analysis pipelines.
- **Prism-ready Excel exports**:cleanly shaped tables (concentration-response / bar-graph
  / heatmap) that paste directly into GraphPad Prism for plotting and curve fitting.

### Key features

- Reproducible BRET processing pipeline
- Reversible, documented exclusions
- Concentration-response plotting: check data quality to decide on exclusions or optionally export plot
- Import & re-use master CSVs:load a previously exported master and continue working
  (exclude, merge, plot, re-export)
- Dataset merging: combine several processed datasets (exported master CSVs and/or raw
  experiment folders) into one working master for pooled analysis

---
### Example data

A small example experiment is provided:

```
examples/
├── 260101_ExampleExperiment/     # one experiment folder (YYMMDD_... naming)
│   ├── protocol.xlsx             #  → 'Protocol' sheet
│   ├── measurement_01.xlsx       #  → 'Table All Cycles' sheet
│   └── measurement_02.xlsx
└── example_master.csv            # the master CSV produced from the folder above
```

Launch the app. Within the **Import & Export Data** tab: "Select folder containing reuslts of experiment" button to 
`examples/260101_ExampleExperiment/`, and "Load Files" — or import `examples/example_master.csv`
directly.
