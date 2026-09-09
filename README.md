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

---

This work is funded by the [Wellcome Trust](https://wellcome.org/research-funding/funding-portfolio/funded-grants/spatiotemporal-assessment-b-arrestin-centred)

---