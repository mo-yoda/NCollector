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
