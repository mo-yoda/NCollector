# Example data

This folder contains a small, self-contained example so you can try the full N Collector
workflow.

The folder name **must** start with a `YYMMDD` date (e.g. `260101_...`). That date links
the protocol to its measurement files; a folder whose name doesn't start with a valid
date is skipped.

---

## How the protocol file is read

N Collector reads everything from the single **`Protocol`** sheet by searching for text
**anchors** and then reading a value at a fixed position relative to each anchor. The
example `protocol.xlsx` is annotated to make this visible:

- 🟩 **green** = the required **anchor** text (the label the app searches for)
- 🟨 **yellow** = the **value(s)** the app extracts for that anchor

Anchors are matched as a **case-sensitive substring**, and where an anchor appears more
than once, the app reads a **specific occurrence** (see the table). Keep the labels spelled
and ordered exactly as below.

### Required fields

| 🟩 Anchor (label on the sheet) | 🟨 Value read from | What it becomes | Format / notes                                                                                                                                                  |
|---|---|---|-----------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `date of measurement` | cell to the right | experiment date | `DD.MM.YY`. Must equal the folder's `YYMMDD` date, or the file is skipped as a mismatch.                                                                        |
| `Cell line layout` | cell directly below | cell-line plate layout | One of `One Line` / `Half` / `Alternating` (see *Dropdown_options* in `protocol.xlsx`).                                                                         |
| `Cell line` (2nd occurrence) | table below | list of cell lines used | The label must appear at least twice on the sheet; the app reads the table anchored on the second `Cell line`.                                                  |
| `DNA` | table below | transfection scheme → Main_Plasmids + per-condition plasmids | Columns = transfection numbers (`1`, `2`, …); rows = plasmid names; a value **> 0** means that plasmid is in that transfection.                                 |
| `Ligand dilution` (1st occurrence) | cell directly below | ligand 1 identity (name) | Mandatory — a missing ligand name will stop extraction.                                                                                                         |
| `final concentration in well (log(M))` (1st occurrence) | table below | ligand 1 concentrations | 8 rows → plate rows A–H, top to bottom.                                                                                                                         |
| `n =` | cell to the right | replicate count | Mandatory to ensure correct documentation, parsed and stored but **not** used for processing (n within app is derived from number of independent measurements). |

### Optional second ligand

| 🟩 Anchor | 🟨 Value read from | What it becomes | Notes |
|---|---|---|---|
| `Ligand dilution` (2nd occurrence) | cell directly below | ligand 2 identity | Presence of a value here switches the plate to two-ligand mode. |
| `final concentration in well (log(M))` (2nd occurrence) | table below | ligand 2 concentrations | Same A–H mapping as ligand 1. |
| `Ligand layout` | cell directly below | ligand plate layout | One of `One Ligand` / `Half` / `Alternating`. Used only with two ligands. |

---

## Details worth knowing when building a protocol
### The `DNA` transfection table

- Columns are the **transfection numbers**; rows are **plasmid names**; a cell **> 0**
  marks that plasmid as part of that transfection.
- Plasmids present in **every** transfection column become the shared **`Main_Plasmids`**
  (the BRET pair). The remaining, per-column plasmids become that condition's
  transfection label.
- Rows named `pcDNA3.1` are ignored (empty-vector filler).


**Row H is always treated as the vehicle control** (no ligand), regardless of what value
is entered there.

---

## Measurement file link to the protocol
The protocol supplies the **layout rules** and the **identity dictionaries**. The actual
per-plate **cell-line names** and **transfection numbers** come from the measurement file's
`Table All Cycles` sheet (`ID2` = cell line, `ID3` = transfections). N Collector combines
the two to label every well:

- **Cell line:** name from `ID2`, placement from the protocol's `Cell line layout`.
- **Transfection / condition:** number from `ID3`, translated to plasmids by the protocol's
  `DNA` table.
- **Ligand identity & concentration:** entirely from the protocol (`Ligand dilution` +
  `final concentration in well (log(M))`), placed by `Ligand layout`.
