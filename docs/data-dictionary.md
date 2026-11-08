# Data dictionary

## `data/raw/logger/`

`VK_Control_C.csv` (Tank 3, control) and `VK_Treatment_T.csv` (Tank 4, heated)
are unmodified Onset HOBO UA-002-64 Pendant exports at 15-minute intervals,
deployed 16:00 on 11 March 2025 and downloaded 26 March.

| Column | Meaning |
|---|---|
| `#` | sequential record number assigned by the logger |
| `Date Time, GMT+10:00` | local Guam time, `MM/DD/YY HH:MM:SS AM/PM` |
| `Temp, °C (…)` | seawater temperature |
| `Intensity, lum/ft² (…)` | illuminance at the logger |
| `Coupler Attached`, `Host Connected`, `Stopped`, `End Of File` | logger events |

Temperature is in **degrees Celsius**, as the CSV column header states. The
`Series:` line of `VK_Control_Details.txt` says °F; that label is wrong and the
values themselves (28–32) confirm it — the series was converted before export.
The two `*_Details.txt` files are the logger metadata panels, recording serial
numbers 21086852 (control) and 21086851 (treatment).

Readings before 16:45 on 11 March were taken while the loggers were still in
air and are excluded by every script.

Illuminance is a **photometric** quantity weighted to human visual
sensitivity. It is not photosynthetically active radiation and cannot be
compared with PAR values reported elsewhere. It is used here only to separate
night from day, and to document that the two loggers were not equally lit.

## `data/raw/pam/`

`VK-1.csv`, `VK-2.csv`, `VK-3.csv` are unmodified Diving-PAM II WinControl
report files, semicolon-delimited, covering the three measurement sessions.
They hold `F`, `Fm'`, `PAR`, `Y(II)` and `ETR` per reading with a timestamp,
but **no species or tank labels** — those were recorded separately and, for
most nights, reconstructed (see the README). A single blank row marks a change
of species; two consecutive blanks indicate failed readings, not a boundary.

`PAR` reads zero throughout, confirming that all measurements were made in
full darkness.

## `data/raw/qpcr/`

`qPCR_Cq_values.csv`, long format, one row per technical replicate.

| Column | Meaning |
|---|---|
| `Sample_ID` | fragment identifier, e.g. `AC1` = *Acropora*, control, fragment 1 |
| `Species` | `Acropora`, `Pavona` or `Sclerophytum` |
| `Treatment` | `Control` or `Heat` |
| `Target` | `Cni_ACTB` (reference), `Endozoicomonas`, `Cladocopium`, `Durusdinium` |
| `Tech_Rep` | technical replicate, 1 or 2 |
| `Cq` | quantification cycle; blank where no amplification occurred |

Biological replicates: *Acropora* n = 4, *Pavona* n = 4, *Sclerophytum* n = 2
per group — two *Sclerophytum* colonies died before sampling. Samples `AC2` and
`AC3` show no *Endozoicomonas* amplification, so *Acropora* control is n = 2
for that target alone. *Acropora* *Endozoicomonas* (Cq 33–35) and *Durusdinium*
(Cq 32–33) sit at or beyond reliable quantification;
`03_qpcr_analysis.R` flags values above Cq 32 and treats values above Cq 36 as
non-detects. No standard curves were run, so amplification efficiency is
assumed equal across targets — a limitation stated in the manuscript.

## `data/processed/`

Four files, one per measurement block, with identical columns. These are the
PAM readings after species, tank and phase were assigned.

| File | Date | Readings |
|---|---|---|
| `PAM_BASELINE_10Mar.csv` | 10 March | 21 |
| `PAM_BASELINE.csv` | 11 March | 69 |
| `PAM_HEAT.csv` | 12–16 March | 326 |
| `PAM_RECOVERY.csv` | 17 March | 41 |

| Column | Meaning |
|---|---|
| `Date` | `YYYY-MM-DD` |
| `Time` | reading time from the instrument; `NA` on 10 March, transcribed by hand |
| `Phase` | as originally labelled; the scripts re-derive phase from date |
| `Elevation_C` | night-time mean elevation of the heated tank, from the logger record |
| `Species` | `Acropora`, `Pavona` or `Sclerophytum` |
| `Tank` | `Control` or `Heat` |
| `Tank_basis` | **how this reading's tank was assigned** — see below |
| `Reading_no` | position of this reading within its night's species block |
| `n_species` | number of readings in that species block |
| `split_at` | reading number at which the block was divided between the two tanks |
| `Fo` | minimum fluorescence; `NA` on 10 March |
| `Fm` | maximum fluorescence; `NA` on 10 March |
| `FvFm` | dark-adapted maximum quantum yield, (Fm − Fo) / Fm |

`Tank_basis` takes three kinds of value:

- `recorded by hand` — tank written down at the time of measurement. 10 March only.
- `your marker HH:MM:SS` — the operator's own recollection, supplied afterwards,
  of the timestamp at which measurement moved from one tank to the next.
- `equal split` — no boundary information; the species block was divided evenly
  in measurement order.

Readings assigned by `equal split` are the weakest link in the dataset. The
`"10only"` sensitivity analysis exists to test whether they matter. They do not
change the conclusion.

## `results/`

Files suffixed `_PUBLISHED` are the outputs exactly as they stand in the
manuscript, kept so that a re-run can be diffed against them. Everything else
in `results/` is regenerated by the scripts.

| File | Corresponds to |
|---|---|
| `tables/TableS1_nightly_temperature.csv` | Supplementary Table S1 |
| `tables/Table_PAM_descriptive.csv` | Supplementary Table S2 (per night; the manuscript aggregates to phase) |
| `tables/Table_PAM_effect_sizes.csv` | Supplementary Table S5 |
| `tables/Table1_qPCR_foldchange.csv` | Table 1 |
| `tables/dose_by_night.csv` | cumulative excess heat at 19:00 each night |
| `reports/PAM_statistics_report.txt` | Supplementary Tables S3, S4 and the Results statistics |
| `reports/qPCR_statistics_report.txt` | the qPCR statistics in Results |
| `reports/logger_report.txt` | the thermal values in Results |
