# Thermal stress experiment on three Guam corals, March 2025

Data and analysis code for the manuscript:

> *Endozoicomonas and Symbiodiniaceae abundance during short-term experimental warming in two hard corals and a soft coral* Shimpi GG, Kim V, Bentlage B. Submitted to Access Microbiology.

Two hard corals, *Acropora* cf. *manni* and *Pavona decussata*, and the soft coral
*Sclerophytum* sp. were warmed together in a single outdoor flow-through tank
at the University of Guam Marine Laboratory, alongside an unheated control
tank. Dark-adapted photosynthetic efficiency was measured nightly, seawater
temperature and illuminance were logged continuously, and symbiont abundance
was quantified by qPCR at the end of warming.

---

## Read this before using the data

**One tank per treatment.** All fragments of a given treatment shared a single
tank, so they are not independent with respect to temperature.
Control-versus-heated differences are reported throughout as **effect sizes
only** and are never tested inferentially. The comparison that *is* valid is
**among species within the heated tank**: all three shared the same water, so
tank identity is held constant and a difference between species cannot be a
tank artefact. Every script enforces this distinction.

**The experimental phases are not the ones originally planned.** The logger
record places heater switch-on at approximately 19:30 on 11 March — the two
tanks were identical at 19:15 and differed by 2.0 °C ninety minutes later.
That night's PAM session ran 19:04–19:17, so 11 March remains a valid
baseline. Phases are therefore baseline 10–11 March, warming 12–16 March,
recovery 17 March. `01_logger_analysis.R` prints the switch-on window.

**Tank assignment for 11–17 March was reconstructed, not recorded.** The
Diving-PAM II exports a continuous series with no tank or species labels. Only
10 March has tank labels written down at the time. For the other nights,
species boundaries were identified from the single blank rows the operator
inserted when moving between species, and readings were then split evenly
between the two tanks within each species block, in measurement order. The
`Tank_basis` column in `data/processed/` records how **every individual
reading** was assigned, so a reader can see which assignments are recorded
fact and which are reconstruction. Because that reconstruction could in
principle bias the result, the whole analysis is repeated using only the
hand-labelled night as baseline; see *Sensitivity analysis* below. The
conclusions do not change.

---

## Layout

```
data/raw/logger/   HOBO UA-002-64 exports, control (Tank 3) and treatment (Tank 4)
data/raw/pam/      Diving-PAM II WinControl exports, three sessions, unlabelled
data/raw/qpcr/     Cq values, long format, two technical replicates per reaction
data/processed/    PAM readings with species, tank and phase assigned
scripts/           00_common.R and three numbered analysis scripts
results/tables/    summary tables as published
results/reports/   full statistical output
results/figures/   written when the scripts are run (not tracked)
docs/              data dictionary
```

## Running the analysis

Requires R ≥ 4.1 (the scripts use the native `|>` pipe) and the packages
`readr`, `dplyr`, `tidyr`, `ggplot2`, `ggbeeswarm`, `purrr` and `scales`.
Anything missing is installed on first run.

From the `scripts/` directory, in order:

```
Rscript 01_logger_analysis.R
Rscript 02_pam_analysis.R
Rscript 03_qpcr_analysis.R
```

`00_common.R` holds the shared settings and helpers and is sourced by the
other three; it is not run on its own. The repository root is resolved
automatically when the scripts are run from `scripts/`; otherwise set
`REPO <- "/path/to/guam-coral-thermal-stress"` before sourcing.

| Script | Reproduces |
|---|---|
| `01_logger_analysis.R` | Supplementary Table S1, the thermal values in Results, Figure 1, Supplementary Figure S3 |
| `02_pam_analysis.R` | all Fv/Fm statistics, Supplementary Tables S2–S6, Figures 2–4, Supplementary Figure S1 |
| `03_qpcr_analysis.R` | Table 1, all qPCR statistics, Figure 5 |

Supplementary Figure S2 (logger illuminance) is not produced here: it is a
plot of the `Light` column of the two raw logger files and carries no
analysis.

Figures are written as PNG, SVG and PDF. The SVG is written with `svglite`
when that package is installed, because svglite emits real `<text>`
elements and the labels stay editable in Inkscape or Illustrator. Without
it the scripts fall back to R's built-in `grDevices::svg()`, which is
cairo-backed and converts every character to an outline path: the figure
looks the same but the text becomes groups of letter-shaped paths that
cannot be retyped. If `svglite` fails to load because of a version clash
with `systemfonts`, reinstall the pair from source together:

```r
install.packages(c("systemfonts", "svglite"), type = "source")
```

Species names are drawn as plotmath (`sp_expr` and `sp_labeller` in
`00_common.R`) so that *Acropora* cf. *manni* and *Sclerophytum* sp. keep
"cf." and "sp." upright while the names themselves are italic. Anywhere
those labellers are used, do not also set `face = "italic"` on the axis, or
the connecting terms are slanted too.

Figure 1 and Supplementary Figure S3 are the same thermal record at two
spans, drawn by one function in `01_logger_analysis.R` so that they cannot
drift apart: Figure 1 ends on 19 March, the experiment plus one day of
convergence, and Figure S3 runs to the end of logging on 27 March. Both put
cumulative excess heat on a right-hand axis. That axis is anchored to the
observed temperature range, which places dose zero on the panel floor, below
both tanks' starting temperatures; it is axis arithmetic and not a
measurement. `PANEL_LAST` in `00_common.R` sets where Figure 1 ends.

Before saving either figure the script tests that the shaded phase bands
actually sit where the heat is — no difference between the tanks before the
band starts, at least 1.5 °C inside it, none again after the experiment — and
stops rather than write a figure whose shading is in the wrong place. The
bands are anchored to heater switch-on at 19:30 on 11 March, which
`results/reports/logger_report.txt` establishes from the paired 15-minute
records.

### Each script checks itself

Every script ends by comparing its own output against the values published in
the manuscript and **stops with an error if they disagree**. The targets are
listed in `PUB` at the top of `00_common.R`. If a script halts there, the
output no longer reproduces the paper and should not be used. The outputs as
published are kept alongside, suffixed `_PUBLISHED`, so a re-run can be
diffed against them directly.

### Sensitivity analysis

Supplementary Table S6 compares the principal between-species result under two
baseline definitions. To reproduce the alternative, set

```r
BASELINE_NIGHTS <- "10only"
```

in `00_common.R` and re-run `02_pam_analysis.R`. That restricts normalisation
to 10 March, the only night whose tank assignment was recorded at the time.
The contrasts come out larger, because that night supplies only three to four
readings per species per tank and the baseline means are correspondingly less
precise, but the conclusion is unchanged. The 10–11 March baseline is used
throughout the manuscript. Note that the verification checks are calibrated to
the primary analysis and will fail under `"10only"`; that is expected.

## Statistical approach

Fv/Fm readings below 0.15 or above 0.80 are removed as measurement artefacts,
followed by an interquartile-range screen within each date, species and tank.
Remaining values are normalised to each species-by-tank group's own baseline,
because species differ substantially in absolute Fv/Fm
(F<sub>2,83</sub> = 59.27, p < 0.001) for reasons unrelated to condition, while
the two fragment sets of a given species do not differ (all p > 0.12).

The headline test is a one-way ANOVA on percentage of baseline across species
within the heated tank during warming, with Tukey HSD contrasts and a
Kruskal–Wallis test as non-parametric confirmation.

Heat exposure is expressed as **cumulative excess heat in °C-days**, obtained
by integrating the temperature difference between the two tanks over the
paired 15-minute logger records. This is referenced to the concurrent control,
not to a regional climatological maximum monthly mean, so it is **not** Degree
Heating Weeks and must not be compared with published DHW values. Effect sizes
are regressed on dose across 11–16 March: the pre-warming night is retained
because it anchors the relationship at zero exposure, and the recovery night
is excluded because cumulative exposure keeps rising once heating has stopped.

qPCR follows Livak's 2<sup>−ΔΔCq</sup> method with `Cni_ACTB` as reference and
the species-matched control as calibrator. All tests are run on ΔΔCq, which is
on a log<sub>2</sub> scale, never on fold change, which is log-normally
distributed. *Cladocopium* serves as a specificity control and is expected to
show no effect.

Fragments were not individually identified between nights, so readings are
treated as independent within each night and no repeated-measures structure is
fitted. This is stated in the manuscript as a limitation.

## Still to be added

Two raw datasets are not yet in this repository and will be deposited before
publication: the **AriaMx instrument export** (the Cq table derived from it is
present, but not the run files or melt curves) and the original **`.hobo`
logger files** (the CSV exports are present). See `docs/TODO.md`.

## Licence

Code is released under the MIT Licence (`LICENSE-CODE`); data under CC BY 4.0
(`LICENSE-DATA`).

## Contact

Gaurav G. Shimpi — nature.gaurav@gmail.com; marinerealm@gmail.com
