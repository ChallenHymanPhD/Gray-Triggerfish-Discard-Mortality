# Gray Triggerfish Post-Release Mortality Analyses

This repository contains the R code and Stan models used to generate posterior predictions, figures, and tables for:

> Hyman AC, Ramsay C, Wilms S, Frazer TK. 2026. Depth and water temperature drive elevated post-release mortality of gray triggerfish (*Balistes capriscus*). *Frontiers in Marine Science* 13:1743605. https://doi.org/10.3389/fmars.2026.1743605

All analyses are based on posterior samples from a modified Cormack–Jolly–Seber (CJS) recapture model with nonlinear covariate effects implemented using a modified softplus transformation.

---

## Repository contents

| File | Description |
|---|---|
| `Hyman et al 2026 Gray Triggerfish Analysis.R` | Main script: data preparation, model fitting, model selection, posterior predictions, figures, and tables |
| `CJS_GTF_0.stan`, `CJS_GTF_1.stan`, `CJS_GTF_2.stan` | Candidate Stan models compared via PSIS-LOO |
| `GTF_Data_for_CJS.csv` | Tag-release and recapture data used to fit the models |
| `Rec_GTF.csv` | Private recreational fishing effort by time step and location |
| `FHO_GTF.csv` | For-hire fishing effort by time step and location |
| `GT Selectivity.csv` | Circle-hook size-selectivity inputs (after Garner et al. 2017) |
| `Longterm GTF.csv` | Long-term observer data used for retrospective mortality predictions |

## Requirements

- R (≥ [version]) and Stan via `rstan` (≥ [version])
- R packages: `tidyverse`, `lubridate`, `readxl`, `sf`, `sp`, `loo`, `bayesplot`, `cowplot`, `ggpubr`, `ggnewscale`, `ggdist`, `scales`, `xtable`

## How to run

1. Download all files into a single directory.
2. Set that directory as your working directory, or open it as an RStudio project.
3. Run `Hyman et al 2026 Gray Triggerfish Analysis.R` from top to bottom. Model fitting may take substantial time depending on hardware.

---

## Model overview

Post-release mortality is estimated using posterior simulations from the fitted recapture model, which:

* Accounts for variable re-encounter probability based on spatiotemporal fishing effort and size selectivity
* Models depth and temperature effects using a softplus transformation
* Selects among candidate models using PSIS-LOO
* Propagates posterior uncertainty by simulating across N = 1,000 posterior draws
* Reports mortality as the posterior median with 10th–90th percentile credible intervals

Baseline survival scenarios (Table 4) explore sensitivity to assumed post-release survival under best-case release conditions:

* S = 1.00
* S = 0.925
* S = 0.85

---

## Figures

**Figure 1: Conceptual diagram.** Not produced here.

**Figure 2: Study site map.** Not produced here.

**Figure 3: Observed data.** Proportions of gray triggerfish recaptured as a function of:
* Panel A: depth (binned to 3-m intervals)
* Panel B: SST (binned to 1 °C intervals)

**Figure 4: Conditional effects.** Predicted post-release mortality as a function of depth, conditional on three representative SST values (15, 22.5, and 30 °C). Shaded ribbons represent 80% credible intervals.

**Figure 5: Projected mortality.**
* Panel A: Annual posterior distributions (violin plots)
* Panel B: Monthly median mortality with credible intervals

**Figure 6: Long-term trends.** Annual median mortality and credible intervals across the full time series.

## Tables

**Table 1: Model equations.** Conceptual; not produced here.

**Table 2: Model selection.** PSIS-LOO results computed in the R script.

**Table 3: Model results.** Maximum a posteriori
