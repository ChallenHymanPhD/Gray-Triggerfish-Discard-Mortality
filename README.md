# Gray Triggerfish Post-Release Mortality Analyses

This repository contains the R code used to generate posterior predictions, figures, and tables for an accepted manuscript examining post-release mortality in gray triggerfish (Balistes capriscus) as a function of capture depth and sea surface temperature (SST).
**NOTE**: Users must download R code, CJS Stan Models, and relevant CSV files to run on home file. Pre-run model outputs have been uploaded for convenience.

All analyses are based on posterior samples from a Cormack–Jolly–Seber (CJS) survival model with nonlinear covariate effects implemented using a modified softplus transformation.
---
## Model Overview

Post-release mortality is estimated using posterior simulations from a fitted CJS model:
*  Accounts for variable detection based on spatiotemporal effort and size-selectivity
*  Model selection based on LOO-IC
*  Depth and temperature effects are modeled using a softplus transformation
*  Posterior uncertainty is propagated by simulating across N = 1000 posterior draws
*  Mortality is reported as:
    - Median (50%)
    - 10th–90th percentile credible interval

Baseline survival scenarios (Table 4) explore sensitivity to assumed post-release survival

*  S = 1.00
*  S = 0.925
*  S = 0.85
---
## Figures
### Figure 1 — Conceptual diagram
Not produced here

### Figure 2 — Study site map
Not produced here

### Figure 3 — Observed data
Proportions of gray triggerfish recaptured as a function of: 
*  Panel A: depth (binned to 3-m intervals)
*  Panel B: SST (binned to 1-degree C intervals)

### Figure 4 — Conditional Effects
Shows predicted post-release mortality as a function of depth, conditional on three representative SST values (15, 22.5, 30 °C). Shaded ribbons represent 80% credible intervals.

### Figure 5 — Projected Mortality

*  Panel A: Annual posterior distributions (violin plots)
*  Panel B: Monthly median mortality with credible intervals

### Figure 6 — Long-Term Trends
Displays annual median mortality and credible intervals across the full time series.

---
## Tables

### Table 1 — Model equations
Conceptual - not produced here

### Table 2 — Model selection 
Results based on PSIS-LOO-IC conducted in R script

### Table 3 — Model results
Maximum a posteriori and 80% Bayesian credible intervals of posterior parameter distributions.

### Table 4 — Conditional Effects
Shows predicted post-release mortality as a function of depth, conditional on three representative SST values (15, 22.5, 30 °C). Maximum a posteriori and 80% Bayesian credible intervals are reported.

---
## Supplemental figures
### Figure S1 — Conceptual of size selectivity from Garner et al (2017)
Selectivity curve based on general exponential logistic selectivity for circle hooks

### Figure S2 — Conceptual diagram of softplus transformation
Conceptual diagram demonstrating the effect of depth on relative post-release survival using a modified softplus function and log link with differing effects sizes relating depth to post-release survival. 

### Figure S3 — Model diagnostics
Model fit based on comparisons of simulated recapture times to observed values
*  Panel A: Histogram of distribution of simulated recapture proportions with the observed recapture proportion indicated by the vertical line
*  Panel B: Posterior predictive calibration plot showing agreement between predicted and observed recapture probabilities. Predicted probabilities were grouped into evenly spaced bins, with the mean posterior predicted probability (x-axis) compared against the observed proportion of recaptures (y-axis) in each bin
*  Panel C: A comparison of the empirical distribution of observed recapture times to the distributions of 1,000 scans from the posterior predictive distribution.
*  Panel D: Empirical cumulative distribution function (ECDF) of observed recapture times compared with ECDFs simulated from posterior predictive draws.

### Figure S4 — Model diagnostics
Posterior predictive checks for best-fitting model. Predicted vs observed proportions of fish recaptured as a function of:
*  Panel A: Depth
*  Panel B: SST
*  Panel C: Fork length
