# Drosophila-climbing-assay

Reproducible R code and data accompanying the manuscript:

> **An individual-tube climbing assay with transparent analysis and quality control for Drosophila locomotor phenotyping**  
> *(Submitted to Animal Behaviour)*

## Overview

This repository contains the data and R scripts used to analyse an individual-tube negative geotaxis (climbing) assay in *Drosophila melanogaster*.

The repository reproduces the principal analyses presented in the manuscript, including:

- Longitudinal age-profile analysis and Day 21 readout-window selection
- Day 21 vehicle-control reproducibility analysis
- Salbutamol worked pharmacological application
- Round-level weighted AUC sensitivity analysis

The scripts generate the statistical models, estimated marginal means, contrasts, diagnostic summaries, supplementary tables and **ggplot2** objects underlying the manuscript figures. Final figure assembly and SVG export were performed manually during manuscript preparation.

---

## Repository structure

```
Drosophila-climbing-assay/
│
├── data/
│   ├── longitudinal.xlsx
│   ├── day21_vehicle_only.xlsx
│   └── salbutamol_application.xlsx
│
├── scripts/
│   ├── 01_longitudinal_analysis.R
│   ├── 02_vehicle_control_analysis.R
│   └── 03_salbutamol_application.R
│
├── output/
│   ├── figures/
│   └── tables/
│
├── README.md
├── LICENSE
└── .gitignore
```

---

## Requirements

The analyses were performed in **R (version 4.x)** using the following packages:

- tidyverse
- readxl
- lme4
- emmeans
- DHARMa
- performance
- patchwork
- scales
- writexl

---

## Running the analyses

Each script is independent and can be run separately.

1. `01_longitudinal_analysis.R`
2. `02_vehicle_control_analysis.R`
3. `03_salbutamol_application.R`

The scripts assume that the working directory is the repository root.

---

## Outputs

The scripts generate:

- Generalised linear mixed models (GLMMs)
- Estimated marginal means
- Treatment and genotype contrasts
- Model diagnostics
- Supplementary tables
- **ggplot2** objects corresponding to the manuscript figures

Supplementary tables are written to the `output/tables/` directory.

---

## Data availability

The datasets used in the manuscript are included in the `data/` directory.

---

## Code availability

All analyses were performed in R and are fully reproducible using the scripts provided in this repository.

---

## Citation

If you use this repository, please cite the associated manuscript once published.