# Drosophila-climbing-assay



Reproducible R code and data accompanying the manuscript:



> \\\\\\\*\\\\\\\*An individual-tube climbing assay with transparent analysis and quality control for Drosophila locomotor phenotyping\\\\\\\*\\\\\\\*  

> (Submitted to \\\\\\\*Animal Behaviour\\\\\\\*)



\## Overview



This repository contains the data and R scripts used to analyse an individual-tube negative geotaxis assay in \*Drosophila melanogaster\*.



The repository reproduces the principal analyses presented in the manuscript, including:



\- Longitudinal age-profile analysis and Day 21 readout-window selection.

\- Day 21 vehicle-control reproducibility analysis.

\- Salbutamol worked pharmacological application.

\- Supplementary sensitivity analyses based on round-level weighted AUC.



The scripts generate the statistical models, estimated marginal means, contrasts, diagnostic summaries, supplementary tables and `ggplot2` objects underlying the manuscript figures. Final figure layout and SVG export were performed manually during manuscript preparation.



\---



\## Repository structure



```

Drosophila-climbing-assay/

│

├── data/

│   ├── longitudinal.xlsx

│   ├── day21\\\\\\\_vehicle\\\\\\\_only.xlsx

│   └── salbutamol\\\\\\\_application.xlsx

│

├── scripts/

│   ├── 01\\\\\\\_longitudinal\\\\\\\_analysis.R

│   ├── 02\\\\\\\_vehicle\\\\\\\_control\\\\\\\_analysis.R

│   └── 03\\\\\\\_salbutamol\\\\\\\_application.R

│

├── output/

│   ├── figures/

│   └── tables/

│

├── README.md

├── LICENSE

└── .gitignore

```



\---



\## Requirements



The analyses were performed in \*\*R (version 4.x)\*\* using the following packages:



\- tidyverse

\- readxl

\- lme4

\- emmeans

\- DHARMa

\- performance

\- patchwork

\- scales

\- writexl



\---



\## Running the analyses



Each script is independent and can be run separately.



1\. Run `01\\\\\\\_longitudinal\\\\\\\_analysis.R`

2\. Run `02\\\\\\\_vehicle\\\\\\\_control\\\\\\\_analysis.R`

3\. Run `03\\\\\\\_salbutamol\\\\\\\_application.R`



The scripts assume that the working directory is the repository root.



\---



\## Outputs



The scripts generate:



\- fitted GLMM objects

\- estimated marginal means

\- treatment and genotype contrasts

\- model diagnostics

\- supplementary tables

\- `ggplot2` objects corresponding to the manuscript figures



Supplementary tables are written to the `output/tables/` directory.



\---



\## Citation



If you use this repository, please cite the associated manuscript once published.

