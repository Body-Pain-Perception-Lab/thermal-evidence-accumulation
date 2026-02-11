# Evidence accumulation in temperature perception: Distinct effects of spatial summation and lateral inhibition

## Overview

This repository contains the analysis code and data for the paper "Evidence accumulation in temperature perception: Distinct effects of spatial summation and lateral inhibition". Which investigates spatial summation and lateral inhibition effects under different temporal conditions (dynamic thermal changes vs. fixed target temperatures).

Larger model files are available via the corresponding OSF repository.

Publication: Link forthcoming.

[OSF repository](https://osf.io/ev2mk/overview)


## Repository Structure

```
thermal-evidence-accumulation/
├── CITATION.cff               # Citation metadata for the repository
├── LICENSE                    # License information (MIT for code, CC BY 4.0 for text/figures)
├── Analysis/                   # Main frequentist analysis scripts and functions
│   ├── Analysis.Rmd            # Primary statistical analysis
│   ├── Plots.Rmd               # Data visualization
│   ├── functions/              # Custom R functions
│   │   ├── clean_data.R        # Data cleaning utilities
│   │   ├── load_data.R         # data loading
│   │   └── utils.R             # Analysis helper functions
│   ├── tables/                 # Table generation scripts
│   └── Workspaces/             # Saved R workspace files
├── Data/                       # Experimental data files
├── STAN/                      # Bayesian analysis with Stan
│   ├── scripts/               # Stan analysis scripts
│   ├── stanmodels/            # Compiled Stan models
│   └── taskmodels/            # Task-specific model files (on OSF)
├── Figures/                   # Generated plots and visualizations
│   ├── png/                   # PNG format figures
│   ├── tiff/                  # High-resolution TIFF figures
│   ├── parameter recovery/    # Model validation plots
│   └── Predictive_checks/     # Bayesian model diagnostics
```
