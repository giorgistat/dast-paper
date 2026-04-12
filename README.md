# 🧮 Decay-Adjusted Spatio-Temporal (DAST) Model
### Repository for the paper
**“A decay-adjusted spatio-temporal model to account for the impact of mass drug administration on neglected tropical disease prevalence”** by Emanuele Giorgi, Claudio Fronterre and Peter Diggle. A preprint of the paper is available on arXiv: https://arxiv.org/abs/2512.03760.

---

## 📂 Contents

This repository contains fully reproducible R scripts corresponding to the analyses and figures presented in the paper:

| Script | Description |
|---------|--------------|
| `01_scenario1_simulation.R` | Simulation to evaluate the predictive performance of the DAST model (Scenario 1). |
| `01_scenario2_simulation.R` | Simulation to evaluate the predictive performance of the DAST model (Scenario 2). |
| `simulation_functions.R` | Auxiliary functions used in the simulation study. |
| `02_sth_kenya_application.R` | Application of the DAST model to **soil-transmitted helminth (STH)** prevalence data from Kenya. |
| `03_lf_madagascar_application.R` | Application of the DAST model to **lymphatic filariasis (LF)** data from Madagascar, including comparison with GLM and GLMM models and projection of MDA rounds needed for elimination. |
| `04_figures.R` | Generates the figures reported in the paper using the analysis outputs and saves them to `figures/`. |

Each script can be executed independently and reproduces the corresponding figures, tables, and results reported in the paper.

---

## ▶️ Usage

All scripts are written for R (≥ 4.2). You can run each one directly from the command line or within R / RStudio.

### Command line
```bash
Rscript 01_scenario1_simulation.R
Rscript 01_scenario2_simulation.R
Rscript 02_sth_kenya_application.R
Rscript 03_lf_madagascar_application.R
Rscript 04_figures.R
```

Scripts assume the required input data files are available in the `data/` directory (not included here because of data-sharing agreements; please contact the authors for access). Output figures are saved to the `figures/` directory.
