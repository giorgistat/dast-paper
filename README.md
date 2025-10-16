# 🧮 Decay-Adjusted Spatio-Temporal (DAST) Model
### Repository for the paper  
**“A decay-adjusted spatio-temporal model to account for the impact of mass drug administration on neglected tropical disease prevalence”** by Emanuele Giorgi, Claudio Fronterre and Peter Diggle. 

---

## 📂 Contents

This repository contains three fully reproducible R scripts corresponding to the analyses presented in the paper:

| Script | Description |
|---------|--------------|
| `01_simulation_study.R` | Simulation experiments used to evaluate the statistical performance and robustness of the DAST model under controlled conditions. |
| `02_sth_kenya_application.R` | Application of the DAST model to **soil-transmitted helminth (STH)** prevalence data from Kenya. |
| `03_lf_madagascar_application.R` | Application of the DAST model to **lymphatic filariasis (LF)** data from Madagascar, including comparison with GLM and GLMM models and projection of MDA rounds needed for elimination. |

Each script can be executed independently and reproduces the corresponding figures, tables, and results reported in the paper.

---

## ▶️ Usage

All scripts are written for R (≥ 4.2).  
You can run each one directly from the command line or within R / RStudio.

### Command line
```bash
Rscript 01_simulation_study.R
Rscript 02_sth_kenya_application.R
Rscript 03_lf_madagascar_application.R
