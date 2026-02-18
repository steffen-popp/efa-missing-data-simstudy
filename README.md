# Exploratory factor analysis with multiply imputed data under Missing at Random
## Steffen Popp, Clarissa Spiegler

Contents:
Markdown
# Online Material: Missing Data Simulation & Reproduction

This repository contains the complete code and documentation for the simulation, reproduction, and analysis of missing data handling methods, with a particular focus on the performance of confidence intervals.

## 📁 Project Structure

### 📄 Documentation (HTML Reports)
The rendered results of the analysis:
* **`doc_reproduction.html`**: A detailed report covering the reproduction of the findings from the original study.
* **`doc_confidence_intervals.html`**: Analysis and comparison of different confidence interval methods (Bootstrap vs. Fieller).

### 💻 R Scripts (`/R`)
The source code is organized into sequential phases:
* **`01_data_generation.R`**: Scripts for initializing and generating the simulation data.
* **`02_par_calc_sim_.R`**: Parallelized computation of the main simulation.
* **`02_par_process_results.R`**: Data cleaning and merging of the parallelized output.
* **`03_analysis_results.R`**: Statistical evaluation of the simulation outcomes.
* **`04_plots.R`** & **`04_plots_poster.R`**: Visualization scripts for the final reports and poster presentations.

**Core Functions & Configuration:**
* **`header.R`**: Central setup file containing required packages, global paths, and variables.
* **`function_covariances.R`**: Custom functions for estimating various covariance matrices.
* **`function_mahalanobis_measures.R`**: Implementation of distance metrics (including Mahalanobis distance and Frobenius norm).

### 📊 Data & Results (`/R/data_and_results`)
* **`1000_sim_data.RData`**: Pre-computed simulation results based on 1000 iterations for immediate analysis.

---

## 🚀 Getting Started

1.  **Execution Order:** To replicate the study from scratch, run the scripts in chronological order (`01` to `04`).
2.  **Dependencies:** Ensure that all R packages specified in `header.R` are installed on your system.
3.  **Working Directory:** All scripts assume the project root (`online_material_missingdata/`) as the working directory.

## 📝 Methodological Notes

* **Reproduction of Distances:** During the replication of Mahalanobis-based distances, it was found that covariance matrices for Listwise and Pairwise deletion methods were frequently not **positive definite**. This is consistent with the original paper's findings (reporting 0% and 0.2% positive definiteness for these methods, respectively).
* **Alternative Metric:** Due to these constraints, the **Frobenius norm** was implemented as an alternative distance measure. While the numerical scale differs from the original paper, the resulting patterns and structural insights remain highly consistent.

---
*Last updated: February 2026*
