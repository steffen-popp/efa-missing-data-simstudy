# Exploratory factor analysis with multiply imputed data under Missing at Random
**Authors:** Steffen Popp, Clarissa Spiegler


# Online Material

## 📁 Project Structure

### 📄 References and Documentation (HTML Reports)
* **`doc_reproduction.html`**: A report with the reproduction results from the original study (Reproduced settings: MCAR small and MCAR large).
* **`doc_confidence_intervals.html`**: Analysis and comparison of different confidence interval methods (Bootstrap vs. Fieller).
* **`references.pdf`**: Contains the References.

### 💻 R Scripts (`/R`)
The source code is organized into sequential phases:
* **`01_data_generation.R`**: Scripts for generating the simulation data.
* **`02_par_calc_sim_.R`**: Parallelized computation of the different simulation settings (reproduction and MAR setting).
* **`02_par_process_results.R`**: Data preparation for further analysis.
* **`03_analysis_results.R`**: Analysis of simulation results.
* **`04_plots.R`** & **`04_plots_poster.R`**: Scripts for visualization of simulation results.

**Core Functions & Configuration:**
* **`header.R`**: Script containing required packages.
* **`function_covariances.R`**: Custom functions for estimating the covariance matrices with different methods.
* **`function_mahalanobis_measures.R`**: Custom functions for calculation of mahalanobis-based distances.

### 📊 Data & Results (`/R/data_and_results`)
* **`1000_sim_data.RData`**: Generated simulation data for 1000 iterations (generated with **`01_data_generation.R`**).

## AI Disclosure

**Note on the use of Artificial Intelligence (AI):**
AI tools were employed for technical assistance in coding (R), particularly in the creation of data visualizations (plots) and optimization of the code. This included enhancing computational efficiency through the parallelization of simulation runs.
---

