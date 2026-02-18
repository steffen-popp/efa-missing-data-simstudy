if (!requireNamespace("pacman")) {install.packages("pacman")}
if (!requireNamespace("eigeninv")) {
  install.packages(
    "https://cran.r-project.org/src/contrib/Archive/eigeninv/eigeninv_2011.8-1.tar.gz",
    repos = NULL,
    type = "source")
  }

pacman::p_load(
  tidyverse,       # data manipulation and visualization
  patchwork,       # combining ggplots
  eigeninv,        # create data from eigenvalues
  mice,            # imputation
  missMDA,         # Handling missing values multivariate data analysis (pca)
  mifa,            # Get covariance matrix of incomplete data using MI
  VIM,             # Visualization and Imputation of Missing Values
  furrr,           # parallel processing with purrr syntax
  matrixcalc,      # matrix calculation
  tictoc,          # time measurement
  doParallel,      # parallel calculation
  foreach,         # parallel calculation
  future,          # parallel calculation
  progressr,       # progress bars for parallel processing
  beepr,           # sound when calculation is done
  ggh4x,           # for individual y scales in facet wrap plots
  viridis,         # plot color palette
  gt,              # for nice tables
  gridExtra,       # for arranging gt tables
  extrafont,       # for nice fonts in plots
  scales           # for percents at axis labels
  )

set.seed(12965)
