# Plotting results -------------------------------------------------------------
rm(list = ls())
source("./R/header.R")
source("./R/02_par_process_results.R")
load("./R/data_and_results/par_eigenvalues.RData")
load("./R/data_and_results/par_mahala_nobis.RData")
load("./R/data_and_results/par_frob_distances.RData")
load("./R/data_and_results/par_coverage.RData")
load("./R/data_and_results/par_averaged_ci.RData")

font_import()
loadfonts()

# 1. table of proportion of times each method leads to pd covariance -----------
pd_proportion_df <- 1 - data.frame(
  eigenvalues_all$mcar_small$proportion_pd_matrices,
  eigenvalues_all$mar_small$proportion_pd_matrices,
  eigenvalues_all$mcar_large$proportion_pd_matrices,
  eigenvalues_all$mar_large$proportion_pd_matrices) %>%
  slice(-1)

names(pd_proportion_df) <- c("mcar_small", "mar_small", 
                             "mcar_large", "mar_large")
idx <- pd_proportion_df %>% is.na()
pd_proportion_df[idx] <- 0

tab_prop_pd <- knitr::kable(
  pd_proportion_df %>% 
    mutate(across(everything(), ~ sprintf("%.3f", .x))), digits = 3  
  )

# 2. table summary -------------------------------------------------------------
data_by_setting <- list_transpose(data)
settings <- names(data_by_setting)


missing_df <- imap_dfr(data_by_setting, function(sim_list, setting_name) {
  imap_dfr(sim_list, function(mat, simrun_id) {
    tibble(setting = setting_name,
           simrun  = as.integer(simrun_id),
           n_missings = mean(is.na(mat)))
      })}) %>% 
  filter(setting != "X") %>% 
  mutate(setting = factor(setting,
                          levels = c("X_mcar_small", "X_mar_small", 
                                     "X_mcar_large", "X_mar_large"),
                          labels = c("MCAR small", "MAR small", 
                                     "MCAR large", "MAR large"))) %>% 
  group_by(setting) %>% 
  summarize(Mean = 1-mean(n_missings),
            Median = 1-median(n_missings),
            SD = sd(n_missings)) %>% 
  as.data.frame() %>% 
  rename(`Mechanism` = setting) %>% 
  mutate(Missing = c("Small", "Small", "Large", "Large"),
         Mechanism = c("MCAR", "MAR", "MCAR", "MAR"), .before = 1)

# proportion of observed
missing_df
tab_missing_summary <- knitr::kable(missing_df, digits = 3)


# 3a. mahala Nobis Distance Boxplot ---------------------------------------------
mahala_nobis_df <- mahala_nobis_distances_all %>% 
  mutate(method = factor(method,
                         levels = c("listwise", "pairwise", "ipca_reg", "mi"),
                         labels = c("Listwise", "Pairwise", "IPCA-reg", "MI")),
  setting = factor(setting, 
                   levels = c("mcar_small", "mar_small", 
                              "mcar_large", "mar_large"),
                   labels = c("MCAR small", "MAR small", 
                              "MCAR large", "MAR large"))) %>% 
  drop_na()


(mahala_nobis_distanz <- ggplot(mahala_nobis_df, 
                                aes(x = as.factor(method), y = value)) +
  geom_boxplot(show.legend = TRUE) +
  facet_wrap(~ setting, axis.labels = "all_x", scales = "free_y") +
    facetted_pos_scales(
      y = list(setting == "MCAR small" ~ scale_y_continuous(limits = c(0,0.5)),
               setting == "MCAR large" ~ scale_y_continuous(limits = c(0,0.3)),
        TRUE ~ scale_y_continuous())) +
  labs(x = NULL,
       y = expression(d(S[italic(comp.)], S[italic(miss.)]))) +
  #theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))
)

(mahala_nobis_distanz_rep <- ggplot(mahala_nobis_df %>% filter(setting %in% c("MCAR small", "MCAR large")), 
                                aes(x = as.factor(method), y = value)) +
    geom_boxplot(show.legend = TRUE) +
    facet_wrap(~ setting, axis.labels = "all_x", scales = "free_y") +
    facetted_pos_scales(
      y = list(setting == "MCAR small" ~ scale_y_continuous(limits = c(0,0.5)),
               setting == "MCAR large" ~ scale_y_continuous(limits = c(0,0.3)),
               TRUE ~ scale_y_continuous())) +
    labs(x = NULL,
         y = expression(d(S[italic(comp.)], S[italic(miss.)]))) +
    #theme_minimal() +
    theme(axis.text.x = element_text(angle = 90))
)

# 3b. frobenius norm ----------------------------------------------------------
# Matrix-Abweichung: Frobenius-Norm (Robust)
(frobenius_norm <- frob_distances_all %>%
  mutate(method = factor(method, 
                         levels = c("listwise", "pairwise", "ipca_reg", "mi"),
                         labels = c("Listwise", "Pairwise", "IPCA-reg", "MI")),
         setting = factor(setting, 
                          levels = c("mcar_small", "mar_small", 
                                     "mcar_large", "mar_large"),
                          labels = c("MCAR small", "MAR small", 
                                     "MCAR large", "MAR large"))) %>%
  #mutate(value = value/280) %>% 
  drop_na(value) %>%
  ggplot(aes(x = method, y = value)) +
  geom_boxplot() +
  facet_wrap(~ setting, scales = "free_y") +
  coord_cartesian(ylim = c(0, 100)) +
  labs(y = "Frobenius Norm",
       x = NULL) +
  theme(axis.text.x = element_text(angle = 90))
)

(frobenius_norm_rep <- frob_distances_all %>%
    mutate(method = factor(method, 
                           levels = c("listwise", "pairwise", "ipca_reg", "mi"),
                           labels = c("Listwise", "Pairwise", "IPCA-reg", "MI")),
           setting = factor(setting, 
                            levels = c("mcar_small", "mar_small", 
                                       "mcar_large", "mar_large"),
                            labels = c("MCAR small", "MAR small", 
                                       "MCAR large", "MAR large"))) %>%
    filter(setting %in% c("MCAR small", "MCAR large")) %>%
    #mutate(value = value/280) %>% 
    drop_na(value) %>%
    ggplot(aes(x = method, y = value)) +
    geom_boxplot() +
    facet_wrap(~ setting, scales = "free_y") +
    coord_cartesian(ylim = c(0, 100)) +
    labs(y = "Frobenius Norm",
         x = NULL) +
    theme(axis.text.x = element_text(angle = 90))
)

# 4. eigenvalue difference -----------------------------------------------------
plot_mcar_small <- ggplot(
  eigenvalues_all$mcar_small$eigenvalue_share_diff_long, 
  aes(x = as.factor(k), y = gamma_diff)) +
   geom_boxplot(show.legend = TRUE) +
   geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
   facet_wrap(~ method, axis.labels = "all_x", scales = "free_y") +
   facetted_pos_scales(
     y = list(method == "Listwise" ~ scale_y_continuous(limits = c(-0.06,0.15)),
              method == "Pairwise" ~ scale_y_continuous(limits = c(-0.03,0.03)),
              method == "IPCA-reg" ~ scale_y_continuous(limits = c(-0.03,0.04)),
              method == "MI" ~ scale_y_continuous(limits = c(-0.03,0.04)),
              TRUE ~ scale_y_continuous())) +
   labs(x = "Number of factors k",
        y = expression(hat(gamma)[italic(k)[italic(miss.)]] - 
                         hat(gamma)[italic(k)[italic(comp.)]]))


plot_mcar_large <- ggplot(eigenvalues_all$mcar_large$eigenvalue_share_diff_long %>% drop_na(gamma_diff), 
                           aes(x = as.factor(k), y = gamma_diff)) +
    geom_boxplot(show.legend = TRUE) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    facet_wrap(~ method, axis.labels = "all_x", scales = "free_y") +
    facetted_pos_scales(
      y = list(
        method == "Listwise" ~ scale_y_continuous(limits = c(-0.06,0.15)),
        method == "Pairwise" ~ scale_y_continuous(limits = c(-0.03,0.04)),
        method == "IPCA-reg" ~ scale_y_continuous(limits = c(-0.1,0.2)),
        method == "MI" ~ scale_y_continuous(limits = c(-0.05,0.1)),
        TRUE ~ scale_y_continuous())) +
    labs(x = "Number of factors k",
         y = expression(hat(gamma)[italic(k)[italic(miss.)]] - 
                          hat(gamma)[italic(k)[italic(comp.)]]))


plot_mar_small <- ggplot(eigenvalues_all$mar_small$eigenvalue_share_diff_long, 
                           aes(x = as.factor(k), y = gamma_diff)) +
    geom_boxplot(show.legend = TRUE) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    facet_wrap(~ method, axis.labels = "all_x", scales = "free_y") +
    labs(x = "Number of factors k",
         y = expression(hat(gamma)[italic(k)[italic(miss.)]] - 
                          hat(gamma)[italic(k)[italic(comp.)]]))


plot_mar_large <- ggplot(eigenvalues_all$mar_large$eigenvalue_share_diff_long, 
                          aes(x = as.factor(k), y = gamma_diff)) +
    geom_boxplot(show.legend = TRUE) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    facet_wrap(~ method, axis.labels = "all_x", scales = "free_y") +
    labs(x = "Number of factors k",
         y = expression(hat(gamma)[italic(k)[italic(miss.)]] - 
                          hat(gamma)[italic(k)[italic(comp.)]]))



# 5. coverage and ci -----------------------------------------------------------
df_ci_coverage <- conf_int_averaged %>% 
  left_join(coverage, by = c("condition", "n_pc"))

  
df_plot_tidy <- df_ci_coverage %>%
  pivot_longer(
    cols = c(coverage_boot, coverage_fieller, width_boot, width_fieller),
    names_to = c(".value", "method"),
    names_sep = "_"
  ) %>%
  mutate(method = recode(method, boot = "Bootstrap", fieller = "Fieller"))


# Coverage (x) vs. RRMSE (y)
plot_ci1 <- ggplot(df_plot_tidy, aes(x = coverage, y = rrmse, color = method)) +
  geom_vline(xintercept = 0.95, linetype = "dashed", 
             color = "red", linewidth = 0.8) +
  geom_point(aes(size = n_pc, alpha = n_pc)) +
  facet_wrap(~condition) +
  scale_x_continuous(labels = scales::percent_format(), limits = c(NA, 1)) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_size_continuous(range = c(2, 10), name = "PC Index (k)") +
  labs(
    title = "Analysis of Estimation Precision and Reliability",
    subtitle = "Relative RMSE vs. Confidence Interval Coverage",
    x = "Coverage Probability (Target = 95%)",
    y = "Relative RMSE (Prediction Error)",
    color = "CI Method"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

plot_ci2 <- ggplot(df_plot_tidy, aes(x = coverage, y = rrmse, color = method)) +
  geom_vline(xintercept = 0.95, linetype = "dashed", color = "red") +
  facet_wrap(~condition) +
  geom_point(aes(size = width, alpha = n_pc)) +
  scale_size_continuous(range = c(2, 10), name = "CI Width (Avg)") +
  scale_alpha_continuous(range = c(0.3, 1), name = "PC Index (k)") +
  scale_x_continuous(labels = scales::percent_format(), limits = c(NA, 1)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Accuracy vs. Coverage vs. Precision (Width)",
    subtitle = "Larger points indicate wider (less precise) intervals",
    x = "Coverage Probability",
    y = "Relative RMSE"
  ) +
  theme_minimal()



plot_ci3 <- ggplot(df_plot_tidy, aes(x = coverage, y = rrmse)) +
  annotate("rect", xmin = 0.90, xmax = 1.01, ymin = -Inf, ymax = Inf, 
           fill = "grey90", alpha = 0.5) +
  geom_vline(xintercept = 0.95, linetype = "dashed", color = "firebrick", size = 0.7) +
  geom_point(aes(fill = method, size = width, alpha = n_pc), 
             shape = 21, color = "white", stroke = 0.5) +
  facet_wrap(~condition, labeller = as_labeller(c(mar_small = "MAR: Small Missingness", 
                                                  mar_large = "MAR: Large Missingness"))) +
  scale_fill_viridis_d(option = "plasma", end = 0.8, name = "Method") +
  scale_size_continuous(range = c(2, 12), name = "CI Width") +
  scale_alpha_continuous(range = c(0.4, 1), name = "PC Index (k)") +
  scale_x_continuous(labels = scales::percent_format(), expand = c(0.02, 0)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Simulation Results: Principal Component Variance Shares",
    subtitle = "Trade-off between Estimation Error (RRMSE) and Interval Reliability (Coverage)",
    x = "Coverage Probability",
    y = "Relative RMSE",
    caption = "Note: Red line indicates 95% nominal coverage. Shaded area represents the 90-100% reliability zone."
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 5, alpha = 1)), # Macht Methoden-Punkte in Legende groß & sichtbar
    size = guide_legend(override.aes = list(shape = 21, fill = "grey50", color = "white")), # Zeigt Kreise für die Breite
    alpha = guide_legend(override.aes = list(shape = 21, fill = "grey50", size = 4)) # Zeigt Kreise für den PC-Index
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey30")
  )

# save plots -------------------------------------------------------------------
# save(tab_missing_summary,
#      tab_prop_pd,
#      mahala_nobis_distanz, mahala_nobis_distanz_rep,
#      frobenius_norm, frobenius_norm_rep,
#      plot_mcar_small, plot_mcar_large,
#      plot_mar_small, plot_mar_large,
#      plot_ci1, plot_ci2, plot_ci3,
#      df_ci_coverage,
#      file = "./R/plots/tables_and_plots.RData")
# 
# ggsave("./R/plots/mahala_nobis_distance_boxplot.pdf", mahala_nobis_distanz,
#        width = 10, height = 10)
# ggsave("./R/plots/frobenius_norm_boxplot.pdf", frobenius_norm,
#        width = 10, height = 10)
# ggsave("./R/plots/eigenvalue_diff_mcar_small.pdf", plot_mcar_small,
#        width = 10, height = 6)
# ggsave("./R/plots/eigenvalue_diff_mcar_large.pdf", plot_mcar_large,
#        width = 10, height = 3)
# ggsave("./R/plots/eigenvalue_diff_mar_small.pdf", plot_mar_small,
#        width = 10, height = 6)
# ggsave("./R/plots/eigenvalue_diff_mar_large.pdf", plot_mar_large,
#        width = 10, height = 6)
# ggsave("./R/plots/ci_coverage_rrmse.pdf", plot_ci1,
#        width = 10, height = 6)
# ggsave("./R/plots/ci_coverage_rrmse_width.pdf", plot_ci2,
#        width = 10, height = 6)
# ggsave("./R/plots/ci_coverage_rrmse_width_nice.pdf", plot_ci3,
#        width = 10, height = 6)
