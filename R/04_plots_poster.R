# load data and results -----------------------------------------------------
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


# color palette for method comparison
pal=c(#"#003f5c",
      "#00457D", # uni blau
      #"#2f4b7c",
      "#665191",
      #"#a05195",
      "#d45087",
      #"#f95d6a",
      "#ff7c43"
      #"#ffa600"
      )

# color palettes for confidence interval plot
pal2 = c("#004e64", "#9a031e")

# Textsizes
ax_txt = 15


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

plot_data <- pd_proportion_df %>%
  mutate(Method = rownames(.)) %>%
  pivot_longer(cols = -Method, 
               names_to = "Scenario", 
               values_to = "Proportion") %>%
  mutate(Scenario = factor(Scenario, 
                           levels = c("mcar_small", "mar_small", 
                                      "mcar_large", "mar_large"),
                           labels = c("MCAR small", "MAR small", 
                                      "MCAR large", "MAR large"))) %>% 
  mutate(Method = factor(Method,
                          levels = c("listwise", "pairwise", 
                                     "ipca_reg", "mi"),
                          labels = c("Listwise", "Pairwise", 
                                     "IPCA-reg", "MI")))

# barplot of proportion of pd matrices for MAR large
plot_bar_pd <- ggplot(plot_data %>% filter(Scenario == "MAR large"), 
                      aes(x = Method, y = Proportion, fill = Method)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values=pal) +
  scale_color_manual(values=pal) +
  labs(x = "Method",
       y = "Proportion",
       fill = "Method"
  ) +
  theme(legend.position = "none",
      text = element_text(family = "Times New Roman"),
      title = element_text(size = 18),
      axis.ticks = element_blank(),
      axis.text = element_text(size = 20),
      axis.title.x = element_text(margin = margin(t = 10), size = 25),
      axis.title.y = element_text(margin = margin(r = 10), size = 25),
      panel.grid.major.x = element_blank(),
      plot.background = element_blank(),
      panel.background = element_blank())


# 2. frobenius norm -----------------------------------------------------
frobenius_plot_df <- frob_distances_all %>%
  mutate(method = factor(method, 
                         levels = c("listwise", "pairwise", "ipca_reg", "mi"),
                         labels = c("Listwise", "Pairwise", "IPCA-reg", "MI")),
         setting = factor(setting, 
                          levels = c("mcar_small", "mar_small", 
                                     "mcar_large", "mar_large"))) %>%
  filter(setting == "mar_large") %>% 
  drop_na(value)

(frobenius_norm <- ggplot(frobenius_plot_df, 
                          aes(x = method, y = value, fill = method)) +
  geom_boxplot(alpha = 0.9) +
  scale_fill_manual(values=pal) +
  scale_color_manual(values=pal) +
  labs(y = "Frobenius Norm",
       x = "Method", fill = "Method") +
  theme(legend.position = "none",
         text = element_text(family = "Times New Roman"),
         title = element_text(size = 18),
         axis.ticks = element_blank(),
         axis.text = element_text(size = 20),
         axis.title.x = element_text(margin = margin(t = 10), size = 25),
         axis.title.y = element_text(margin = margin(r = 10), size = 25),
         plot.background = element_blank(),
         panel.background = element_blank()))




# 4. eigenvalue difference -----------------------------------------------------
plot_mar_large <- ggplot(eigenvalues_all$mar_large$eigenvalue_share_diff_long, 
                          aes(x = as.factor(k), y = gamma_diff, fill = method))+
    geom_boxplot(show.legend = TRUE, alpha = 0.8) +
    scale_fill_manual(values=pal) +
    scale_color_manual(values=pal) +
    geom_hline(yintercept = 0, color = "red", 
               linetype = "dashed", linewidth = 0.8) +
    facet_wrap(~ method, axis.labels = "all_x", 
               #scales = "free_y", 
               ncol = 2) +
    labs(x = "Number of Factors k",
         y = expression(Delta[gamma * "," * italic(k)])) +
    theme(legend.position = "none",
          text = element_text(family = "Times New Roman"),
          title = element_text(size = 18),
          axis.ticks = element_blank(),
          axis.text = element_text(size = ax_txt),
          axis.title.x = element_text(margin = margin(t = 10), size = 18),
          axis.title.y = element_text(margin = margin(r = 10), size = 18),
          strip.text = element_text(size = 18, face = "bold"),
          plot.background = element_blank(),
          panel.background = element_blank())


# plot boxplots and only one row
eig_mar_large_croped <- eigenvalues_all$mar_large$eigenvalue_share_diff_long %>% 
  filter(k <=10)
  
plot_mar_large_croped <- ggplot(eig_mar_large_croped, 
                          aes(x = as.factor(k), y = gamma_diff, fill = method))+
    geom_boxplot(show.legend = TRUE, alpha = 0.8) +
    scale_fill_manual(values=pal) +
    scale_color_manual(values=pal) +
    geom_hline(yintercept = 0, color = "red", 
               linetype = "dashed", linewidth = 0.8) +
    facet_wrap(~ method, axis.labels = "all_x", 
               #scales = "free_y", 
               ncol = 4) +
    coord_cartesian(ylim = c(-0.05, 0.15)) +
    labs(x = "Number of Factors k",
         y = expression(Delta[gamma * "," * italic(k)])) +
    theme(legend.position = "none",
          text = element_text(family = "Times New Roman"),
          title = element_text(size = 18),
          axis.ticks = element_blank(),
          axis.text = element_text(size = ax_txt),
          axis.title.x = element_text(margin = margin(t = 10), size = 18),
          axis.title.y = element_text(margin = margin(r = 5), size = 18),
          strip.text = element_text(size = 18, face = "bold"),
          plot.background = element_blank(),
          panel.background = element_blank())


# 5. coverage and ci ---------------------------------------------------
df_ci_coverage <- conf_int_averaged %>% 
  left_join(coverage, by = c("condition", "n_pc"))

df_plot_tidy <- df_ci_coverage %>%
  pivot_longer(
    cols = c(coverage_boot, coverage_fieller, width_boot, width_fieller),
    names_to = c(".value", "method"),
    names_sep = "_"
  ) %>%
  mutate(method = recode(method, boot = "Bootstrap", fieller = "Fieller")) %>% 
  filter(condition == "mar_large")


plot_ci <- ggplot(df_plot_tidy, aes(x = coverage, y = rrmse, fill = method)) +
  geom_vline(xintercept = 0.95, linetype = "dashed", 
             color = "firebrick", linewidth = 0.7) +
  geom_point(aes(fill = method, size = width, alpha = n_pc), 
             shape = 21, color = "white", stroke = 0.5) +
  geom_point(aes(fill = method), shape = 16, size = 0.8, color = "black",
             show.legend = FALSE) +
  scale_fill_manual(values=pal2) +
  scale_color_manual(values=pal2) +
  scale_size_continuous(range = c(2, 12), name = "CI Width") +
  scale_alpha_continuous(range = c(0.4, 1), name = "Nr. of Factors k") +
  scale_x_continuous(labels = scales::percent_format(), expand = c(0.02, 0)) +
  scale_y_continuous(labels = scales::percent_format()) +
  coord_cartesian(xlim = c(0, 1.05)) +
  labs(x = "Coverage", y = "RRMSE", fill = "Method") +
  guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1)), 
         size = guide_legend(override.aes = list(shape = 21, fill = "grey50", 
                                                 color = "white")),
         alpha = guide_legend(override.aes = list(shape = 21, fill = "grey30", 
                                                  size = 4))) +
  theme(legend.position = "right",
        text = element_text(family = "Times New Roman"),
        title = element_text(size = 18),
        axis.ticks = element_blank(),
        axis.text = element_text(size = ax_txt),
        axis.title.x = element_text(margin = margin(t = 10), size = 18),
        axis.title.y = element_text(margin = margin(r = 10), size = 18),
        strip.background = element_rect(fill = "grey95", color = NA),
        strip.text = element_text(size = 18, face = "bold"),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(), 
        legend.background = element_blank())


# Save plots -------------------------------------------------------------------
# ggsave("./Poster/images/proportion_pd_barplot.pdf", plot_bar_pd,
#        width = 6.92, height = 4.51, device = cairo_pdf)
# ggsave("./Poster/images/frobenius_norm_boxplot.pdf", frobenius_norm,
#        width = 6.92, height = 4.51, device = cairo_pdf)
# ggsave("./Poster/images/eigenvalue_diff_mar_large.pdf", plot_mar_large,
#        width = 12, height = 8, device = cairo_pdf)
# ggsave("./Poster/images/eigenvalue_diff_mar_large_croped.pdf",
#        plot_mar_large_croped, width = 12, height = 4, device = cairo_pdf)
# ggsave("./Poster/images/ci_coverage_rrmse_width.pdf", plot_ci,
#        width = 10, height = 4.5, device = cairo_pdf)