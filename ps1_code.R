library(haven)
library(tidyverse)
library(here)
library(rlang)
library(kableExtra)

# Main analysis ----------------------------------------------------------------
# Load helper scripts ---- 
source(here("ps1_functions.R"))

# Declare useful paths ----
output <- here("output")
plots <- here(output, "plots")

# Load the VV data ---- 
df_vv <- read_dta(here("voigtlaender-voth", "data_vv_clean.dta")) %>% 
  as.data.frame()
df_vv <- df_vv %>% filter(exist1349 == 1)

# 2b: Panel A: OLS ----
allvars_vv <- c("pog1349", "logpop25c", "perc_JEW25", "perc_PROT25")
panelA_vv <- ols_run(df_vv, y = "pog20s", covars = allvars_vv)
se_panelA_vv <- se_clust(df_vv, panelA_vv, clust_var = "kreis_nr")

# 2b: Panel B: Matching on Covariates ---- 
covars_vv <- c("logpop25c", "perc_JEW25", "perc_PROT25")
panelB_vv <- find_treated_nn(df_vv, y = "pog20s", covars = covars_vv,
                             treat_var = "pog1349", k = 4, estimand = "att")

se_panelB_vv <- se_boot(panelB_vv, B = 100, var_by = "pog1349")

# 2b: Panel C: Matching using Geography ----
covars_geo_vv <- c("Longitude", "Latitude")
panelC_vv <- find_treated_nn(df_vv, y = "pog20s", covars = covars_geo_vv,
                             treat_var = "pog1349", k = 2)
se_panelC_vv <- se_boot(panelC_vv, B = 100, var_by = "pog1349")

# 2e: Pscore Matching ----
nbins <- 6
pscore_mod_vv <- logit_run(df_vv, y = "pog1349", covars = covars_vv)

df_vv_block <- get_df_blocked(df_vv, y = "pog20s", treat_var = "pog1349", 
                              pscore_mod_vv, nbins = nbins)

# Inspect balance on covariates within blocks ----
# The covariates don't look very balanced, so I choose to use nearest neighbor
# matching based on the propensity score.
# Bin size is 6 because anything else leads to the covariates having missing 
# values.
df_vv_counts <- df_vv_block %>%
  group_by(block, pog1349) %>%
  summarise(n = n())

covars_vv %>% 
  walk(function(x) {
    plot_df_blocked(df_vv_block, y = x, treat_var = "pog1349")
    ggplot2::ggsave(here(plots, paste0(x, ".png")))
  })

# NN matching based on the propensity score ---- 
ps_att_vv <- find_treated_nn(
  df_vv_block, 
  y = "pog20s", 
  covars = "px",
  treat_var = "pog1349", 
  k = 1,
  estimand = "att"
)

se_ps_att_vv <- se_boot(ps_att_vv, B = 100, var_by = "pog1349")


ps_atu_vv <- find_treated_nn(
  df_vv_block, 
  y = "pog20s", 
  covars = "px",
  treat_var = "pog1349", 
  k = 1,
  estimand = "atu"
)

se_ps_atu_vv <- se_boot(ps_atu_vv, B = 100, var_by = "pog1349")

ps_ate_vv <- find_treated_nn(
  df_vv_block, 
  y = "pog20s", 
  covars = "px",
  treat_var = "pog1349", 
  k = 1,
  estimand = "ate"
)

se_ps_ate_vv <- se_boot(ps_ate_vv, B = 100, var_by = "pog1349")

# 7: Monte Carlo ---------------------------------------------------------------
# Helper functions ----
# Function to calculate the weights for a given cohort in a given year
weight_calc <- function(df_resid, ref_cohort) {
  props <- df_resid %>% 
    group_by(cohort) %>% 
    summarise(
      prop = n()/nrow(df_resid),
      .groups = "keep"
    )
  
  df_sum <- df_resid %>% 
    group_by(cohort, year) %>% 
    summarise(
      mean_resid = mean(resid),
      .groups = "keep"
    ) 
  
  df_out <- df_sum %>% 
    left_join(props, by = "cohort") %>% 
    ungroup() %>% 
    mutate(wt = mean_resid * prop)
  
  df_out <- df_out %>% 
    filter(cohort != ref_cohort & year != 1)
}

# Relative weights for cohort e at relative time g, looping over 
# different relative times r 
get_wts_e_r <- function(df_did, e, r, g, rel_times_exclude) {
  df_reg <- df_did %>% 
    select(cohort, starts_with("rel_time"), starts_with("unit_"),
           starts_with("year_")) %>% 
    mutate(tag_obs = if_else(cohort == e & rel_time == r, 1, 0)) %>% 
    mutate(tag_obs = replace_na(tag_obs, 0))
  
  # Exclude unit, year, and relative time dummies
  unit_vars <- grep("unit_", names(df_reg), value = TRUE)[-1]
  year_vars <- grep("year_", names(df_reg), value = TRUE)[-1]
  rel_times <- grep("rel_time_", names(df_reg), value = TRUE)
  
  rel_times_include <- rel_times[!(rel_times %in% paste0("rel_time_", 
                                                         rel_times_exclude))]

  
  # TWFE to get the weights: indicator for cohort e and relative time r 
  # on unit and time FEs + relative time indicators
  wts_reg <- ols_run(df_reg, y = "tag_obs", covars = c(unit_vars,
                                                       year_vars,
                                                       rel_times_include))
  
  # Grab the correct coefficients 
  # +1 to account for the constant
  coefs_tag <- which(wts_reg[["var_names"]] == paste0("rel_time_", g)) + 1
  coefs_e_r <- wts_reg[["estimates"]][coefs_tag][1]
}

# Function that grabs all the weights used in estimating the dynamic 
# TWFE coefficient on a single relative time indicator
# Returns a list of the weights on each cohort e at relative time r
# g = the coefficient to estimate 
get_wts_dynamic <- function(df_did, g, rel_times_exclude) {
  cohorts_times <- df_did %>% 
    distinct(cohort, rel_time) %>% 
    arrange(cohort, rel_time)
  
  wts <- map2(
    cohorts_times$cohort,
    cohorts_times$rel_time,
    ~ as.data.frame(get_wts_e_r(df_did, e = .x, r = .y, g = g,
                                rel_times_exclude = rel_times_exclude))
  ) %>% 
    list_rbind() %>% 
    cbind(cohorts_times$cohort) %>% 
    cbind(cohorts_times$rel_time)
  
  names(wts) <- c("wt", "cohort", "rel_time")
  
  return(wts)
}

# Calculates the implied coefficient on a certain relative time in the TWFE
# dynamic specification
coef_calculate <- function(df_wts) {
  df_wts <- df_wts %>% 
    mutate(att_contrib = att * wt) %>% 
    mutate(att_homog_contrib = att_homog * wt) %>% 
    summarise(b_att = sum(att_contrib, na.rm = TRUE),
              b_att_homog = sum(att_homog_contrib, na.rm = TRUE),
              .groups = "keep")
  
}

plot_wt_rel_time <- function(df_wts, r) {
  ggplot(df_wts) + 
    geom_point(aes(x = rel_time, y = wt, 
                   color = as.factor(cohort))) + 
    geom_line(aes(x = rel_time, y = wt, group = cohort,
                  color = as.factor(cohort))) + 
    labs(x = "Relative Time", y = "Weight", color = "Cohort") + 
    ggtitle(paste0("Weights on Two-Way Fixed Effects Coefficient"), 
            subtitle = paste0("Dynamic Specification, Relative Time ", r)) + 
    theme_minimal()
}

coefs_plot <- function(df_coefs, g, M) {
  df_plot <- df_coefs %>% 
    mutate(var = if_else(var == "b_att_mean", "Heterogeneous ATT's", 
                          "Homogeneous ATT's"))
  ggplot(df_plot, aes(x = var, y = estim)) + 
    geom_point(color = "dodgerblue") + 
    geom_errorbar(aes(ymin = ci_95lo, ymax = ci_95up), width = 0.1,
                  color = "dodgerblue") + 
    labs(x = "", y = "Estimate") + 
    ggtitle(paste0("Coefficient Plot, Relative Time ", g),
            subtitle = paste0("Average Coefficient from ", M, " simulations")) + 
    theme_minimal() + 
    theme(legend.position = "none")
}

dynamic_wts_df_process <- function(df, trt_df, trt_df_homog) {
  df_out <- df %>% 
    group_by(cohort, rel_time) %>% 
    summarise(wt = mean(wt), .groups = "keep") 
  
  df_out <- df_out %>% 
    left_join(trt_df, by = c("cohort", "rel_time")) %>% 
    select(-year) %>% 
    left_join(trt_df_homog, by = c("cohort", "rel_time")) %>% 
    select(-year)
}

dynamic_coef_df_process <- function(df, M) {
  df_out <- df %>% 
    summarise(b_att_mean = mean(b_att),
              b_att_sd = sd(b_att),
              b_att_homog_mean = mean(b_att_homog),
              b_att_homog_sd = sd(b_att_homog),
              .groups = "keep") %>% 
    pivot_longer(cols = contains("_mean"), names_to = "var",
                 values_to = "estim") %>% 
    mutate(sd = if_else(var == "b_att_mean", b_att_sd, b_att_homog_sd)) %>% 
    select(-contains("_sd"))
  
  # Confidence intervals
  df_out <- df_out %>% 
    mutate(ci_95up = estim + qt(.975, M-1)*sd,
           ci_95lo = estim - qt(.975, M-1)*sd)
}

# Cohort-weighted ATT 
weighted_att_calc <- function(df, r) {
  df_att <- df %>% 
    filter(rel_time == r) %>% 
    group_by(cohort) %>% 
    summarise(size = n(),
              att = att,
              .groups = "keep") %>% 
    slice(1) %>% 
    ungroup() %>% 
    mutate(pr = size/sum(size))
  
  df_att <- df_att %>% 
    mutate(att = replace_na(att, 0)) %>% 
    mutate(wt_att_contrib = att * pr)
  
  att <- sum(df_att$wt_att_contrib)
}

# CS approach 
# Function that calculates ATT(g, t)
cs_gt_calc <- function(df, e, r) {
  # Setup for r >= 0 
  if (r >= 0) {
    df_chrt <- df %>% 
      filter(cohort == e) %>% 
      filter(year == e + r | year == e - 1) %>% 
      mutate(pre = if_else(year == e - 1, 1, 0))
    
    tag_not_yet <- df %>% 
      filter(year == e + r & d == 0) 
    
    chrts_not_yet <- unique(tag_not_yet$cohort)
    
    df_not_yet <- df %>% 
      filter(cohort %in% chrts_not_yet) %>% 
      filter(year == e + r | year == e - 1) %>% 
      mutate(pre = if_else(year == e - 1, 1, 0))
    
    if (nrow(df_not_yet) == 0) {
      return(NA_real_)
    } 
  } else {
    df_chrt <- df %>% 
      filter(cohort == e) %>% 
      filter(year == e + r | year == e + r - 1) %>% 
      mutate(pre = if_else(year == e +r - 1, 1, 0))
    
    tag_not_yet <- df %>% 
      filter(year == e + r & d == 0) 
    
    chrts_not_yet <- unique(tag_not_yet$cohort)
    
    df_not_yet <- df %>% 
      filter(cohort %in% chrts_not_yet) %>% 
      filter(year == e + r | year == e + r - 1) %>% 
      mutate(pre = if_else(year == e + r - 1, 1, 0))
    
    if (nrow(df_not_yet) == 0) {
      return(NA_real_)
    }
  }
  
  # Quantities 
  chrt_post <- mean(df_chrt[df_chrt[, "pre"] == 0, ]$y, na.rm = TRUE)
  chrt_pre <- mean(df_chrt[df_chrt[, "pre"] == 1, ]$y, na.rm = TRUE)
  not_yet_post <- mean(df_not_yet[df_not_yet[, "pre"] == 0, ]$y, na.rm = TRUE)
  not_yet_pre <- mean(df_not_yet[df_not_yet[, "pre"] == 1, ]$y, na.rm = TRUE)
  
  att_gt <- (chrt_post - chrt_pre) - (not_yet_post - not_yet_pre)
}

# Function that aggregates ATT(g, t) for a relative time 
cs_did <- function(df, r) {
  
  df_att <- df %>% 
    filter(rel_time == r) %>% 
    group_by(cohort) %>% 
    summarise(size = n(), .groups = "keep") 
  
  chts <- unique(df_att$cohort)
  
  # Run the main CS estimator
  att_gt <- chts %>% 
    map(~ cs_gt_calc(df, e = .x, r = r))
  
  df_att_gt <- tibble(
    cohort = chts,
    att_gt = as.numeric(att_gt)
  )
  
  # Left bind to the main df and calculate final estimate, dropping groups 
  # with no control comparisons
  df_att <- df_att %>% 
    left_join(df_att_gt, by = "cohort") %>% 
    ungroup() %>% 
    mutate(pr = size/sum(size)) %>% 
    mutate(att_contrib = pr * att_gt)
  
  att_e <- sum(df_att$att_contrib, na.rm = TRUE)
}

# Imputation approach 
impute_did <- function(df, r, rel_times_exclude = -1) {
  
  # Select untreated subsample and regress on cohort + time FEs
  df_d0 <- df %>% 
    filter(d == 0)
  
  cohorts_include <- paste0("cohort_", unique(df_d0$cohort)[-1])
  years_include <- paste0("year_", unique(df_d0$year)[-1])
  
  mod_d0 <- ols_run(df = df_d0, y = "y", 
                    covars = c(cohorts_include, years_include))
  
  estims_d0 <- mod_d0[["estimates"]]
  
  mat_predict <- matrixPrep(df, cols = c(cohorts_include, years_include),
                            constant = TRUE)
  
  # Get fitted values and subtract them from each actual outcome to construct
  # Y dot
  y0_fit <- mat_predict %*% estims_d0
  df <- df %>% 
    cbind(y0_fit = y0_fit) %>% 
    mutate(y_dot = y - y0_fit)
  
  # Regress Y dot on the relative time indicators.
  rel_times <- grep("rel_time_", names(df), value = TRUE)
  
  rel_times_include <- rel_times[!(rel_times %in% paste0("rel_time_", 
                                                         rel_times_exclude))]
  
  impute_mod <- ols_run(df, y = "y_dot", 
                        covars = rel_times_include)
  
  # Grab the right coefficient 
  coef_tag <- which(impute_mod[["var_names"]] == paste0("rel_time_", r)) + 1
  att_impute <- impute_mod[["estimates"]][coef_tag]
}

# Routine ----
# N units, M simulations 
did_mc <- function(N, M) {
  # Declare empty list for MC results 
  results_7a <- data.frame() 
  results_7a_log <- data.frame()
  results_7b <- data.frame()
  results_7c <- data.frame()
  results_7c_coefs <- data.frame()
  results_7d <- data.frame() 
  results_7d_coefs <- data.frame() 
  results_7e <- data.frame() 
  
  for (m in 1:M) {
    if (m == 1) {
      print("Iteration 1........................................")
    }
    if (m %% 10 == 0) {
      print(paste0("Iteration ", m, "........................................"))
    }
    # Setup ----
    # DID components 
    years <- 1:6
    cohorts <- c(3:5)
    probs <- rep(1/length(cohorts), times = length(cohorts))
    units <- 1:N
    
    unit_vec <- rep(units, times = length(years))
    years_vec <- rep(years, each = length(units))
    
    df_mc <- tibble(
      unit = unit_vec,
      year = years_vec
    )
    
    # Create the data ----
    # Assign units to cohorts 
    df_mc <- df_mc %>% 
      group_by(unit) %>% 
      mutate(cohort = sample(cohorts, size = 1, prob = probs)) %>% 
      ungroup() 
    
    # Cohort FEs
    # Gaussians are centered at positive numbers to play nice with log 
    df_mc <- df_mc %>% 
      mutate(alpha = if_else(cohort == 3, rnorm(1, mean = 10, sd = 3), NA_real_)) %>% 
      mutate(alpha = if_else(cohort == 4, rnorm(1, mean = 100, sd = 3), alpha)) %>% 
      mutate(alpha = if_else(cohort == 5, rnorm(1, mean = 20, sd = 3), alpha)) %>% 
      mutate(alpha = if_else(is.na(cohort), rnorm(1, mean = 5, sd = 3), alpha))
    
    # Time FEs: N(0, 2)
    df_mc <- df_mc %>% 
      group_by(year) %>% 
      mutate(lambda = rnorm(1, mean = 0, sd = 2)) %>% 
      ungroup() 
    
    # Unit shocks that vary by unit and time: N(0, 1)
    df_mc <- df_mc %>% 
      group_by(unit, year) %>% 
      mutate(eps = rnorm(1)) %>% 
      ungroup()
    
    # Treatment effects
    # For cohort 3: c(4, 8, 10, 10)
    chrt3_att <- tibble(
      year = c(3, 4, 5, 6),
      att = c(4, 8, 10, 10),
      cohort = 3
    )
    
    chrt3_att_homog <- tibble(
      year = c(3, 4, 5, 6),
      att_homog = c(4, 8, 10, 10),
      cohort = 3
    )
    
    # For cohort 4: c(NA_real_, 15, 20, 21)
    chrt4_att <- tibble(
      year = c(3, 4, 5, 6),
      att = c(NA_real_, 15, 20, 21),
      cohort = 4
    )
    
    chrt4_att_homog <- tibble(
      year = c(3, 4, 5, 6),
      att_homog = c(NA_real_, 4, 8, 10),
      cohort = 4
    )
    
    # For cohort 5
    chrt5_att <- tibble(
      year = c(3, 4, 5, 6),
      att = c(NA_real_, NA_real_, 0.5, 18),
      cohort = 5
    )
    
    chrt5_att_homog <- tibble(
      year = c(3, 4, 5, 6),
      att_homog = c(NA_real_, NA_real_, 4, 8),
      cohort = 5
    )
    
    trt_df <- rbind(chrt3_att, chrt4_att, chrt5_att) %>% 
      mutate(rel_time = year - cohort)
    trt_df_homog <- rbind(chrt3_att_homog, chrt4_att_homog, 
                          chrt5_att_homog) %>% 
      mutate(rel_time = year - cohort)
    
    # Bind to the main data 
    df_mc <- df_mc %>% 
      left_join(trt_df, by = c("cohort", "year")) %>% 
      select(-rel_time) %>% 
      left_join(trt_df_homog, by = c("cohort", "year")) %>% 
      select(-rel_time)
    
    # Observed Y_{it} conditional on cohort = cohort FE + time FE + treatment + 
    # idiosyncratic shock.
    # Also, generate dummies and relative time.
    df_mc <- df_mc %>% 
      mutate(d = if_else(year >= cohort, 1, 0)) %>% 
      mutate(d = if_else(is.na(cohort), 0, d)) %>% 
      mutate(y = if_else(!is.na(att), alpha + lambda + d * att + eps,
                         alpha + lambda + eps)) %>%
      # outcome with homogeneous relative time treatment effects 
      mutate(y_homog = if_else(!is.na(att), alpha + lambda + d * att_homog + eps,
                               alpha + lambda + eps)) %>% 
      mutate(logy = log(y)) %>% 
      mutate(rel_time = year - cohort) %>% 
      fastDummies::dummy_cols(select_columns = c("unit", "year", "cohort", 
                                                 "rel_time"))
    
    # Part 7a ---- 
    # Common trends: compare cohorts 3 and 4 from periods 1 to 2
    df_sum_mc <- df_mc %>% 
      filter(year == 1 | year == 2) %>% 
      group_by(cohort, year) %>% 
      summarise(y0 = mean(y), .groups = "keep")
    
    results_7a <- results_7a %>% 
      rbind(df_sum_mc)
    
    df_sum_log_mc <- df_mc %>% 
      filter(year == 1 | year == 2) %>% 
      group_by(cohort, year) %>% 
      summarise(logy0 = mean(logy, na.rm = TRUE), .groups = "keep")
    
    results_7a_log <- results_7a_log %>% 
      rbind(df_sum_log_mc) 
    
    # Part 7b ----
    unit_vars <- grep("unit_", names(df_mc), value = TRUE)[-1]
    year_vars <- grep("year_", names(df_mc), value = TRUE)[-1]
    d_resid_mod <- ols_run(df_mc, y = "d", covars = c(unit_vars, year_vars))
    d_resid <- as.vector(d_resid_mod[["residuals"]])
    
    
    df_mc <- df_mc %>% 
      mutate(resid = d_resid)
    
    results_7b <- results_7b %>% 
      rbind(weight_calc(df_mc, ref_cohort = 3))
    
    # Part 7c ----
    if (min(df_mc$rel_time) != -1) {
      rel_times_exclude <- c(-1, min(df_mc$rel_time))
    } else {
      rel_times_exclude <- c(-1, 0)
    }
    
    # Looking at the coefficient 2 periods before treatment 
    # Excluding relative times r = -4, -1
    # AS's results mean that the weights on -4, -1 (excluded times) should sum to 
    # -1, the weights on -2 should sum to 1, and the weights on all other relative
    # times should sum to zero.
    g <- -2
    rel_times_exclude <- c(-4, -1)
    
    wts_neg2 <- get_wts_dynamic(df_mc, g = g, 
                                rel_times_exclude = rel_times_exclude)
    
    results_7c <- results_7c %>% rbind(wts_neg2)
    wts_neg2 <- wts_neg2 %>% 
      left_join(trt_df, by = c("cohort", "rel_time")) %>% 
      select(-year) %>% 
      left_join(trt_df_homog, by = c("cohort", "rel_time")) %>% 
      select(-year)
    
    estims <- coef_calculate(wts_neg2)
    results_7c_coefs <- results_7c_coefs %>% rbind(estims)
    
    # Part 7d ----
    rel_times_exclude_new <- c(-4, 0)
    
    wts_neg2_new <- get_wts_dynamic(df_mc, g = g, 
                                    rel_times_exclude = rel_times_exclude_new)
    
    results_7d <- results_7d %>% rbind(wts_neg2_new)
    wts_neg2_new <- wts_neg2_new %>% 
      left_join(trt_df, by = c("cohort", "rel_time")) %>% 
      select(-year) %>% 
      left_join(trt_df_homog, by = c("cohort", "rel_time")) %>% 
      select(-year)
    
    estims_new <- coef_calculate(wts_neg2_new)
    results_7d_coefs <- results_7d_coefs %>% rbind(estims_new)
    
    # Part 7e ----
    att_r0 <- weighted_att_calc(df_mc, r = 0)
    cs_did_r0 <- cs_did(df_mc, r = 0)
    impute_did_r0 <- impute_did(df_mc, r = 0)
    
    df_estims_7e <- tibble(
      att_r0 = att_r0,
      cs_did_r0 = cs_did_r0,
      impute_did_r0 = impute_did_r0
    )
    
    results_7e <- results_7e %>% rbind(df_estims_7e)
  }
  
  # Plot results from part (a) ----
  # Nonlog 
  df_7a_plot <- results_7a %>% 
    group_by(cohort, year) %>% 
    summarise(y0_mean = mean(y0),
              y0_sd = sd(y0),
              y0_ci95up = y0_mean + qt(.975, M - 1) * y0_sd,
              y0_ci95low = y0_mean - qt(.975, M - 1) * y0_sd,
              .groups = "keep")
  
  ggplot(df_7a_plot) + 
    geom_point(aes(x = as.factor(year), y = y0_mean, 
                   color = as.factor(cohort))) + 
    geom_line(aes(group = cohort, 
                  x = as.factor(year), 
                  y = y0_mean,
                  color = as.factor(cohort))) + 
    labs(x = "Year", y = "Avg. Untreated Outcome", color = "Cohort") + 
    ggtitle("Parallel Trends for Y", subtitle = paste0(M, " Simulations")) + 
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave(here("output", "plots", paste0("plt_ct_n", N, "_m", M, ".png")), 
         bg = "white", width = 6, height = 4)
  
  # Log 
  df_7a_plot_log <- results_7a_log %>% 
    group_by(cohort, year) %>% 
    summarise(y0_mean = mean(logy0, na.rm = TRUE),
              y0_sd = sd(logy0),
              y0_ci95up = y0_mean + qt(.975, M - 1) * y0_sd,
              y0_ci95low = y0_mean - qt(.975, M - 1) * y0_sd,
              .groups = "keep")
  
  ggplot(df_7a_plot_log) + 
    geom_point(aes(x = as.factor(year), y = y0_mean, 
                   color = as.factor(cohort))) + 
    geom_line(aes(group = cohort, 
                  x = as.factor(year), 
                  y = y0_mean,
                  color = as.factor(cohort))) + 
    labs(x = "Year", y = "Avg. Untreated Outcome", color = "Cohort") + 
    ggtitle("Parallel Trends for log(Y)", subtitle = paste0(M, " Simulations")) + 
    theme_minimal() + 
    theme(legend.position = "bottom")
  
  ggsave(here("output", "plots", paste0("plt_ct_log_n", N, "_m", M, ".png")), 
         bg = "white", width = 6, height = 4)
  
  # Plot results from part b ----
  ggplot(results_7b, aes(x = wt, color = as.factor(year))) + 
    geom_density() + 
    facet_wrap(vars(cohort), nrow = 2) + 
    labs(x = "Weight", y = "Density", color = "Year") + 
    ggtitle("Cohort Weight Contributions, Static TWFE Coefficient",
            subtitle = paste0(M, " Simulations")) + 
    theme_minimal()
  
  ggsave(here("output", "plots", paste0("plt_twfe_static_n", N, "_m", M, 
                                        ".png")), bg = "white")
  
  # Plot weights in part c and show pretrend testing 
  # is problematic ----
  df_7c_final <- dynamic_wts_df_process(results_7c, trt_df, trt_df_homog)
  plot_wt_rel_time(df_7c_final, g)
  ggsave(here("output", "plots", paste0("plt_wts_7c_n", N, "_m", M, ".png")), 
         width = 6, height = 4, bg = "white")
  
  results_7c_coefs_final <- dynamic_coef_df_process(results_7c_coefs, M = M)
  coefs_plot(results_7c_coefs_final, g = g, M = M)
  ggsave(here("output", "plots", paste0("plt_coefs_7c_n", N, "_m", M, ".png")), 
         width = 6, height = 4, bg = "white")
  
  # Part 7d ----
  # Throws a matrix not invertible error when wrapped in purrr:safely()
  rel_times_exclude_one <- -1
  safe_get_wts <- safely(get_wts_dynamic)
  wts_error <- safe_get_wts(df_mc, g = g, 
                            rel_times_exclude = rel_times_exclude_one)
  print(wts_error$error)
  
  df_7d_final <- dynamic_wts_df_process(results_7d, trt_df, trt_df_homog)
  plot_wt_rel_time(df_7d_final, g)
  ggsave(here("output", "plots", paste0("plt_wts_7d_n", N, "_m", M, ".png")), 
         width = 7, height = 5, bg = "white")
  
  results_7d_coefs_final <- dynamic_coef_df_process(results_7d_coefs, M = M)
  coefs_plot(results_7d_coefs_final, g = g, M = M)
  ggsave(here("output", "plots", paste0("plt_coefs_7d_n", N, "_m", M, ".png")), 
         width = 7, height = 5, bg = "white")
  
  # Part 7e ----
  df_plot_7e <- results_7e %>% 
    pivot_longer(cols = c("cs_did_r0", "impute_did_r0"), values_to = "estim", 
                 names_to = "estimator") %>% 
    mutate(bias = estim - att_r0) %>% 
    mutate(estimator = if_else(estimator == "cs_did_r0", "Direct (CS)",
                               "Imputation"))
  
  df_biasvar_7e <- df_plot_7e %>% 
    group_by(estimator) %>% 
    summarise(bias_avg = mean(bias),
              sd = sd(estim), .groups = "keep") %>% 
    mutate(across(where(is.numeric), ~ round(.x, 3)))
  
  df_biasvar_7e %>% 
    kbl(format = "latex", caption = "Estimator Performance", 
        booktabs = TRUE) %>%
    kable_styling(latex_options = c("striped", "hold_position")) %>% 
    save_kable(here("output", paste0("tbl_estims_7e_n", N, "_m", M, ".png")))
  
  
  # Plot 
  ggplot(df_plot_7e) + 
    geom_histogram(aes(x = bias, fill = as.factor(estimator)), binwidth = .05) +
    labs(x = "Bias (Bins of 0.1)", y = "Count", fill = "Estimator") + 
    ggtitle("Bias of Estimated Cohort-Weighted ATT", 
            subtitle = paste0(M, " Simulations")) + 
    theme_minimal() 
  
  ggsave(here("output", "plots", 
              paste0("plt_estim_perf_n", N, "_m", M, ".png")), 
         bg = "white", width = 7, height = 5)
}

# Test -------------------------------------------------------------------------
test <- did_mc(N = 100, M = 200)

# Code -------------------------------------------------------------------------
# Setup ----
years <- 1:6
cohorts <- c(3:5)
probs <- rep(1/length(cohorts), times = length(cohorts))
units <- 1:100

unit_vec <- rep(units, times = length(years))
years_vec <- rep(years, each = length(units))

df_mc <- tibble(
  unit = unit_vec,
  year = years_vec
)

# Create the data ----
# Assign units to cohorts 
df_mc <- df_mc %>% 
  group_by(unit) %>% 
  mutate(cohort = sample(cohorts, size = 1, prob = probs)) %>% 
  ungroup() 

# Cohort FEs
# Gaussians are centered at positive numbers to play nice with log 
df_mc <- df_mc %>% 
  mutate(alpha = if_else(cohort == 3, rnorm(1, mean = 10, sd = 3), NA_real_)) %>% 
  mutate(alpha = if_else(cohort == 4, rnorm(1, mean = 30, sd = 3), alpha)) %>% 
  mutate(alpha = if_else(cohort == 5, rnorm(1, mean = 20, sd = 3), alpha)) %>% 
  mutate(alpha = if_else(is.na(cohort), rnorm(1, mean = 5, sd = 3), alpha))

# Time FEs: N(0, 2)
df_mc <- df_mc %>% 
  group_by(year) %>% 
  mutate(lambda = rnorm(1, mean = 0, sd = 2)) %>% 
  ungroup() 

# Unit shocks that vary by unit and time: N(0, 1)
df_mc <- df_mc %>% 
  group_by(unit, year) %>% 
  mutate(eps = rnorm(1)) %>% 
  ungroup()

# Treatment effects
# For cohort 3: c(4, 8, 10, 10)
chrt3_att <- tibble(
  year = c(3, 4, 5, 6),
  att = c(4, 8, 10, 10),
  cohort = 3
)

chrt3_att_homog <- tibble(
  year = c(3, 4, 5, 6),
  att_homog = c(4, 8, 10, 10),
  cohort = 3
)

# For cohort 4: c(NA_real_, 15, 20, 21)
chrt4_att <- tibble(
  year = c(3, 4, 5, 6),
  att = c(NA_real_, 15, 20, 21),
  cohort = 4
)

chrt4_att_homog <- tibble(
  year = c(3, 4, 5, 6),
  att_homog = c(NA_real_, 4, 8, 10),
  cohort = 4
)

# For cohort 5
chrt5_att <- tibble(
  year = c(3, 4, 5, 6),
  att = c(NA_real_, NA_real_, 0.5, 18),
  cohort = 5
)

chrt5_att_homog <- tibble(
  year = c(3, 4, 5, 6),
  att_homog = c(NA_real_, NA_real_, 4, 8),
  cohort = 5
)

trt_df <- rbind(chrt3_att, chrt4_att, chrt5_att) %>% 
  mutate(rel_time = year - cohort)
trt_df_homog <- rbind(chrt3_att_homog, chrt4_att_homog, 
                      chrt5_att_homog) %>% 
  mutate(rel_time = year - cohort)

# Bind to the main data 
df_mc <- df_mc %>% 
  left_join(trt_df, by = c("cohort", "year")) %>% 
  select(-rel_time) %>% 
  left_join(trt_df_homog, by = c("cohort", "year")) %>% 
  select(-rel_time)

# Observed Y_{it} conditional on cohort = cohort FE + time FE + treatment + 
# idiosyncratic shock.
# Also, generate dummies and relative time.
df_mc <- df_mc %>% 
  mutate(d = if_else(year >= cohort, 1, 0)) %>% 
  mutate(y = if_else(!is.na(att), alpha + lambda + d * att + eps,
                     alpha + lambda + eps)) %>%
  # outcome with homogeneous relative time treatment effects 
  mutate(y_homog = if_else(!is.na(att), alpha + lambda + d * att_homog + eps,
                           alpha + lambda + eps)) %>% 
  mutate(logy = log(y)) %>% 
  mutate(rel_time = year - cohort) %>% 
  fastDummies::dummy_cols(select_columns = c("unit", "year", "cohort", 
                                             "rel_time"))

# Part 7a ---- 
# Common trends: compare cohorts 3 and 4 from periods 1 to 2
df_sum_mc <- df_mc %>% 
  filter(year == 1 | year == 2) %>% 
  group_by(cohort, year) %>% 
  summarise(y0 = mean(y))

df_sum_log_mc <- df_mc %>% 
  filter(year == 1 | year == 2) %>% 
  group_by(cohort, year) %>% 
  summarise(logy0 = mean(logy, na.rm = TRUE))

# Plots: no common trends for log 
ggplot(df_sum_mc) + 
  geom_point(aes(x = as.factor(year), y = y0, color = as.factor(cohort))) + 
  geom_line(aes(group = cohort, 
                x = as.factor(year), 
                y = y0,
                color = as.factor(cohort)))

ggplot(df_sum_log_mc) + 
  geom_point(aes(x = as.factor(year), y = logy0, color = as.factor(cohort))) + 
  geom_line(aes(group = cohort, 
                x = as.factor(year), 
                y = logy0,
                color = as.factor(cohort)))

# Part 7b ----
unit_vars <- grep("unit_", names(df_mc), value = TRUE)[-1]
year_vars <- grep("year_", names(df_mc), value = TRUE)[-1]
d_resid_mod <- ols_run(df_mc, y = "d", covars = c(unit_vars, year_vars))
d_resid <- as.vector(d_resid_mod[["residuals"]])
  

df_mc <- df_mc %>% 
  mutate(resid = d_resid)

# Function to calculate the weights for a given cohort in a given year
weight_calc <- function(df_resid, ref_cohort) {
  props <- df_resid %>% 
    group_by(cohort) %>% 
    summarise(
      prop = n()/nrow(df_resid)
    )
  
  df_sum <- df_resid %>% 
    group_by(cohort, year) %>% 
    summarise(
      mean_resid = mean(resid)
    ) 
  
  df_out <- df_sum %>% 
    left_join(props, by = "cohort") %>% 
    ungroup() %>% 
    mutate(wt = mean_resid * prop)
  
  df_out <- df_out %>% 
    filter(cohort != ref_cohort & year != 1)
}

wts_twfe_stat <- weight_calc(df_mc, ref_cohort = 3)

# Part 7c ----
get_wts_e_r <- function(df_did, e, r, g, rel_times_exclude) {
  df_reg <- df_did %>% 
    select(cohort, starts_with("rel_time"), starts_with("unit_"),
           starts_with("year_")) %>% 
    mutate(tag_obs = if_else(cohort == e & rel_time == r, 1, 0)) %>% 
    mutate(tag_obs = replace_na(tag_obs, 0))
  
  # Exclude unit, year, and relative time dummies
  unit_vars <- grep("unit_", names(df_reg), value = TRUE)[-1]
  year_vars <- grep("year_", names(df_reg), value = TRUE)[-1]
  rel_times <- grep("rel_time_", names(df_reg), value = TRUE)
  
  rel_times_include <- rel_times[!(rel_times %in% paste0("rel_time_", 
                                                         rel_times_exclude))]
  
  print(rel_times_include)
  
  print(paste0("Cohort: ", e))
  print(paste0("Relative time: ", r))
  
  # TWFE to get the weights: indicator for cohort e and relative time r 
  # on unit and time FEs + relative time indicators
  wts_reg <- ols_run(df_reg, y = "tag_obs", covars = c(unit_vars,
                                                       year_vars,
                                                       rel_times_include))
  
  # Grab the correct coefficients 
  # +1 to account for the constant
  coefs_tag <- which(wts_reg[["var_names"]] == paste0("rel_time_", g)) + 1
  coefs_e_r <- wts_reg[["estimates"]][coefs_tag][1]
  print(coefs_e_r)
}

if (min(df_mc$rel_time) != -1) {
  rel_times_exclude <- c(-1, min(df_mc$rel_time))
} else {
  rel_times_exclude <- c(-1, 0)
}

# Function that grabs all the weights used in estimating the dynamic 
# TWFE coefficient on a single relative time indicator
# Returns a list of the weights on each cohort e at relative time r
# g = the coefficient to estimate 
get_wts_dynamic <- function(df_did, g, rel_times_exclude) {
  cohorts_times <- df_did %>% 
    distinct(cohort, rel_time) %>% 
    arrange(cohort, rel_time)
  
  print(cohorts_times)
  
  wts <- map2(
    cohorts_times$cohort,
    cohorts_times$rel_time,
    ~ as.data.frame(get_wts_e_r(df_did, e = .x, r = .y, g = g,
                  rel_times_exclude = rel_times_exclude))
  ) %>% 
    list_rbind() %>% 
    cbind(cohorts_times$cohort) %>% 
    cbind(cohorts_times$rel_time)
  
  names(wts) <- c("wt", "cohort", "rel_time")
  
  return(wts)
}

plot_wt_rel_time <- function(df_wts, r) {
  ggplot(df_wts) + 
    geom_point(aes(x = rel_time, y = wt, 
                   color = as.factor(cohort))) + 
    geom_line(aes(x = rel_time, y = wt, group = cohort,
                  color = as.factor(cohort))) + 
    labs(x = "Relative Time", y = "Weight", color = "Cohort") + 
    ggtitle(paste0("Weights on TWFE Coef"), 
            subtitle = paste0("Relative Time ", r)) + 
    theme_minimal()
}

coef_calculate <- function(df_wts) {
  df_wts <- df_wts %>% 
    mutate(att_contrib = att * wt) %>% 
    mutate(att_homog_contrib = att_homog * wt) %>% 
    summarise(b_att = sum(att_contrib, na.rm = TRUE),
              b_att_homog = sum(att_homog_contrib, na.rm = TRUE))
  
}

coefs_plot <- function(df_coefs) {
  df_plot <- df_coefs %>% 
    pivot_longer(cols = c("b_att", "b_att_homog"), names_to = "var",
                 values_to = "coef")
  
  ggplot(df_plot) + 
    geom_point(aes(x = var, y = coef), color = "red") + 
    theme_minimal() + 
    theme(legend.position = "none")
}

# Looking at the coefficient 2 periods before treatment 
# Excluding relative times r = -4, -1
# AS's results mean that the weights on -4, -1 (excluded times) should sum to 
# -1, the weights on -2 should sum to 1, and the weights on all other relative
# times should sum to zero.
g <- -2
rel_times_exclude <- c(-4, -1)

wts_neg2 <- get_wts_dynamic(df_mc, g = g, rel_times_exclude = rel_times_exclude)
wts_neg2 <- wts_neg2 %>% 
  left_join(trt_df, by = c("cohort", "rel_time")) %>% 
  select(-year) %>% 
  left_join(trt_df_homog, by = c("cohort", "rel_time")) %>% 
  select(-year)

plot_wt_rel_time(wts_neg2, g)
estims <- coef_calculate(wts_neg2)
coefs_plot(estims)

# Part 7d ----
# Throws a matrix not invertible error when wrapped in purrr:safely()
rel_times_exclude_one <- -1
safe_get_wts <- safely(get_wts_dynamic)
wts_error <- safe_get_wts(df_mc, g = g, 
                         rel_times_exclude = rel_times_exclude_one)
print(wts_error$error)

# Now, try leaving out a period with non-zero treatment effects, e.g.
# relative time -4 combined with 0
rel_times_new <- c(-4, 0)

wts_neg2_new <- get_wts_dynamic(df_mc, g = g, 
                                rel_times_exclude = rel_times_new)
wts_neg2_new <- wts_neg2_new %>% 
  left_join(trt_df, by = c("cohort", "rel_time")) %>% 
  select(-year) %>% 
  left_join(trt_df_homog, by = c("cohort", "rel_time")) %>% 
  select(-year)

plot_wt_rel_time(wts_neg2_new, g)
estims_new <- coef_calculate(wts_neg2_new)
coefs_plot(estims_new)

# Part 7e ----
# Cohort-weighted ATT 
weighted_att_calc <- function(df, r) {
  df_att <- df %>% 
    filter(rel_time == r) %>% 
    group_by(cohort) %>% 
    summarise(size = n(),
              att = att) %>% 
    slice(1) %>% 
    ungroup() %>% 
    mutate(pr = size/sum(size))
  
  df_att <- df_att %>% 
    mutate(att = replace_na(att, 0)) %>% 
    mutate(wt_att_contrib = att * pr)

  att <- sum(df_att$wt_att_contrib)
}

att_r0 <- weighted_att_calc(df_mc, r = 0)

# CS approach 
# Function that calculates ATT(g, t)
cs_gt_calc <- function(df, e, r) {
  # Setup for r >= 0 
  if (r >= 0) {
    df_chrt <- df %>% 
      filter(cohort == e) %>% 
      filter(year == e + r | year == e - 1) %>% 
      mutate(pre = if_else(year == e - 1, 1, 0))
    
    tag_not_yet <- df %>% 
      filter(year == e + r & d == 0) 
    
    chrts_not_yet <- unique(tag_not_yet$cohort)
    
    df_not_yet <- df %>% 
      filter(cohort %in% chrts_not_yet) %>% 
      filter(year == e + r | year == e - 1) %>% 
      mutate(pre = if_else(year == e - 1, 1, 0))
    
    if (nrow(df_not_yet) == 0) {
      return(NA_real_)
    } 
  } else {
    df_chrt <- df %>% 
      filter(cohort == e) %>% 
      filter(year == e + r | year == e + r - 1) %>% 
      mutate(pre = if_else(year == e +r - 1, 1, 0))
    
    tag_not_yet <- df %>% 
      filter(year == e + r & d == 0) 
    
    chrts_not_yet <- unique(tag_not_yet$cohort)
    
    df_not_yet <- df %>% 
      filter(cohort %in% chrts_not_yet) %>% 
      filter(year == e + r | year == e + r - 1) %>% 
      mutate(pre = if_else(year == e + r - 1, 1, 0))
    
    if (nrow(df_not_yet) == 0) {
      return(NA_real_)
    }
  }
  
  # Quantities 
  chrt_post <- mean(df_chrt[df_chrt[, "pre"] == 0, ]$y, na.rm = TRUE)
  chrt_pre <- mean(df_chrt[df_chrt[, "pre"] == 1, ]$y, na.rm = TRUE)
  not_yet_post <- mean(df_not_yet[df_not_yet[, "pre"] == 0, ]$y, na.rm = TRUE)
  not_yet_pre <- mean(df_not_yet[df_not_yet[, "pre"] == 1, ]$y, na.rm = TRUE)
  
  att_gt <- (chrt_post - chrt_pre) - (not_yet_post - not_yet_pre)
}

# Function that aggregates ATT(g, t) for a relative time 
cs_did <- function(df, r) {
  
  df_att <- df %>% 
    filter(rel_time == r) %>% 
    group_by(cohort) %>% 
    summarise(size = n()) 
  
  chts <- unique(df_att$cohort)
  
  # Run the main CS estimator
  att_gt <- chts %>% 
    map(~ cs_gt_calc(df, e = .x, r = r))
  
  df_att_gt <- tibble(
    cohort = chts,
    att_gt = as.numeric(att_gt)
  )
  
  # Left bind to the main df and calculate final estimate, dropping groups 
  # with no control comparisons
  df_att <- df_att %>% 
    left_join(df_att_gt, by = "cohort") %>% 
    mutate(pr = size/sum(size)) %>% 
    mutate(att_contrib = pr * att_gt)
  
  att_e <- sum(df_att$att_contrib, na.rm = TRUE)
}


cs_did_r0 <- cs_did(df_mc, r = 0)

# Imputation approach 
impute_did <- function(df, r, rel_times_exclude = -1) {
  
  # Select untreated subsample and regress on cohort + time FEs
  df_d0 <- df %>% 
    filter(d == 0)
  
  cohorts_include <- paste0("cohort_", unique(df_d0$cohort)[-1])
  years_include <- paste0("year_", unique(df_d0$year)[-1])
  
  mod_d0 <- ols_run(df = df_d0, y = "y", 
                    covars = c(cohorts_include, years_include))
  
  estims_d0 <- mod_d0[["estimates"]]
  print(estims_d0)
  
  mat_predict <- matrixPrep(df, cols = c(cohorts_include, years_include),
                            constant = TRUE)
  
  # Get fitted values and subtract them from each actual outcome to construct
  # Y dot
  y0_fit <- mat_predict %*% estims_d0
  df <- df %>% 
    cbind(y0_fit = y0_fit) %>% 
    mutate(y_dot = y - y0_fit)
  
  # Regress Y dot on the relative time indicators.
  rel_times <- grep("rel_time_", names(df), value = TRUE)
  
  rel_times_include <- rel_times[!(rel_times %in% paste0("rel_time_", 
                                                         rel_times_exclude))]
  
  impute_mod <- ols_run(df, y = "y_dot", 
                        covars = rel_times_include)
  print(impute_mod[["var_names"]])
  print(impute_mod[["estimates"]])
  
  # Grab the right coefficient 
  coef_tag <- which(impute_mod[["var_names"]] == paste0("rel_time_", r)) + 1
  att_impute <- impute_mod[["estimates"]][coef_tag]
}

test <- impute_did(df_mc, r = 0)
