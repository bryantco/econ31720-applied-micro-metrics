# Auxiliary functions ----------------------------------------------------------
# Auxiliary ----
matrixPrep <- function(df, cols, constant = FALSE) {
  df_out <- df %>% select(all_of(cols))
  
  mat_out <- as.matrix(df_out)
  colnames(mat_out) <- NULL 
  rownames(mat_out) <- NULL 
  
  if (constant) {
    mat_out <- addConstant(mat_out)
  }
  
  mat_out
}

addConstant <- function(mat) {
  ones <- matrix(rep(1, nrow(mat)), nrow = nrow(mat))
  mat <- cbind(ones, mat)
}

# Computational functions ------------------------------------------------------
# Estimation ----
# Least squares 
ols_fit <- function(X, y) {
  Xt <- t(X)
  XtX <- Xt %*% X
  Xty <- Xt %*% y
  b <- solve(XtX) %*% Xty
}

# Logit helpers ----
# Sigmoid function
sigmoid <- function(z) {
  (exp(z))/(1 + exp(z))
}

# Logit log likelihood 
# beta is the k x 1 parameter vector 
# X is assumed to be N x k
logit.lik <- function(beta, y, X) {
  sig_Xb <- sigmoid(X %*% beta)
  ll_vec <- y * log(sig_Xb) + (1 - y) * log(1 - sig_Xb)
  nll <- -sum(ll_vec) # because stats::optim() performs minimization
}

# Matching helpers ---- 
# Calculates inverse variance weighting matrix for input matrix X
v_calc <- function(X) {
  # Handle the case where X is only one column 
  if (dim(X)[2] == 1) {
    V <- 1/diag(var(X))
  } else {
    V <- diag(1/diag(var(X))) 
  }
  
  return(V)
}

# Calculates the distance based on a supplied matrix V 
# x, z assumed to be vectors (that as.matrix() converts to column vectors)
distance_calc <- function(x, z, V) {
  x <- as.matrix(x)
  z <- as.matrix(z)
  dist <- (t(x-z) %*% V %*% (x-z))^(1/2)
  dist <- dist[1, 1]
}


# Inference ---- 
# Clustered SEs
se_clust <- function(df, ols_fit, clust_var) {
  # Grab useful values 
  df <- df[ols_fit$tag_obs, ] # grab obs. used in estimation
  clust_var_values <- unique(df[, clust_var])
  m <- length(clust_var_values) # number of clusters 
  n <- ols_fit$n_obs
  k <- ols_fit$k
  
  meat_list <- lapply(clust_var_values, function(i) {
    indices <- which(df[, clust_var] == i)
    
    # Select the X matrix and residuals 
    Xi <- matrixPrep(df[indices, ], cols = ols_fit$var_names,
                     constant = ols_fit$constant)
    
    eta_i <- as.matrix(ols_fit$residuals[indices], ncol = 1)
    
    # Calculate the meat for each cluster
    XiTXi <- t(Xi) %*% eta_i %*% t(eta_i) %*% Xi
  })
  
  meat <- Reduce("+", meat_list)
  X <- ols_fit$design_mat 
  bread <- solve(t(X) %*% X)
  
  # Final result: only output the SE, making Stata finite sample correction 
  vcov <- (m/(m-1)) * ((n-1)/(n-k)) * (bread %*% meat %*% bread)
  se_out <- sqrt(diag(vcov))
}

# Bootstrap 
se_boot <- function(fun_results, seed = 12345, B = 50, var_by) {
  # Extract the call, data used, and number of obs. from the model object
  call <- fun_results[["call"]]
  tag_obs <- fun_results[["tag_obs"]]
  df <- rlang::eval_tidy(call[["df"]])[tag_obs, ]
  n <- fun_results[["n_obs"]]
  
  # Handles proportions in creating the bootstrap sample 
  levels <- df %>% 
    dplyr::select(!!sym(var_by)) %>% 
    dplyr::distinct() %>% 
    dplyr::pull() 
  
  var_values <- df %>% 
    dplyr::select(!!sym(var_by)) %>% 
    dplyr::pull() 
  
  obs_labels <- match(var_values, levels)
  probs <- levels %>% 
    purrr::map_dbl(~ sum(var_values == .x))
  
  probs <- probs/length(var_values)
  
  obs_probs <- obs_labels %>% 
    purrr::map_dbl(~ probs[.x])
  
  idx <- 1:nrow(df) # indices 
  
  # Bootstrap the estimates 
  set.seed(seed)
  
  estims_bstrap <- 1:B %>% purrr::map(function(x) {
    # Sample rows 
    idx_sample <- sample(idx, size = n, replace = TRUE, prob = obs_probs)
    # df_boot <- dplyr::slice_sample(df, n = n, replace = TRUE)
    df_boot <- df[idx_sample, ]
    call_new <- call
    call_new$df <- df_boot
    fn_name <- rlang::as_string(call[[1]]) # function name 
    estim <- do.call(fn_name, as.list(call_new)[-1])$estimates
  })
  
  estims_bstrap <- estims_bstrap %>% 
    purrr::reduce(rbind)
  
  se_bstrap <- sqrt(diag(var(estims_bstrap)))
}

# Routines ---------------------------------------------------------------------
# Function arguments are consistently coded to play nice with other functions 
# eg bootstrapping 
ols_run <- function(df, y, covars, constant = TRUE) {
  
  # Select nonempty rows and prep 
  df <- df[, c(covars, y)]
  tag_obs <- which(rowSums(is.na(df))==0) # tag nonmissing obs.
  df <- na.omit(df)
  n <- nrow(df)
  k <- length(covars) + constant
  X <- matrixPrep(df, covars, constant = constant)
  y <- matrixPrep(df, y)
  
  
  # Fit the model 
  b <- ols_fit(X, y)
  coefs <- c(t(b[, 1])) # Convert matrix coefs to vector
  
  # Get the fitted values: Xb
  fitted_values <- X %*% b
  
  # Get the residuals: y - Xb 
  residuals <- y - fitted_values
  
  # R-squared and adjusted R-squared
  sst <- sum((y - mean(y))^2)
  ssr <- sum(residuals^2)
  r_sq <- 1 - ssr/sst
  r_sq_adjust <- 1 - (1 - r_sq) * ((n - 1)/(n - k))
  
  # Final output
  ols_out <- list(
    design_mat = X,
    call = match.call(),
    estimates = coefs,
    constant = constant,
    n_obs = n,
    k = k,
    tag_obs = tag_obs,
    var_names = covars,
    fitted_values = fitted_values,
    residuals = residuals,
    r_sq = r_sq,
    r_sq_adjust = r_sq_adjust
  )
}

# NN Matching for the ATT
find_treated_nn <- function(df, y, k, covars, treat_var, 
                            estimand = "att") {
  
  # Setup: grab matrix of covariates, both in its entirety and separated by 
  # treated/untreated, and treated/untreated outcomes 
  estim <- NULL 
  
  df <- df %>% dplyr::select(!!sym(y), !!sym(treat_var), all_of(covars))
  tag_obs <- which(rowSums(is.na(df))==0) # tag nonmissing obs.
  df <- na.omit(df)
  
  df_covars <- df %>% dplyr::select(all_of(covars))
  
  df_trt <- df %>% 
    dplyr::filter(!!sym(treat_var) == 1)
  
  df_cntrl <- df %>% 
    dplyr::filter(!!sym(treat_var) == 0) 
  
  y1 <- df_trt %>% 
    dplyr::select(!!sym(y)) %>% 
    unlist() %>% 
    unname()
  
  y0 <- df_cntrl %>% 
    dplyr::select(!!sym(y)) %>% 
    unlist() %>% 
    unname() 
  
  x1 <- df_trt %>% dplyr::select(all_of(covars)) 
  x0 <- df_cntrl %>% dplyr::select(all_of(covars)) 
  
  # Calculation begins 
  # First, find the inverse variance weighting matrix. Then, if the 
  # estimand is the ATT, calculate the distance matrix as an N1 x N0 matrix of 
  # matches. 
  V <- v_calc(as.matrix(df_covars))
  
  if (estimand == "att" | estimand == "ate") {
    # Match treated obs. to untreated obs.
    distances <- x0 %>%
      apply(., 1, function(x) {
        apply(x1, 1, function(y) distance_calc(y, x, V))
      })
    
    # Find the top k values in each row of the distance matrix 
    top_k <- list() 
    
    for (row in 1:nrow(distances)) {
      top_k[[row]] <- sort(unique(distances[row, ]))
    }
    
    top_k <- top_k %>% 
      map(function(x) x[1:k])
    
    top_k_indx <- 1:nrow(distances) %>% 
      purrr::map(., ~ which(distances[.x, ] %in% top_k[[.x]]))
    
    # Now, compute averages of the untreated outcomes based on the indices
    y0_use <- top_k_indx %>% 
      map(., ~ y0[.x]) %>% 
      map_dbl(., ~ mean(.x))
    
  } 
  
  if (estimand == "atu" | estimand == "ate") {
    # Match untreated obs. to treated obs.
    distances <- x1 %>%
      apply(., 1, function(x) {
        apply(x0, 1, function(y) distance_calc(y, x, V))
      })
    
    # Find the top k values in each row of the distance matrix 
    top_k <- list() 
    
    for (row in 1:nrow(distances)) {
      top_k[[row]] <- sort(unique(distances[row, ]))
    }
    
    top_k <- top_k %>% 
      map(function(x) x[1:k])
    
    top_k_indx <- 1:nrow(distances) %>% 
      purrr::map(., ~ which(distances[.x, ] %in% top_k[[.x]]))
    
    # Now, compute averages of the treated outcomes based on the indices
    y1_use <- top_k_indx %>% 
      map(., ~ y1[.x]) %>% 
      map_dbl(., ~ mean(.x))
  }
  
  # Calculate treatment effects 
  # For the ATT, only impute y0. For the ATU, only impute y1. For the ATE,
  # combine imputations. 
  if (estimand == "att") {
    estim <- mean(y1 - y0_use)
  } else if (estimand == "atu") {
    estim <- mean(y0 - y1_use)
  } else if (estimand == "ate") {
    y1_use <- c(y1, y1_use)
    y0_use <- c(y0, y0_use)
    
    estim <- mean(y1_use - y0_use)
  }
  
  n <- length(y0) + length(y1)
  
  nn_out <- list(
    estimates = estim,
    call = match.call(),
    n_obs = n,
    tag_obs = tag_obs
  )
}

# Logit routine
logit_run <- function(df, y, covars, constant = TRUE) {
  
  # Select nonempty rows and prep 
  df <- df[, c(covars, y)]
  tag_obs <- which(rowSums(is.na(df))==0) # tag nonmissing obs.
  df <- na.omit(df)
  n <- nrow(df)
  k <- length(covars) + constant
  X <- matrixPrep(df, covars, constant = constant)
  y <- matrixPrep(df, y)
  
  
  # Fit the model 
  b0 <- rep(0, times = k) # initialize starting values as zero 
  b <- optim(b0, logit.lik, y = y, X = X, method = "BFGS")[["par"]]
  
  # Get the fitted values: Xb
  fitted_values <- sigmoid(X %*% b)
  
  # Get the residuals: y - Xb 
  residuals <- y - fitted_values
  
  # Final output
  logit_out <- list(
    call = match.call(),
    estimates = b,
    constant = constant,
    n_obs = n,
    k = k,
    tag_obs = tag_obs,
    var_names = covars,
    fitted_values = fitted_values,
    residuals = residuals
  )
}

# Takes in a df and a logit model, the result of logit_run(), 
# and returns the same df but with the fitted values (propensity scores)
# and the fitted values split into nbins number of bins
get_df_blocked <- function(df, y, treat_var, logit_mod, nbins) {
  # Pull out useful values
  tag_obs <- logit_mod[["tag_obs"]]
  covars <- logit_mod[["var_names"]]
  yhat <- logit_mod[["fitted_values"]]
  
  df <- df[tag_obs, ][, c(y, treat_var, covars)]
  block <- cut(yhat, 
               breaks = nbins, 
               labels = FALSE)
  
  df_out <- cbind(df, yhat, block) %>% 
    dplyr::rename(px = yhat)
}

plot_df_blocked <- function(df, y, treat_var) {
  ggplot2::ggplot(df) +
    ggplot2::geom_density(
      aes(
        x = !!sym(y),
        fill = as.factor(!!sym(treat_var)),
        color = as.factor(!!sym(treat_var))
      ),
      position = "identity",
      alpha = 0.5
    ) + 
    labs(x = y, y = "Density", fill = treat_var) + 
    ggplot2::facet_wrap(vars(block)) + 
    ggplot2::guides(color = "none")
}

# DID Monte Carlo helper functions ----
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