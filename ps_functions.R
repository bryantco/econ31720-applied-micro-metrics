################################################################################
# AUXILIARY FUNCTIONS 
################################################################################
matrixPrep <- function(df, cols, constant = FALSE) {
  df_out <- df %>% dplyr::select(all_of(cols))
  
  mat_out <- as.matrix(df_out)
  colnames(mat_out) <- NULL 
  rownames(mat_out) <- NULL 
  
  if (constant) {
    mat_out <- addConstant(mat_out)
  }

  mat_out
}

addConstant <- function(mat) {
  mat <- cbind(matrix(rep(1, nrow(mat)), nrow = nrow(mat)), mat)
}

################################################################################
# COMPUTATIONAL FUNCTIONS 
################################################################################
# Estimation ----
# Least squares 
ols_fit <- function(X, y) {
  XtX <- crossprod(X, X)
  Xty <- crossprod(X, y)
  b <- solve(XtX, Xty, tol = 1e-19)
}

# IV estimator 
iv_fit <- function(Z, X, y) {
  ZtX <- crossprod(Z, X)
  Zty <- crossprod(Z, y)
  b <- solve(ZtX, Zty, tol = 1e-19)
}

# TSLS
tsls_fit <- function(ZtX, Zty, ZtZ) {
  denom <- t(ZtX) %*% solve(ZtZ) %*% ZtX
  num <- t(ZtX) %*% solve(ZtZ) %*% Zty
  b <- solve(denom, num, tol = 1e-19)
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
# F-test of all the slope coefficients from OLS 
test_f_all <- function(ols_fit, df, ols_vcov) {
  J <- ols_fit$k - 1
  zeros <- matrix(rep(0, J), nrow = J)
  R <- cbind(zeros, diag(J))
  b2 <- as.matrix(ols_fit$estimates[-1], nrow = J)
  
  # Get the F-stat
  f <- (t(b2) %*% solve(R %*% ols_vcov %*% t(R), tol = 1e-30) %*% b2)/J
  f <- as.numeric(f)
  
  # Conduct the F-test 
  # The statistic JF is asymptotic chi-squared with J degrees of freedom 
  p <- pchisq(f, df = J, lower.tail = FALSE)
  
  # Return output 
  f_out <- list(
    f = f,
    p = p
  )
}

# Normal OLS SEs 
se_ols <- function(ols_fit) {
  s2 <- ols_fit$ssr/(ols_fit$n_obs - ols_fit$k)
  X <- ols_fit$design_mat
  
  vcov <- s2 * solve(crossprod(X, X))
  se <- sqrt(diag(vcov))
  
  return(se)
}

# Clustered SEs
se_clust <- function(df, lin_fit, clust_var, fe_list = NULL,
                     vcov_return = FALSE, ssc = TRUE) {
  
  # Grab useful values 
  df <- df[lin_fit$tag_obs, ] # grab obs. used in estimation
  clust_var_values <- unique(df[, clust_var])
  m <- length(clust_var_values) # number of clusters 
  n <- lin_fit$n_obs
  k <- lin_fit$k
  meat <- matrix(0, nrow = k, ncol = k)
  j <- 1
  
  for (i in clust_var_values) {
    indices <- which(df[, clust_var] == i)
    Xi <- matrix(lin_fit[["design_mat"]][indices, ], nrow = length(indices))
      
    eta_i <- as.matrix(lin_fit$residuals[indices], ncol = 1)
    XiTXi <- t(Xi) %*% eta_i %*% t(eta_i) %*% Xi
      
    meat <- meat + XiTXi
    j <- j + 1
  }
  
  
  X <- lin_fit$design_mat 
  bread <- solve(t(X) %*% X)

  # Final result: only output the SE, making Stata finite sample correction
  # In the finite sample correction, add back any nested fixed effects to 
  # prevent them from being subtracted off. 
  k_eff <- k
  if (!is.null(fe_list)) {
    k_eff <- k - length(fe_list) - lin_fit$constant - 1
  }
  
  if (ssc) {
    vcov <- (m/(m-1)) * ((n-1)/(n-k_eff)) * (bread %*% meat %*% bread)
  } else {
    vcov <- bread %*% meat %*% bread
  }
  
  # If user requested full variance-covariance matrix, return it. Else,
  # only return the standard errors
  if (vcov_return) {
    return(vcov)
  } else {
    return(sqrt(diag(vcov)))
  }
}

# Bootstrap ----
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

se_boot_block <- function(fun_results, seed = 12345, B = 50, block) {
  # Extract the call, data used, and number of obs. from the model object
  call <- fun_results[["call"]]
  tag_obs <- fun_results[["tag_obs"]]
  df <- rlang::eval_tidy(call[["df"]])[tag_obs, ]
  
  # Count block sizes 
  levels <- df %>% 
    dplyr::select(!!sym(block)) %>% 
    dplyr::distinct() %>% 
    dplyr::pull() 
  
  var_values <- df %>% 
    dplyr::select(!!sym(block)) %>% 
    dplyr::pull() 
  
  obs_labels <- match(var_values, levels)
  
  probs <- levels %>% 
    purrr::map_dbl(~ sum(var_values == .x))
  
  probs <- probs/length(var_values)
  
  df <- df %>% 
    cbind(obs_labels)
  
  print(paste0("Number of blocks: ", length(levels)))
  idx <- 1:length(levels) # indices: number of blocks 
  
  # Bootstrap the estimates 
  set.seed(seed)
  
  estims_bstrap <- 1:B %>% purrr::map(function(x) {
    # Sample blocks 
    idx_sample <- sample(idx, size = length(levels), replace = TRUE, 
                         prob = probs)
    df_boot <- df[df$obs_labels %in% idx_sample, ]
    call_new <- call
    call_new$df <- df_boot
    fn_name <- rlang::as_string(call[[1]]) # function name 
    estim <- do.call(fn_name, as.list(call_new)[-1])$estimates
  })
  
  estims_bstrap <- estims_bstrap %>% 
    purrr::reduce(rbind)
  
  se_bstrap <- sqrt(diag(var(estims_bstrap, na.rm = TRUE)))
}

# Anderson-Rubin ----
# Two-tailed t-test
test_t <- function(ols_fit, var_name, b0 = 0, se_vec) {
  if(!var_name %in% ols_fit$var_names) {
    stop("Error! Variable not found in the regression.")
  }
  
  if (var_name == "constant") {
    pos <- 1
  } else {
    pos <- ols_fit$constant + which(ols_fit$var_names == var_name)
  }
  
  # Grab the correct standard error and estimate 
  se <- se_vec[pos]
  b <- ols_fit$estimates[pos]
  
  # Construct the t-statistic 
  t <- (b - b0)/se 
  df <- ols_fit$n_obs - ols_fit$k
  p_val <- 2 * stats::pt(abs(t), lower.tail = FALSE, df = df)
  
  t_out <- list(
    t = t,
    p_val = p_val
  )
}

beta_grid <- seq(-5, 5, by = .01)

test_ar <- function(tsls_fit, df, d_var, z_var,
                    beta_grid, se_meth, se_opts, alpha) {
  # Set up the confidence set 
  conf_set <- c()
  p_vals <- c()
  
  # Grab needed quantities 
  y <- tsls_fit$y
  X <- tsls_fit$X_mat
  D <- df[tsls_fit$tag_obs, ] %>% 
    select(!!sym(d_var)) %>% 
    as.matrix(ncol = 1) 
  
  # Loop over the grid of betas 
  for (b in beta_grid) {
    print(paste0("Testing the beta value ... ", b))
    # Grab the residuals 
    residuals_null <- y - D*b
    
    # Fit the AR regression: regress the residuals onto the first stage
    df_temp <- df[tsls_fit$tag_obs, ] %>% 
      cbind(residuals_null) %>% 
      as.data.frame()
    
    reg_ar <- ols_run(df_temp, y = "residuals_null", 
                      covars = tsls_fit$var_names_fs)
    
    if (se_meth == "clust") {
      se_args <- list(
        df = df_temp,
        lin_fit = reg_ar
      )
      
      se_args <- c(se_args, se_opts)
      reg_ar_se <- do.call("se_clust", args = se_args)
    }
    
    # Compute the t-test on the scalar instrument for the scalar endogeneous
    # variable 
    # If we fail to reject the null hypothesis that beta = the current value,
    # then add the current value to the AR confidence set.
    t <- test_t(reg_ar, var_name = z_var, se_vec = reg_ar_se)
    p_vals <- c(p_vals, t$p_val)
    
    if (t$p_val <= alpha) {
      conf_set <- c(conf_set, "out")
    } else {
      conf_set <- c(conf_set, "in")
    }
  }
  
  results <- tibble(
    beta = beta_grid,
    conf_set = conf_set,
    p_vals = p_vals
  )
}


################################################################################
# ROUTINES 
################################################################################
# Function arguments are consistently coded to play nice with other functions 
# eg bootstrapping 
# If doing a constant-only regression, specify y and covars to be the same
# outcome variable string, and set constant = FALSE.
ols_run <- function(df, y, covars, constant = TRUE, wt = NULL) {
  
  # Select nonempty rows and prep 
  if (!is.null(wt)) {
    df <- df[, c(covars, y, wt)]
  } else {
    df <- df[, c(covars, y)]
  }

  tag_obs <- which(rowSums(is.na(df))==0) # tag nonmissing obs.
  df <- na.omit(df)
  
  if (!is.null(wt)) {
    wt <- matrixPrep(df, cols = wt) %>% as.vector()
  } 
  
  df <- df[, c(covars, y)]
  n <- nrow(df)
  k <- length(covars) + constant
  
  # Prep X
  # Handle the case where a regression only includes an intercept 
  if (length(unique(c(y, covars))) > 1) {
    X_unweight <- matrixPrep(df, covars, constant = constant) 
  } else {
    # only matrix of ones if user requested a constant-only regression 
    X_unweight <- matrix(rep(1, n), nrow = n)
    
  }
  y_unweight <- matrixPrep(df, y) 
  
  X <- X_unweight
  y <- y_unweight
  
  if (!is.null(wt)) {
    X <- X_unweight * sqrt(wt)
    y <- y_unweight * sqrt(wt)
  } 
  
  # Fit the model 
  b <- ols_fit(X, y)
  coefs <- c(t(b[, 1])) # Convert matrix coefs to vector
  
  # Get the fitted values: Xb
  fitted_values <- X %*% b
  fitted_values_unweight <- X_unweight %*% b
  
  # Get the residuals: y - Xb 
  residuals <- y - fitted_values
  residuals_unweight <- y_unweight - fitted_values_unweight
  
  # Sum of square residuals 
  ssr <- sum(residuals^2)
  
  # R-squared and adjusted R-squared
  sst <- sum((y - mean(y))^2)
  ssr <- sum(residuals^2)
  r_sq <- 1 - ssr/sst
  r_sq_adjust <- 1 - (1 - r_sq) * ((n - 1)/(n - k))
  
  # Final output
  ols_out <- list(
    design_mat = X,
    X_unweight = X_unweight,
    y_unweight = y_unweight,
    call = match.call(),
    estimates = coefs,
    constant = constant,
    n_obs = n,
    k = k,
    tag_obs = tag_obs,
    var_names = covars,
    fitted_values = fitted_values,
    residuals = residuals,
    residuals_unweight = residuals_unweight,
    ssr = ssr,
    r_sq = r_sq,
    r_sq_adjust = r_sq_adjust
  )
}

# TSLS
tsls_run <- function(df, y, x_vars, z_list, constant = TRUE) {
  # Select nonempty rows and prep 
  vars_endo <- intersect(x_vars, names(z_list)) # endogeneous variables
  instruments <- z_list %>% unlist() %>% unname() # instruments 
  z_vars <- c(x_vars[!x_vars %in% vars_endo], instruments)
  all_vars <- union(x_vars, instruments)
  
  # Grab names 
  vars_non_endo_names <- setdiff(x_vars, vars_endo)
  var_names_fs <- c(vars_non_endo_names, instruments)
  
  # Select nonempty rows and prep 
  df <- df[, c(all_vars, y)]
  
  tag_obs <- which(rowSums(is.na(df))==0) # tag nonmissing obs.
  df <- na.omit(df)
  n <- nrow(df)
  k <- length(x_vars) + constant
  
  # Prep X, Z, and y
  X <- matrixPrep(df, x_vars, constant = constant) 
  Z <- matrixPrep(df, z_vars, constant = constant)
  y <- matrixPrep(df, y) 
  
  # Calculate the ingredients of the TSLS fit 
  ZtX <- crossprod(Z, X)
  Zty <- crossprod(Z, y)
  ZtZ <- crossprod(Z, Z)
  
  # Save the design matrix for later: PzX = Z(Z'Z)^{-1}Z'X
  design_mat <- Z %*% solve(ZtZ) %*% ZtX
  
  # Fit the model 
  b <- tsls_fit(ZtX, Zty, ZtZ)
  coefs <- c(t(b[, 1])) # Convert matrix coefs to vector
  
  # Get the residuals as y - Xb, NOT y - Xhat %*% b
  residuals <- y - X %*% b
  
  # Final output
  tsls_out <- list(
    design_mat = design_mat,
    y = y,
    Z_mat = Z,
    X_mat = X,
    call = match.call(),
    estimates = coefs,
    constant = constant,
    n_obs = n,
    k = k,
    tag_obs = tag_obs,
    var_names = x_vars,
    var_names_fs = var_names_fs,
    residuals = residuals
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
    estim <- mean(y1_use - y0)
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

# cs_did routine 
cs_did_run <- function(df, rel_time, covars = NULL) {
  if (!is.null(covars)) {
    df <- df[, c(covars, "y", "cohort", "rel_time", "year", "d")]
    df <- na.omit(df)
    n <- nrow(df)
    tag_obs <- which(rowSums(is.na(df))==0) # tag nonmissing obs.
  } else {
    df <- df[, c("y", "cohort", "rel_time", "year", "d")]
    df <- na.omit(df)
    n <- nrow(df)
    tag_obs <- which(rowSums(is.na(df))==0) # tag nonmissing obs.
  } 
  
  cs_did_results <- cs_did(df, rel_time, covars)
  
  # Final output
  cs_did_out <- list(
    call = match.call(),
    estimates = cs_did_results,
    n_obs = n,
    tag_obs = tag_obs
  )
}

# DID imputation routine 
impute_did_run <- function(df, r, rel_times_exclude = -1, covars = NULL) {
  if (!is.null(covars)) {
    df <- df %>% 
      dplyr::select(y, d, starts_with("cohort_"), starts_with("year_"),
             starts_with("rel_time_"), cohort, year, all_of(covars))
   df <- na.omit(df)
   n <- nrow(df)
   tag_obs <- which(rowSums(is.na(df))==0) # tag nonmissing obs.
  } else {
    df <- df %>% 
      dplyr::select(y, d, starts_with("cohort_"), starts_with("year_"),
             starts_with("rel_time_"), cohort, year)
    df <- na.omit(df)
    n <- nrow(df)
    tag_obs <- which(rowSums(is.na(df))==0) # tag nonmissing obs.
  }

  impute_did_results <- impute_did(df, r, covars = covars)
  
  # Final output
  impute_did_out <- list(
    call = match.call(),
    estimates = impute_did_results,
    n_obs = n,
    tag_obs = tag_obs
  )
}

# Synthetic control routines ----
# Returns x1 as a K x 1 matrix and x0 as a K x N0 matrix 
data_prep_sc <- function(df, treat_var, treat_unit, treat_time_var, 
                         treat_time, covars, covars_lag = NULL,
                         outcome_var) {
  
  # Complete the data 
  df <- df %>% 
    tidyr::complete(!!sym(treat_var), !!sym(treat_time_var))
  
  # Handle lags. If no lags are requested, take the mean of every covariate 
  # across the pre-period. 
  if (is.null(covars_lag)) {
    # Grab covariates 
    x1_mat <- df %>% 
      dplyr::filter(!!sym(treat_var) == treat_unit) %>% 
      dplyr::filter(!!sym(treat_time_var) < treat_time) %>% 
      dplyr::summarise(across(all_of(covars), ~ mean(.x, na.rm = TRUE))) %>% 
      as.matrix() %>% 
      t()
    
    x0_mat <- df %>% 
      dplyr::filter(!!sym(treat_var) != treat_unit) %>% 
      dplyr::filter(!!sym(treat_time_var) < treat_time) %>% 
      dplyr::group_by(!!sym(treat_var)) %>% 
      dplyr::summarise(across(all_of(covars), ~ mean(.x, na.rm = TRUE))) 
    
    tag_obs_x0 <- which(rowSums(is.na(x0_mat))==0) # tag nonmissing obs.
    tag_obs_x0_names <- x0_mat %>% 
      na.omit() %>% 
      dplyr::distinct(!!sym(treat_var)) %>% 
      dplyr::pull() 
    
    x0_mat <- x0_mat %>% 
      na.omit() %>% 
      as.matrix() %>% 
      t() 
    
    # Grab outcomes 
    z1_mat <- df %>% 
      dplyr::filter(!!sym(treat_var) == treat_unit) %>% 
      dplyr::filter(!!sym(treat_time_var) < treat_time) %>% 
      dplyr::select(!!sym(outcome_var)) %>% 
      as.matrix() 
    
    z0_mat <- df %>% 
      dplyr::filter(!!sym(treat_var) != treat_unit) %>% 
      dplyr::filter(!!sym(treat_time_var) < treat_time) %>% 
      dplyr::filter(!!sym(treat_var) %in% tag_obs_x0_names) %>% 
      dplyr::select(!!sym(treat_var), !!sym(treat_time_var), 
                    !!sym(outcome_var)) %>% 
      tidyr::pivot_wider(names_from = !!sym(treat_time_var), 
                         values_from = !!sym(outcome_var)) %>% 
      dplyr::select(-!!sym(treat_var)) %>% 
      as.matrix() %>% 
      t() 
    
    return(list(x1 = x1_mat, x0 = x0_mat, x0_obs = tag_obs_x0, 
                x0_units = tag_obs_x0_names,
                z1_mat = z1_mat, z0_mat = z0_mat))
    
  } else {
    x1 <- df %>% 
      dplyr::filter(!!sym(treat_var) == treat_unit)
    x0 <- df %>% 
      dplyr::filter(!!sym(treat_var) != treat_unit)
    
    x0_n <- nrow(x0 %>% dplyr::distinct(!!sym(treat_var)))
    
    # Declare empty matrices 
    x1_mat <- matrix(nrow = 1, ncol = 0)
    x0_mat <- matrix(nrow = x0_n, ncol = 0)
    
    # Bind unit names 
    x0_names <- x0 %>% dplyr::distinct(!!sym(treat_var)) %>% dplyr::pull() 
    x0_mat <- x0_mat %>% cbind(x0_names)
    
    # First, handle variables not specified in the lag 
    vars_pre <- covars[!covars %in% names(covars_lag)]
    
    if (length(vars_pre) > 0) {
      x1_pre <- x1 %>% 
        dplyr::filter(!!sym(treat_var) == treat_unit) %>% 
        dplyr::summarise(across(all_of(vars_pre), ~ mean(.x, na.rm = TRUE))) %>% 
        as.matrix()
      
      x0_pre <- x0 %>% 
        dplyr::filter(!!sym(treat_var) != treat_unit) %>% 
        dplyr::group_by(!!sym(treat_var)) %>% 
        dplyr::summarise(across(all_of(vars_pre), ~ mean(.x, na.rm = TRUE))) %>% 
        dplyr::select(-!!sym(treat_var)) %>% 
        as.matrix()
      
      
      x1_mat <- x1_mat %>% cbind(x1_pre)
      x0_mat <- x0_mat %>% cbind(x0_pre)
    }
    
    
    for (var in names(covars_lag)) {
      for (lag in covars_lag[[var]]) {
        # Handle treated covariates 
        var_lag <- df %>%
          dplyr::select(!!sym(treat_time_var), !!sym(var), !!sym(treat_var)) %>% 
          dplyr::filter(!!sym(treat_time_var) == treat_time - lag)
        
        x1_mat <- x1_mat %>% 
          cbind(var_lag %>%
                  dplyr::filter(!!sym(treat_var) == treat_unit) %>% 
                  dplyr::pull(!!sym(var))
          )
        
        x0_mat <- x0_mat %>% 
          cbind(var_lag %>%
                  dplyr::filter(!!sym(treat_var) != treat_unit) %>% 
                  dplyr::pull(!!sym(var))
          )
      }
    }
    
    names <- covars_lag %>% 
      purrr::imap(~ paste0(.y, "_l_", .x)) %>% 
      unlist() %>% 
      unname()
    
    if (length(vars_pre) > 0) {
      names <- c(vars_pre, names)
    }
    
    tag_obs_x0 <- which(rowSums(is.na(x0_mat))==0) # tag nonmissing obs.
    x0_mat_omit <- x0_mat %>% na.omit() 
    tag_obs_x0_names <- unique(x0_mat_omit[, 1])
    x0_mat <- matrix(as.numeric(x0_mat[, -1] %>% na.omit()), 
                     ncol = ncol(x0_mat[, -1]))
    
    colnames(x0_mat) <- names
    colnames(x1_mat) <- names 
    
    if (nrow(na.omit(x1_mat)) == 0) {
      stop("Error! Treated unit is missing some values to match the pre-period")
    }
    
    # Transpose to get the correct dimensions 
    x1_mat <- x1_mat %>% t() 
    x0_mat <- x0_mat %>% t()
    
    # Grab outcomes 
    # z1: length of pre-period x 1
    # z0: length of pre-period x 
    z1_mat <- df %>% 
      dplyr::filter(!!sym(treat_var) == treat_unit) %>% 
      dplyr::filter(!!sym(treat_time_var) < treat_time) %>% 
      dplyr::select(!!sym(outcome_var)) %>% 
      as.matrix() 
    
    z0_mat <- df %>% 
      dplyr::filter(!!sym(treat_var) != treat_unit) %>% 
      dplyr::filter(!!sym(treat_time_var) < treat_time) %>% 
      dplyr::filter(!!sym(treat_var) %in% tag_obs_x0_names) %>% 
      dplyr::select(!!sym(treat_var), !!sym(treat_time_var), 
                    !!sym(outcome_var)) %>% 
      tidyr::pivot_wider(names_from = !!sym(treat_time_var), 
                         values_from = !!sym(outcome_var)) %>% 
      dplyr::select(-!!sym(treat_var)) %>% 
      as.matrix() %>% 
      t() 
    
    # x0: number of predictors x number of control donors
    # x1: number of predictors x 1
    out <- list(
      x1 = x1_mat, 
      x0 = x0_mat, 
      x0_obs = tag_obs_x0,
      x0_units = tag_obs_x0_names,
      z0_mat = z0_mat,
      z1_mat = z1_mat
    )
  }
  
}

scRun <- function(mat_list, params, V_method = NULL) {
  # Check correct weighting matrix methods 
  v_allowed <- c("identity", "mahalanobis")
  if (!is.null(V_method)) {
    if (!V_method %in% v_allowed) {
      stop("Error! You have specified a weighting option that is not allowed.") 
    }
  }
  
  x0 <- mat_list$x0
  x1 <- mat_list$x1 
  z0 <- mat_list$z0
  z1 <- mat_list$z1
  
  # Normalize x0 and x1 
  x0_n <- ncol(x0)
  temp <- cbind(x0, x1)
  # Divide by the standard deviation => everything has unit variance 
  scaled <- t(scale(t(temp), center = FALSE, 
                    scale = apply(temp, 1, sd, na.rm = TRUE)))
  x0_norm <- scaled[, 1:x0_n] %>% as.matrix()
  x1_norm <- scaled[, ncol(scaled)] %>% as.matrix()
  
  # Set up the weight matrix V which is square and has dimensions 
  # no. of predictors x no. of predictors
  if (is.null(V_method)) {
    k <- nrow(x0_norm)
    
    # First, try a diagonal matrix with equal weights on the covariates 
    v_diag <- rep(1/k, k)
    vstar_diag <- optim(par = v_diag, fn = scLoss, x0_norm = x0_norm,
                        x1_norm = x1_norm, z0 = z0, z1 = z1, 
                        params = params, method = "BFGS")$par
    
    # Set up a regression of Z onto X
    X <- cbind(x1_norm, x0_norm)
    X <- cbind(rep(1, ncol(X)),t(X))
    Z <- cbind(z1, z0)
    beta <- try(solve(crossprod(X, X)) %*% t(X) %*% t(Z), silent=TRUE)
    
    if(inherits(beta,"try-error")) {
      V <- vstar_diag
    } else {
      beta <- beta[-1, ] # drop the constants 
      # Grab the diagonals of the outer product 
      vstar_beta <- diag(beta %*% t(beta))
      V <- vstar_beta/sum(vstar_beta) 
    }
  } else if (V_method == "identity") {
    k <- nrow(x0_norm)
    V <- rep(1/k, k)
  } else if (V_method == "mahalanobis") {
    X <- cbind(x1, x0)
    V <- diag(v_calc(t(X)))
  }
  
  
  # Now, optimize over W
  Vstar <- diag(x = as.numeric(abs(V)/sum(abs(V))),
                nrow = length(V), ncol = length(V))
  
  
  mod.w <- list() 
  
  mod.w$A <- matrix(rep(1, times = ncol(x0_norm)), nrow = 1) # weights sum to 1 
  mod.w$Q <- 0.5 * t(x0_norm) %*% Vstar %*% x0_norm
  mod.w$obj <- as.vector(-t(t(x0_norm) %*% Vstar %*% x1_norm))
  mod.w$rhs <- c(1)
  mod.w$sense <- c("=")
  
  # Return solution as a vector along with weights 
  w <- round(gurobi::gurobi(mod.w, params = params)$x, 3)
}

scLoss <- function(v, x0_norm, x1_norm, z0, z1, params) {
  # Normalize V
  V <- diag(x = as.numeric(abs(v)/sum(abs(v))),
            nrow = length(v), ncol = length(v))
  
  # Find the optimal W by solving a quadratic program
  mod.w <- list()
  
  mod.w$A <- matrix(rep(1, times = ncol(x0_norm)), nrow = 1) # weights sum to 1
  mod.w$Q <- 0.5 * t(x0_norm) %*% V %*% x0_norm
  mod.w$obj <- as.vector(-t(t(x0_norm) %*% V %*% x1_norm))
  mod.w$rhs <- c(1)
  mod.w$sense <- c("=")
  
  w <- round(gurobi::gurobi(mod.w, params = params)$x, 3)
  
  
  # Calculate MSE of the SC estimator 
  # z1 is pre-period x 1, z0 is pre-period x N0, W is N0 x 1
  loss <- crossprod(z1 - z0 %*% w, 
                    z1 - z0 %*% w) %>% as.numeric()
  loss <- loss/nrow(z0)
  return(loss)
}

sc_att_calc <- function(df, treat_var, treat_unit, outcome_var, W,
                        treat_time_var, dataprep_df, treat_time) {
  
  # Grab outcomes for treated unit 
  # Not including year and treatment identifier, this will be a T x 1 dataframe 
  # where T = number of periods in the data 
  z1_full <- df %>% 
    dplyr::filter(!!sym(treat_var) == treat_unit) %>% 
    dplyr::select(!!sym(treat_time_var), !!sym(outcome_var), !!sym(treat_var))
  
  z1_names <- z1_full %>% dplyr::select(!!sym(treat_var)) %>% dplyr::pull() 
  z1_times <- z1_full %>% 
    dplyr::select(!!sym(treat_time_var)) %>% 
    dplyr::pull() 
  
  z1_full <- z1_full %>% 
    dplyr::select(-!!sym(treat_time_var), -!!sym(treat_var))
  
  # Grab outcomes for control units 
  z0_full <- df %>% 
    dplyr::filter(!!sym(treat_var) %in% dataprep_df$x0_units) %>% 
    dplyr::select(!!sym(treat_time_var), !!sym(outcome_var), !!sym(treat_var))
  
  z0_full <- z0_full %>% 
    tidyr::pivot_wider(names_from = !!sym(treat_time_var), 
                       values_from = !!sym(outcome_var)) %>% 
    dplyr::select(-!!sym(treat_var))
  
  # Convert to matrix and calculate ATTs 
  z1_full <- as.matrix(z1_full)
  z0_full <- as.matrix(z0_full) %>% t()
  W <- as.matrix(W, nrow = length(w))
  
  att_mat <- z1_full - z0_full %*% W
  
  att_df <- cbind(z1_names, z1_times, as.data.frame(att_mat))
  
  colnames(att_df) <- c(treat_var, treat_time_var, "att")
  
  att_df <- att_df %>% 
    dplyr::mutate(event_time = !!sym(treat_time_var) - treat_time)
  
  return(att_df)
}

scRoutine <- function(df, treat_var, treat_unit, treat_time_var, treat_time,
                      covars, outcome_var, covars_lag, params,
                      outcome_var_att, V_method = NULL) {
  
  df_in <- data_prep_sc(df, treat_var = treat_var, treat_unit = treat_unit, 
                        treat_time_var = treat_time_var, 
                        treat_time = treat_time,
                        covars = covars, outcome_var = outcome_var, 
                        covars_lag = covars_lag)
  
  w <- scRun(df_in, params = params, V_method = V_method)
  names <- df_in$x0_units
  att <- sc_att_calc(df, treat_var, treat_unit, outcome_var = outcome_var_att,
                     W = w, treat_time_var, dataprep_df = df_in,
                     treat_time = treat_time)
  
  return(list(w = w, names = names, att = att))
}

# Marginal Treatment Effect routines ----
# Function to calculate the MTE with binary treatment 
mte_calc <- function(df, y_var, treat_var) {
  # Grab the data with D = 1 and D = 0
  df_treat <- df %>% 
    filter(!!sym(treat_var) == 1)
  
  df_untreat <- df %>% 
    filter(!!sym(treat_var) == 0)
  
  # Separate regressions 
  reg_d1 <- ols_run(df_treat, y = y_var, covars = "phat")
  reg_d0 <- ols_run(df_untreat, y = y_var, covars = "phat")
  
  # Calculate the MTE 
  b0 <- reg_d0$estimates[2]*2
  a0 <- reg_d0$estimates[1] - b0/2
  a1 <- reg_d1$estimates[1]
  b1 <- reg_d1$estimates[2]*2
  
  mte_vec <- (a1 - a0) + (b1 - b0) * u_grid  
  
  mte_out <- list(
    mte_vec = mte_vec,
    y_var = y_var,
    treat_var = treat_var,
    call = match.call()
  )
}

mte_bootstrap <- function(mte_results, u_grid, reps = 100) {
  # Setup for bootstrap: declare matrix with rows equal to the size of the 
  # U grid. Each column will contain the results from one iteration.
  mte_bstrap_mat <- matrix(nrow = length(u_grid), ncol = 1,
                           data = u_grid)  
  
  # Grab the data used in estimation 
  # This assumes a df has all complete observations for now 
  call <- mte_results[["call"]]
  df <- rlang::eval_tidy(call[["df"]])
  treat_var <- mte_results[["treat_var"]]
  y_var <- mte_results[["y_var"]]
  
  # Handles proportions in creating the bootstrap sample 
  levels <- df %>% 
    dplyr::select(!!sym(treat_var)) %>% 
    dplyr::distinct() %>% 
    dplyr::pull() 
  
  var_values <- df %>% 
    dplyr::select(!!sym(treat_var)) %>% 
    dplyr::pull() 
  
  obs_labels <- match(var_values, levels)
  probs <- levels %>% 
    purrr::map_dbl(~ sum(var_values == .x))
  
  probs <- probs/length(var_values)
  
  obs_probs <- obs_labels %>% 
    purrr::map_dbl(~ probs[.x])
  
  idx <- 1:nrow(df) # indices 
  
  for (i in 1:reps) {
    # Create bootstrap treated and untreated observations 
    idx_sample <- sample(idx, size = nrow(df), replace = TRUE, prob = obs_probs)
    df_boot <- df[idx_sample, ]
    
    # Calculate the MTE over the vector of U's
    mte_boot <- mte_calc(df = df_boot, y_var = y_var, 
                         treat_var = treat_var)$mte_vec
    
    mte_bstrap_mat <- mte_bstrap_mat %>% cbind(mte_boot)
  }
  
  # Compute the bootstrap standard deviation if requested.
  mte_bstrap_se <- mte_bstrap_mat[, -1] %>% 
    apply(., 1, sd)
  
  return(mte_bstrap_se)
}

# Function to calculate the ATT weights 
att_mte_wt <- function(df, treat_var, u_vec) {
  pr_treat <- nrow(df %>% filter(!!sym(treat_var) == 1))/nrow(df)
  
  p_vec <- eval(rlang::parse_expr(paste0("df$", "phat")))
  
  num <- u_vec %>% 
    purrr::map_dbl(~ sum(p_vec >= .x)/length(p_vec))
  
  wt <- num/pr_treat
}

# Function to calculate the ATU weights 
atu_mte_wt <- function(df, treat_var, u_vec) {
  pr_untreat <- nrow(df %>% filter(!!sym(treat_var) == 0))/nrow(df)
  
  p_vec <- eval(rlang::parse_expr(paste0("df$", "phat")))
  
  num <- u_vec %>% 
    purrr::map_dbl(~ sum(p_vec < .x)/length(p_vec))
  
  wt <- num/pr_untreat
}

# Function to bootstrap target parameters from MTE 
target_param_bstrap <- function(mte_results, reps, target, u_grid) {
  
  # Setup for bootstrap: declare matrix with columns equal to the number of 
  # iterations 
  target_bstrap_mat <- matrix(nrow = 1, ncol = reps)
  
  # Grab the data used in estimation 
  # This assumes a df has all complete observations for now 
  call <- mte_results[["call"]]
  df <- rlang::eval_tidy(call[["df"]])
  treat_var <- mte_results[["treat_var"]]
  y_var <- mte_results[["y_var"]]
  
  # Handles proportions in creating the bootstrap sample 
  levels <- df %>% 
    dplyr::select(!!sym(treat_var)) %>% 
    dplyr::distinct() %>% 
    dplyr::pull() 
  
  var_values <- df %>% 
    dplyr::select(!!sym(treat_var)) %>% 
    dplyr::pull() 
  
  obs_labels <- match(var_values, levels)
  probs <- levels %>% 
    purrr::map_dbl(~ sum(var_values == .x))
  
  probs <- probs/length(var_values)
  
  obs_probs <- obs_labels %>% 
    purrr::map_dbl(~ probs[.x])
  
  idx <- 1:nrow(df) # indices 
  
  for (i in 1:reps) {
    # Create bootstrap treated and untreated observations 
    idx_sample <- sample(idx, size = nrow(df), replace = TRUE, prob = obs_probs)
    df_boot <- df[idx_sample, ]
    
    # Calculate the MTE over the vector of U's
    mte_boot <- mte_calc(df = df_boot, y_var = y_var, 
                         treat_var = treat_var)$mte_vec
    
    # Gather results
    mte_df <- tibble(
      u = u_grid,
      mte = mte_boot 
    )
    
    # Compute the target parameter 
    # ATT 
    if (target == "att") {
      mte_df <- mte_df %>% 
        mutate(att_wt = att_mte_wt(df_boot, treat_var = "jikokoa", u = u)) 
      
      att <- sum(mte_df$mte * mte_df$att_wt)/length(mte_df$u)
      
      target_bstrap_mat[, i] <- att
      
    } else if (target == "atu") {
      mte_df <- mte_df %>% 
        mutate(atu_wt = atu_mte_wt(df_boot, treat_var = "jikokoa", u = u))
      
      atu <- sum(mte_df$mte * mte_df$atu_wt)/length(mte_df$u)
      
      target_bstrap_mat[, i] <- atu
    } else if (target == "ate") {
      ate <- sum(mte_df$mte)/length(mte_df$u)
      
      target_bstrap_mat[, i] <- ate
    }
  }
  
  se <- sd(as.vector(target_bstrap_mat))
  
  return(se)
}

################################################################################
# DID MONTE CARLO FUNCTIONS 
################################################################################
# Helper functions ----
# Function to calculate all the weights on the static TWFE estimator, based on 
# the residualized treatment indicator.
weight_calc <- function(df_resid, ref_cohort) {
  props <- df_resid %>% 
    dplyr::group_by(cohort) %>% 
    dplyr::summarise(
      prop = n()/nrow(df_resid),
      .groups = "keep"
    )
  
  df_sum <- df_resid %>% 
    dplyr::group_by(cohort, year) %>% 
    dplyr::summarise(
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
    dplyr::select(cohort, starts_with("rel_time"), starts_with("unit_"),
           starts_with("year_")) %>% 
    dplyr::mutate(tag_obs = if_else(cohort == e & rel_time == r, 1, 0)) %>% 
    dplyr::mutate(tag_obs = replace_na(tag_obs, 0))
  
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
    dplyr::distinct(cohort, rel_time) %>% 
    dplyr::arrange(cohort, rel_time)
  
  wts <- map2(
    cohorts_times$cohort,
    cohorts_times$rel_time,
    ~ as.data.frame(get_wts_e_r(df_did, e = .x, r = .y, g = g,
                                rel_times_exclude = rel_times_exclude))
  ) %>% 
    purrr::list_rbind() %>% 
    cbind(cohorts_times$cohort) %>% 
    cbind(cohorts_times$rel_time)
  
  names(wts) <- c("wt", "cohort", "rel_time")
  
  return(wts)
}

# Calculates the implied coefficient on a certain relative time in the TWFE
# dynamic specification
coef_calculate <- function(df_wts) {
  df_wts <- df_wts %>% 
    dplyr::mutate(att_contrib = att * wt) %>% 
    dplyr::mutate(att_homog_contrib = att_homog * wt) %>% 
    dplyr::summarise(b_att = sum(att_contrib, na.rm = TRUE),
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
    dplyr::mutate(var = if_else(var == "b_att_mean", "Heterogeneous ATT's", 
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
    dplyr::group_by(cohort, rel_time) %>% 
    dplyr::summarise(wt = mean(wt), .groups = "keep") 
  
  df_out <- df_out %>% 
    left_join(trt_df, by = c("cohort", "rel_time")) %>% 
    dplyr::select(-year) %>% 
    left_join(trt_df_homog, by = c("cohort", "rel_time")) %>% 
    dplyr::select(-year)
}

dynamic_coef_df_process <- function(df, M) {
  df_out <- df %>% 
    dplyr::summarise(b_att_mean = mean(b_att),
              b_att_sd = sd(b_att),
              b_att_homog_mean = mean(b_att_homog),
              b_att_homog_sd = sd(b_att_homog),
              .groups = "keep") %>% 
    pivot_longer(cols = contains("_mean"), names_to = "var",
                 values_to = "estim") %>% 
    dplyr::mutate(sd = if_else(var == "b_att_mean", b_att_sd, b_att_homog_sd)) %>% 
    dplyr::select(-contains("_sd"))
  
  # Confidence intervals
  df_out <- df_out %>% 
    dplyr::mutate(ci_95up = estim + qt(.975, M-1)*sd,
           ci_95lo = estim - qt(.975, M-1)*sd)
}

# Cohort-weighted ATT 
weighted_att_calc <- function(df, r) {
  df_att <- df %>% 
    dplyr::filter(rel_time == r) %>% 
    dplyr::group_by(cohort) %>% 
    dplyr::summarise(size = n(),
              att = att,
              .groups = "keep") %>% 
    dplyr::slice(1) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(pr = size/sum(size))
  
  df_att <- df_att %>% 
    dplyr::mutate(att = replace_na(att, 0)) %>% 
    dplyr::mutate(wt_att_contrib = att * pr)
  
  att <- sum(df_att$wt_att_contrib)
}

# CS approach 
# Function that calculates ATT(g, t)
cs_gt_calc <- function(df, e, r, covar_str = NULL) {
  if (!is.null(covar_str)) {
    df <- df %>% 
      dplyr::mutate(across(everything(), ~ as.character(.x))) %>% 
      dplyr::filter(eval(rlang::parse_expr(covar_str)))
  }
  # Setup for r >= 0 
  if (r >= 0) {
    df_chrt <- df %>% 
      dplyr::filter(cohort == e) %>% 
      dplyr::filter(year == e + r | year == e - 1) %>% 
      dplyr::mutate(pre = if_else(year == e - 1, 1, 0))
    
    tag_not_yet <- df %>% 
      dplyr::filter(year == e + r & d == 0) 
    
    chrts_not_yet <- unique(tag_not_yet$cohort)
    
    df_not_yet <- df %>% 
      dplyr::filter(cohort %in% chrts_not_yet) %>% 
      dplyr::filter(year == e + r | year == e - 1) %>% 
      dplyr::mutate(pre = if_else(year == e - 1, 1, 0))
    
    if (nrow(df_not_yet) == 0) {
      return(NA_real_)
    } 
  } else {
    df_chrt <- df %>% 
      dplyr::filter(cohort == e) %>% 
      dplyr::filter(year == e + r | year == e + r - 1) %>% 
      dplyr::mutate(pre = if_else(year == e +r - 1, 1, 0))
    
    tag_not_yet <- df %>% 
      dplyr::filter(year == e + r & d == 0) 
    
    chrts_not_yet <- unique(tag_not_yet$cohort)
    
    df_not_yet <- df %>% 
      dplyr::filter(cohort %in% chrts_not_yet) %>% 
      dplyr::filter(year == e + r | year == e + r - 1) %>% 
      dplyr::mutate(pre = if_else(year == e + r - 1, 1, 0))
    
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
cs_did <- function(df, r, covars = NULL) {
  df <- na.omit(df)
  
  if (!is.null(covars)) {
    # needs a different grouping structure 
    df_att <- df %>% 
      dplyr::filter(rel_time == r) %>% 
      dplyr::mutate(across(all_of(covars), ~ as.character(.x))) %>% 
      dplyr::group_by(across(all_of(c("cohort", covars)))) %>% 
      dplyr::summarise(size = n(), .groups = "keep") 
    
  } else {
  df_att <- df %>% 
    dplyr::filter(rel_time == r) %>% 
    dplyr::group_by(cohort) %>% 
    dplyr::summarise(size = n(), .groups = "keep") 
  } 
  
  chts <- unique(df_att$cohort)
  att_gt <- c()
  
  if (!is.null(covars)) {
    cells <- unique(df[covars])
    df_loop <- df %>% 
      dplyr::select(all_of(covars)) %>% 
      dplyr::distinct(across(all_of(covars))) %>% 
      dplyr::mutate(across(all_of(covars), 
                           ~ paste0(dplyr::cur_column(), "==", 
                                    paste0("\"", .x, "\""))))
    
    covar_strs <- df_loop %>% 
      tidyr::unite(col = "var", sep = " & ") %>% 
      dplyr::pull() 
    
    for (cht in chts) {
      print(paste0("Cohort: ", cht))
      for (str in covar_strs) {
        att_gt <- append(att_gt, 
                         cs_gt_calc(df, e = cht, r = r, covar_str = str))
      }
    }
    
    df_att_gt <- tibble(
      cohort = rep(chts, each = length(covar_strs)),
      att_gt = att_gt,
      covar_str = rep(covar_strs, each = length(chts))
    ) %>% 
      tidyr::separate(covar_str, sep = "&", into = covars)
    
    #names(df_att_gt) <- c("cohort", "att_gt", covars)
    df_att_gt %>% glimpse()
    
    # Left bind to the main df and calculate final estimate, dropping groups 
    # with no control comparisons
    df_att <- df_att %>% 
      left_join(df_att_gt, by = c("cohort", covars)) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(pr = size/sum(size)) %>% 
      dplyr::mutate(att_contrib = pr * att_gt)
    
    att_e <- sum(df_att$att_contrib, na.rm = TRUE)
    
  } else {
  # Run the main CS estimator
  att_gt <- chts %>% 
    purrr::map(~ cs_gt_calc(df, e = .x, r = r))
  
  df_att_gt <- tibble(
    cohort = chts,
    att_gt = as.numeric(att_gt)
  )
  
  # Left bind to the main df and calculate final estimate, dropping groups 
  # with no control comparisons
  df_att <- df_att %>% 
    left_join(df_att_gt, by = "cohort") %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(pr = size/sum(size)) %>% 
    dplyr::mutate(att_contrib = pr * att_gt)
  
  att_e <- sum(df_att$att_contrib, na.rm = TRUE)
  }
}

# Imputation approach 
impute_did <- function(df, r, rel_times_exclude = -1, covars = NULL) {
  
  # Select untreated subsample and regress on cohort + time FEs
  df_d0 <- df %>% 
    dplyr::filter(d == 0)
  
  cohorts_include <- paste0("cohort_", unique(df_d0$cohort)[-1])
  years_include <- paste0("year_", unique(df_d0$year)[-1])
  
  covars_include <- c(cohorts_include, years_include)
  
  if (!is.null(covars)) {
    covars <- covars[!grepl("year_", covars) | 
                       grepl(paste(years_include, collapse = "|"), covars)]
    covars_include <- c(covars_include, covars)
  }
  
  mod_d0 <- ols_run(df = df_d0, y = "y", 
                    covars = covars_include)
  
  estims_d0 <- mod_d0[["estimates"]]
  
  mat_predict <- matrixPrep(df, cols = covars_include, constant = TRUE)
  
  # Get fitted values and subtract them from each actual outcome to construct
  # Y dot
  y0_fit <- mat_predict %*% estims_d0
  df <- df %>% 
    cbind(y0_fit = y0_fit) %>% 
    dplyr::mutate(y_dot = y - y0_fit)
  
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