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
  print(meat)
  X <- ols_fit$design_mat 
  bread <- solve(t(X) %*% X)
  print(m)
  
  # Final result: only output the SE, making Stata finite sample correction 
  vcov <- (m/(m-1)) * ((n-1)/(n-k)) * (bread %*% meat %*% bread)
  se_out <- sqrt(diag(vcov))
}

# Bootstrap 
se_boot <- function(fun_results, seed = 12345, B = 50) {
  # Extract the call, data used, and number of obs. from the model object
  call <- fun_results[["call"]]
  df <- rlang::eval_tidy(call[["df"]])
  n <- fun_results[["n_obs"]]
  
  # Bootstrap the estimates 
  set.seed(seed)
  
  estims_bstrap <- 1:B %>% purrr::map(function(x) {
    df_boot <- dplyr::slice_sample(df, n = n, replace = FALSE)
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
    
    print(y0_use)
    
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