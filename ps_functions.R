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

# Matching helpers ---- 
# Calculates inverse variance weighting matrix for input matrix X
v_calc <- function(X) {
  n_vec <- apply(X, 2, function(x) {sum(!is.na(x))}) # nonmissing obs.
  colmeans <- colMeans(X, na.rm = TRUE)
  # R recycles, and doesn't broadcast like Python
  colmeans <- rep(colmeans, rep(nrow(X), ncol(X))) 
  X_out <- colSums((X - colmeans)^2)/(n_vec - 1) # k x 1 for each covariate
  
  V <- diag(1/X_out)
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
    intercept = constant,
    fitted_values = fitted_values,
    residuals = residuals,
    r_sq = r_sq,
    r_sq_adjust = r_sq_adjust
  )
}

# NN Matching for the ATT
find_treated_nn <- function(df, y, k, covars, treat_var) {
  # Setup: grab matrix of covariates, both in its entirety and separated by 
  # treated/untreated, and treated/untreated outcomes 
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
  # First, find the inverse variance weighting matrix. Then, calculate 
  # the distance matrix as an N1 x N0 matrix of matches 
  V <- v_calc(as.matrix(df_covars))
  
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
  
  # Calculate final estimate of the ATT, no. of obs, and output 
  att <- mean(y1 - y0_use)
  n <- length(y0) + length(y1)
  
  att_out <- list(
    estimates = att,
    call = match.call(),
    n_obs = n,
    tag_obs = tag_obs
  )
}