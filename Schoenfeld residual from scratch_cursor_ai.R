# Compute Schoenfeld residuals from scratch (Cox PH, Breslow ties)

ensure_matrix <- function(X) {
  if (is.vector(X)) {
    X <- matrix(X, ncol = 1)
  }
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  X
}

safe_sort_index <- function(time, decreasing = FALSE) {
  # Stable order by time; ties keep input order
  order(seq_along(time), time, decreasing = decreasing)
}

cox_fit_nr <- function(time, status, X, max_iter = 50, tol = 1e-8) {
  # Fit Cox PH via Newton-Raphson with Breslow approximation for ties
  # Inputs:
  # - time: numeric vector of follow-up times
  # - status: 1 = event, 0 = censored
  # - X: covariate matrix (n x p)
  # Returns: list(beta, eta, loglik)

  n <- length(time)
  if (length(status) != n) stop("status length must equal time length")
  X <- ensure_matrix(X)
  if (nrow(X) != n) stop("nrow(X) must equal length(time)")
  p <- ncol(X)

  # Sort by increasing time; within ties keep order
  ord <- order(time, seq_along(time))
  time <- time[ord]
  status <- status[ord]
  X <- X[ord, , drop = FALSE]

  beta <- rep(0, p)
  loglik_prev <- -Inf

  for (iter in seq_len(max_iter)) {
    eta <- as.vector(X %*% beta)
    r <- exp(eta)

    # Compute cumulative risk set sums efficiently by traversing from last to first
    # S0[i] = sum_{k: time_k >= time_i} r_k
    # S1[i,] = sum X_k r_k over same risk set
    # S2 aggregate is used for Hessian (sum outer(X, X) r) over risk sets at event times
    S0 <- numeric(n)
    S1 <- matrix(0, nrow = n, ncol = p)

    cum_r <- 0
    cum_Xr <- rep(0, p)
    for (i in n:1) {
      cum_r <- cum_r + r[i]
      cum_Xr <- cum_Xr + r[i] * X[i, ]
      S0[i] <- cum_r
      S1[i, ] <- cum_Xr
    }

    # Group indices by unique event times for events only (Breslow)
    unique_times <- unique(time[status == 1])
    score <- rep(0, p)
    H <- matrix(0, nrow = p, ncol = p)
    loglik <- 0

    for (t0 in unique_times) {
      idx_event <- which(time == t0 & status == 1)
      # Risk set is everyone with time >= t0, i = first index with time == t0
      i0 <- min(which(time == t0))
      d <- length(idx_event)

      S0_t <- S0[i0]
      S1_t <- S1[i0, ]

      # For Hessian we also need S2 = sum r x x^T over risk set
      # Compute naively for clarity; acceptable for moderate n
      idx_risk <- which(time >= t0)
      X_risk <- X[idx_risk, , drop = FALSE]
      r_risk <- r[idx_risk]
      # S2 = sum r_i x_i x_i^T
      S2_t <- matrix(0, nrow = p, ncol = p)
      for (j in seq_along(r_risk)) {
        xj <- X_risk[j, ]
        S2_t <- S2_t + r_risk[j] * tcrossprod(xj)
      }

      # Contribution of this tied time (Breslow): d * (S1/S0) etc.
      x_bar <- S1_t / S0_t
      score <- score + colSums(X[idx_event, , drop = FALSE]) - d * x_bar

      H <- H - (S2_t / S0_t - tcrossprod(S1_t, S1_t) / (S0_t^2)) * d

      # Partial log-likelihood contribution
      loglik <- loglik + sum(eta[idx_event]) - d * log(S0_t)
    }

    # Newton step
    step <- tryCatch(solve(H, score), error = function(e) NULL)
    if (is.null(step)) {
      warning("Hessian is singular; stopping early.")
      break
    }

    beta <- beta - step

    if (abs(loglik - loglik_prev) < tol) {
      break
    }
    loglik_prev <- loglik
  }

  list(beta = beta, eta = as.vector(X %*% beta), loglik = loglik, order = ord)
}

schoenfeld_residuals <- function(time, status, X, beta = NULL, fit = NULL) {
  # Compute Schoenfeld residuals for each event time
  # If beta is NULL, fit Cox model first.

  n <- length(time)
  X <- ensure_matrix(X)
  if (nrow(X) != n) stop("nrow(X) must equal length(time)")

  if (is.null(beta)) {
    fit <- cox_fit_nr(time, status, X)
    beta <- fit$beta
  }

  # Sort by time to reuse risk set logic
  ord <- order(time, seq_along(time))
  time_s <- time[ord]
  status_s <- status[ord]
  X_s <- X[ord, , drop = FALSE]

  eta <- as.vector(X_s %*% beta)
  r <- exp(eta)

  n <- length(time_s)
  p <- ncol(X_s)
  # cumulative sums from the end for risk sets
  S0 <- numeric(n)
  S1 <- matrix(0, nrow = n, ncol = p)

  cum_r <- 0
  cum_Xr <- rep(0, p)
  for (i in n:1) {
    cum_r <- cum_r + r[i]
    cum_Xr <- cum_Xr + r[i] * X_s[i, ]
    S0[i] <- cum_r
    S1[i, ] <- cum_Xr
  }

  # For each event i, residual = x_i - E[X | risk set at time_i]
  # Using Breslow: expectation uses full risk set at that time
  idx_events <- which(status_s == 1)
  res <- matrix(NA_real_, nrow = length(idx_events), ncol = p)
  rownames(res) <- idx_events
  colnames(res) <- colnames(X)

  for (k in seq_along(idx_events)) {
    i <- idx_events[k]
    # Find first index with same time (risk set anchor)
    i0 <- min(which(time_s == time_s[i]))
    mu_t <- S1[i0, ] / S0[i0]
    res[k, ] <- X_s[i, ] - mu_t
  }

  # Return in original subject order for events, include mapping
  original_indices <- ord[idx_events]
  rownames(res) <- original_indices
  list(residuals = res, beta = beta, event_rows = original_indices)
}

#############################
# Minimal reproducible example
#############################
if (identical(environment(), globalenv())) {
  set.seed(123)
  n <- 100
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  beta_true <- c(log(1.5), -0.7)
  linpred <- beta_true[1] * x1 + beta_true[2] * x2
  # Exponential baseline hazard
  u <- runif(n)
  time <- -log(u) / exp(linpred)
  # Independent censoring
  cens <- rexp(n, rate = 0.1)
  status <- as.integer(time <= cens)
  time <- pmin(time, cens)

  X <- cbind(x1 = x1, x2 = x2)

  fit <- cox_fit_nr(time, status, X)
  cat("Estimated beta:\n")
  print(fit$beta)

  sch <- schoenfeld_residuals(time, status, X, beta = fit$beta)
  cat("Schoenfeld residuals (first 5 events):\n")
  print(head(sch$residuals, 5))
}