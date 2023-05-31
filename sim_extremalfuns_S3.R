library(tidyverse)
library(fields)

# For callback/tracer
library(CSwR)

library(Rcpp)
sourceCpp("fBs_cov.cpp")

sim_extremalfuns <- function(model, x = NULL, quietly = FALSE, cb = NULL, cb_nsim = NULL) {
  # If model object contains $x, ignore argument x.
  if (!is.null(model$x)) x <- model$x
  if (is.null(x)) stop("Please provide points!")
  
  if (!is.null(dim(x)))
    N <- dim(x)[1]
  else
    N <- length(x)
  
  
  Gamma <- rexp(1)
  Y <- rPx(model, x, k = 1)
  nsim <- 1
  
  if (isFALSE(all.equal(Y[1], 1, check.names = FALSE)))
    warning("Y(x1) != 1. Make sure you're simulating under Px.")
  
  Z <- zY <- Gamma^(-1) * Y

  accept <- TRUE; k <- 1
  if(!is.null(cb)) cb()
  
  if (N == 1) return(Z)
  
  for (k in 2:N) {
    if(!quietly) cat("\r n:", k, "/", N)
    
    Y <- NA; zY <- NA; Gamma <- NA; accept <- NA
    if(!is.null(cb)) cb()
    
    Gamma <- rexp(1)
    while (Gamma^(-1) > Z[k]) {
      Y <- rPx(model, x, k)
      zY <- Gamma^(-1) * Y
      
      nsim <- nsim + 1
      
      accept <- FALSE
      if (all(zY[1:(k-1)] < Z[1:(k-1)])) {
        accept <- TRUE
        Z <- pmax(Z, zY)
      }
      
      if(!is.null(cb)) cb()
      Gamma <- Gamma + rexp(1)
    }
    Y <- NA; zY <- NA; accept <- NA
    if(!is.null(cb)) cb()
  }
  if(!quietly) cat("\n")
  if(!is.null(cb_nsim)) cb_nsim()
  Z
}

sim_sumnorm <- function(model, x = NULL, quietly = FALSE, cb = NULL, cb_nsim = NULL) {
  # If model object contains $x, ignore argument x.
  if (!is.null(model$x)) x <- model$x
  if (is.null(x)) stop("Please provide points!")
  
  if (!is.null(dim(x)))
    N <- dim(x)[1]
  else
    N <- length(x)
  
  Z <- numeric(N)
  Gamma <- rexp(1)
  
  nsim <- 0
  
  while (Gamma^(-1) * N > min(Z)) {
    n <- sample(N, 1)
    Y <- rPx(model, x, n)
    Z <- pmax(Z, Gamma^(-1) * N * Y / sum(Y))
    Gamma <- Gamma + rexp(1)
    
    nsim <- nsim + 1
  }
  if(!is.null(cb_nsim)) cb_nsim()
  Z
}

## rPx methods
## rPx() function calls rPx.<class> based on model object.
rPx <- function(model, x, k) {
  UseMethod("rPx", model)
}

rPx.BR <- function(model, x, k) {
  x_zero <- (x == 0)
  
  Z <- rnorm(ncol(model$L))
  
  W <- numeric(length(x))
  W[!x_zero] <- c(Z %*% model$L)
  
  exp(W - W[k] - abs(x - x[k])^model$alpha / 2)
}

rPx.BR_sheet <- function(model, x, k) {
  Z <- rnorm(ncol(model$L))
  W <- c(Z %*% model$L)

  ## Rcpp method to quickly calculate ||(x-x_n)/s||^α.
  sig <- rPx_BR_sheet_cov(model$x, k, model$alpha, model$s)

  exp(W - W[k] - sig / 2)
}

rPx.MM <- function(model, x, k) {
  u <- model$rh()
  model$h(x + u - x[k]) / model$h(u)
}
rPx.MM_sheet <- function(model, x, k) {
  u <- model$rh()
  y <- x + matrix(u-t(x)[,k], ncol = ncol(x), nrow = nrow(x), byrow = TRUE)
  model$h(y) / model$h(u)
}

# Max-stable model factories
BR_factory <- function(x, alpha = 1, s = 1) {
  x_mat <- as.matrix(x)
  L <- fBs_cov(x_mat, alpha, s) %>% chol()
  
  structure(list(x = x_mat, alpha = alpha, s = s, L = L), class = "BR_sheet")
}
BR_factory_univariate <- function(x, alpha = 1, s = 1) {
  x0 <- x[x != 0]
  fBm_cov <- function(t, s) 0.5 * (abs(t)^alpha + abs(s)^alpha - abs(t-s)^alpha)
  Sigma <- outer(x0,x0,fBm_cov)
  L <- chol(Sigma)
  
  structure(list(x = x, alpha = alpha, s = s, L = L), class = "BR")
}

smith_factory <- function(type = c("smith", "smith_sheet"), sigma = 1) {
  if (type == "smith") {
    h <- function(x) dnorm(x, sd = sigma)
    rh <- function() rnorm(1, sd = sigma)
    class <- c("MM", "smith")
  } else if (type == "smith_sheet") {
    h <- function(x) dnorm(x[,1], sd = sigma) * dnorm(x[,2], sd = sigma)
    rh <- function() matrix(rnorm(2, sd = sigma), nrow = 1) # k doesn't matter, we only call with k = 1
    class <- c("MM_sheet", "smith")
  } else {
    stop("type must be one of \"smith\" or \"smith_sheet\".")
  }

  structure(
    list(h = h, rh = rh),
    class = class
  )
}

interpolate_sheet <- function(df, length.out = 100) {
  require(fields)
  
  df <- df[order(df[[1]], df[[2]]),]
  
  obj <- list(
    x = unique(df[[1]]),
    y = unique(df[[2]]),
    z = matrix(
      df[[3]],
      nrow = length(unique(df[[1]]))
    )
  )
  loc <- list(
    x = seq(min(obj$x), max(obj$x), length.out = length.out),
    y = seq(min(obj$y), max(obj$y), length.out = length.out)
  )
  interpolated <- interp.surface.grid(obj, loc)
  
  out <- expand_grid(
    x = interpolated$x,
    y = interpolated$y
  ) %>% 
    mutate(Z = c(interpolated$z)) %>% 
    set_names(colnames(df)[1:3])
}

extcoeff_BR <- function(h, s, alpha) {
  2 * pnorm(sqrt((h / s)^alpha/4))
}