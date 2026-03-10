#' Chenyi Shirley He
#' Likelihood-Based Inference & Monte Carlo Validation


# =============================================================================
# PART 1: Exponential Regression via MLE
# =============================================================================

#' Negative log-likelihood for exponential regression
#'
#' @param par numeric vector of length 2: (beta0, beta1)
#' @param y   numeric vector of strictly positive responses
#' @param x   numeric vector of covariates
#' @return numeric scalar: negative log-likelihood
negloglik_expreg <- function(par, y, x) {
  beta0  <- par[1]
  beta1  <- par[2]
  lambda <- exp(beta0 + beta1 * x)
  sum(-log(lambda) + lambda * y)
}


#' Fit exponential regression by MLE
#'
#' @param y numeric vector of strictly positive responses
#' @param x numeric vector of covariates
#' @return list: par_hat, hessian, se
fit_expreg <- function(y, x) {
  opt <- optim(
    par     = c(0, 0),
    fn      = negloglik_expreg,
    y       = y,
    x       = x,
    hessian = TRUE
  )
  H_mat   <- opt$hessian
  cov_mat <- solve(H_mat)
  se      <- sqrt(diag(cov_mat))
  list(par_hat = opt$par, hessian = H_mat, se = se)
}


#' Wald confidence intervals
#'
#' @param y     numeric vector of responses
#' @param x     numeric vector of covariates
#' @param alpha significance level (default 0.05)
#' @return data.frame with rows beta0, beta1 and columns hat, lo, hi, se
CI <- function(y, x, alpha = 0.05) {
  fit     <- fit_expreg(y, x)
  par_hat <- fit$par_hat
  se      <- fit$se
  z       <- qnorm(1 - alpha / 2)
  data.frame(
    hat = par_hat,
    lo  = par_hat - z * se,
    hi  = par_hat + z * se,
    se  = se,
    row.names = c("beta0", "beta1")
  )
}


#' beta0 as closed-form function of beta1 (dimension reduction)
#'
#' @param beta1 numeric scalar
#' @param y     numeric vector of responses
#' @param x     numeric vector of covariates
#' @return numeric scalar
beta0_hat <- function(beta1, y, x) {
  log(length(y)) - log(sum(y * exp(beta1 * x)))
}


#' 1D profile log-likelihood as a function of beta1
#'
#' @param beta1 numeric scalar or vector
#' @param y     numeric vector of responses
#' @param x     numeric vector of covariates
#' @return numeric: log-likelihood value(s)
loglik_1d <- function(beta1, y, x) {
  sapply(beta1, function(b1) {
    b0     <- beta0_hat(b1, y, x)
    lambda <- exp(b0 + b1 * x)
    sum(log(lambda) - lambda * y)
  })
}


#' Monte Carlo coverage simulation for 95% Wald CI on beta1
#'
#' @param n     sample size
#' @param beta0 true beta0
#' @param beta1 true beta1
#' @param x     observed covariates (resampled from)
#' @param N     number of replications (default 10000)
#' @param alpha significance level (default 0.05)
#' @param seed  optional random seed
#' @return numeric: empirical coverage probability
estimate_coverage <- function(n, beta0, beta1, x, N = 10000,
                              alpha = 0.05, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  z        <- qnorm(1 - alpha / 2)
  coverage <- logical(N)
  for (i in seq_len(N)) {
    x_sam      <- sample(x, n, replace = TRUE)
    y_sam      <- rexp(n, rate = exp(beta0 + beta1 * x_sam))
    fit        <- fit_expreg(y_sam, x_sam)
    b1_hat     <- fit$par_hat[2]
    se_b1      <- fit$se[2]
    coverage[i] <- (beta1 >= b1_hat - z * se_b1) &&
                   (beta1 <= b1_hat + z * se_b1)
  }
  mean(coverage)
}


# =============================================================================
# PART 2: Latent Variable Poisson Model
# =============================================================================

#' Crude Monte Carlo estimator for the latent-variable Poisson PMF
#'
#' Model: X ~ N(0, sigma^2), Y | X=x ~ Poisson(exp(mu + x))
#'
#' @param m     integer vector of count values
#' @param mu    log-mean parameter
#' @param sigma standard deviation of latent variable
#' @param K     number of Monte Carlo samples
#' @return data.frame: m, pmf, se
dlatentpois_mc <- function(m, mu, sigma, K) {
  x      <- rnorm(K, mean = 0, sd = sigma)
  result <- sapply(m, function(mi) {
    g <- exp(-exp(mu + x) + mi * (mu + x)) / factorial(mi)
    c(pmf = mean(g), se = sd(g) / sqrt(K))
  })
  result <- matrix(result, nrow = 2,
                   dimnames = list(c("pmf", "se"), NULL))
  data.frame(m = m, pmf = result["pmf", ], se = result["se", ])
}


#' Importance sampling estimator for the latent-variable Poisson PMF
#'
#' Uses proposal N(a, sigma^2) with recommended a = log(m) - mu
#'
#' @param m     integer vector of count values
#' @param mu    log-mean parameter
#' @param sigma standard deviation of latent variable
#' @param a     proposal mean shift
#' @param K     number of Monte Carlo samples
#' @return data.frame: m, pmf, se
dlatentpois_is <- function(m, mu, sigma, a, K) {
  x      <- rnorm(K, mean = a, sd = sigma)
  w      <- dnorm(x, mean = 0, sd = sigma) / dnorm(x, mean = a, sd = sigma)
  result <- sapply(m, function(mi) {
    g <- w * exp(-exp(mu + x) + mi * (mu + x)) / factorial(mi)
    c(pmf = mean(g), se = sd(g) / sqrt(K))
  })
  result <- matrix(result, nrow = 2,
                   dimnames = list(c("pmf", "se"), NULL))
  data.frame(m = m, pmf = result["pmf", ], se = result["se", ])
}
