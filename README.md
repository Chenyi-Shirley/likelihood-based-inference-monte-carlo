# Likelihood-Based Inference & Monte Carlo Validation

**R | Maximum Likelihood Estimation | Monte Carlo Simulation | Importance Sampling**

---

## Overview

This project implements and evaluates two core statistical computing techniques:

1. **Exponential Regression via MLE** — Fitting a parametric survival/hazard model using maximum likelihood estimation, with analytical score functions and Wald confidence intervals.
2. **Latent Variable Poisson Model** — Approximating an intractable marginal PMF using crude Monte Carlo and importance sampling, with variance reduction analysis.

---

## Project Structure

```
.
├── report.html          # Rendered report (open in browser)
├── functions.R          # All custom R functions (documented)
├── README.md            # This file
└── data/
    └── expdata.csv      # Dataset used for exponential regression
```

---

## Part 1 — Exponential Regression via MLE

### Model

Response variable $Y_i \sim \text{Exponential}(\lambda_i)$ with log-linear mean:

$$\lambda_i = e^{\beta_0 + \beta_1 x_i}$$

### Key Results

| Parameter | MLE Estimate | Std. Error | 95% Wald CI |
|-----------|-------------|------------|-------------|
| $\beta_0$ | 0.6499 | 0.1191 | (0.4165, 0.8832) |
| $\beta_1$ | 0.5453 | 0.2084 | (0.1369, 0.9538) |
| $e^{\beta_1}$ (hazard ratio) | 1.7252 | 0.3595 | (1.1467, 2.5955) |
| $e^{-\beta_1}$ (waiting time multiplier) | 0.5796 | 0.1208 | (0.3853, 0.8720) |

### Coverage Simulation

Empirical coverage of the 95% Wald CI for $\beta_1$ ($N = 10{,}000$ replications):

| $n$ | Coverage | MCSE |
|-----|----------|------|
| 5 | 0.9257 | 0.0026 |
| 10 | 0.9413 | 0.0024 |
| 20 | 0.9490 | 0.0022 |
| 50 | 0.9495 | 0.0022 |
| 100 | 0.9514 | 0.0022 |
| 200 | 0.9486 | 0.0022 |
| 500 | 0.9494 | 0.0022 |

Coverage is near-nominal even at $n = 20$, consistent with asymptotic normality of the MLE.

---

## Part 2 — Latent Variable Poisson Model

### Model

$$X_i \sim N(0, \sigma^2), \quad Y_i \mid X_i = x \sim \text{Poisson}(e^{\mu + x})$$

The marginal PMF $p(m \mid \mu, \sigma)$ is intractable in closed form.

### Monte Carlo vs Importance Sampling (at $m = 10$)

| Method | Estimate | Std. Error |
|--------|----------|------------|
| Crude MC | 0.002880 | 0.000191 |
| Importance Sampling | 0.002845 | 0.000038 |

Importance sampling achieves a **~5× reduction in standard error**.

---

## Functions (`functions.R`)

| Function | Description |
|----------|-------------|
| `negloglik_expreg(par, y, x)` | Negative log-likelihood for exponential regression |
| `fit_expreg(y, x)` | MLE via `optim()`, returns estimates, Hessian, SE |
| `CI(y, x, alpha)` | Wald confidence intervals using delta method |
| `estimate_coverage(n, beta0, beta1, x, N, alpha)` | Monte Carlo coverage simulation |
| `dlatentpois_mc(m, mu, sigma, K)` | Crude MC estimator for latent Poisson PMF |
| `dlatentpois_is(m, mu, sigma, a, K)` | Importance sampling estimator |

---

## How to Run

```r
source("functions.R")

data <- read.csv("data/expdata.csv")
y <- data$y
x <- data$x

# Fit model
fit <- fit_expreg(y, x)

# Confidence intervals
CI(y, x, alpha = 0.05)

# Monte Carlo coverage
estimate_coverage(n = 100, beta0 = fit$par_hat[1], beta1 = fit$par_hat[2], x = x)

# Latent Poisson PMF
dlatentpois_mc(m = 0:15, mu = 0.2, sigma = 0.8, K = 10000)
dlatentpois_is(m = 10, mu = 0.2, sigma = 0.8, a = log(10) - 0.2, K = 5000)
```

---

## Dependencies

- **R** ≥ 4.0  
- Base R only: `stats`, `knitr`

---

## Author

**Chenyi Shirley He**
