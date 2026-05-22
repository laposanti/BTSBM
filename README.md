# BTSBM

[![R-CMD-check](https://github.com/laposanti/BTSBM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/laposanti/BTSBM/actions/workflows/R-CMD-check.yaml)

Bayesian inference for Bradley-Terry stochastic block models for paired comparison data.

## Installation

```r
# After CRAN release
# install.packages("BTSBM")

# Development version
remotes::install_github("laposanti/BTSBM")
```

## Minimal example

```r
library(BTSBM)

w_ij <- ATP_2000_2025$`2017`$Y_ij

fit <- gibbs_bt_sbm(
  w_ij = w_ij,
  a = 4,
  prior = "GN",
  gamma_GN = 0.8,
  T_iter = 500,
  T_burn = 250,
  verbose = FALSE
)

post <- relabel_by_lambda(fit$x_samples, fit$lambda_samples)
plot_block_adjacency(fit = post, w_ij = w_ij)
```

## Documentation

- Package website: https://laposanti.github.io/BTSBM
- Function reference: https://laposanti.github.io/BTSBM/reference/index.html
- Getting started vignette: https://laposanti.github.io/BTSBM/articles/getting-started.html
