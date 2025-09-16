---
title: "BTSBM: Bayesian Bradley–Terry SBM"
---

# Description

Implements Bayesian inference for Bradley–Terry models with a 
    Stochastic Block Model (SBM) prior on the players or items being compared. 
    Includes Gibbs-type MCMC samplers, posterior relabeling, model diagnostics, 
    and visualization utilities for inferred skill clusters and block 
    interaction strengths. Suitable for applications such as sports analytics, 
    psychometrics, and ranking problems with latent group structure.



## Installation

Install the development version from GitHub:

```r


# install.packages("devtools")
devtools::install_github("laposanti/BTSBM")

Example
library(BTSBM)

set.seed(123)
K <- 6
n <- matrix(0, K, K)
n[upper.tri(n)] <- sample(0:5, sum(upper.tri(n)), TRUE)
n <- n + t(n); diag(n) <- 0

w <- matrix(0, K, K)
w[upper.tri(w)] <- rbinom(sum(upper.tri(w)), n[upper.tri(n)], 0.5)
w <- w + t(n - w); diag(w) <- 0

fit <- gibbs_bt_sbm(w, n, a = 1, b = 1, prior = "GN",
                    n_iter = 500, burnin = 250, verbose = FALSE)

```


