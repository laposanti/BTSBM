---
title: "BTSBM: Bayesian Bradley–Terry Stochastic Block Model"
---

# Overview

**BTSBM** implements Bayesian inference for the *Bradley–Terry Stochastic Block Model* (BT–SBM),  
a model that combines **pairwise comparison data** with **block clustering** of items.

The package provides:
- Gibbs-type MCMC samplers for inference on the posterior distribution;
- Posterior relabeling and model diagnostics;
- Visualization tools for cluster assignments and block interaction strengths.

Applications include **sports analytics**, **psychometrics**, and **ranking problems** with hidden group structure.

➡️ Jump directly to:
- [Function Reference](reference/index.html)
- [Getting Started Vignette](articles/getting-started.html)

---

## Model Summary

The BT–SBM assumes that each item belongs to a latent cluster,  
and that the probability of one item defeating another depends on both  
their cluster-level *interaction strength* and their individual *skill parameter*.

Formally, for items \\( i, j \\):
\\[
\Pr(i \text{ beats } j) = \frac{\exp(\lambda_{x_i}}{\lambda_{x_i}+\lambda_{x_j}}
\\]
where:
- \\( \lambda_{x_i} \\) is the individual skill of item \\( i \\);
- \\( x_i \\in \{1, \dots, K\} \\) is its latent cluster label;

This implies that all items in the same cluster \\(i: x_i = k \\) share the same strength \\( \lambda_i = k \\)
---

## Required Inputs

To fit the model, you need **aggregated pairwise comparison data** in the form of a square matrix:

- `w`, the  matrix of *wins*, where `w[i, j]` is the number of times item `i` succeeds over item `j`; 
- Diagonal entries must be 0

---

## Installation

Install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("laposanti/BTSBM")
library(BTSBM)
```

---

## Example

Below is a minimal working example simulating toy data and fitting the BT–SBM with a Gnedin prior:

```r
set.seed(123)
K <- 6

# number of matches per block pair
n <- matrix(0, K, K)
n[upper.tri(n)] <- sample(0:5, sum(upper.tri(n)), TRUE)
n <- n + t(n); diag(n) <- 0

# number of wins per block pair
w <- matrix(0, K, K)
w[upper.tri(w)] <- rbinom(sum(upper.tri(w)), n[upper.tri(n)], 0.5)
w <- w + t(n - w); diag(w) <- 0

# fit the model
fit <- gibbs_bt_sbm(w, n, a = 1, b = 1, prior = "GN",
                    n_iter = 500, burnin = 250, verbose = FALSE)



