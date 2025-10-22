---
title: 'BTSBM: Bradley–Terry Stochastic Block Model'
---

# Overview

**BTSBM** implements Bayesian inference for the *Bradley–Terry Stochastic Block Model* (BT–SBM), a model that combines **pairwise comparison data** with **block clustering** of items.

The package provides:

1) Gibbs-type MCMC samplers for inference on the posterior distribution;

2) Posterior relabeling and model diagnostics;

3) Visualization tools for cluster assignments and block interaction strengths.

Applications include **sports analytics**, **psychometrics**, and **ranking problems** with hidden group structure.

➡️ Jump directly to:


- The complete list of functions provided : [Function Reference](https://laposanti.github.io/BTSBM/reference/index.html)

- A vignette where a minimal workflow is provided: [Getting Started Vignette](https://laposanti.github.io/BTSBM/articles/getting-started.html)

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
````

Then load the package:

```r
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
fit <- gibbs_bt_sbm(w, a = 4, prior = "GN",
                    n_iter = 500, burnin = 250, verbose = FALSE)
```

---

## Workflow at a Glance

1. **Prepare input matrix** (`w`) from your pairwise data.
2. **Fit the model** with [`gibbs_bt_sbm()`](https://laposanti.github.io/BTSBM/reference/gibbs_bt_sbm.html).
3. **Relabel the output** with [`relabel_by_lambda()`](https://laposanti.github.io/BTSBM/reference/relabel_by_lambda.html).
3. **Visualize clustering structure** using [`plot_block_adjacency()`](https://laposanti.github.io/BTSBM/reference/plot_block_adjacency.html).
4. **Compare priors** (`DM`, `GN`, `PY`) via Leave-One-Out Information Criterion or other clustering metrics. [`make_bt_cluster_loo()`](https://laposanti.github.io/BTSBM/reference/make_bt_cluster_loo.html).

---

## Learn More

* 📘 [Function Reference](https://laposanti.github.io/BTSBM/reference/index.html): Complete list of functions and documentation.
* 📄 [Getting Started Vignette](https://laposanti.github.io/BTSBM/articles/btsbm.html): Conceptual background and reproducible examples.

---

## Citation

Santi, L., Friel, N. (2025). *The Bradley–Terry Stochastic Block Model.*
(Working paper, University College Dublin)

