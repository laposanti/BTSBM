---
title: 'BTSBM: a package for Bradley–Terry Stochastic Block Model'
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

## When this applies

Use this model whenever outcomes can be summarized as pairwise wins/losses or preferences:

* Sports & games: player/team `i` beats `j`.
* Information retrieval / ranking: document `i` preferred to `j` by a judge/user.
* A/B testing at scale: variant `i` preferred over `j` in head-to-head trials.
* Psychometrics & sensory studies: stimulus `i` chosen over `j`.
* Model selection by humans: method `i` judged better than `j`.

The model only requires that each comparison yields a binary outcome (“`i` over `j`”), possibly repeated and aggregated into counts.

---

## Required Inputs

Required inputs

To fit the model, provide aggregated pairwise comparison data as a square matrix:

* `w`: an (`n \times n`) matrix of wins/preferences, where `w[i, j]` is the number of times item (`i`) is preferred to item (`j`).

* Diagonal entries must be 0: `w[i, i] = 0` for all (`i`).

It follows that, for each unordered pair (`{i,j}`), the total number of comparisons is `n_{ij} ;=; w_{ij} + w_{ji}`, with `n_{ij} \in N` and `n{ii}=0`.
Throughout, we use the language of competition (e.g., “player (`i`) beats player (`j`)”), but the same structure applies to any pairwise comparison task (products, algorithms, stimuli, judges’ preferences, etc.).

What entries mean

* `w[i, j] = k` means item `i` beat/was preferred to `j` exactly `k` times.

* Binary data from single encounters are a special case: `w[i, j] ∈ {0,1}` and `w[i, j] + w[j, i] = 1` if they were compared once (or 0 if never compared).

Examples of `w`

(A) Aggregated counts (general case)

```r

# items: A, B, C, D
# w[i, j] = number of times i beat j
w <- matrix(
  c( 0, 3, 0, 2,
     1, 0, 4, 0,
     2, 0, 0, 1,
     0, 1, 3, 0 ),
  nrow = 4, byrow = TRUE,
  dimnames = list(c("A","B","C","D"), c("A","B","C","D"))
)
```

(B) Binary outcomes (single round)

```r

# one comparison per pair observed (0/1 wins)
w_bin <- matrix(
  c( 0, 1, 0, 1,
     0, 0, 1, 0,
     1, 0, 0, 0,
     0, 1, 1, 0 ),
  nrow = 4, byrow = TRUE,
  dimnames = list(c("A","B","C","D"), c("A","B","C","D"))
)
```

Attention: the w matrix must not be symmetric!
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
2. **Check the Gnedin prior hyperparameters** using `gnedin_K_mean` and `gnedin_K_var`
3. **Fit the model** with [`gibbs_bt_sbm()`](https://laposanti.github.io/BTSBM/reference/gibbs_bt_sbm.html).
4. **Relabel the output** with [`relabel_by_lambda()`](https://laposanti.github.io/BTSBM/reference/relabel_by_lambda.html).
5. **Visualize clustering structure** using [`plot_block_adjacency()`](https://laposanti.github.io/BTSBM/reference/plot_block_adjacency.html).
---

## Learn More

* 📘 [Function Reference](https://laposanti.github.io/BTSBM/reference/index.html): Complete list of functions and documentation.
* 📄 [Getting Started Vignette](https://laposanti.github.io/BTSBM/articles/getting-started.html): Conceptual background and reproducible examples.

---

## Citation

Santi, L., Friel, N. (2025). *The Bradley–Terry Stochastic Block Model.*
(Working paper, University College Dublin)

