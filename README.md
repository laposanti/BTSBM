# BTSBM

[![R-CMD-check](https://github.com/laposanti/BTSBM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/laposanti/BTSBM/actions/workflows/R-CMD-check.yaml)

Bayesian inference for Bradley-Terry stochastic block models for paired comparison data.

The development version also includes an initial Plackett--Luce interface for
strict complete and top-`m` rankings: validated ranking data, simple PL, and
item-block PL--SBM fitting, ranking mixtures, and joint ranking/item latent
block models. The mixture and LBM currently use clear R reference kernels;
the PLuce C++ kernels remain a planned performance upgrade.

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

## Plackett--Luce ranking example

```r
library(BTSBM)

# Synthetic teaching data with Sushi Preference Dataset-style item labels.
# It contains no rows from the original, non-redistributable benchmark.
rankings <- sushi_toy(n_rankings = 40, rank_length = 5)

fit_pl <- fit_btsbm(
  rankings,
  pl_sbm_model(
    latent_strength = latent_strength(shape = 1, rate = 1),
    item_clustering = gnedin_prior(gnedin_hyperparameter = 0.5)
  ),
  mcmc_control(iter = 600, warmup = 300, seed = 1)
)

posterior_strength(fit_pl, summary = "mean")
posterior_similarity(fit_pl, target = "item")
```

`item_clustering` can also use `dirichlet_process_prior()`,
`pitman_yor_prior()`, or `finite_partition(max_clusters = ...)`. See
`?clustering_prior` for the common parameterisation.

See `vignette("pl-ranking-models", package = "BTSBM")` for the ranking-data
format, model scope, references, and the original Sushi benchmark's licence.
See `vignette("toy-data-gallery", package = "BTSBM")` for reproducible tennis,
Eurovision-inspired, Sushi-inspired, and joint-clustering examples.
See [IMPLEMENTATION_STATUS.md](IMPLEMENTATION_STATUS.md) for the current model
matrix and the remaining PLuce-porting work.

## Documentation

- Package website: https://laposanti.github.io/BTSBM
- Function reference: https://laposanti.github.io/BTSBM/reference/index.html
- Getting started vignette: https://laposanti.github.io/BTSBM/articles/getting-started.html
