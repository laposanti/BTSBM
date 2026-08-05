# Implementation status and next work

## Working now

The package builds, installs, and passes its test suite. The common
high-level interface supports:

| Model | Status | Main entry point |
|----|----|----|
| Simple Bradley–Terry | working | [`bt_model()`](https://laposanti.github.io/BTSBM/reference/bt_model.md) + [`fit_btsbm()`](https://laposanti.github.io/BTSBM/reference/fit_btsbm.md) |
| Bradley–Terry item SBM | working | [`bt_sbm_model()`](https://laposanti.github.io/BTSBM/reference/bt_sbm_model.md) + [`fit_btsbm()`](https://laposanti.github.io/BTSBM/reference/fit_btsbm.md) |
| Simple Plackett–Luce | working | [`pl_model()`](https://laposanti.github.io/BTSBM/reference/pl_model.md) + [`fit_btsbm()`](https://laposanti.github.io/BTSBM/reference/fit_btsbm.md) |
| Plackett–Luce item SBM | working | [`pl_sbm_model()`](https://laposanti.github.io/BTSBM/reference/pl_sbm_model.md) + [`fit_btsbm()`](https://laposanti.github.io/BTSBM/reference/fit_btsbm.md) |
| Plackett–Luce ranking mixture | working R reference sampler | [`pl_mixture_model()`](https://laposanti.github.io/BTSBM/reference/pl_mixture_model.md) |
| Plackett–Luce latent block model | working R reference sampler | [`pl_lbm_model()`](https://laposanti.github.io/BTSBM/reference/pl_lbm_model.md) |
| Bradley–Terry mixture/co-clustering | deferred | – |

PL code now uses the public names `rankings`, `latent_strength`,
`item_clustering`, `ranking_clustering`, and `gnedin_hyperparameter`.
The compatibility names
[`as_pl_data()`](https://laposanti.github.io/BTSBM/reference/as_rankings.md),
[`gamma_ability()`](https://laposanti.github.io/BTSBM/reference/latent_strength.md),
and the earlier model arguments remain available, but new examples and
documentation use the descriptive vocabulary.

PL–SBM users can choose one of four clustering priors:

``` r

gnedin_prior(gnedin_hyperparameter = 0.5)
dirichlet_process_prior(concentration = 1)
pitman_yor_prior(concentration = 1, discount = 0.2)
finite_partition(max_clusters = 4, concentration = 1)
```

The current tests exercise each choice with a short PL–SBM fit.

## Required before the next PL release

1.  **Optimise PL mixture.** Bring the tested ranking-clustering kernel
    from `PLuce/PL-SBM/R/plsbm_sample_mixture.R` and its C++
    sufficient-statistic helpers into `R/` and `src/`. Expose
    ranking-cluster draws, mixture weights, ranking-level posterior
    similarity, and ranking-level
    [`log_lik()`](https://laposanti.github.io/BTSBM/reference/log_lik.md).
2.  **Optimise PL–LBM.** Port the joint row/item allocation kernel from
    `PLuce/PL-SBM/R/plsbm_joint_clustering.R`. Store both partitions,
    the `ranking_cluster × item_cluster` latent-strength draws, and
    separate diagnostics/moves for the two partitions.
3.  **Parity tests against PLuce.** Use fixed, small complete and
    top-`m` ranking fixtures to compare augmentation statistics,
    allocation-score calculations, seed-controlled traces, and posterior
    summaries. These tests should be in addition to the package-facing
    tests, not replaced by them.
4.  **Sampler diagnostics.** Add posterior-predictive ranking checks,
    multi-chain convergence helpers, effective sample size summaries,
    and mixture/LBM PSIS-LOO tests once the corresponding likelihood
    extractors exist.
5.  **Performance path.** The shared PL mixture/LBM augmentation and
    exposure kernel is now registered C++ with a fixed-seed parity test
    against its R reference. Move the remaining collapsed allocation
    sweeps and sufficient-statistic aggregation only after equivalent
    parity tests are green.

## Design work to schedule separately

- Extend ranking input to varying offered sets, ties, and long-form
  rankings; the present interface intentionally assumes one shared item
  universe and strict complete or fixed-length top-`m` rankings.
- Decide whether the package name should become model-neutral before a
  broad PL release.
- Resolve legacy vignette packaging: a source check currently has
  warnings because historic vignette sources have no generated
  `inst/doc` outputs. This does not affect installation, tests, or
  vignette rebuilding, but should be cleaned before release.
- Keep BT mixture and BT co-clustering deferred until there is an agreed
  source/observation data contract; they cannot be obtained by merely
  copying the PL row-mixture interface.
