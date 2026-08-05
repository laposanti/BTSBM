# Architecture proposal: a unified package for BT and Plackett–Luce models

## Purpose

`BTSBM` currently fits a Bayesian Bradley--Terry (BT) model and a BT stochastic
block model (BT--SBM) from a square matrix of pairwise win counts.  The adjacent
`PLuce` project provides working Plackett--Luce (PL) inference for the three
PL latent structures that should be brought into this package.  The delivery
scope is therefore:

| Likelihood | Baseline | Item partition | Ranking/sample mixture | Latent block model |
|---|---|---|---|---|
| Bradley--Terry | retain | retain BT--SBM | **deferred** | **deferred** |
| Plackett--Luce | add | add PL--SBM | add PL mixture | add PL--LBM |

Here PL--LBM means the existing PLuce *joint clustering* model: it clusters
ranking rows and items simultaneously.  It is the concrete PL co-clustering
implementation to port.  BT mixture and BT co-clustering remain documented
future models, but are explicitly outside the present implementation plan.

The package should be organised around these two independent decisions:

1. **Observation likelihood**: a pairwise outcome (BT) or an ordered/partial
   ranking (PL).
2. **Latent structure**: none, a mixture over observations, a partition of
   items, or a partition of both sources and items.

This avoids building six unrelated samplers.  Each model becomes a particular
combination of a likelihood module, a structural-prior module, and an inference
kernel.  The existing BT--SBM becomes the first implementation of that
combination, not a special architecture.

## Initial implementation status

The first implementation increment has started. It includes `as_bt_data()`,
`as_rankings()` (with `as_pl_data()` retained as a compatibility alias), shared model/control constructors, the S3 `btsbm_fit` contract,
and pure-R Gibbs samplers for simple PL and PL--SBM.  Tests cover input
validation, fit dimensions, ranking-level log likelihoods, the current BT
wrappers, and PL mixture/PL--LBM reference samplers. The mixture and LBM
samplers use the PLuce exponential-race architecture in readable R; their C++
augmentation/exposure kernel is registered and parity-tested against the R
reference. Further C++ optimisation remains gated on equivalent parity tests.
BT mixture and BT co-clustering remain deferred.

## Package identity

The name `BTSBM` accurately describes the current package but will no longer
describe base PL, mixtures, or co-clustering.  Because the current version is
`0.1.0`, decide before the first broad release whether to retain the name for
continuity or rename the package to a model-neutral name (for example,
`rankblocks`).  If `BTSBM` is retained, its Title and README should say
“Bayesian models for pairwise comparisons and rankings with latent structure”,
and the class/interface should not hard-code `bt_sbm`.  The examples below use
`btsbm_*` only as a conservative migration path; a renamed package should use a
neutral prefix such as `rank_*` instead.

## PLuce is the reference implementation

The source to port is `/Users/lapo_santi/Desktop/Nial/PLuce/PL-SBM`, with the
following division of responsibilities:

| PLuce component | Existing capability | Package action |
|---|---|---|
| `R/SBM_PL.R` | Gnedin item clustering, exponential-race augmentation | Port as the PL--SBM kernel |
| `R/plsbm_sample_mixture.R` | finite and Gnedin clustering of ranking rows; split--merge moves | Port as the PL-mixture kernel |
| `R/plsbm_joint_clustering.R` | joint ranking-row/item clustering, `C × K` ability matrix | Port as the PL--LBM kernel |
| `src/*.cpp` | augmentation, exposure aggregation, collapsed allocation scores | Compile through this package's normal Rcpp registration |
| `R/plsbm_interface.R` | `fit_PL()` facade, partition summaries, simulations | Re-express as package-native constructors, methods, and tests |

`Reproducibility Support PL-LBM` is a valuable frozen comparison target, but
not a second source tree to import.  Select `PL-SBM/R` and `PL-SBM/src` as the
code source of truth, then use the reproducibility version's simulations and
saved results for parity tests.  In particular, do **not** copy PLuce's runtime
`source()` / `Rcpp::sourceCpp()` loader or its use of `.GlobalEnv`; package
loading and generated `RcppExports` should own compilation and registration.

## Recommended model meanings

The names “mixture” and “co-clustering” are not sufficient by themselves; they
need unambiguous generative definitions in the API and documentation.

### Common ability parameterisation

Use positive abilities \(\lambda_i\), or log-abilities
\(\eta_i=\log\lambda_i\), internally.  The public summaries should normally
report item-level \(\eta_i\) (and optionally \(\lambda_i\)) because they are
available for every model.  Apply one identifiability constraint per relevant
ability vector, for example `mean(eta) = 0`.

For BT, an observation between items \(i\) and \(j\) has

\[
  P(i \succ j) = \frac{\lambda_i}{\lambda_i + \lambda_j}.
\]

For a strict partial PL ranking \(\pi=(\pi_1,\ldots,\pi_m)\) from an offered
set \(A\), use

\[
  P(\pi\mid A,\lambda) =
  \prod_{q=1}^{m}
  \frac{\lambda_{\pi_q}}
       {\sum_{j\in A\setminus\{\pi_1,\ldots,\pi_{q-1}\}}\lambda_j}.
\]

For the initial port, preserve PLuce's well-tested data assumption: every row
shares the same universe of \(n\) available items and `rho` may record either a
complete ranking or its top \(m\) prefix.  The denominator is therefore the
full item total minus previously selected items.  Add varying offered sets and
ties only after that path has a regression-tested package implementation.

### The PL model lattice to port

Let \(\rho_{\ell,\cdot}\) be ranking row \(\ell\), let \(z_\ell\) be its
sample/ranking-cluster label, and let \(x_i\) be an item-block label.  The four
PL models are a simple nested lattice:

| Model | Latent labels | Ability used for item \(i\) in row \(\ell\) |
|---|---|---|
| `pl_model()` | none | \(\lambda_i\) |
| `pl_sbm_model()` | \(x_i\) | \(\lambda_{x_i}\) |
| `pl_mixture_model()` | \(z_\ell\) | \(\lambda_{z_\ell,i}\) |
| `pl_lbm_model()` | \(z_\ell, x_i\) | \(\lambda_{z_\ell,x_i}\) |

`pl_sbm_model()` is PLuce's `mode = "item_clustering"`.
`pl_mixture_model()` is `mode = "sample_clustering"`: it clusters whole
ranking rows, not ranking positions or individual choices.  For a finite
mixture it has a fixed number \(C\) of row clusters and Dirichlet weights.  The
Gnedin version infers occupied \(C\), optionally accelerated with split--merge
moves.  The model is a mixture of PL ability vectors, one for each ranking-row
cluster.

`pl_lbm_model()` is PLuce's `mode = "joint_clustering"`, and is the preferred
public name for the PL co-clustering model.  It has a \(C\times K\) ability
matrix: row-cluster \(c\) ranks item-block \(k\) with weight
\(\lambda_{c k}\).  Both partitions can use Gnedin priors; the item partition
may alternatively be fixed at \(K\).  Call a row a `sample` in low-level
objects, but expose `ranking_id` in the user API.  It may represent a rater or
context only when there is one ranking per rater/context.

The exponential-race augmentation used in PLuce gives the sufficient statistics
that make this compositional design practical.  For each ranking position draw
\(Z_{\ell q}\) from an exponential risk-set time, then aggregate wins and
exposures as follows:

| Model | Gamma-update statistics |
|---|---|
| PL | \(w_i, S_i\) |
| PL--SBM | \(w_k, S_k\) |
| PL mixture | \(w_{c i}, S_{c i}\) |
| PL--LBM | \(w_{c k}, S_{c k}\) |

Thus PL--LBM is not a separate likelihood.  It is the same PL augmentation
with two allocation maps and a different aggregation target.

### Identifiability and partition priors

PL likelihoods identify relative row weights only.  Preserve the PLuce
`norm_method` concept in the model/control specification, with
`"logmean0"` as the default: normalise each ability row to geometric mean one
and retain the auxiliary scale only when the exact augmented transition needs
it.  The identified representation can also be recorded as
\(\pi_{c k}=\lambda_{c k}/\sum_k\lambda_{c k}\) and a nuisance scale
\(\tau_c\), as in PLuce's `compute_pi_tau()` implementation.  Never compare
or relabel unnormalised PL ability rows by their arbitrary scale.

The BT--SBM and PL partitioned models share DP, Pitman--Yor,
Dirichlet--multinomial, and Gnedin priors.  They should be reusable
`partition_prior` objects, not several hyperparameters repeated in every
public fitting function.  The first PL port should cover the Gnedin and
finite-prior cases already exercised by PLuce; the other existing BT priors can
be added to PL only after the shared interface is stable.

### Deferred BT mixture and BT co-clustering

Do not implement `bt_mixture_model()` or `bt_coclustering_model()` in this
cycle.  They require a decision about what is clustered: a match, an aggregated
pair count, a rater/session, or another context.  In particular, a mixture per
matrix cell is not generally equivalent to a mixture per contest.  The current
pairwise matrix has no row/source entity, so it cannot support the PL--LBM-style
co-clustering model without a new long-form BT data contract.  Reserve the
names and keep the model interface extensible, but do not expose incomplete
constructors or imply that they are implemented.

## Public interface

### User-facing vocabulary

The PL interface uses names that express statistical roles rather than source
code symbols: `rankings` replaces `rho`, `latent_strength` replaces `lambda`,
and `gnedin_hyperparameter` replaces `gamma_GN`. The underlying samplers may
still use compact mathematical variables, but those are not part of the public
contract. A model therefore reads like its generative story:

```r
pl_sbm_model(
  latent_strength = latent_strength(shape = 1, rate = 1),
  item_clustering = gnedin_prior(gnedin_hyperparameter = 0.5)
)
```

All implemented PL clustering models accept `gnedin_prior()`,
`dirichlet_process_prior()`, `pitman_yor_prior()`, or
`finite_partition(max_clusters = ...)`. `clustering_prior()` exposes the same
choice through one explicit `family` argument. Legacy constructors remain as
compatibility aliases during migration.

### Data constructors

Accept familiar matrices for simple pairwise data, but turn all inputs into
typed data objects before fitting.

```r
# Existing matrix workflow remains valid
pairwise <- as_bt_data(wins = w_ij)

# PLuce-parity input: rows are ranked item ids; top-m is allowed
ranking_data <- as_rankings(rankings, item_count = n)
```

`as_bt_data()` should validate zero/self comparisons, non-negative counts,
unique item ids, and whether aggregation is allowed.  For the first PL release,
`as_rankings(rankings, item_count)` validates a dense integer `L × m` ranking matrix:
every entry is in `1:n_items`, no item is repeated in a row, and every row has
the same full availability universe.  Its `ranking_id` is initially the row
index, but is stored explicitly so later input forms can retain external ids.

An `as_pl_data()` long-form method with `ranking_id`, `item_id`, `rank`, and an
explicit availability table is the correct next input extension.  It should
not delay the PLuce-parity port, and it should not be silently emulated by
treating absent ranked items as unavailable.

Internally, `btsbm_pairwise_data` can contain an edge table and an optional
fast sparse representation.  `btsbm_ranking_data` should initially contain the
validated integer `rho`, `n_items`, `n_rankings`, `rank_length`, and ranking
ids; later it can add choice-set offsets.  Inference code consumes these
prepared forms, not raw matrices or data frames.

### Model and control constructors

Use a single fit entry point plus small explicit constructors:

```r
model <- bt_sbm_model(
  latent_strength = latent_strength(shape = 4),
  item_clustering = gnedin_prior(gnedin_hyperparameter = 0.8)
)

fit <- fit_btsbm(
  pairwise, model,
  control = mcmc_control(iter = 2_000, warmup = 1_000, chains = 4)
)
```

Equivalent convenience constructors are useful discoverability wrappers:

```r
bt_model()                  # simple BT
pl_model()                  # simple PL
bt_sbm_model(item_clustering = gnedin_prior(gnedin_hyperparameter = .8))
pl_sbm_model(item_clustering = gnedin_prior(gnedin_hyperparameter = .8))
pl_mixture_model(ranking_clustering = gnedin_prior(gnedin_hyperparameter = .6))
pl_lbm_model(
  ranking_clustering = gnedin_prior(gnedin_hyperparameter = .6),
  item_clustering = gnedin_prior(gnedin_hyperparameter = .6)
)
```

For a fixed finite PL mixture, the equivalent constructor is
`pl_mixture_model(ranking_clustering = finite_partition(max_clusters = C), ...)`. Do not give
every fitting function a long list of `alpha_PY`, `sigma_PY`, `beta_DM`,
`K_DM`, and `gamma_GN`.  A prior constructor validates only the parameters it
needs, makes defaults visible, and can be used in every implemented partition.

`mcmc_control()` should hold iteration counts, seed, stored state, and
algorithmic moves.  This folds the many PLuce tuning arguments into a clear
structure, for example `partition_moves = partition_moves(split_merge = 2,
pair_proposal = "informed", restricted_scans = 1)` and separate `item_updates`
and `ranking_updates` settings for PL--LBM.

### A stable fit contract

Every call to `fit_btsbm()` should return an S3 `btsbm_fit` object with a
common top-level shape:

```r
list(
  call       = match.call(),
  data       = <validated btsbm_*_data>,
  model      = <btsbm_model>,
  control    = <mcmc_control>,
  draws      = <raw sampler draws>,
  diagnostics = <chain and sampler diagnostics>,
  elapsed    = <timing information>
)
```

The top-level contract is stable, while `draws` follows the model's actual
parameter dimensions.  This mirrors PLuce's useful distinction between raw
allocations, raw/identified `lambda`, and implied item weights.

| Field | Present for | Meaning |
|---|---|---|
| `item_cluster` | BT--SBM, PL--SBM, PL--LBM | `S × n` raw item labels |
| `ranking_cluster` | PL mixture, PL--LBM | `S × L` raw ranking-row labels |
| `n_item_clusters` | item-partition models | length-`S` occupied-`K` trace |
| `n_ranking_clusters` | PL mixture, PL--LBM | length-`S` occupied-`C` trace |
| `ability` | every model | model-native positive ability draws: `S × n`, ragged `K`, `C × n`, or ragged `C × K` |
| `ability_scale` | PL models when needed | auxiliary row scales used by the exact augmented state |
| `ability_simplex` | PL mixture/LBM | identified row-normalised weights `pi` (and optional `tau`) |
| `augmented` | optional | latent exponential-race or BT augmentation draws |

For a PL mixture or LBM, do not manufacture an ambiguous universal `S × n`
`item_ability` field.  Instead provide an explicit accessor such as
`implied_ability(fit, ranking = ..., summary = ...)`, which returns the
ability vector conditional on a selected ranking row/cluster or an explicitly
declared averaging rule.  The accessor may produce the rectangular quantities
needed for a plot, while `log_lik()` and prediction use the model-native
`ability` with its allocations.  This prevents a component-averaged ability
from being mistaken for a likelihood parameter.

Implement S3 methods:

```r
print(fit); summary(fit); plot(fit, type = "ability")
log_lik(fit, unit = "observation")
posterior_membership(fit, target = "item")
posterior_similarity(fit, target = "item")
predict(fit, newdata, type = "probability")
```

`log_lik()` must return `S × R`, where columns are conditionally independent
units: observed pair/count cells or individual contests for BT, and whole
ranking rows for every PL mode.  It should carry a `unit_index` attribute.
This replaces model-specific LOO builders and makes `loo_btsbm(fit)` a thin,
safe wrapper around `loo::loo()`.

## Internal layout

Suggested files (the exact file names are less important than the boundaries):

```text
R/
  data-bt.R                 # as_bt_data(), validation, matrix-to-edge conversion
  data-pl.R                 # as_pl_data(), ranking and choice-set preparation
  model.R                   # btsbm_model and public model constructors
  prior-ability.R           # gamma/log-normal ability priors and constraints
  prior-partition.R         # DP, PY, DM, Gnedin; predictive masses
  fit.R                     # fit_btsbm(), btsbm_fit, print/summary methods
  kernel-bt.R               # BT likelihood, augmentation, updates
  kernel-pl.R               # PL likelihood and exponential-race augmentation
  structure-none.R          # no allocation state
  structure-item-sbm.R      # BT--SBM and PL--SBM item updates
  structure-pl-mixture.R    # ranking-row allocation and finite/Gnedin priors
  structure-pl-lbm.R        # joint ranking-row/item allocation updates
  posterior.R               # relabelling, summaries, membership and similarity
  predictive.R              # log_lik(), predict(), posterior_predict()
  loo.R                     # loo_btsbm(), compare_btsbm()
  plotting.R                # plots dispatching from btsbm_fit
  compat.R                  # deprecated wrappers for the 0.1 API
src/
  bt_augmentation.cpp       # optional, narrow C++ numerical kernels
  pl_augmentation.cpp       # draw exponential race times and exposures
  pl_allocation.cpp          # collapsed allocation candidate scores
```

The fitting loop should depend on a small kernel protocol rather than an
`if`/`switch` branch at every update:

```r
kernel <- make_kernel(data, model)
state  <- initialize_state(data, model, control)

for (iter in seq_len(control$iter)) {
  state <- kernel$augment(state, data, model)
  state <- kernel$update_ability(state, data, model)
  state <- kernel$update_structure(state, data, model)
  state <- kernel$identify(state, model)
  recorder$save(state, iter)
}
```

For the current positive-ability BT sampler, Gamma/exponential augmentation
provides conjugate ability updates.  PLuce confirms that PL has the analogous
exponential-race augmentation, so PL--SBM, PL mixture, and PL--LBM can share a
single risk-set/exposure module.  The mixture and LBM differences belong in
their allocation and aggregation modules, not in three copied PL likelihood
implementations.

C++ should remain an implementation detail behind narrowly named functions
such as `bt_draw_exposure()`, `pl_draw_race_times()`,
`pl_aggregate_exposures()`, and `pl_candidate_log_marginals()`.  It must not
receive a package-wide `fit` list or decide model semantics.  This keeps the R
model layer testable and lets PLuce's optimised kernels enter without importing
its script-oriented runtime design.

## Posterior summaries and label switching

There are three different exchangeabilities:

1. item-block labels in an SBM;
2. ranking-row cluster labels in a PL mixture; and
3. ranking-row and item-block labels in PL--LBM.

The current `relabel_by_lambda()` is appropriate only for a one-dimensional
ordered ability-block model.  Replace it internally with a generic
`canonicalize_draws()` dispatcher:

* **SBM**: compute partition summaries from raw allocations, then optionally
  order item blocks by a stated scalar (for example, posterior block ability).
* **PL mixture**: relabel ranking clusters by an allocation-based method such
  as ECR or Stephens' method, then report component-specific predictive
  summaries.
* **PL--LBM**: resolve ranking-row and item permutations jointly; report both
  posterior similarity matrices and cell-level (\(C\times K\)) summaries.

Canonical labels are presentation aids, not data required by likelihood or LOO
calculations.  Plotting and summary methods should use the fit object and call
the appropriate canonicalisation method; users should not pass separate
`x_samples` and `lambda_samples` by hand.

## Concrete refactor: the current BT--SBM function

Today, `gibbs_bt_sbm()` validates a matrix, chooses an urn with `switch()`,
manages labels and latent variables, runs the Markov chain, and returns an
unnamed list.  It is the right prototype for BT--SBM, but extending it with PL
risk sets, mixture allocations, and two partitions would make it unmaintainable.

The desired public use is:

```r
data <- as_bt_data(wins = w_ij)
model <- bt_sbm_model(
  ability_prior = gamma_ability(shape = 4),
  item_prior = partition_prior("gnedin", gamma = 0.8)
)
fit <- fit_btsbm(
  data, model,
  control = mcmc_control(iter = 500, warmup = 250, seed = 1, store = "minimal")
)

summary(fit)
loo_btsbm(fit)
plot(fit, type = "membership")
```

The compatibility function should be small and preserve the existing workflow:

```r
#' @export
gibbs_bt_sbm <- function(w_ij, a = 4, prior = c("DP", "PY", "DM", "GN"),
                         alpha_PY = NA_real_, sigma_PY = NA_real_,
                         beta_DM = NA_real_, K_DM = NA_integer_,
                         gamma_GN = NA_real_, T_iter = 2000, T_burn = 1000,
                         init_x = NULL, store_z = FALSE, verbose = TRUE) {
  lifecycle::deprecate_soft("0.2.0", "gibbs_bt_sbm()", "fit_btsbm()")

  prior <- match.arg(prior)
  fit <- fit_btsbm(
    data = as_bt_data(wins = w_ij),
    model = bt_sbm_model(
      ability_prior = gamma_ability(shape = a),
      item_prior = partition_prior_from_legacy(
        prior, alpha_PY, sigma_PY, beta_DM, K_DM, gamma_GN
      )
    ),
    control = mcmc_control(
      iter = T_iter, warmup = T_burn, init_item_cluster = init_x,
      store = if (store_z) "augmented" else "minimal", progress = verbose
    )
  )

  as_legacy_bt_sbm_fit(fit)
}
```

The new `fit_btsbm()` owns validation, storage, diagnostics, and the shared
fit contract.  `bt_sbm_model()` owns model selection.  A BT kernel owns only
the augmentation and conditional updates.  `as_legacy_bt_sbm_fit()` returns
the old `x_samples`, `lambda_samples`, `K_per_iter`, and `z_samples` names
during a documented transition period.

The use of `lifecycle` above is illustrative; it can be introduced as an
import, or the first release can issue `warning()` without adding a dependency.

## Concrete refactor: the PLuce interface

PLuce's `fit_PL(rho, mode = ..., ...)` is the right transitional facade, but
its current signature mixes model choice, priors, sampler moves, compilation,
and output normalisation.  Preserve the three modes, while moving their meaning
into data, model, and control objects:

| Current PLuce mode | New constructor | Main raw draws |
|---|---|---|
| `"item_clustering"` | `pl_sbm_model()` | `item_cluster`, `ability[K]` |
| `"sample_clustering"` | `pl_mixture_model()` | `ranking_cluster`, `ability[C × n]` |
| `"joint_clustering"` | `pl_lbm_model()` | `ranking_cluster`, `item_cluster`, `ability[C × K]` |

For example, a PL--LBM fit becomes:

```r
rankings <- as_pl_data(rho, n_items = n_items)
model <- pl_lbm_model(
  ability_prior = gamma_ability(shape = 1),
  ranking_prior = partition_prior("gnedin", gamma = 0.5),
  item_prior = partition_prior("gnedin", gamma = 0.5),
  identifiability = pl_identifiability("logmean0")
)
fit <- fit_btsbm(
  rankings, model,
  control = mcmc_control(
    iter = 4_000, warmup = 1_000, seed = 42,
    partition_moves = partition_moves(
      ranking_split_merge = 1, item_split_merge = 1,
      pair_proposal = "informed", restricted_scans = 1
    )
  )
)

posterior_membership(fit, target = "ranking")
posterior_membership(fit, target = "item")
```

If a package-level `fit_PL()` transition helper is useful, it should only map
the old mode to a model constructor and call the common fitter:

```r
fit_PL <- function(rho, mode = c("item_clustering", "sample_clustering", "joint_clustering"),
                   n_iter = 4000, burn = 1000, C = NULL, K_fixed = NULL,
                   a = 1, b = 1, gamma_gn = .5, gamma_gn_z = .5,
                   n_items = NULL, ...) {
  mode <- match.arg(mode)
  model <- switch(
    mode,
    item_clustering = pl_sbm_model(
      ability_prior = gamma_ability(shape = a, rate = b),
      item_prior = partition_prior("gnedin", gamma = gamma_gn)
    ),
    sample_clustering = pl_mixture_model(
      ability_prior = gamma_ability(shape = a, rate = b),
      ranking_prior = if (is.null(C)) partition_prior("gnedin", gamma = gamma_gn) else finite_partition(C)
    ),
    joint_clustering = pl_lbm_model(
      ability_prior = gamma_ability(shape = a, rate = b),
      ranking_prior = partition_prior("gnedin", gamma = gamma_gn_z),
      item_prior = if (is.null(K_fixed)) partition_prior("gnedin", gamma = gamma_gn) else finite_partition(K_fixed)
    )
  )
  fit_btsbm(as_pl_data(rho, n_items = n_items), model,
             mcmc_control(iter = n_iter, warmup = burn, ...))
}
```

This deliberately does not reproduce the PLuce runtime loader.  Its C++
kernels become registered package code, `...` is validated by
`mcmc_control()`, and each returned fit has the same `btsbm_fit` contract as
the existing BT models.

## Current code issues to resolve before expanding

These are structural concerns visible in the current package, not reasons to
change the statistical model:

* `make_bt_simple_loo()`, `make_bt_cluster_loo()`, and
  `compare_bt_models_loo()` are defined in both `R/utils.R` and
  `R/loo_helper.R`.  There must be one implementation, tested through the
  generic `log_lik()` path.
* `gibbs_bt_sbm()` allocates `L_cap_trace` and documents
  `L_cap_per_iter`, but does not fill or return it.  The implementation and
  documented fit contract need to agree before it becomes a compatibility API.
* Fit results are raw lists, while downstream functions sometimes expect raw
  labels and sometimes relabelled arrays.  A classed fit object removes this
  implicit, fragile contract.
* `relabel_by_lambda()` directly invokes optional `salso` functionality and
  is specialised to one item partition.  Posterior summaries should test for
  optional packages when called and dispatch by model structure.
* The current simulation function sets `seed = 123` instead of calling
  `set.seed(seed)`.  Simulation should share the same input validation and
  reproducibility conventions as fitting.
* The only current tests are plotting tests.  New kernels require likelihood,
  input-validation, sampler-invariant, posterior-predictive, and LOO tests;
  simulation-based calibration should be used for the core samplers.
* PLuce's `rho` kernels assume one dense item universe for every ranking row,
  one-based item ids, and no repeated item in a row.  Enforce those conditions
  in `as_pl_data()` before calling C++; do not leave validation distributed
  across three samplers.
* PLuce currently contains both the active `PL-SBM` tree and a reproducibility
  copy.  Port exactly one implementation into `R/` and `src/`; use the other
  for fixed-seed and posterior-summary comparison tests so they cannot drift.

## Implementation sequence

1. Add regression tests around the existing BT and BT--SBM examples, then
   remove duplicate LOO implementations and fix the documented sampler-output
   mismatch.
2. Introduce `as_bt_data()`, `as_pl_data()`, `btsbm_fit`, `mcmc_control()`,
   and `log_lik()`.  Keep `gibbs_bt_simple()` and `gibbs_bt_sbm()` as wrappers.
3. Port and test the common PL exponential-race augmentation from PLuce using
   fixed `rho` fixtures: complete rankings, top-\(m\) rankings, invalid ids,
   duplicate ids, and seeded posterior-predictive checks.
4. Implement baseline PL and port PL--SBM (`SBM_PL.R`) onto the common fit
   object.  Compare its traces, posterior `K`, and implied item abilities with
   the PLuce reference implementation.
5. Port PL mixture (`plsbm_sample_mixture.R`) for both finite and Gnedin row
   partitions.  Add ranking-row PSM/minVI summaries, finite-mixture weights,
   split--merge diagnostics, and ranking-level `log_lik()`.
6. Port PL--LBM (`plsbm_joint_clustering.R`) last: joint \(C,K\) traces,
   separate item/ranking PSMs, `C × K` identified ability rows, and independent
   split--merge diagnostics for each partition.  Test against PLuce's joint
   simulation/reproducibility artifacts before expanding the model.
7. Publish PL--SBM, PL-mixture, and PL--LBM vignettes plus a model-comparison
   vignette.  Mark the old function-specific BT LOO and relabelling helpers as
   superseded only when the generic methods cover their use cases.
8. **Defer BT mixture and BT co-clustering.**  Revisit them only after a
   long-form BT observation/source data design is agreed; they are not a
   dependency of the PL work above.

This order reuses PLuce's mature algorithms without importing its script
architecture, produces useful PL releases incrementally, and keeps the
unresolved BT extensions from delaying them.
