# Gibbs sampler for the Bradley–Terry Stochastic Block Model (BT–SBM)

Runs a Gibbs sampler for the BT–SBM with optional DP/PY/DM/GN priors on
the partition. Returns raw draws plus minimal summaries; `n_ij` is
computed internally as `w_ij + t(w_ij)`.

## Usage

``` r
gibbs_bt_sbm(
  w_ij,
  a = 4,
  prior = c("DP", "PY", "DM", "GN"),
  alpha_PY = NA_real_,
  sigma_PY = NA_real_,
  beta_DM = NA_real_,
  K_DM = NA_integer_,
  gamma_GN = NA_real_,
  T_iter = 2000,
  T_burn = 1000,
  init_x = NULL,
  store_z = FALSE,
  verbose = TRUE
)
```

## Arguments

- w_ij:

  Integer or numeric square matrix \\n \times n\\ of directed wins (i
  over j). Must be nonnegative with zero diagonal. The function builds
  \\n\_{ij} = w\_{ij} + w\_{ji}\\ internally.

- a:

  Positive shape parameter for the Gamma prior \\\lambda_k \sim
  \mathrm{Gamma}(a,b)\\. The algorithm uses \\b = \exp(\psi(a))\\ so
  that \\\mathbb{E}\[\log \lambda_k\] = 0\\ a priori.

- prior:

  Character scalar, one of `"DP"`, `"PY"`, `"DM"`, `"GN"`.

- alpha_PY, sigma_PY:

  Hyperparameters for Pitman–Yor / Dirichlet Process. For `prior="DP"`
  use `alpha_PY` (with `sigma_PY` ignored). For `prior="PY"` use both
  `alpha_PY` and `sigma_PY in (0,1)`.

- beta_DM, K_DM:

  Hyperparameters for the finite Dirichlet–Multinomial prior. `K_DM` is
  the maximum number of allowed clusters.

- gamma_GN:

  Hyperparameter for the Gnedin prior.

- T_iter, T_burn:

  Integers: total iterations and burn-in. Must satisfy
  `T_burn < T_iter`.

- init_x:

  Optional integer vector of length `n` with initial labels (1-based).

- store_z:

  Logical; if `TRUE`, store latent `Z` draws (heavy).

- verbose:

  Logical; if `TRUE`, prints progress every 1000 iterations.

## Value

A `list` with:

- `x_samples`: integer matrix \\S \times n\\ of raw labels (\\S =
  T\_{\mathrm{iter}}-T\_{\mathrm{burn}}\\).

- `lambda_samples`: list of length \\S\\; each element is a numeric
  vector of length \\L\_{\mathrm{cap}}\\ for that draw, with `NA` at
  empty labels.

- `K_per_iter`: integer vector length \\S\\ (occupied cluster count per
  saved draw).

- `L_cap_per_iter`: integer vector length \\S\\ (label capacity trace).

- `z_samples`: if `store_z=TRUE`, a numeric array \\S \times n \times
  n\\; otherwise `NULL`.

## Details

Row names of `w_ij` (if present) are propagated to item-level outputs;
otherwise items are named `Item_1, ..., Item_n`.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
n <- 6L
w <- matrix(0L, n, n)
w[lower.tri(w)] <- rpois(sum(lower.tri(w)), 2)
diag(w) <- 0
rownames(w) <- colnames(w) <- paste0("P", seq_len(n))
fit <- gibbs_bt_sbm(
  w_ij = w, prior = "GN", gamma_GN = 0.5,
  T_iter = 200, T_burn = 100, verbose = FALSE
)
} # }
```
