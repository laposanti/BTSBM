# Figure 5 plotting function Lambda uncertainty plot (per player)

Forest plot of \\\lambda\\ with uncertainty intervals, using relabeled
draws. Points colored by the hard partition (`fit$estimates$x_hat`).

## Usage

``` r
plot_lambda_uncertainty(
  fit,
  w_ij = NULL,
  players_df = NULL,
  players_id_col = NULL,
  labels = NULL,
  labels_key_col = NULL,
  labels_value_col = NULL,
  clean_display_labels = NULL,
  order_ids = NULL,
  log_base = exp(1),
  max_n_clust = NULL,
  prob = 0.9,
  palette = NULL,
  clean_fun = function(x) x,
  x_hat = NULL,
  filter_lambdas = TRUE,
  conditional = TRUE,
  ...
)
```

## Arguments

- fit:

  Output from `gibbs_BT_SBM()` with `opt_lambda$lambda_item` computed
  (set `keep_lambda=TRUE` when sampling).

- w_ij:

  Optional wins matrix to compute marginal wins for ordering and
  annotation.

- players_df:

  Optional data.frame used to infer display labels.

- players_id_col:

  Optional id column name in `players_df`.

- labels:

  Optional display labels; either named by item ids or aligned with
  items.

- labels_key_col:

  Optional key column in `players_df` used to match item ids.

- labels_value_col:

  Optional value column in `players_df` containing display labels.

- clean_display_labels:

  Logical; if `TRUE`, applies `clean_fun` to resolved display labels.

- order_ids:

  Optional explicit player ordering.

- log_base:

  Base for the x-axis logarithm (10 or e). Defaults to 10.

- max_n_clust:

  Where to filter the mcmc x_t. If not specified we use the modal K

- prob:

  Interval probability for HPD (e.g., 0.90).

- palette:

  Named colors for clusters.

- clean_fun:

  Optional player-name cleaner. Default: identity.

- x_hat:

  Optional hard partition; if `NULL`, inferred from `fit`.

- filter_lambdas:

  Logical; when `conditional = TRUE`, keeps lambda draws aligned to
  filtered `x` draws.

- conditional:

  Logical; if `TRUE`, intervals are computed conditional on each item's
  hard cluster.

- ...:

  Reserved for future extensions.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) { # \dontrun{
# fit with keep_lambda=TRUE
p <- plot_lambda_uncertainty(fit, prob = 0.90)
} # }
```
