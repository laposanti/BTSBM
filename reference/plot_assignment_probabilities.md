# Figure 4 plotting function Assignment probability heatmap (BT-SBM)

Heatmap of item-wise posterior assignment probabilities for clusters
(relabeled so that Cluster 1 is the top block by decreasing
\\\lambda\\). Items are ordered by their most probable cluster and
marginal wins.

## Usage

``` r
plot_assignment_probabilities(
  fit,
  w_ij = NULL,
  max_n_clust = NULL,
  players_df = NULL,
  players_id_col = NULL,
  labels = NULL,
  labels_key_col = NULL,
  labels_value_col = NULL,
  clean_display_labels = NULL,
  clean_fun = clean_players_names,
  x_hat = NULL,
  order_ids = NULL,
  k_show = NULL,
  fill_low = "#FFFFCC",
  fill_high = "#006400"
)
```

## Arguments

- fit:

  Output list from `gibbs_BT_SBM()` (must include
  `relabeled$assign_prob`).

- w_ij:

  Optional wins matrix to compute marginal wins for ordering and
  annotation. If `NULL`, items are ordered by most-probable cluster
  only.

- max_n_clust:

  Where to filter the mcmc x_t. If not specified we use the modal K

- players_df:

  Optional data.frame for label lookup (e.g. season metadata).

- players_id_col:

  Optional column name in `players_df` that matches item IDs.

- labels:

  Optional character vector of display labels. Either named by item ID,
  or aligned with the items.

- labels_key_col:

  Optional column name in `players_df` to match item IDs.

- labels_value_col:

  Optional column name in `players_df` containing display labels.

- clean_display_labels:

  If `TRUE`, applies `clean_fun` to display labels.

- clean_fun:

  Optional function to prettify names. Default: identity.

- k_show:

  Optional integer number of clusters to show (defaults to all columns
  in `assign_prob`).

- fill_low, fill_high:

  Colors for the heatmap gradient low/high.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) { # \dontrun{
p <- plot_assignment_probabilities(fit, w_ij)
} # }
```
