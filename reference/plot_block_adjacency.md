# Figure 3 plotting function Block-ordered adjacency heatmap (BT-SBM)

Plot the adjacency (wins / matches) heatmap with rows/cols ordered by
the hard partition (`fit$estimates$x_hat`) and, within blocks, by
marginal wins. Draws block boundary lines and (optionally) a side color
strip for block IDs if ggside is installed. Optionally relabels blocks
so that Block 1 has the largest mean player strength \\\lambda\\ (mean
computed across players in the block).

## Usage

``` r
plot_block_adjacency(
  fit,
  w_ij,
  x_hat = NULL,
  relabel_blocks = c("none", "avg_lambda"),
  lambda_hat = NULL,
  relabel_decreasing = TRUE,
  clean_fun = clean_players_names,
  players_df = NULL,
  players_id_col = NULL,
  labels = NULL,
  labels_key_col = NULL,
  labels_value_col = NULL,
  clean_display_labels = NULL,
  order_ids = NULL,
  palette = NULL,
  fill_low = "#FFFFCC",
  fill_high = "#006400",
  no_match_label = "No match played",
  no_match_color = "white",
  mark_unplayed = TRUE,
  unplayed_mark = "×",
  unplayed_mark_size = 2.2,
  bw_preview = TRUE
)
```

## Arguments

- fit:

  Output list from `gibbs_BT_SBM()`.

- w_ij:

  Integer matrix of wins (same players & order used in `fit`).

- x_hat:

  partition point estimate. One n-length integer vector.

- clean_fun:

  Optional function to prettify player names. Default: identity.

- palette:

  Named colors for blocks (as character vector). Defaults to a
  Wimbledon-ish palette.

- fill_low, fill_high:

  Colors for the heatmap pixels gradient low/high.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) { # \dontrun{
p <- plot_block_adjacency(fit, w_ij)
} # }
```
