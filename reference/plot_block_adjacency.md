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
  unplayed_mark_size = 2.2
)
```

## Arguments

- fit:

  Output list from `gibbs_BT_SBM()`.

- w_ij:

  Integer matrix of wins (same players & order used in `fit`).

- x_hat:

  partition point estimate. One n-length integer vector.

- relabel_blocks:

  Either `"none"` or `"avg_lambda"`; controls optional relabeling of
  `x_hat`.

- lambda_hat:

  Optional lambda draws/vector used when
  `relabel_blocks = "avg_lambda"`.

- relabel_decreasing:

  Logical; direction used in `avg_lambda` relabeling.

- clean_fun:

  Optional function to prettify player names. Default: identity.

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

- palette:

  Named colors for blocks (as character vector). Defaults to a
  Wimbledon-ish palette.

- fill_low, fill_high:

  Colors for the heatmap pixels gradient low/high.

- no_match_label:

  Label used for pairs with zero observed matches.

- no_match_color:

  Fill colour for pairs with zero observed matches.

- mark_unplayed:

  Logical; if `TRUE`, overlays `unplayed_mark` on unplayed pairs.

- unplayed_mark:

  Character marker used for unplayed pairs.

- unplayed_mark_size:

  Numeric marker size for unplayed pairs.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) { # \dontrun{
p <- plot_block_adjacency(fit, w_ij)
} # }
```
