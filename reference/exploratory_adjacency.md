# Adjacency heatmap (wins per matches)

Exploratory adjacency heatmap showing wins for each ordered pair of
players. Optionally restrict to a subset of players and/or order
rows/cols by an external metric such as year-end ranking.

## Usage

``` r
exploratory_adjacency(
  Y_ij,
  N_ij = NULL,
  names = NULL,
  last_rank = NULL,
  order_by = c("auto", "none", "last_rank"),
  players_df = NULL,
  players_id_col = NULL,
  players_rank_col = NULL,
  labels = NULL,
  labels_key_col = NULL,
  labels_value_col = NULL,
  clean_display_labels = NULL,
  clean_fun = clean_players_names,
  max_wins = 4L,
  show_margins = TRUE,
  legend_title = "Number of wins",
  no_match_label = "No match played",
  no_match_color = "white",
  mark_unplayed = TRUE,
  unplayed_mark = "×",
  unplayed_mark_size = 2.2,
  bw_preview = TRUE
)
```

## Arguments

- Y_ij:

  Integer matrix of wins. Entry (i, j) is wins of i over j.

- N_ij:

  Integer matrix of matches. If `NULL`, computed as `Y_ij + t(Y_ij)`.

- names:

  Optional character vector of player names to include. Can be either
  raw matrix dimnames (slugs) or display names (after cleaning). If
  `NULL`, plots all players.

- last_rank:

  Optional numeric vector of year-end ranks. Can be named (names are
  player identifiers) or an unnamed vector aligned with the players
  being plotted.

- order_by:

  Ordering criterion. `"last_rank"` orders by ascending rank (rank 1 at
  the top). `"none"` keeps matrix order. `"auto"` uses `"last_rank"` if
  ranking info is provided, otherwise `"none"`.

- players_df:

  Optional data.frame with columns `player_slug` and `last_rank` (e.g.,
  season-level metadata). Used when `order_by = "last_rank"`.

- clean_fun:

  Name-cleaning function applied to display labels.

- max_wins:

  Maximum win count shown explicitly in the legend; larger counts are
  bucketed into `paste0(max_wins, "+")`.

- show_margins:

  If `TRUE` and package `ggside` is available, draws a side bar with
  total matches per player.

## Value

A `ggplot` object.
