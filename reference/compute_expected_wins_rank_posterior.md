# Posterior rank summaries from player-strength draws

Convert posterior draws of player-level strength parameters into
posterior summaries for player ranks. Ranks are computed draw by draw
with rank 1 as the strongest player.

## Usage

``` r
compute_expected_wins_rank_posterior(
  lambda_item,
  w_ij = NULL,
  player_names = NULL,
  burn_in = 0,
  thin = 1,
  ci = 0.95,
  ties.method = "average",
  rank_grid_by = 0.5
)
```

## Arguments

- lambda_item:

  Numeric matrix or list of posterior draws for player strengths.

- w_ij:

  Optional wins matrix kept for API compatibility and ignored.

- player_names:

  Optional names for the players.

- burn_in:

  Burn-in as an integer iteration count or a fraction in `[0, 1)`.

- thin:

  Thinning interval.

- ci:

  Posterior interval level.

- ties.method:

  Passed to [`base::rank()`](https://rdrr.io/r/base/rank.html).

- rank_grid_by:

  Spacing for the posterior rank probability grid.

## Value

A list containing rank draws, a posterior rank-probability matrix, and
summary statistics for each player.
