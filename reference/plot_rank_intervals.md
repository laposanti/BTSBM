# Plot posterior rank intervals

Draw a rank-interval plot from posterior rank summaries, posterior
lambda draws, or the relabelled output returned by
[`relabel_by_lambda()`](https://laposanti.github.io/BTSBM/reference/relabel_by_lambda.md).

## Usage

``` r
plot_rank_intervals(
  x,
  max_players = NULL,
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

- x:

  Either a posterior-rank summary, posterior lambda draws, or the output
  of
  [`relabel_by_lambda()`](https://laposanti.github.io/BTSBM/reference/relabel_by_lambda.md).

- max_players:

  Optional maximum number of players to display.

- w_ij:

  Optional wins matrix passed through when rank summaries must be
  computed from lambda draws.

- player_names:

  Optional player names passed through to
  [`compute_expected_wins_rank_posterior()`](https://laposanti.github.io/BTSBM/reference/compute_expected_wins_rank_posterior.md).

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

A `ggplot` object.
