# Create validated Plackett–Luce ranking data

Rows of `rankings` are preference orderings and columns are successive
positions, from most to least preferred. Each row may be a complete
ranking or a top-`m` prefix from one common universe of items.
Internally this is represented in the compact integer format used by the
PLuce reference code, but users never need to work with that internal
name.

## Usage

``` r
as_rankings(
  rankings,
  item_count = NULL,
  ranking_labels = NULL,
  item_labels = NULL
)

as_pl_data(rho, n_items = NULL, ranking_id = NULL, item_labels = NULL)
```

## Arguments

- rankings:

  Integer-like matrix. Row `l` lists item ids from most to least
  preferred in ranking `l`.

- item_count:

  Total number of available items. If `NULL`, it is inferred as
  `max(rankings)`; supply it when a top-`m` dataset has never-ranked
  items.

- ranking_labels:

  Optional unique character vector identifying ranking rows.

- item_labels:

  Optional unique character vector of length `item_count`.

- rho:

  Deprecated name for `rankings`.

- n_items:

  Deprecated name for `item_count`.

- ranking_id:

  Deprecated name for `ranking_labels`.

## Value

An object of class `btsbm_ranking_data`.

## References

Caron, F. and Doucet, A. (2012). Efficient Bayesian inference for
generalized Bradley–Terry models. *Journal of Computational and
Graphical Statistics*, 21(1), 174–196.
[doi:10.1080/10618600.2012.638220](https://doi.org/10.1080/10618600.2012.638220)
