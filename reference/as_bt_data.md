# Create validated Bradley–Terry pairwise data

Converts a square directed win-count matrix into the package's validated
pairwise-data representation. This is the boundary between the legacy
matrix API and the common fitting interface.

## Usage

``` r
as_bt_data(wins, item_labels = NULL)
```

## Arguments

- wins:

  Numeric square matrix. Entry `wins[i, j]` is the number of wins for
  item `i` over item `j`.

- item_labels:

  Optional character vector of item names. By default, matrix row names
  are used when available.

## Value

An object of class `btsbm_pairwise_data`.
