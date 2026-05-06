# Build BT clustered log-likelihood matrix (player-level)

Convenience alternative returning a plain \\S \times D\\ matrix (\\D\\ =
\#pairs with i\<j). Here `x_draws` are labels and `lambda_draws` are
cluster rates per draw.

## Usage

``` r
make_bt_loo_cluster(x_draws, lambda_draws, w_ij)
```

## Arguments

- x_draws:

  integer \\S \times n\\ matrix of cluster labels per draw.

- lambda_draws:

  numeric \\S \times K\\ matrix of cluster rates per draw.

- w_ij:

  \\n \times n\\ pairwise wins matrix (i over j, diag = 0).

## Value

Numeric matrix \\S \times D\\ of log-likelihoods.
