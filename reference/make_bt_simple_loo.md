# Log-likelihood matrix for the *simple* Bradley–Terry model

Builds an T_iter x D matrix of log-likelihood values, where T_iter is
the number of posterior draws and D is the number of observed unordered
pairs (i\<j) with `n_ij > 0`. This is suitable as input to loo.

Builds an \\S \times D\\ matrix of log-likelihood values, where \\S\\ is
the number of posterior draws and \\D\\ is the number of observed
unordered pairs (i\<j) with \\n\_{ij} \> 0\\.

## Usage

``` r
make_bt_simple_loo(w_ij, lambda_samples)

make_bt_simple_loo(w_ij, lambda_samples)
```

## Arguments

- w_ij:

  integer/numeric \\n \times n\\ wins (i over j).

- lambda_samples:

  numeric \\S \times n\\ matrix of player-specific rates \\\lambda_i\\.

## Value

A list with:

- `ll` — T_iter x D matrix of log-likelihoods.

- `obs_idx` — D x 2 matrix of (i,j) indices defining each column.

A list with:

- `ll` — \\S \times D\\ matrix of log-likelihoods.

- `obs_idx` — \\D \times 2\\ matrix of (i,j) indices defining each
  column.
