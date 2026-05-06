# Log-likelihood matrix for the BT–SBM (clustered) model

Builds an T_iter x D matrix of log-likelihood values using cluster
labels \\x_i\\ and cluster rates \\\lambda_k\\. Assumes `x_samples` and
`lambda_samples` are *relabelled consistently* (e.g. via
`inference_helper`).

Builds an \\S \times D\\ matrix using cluster labels \\x_i\\ and cluster
rates \\\lambda_k\\. Assumes inputs are *relabelled consistently* or
that cluster ids in `x_samples[s, ]` are in `1..K` where
`K = ncol(lambda_samples)`.

## Usage

``` r
make_bt_cluster_loo(w_ij, lambda_samples, x_samples)

make_bt_cluster_loo(w_ij, lambda_samples, x_samples)
```

## Arguments

- w_ij:

  integer/numeric \\n \times n\\ wins (i over j).

- lambda_samples:

  numeric \\S \times K\\ matrix of cluster rates \\\lambda_k\\.

- x_samples:

  integer \\S \times n\\ matrix of cluster labels for each item.

## Value

A list with:

- `ll` — T_iter x D matrix of log-likelihoods.

- `obs_idx` — D x 2 matrix of (i,j) indices defining each column.

A list with:

- `ll` — \\S \times D\\ matrix of log-likelihoods.

- `obs_idx` — \\D \times 2\\ matrix of (i,j) indices defining each
  column.

## Examples

``` r
if (FALSE) { # \dontrun{
# After running your clustered sampler and relabeling:
# ll_obj <- make_bt_cluster_loo(w, n, out$lambda_samples_relabel, out$x_samples_relabel)
} # }
```
