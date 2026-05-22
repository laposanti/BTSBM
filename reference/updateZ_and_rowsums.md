# Update cluster assignments and row sums

Internal helper used by the Gibbs sampler to update latent cluster
assignments and the corresponding row-sum summaries.

## Usage

``` r
updateZ_and_rowsums(n_ij, x, lambda)
```

## Arguments

- n_ij:

  Integer matrix of pairwise match counts.

- x:

  Integer vector of current cluster assignments.

- lambda:

  Numeric vector of current item strengths.
