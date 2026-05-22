# Relabel partitions by decreasing \\\lambda\\, compute point estimates, and quantify uncertainty via credible balls

Given posterior samples of labels `x_samples` (S x N) and corresponding
cluster-level intensities `lambda_samples` per iteration, this function:
(i) relabels each draw so that cluster `1` has the largest \\\lambda\\,
cluster `2` the second largest, etc.; (ii) computes posterior similarity
matrix (PSM) and partition point estimates (minVI and Binder); and (iii)
constructs a VI-credible ball around the minVI partition and returns the
*extremal* partitions on its surface (lower/upper), relabelled by
decreasing strength for interpretability. It also returns per-item
assignment probabilities and the posterior distribution of the number of
clusters.

## Usage

``` r
relabel_by_lambda(x_samples, lambda_samples)
```

## Arguments

- x_samples:

  Integer matrix (S x N) of sampled partitions; rows index MCMC
  iterations and columns index items. Cluster labels are arbitrary
  across iterations.

- lambda_samples:

  Either:

  - a **list** of length S; element `[[s]]` is a numeric vector of
    cluster intensities \\\lambda\_\ell\\ indexed by the *raw* label id
    used at iteration `s` (NAs allowed for non-occupied ids), or

  - a **matrix** (S x L) whose row `s` gives \\\lambda\_\ell\\ for label
    \\\ell\\ at iteration `s` (sparse columns with NAs permitted).

## Value

A named list with components:

- `x_samples_relabel`:

  Integer matrix (S x N) of relabelled draws (labels `1..K_s` per
  iteration `s`, ordered by decreasing \\\lambda\\).

- `lambda_samples_relabel`:

  Numeric matrix (S x N) assigning each item its cluster's \\\lambda\\
  after relabelling.

- `co_clustering`:

  Posterior similarity matrix (N x N).

- `minVI_partition`:

  Partition estimated by minimizing posterior expected VI using salso.

- `partition_binder`:

  Partition estimated by Binder's loss using salso.

- `n_clusters_each_iter`:

  Integer vector (length S) of occupied cluster counts per iteration.

- `block_count_distribution`:

  Data frame with columns `num_blocks`, `count`, `prob` summarizing the
  posterior of the number of clusters.

- `item_cluster_assignment_probs`:

  Data frame (N x Kmax) of per-item marginal assignment probabilities
  (columns `Cluster_1`, `Cluster_2`, ...).

- `avg_top_block_count`:

  Average size of the top-\\\lambda\\ cluster across iterations.

- `top_block_count_per_iter`:

  Integer vector (length S) with the size of the top-\\\lambda\\ cluster
  per iteration.

- `top_block_size_distribution`:

  Data frame with columns `top_block_size`, `count`, `prob` summarizing
  the posterior distribution of the strongest block size.

- `cluster_lambda_ordered`:

  List of length S; each element is the vector of cluster \\\lambda\\
  values for that iteration, ordered decreasingly.

## Details

**Relabelling.** For each iteration, occupied labels are compacted to
`1..K_s` and reordered by decreasing \\\lambda\\, producing a canonical
“1 = strongest” labelling.

**Point estimation.** The posterior similarity matrix is computed from
relabelled draws; the minVI and Binder partitions are estimated with
salso using the generalized variation of information and Binder losses.

## Input requirements

- `x_samples` must be integer-valued with no missing items per row.

- `lambda_samples` may be sparse (NAs for non-occupied labels).

- salso must be available for the point-estimate summaries.

## References

Wade, S., 2023. Bayesian cluster analysis. Philosophical Transactions of
the Royal Society A: Mathematical, Physical and Engineering Sciences
381, 20220149. https://doi.org/10.1098/rsta.2022.0149

## See also

[`salso`](https://rdrr.io/pkg/salso/man/salso.html),
[`VI`](https://rdrr.io/pkg/salso/man/partition.loss.html),
[`binder`](https://rdrr.io/pkg/salso/man/partition.loss.html),
[`comp.psm`](https://rdrr.io/pkg/mcclust/man/comp.psm.html)

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
S <- 50; N <- 15
x_samps <- matrix(sample(1:4, S*N, TRUE), S, N)
# Sparse lambda per-iter: labels up to 6, only first 4 occupied typically
lam_list <- replicate(S, { v <- rep(NA_real_, 6); v[1:4] <- rexp(4, 1); v }, simplify = FALSE)
out <- relabel_by_lambda(x_samps, lam_list)
out$minVI_partition[1:10]
out$K_VI_lower; out$K_VI_upper
} # }
```
