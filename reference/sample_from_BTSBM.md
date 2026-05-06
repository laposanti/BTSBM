# Sample a Bradley–Terry Stochastic Block Model (BT-SBM) tournament

Generates a sparse undirected match matrix \\N\\ and corresponding
directed win counts matrix \\W\\ under the BTSBM generating process.
Block strengths are specified via a positive vector \\\lambda\\, mapped
to pairwise win probabilities \\\theta\_{ab} = \lambda_a / (\lambda_a +
\lambda_b)\\.

## Usage

``` r
sample_from_BTSBM(
  n_players,
  K,
  x = NULL,
  mean_matches = 5,
  p_edge = 0.5,
  lambda = NULL,
  lambda_base = 0.08,
  lambda_ratio = 2.3,
  reverse_lambda = TRUE,
  seed = NULL
)
```

## Arguments

- n_players:

  Integer. Number of players \\n\\.

- K:

  Integer. Number of blocks.

- x:

  Optional integer vector of length `n_players` with values in `1:K`. If
  `NULL`, uses `rep(seq_len(K), length.out = n_players)`.

- mean_matches:

  Positive numeric. Poisson mean for match counts on present edges.

- p_edge:

  Numeric in \[0,1\]. Probability that an unordered pair is present in
  the topology.

- lambda:

  Optional positive numeric vector of length `K`. If `NULL`, a geometric
  sequence is used: `lambda_base * lambda_ratio^(0:(K-1))`.

- lambda_base:

  Positive numeric. Base of the geometric schedule when `lambda` is not
  provided. Default `0.08` (as in your snippet).

- lambda_ratio:

  Positive numeric. Ratio of the geometric schedule when `lambda` is not
  provided. Default `2.3` (as in your snippet).

- reverse_lambda:

  Logical. If `TRUE`, uses `rev(lambda)` when computing `theta`,
  matching your `theta_star <- lambda_to_theta(rev(lambda_star))`.

- seed:

  Optional integer. If provided, sets the RNG seed for reproducibility
  *within the function call*.

## Value

A list with components:

- `N` — \\n \times n\\ integer matrix of symmetric match counts with
  zero diagonal.

- `W` — \\n \times n\\ integer matrix of directed win counts with
  `w[i, j] + w[j, i] = N[i, j]`.

- `x` — integer vector of block labels in `1:K`.

- `lambda` — numeric vector of length `K` used to build `theta`.

- `theta` — \\K \times K\\ matrix with entries \\\theta\_{ab} =
  \lambda_a / (\lambda_a + \lambda_b)\\.

## Details

The workflow is:

1.  Sample a symmetric match-count matrix \\N\\ by first drawing a
    random subset of unordered pairs with probability `p_edge`, then
    assigning each selected pair a Poisson\\(mean\\matches)\\ number of
    matches.

2.  Compute the block-to-block Bradley–Terry matrix \\\theta = \lambda
    \mathbf{1}^\top / (\lambda \mathbf{1}^\top +
    \mathbf{1}\lambda^\top)\\.

3.  For each unordered pair with \\N\_{ij}\>0\\, draw \\w\_{ij} \sim
    \text{Binomial}(N\_{ij}, \theta\_{x_i, x_j})\\ and set \\w\_{ji} =
    N\_{ij} - w\_{ij}\\.

By default, block labels are assigned deterministically as
`rep(seq_len(K), length.out = n_players)` to mirror your example, but
you can pass a custom label vector via `x`.
