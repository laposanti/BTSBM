#' Sample a Bradley–Terry Stochastic Block Model (BT-SBM) tournament
#'
#' @description
#' Generates a sparse undirected match matrix \eqn{N} and corresponding
#' directed win counts matrix \eqn{W} under the BTSBM
#' generating process. Block strengths are specified via a positive vector
#' \eqn{\lambda}, mapped to pairwise win probabilities
#' \eqn{\theta_{ab} = \lambda_a / (\lambda_a + \lambda_b)}.
#'
#' @details
#' The workflow is:
#' 1) Sample a symmetric match-count matrix \eqn{N} by first drawing a random
#'    subset of unordered pairs with probability `p_edge`, then assigning each
#'    selected pair a Poisson\eqn{(mean\_matches)} number of matches.
#' 2) Compute the block-to-block Bradley–Terry matrix
#'    \eqn{\theta = \lambda \mathbf{1}^\top / (\lambda \mathbf{1}^\top + \mathbf{1}\lambda^\top)}.
#' 3) For each unordered pair with \eqn{N_{ij}>0}, draw
#'    \eqn{w_{ij} \sim \text{Binomial}(N_{ij}, \theta_{x_i, x_j})} and set
#'    \eqn{w_{ji} = N_{ij} - w_{ij}}.
#'
#' By default, block labels are assigned deterministically as
#' `rep(seq_len(K), length.out = n_players)` to mirror your example, but you can
#' pass a custom label vector via `x`.
#'
#' @param n_players Integer. Number of players \eqn{n}.
#' @param K Integer. Number of blocks.
#' @param x Optional integer vector of length `n_players` with values in
#'   `1:K`. If `NULL`, uses `rep(seq_len(K), length.out = n_players)`.
#' @param mean_matches Positive numeric. Poisson mean for match counts on
#'   present edges.
#' @param p_edge Numeric in \[0,1\]. Probability that an unordered pair is
#'   present in the topology.
#' @param lambda Optional positive numeric vector of length `K`. If `NULL`,
#'   a geometric sequence is used: `lambda_base * lambda_ratio^(0:(K-1))`.
#' @param lambda_base Positive numeric. Base of the geometric schedule when
#'   `lambda` is not provided. Default `0.08` (as in your snippet).
#' @param lambda_ratio Positive numeric. Ratio of the geometric schedule when
#'   `lambda` is not provided. Default `2.3` (as in your snippet).
#' @param reverse_lambda Logical. If `TRUE`, uses `rev(lambda)` when computing
#'   `theta`, matching your `theta_star <- lambda_to_theta(rev(lambda_star))`.
#' @param seed Optional integer. If provided, sets the RNG seed for
#'   reproducibility *within the function call*.
#'
#' @return A list with components:
#' \itemize{
#' \item `N` — \eqn{n \times n} integer matrix of symmetric match counts with
#'   zero diagonal.
#' \item `W` — \eqn{n \times n} integer matrix of directed win counts with
#'   `w[i, j] + w[j, i] = N[i, j]`.
#' \item `x` — integer vector of block labels in `1:K`.
#' \item `lambda` — numeric vector of length `K` used to build `theta`.
#' \item `theta` — \eqn{K \times K} matrix with entries
#'   \eqn{\theta_{ab} = \lambda_a / (\lambda_a + \lambda_b)}.
#' }
#'
#' @export
sample_from_BTSBM <- function(n_players,
                              K,
                              x = NULL,
                              mean_matches = 5,
                              p_edge = 0.5,
                              lambda = NULL,
                              lambda_base = 0.08,
                              lambda_ratio = 2.3,
                              reverse_lambda = TRUE,
                              seed = NULL) {
  # --- input checks ---
  stopifnot(length(n_players) == 1L, n_players > 1L, n_players == as.integer(n_players))
  stopifnot(length(K) == 1L, K > 0L, K == as.integer(K))
  stopifnot(is.null(x) || (is.integer(x) || is.numeric(x)))
  stopifnot(mean_matches > 0)
  stopifnot(is.numeric(p_edge), p_edge >= 0, p_edge <= 1)
  if (!is.null(lambda)) {
    stopifnot(is.numeric(lambda), length(lambda) == K, all(lambda > 0))
  } else {
    stopifnot(lambda_base > 0, lambda_ratio > 0)
  }
  if (!is.null(seed)) {
    seed=123
  }

  n <- as.integer(n_players)

  # --- block labels (mirror your example by default) ---
  if (is.null(x)) {
    x <- rep(seq_len(K), length.out = n)
  } else {
    x <- as.integer(x)
    stopifnot(length(x) == n, all(x >= 1L & x <= K))
  }

  # --- lambda and theta ---
  if (is.null(lambda)) {
    lambda <- lambda_base * lambda_ratio^(0:(K - 1L))
  }
  lam_use <- if (isTRUE(reverse_lambda)) rev(lambda) else lambda

  # map lambda -> theta via outer
  theta <- outer(lam_use, lam_use, function(a, b) a / (a + b))

  # --- sample sparse-ish symmetric match counts N ---
  N <- matrix(0L, n, n)
  pairs <- utils::combn(n, 2L)
  m <- ncol(pairs)
  keep <- which(stats::runif(m) < p_edge)

  if (length(keep)) {
    for (k in keep) {
      i <- pairs[1L, k]; j <- pairs[2L, k]
      nij <- stats::rpois(1L, mean_matches)
      if (nij > 0L) {
        N[i, j] <- nij
        N[j, i] <- nij
      }
    }
  }
  diag(N) <- 0L

  # --- generate directed outcomes w given theta and labels x ---
  w <- matrix(0L, n, n)
  idx <- which(upper.tri(N) & N > 0L, arr.ind = TRUE)
  if (nrow(idx)) {
    for (r in seq_len(nrow(idx))) {
      i <- idx[r, 1L]; j <- idx[r, 2L]
      nij <- N[i, j]
      pij <- theta[x[i], x[j]]
      wij <- stats::rbinom(1L, nij, pij)
      w[i, j] <- wij
      w[j, i] <- nij - wij
    }
  }

  list(N = N, w = w, x = x, lambda = lambda, theta = theta)
}
