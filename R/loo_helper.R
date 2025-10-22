#' Log-likelihood matrix for the *simple* Bradley–Terry model
#'
#' Builds an T_iter x D matrix of log-likelihood values, where T_iter is the number of
#' posterior draws and D is the number of observed unordered pairs (i<j) with
#' \code{n_ij > 0}. This is suitable as input to \pkg{loo}.
#'
#' @param w_ij integer/numeric n x n wins (i over j).
#' @param lambda_samples numeric T_iter x n matrix of player-specific rates \eqn{\lambda_i}.
#'
#' @return A list with:
#' \itemize{
#' \item \code{ll} — T_iter x D matrix of log-likelihoods.
#' \item \code{obs_idx} — D x 2 matrix of (i,j) indices defining each column.
#' }
#'
#' @export
#'
make_bt_simple_loo <- function(w_ij, lambda_samples) {
  n_ij = w_ij + t(w_ij)
  stopifnot(is.matrix(w_ij), is.matrix(n_ij), is.matrix(lambda_samples))
  K <- nrow(w_ij)
  if (!identical(dim(w_ij), dim(n_ij))) stop("w_ij and n_ij must have same dims.")
  if (!isTRUE(all.equal(n_ij, t(n_ij)))) stop("n_ij must be symmetric.")
  if (any(diag(n_ij) != 0)) stop("Diagonal of n_ij must be zero.")
  T_iter <- nrow(lambda_samples)
  if (ncol(lambda_samples) != K) stop("lambda_samples must have K columns.")

  idx <- which(upper.tri(n_ij) & n_ij > 0, arr.ind = TRUE)
  D <- nrow(idx)
  ll <- matrix(NA_real_, T_iter, D)

  for (s in seq_len(T_iter)) {
    lam <- lambda_samples[s, ]
    for (d in seq_len(D)) {
      i <- idx[d, 1]; j <- idx[d, 2]
      n <- n_ij[i, j]; w <- w_ij[i, j]
      p <- lam[i] / (lam[i] + lam[j])
      ll[s, d] <- stats::dbinom(w, size = n, prob = p, log = TRUE)
    }
  }
  list(ll = ll, obs_idx = idx)
}

#' Log-likelihood matrix for the BT–SBM (clustered) model
#'
#' Builds an T_iter x D matrix of log-likelihood values using cluster labels \eqn{x_i}
#' and cluster rates \eqn{\lambda_k}. Assumes \code{x_samples} and
#' \code{lambda_samples} are *relabelled consistently* (e.g. via \code{inference_helper}).
#'
#' @param w_ij integer/numeric n x n wins (i over j).
#' @param lambda_samples numeric T_iter x n matrix of cluster rates \eqn{\lambda_k}.
#' @param x_samples integer T_iter x n matrix of cluster labels for each player.
#'
#' @return A list with:
#' \itemize{
#' \item \code{ll} — T_iter x D matrix of log-likelihoods.
#' \item \code{obs_idx} — D x 2 matrix of (i,j) indices defining each column.
#' }
#'
#' @examples
#' \dontrun{
#' # After running your clustered sampler and relabeling:
#' # ll_obj <- make_bt_cluster_loo(w, n, out$lambda_samples_relabel, out$x_samples_relabel)
#' }
#' @export
# Drop-in replacement: accepts lambda_samples as either
#  - matrix S x Kmax (with NA padding), or
#  - list length S with numeric vectors of length K_s
make_bt_cluster_loo <- function(w_ij, lambda_samples, x_samples) {
  stopifnot(is.matrix(w_ij), is.matrix(x_samples))
  n <- nrow(w_ij); stopifnot(n == ncol(w_ij))
  S <- nrow(x_samples); stopifnot(S > 0)

  # Accessor for lambda at iteration s (length K_s)
  if (is.list(lambda_samples)) {
    get_lam <- function(s) {
      lam <- lambda_samples[[s]]
      if (!is.numeric(lam)) stop("lambda_samples[[s]] must be numeric.")
      lam
    }
    Kmax <- max(vapply(lambda_samples, length, 1L))
  } else if (is.matrix(lambda_samples)) {
    stopifnot(nrow(lambda_samples) == S)
    get_lam <- function(s) {
      lam <- lambda_samples[s, ]
      lam[is.finite(lam)]  # allow NA padding on the right
    }
    Kmax <- ncol(lambda_samples)
  } else stop("lambda_samples must be list(S) or matrix S x Kmax.")

  n_ij <- w_ij + t(w_ij)
  if (!isTRUE(all.equal(n_ij, t(n_ij)))) stop("n_ij must be symmetric.")
  if (any(diag(n_ij) != 0)) stop("Diagonal of n_ij must be zero.")

  idx <- which(upper.tri(n_ij) & n_ij > 0, arr.ind = TRUE)   # MATCH simple builder
  D <- nrow(idx)
  ll <- matrix(NA_real_, S, D)

  for (s in seq_len(S)) {
    xs <- as.integer(x_samples[s, ])
    lam <- get_lam(s)                   # length K_s
    Ks  <- length(lam)
    if (any(xs < 1 | xs > Ks, na.rm = TRUE))
      stop(sprintf("x_samples row %d has labels outside 1..K_s (K_s=%d).", s, Ks))
    lam_player <- lam[xs]               # length n

    for (d in seq_len(D)) {
      i <- idx[d, 1]; j <- idx[d, 2]
      nij <- n_ij[i, j]; wij <- w_ij[i, j]
      pij <- lam_player[i] / (lam_player[i] + lam_player[j])
      ll[s, d] <- stats::dbinom(wij, size = nij, prob = pij, log = TRUE)
    }
  }
  list(ll = ll, obs_idx = idx)
}


#' Compare BT models with Pareto-smoothed importance sampling LOO
#'
#' Convenience wrapper around \pkg{loo} for comparing two log-likelihood
#' matrices (simple vs clustered BT).
#'
#' @param simple_llo list returned by \code{make_bt_simple_loo()}.
#' @param cluster_llo list returned by \code{make_bt_cluster_loo()}.
#'
#' @return A list with:
#' \itemize{
#' \item \code{simple} — \code{loo} object for the simple BT.
#' \item \code{cluster} — \code{loo} object for the clustered BT–SBM.
#' \item \code{comparison} — result of \code{loo::compare_models()}.
#' }
#' @export
#'
compare_bt_models_loo <- function(simple_llo, cluster_llo) {
  if (!requireNamespace("loo", quietly = TRUE))
    stop("Package 'loo' not installed.")
  for (x in list(simple_llo, cluster_llo)) {
    if (!is.list(x) || is.null(x$ll)) stop("Each input must be a list with $ll.")
  }
  loo_simple  <- loo::loo(simple_llo$ll)
  loo_cluster <- loo::loo(cluster_llo$ll)
  comp <- loo::loo_compare(loo_simple, loo_cluster)
  list(simple = loo_simple, cluster = loo_cluster, comparison = comp)
}

