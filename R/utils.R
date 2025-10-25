# R/utils-internal.R

#' @keywords internal
.safe_log <- function(x) ifelse(x > 0, log(x), -Inf)

#' @keywords internal
.logsumexp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

#' Slice sampler on log-scale (internal)
#'
#' Univariate slice sampling on \eqn{\log a}. Used internally for hyperparameter updates.
#'
#' @param logpost Function taking `loga` and returning log-posterior up to a constant.
#' @param loga0 Numeric; current log value.
#' @param w Step-out width (default 1).
#' @param m Max step-out steps (default 20).
#' @param lower,upper Hard bounds on `loga`.
#' @return A new `loga` draw.
#' @keywords internal
.slice_on_log <- function(
    logpost, loga0, w = 1.0, m = 20L, lower = -Inf, upper = Inf
) {
  f0 <- logpost(loga0)
  if (!is.finite(f0)) stop("Initial loga has -Inf logpost.")
  y <- f0 - stats::rexp(1)
  L <- loga0 - stats::runif(1, 0, w)
  R <- L + w
  J <- floor(stats::runif(1) * m); K <- (m - 1L) - J
  while (J > 0L && L > lower && logpost(L) > y) { L <- L - w; J <- J - 1L }
  while (K > 0L && R < upper && logpost(R) > y) { R <- R + w; K <- K - 1L }
  for (iter in 1:100) {
    loga1 <- stats::runif(1, max(L, lower), min(R, upper))
    f1 <- logpost(loga1)
    if (is.finite(f1) && f1 >= y) return(loga1)
    if (loga1 < loga0) L <- loga1 else R <- loga1
  }
  warning("Slice sampler reached max shrink steps; returning last draw.")
  loga1
}

#' Gamma-Poisson integrated log-likelihood (new block)
#' @keywords internal
.new_cluster_integral_log <- function(wi, Zi, a_now, b_now) {
  (a_now * log(b_now)) + lgamma(a_now + wi) - lgamma(a_now) -
    (a_now + wi) * log(b_now + Zi)
}

#' Anchored version of the integrated log-score (drop a-only constant)
#' @keywords internal
.new_cluster_integral_log_anchored <- function(wi, Zi, a_now, b_now) {
  lgamma(a_now + wi) - (a_now + wi) * log(b_now + Zi)
}

#' Sample index from log-weights
#' @keywords internal
.sample_from_logweights <- function(logw) {
  if (all(!is.finite(logw))) return(sample.int(length(logw), 1L))
  lse <- .logsumexp(logw)
  if (!is.finite(lse)) return(sample.int(length(logw), 1L))
  p <- exp(logw - lse); p[!is.finite(p)] <- 0
  s <- sum(p)
  if (s <= 0) sample.int(length(p), 1L) else sample.int(length(p), 1L, prob = p/s)
}

# ---------- Naming helpers ----------

#' Format player/item names as "Surname F."
#'
#' @param name character scalar or vector (e.g., "Roger Federer").
#' @return Character vector with formatted names (e.g., "Federer R.").
#' @examples
#' clean_players_names(c("Roger Federer", "Nadal"))
#' @export
clean_players_names <- function(name) {
  name_parts <- strsplit(name, " ")
  vapply(name_parts, function(parts) {
    parts <- unlist(parts)
    if (length(parts) == 1L) return(parts)
    paste(parts[length(parts)], substr(parts[1], 1, 1), ".", sep = " ")
  }, FUN.VALUE = character(1))
}

# ---------- Priors & urn weights ----------

#' Gnedin prior normalization weight K(n,k)
#'
#' Computes the unnormalized mass function term used in Gnedin-type priors.
#'
#' @param n integer(1) total items.
#' @param k integer vector of block counts.
#' @param gamma numeric(1) parameter \eqn{\gamma > 0}.
#' @return Numeric vector of weights \eqn{K(n,k)}.
#' @keywords internal
#' @references
#' Legramanti, S., Rigon, T., Durante, D., Dunson, D.B., 2022. Extended stochastic block models with application to criminal networks. The Annals of Applied Statistics 16. https://doi.org/10.1214/21-aoas1595

HGnedin <- function(n, k, gamma = 0.5) {
  exp(lchoose(n, k) + lgamma(k - gamma) - lgamma(1 - gamma) +
        log(gamma) + lgamma(n + gamma - k) - lgamma(n + gamma))
}

#' Urn weight: Dirichlet–Multinomial with finite cap K_DM
#' @param v_minus integer sizes of occupied clusters (excluding focal item).
#' @param beta_DM numeric concentration (>0).
#' @param K_DM integer maximum number of clusters (>=1).
#' @return Numeric vector length \eqn{K+1}: existing weights and new-cluster mass.
#' #' @references
#' Legramanti, S., Rigon, T., Durante, D., Dunson, D.B., 2022. Extended stochastic block models with application to criminal networks. The Annals of Applied Statistics 16. https://doi.org/10.1214/21-aoas1595

#' @keywords internal
urn_DM <- function(v_minus, beta_DM, K_DM) {
  K <- length(v_minus)
  c(v_minus + beta_DM, beta_DM * (K_DM - K) * (K_DM > K))
}

#' Urn weight: Dirichlet Process
#' @param v_minus integer sizes of occupied clusters (excluding focal item).
#' @param alpha_PY numeric concentration (>0).
#' @return Numeric vector length \eqn{K+1}: existing weights and new-cluster mass.
#' #' @references
#' Legramanti, S., Rigon, T., Durante, D., Dunson, D.B., 2022. Extended stochastic block models with application to criminal networks. The Annals of Applied Statistics 16. https://doi.org/10.1214/21-aoas1595

#' @keywords internal
urn_DP <- function(v_minus, alpha_PY) {
  c(v_minus, alpha_PY)
}

#' Urn weight: Gnedin heuristic
#' @param v_minus integer sizes of occupied clusters (excluding focal item).
#' @param gamma numeric parameter (>0).
#' @return Numeric vector length \eqn{K+1}.
#' #' @references
#' Legramanti, S., Rigon, T., Durante, D., Dunson, D.B., 2022. Extended stochastic block models with application to criminal networks. The Annals of Applied Statistics 16. https://doi.org/10.1214/21-aoas1595

#' @keywords internal
urn_GN <- function(v_minus, gamma) {
  K  <- length(v_minus)
  n_ <- sum(v_minus) + 1L               # +1 for the focal item
  c((v_minus + 1) * (n_ - K + gamma),
    K^2 - K * gamma)
}

#' Urn weight: Pitman–Yor process
#' @param v_minus integer sizes of occupied clusters (excluding focal item).
#' @param alpha_PY numeric concentration (> -sigma_PY).
#' @param sigma_PY numeric discount in [0,1).
#' @return Numeric vector length \eqn{K+1}.
#' #' @references
#' Legramanti, S., Rigon, T., Durante, D., Dunson, D.B., 2022. Extended stochastic block models with application to criminal networks. The Annals of Applied Statistics 16. https://doi.org/10.1214/21-aoas1595

#' @keywords internal
urn_PY <- function(v_minus, alpha_PY, sigma_PY) {
  K <- length(v_minus)
  c(v_minus - sigma_PY, alpha_PY + K * sigma_PY)
}

#' Expected number of clusters under finite/inf PY/DM-like priors (helper)
#'
#' @param n integer(1) number of items.
#' @param sigma numeric(1) discount (0 for DP/DM).
#' @param theta numeric(1) concentration parameter.
#' @param K_DM integer(1) maximum number of clusters, or \code{Inf}.
#' @return Numeric(1) expected number of clusters.
#' #' @references
#' Legramanti, S., Rigon, T., Durante, D., Dunson, D.B., 2022. Extended stochastic block models with application to criminal networks. The Annals of Applied Statistics 16. https://doi.org/10.1214/21-aoas1595

#' @keywords internal
expected_cl_py <- function(n, sigma, theta, K_DM) {
  n <- as.integer(n)
  stopifnot(sigma >= 0, sigma < 1, theta > - sigma, n > 0, K_DM > 1)
  if (is.infinite(K_DM)) {
    if (sigma == 0) {
      out <- theta * sum(1 / (theta - 1 + 1:n))
    } else {
      out <- (1 / sigma) * exp(lgamma(theta + sigma + n) - lgamma(theta + sigma) -
                                 lgamma(theta + n) + lgamma(theta + 1)) - theta / sigma
    }
  } else {
    if (sigma == 0) {
      index <- 0:(n - 1L)
      out <- K_DM - K_DM * exp(sum(log(index + theta * (1 - 1 / K_DM)) - log(theta + index)))
    } else {
      out <- NA_real_  # finite-cap PY not provided here
    }
  }
  out
}

# ---------- Simple BT sampler (no clustering) ----------

#' Simple Bradley–Terry Gibbs sampler (no clustering)
#'
#' @param w_ij integer/numeric \eqn{n \times n} wins from i over j (diag = 0, nonnegative).
#' @param a,b numeric(1) Gamma(a,b) prior on each \eqn{\lambda_i}.
#' @param T_iter,T_burn integers; total iterations and burn-in. Require \code{T_burn < T_iter}.
#' @param verbose logical; print progress every 1000 iterations.
#'
#' @return A list with \code{lambda_samples} (matrix of size \eqn{S \times n}, \eqn{S=T_{\text{iter}}-T_{\text{burn}}}).
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 6L
#' w <- matrix(0L, n, n)
#' w[upper.tri(w)] <- rpois(sum(upper.tri(w)), 2)
#' w <- w + t(w) - diag(diag(w))
#' fit <- gibbs_bt_simple(w, a = 1, b = 1, T_iter = 500, T_burn = 100, verbose = FALSE)
#' }
#' @importFrom stats rgamma
#' @export
gibbs_bt_simple <- function(
    w_ij,
    a = 0.01, b = 0.1,
    T_iter = 5000, T_burn = 1000,
    verbose = TRUE
) {
  stopifnot(is.matrix(w_ij))
  n <- nrow(w_ij)
  stopifnot(n == ncol(w_ij))
  if (any(w_ij < 0)) stop("w_ij: negative entries not allowed.")
  if (any(diag(w_ij) != 0)) stop("w_ij: diagonal must be zero.")
  if (T_burn >= T_iter) stop("Require T_burn < T_iter.")

  n_ij <- w_ij + t(w_ij)

  w_i <- rowSums(w_ij)
  lambda <- stats::rgamma(n, a, b)
  Z      <- matrix(0, n, n)

  S <- T_iter - T_burn
  lambda_store <- matrix(NA_real_, S, n)
  keep_idx <- 0L

  for (iter in seq_len(T_iter)) {
    # Polya–Gamma–like gamma trick for BT via latent exponential/gamma:
    for (i in 1:(n - 1L)) {
      for (j in (i + 1L):n) {
        nij <- n_ij[i, j]
        if (nij > 0) {
          rate   <- lambda[i] + lambda[j]
          Z_ij   <- stats::rgamma(1, nij, rate)
          Z[i, j] <- Z[j, i] <- Z_ij
        } else {
          Z[i, j] <- Z[j, i] <- 0
        }
      }
    }

    for (i in seq_len(n)) {
      shape <- a + w_i[i]
      rate  <- b + sum(Z[i, ])
      lambda[i] <- stats::rgamma(1, shape, rate)
    }

    if (iter > T_burn) {
      keep_idx <- keep_idx + 1L
      lambda_store[keep_idx, ] <- lambda
    }

    if (verbose && iter %% 1000L == 0L)
      cat("simple BT | iter", iter, "\n")
  }

  invisible(list(lambda_samples = lambda_store))
}

# ---------- LOO builders & model comparison ----------

#' Log-likelihood matrix for the *simple* Bradley–Terry model
#'
#' Builds an \eqn{S \times D} matrix of log-likelihood values, where \eqn{S}
#' is the number of posterior draws and \eqn{D} is the number of observed
#' unordered pairs (i<j) with \eqn{n_{ij} > 0}.
#'
#' @param w_ij integer/numeric \eqn{n \times n} wins (i over j).
#' @param lambda_samples numeric \eqn{S \times n} matrix of player-specific rates \eqn{\lambda_i}.
#'
#' @return A list with:
#' \itemize{
#' \item \code{ll} — \eqn{S \times D} matrix of log-likelihoods.
#' \item \code{obs_idx} — \eqn{D \times 2} matrix of (i,j) indices defining each column.
#' }
#' @importFrom stats dbinom
#' @export
make_bt_simple_loo <- function(w_ij, lambda_samples) {
  stopifnot(is.matrix(w_ij), is.matrix(lambda_samples))
  n <- nrow(w_ij); stopifnot(n == ncol(w_ij))
  if (any(diag(w_ij) != 0, na.rm = TRUE)) stop("Diagonal of w_ij must be zero (ignoring NAs).")
  S <- nrow(lambda_samples)
  if (ncol(lambda_samples) != n) stop("lambda_samples must have n columns.")

  # Symmetric match counts; may contain NAs if w_ij has NAs
  n_ij <- w_ij + t(w_ij)

  # Only keep pairs with finite, positive counts
  idx <- which(upper.tri(n_ij) & is.finite(n_ij) & (n_ij > 0), arr.ind = TRUE)
  D <- nrow(idx)
  ll <- matrix(NA_real_, S, D)

  for (s in seq_len(S)) {
    lam <- lambda_samples[s, ]
    for (d in seq_len(D)) {
      i <- idx[d, 1]; j <- idx[d, 2]
      nij <- n_ij[i, j]
      wij <- w_ij[i, j]

      # validity checks (NA-safe)
      if (!is.finite(nij) || !is.finite(wij)) {
        ll[s, d] <- NA_real_; next
      }
      li <- lam[i]; lj <- lam[j]
      if (!is.finite(li) || !is.finite(lj)) {
        ll[s, d] <- NA_real_; next
      }
      denom <- li + lj
      if (!is.finite(denom) || denom <= 0) {
        # undefined p when both zero or negative/NaN
        ll[s, d] <- NA_real_; next
      }

      p <- li / denom
      # ensure p is in [0,1]; tiny numeric noise guard
      if (!is.finite(p) || p < 0 || p > 1) {
        ll[s, d] <- NA_real_; next
      }

      ll[s, d] <- stats::dbinom(wij, size = nij, prob = p, log = TRUE)
    }
  }
  list(ll = ll, obs_idx = idx)
}

#' Log-likelihood matrix for the BT–SBM (clustered) model
#'
#' Builds an \eqn{S \times D} matrix using cluster labels \eqn{x_i}
#' and cluster rates \eqn{\lambda_k}. Assumes inputs are *relabelled consistently*
#' or that cluster ids in \code{x_samples[s, ]} are in \code{1..K} where
#' \code{K = ncol(lambda_samples)}.
#'
#' @param w_ij integer/numeric \eqn{n \times n} wins (i over j).
#' @param lambda_samples numeric \eqn{S \times K} matrix of cluster rates \eqn{\lambda_k}.
#' @param x_samples integer \eqn{S \times n} matrix of cluster labels for each item.
#'
#' @return A list with:
#' \itemize{
#' \item \code{ll} — \eqn{S \times D} matrix of log-likelihoods.
#' \item \code{obs_idx} — \eqn{D \times 2} matrix of (i,j) indices defining each column.
#' }
#' @importFrom stats dbinom
#' @export
make_bt_cluster_loo <- function(w_ij, lambda_samples, x_samples) {
  stopifnot(is.matrix(w_ij), is.matrix(lambda_samples), is.matrix(x_samples))
  n <- nrow(w_ij); stopifnot(n == ncol(w_ij))
  S <- nrow(x_samples)
  if (nrow(lambda_samples) != S) stop("Row count mismatch: nrow(lambda_samples) must equal nrow(x_samples) (S).")
  if (ncol(x_samples) != n) stop("x_samples must have n columns (items).")

  n_ij <- w_ij + t(w_ij)
  if (!isTRUE(all.equal(n_ij, t(n_ij)))) stop("n_ij must be symmetric.")
  if (any(diag(n_ij) != 0)) stop("Diagonal of n_ij must be zero.")

  idx <- which(upper.tri(n_ij) & n_ij > 0, arr.ind = TRUE)
  D <- nrow(idx)
  ll <- matrix(NA_real_, S, D)

  K <- ncol(lambda_samples)
  # sanity: labels must be in 1..K (allow NAs if present)
  if (any(x_samples < 1 | x_samples > K, na.rm = TRUE))
    stop("x_samples contains labels outside 1..K, with K = ncol(lambda_samples).")

  for (s in seq_len(S)) {
    lam <- lambda_samples[s, ]      # lambda_1..lambda_K for draw s
    xs  <- x_samples[s, ]           # labels in {1..K}
    lam_player <- lam[xs]           # length n
    for (d in seq_len(D)) {
      i <- idx[d, 1]; j <- idx[d, 2]
      nij <- n_ij[i, j]
      if (nij == 0) { ll[s, d] <- 0; next }
      wij <- w_ij[i, j]
      p <- lam_player[i] / (lam_player[i] + lam_player[j])
      ll[s, d] <- stats::dbinom(wij, size = nij, prob = p, log = TRUE)
    }
  }
  list(ll = ll, obs_idx = idx)
}

#' Compare BT models with Pareto-smoothed importance sampling LOO
#'
#' @param simple_llo list returned by \code{make_bt_simple_loo()}.
#' @param cluster_llo list returned by \code{make_bt_cluster_loo()}.
#' @return A list with:
#' \itemize{
#' \item \code{simple} — \code{loo} object for the simple BT.
#' \item \code{cluster} — \code{loo} object for the clustered BT–SBM.
#' \item \code{comparison} — result of \code{loo::compare_models()}.
#' }
#' @export
compare_bt_models_loo <- function(simple_llo, cluster_llo) {
  if (!requireNamespace("loo", quietly = TRUE))
    stop("Package 'loo' not installed. Install it to run LOO.")
  if (!is.list(simple_llo) || !is.list(cluster_llo) ||
      is.null(simple_llo$ll) || is.null(cluster_llo$ll))
    stop("Inputs must be lists with an 'll' matrix component.")

  loo_simple  <- loo::loo(simple_llo$ll)
  loo_cluster <- loo::loo(cluster_llo$ll)
  comp <- loo::compare_models(loo_simple, loo_cluster)
  list(simple = loo_simple, cluster = loo_cluster, comparison = comp)
}

#' Build BT clustered log-likelihood matrix (player-level)
#'
#' Convenience alternative returning a plain \eqn{S \times D} matrix
#' (\eqn{D} = #pairs with i<j). Here \code{x_draws} are labels and
#' \code{lambda_draws} are cluster rates per draw.
#'
#' @param x_draws integer \eqn{S \times n} matrix of cluster labels per draw.
#' @param lambda_draws numeric \eqn{S \times K} matrix of cluster rates per draw.
#' @param w_ij \eqn{n \times n} pairwise wins matrix (i over j, diag = 0).
#' @return Numeric matrix \eqn{S \times D} of log-likelihoods.
#' @importFrom stats dbinom
#' @keywords internal
make_bt_loo_cluster <- function(x_draws, lambda_draws, w_ij) {
  stopifnot(is.matrix(x_draws), is.matrix(lambda_draws), is.matrix(w_ij))
  n <- nrow(w_ij); stopifnot(n == ncol(w_ij))
  S <- nrow(x_draws)
  if (nrow(lambda_draws) != S) stop("Row mismatch: lambda_draws and x_draws must have same S.")
  n_ij <- w_ij + t(w_ij)

  pair_idx <- which(upper.tri(n_ij), arr.ind = TRUE)
  D <- nrow(pair_idx)
  K <- ncol(lambda_draws)

  # labels must not exceed K
  if (any(x_draws < 1 | x_draws > K, na.rm = TRUE))
    stop("x_draws contains labels outside 1..K, with K = ncol(lambda_draws).")

  log_lik <- matrix(0, nrow = S, ncol = D)

  for (s in seq_len(S)) {
    xs   <- x_draws[s, ]
    lam  <- lambda_draws[s, ]
    lam_player <- lam[xs]

    for (d in seq_len(D)) {
      i <- pair_idx[d, 1]; j <- pair_idx[d, 2]
      nij <- n_ij[i, j]
      if (nij == 0) {
        log_lik[s, d] <- 0
      } else {
        wij <- w_ij[i, j]
        pij <- lam_player[i] / (lam_player[i] + lam_player[j])
        log_lik[s, d] <- stats::dbinom(wij, size = nij, prob = pij, log = TRUE)
      }
    }
  }

  colnames(log_lik) <- apply(pair_idx, 1, paste, collapse = ",")
  log_lik
}

# ---------- Misc ----------

#' Shannon entropy (nats)
#' @param p numeric vector of nonnegative masses.
#' @return numeric(1) entropy in nats.
#' @keywords internal
shannon_entropy <- function(p) {
  probs <- p / sum(p)
  probs <- probs[probs > 0]  # remove zeros to avoid log(0)
  -sum(probs * log(probs))
}
#' Relabel partitions by decreasing \eqn{\lambda}, compute point estimates, and quantify uncertainty via credible balls
#'
#' Given posterior samples of labels \code{x_samples} (S x N) and corresponding
#' cluster-level intensities \code{lambda_samples} per iteration, this function:
#' (i) relabels each draw so that cluster \code{1} has the largest \eqn{\lambda},
#' cluster \code{2} the second largest, etc.; (ii) computes posterior similarity
#' matrix (PSM) and partition point estimates (minVI and Binder); and
#' (iii) constructs a VI-credible ball around the minVI partition and returns the
#' \emph{extremal} partitions on its surface (lower/upper), relabelled by decreasing
#' strength for interpretability. It also returns per-item assignment probabilities
#' and the posterior distribution of the number of clusters.
#'
#' @param x_samples Integer matrix (S x N) of sampled partitions; rows index MCMC
#'   iterations and columns index items. Cluster labels are arbitrary across iterations.
#' @param lambda_samples Either:
#'   \itemize{
#'     \item a \strong{list} of length S; element \code{[[s]]} is a numeric vector
#'           of cluster intensities \eqn{\lambda_\ell} indexed by the \emph{raw}
#'           label id used at iteration \code{s} (NAs allowed for non-occupied ids), or
#'     \item a \strong{matrix} (S x L) whose row \code{s} gives \eqn{\lambda_\ell}
#'           for label \eqn{\ell} at iteration \code{s} (sparse columns with NAs permitted).
#'   }
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{x_samples_relabel}}{Integer matrix (S x N) of relabelled draws
#'         (labels \code{1..K_s} per iteration \code{s}, ordered by decreasing \eqn{\lambda}).}
#'   \item{\code{lambda_samples_relabel}}{Numeric matrix (S x N) assigning each item its
#'         cluster's \eqn{\lambda} after relabelling.}
#'   \item{\code{co_clustering}}{Posterior similarity matrix (N x N).}
#'   \item{\code{minVI_partition}}{Partition estimated by minimizing posterior expected VI
#'         (first solution returned by \code{mcclust.ext::minVI}).}
#'   \item{\code{partition_binder}}{Partition estimated by Binder's loss
#'         (\code{mcclust.ext::minbinder.ext}).}
#'   \item{\code{credible_ball_lower_partition}}{Partition on the \emph{surface} of the
#'         95\% VI-credible ball that attains one extremum (relabelled by decreasing
#'         posterior-mean item strength).}
#'   \item{\code{credible_ball_upper_partition}}{Analogous extremal partition on the
#'         credible-ball surface (relabelled).}
#'   \item{\code{K_VI_lower}}{Number of clusters in \code{credible_ball_lower_partition}.}
#'   \item{\code{K_VI_upper}}{Number of clusters in \code{credible_ball_upper_partition}.}
#'   \item{\code{n_clusters_each_iter}}{Integer vector (length S) of occupied cluster counts per iteration.}
#'   \item{\code{block_count_distribution}}{Data frame with columns
#'         \code{num_blocks}, \code{count}, \code{prob} summarizing the posterior of the
#'         number of clusters.}
#'   \item{\code{item_cluster_assignment_probs}}{Data frame (N x Kmax) of per-item
#'         marginal assignment probabilities (columns \code{Cluster_1}, \code{Cluster_2}, ...).}
#'   \item{\code{avg_top_block_count}}{Average size of the top-\eqn{\lambda} cluster across iterations.}
#'   \item{\code{top_block_count_per_iter}}{Integer vector (length S) with the size
#'         of the top-\eqn{\lambda} cluster per iteration.}
#'   \item{\code{cluster_lambda_ordered}}{List of length S; each element is the vector of
#'         cluster \eqn{\lambda} values for that iteration, ordered decreasingly.}
#' }
#'
#' @details
#' \strong{Relabelling.} For each iteration, occupied labels are compacted to \code{1..K_s}
#' and reordered by decreasing \eqn{\lambda}, producing a canonical “1 = strongest” labelling.
#'
#' \strong{Point estimation.} The posterior similarity matrix is computed from relabelled
#' draws; minVI and Binder partitions are obtained via \pkg{mcclust.ext}.
#'
#' \strong{Credible ball and extremal partitions.} A 95\% credible ball in partition space
#' (under VI) is constructed around the minVI partition. We report the \emph{extreme}
#' partitions on the ball's surface (in the sense of maximal VI distance from the centre),
#' as returned by \code{mcclust.ext::credibleball}. These are then relabelled by decreasing
#' posterior-mean item strength to ensure a consistent “strength ordering” across summaries.
#' The associated cluster counts \code{K_VI_lower} and \code{K_VI_upper} characterize the
#' local structural uncertainty around the point estimate; they are \emph{not} marginal
#' posterior quantiles of \eqn{K}.
#'
#' @section Input requirements:
#' \itemize{
#'   \item \code{x_samples} must be integer-valued with no missing items per row.
#'   \item \code{lambda_samples} may be sparse (NAs for non-occupied labels).
#'   \item \pkg{mcclust} and \pkg{mcclust.ext} must be available; \code{credibleball}
#'         is expected to return lower/upper partitions in either \code{c.lower/c.upper}
#'         or \code{c.lowervert/c.uppervert}.
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' S <- 50; N <- 15
#' x_samps <- matrix(sample(1:4, S*N, TRUE), S, N)
#' # Sparse lambda per-iter: labels up to 6, only first 4 occupied typically
#' lam_list <- replicate(S, { v <- rep(NA_real_, 6); v[1:4] <- rexp(4, 1); v }, simplify = FALSE)
#' out <- relabel_by_lambda(x_samps, lam_list)
#' out$minVI_partition[1:10]
#' out$K_VI_lower; out$K_VI_upper
#' }
#'
#' @references
#' Wade, S., 2023. Bayesian cluster analysis. Philosophical Transactions of the Royal Society A: Mathematical, Physical and Engineering Sciences 381, 20220149. https://doi.org/10.1098/rsta.2022.0149
#'
#'
#' @seealso \code{\link[mcclust.ext]{minVI}}, \code{\link[mcclust.ext]{minbinder.ext}},
#'   \code{\link[mcclust.ext]{credibleball}}, \code{\link[mcclust]{comp.psm}}
#'
#' @import mcclust
#' @import mcclust.ext
#' @importFrom stats rexp
#' @export

relabel_by_lambda <- function(x_samples, lambda_samples) {
  stopifnot(is.matrix(x_samples))
  S <- nrow(x_samples); N <- ncol(x_samples)

  # Accept lambda as (i) list length S of numeric vectors indexed by raw label ID (sparse with NAs),
  # or (ii) matrix S x L (sparse columns with NAs).
  is_list_format <- is.list(lambda_samples)
  get_lambda_vec <- function(iter) {
    if (is_list_format) {
      v <- lambda_samples[[iter]]
      if (!is.numeric(v)) stop("lambda_samples[[iter]] must be numeric.")
      v
    } else {
      lambda_samples[iter, ]
    }
  }

  # Outputs
  x_relabeled              <- matrix(NA_integer_, S, N)
  lambda_per_item          <- matrix(NA_real_,    S, N)   # per-item lambda after relabel (S x N)
  cluster_lambda_ordered   <- vector("list", S)
  n_clusters_each_iter     <- integer(S)
  top_block_count_per_iter <- integer(S)

  # --- Per-iteration relabel by decreasing lambda ---
  for (iter in seq_len(S)) {
    xi <- as.integer(x_samples[iter, ])
    occ_raw <- sort(unique(xi))
    xi_seq  <- match(xi, occ_raw)           # 1..K
    K       <- max(xi_seq)

    lam_vec_full <- get_lambda_vec(iter)     # sparse per-label lambda; NA at non-occupied ids
    lam_occ <- rep(NA_real_, length(occ_raw))
    ok_idx  <- occ_raw <= length(lam_vec_full)
    lam_occ[ok_idx] <- lam_vec_full[occ_raw[ok_idx]]

    ord <- order(lam_occ, decreasing = TRUE, na.last = TRUE)
    occ_ord <- occ_raw[ord]
    lam_ord <- lam_occ[ord]
    if (anyNA(lam_ord)) lam_ord[is.na(lam_ord)] <- .Machine$double.xmin

    # map current labels -> rank by lambda (1 strongest)
    raw_to_ord_id <- integer(max(occ_ord))
    raw_to_ord_id[occ_ord] <- seq_len(K)
    xi_new <- raw_to_ord_id[ occ_raw[ xi_seq ] ]

    x_relabeled[iter, ]    <- xi_new
    lambda_per_item[iter, ] <- lam_ord[ xi_new ]
    cluster_lambda_ordered[[iter]] <- lam_ord
    n_clusters_each_iter[iter]     <- K
    top_block_count_per_iter[iter] <- sum(xi_new == 1L)
  }

  # modal K from relabelled draws
  modal_K <- as.integer(names(which.max(table(n_clusters_each_iter))))

  # --- Posterior similarity & point estimates (on relabelled x) ---
  psm <- mcclust::comp.psm(x_relabeled)
  partition_binder <- mcclust.ext::minbinder.ext(
    psm, cls.draw = x_relabeled, method = "avg", max.k = modal_K
  )$cl
  partition_minVI  <- mcclust.ext::minVI(
    psm, cls.draw = x_relabeled, method = "all", max.k = modal_K
  )$cl[1, ]

  # --- Credible ball (VI) around minVI ---
  x_ball <- mcclust.ext::credibleball(
    c.star = partition_minVI,
    cls.draw = x_relabeled,
    c.dist = "VI"
  )

  # Helper: relabel ANY partition by decreasing posterior-mean item strength
  # (gives a consistent ordering for minVI and credible-ball extrema)
  relabel_partition_by_item_mean_lambda <- function(z, lambda_item_mean) {
    stopifnot(length(z) == length(lambda_item_mean))
    z <- as.integer(z)
    labs <- sort(unique(z))
    # cluster mean strength from item-level posterior means
    cl_means <- vapply(labs, function(k) mean(lambda_item_mean[z == k], na.rm = TRUE), numeric(1))
    # order clusters by decreasing mean lambda, assign new ids 1..K
    ord <- order(cl_means, decreasing = TRUE)
    new_ids <- seq_along(labs)
    names(new_ids) <- labs[ord]
    z_new <- new_ids[as.character(z)]
    as.integer(z_new)
  }

  # Posterior-mean item strengths (S-averaged, after per-iter relabel)
  lambda_item_mean <- colMeans(lambda_per_item, na.rm = TRUE)

  # Extract lower/upper partitions from the credible ball.
  # Different mcclust.ext versions use c.lower/c.upper or c.lowervert/c.uppervert;
  # be robust to either.
  get_part <- function(obj, name1, name2) {
    if (!is.null(obj[[name1]])) obj[[name1]] else obj[[name2]]
  }
  c_lower_raw <- get_part(x_ball, "c.lower", "c.lowervert")
  c_upper_raw <- get_part(x_ball, "c.upper", "c.uppervert")

  if (is.null(c_lower_raw) || is.null(c_upper_raw)) {
    stop("credibleball() did not return lower/upper partitions in expected slots.")
  }

  # Relabel lower/upper by decreasing item-mean lambda
  c_lower_rl <- relabel_partition_by_item_mean_lambda(c_lower_raw, lambda_item_mean)
  c_upper_rl <- relabel_partition_by_item_mean_lambda(c_upper_raw, lambda_item_mean)

  K_VI_lower <- length(unique(c_lower_rl))
  K_VI_upper <- length(unique(c_upper_rl))

  # --- Assignment probabilities (up to Kmax = N, sparse beyond occupied levels) ---
  Kmax <- N
  assignment_probs <- matrix(0, nrow = N, ncol = Kmax)
  for (k in seq_len(Kmax)) {
    assignment_probs[, k] <- colMeans(x_relabeled == k, na.rm = TRUE)
  }
  colnames(assignment_probs) <- paste0("Cluster_", seq_len(Kmax))
  rownames(assignment_probs) <- paste0("Item_", seq_len(N))
  assignment_probs_df <- as.data.frame(assignment_probs)

  # --- Block-count distribution across iterations ---
  bc_tab <- table(n_clusters_each_iter)
  block_count_df <- data.frame(
    num_blocks = as.integer(names(bc_tab)),
    count      = as.vector(bc_tab),
    prob       = as.vector(bc_tab) / sum(bc_tab)
  )

  list(
    # relabelled draws and lambdas
    x_samples_relabel             = x_relabeled,             # S x N
    lambda_samples_relabel        = lambda_per_item,         # S x N
    cluster_lambda_ordered        = cluster_lambda_ordered,  # list length S

    # summaries
    co_clustering                 = psm,
    minVI_partition               = partition_minVI,
    partition_binder              = partition_binder,
    n_clusters_each_iter          = n_clusters_each_iter,
    block_count_distribution      = block_count_df,
    item_cluster_assignment_probs = assignment_probs_df,
    avg_top_block_count           = mean(top_block_count_per_iter),
    top_block_count_per_iter      = top_block_count_per_iter,

    # credible-ball extrema (RE-LABELLED by decreasing strength)
    credible_ball_lower_partition = c_lower_rl,
    credible_ball_upper_partition = c_upper_rl,
    K_VI_lower                    = K_VI_lower,
    K_VI_upper                    = K_VI_upper
  )
}

#' Map \eqn{\lambda} to Bradley–Terry \eqn{\theta = \lambda_i / (\lambda_i + \lambda_j)}
#'
#' @param lambda Numeric vector of positive rates \eqn{\lambda_i}.
#' @return A matrix with entries \eqn{\theta_{ij} = \lambda_i / (\lambda_i + \lambda_j)}.
#' Diagonal is 1/2 by the formula; you may overwrite if you prefer NA on the diagonal.
#' @examples
#' lambda_to_theta(c(1,2,4))
#' @export
lambda_to_theta <- function(lambda) {
  outer(lambda, lambda, function(a, b) a / (a + b))
}




#' Format a block-count posterior distribution as a pretty- table
#'
#' Given the output of relabel_by_lambda, this function produces a pretty table
#' with the posterior distribution of K, the number of blocks induced by the partition x.
#' It formats the result with `knitr::kable` and `kableExtra::kable_styling`.
#'
#' @param post The output of relabel_by_lambda()
#' @param p_min Numeric scalar. Probabilities \eqn{\le} this threshold are
#'   dropped before pivoting. Defaults to \code{0.05}.
#' @param format Character string. The format of the output table.
#' Default is "html". Switch to "latex" for a document-ready table.
#' @param digits Integer passed to \code{knitr::kable} for rounding. Defaults to \code{3}.
#' @param latex_options Character vector of LaTeX options passed to
#'   \code{kableExtra::kable_styling}. Defaults to \code{c("striped","hold_position")}.
#'
#' @return A \code{kableExtra} / \code{knitr_kable} object (suitable for LaTeX or HTML,
#'   depending on \code{format}).
#'
#' @details
#' The function produces a pretty table with the posterior distribution of K,
#' the number of blocks induced by the partition x.
#' it filters the clusters that have a p(K > \code{p_min}).
#'
#'
#' @importFrom dplyr select filter
#' @importFrom tidyr pivot_wider
#' @importFrom tidyselect all_of
#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling
#' @export
pretty_table_K_distribution <- function(
    post,
    format = 'html',
    p_min   = 0.05,
    digits     = 3,
    latex_options = c("striped", "hold_position")
) {
  # sanity checks

  df_wide <-
    post$block_count_distribution |>
    dplyr::select(
      num_blocks = num_blocks,
      prob       = prob
    ) |>
    dplyr::filter(prob > p_min) |>
    tidyr::pivot_wider(
      names_from  = num_blocks,
      values_from = prob
    )

  knitr::kable(df_wide, digits = digits,format = format) |>
    kableExtra::kable_styling(latex_options = latex_options)
}
