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
HGnedin <- function(n, k, gamma = 0.5) {
  exp(lchoose(n, k) + lgamma(k - gamma) - lgamma(1 - gamma) +
        log(gamma) + lgamma(n + gamma - k) - lgamma(n + gamma))
}

#' Urn weight: Dirichlet–Multinomial with finite cap K_DM
#' @param v_minus integer sizes of occupied clusters (excluding focal item).
#' @param beta_DM numeric concentration (>0).
#' @param K_DM integer maximum number of clusters (>=1).
#' @return Numeric vector length \eqn{K+1}: existing weights and new-cluster mass.
#' @keywords internal
urn_DM <- function(v_minus, beta_DM, K_DM) {
  K <- length(v_minus)
  c(v_minus + beta_DM, beta_DM * (K_DM - K) * (K_DM > K))
}

#' Urn weight: Dirichlet Process
#' @param v_minus integer sizes of occupied clusters (excluding focal item).
#' @param alpha_PY numeric concentration (>0).
#' @return Numeric vector length \eqn{K+1}: existing weights and new-cluster mass.
#' @keywords internal
urn_DP <- function(v_minus, alpha_PY) {
  c(v_minus, alpha_PY)
}

#' Urn weight: Gnedin heuristic
#' @param v_minus integer sizes of occupied clusters (excluding focal item).
#' @param gamma numeric parameter (>0).
#' @return Numeric vector length \eqn{K+1}.
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
    a = 0.01, b = 1,
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
  if (any(diag(w_ij) != 0)) stop("Diagonal of w_ij must be zero.")
  S <- nrow(lambda_samples)
  if (ncol(lambda_samples) != n) stop("lambda_samples must have n columns.")
  n_ij <- w_ij + t(w_ij)

  idx <- which(upper.tri(n_ij) & n_ij > 0, arr.ind = TRUE)
  D <- nrow(idx)
  ll <- matrix(NA_real_, S, D)

  for (s in seq_len(S)) {
    lam <- lambda_samples[s, ]
    for (d in seq_len(D)) {
      i <- idx[d, 1]; j <- idx[d, 2]
      nij <- n_ij[i, j]; wij <- w_ij[i, j]
      p <- lam[i] / (lam[i] + lam[j])
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

# ---------- Relabeling ----------

#' Relabel cluster assignments by decreasing lambda
#'
#' Given posterior samples of labels `x_samples` (S x N) and corresponding
#' cluster-level intensities `lambda_samples` per iteration, produce a
#' canonical relabeling where cluster 1 has the largest \eqn{\lambda},
#' cluster 2 the second largest, etc. Also computes co-clustering summaries
#' and assignment probabilities.
#'
#' @param x_samples Integer matrix S x N of sampled labels (arbitrary ids per iter).
#' @param lambda_samples Either:
#' \itemize{
#'   \item a **list** of length S with numeric vectors indexed by raw label id
#'         (NAs allowed for non-occupied ids), or
#'   \item a **matrix** S x L giving per-iteration \eqn{\lambda_\ell} for label \eqn{\ell}.
#' }
#'
#' @return A list with components:
#' \describe{
#'   \item{x_samples_relabel}{S x N integer matrix of relabeled draws (1..K per iter, ordered by \eqn{\lambda}).}
#'   \item{lambda_samples_relabel}{S x N numeric matrix assigning each item its cluster's \eqn{\lambda} after relabeling.}
#'   \item{item_cluster_assignment_probs}{N x Kmax data frame of marginal assignment probabilities.}
#'   \item{block_count_distribution}{Data frame of the distribution of the number of blocks across iterations.}
#'   \item{avg_top_block_count}{Average size of the top-\eqn{\lambda} block.}
#'   \item{co_clustering}{Posterior similarity matrix (N x N).}
#'   \item{minVI_partition}{Hard partition via minVI.}
#'   \item{partition_binder}{Hard partition via Binder's loss.}
#'   \item{n_clusters_each_iter}{Integer vector length S with number of blocks per iteration.}
#'   \item{top_block_count_per_iter}{Integer vector length S with size of the top block per iteration.}
#'   \item{cluster_lambda_ordered}{List length S of ordered \eqn{\lambda} vectors (length K per iter).}
#' }
#'
#' @details
#' The function first compacts raw labels to 1..K within each iteration,
#' then orders the occupied labels by decreasing \eqn{\lambda}, producing a
#' canonical labeling. Co-clustering summaries use \pkg{mcclust} and \pkg{mcclust.ext}.
#'
#' @examples
#' \dontrun{
#' S <- 100; N <- 20
#' set.seed(42)
#' x_samps <- matrix(sample(1:3, S*N, TRUE), S, N)
#' lam_list <- replicate(S, { v <- rep(NA_real_, 5); v[1:3] <- runif(3, 0.5, 2); v }, simplify=FALSE)
#' out <- relabel_by_lambda(x_samps, lam_list)
#' table(out$block_count_distribution$num_blocks)
#' }
#' @importFrom stats runif rexp
#' @import mcclust
#' @import mcclust.ext
#' @export
relabel_by_lambda <- function(x_samples, lambda_samples) {
  stopifnot(is.matrix(x_samples))
  S <- nrow(x_samples); N <- ncol(x_samples)

  # Accept (i) list length S of numeric vectors indexed by label ID (sparse with NAs),
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
  lambda_per_item          <- matrix(NA_real_,    S, N)     # per-item lambda after relabel (S x N)
  cluster_lambda_ordered   <- vector("list", S)             # per-iter numeric vector (length = K_iter)
  n_clusters_each_iter     <- integer(S)
  top_block_count_per_iter <- integer(S)

  for (iter in seq_len(S)) {
    xi <- as.integer(x_samples[iter, ])
    # compact to 1..K
    occ_raw <- sort(unique(xi))
    xi_seq <- match(xi, occ_raw)   # 1..K
    K <- max(xi_seq)

    lam_vec_full <- get_lambda_vec(iter)   # sparse per-label lambda with NAs at non-occupied ids

    # pull lambda for the occupied *raw* labels, keep order aligned with occ_raw
    lam_occ <- rep(NA_real_, length(occ_raw))
    ok_idx  <- occ_raw <= length(lam_vec_full)
    lam_occ[ok_idx] <- lam_vec_full[occ_raw[ok_idx]]

    # order clusters by decreasing lambda (NAs last)
    ord <- order(lam_occ, decreasing = TRUE, na.last = TRUE)
    occ_ord <- occ_raw[ord]
    lam_ord <- lam_occ[ord]
    if (anyNA(lam_ord)) lam_ord[is.na(lam_ord)] <- .Machine$double.xmin

    # map seq labels 1..K to new ids by lambda order: new ids = 1..K
    seq_to_raw <- occ_raw
    raw_to_ord_id <- integer(max(occ_ord))
    raw_to_ord_id[occ_ord] <- seq_len(K)
    xi_new <- raw_to_ord_id[ seq_to_raw[ xi_seq ] ]
    x_relabeled[iter, ] <- xi_new

    # per-item lambda after relabel
    lambda_per_item[iter, ] <- lam_ord[ xi_new ]

    # store per-iter ordered cluster lambda vector (length = K)
    cluster_lambda_ordered[[iter]] <- lam_ord

    n_clusters_each_iter[iter]     <- K
    top_block_count_per_iter[iter] <- sum(xi_new == 1L)
  }

  # Posterior similarity & hard partitions (on relabeled x)
  psm <- mcclust::comp.psm(x_relabeled)
  partition_binder <- mcclust.ext::minbinder.ext(psm, cls.draw = x_relabeled, method = "avg")$cl
  partition_minVI  <- mcclust.ext::minVI(psm, cls.draw = x_relabeled, method = "all")$cl[1, ]

  # Per-item assignment probabilities P(item i -> k), with k up to Kmax = N
  Kmax <- N
  assignment_probs <- matrix(0, nrow = N, ncol = Kmax)
  for (k in seq_len(Kmax)) {
    assignment_probs[, k] <- colMeans(x_relabeled == k, na.rm = TRUE)
  }
  colnames(assignment_probs) <- paste0("Cluster_", seq_len(Kmax))
  rownames(assignment_probs) <- paste0("Item_", seq_len(N))
  assignment_probs_df <- as.data.frame(assignment_probs)

  # Block-count distribution across iterations
  bc_tab <- table(n_clusters_each_iter)
  block_count_df <- data.frame(
    num_blocks = as.integer(names(bc_tab)),
    count      = as.vector(bc_tab),
    prob       = as.vector(bc_tab) / sum(bc_tab)
  )

  list(
    x_samples_relabel             = x_relabeled,             # S x N
    lambda_samples_relabel        = lambda_per_item,         # S x N
    item_cluster_assignment_probs = assignment_probs_df,     # N x Kmax
    block_count_distribution      = block_count_df,
    avg_top_block_count           = mean(top_block_count_per_iter),
    co_clustering                 = psm,
    minVI_partition               = partition_minVI,
    partition_binder              = partition_binder,
    n_clusters_each_iter          = n_clusters_each_iter,
    top_block_count_per_iter      = top_block_count_per_iter,
    cluster_lambda_ordered        = cluster_lambda_ordered
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
