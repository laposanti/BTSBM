#' Format player names as "Surname F."
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

#' Post-processing for BT–SBM posterior draws
#'
#' Relabels each MCMC draw by descending block rate \eqn{\lambda}, computes the
#' posterior similarity matrix (PSM), Binder and minVI point partitions, cluster
#' count summaries, and per-player assignment probabilities.
#'
#' @param x_samples integer matrix (M x K). Cluster labels per draw (each row = one draw).
#' @param lambda_samples numeric matrix (M x K). Block rates per draw (aligned with labels).
#'
#' @return A list with components:
#' \itemize{
#' \item \code{x_samples_relabel}, \code{lambda_samples_relabel} — relabeled draws.
#' \item \code{co_clustering} — PSM (K x K).
#' \item \code{partition_binder}, \code{minVI_partition}, \code{partition_expected}.
#' \item \code{n_clusters_each_iter}, \code{avg_n_clusters}.
#' \item \code{player_block_assignment_probs} (K x K data.frame).
#' \item \code{block_count_distribution} (data.frame with \code{num_blocks}, \code{count}, \code{prob}).
#' \item \code{top_block_count_per_iter}, \code{avg_top_block_count}.
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' M <- 100; K <- 6
#' x  <- matrix(sample.int(3, M*K, TRUE), M, K)
#' lam <- matrix(rexp(M*K), M, K)
#' out <- inference_helper(x, lam)
#' str(out)
#' }
#' @importFrom mcclust comp.psm
#' @importFrom mcclust.ext minbinder.ext minVI
#' @importFrom stats median
#' @export
#'
inference_helper <- function(x_samples, lambda_samples){
  stopifnot(is.matrix(x_samples), is.matrix(lambda_samples))
  n_iter <- nrow(x_samples)
  K <- ncol(x_samples)
  stopifnot(nrow(lambda_samples) == n_iter, ncol(lambda_samples) == K)

  new_x_samples <- matrix(NA_integer_, n_iter, K)
  new_lambda_samples <- matrix(NA_real_, n_iter, K)

  for (iter in seq_len(n_iter)) {
    occupied <- sort(unique(x_samples[iter, ]))
    current_lambda <- lambda_samples[iter, occupied]
    ord <- order(current_lambda, decreasing = TRUE)
    sorted_lambda <- current_lambda[ord]

    label_map <- rep(NA_integer_, K)
    for (r in seq_along(ord)) {
      old_label <- occupied[ord[r]]
      label_map[old_label] <- r
    }

    new_x_samples[iter, ] <- label_map[x_samples[iter, ]]
    for (i in seq_len(K)) {
      new_label_i <- new_x_samples[iter, i]
      new_lambda_samples[iter, i] <-
        if (!is.na(new_label_i) && new_label_i >= 1L && new_label_i <= length(sorted_lambda)) {
          sorted_lambda[new_label_i]
        } else {
          NA_real_
        }
    }
  }

  co_clust <- mcclust::comp.psm(new_x_samples)
  partition_binder <- mcclust.ext::minbinder.ext(psm = co_clust)$cl
  partition_minVI  <- mcclust.ext::minVI(psm = co_clust)$cl
  partition_expected <- apply(new_x_samples, 2, stats::median, na.rm = TRUE)

  lambda_means   <- colMeans(new_lambda_samples, na.rm = TRUE)
  lambda_medians <- apply(new_lambda_samples, 2, stats::median, na.rm = TRUE)

  n_clust_vec <- apply(new_x_samples, 1, function(z) length(unique(z[!is.na(z)])))
  avg_n_clust <- mean(n_clust_vec)

  assignment_probs <- matrix(0, nrow = K, ncol = K)
  for (i in seq_len(K)) {
    for (k in seq_len(K)) {
      assignment_probs[i, k] <- mean(new_x_samples[, i] == k, na.rm = TRUE)
    }
  }
  colnames(assignment_probs) <- paste0("Cluster_", seq_len(K))
  rownames(assignment_probs) <- paste0("Player_", seq_len(K))
  assignment_probs_df <- as.data.frame(assignment_probs)

  block_count_freq <- table(n_clust_vec)
  block_count_df <- data.frame(
    num_blocks = as.numeric(names(block_count_freq)),
    count = as.vector(block_count_freq),
    prob  = as.vector(block_count_freq) / sum(block_count_freq)
  )

  top_block_count_per_iter <- rowSums(new_x_samples == 1L, na.rm = TRUE)
  avg_top_block_count <- mean(top_block_count_per_iter)

  list(
    partition_binder              = partition_binder,
    partition_expected            = partition_expected,
    x_samples_relabel             = new_x_samples,
    lambda_samples_relabel        = new_lambda_samples,
    co_clustering                 = co_clust,
    minVI_partition               = partition_minVI,
    n_clusters_each_iter          = n_clust_vec,
    avg_n_clusters                = avg_n_clust,
    lambda_means                  = lambda_means,
    lambda_medians                = lambda_medians,
    player_block_assignment_probs = assignment_probs_df,
    block_count_distribution      = block_count_df,
    top_block_count_per_iter      = top_block_count_per_iter,
    avg_top_block_count           = avg_top_block_count
  )
}

#' Gnedin prior normalization weight H(V,h)
#'
#' Computes the unnormalized mass function term used in Gnedin-type priors.
#'
#' @param V integer(1) total items.
#' @param h integer vector of block counts.
#' @param gamma numeric(1) parameter \eqn{\gamma > 0}.
#' @return Numeric vector of weights \eqn{H(V,h)}.
#' @keywords internal
HGnedin <- function(V, h, gamma = 0.5){
  exp(lchoose(V, h) + lgamma(h - gamma) - lgamma(1 - gamma) +
        log(gamma) + lgamma(V + gamma - h) - lgamma(V + gamma))
}

#' Urn weight: Dirichlet–Multinomial with max H_DM clusters
#' @param v_minus integer sizes of occupied clusters (excluding current item).
#' @param beta_DM numeric concentration.
#' @param H_DM integer maximum number of clusters.
#' @return Numeric vector of length \eqn{H+1}: existing weights and new-cluster mass.
#' @keywords internal
urn_DM <- function(v_minus, beta_DM, H_DM){
  H <- length(v_minus)
  c(v_minus + beta_DM, beta_DM * (H_DM - H) * (H_DM > H))
}

#' Urn weight: Dirichlet Process
#' @param v_minus integer sizes of occupied clusters (excluding current item).
#' @param alpha_PY numeric concentration.
#' @return Numeric vector of length \eqn{H+1}: existing weights and new-cluster mass.
#' @keywords internal
urn_DP <- function(v_minus, alpha_PY){
  c(v_minus, alpha_PY)
}

#' Urn weight: Gnedin heuristic
#' @param v_minus integer sizes of occupied clusters (excluding current item).
#' @param gamma numeric parameter.
#' @return Numeric vector of length \eqn{H+1}.
#' @keywords internal
urn_GN <- function(v_minus, gamma) {
  H  <- length(v_minus)
  n_ <- sum(v_minus) + 1L               # +1 for the focal item
  c((v_minus + 1) * (n_ - H + gamma),
    H^2 - H * gamma)
}

#' Urn weight: Pitman–Yor process
#' @param v_minus integer sizes of occupied clusters (excluding current item).
#' @param alpha_PY numeric concentration.
#' @param sigma_PY numeric discount in [0,1).
#' @return Numeric vector of length \eqn{H+1}.
#' @keywords internal
urn_PY <- function(v_minus, alpha_PY, sigma_PY){
  H <- length(v_minus)
  c(v_minus - sigma_PY, alpha_PY + H * sigma_PY)
}

#' Expected number of clusters under finite/inf PY/DM-like priors (helper)
#'
#' @param n integer(1) number of items.
#' @param sigma numeric(1) discount (0 for DP/DM).
#' @param theta numeric(1) concentration parameter.
#' @param H integer(1) maximum number of clusters, or \code{Inf}.
#' @return Numeric(1) expected number of clusters.
#' @keywords internal
expected_cl_py <- function(n, sigma, theta, H){
  n <- as.integer(n)
  stopifnot(sigma >= 0, sigma < 1, theta > - sigma, n > 0, H > 1)
  if (is.infinite(H)) {
    if (sigma == 0) {
      out <- theta * sum(1 / (theta - 1 + 1:n))
    } else {
      out <- (1 / sigma) * exp(lgamma(theta + sigma + n) - lgamma(theta + sigma) -
                                 lgamma(theta + n) + lgamma(theta + 1)) - theta / sigma
    }
  } else {
    if (sigma == 0) {
      index <- 0:(n - 1L)
      out <- H - H * exp(sum(log(index + theta * (1 - 1 / H)) - log(theta + index)))
    } else {
      out <- NA_real_  # left undefined here
    }
  }
  out
}

#' Simple Bradley–Terry Gibbs sampler (no clustering)
#'
#' @param w_ij integer/numeric K x K wins from i over j.
#' @param n_ij integer/numeric K x K total matches (symmetric, diag = 0).
#' @param a,b numeric(1) Gamma(a,b) prior on each \eqn{\lambda_i}.
#' @param n_iter,burnin integers; total iterations and burn-in.
#' @param verbose logical; print progress.
#'
#' @return A list with \code{lambda_samples} (matrix of size (n_iter-burnin) x K).
#' @examples
#' \dontrun{
#' set.seed(1)
#' K <- 6
#' n <- matrix(0, K, K); n[upper.tri(n)] <- sample(0:3, sum(upper.tri(n)), TRUE)
#' n <- n + t(n); diag(n) <- 0
#' w <- matrix(0, K, K); w[upper.tri(w)] <- rbinom(sum(upper.tri(n)), n[upper.tri(n)], 0.5)
#' w <- w + t(n - w); diag(w) <- 0
#' fit <- gibbs_bt_simple(w, n, a = 1, b = 1, n_iter = 500, burnin = 100, verbose = FALSE)
#' }
#' @importFrom stats rgamma
#' @export
gibbs_bt_simple <- function(w_ij, n_ij,
                            a = 0.01, b = 1,
                            n_iter = 5000, burnin = 1000,
                            verbose = TRUE) {
  K <- nrow(w_ij)
  stopifnot(K == ncol(w_ij),
            K == nrow(n_ij), K == ncol(n_ij))

  w_i <- rowSums(w_ij)
  lambda <- stats::rgamma(K, a, b)
  Z      <- matrix(0, K, K)

  n_keep <- n_iter - burnin
  lambda_store <- matrix(NA_real_, n_keep, K)
  keep_idx <- 0L

  for (iter in seq_len(n_iter)) {
    for (i in 1:(K-1)) {
      for (j in (i+1):K) {
        nij <- n_ij[i, j]
        if (nij > 0) {
          rate  <- lambda[i] + lambda[j]
          Z_ij  <- stats::rgamma(1, nij, rate)
          Z[i,j] <- Z[j,i] <- Z_ij
        } else {
          Z[i,j] <- Z[j,i] <- 0
        }
      }
    }

    for (i in seq_len(K)) {
      shape <- a + w_i[i]
      rate  <- b + sum(Z[i, ])
      lambda[i] <- stats::rgamma(1, shape, rate)
    }

    if (iter > burnin) {
      keep_idx <- keep_idx + 1L
      lambda_store[keep_idx, ] <- lambda
    }

    if (verbose && iter %% 1000 == 0L)
      cat("simple BT | iter", iter, "\n")
  }

  invisible(list(lambda_samples = lambda_store))
}

#' Log-likelihood matrix for the *simple* Bradley–Terry model
#'
#' Builds an S x D matrix of log-likelihood values, where S is the number of
#' posterior draws and D is the number of observed unordered pairs (i<j) with
#' \code{n_ij > 0}.
#'
#' @param w_ij integer/numeric K x K wins (i over j).
#' @param n_ij integer/numeric K x K total matches (symmetric, diag=0).
#' @param lambda_samples numeric S x K matrix of player-specific rates \eqn{\lambda_i}.
#'
#' @return A list with:
#' \itemize{
#' \item \code{ll} — S x D matrix of log-likelihoods.
#' \item \code{obs_idx} — D x 2 matrix of (i,j) indices defining each column.
#' }
#' @importFrom stats dbinom
#' @export
make_bt_simple_loo <- function(w_ij, n_ij, lambda_samples) {
  stopifnot(is.matrix(w_ij), is.matrix(n_ij), is.matrix(lambda_samples))
  K <- nrow(w_ij)
  if (!identical(dim(w_ij), dim(n_ij))) stop("w_ij and n_ij must have same dims.")
  if (!isTRUE(all.equal(n_ij, t(n_ij)))) stop("n_ij must be symmetric.")
  if (any(diag(n_ij) != 0)) stop("Diagonal of n_ij must be zero.")
  S <- nrow(lambda_samples)
  if (ncol(lambda_samples) != K) stop("lambda_samples must have K columns.")

  idx <- which(upper.tri(n_ij) & n_ij > 0, arr.ind = TRUE)
  D <- nrow(idx)
  ll <- matrix(NA_real_, S, D)

  for (s in seq_len(S)) {
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
#' Builds an S x D matrix of log-likelihood values using cluster labels \eqn{x_i}
#' and cluster rates \eqn{\lambda_k}. Assumes inputs are *relabelled consistently*.
#'
#' @param w_ij integer/numeric K x K wins (i over j).
#' @param n_ij integer/numeric K x K total matches (symmetric, diag=0).
#' @param lambda_samples numeric S x K matrix of cluster rates \eqn{\lambda_k}.
#' @param x_samples integer S x K matrix of cluster labels for each player.
#'
#' @return A list with:
#' \itemize{
#' \item \code{ll} — S x D matrix of log-likelihoods.
#' \item \code{obs_idx} — D x 2 matrix of (i,j) indices defining each column.
#' }
#' @importFrom stats dbinom
#' @export
make_bt_cluster_loo <- function(w_ij, n_ij, lambda_samples, x_samples) {
  stopifnot(is.matrix(w_ij), is.matrix(n_ij),
            is.matrix(lambda_samples), is.matrix(x_samples))
  K <- nrow(w_ij)
  if (!identical(dim(w_ij), dim(n_ij))) stop("w_ij and n_ij must have same dims.")
  if (!isTRUE(all.equal(n_ij, t(n_ij)))) stop("n_ij must be symmetric.")
  if (any(diag(n_ij) != 0)) stop("Diagonal of n_ij must be zero.")
  if (!identical(dim(lambda_samples), dim(x_samples)))
    stop("lambda_samples and x_samples must have identical dims (S x K).")

  S <- nrow(x_samples)
  idx <- which(upper.tri(n_ij) & n_ij > 0, arr.ind = TRUE)
  D <- nrow(idx)
  ll <- matrix(NA_real_, S, D)

  for (s in seq_len(S)) {
    lam <- lambda_samples[s, ]  # λ_1..λ_K
    xs  <- x_samples[s, ]       # labels in {1..K}
    lam_player <- lam[xs]
    for (d in seq_len(D)) {
      i <- idx[d, 1]; j <- idx[d, 2]
      n <- n_ij[i, j]
      if (n == 0) { ll[s, d] <- 0; next }
      w <- w_ij[i, j]
      p <- lam_player[i] / (lam_player[i] + lam_player[j])
      ll[s, d] <- stats::dbinom(w, size = n, prob = p, log = TRUE)
    }
  }
  list(ll = ll, obs_idx = idx)
}

#' Compare BT models with PSIS-LOO
#'
#' @param simple_llo list returned by \code{make_bt_simple_loo()}.
#' @param cluster_llo list returned by \code{make_bt_cluster_loo()}.
#' @return A list with:
#' \itemize{
#' \item \code{simple} — \code{loo} object for the simple BT.
#' \item \code{cluster} — \code{loo} object for the clustered BT–SBM.
#' \item \code{comparison} — result of \code{loo::compare()}.
#' }
#' @export
#'
#'
#'
compare_bt_models_loo <- function(simple_llo, cluster_llo) {
  if (!requireNamespace("loo", quietly = TRUE))
    stop("Package 'loo' not installed. Install it to run LOO.")
  if (!is.list(simple_llo) || !is.list(cluster_llo) ||
      is.null(simple_llo$ll) || is.null(cluster_llo$ll))
    stop("Inputs must be lists with an 'll' matrix component.")

  loo_simple  <- loo::loo(simple_llo$ll)
  loo_cluster <- loo::loo(cluster_llo$ll)
  comp <- loo::compare(loo_simple, loo_cluster)
  list(simple = loo_simple, cluster = loo_cluster, comparison = comp)
}

#' Build BT clustered log-likelihood matrix (player-level)
#'
#' Convenience alternative returning a plain S x P matrix (P = #pairs).
#' @param x_draws integer S x K matrix of cluster labels per draw.
#' @param lambda_draws numeric S x K matrix of cluster rates per draw.
#' @param w_ij, n_ij K x K wins and matches.
#' @return Numeric matrix S x P of log-likelihoods.
#' @importFrom stats dbinom
#' @keywords internal
make_bt_loo_cluster <- function(x_draws,
                                lambda_draws,
                                w_ij,
                                n_ij) {
  stopifnot(dim(x_draws) == dim(lambda_draws),
            nrow(w_ij)  == ncol(w_ij),
            nrow(w_ij)  == nrow(n_ij),
            ncol(w_ij)  == ncol(n_ij))

  M <- nrow(x_draws)
  K <- ncol(x_draws)

  pair_idx <- which(upper.tri(n_ij), arr.ind = TRUE)
  P <- nrow(pair_idx)

  log_lik <- matrix(0, nrow = M, ncol = P)

  for (m in seq_len(M)) {
    x_m   <- x_draws[m, ]
    lam_m <- lambda_draws[m, ]
    lambda_player <- lam_m[x_m]

    for (p in seq_len(P)) {
      i <- pair_idx[p, 1]
      j <- pair_idx[p, 2]
      nij <- n_ij[i, j]
      if (nij == 0) {
        log_lik[m, p] <- 0
      } else {
        wij <- w_ij[i, j]
        pij <- lambda_player[i] / (lambda_player[i] + lambda_player[j])
        log_lik[m, p] <- stats::dbinom(wij, size = nij, prob = pij, log = TRUE)
      }
    }
  }

  colnames(log_lik) <- apply(pair_idx, 1, paste, collapse = ",")
  log_lik
}

#' Shannon entropy (nats)
#' @param p numeric vector of nonnegative masses.
#' @return numeric(1) entropy in nats.
#' @keywords internal
shannon_entropy <- function(p) {
  probs <- p / sum(p)
  probs <- probs[probs > 0]  # remove zeros to avoid log(0)
  -sum(probs * log(probs))
}
