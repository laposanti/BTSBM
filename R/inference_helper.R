#' Post-processing for BT–SBM posterior draws
#'
#' Relabels each MCMC draw by descending block rate \eqn{\lambda}, computes the
#' posterior similarity matrix (PSM), point partitions (Binder and minVI),
#' summary stats on the number of clusters, and per-player assignment probabilities.
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
#' @details
#' Relabeling is done *within* each draw by sorting occupied blocks by
#' decreasing \eqn{\lambda}. Block 1 is the “top” block in that draw.
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
#' @export
#'
inference_helper <- function(x_samples, lambda_samples) {
  stopifnot(is.matrix(x_samples), is.matrix(lambda_samples))
  M <- nrow(x_samples); K <- ncol(x_samples)
  if (!identical(dim(x_samples), dim(lambda_samples)))
    stop("x_samples and lambda_samples must have the same dimensions (M x K).")

  # ---------- 1) Relabel by descending lambda within each draw ----------
  new_x  <- matrix(NA_integer_, M, K)
  new_lam <- matrix(NA_real_, M, K)

  for (m in seq_len(M)) {
    occ <- sort(unique(x_samples[m, ]))
    lam_occ <- lambda_samples[m, occ]
    ord <- order(lam_occ, decreasing = TRUE)
    sorted_lambda <- lam_occ[ord]

    # map old -> new labels (1..H by decreasing lambda)
    label_map <- rep(NA_integer_, K)
    for (r in seq_along(ord)) {
      old_label <- occ[ord[r]]
      label_map[old_label] <- r
    }

    new_x[m, ] <- label_map[x_samples[m, ]]
    # assign lambda to players via relabeled cluster id
    for (i in seq_len(K)) {
      ki <- new_x[m, i]
      new_lam[m, i] <- if (is.na(ki)) NA_real_ else sorted_lambda[ki]
    }
  }

  # ---------- 2) PSM and point partitions ----------
  psm <- mcclust::comp.psm(new_x)
  binder <- mcclust.ext::minbinder.ext(psm = psm)$cl
  minvi  <- mcclust.ext::minVI(psm = psm)$cl
  part_expected <- apply(new_x, 2, stats::median, na.rm = TRUE)

  # ---------- 3) Cluster counts ----------
  n_clusters_each_iter <- apply(new_x, 1, function(z) length(unique(z[!is.na(z)])))
  avg_n_clusters <- mean(n_clusters_each_iter)

  # posterior over #blocks
  bc_tab <- table(n_clusters_each_iter)
  block_count_df <- data.frame(
    num_blocks = as.numeric(names(bc_tab)),
    count      = as.vector(bc_tab),
    prob       = as.vector(bc_tab) / sum(bc_tab)
  )

  # ---------- 4) Assignment probabilities P(player i in cluster k) ----------
  assign_probs <- matrix(0, nrow = K, ncol = K)
  for (i in seq_len(K)) {
    for (k in seq_len(K)) assign_probs[i, k] <- mean(new_x[, i] == k, na.rm = TRUE)
  }
  colnames(assign_probs) <- paste0("Cluster_", seq_len(K))
  rownames(assign_probs) <- paste0("Player_", seq_len(K))
  assign_probs_df <- as.data.frame(assign_probs)

  # ---------- 5) Size of the top block (label 1) ----------
  top_block_count_per_iter <- rowSums(new_x == 1L, na.rm = TRUE)
  avg_top_block_count <- mean(top_block_count_per_iter)

  list(
    x_samples_relabel             = new_x,
    lambda_samples_relabel        = new_lam,
    co_clustering                 = psm,
    partition_binder              = binder,
    minVI_partition               = minvi,
    partition_expected            = part_expected,
    n_clusters_each_iter          = n_clusters_each_iter,
    avg_n_clusters                = avg_n_clusters,
    player_block_assignment_probs = assign_probs_df,
    block_count_distribution      = block_count_df,
    top_block_count_per_iter      = top_block_count_per_iter,
    avg_top_block_count           = avg_top_block_count
  )
}

#' Shannon entropy of a discrete distribution
#' @param p numeric vector of nonnegative masses.
#' @return numeric(1) entropy in nats.
#' @keywords internal
shannon_entropy <- function(p) {
  probs <- p / sum(p)
  probs <- probs[probs > 0]
  -sum(probs * log(probs))
}
