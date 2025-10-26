#' Block-ordered adjacency heatmap (BT-SBM)
#'
#' Plot the adjacency (wins / matches) heatmap with rows/cols ordered by
#' the hard partition (\code{fit$estimates$x_hat}) and, within blocks,
#' by marginal wins. Draws block boundary lines and (optionally) a side color strip
#' for block IDs if \pkg{ggside} is installed.
#'
#' @param fit Output list from \code{gibbs_BT_SBM()}.
#' @param w_ij Integer matrix of wins (same players & order used in \code{fit}).
#' @param x_hat partition point estimate. One n-length integer vector.
#' @param clean_fun Optional function to prettify player names. Default: identity.
#' @param palette Named colors for blocks (as character vector). Defaults to a
#'   Wimbledon-ish palette.
#' @param fill_low,fill_high Colors for the heatmap pixels gradient low/high.
#' @return A \code{ggplot} object.
#' @examples
#' \dontrun{
#' p <- plot_block_adjacency(fit, w_ij)
#' }
#' @export
#'
plot_block_adjacency <- function(
    fit,
    w_ij,
    x_hat = NULL,
    clean_fun = clean_players_names,
    palette = NULL,
    fill_low = "#FFFFCC",
    fill_high = "#006400"
){
  stopifnot(is.matrix(w_ij))
  N_ij = w_ij + t(w_ij)
  if(is.null(x_hat))  x_hat <- fit$minVI_partition
  stopifnot(length(x_hat) == nrow(w_ij))
  # Input names
  pl_names <- rownames(w_ij)

  make_palette <- function(
    n,
    low  = "#00441B",                       # was "#FFFFFF"
    mids = c("#CDEB8B", "#78AB46", "#FFD700"),
    high = "#FF8C00",
    bias = 1
  ) {
    grDevices::colorRampPalette(c(low, mids, high), bias = bias)(n)
  }

  K= length(unique(x_hat))
  if(is.null(palette)) palette <- make_palette(K)
  stopifnot(length(palette) == K)
  if (is.null(pl_names)) pl_names <- paste0("Item_", seq_len(nrow(w_ij)))
  colnames(w_ij) <- pl_names
  rownames(w_ij) <- pl_names


  # Add cluster info to players
  df_cl <- data.frame(
    players = pl_names,
    cl = x_hat,
    marginal_win = rowSums(w_ij)
  )

  # Summarize block compositions
  block_players <- aggregate(players ~ cl, data = df_cl, FUN = \(x) paste(x, collapse = ", "))
  block_players$Block_Size <- sapply(strsplit(block_players$players, ", "), length)

  # Prepare long format for adjacency matrix
  Y_long <- reshape2::melt(w_ij)
  colnames(Y_long) <- c("Winner", "Loser", "Win_Count")
  Y_long$Matches_Count <- reshape2::melt(N_ij)$value

  K <- length(unique(x_hat))

  Y_long_plot <- Y_long %>%
    mutate(perc_success = Win_Count / Matches_Count) %>%
    left_join(df_cl, by = c("Loser" = "players")) %>%
    rename(row_cl = cl, marginal_win_row = marginal_win) %>%
    left_join(df_cl, by = c("Winner" = "players")) %>%
    rename(col_cl = cl, marginal_win_col = marginal_win) %>%
    mutate(
      Winner = sapply(stringr::str_to_title(gsub("_", " ", Winner)), clean_fun),
      Loser  = sapply(stringr::str_to_title(gsub("_", " ", Loser)),  clean_fun)
    )

  # Reorder factors by cluster and marginal wins
  Y_long_plot <- Y_long_plot %>%
    mutate(
      Winner = factor(Winner, levels = unique(Winner[order(col_cl, -marginal_win_col,decreasing = T)])),
      Loser  = factor(Loser,  levels = unique(Loser[order(row_cl, -marginal_win_row)])),
      col_cl = factor(col_cl, ordered = TRUE)
    )
  # 1. Get block boundaries for the x-axis (i.e. losers) using row_cl:
  v_lines_list <- Y_long_plot %>%
    group_by(row_cl) %>%
    summarize(x_break = max(as.numeric(Loser)), .groups = "drop") %>%
    pull(x_break)

  # Remove the last boundary so we don't draw a line at the extreme edge:
  v_lines_list <- v_lines_list[-length(v_lines_list)]

  # 2. Get block boundaries for the y-axis (i.e. winners) using col_cl:
  h_lines_list <- Y_long_plot %>%
    group_by(col_cl) %>%
    summarize(y_break = min(as.numeric(Winner)), .groups = "drop") %>%
    pull(y_break)

  # Remove the last boundary similarly:
  h_lines_list <- h_lines_list[-length(h_lines_list)]

  geom_adjacency_fixed <- ggplot(Y_long_plot, aes(x = Loser, y = Winner)) +
    geom_tile(aes(fill = perc_success), color = 'grey40') +
    scale_fill_gradient(low = "#FFFFCC", high = "#006400", na.value = "#009680") +
    ggside::geom_ysidecol(aes(color = factor(col_cl))) +
    scale_color_manual(values = palette) +
    geom_vline(xintercept = unlist(v_lines_list) + 0.5, color = 'black', linewidth = 0.3) +
    geom_hline(yintercept = unlist(h_lines_list) - 0.5, color = 'black', linewidth = 0.3) +
    labs(
      x = "Items (ordered by block)",
      y = "Items (ordered by block)",
      fill = "% success",
      color = "Block"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.y = element_text(size = 7),
      axis.text.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "left"
    ) +
    ggside::theme_ggside_void() +
    scale_y_discrete(guide = guide_axis(n.dodge = 2)) +
    coord_fixed(ratio = 1)

  print(geom_adjacency_fixed)
}


#' Assignment probability heatmap (BT-SBM)
#'
#' Heatmap of item-wise posterior assignment probabilities for clusters
#' (relabeled so that Cluster 1 is the top block by decreasing \eqn{\lambda}).
#' Items are ordered by their most probable cluster and marginal wins.
#'
#' @param fit Output list from \code{gibbs_BT_SBM()} (must include \code{relabeled$assign_prob}).
#' @param w_ij Optional wins matrix to compute marginal wins for ordering and annotation.
#'   If \code{NULL}, items are ordered by most-probable cluster only.
#' @param clean_fun Optional function to prettify names. Default: identity.
#' @param k_show Optional integer number of clusters to show (defaults to all columns in \code{assign_prob}).
#' @param max_n_clust Where to filter the mcmc x_t. If not specified we use the modal K
#' @param fill_low,fill_high Colors for the heatmap gradient low/high.
#' @return A \code{ggplot} object.
#' @examples
#' \dontrun{
#' p <- plot_assignment_probabilities(fit, w_ij)
#' }
#' @export
plot_assignment_probabilities <- function(
    fit,
    w_ij = NULL,
    max_n_clust= NULL,
    clean_fun = clean_players_names,
    k_show = NULL,
    fill_low = "#FFFFCC",
    fill_high = "#006400"
){


  item_names <- rownames(w_ij)
  if (is.null(item_names)) item_names <- paste0("Item_", seq_len(nrow(w_ij)))


  count_cl = function(x){
    length(unique(x))
  }

  x_samples       <- fit$x_samples
  lambda_samples  <- fit$lambda_samples

  unique_count <- apply(x_samples, 1, function(z) length(unique(z)))
  modal_K = as.numeric(names(which.max(table(unique_count))))
  if(is.null(max_n_clust)){max_n_clust <- modal_K}

  x_samples_sub <- x_samples[which(unique_count == max_n_clust), ]
  lambda_samples_sub <- lambda_samples[which(unique_count == max_n_clust)]

  block_prob <- fit$item_cluster_assignment_probs[, 1:max_n_clust]
  block_prob <- as.data.frame(block_prob)
  block_prob <- cbind(block_prob, pl_name = item_names)  # <-- boom if nrow(block_prob)==0


  assignment_probs_long <- block_prob %>%
    pivot_longer(cols = 1:max_n_clust, names_to = "Cluster", values_to = "prob") %>%
    mutate(
      Cluster = gsub(x = Cluster, replacement = " ", pattern = "_"),
      pl_name = gsub(x = pl_name, replacement = " ", pattern = "_")
    )

  # Cluster with highest assignment probability
  max_prob_clusters <- assignment_probs_long %>%
    group_by(pl_name) %>%
    summarize(Cl_ass = Cluster[which.max(prob)], .groups = "drop")

  # Marginal win percentages
  marg_pro_win <- data.frame(
    pl_name = item_names,
    marg_pro_win  = rowSums(w_ij),
    marg_pro_loss = colSums(w_ij)
  ) %>%
    mutate(
      pct_win = marg_pro_win / (marg_pro_loss + marg_pro_win),
      pl_name = gsub(x = pl_name, replacement = " ", pattern = "_")
    )

  assignment_probs_long_plot <- assignment_probs_long %>%
    left_join(max_prob_clusters, by = "pl_name") %>%
    left_join(marg_pro_win, by = "pl_name") %>%
    ungroup()|>
    mutate(
      pl_name = sapply(stringr::str_to_title(gsub("_", " ", pl_name)), clean_fun),
    )|>
    mutate(pl_name = factor(pl_name,
                            levels = unique(pl_name[order(Cl_ass, -marg_pro_win,decreasing = T)]),
                            ordered = TRUE))

  ass_prob_plot <- ggplot(assignment_probs_long_plot) +
    geom_tile(aes(x = Cluster, y = pl_name, fill = prob)) +
    scale_fill_gradient(low = "#FFFFCC", high = "#006400", na.value = "#009680") +
    labs(x = "", y = "", fill = "Assign. Prob.") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    scale_y_discrete(guide = guide_axis(n.dodge = 2))

  ass_prob_plot
}


#' Lambda uncertainty plot (per player)
#'
#' Forest plot of per-player posterior \eqn{\lambda} with uncertainty intervals,
#' using relabeled draws. Points colored by the hard partition (\code{fit$estimates$x_hat}).
#'
#' @param fit Output from \code{gibbs_BT_SBM()} with \code{opt_lambda$lambda_item} computed
#'   (set \code{keep_lambda=TRUE} when sampling).
#' @param w_ij Optional wins matrix to compute marginal wins for ordering and annotation.
#' @param log_base Base for the x-axis logarithm (10 or e). Defaults to 10.
#' @param prob Interval probability for HPD (e.g., 0.90).
#' @param max_n_clust Where to filter the mcmc x_t. If not specified we use the modal K
#' @param palette Named colors for clusters.
#' @param clean_fun Optional player-name cleaner. Default: identity.
#' @return A \code{ggplot} object.
#' @examples
#' \dontrun{
#' # fit with keep_lambda=TRUE
#' p <- plot_lambda_uncertainty(fit, prob = 0.90)
#' }
#' @export
#'
plot_lambda_uncertainty <- function(
    fit,
    w_ij,
    log_base = 2.718,
    max_n_clust = NULL,
    prob = 0.90,
    palette = NULL,
    clean_fun = clean_players_names,
    x_hat = NULL
){

  # Expected shapes:
  #   lambda_samples_relabel: S times N   (per-iteration, per-item lambda)
  #   x_samples_relabel      : S times N   (per-iteration, per-item cluster labels)
  #   partition_binder       : N       (final hard clustering per item)
  lambda_item <- fit$lambda_samples_relabel
  x_samples   <- fit$x_samples_relabel
  if(is.null(x_hat)) x_hat = fit$minVI_partition
  cluster_hat <- x_hat
  make_palette <- function(
    n,
    low  = "#00441B",                       # was "#FFFFFF"
    mids = c("#CDEB8B", "#78AB46", "#FFD700"),
    high = "#FF8C00",
    bias = 1
  ) {
    grDevices::colorRampPalette(c(low, mids, high), bias = bias)(n)
  }

  K= length(unique(x_hat))
  if(is.null(palette)) palette <- make_palette(K)

  if (is.null(lambda_item) || is.null(x_samples))
    stop("Expected fit$lambda_samples_relabel and fit$x_samples_relabel to be present.")

  # Item names ---------------------------------------------------------------
  item_names <- colnames(w_ij)
  if (is.null(item_names)) item_names <- paste0("Item_", seq_len(ncol(lambda_item)))
  colnames(lambda_item) <- item_names

  # Optionally restrict to iterations with exactly max_n_clust clusters -------
  if (!is.null(x_samples) && nrow(x_samples) == nrow(lambda_item)) {
    unique_count <- apply(x_samples, 1, function(z) length(unique(z)))
    modal_K = as.numeric(names(which.max(table(unique_count))))
    if(is.null(max_n_clust)){max_n_clust <- modal_K}
    keep_idx <- which(unique_count == max_n_clust)
    if (length(keep_idx) > 0L) {
      lambda_player <- lambda_item[keep_idx, , drop = FALSE]
    } else {
      # fallback: use all iterations if no iteration matches the requested count
      lambda_player <- lambda_item
      warning("No iterations with exactly max_n_clust; using all iterations.")
    }
  } else {
    lambda_player <- lambda_item
  }

  # Numerical guard for log ---------------------------------------------------
  eps <- .Machine$double.eps
  if (any(lambda_player <= 0, na.rm = TRUE)) {
    warning("Non-positive lambda encountered; adding a small epsilon before log.")
  }

  # Log helpers honoring the chosen base
  logb <- function(x) log(pmax(x, eps), base = log_base)
  powb <- function(x) log_base^x

  # HPD on log scale, then back-transform
  log_lp <- apply(lambda_player, 2, logb)            # S' x N  (or matrix if S'>1)
  # Ensure matrix
  if (is.vector(log_lp)) log_lp <- matrix(log_lp, ncol = 1)

  hpd_mat <- t(apply(log_lp, 2, function(v)
    coda::HPDinterval(coda::as.mcmc(v), prob = prob)[1, ]))

  player_summ <- data.frame(
    player  = item_names,
    mean    = powb(colMeans(log_lp)),    # geometric mean
    low     = powb(hpd_mat[, 1]),
    up      = powb(hpd_mat[, 2]),
    cluster = as.integer(cluster_hat),
    row.names = NULL
  )

  # Plot ----------------------------------------------------------------------
  # colour mapping limited to present clusters to avoid unused palette entries
  present_clusters <- sort(unique(player_summ$cluster))

  plot_lambda <- player_summ |>
    dplyr::arrange(cluster, dplyr::desc(mean)) |>
    dplyr::mutate(
      player = sapply(stringr::str_to_title(gsub("_", " ", player)), clean_fun),
      player = factor(player, levels = rev(player)),
      surname = gsub(".*[ _]", "", player)
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = mean, y = player, colour = factor(cluster))) +
    ggplot2::geom_pointrange(ggplot2::aes(xmin = low, xmax = up), size = 0.4, fatten = 0.6) +
    scale_x_continuous(
      trans  = scales::log_trans(base = log_base),
      labels = scales::label_number(accuracy = 0.01),  # 2 decimals
      breaks = scales::breaks_log(base = log_base)     # optional: nice log breaks
    )+
    ggplot2::scale_colour_manual(values = palette) +
    ggplot2::scale_y_discrete(guide = ggplot2::guide_axis(n.dodge = 2)) +
    ggplot2::labs(
      x = bquote(lambda~"(posterior mean)- x-axis in log scale"),
      y = NULL,
      colour = "Cluster"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 6),
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank()
    )

  return(plot_lambda)
}

