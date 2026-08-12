# Beginner-oriented plots for the common BTSBM fit interface.
#
# The older plotting functions remain available for the legacy sampler output.
# These functions accept a btsbm_fit object directly and avoid exposing
# internal MCMC labels as though they were stable scientific quantities.

.require_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
  }
}

#' Plot posterior item strengths and their uncertainty
#'
#' This forest plot answers a first analysis question: which items appear
#' stronger, and how uncertain is that ordering? Overlapping intervals mean
#' that the data do not cleanly distinguish the corresponding strengths.
#'
#' @param fit A simple BT/PL or item-SBM btsbm_fit object.
#' @param credible_mass Width of the central posterior interval, in (0, 1).
#'
#' @return A ggplot object.
#' @export
plot_strength_summary <- function(fit, credible_mass = 0.9) {
  .require_ggplot2()
  summary_data <- strength_summary(fit, credible_mass = credible_mass)
  summary_data$item <- factor(summary_data$item, levels = rev(summary_data$item))
  likelihood_label <- if (fit$model$likelihood == "pl") {
    "Relative item strength (geometric mean = 1)"
  } else {
    "Relative item strength"
  }

  ggplot2::ggplot(summary_data, ggplot2::aes(x = .data$mean, y = .data$item)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
      width = 0, linewidth = 0.5, colour = "#4D4D4D"
    ) +
    ggplot2::geom_point(size = 2.2, colour = "#0072B2") +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "#7F7F7F") +
    ggplot2::labs(
      x = likelihood_label,
      y = NULL,
      title = "Posterior item strengths",
      subtitle = paste0(round(credible_mass * 100), "% credible intervals"),
      caption = "A larger strength indicates a more preferred item (PL) or a stronger competitor (BT)."
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

#' Plot the posterior similarity of items or rankings
#'
#' A dark cell means that two entities were often assigned to the same latent
#' group across posterior draws. The plot is label-invariant, which makes it
#' safer to interpret than the raw MCMC cluster labels.
#'
#' @param fit A btsbm_fit object containing an item or ranking partition.
#' @param target Whether to plot latent groups of "item" or "ranking".
#' @param order_rows Whether to cluster the rows and columns for display.
#'
#' @return A ggplot object.
#' @export
plot_posterior_similarity <- function(fit, target = c("item", "ranking"),
                                      order_rows = TRUE) {
  .require_ggplot2()
  if (!inherits(fit, "btsbm_fit")) {
    stop("'fit' must be a btsbm_fit object.", call. = FALSE)
  }
  target <- match.arg(target)
  similarity <- posterior_similarity(fit, target = target)
  labels <- if (target == "item") fit$data$item_labels else fit$data$ranking_labels
  display_order <- seq_along(labels)
  if (isTRUE(order_rows) && length(labels) > 1L) {
    dissimilarity_matrix <- 1 - similarity
    dissimilarity_matrix[dissimilarity_matrix < 0] <- 0
    diag(dissimilarity_matrix) <- 0
    dissimilarity <- stats::as.dist(dissimilarity_matrix)
    display_order <- stats::hclust(dissimilarity, method = "average")$order
  }
  similarity <- similarity[display_order, display_order, drop = FALSE]
  labels <- labels[display_order]
  plot_data <- expand.grid(
    row = seq_along(labels), column = seq_along(labels),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  plot_data$probability <- as.vector(similarity)
  plot_data$row_label <- factor(labels[plot_data$row], levels = rev(labels))
  plot_data$column_label <- factor(labels[plot_data$column], levels = labels)

  ggplot2::ggplot(
    plot_data, ggplot2::aes(x = .data$column_label, y = .data$row_label)
  ) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$probability), colour = "white", linewidth = 0.15) +
    ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#08519C", limits = c(0, 1)) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      fill = "Probability",
      title = paste("Posterior similarity of", if (target == "item") "items" else "rankings"),
      subtitle = "Probability of belonging to the same latent group"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )
}

#' Plot saved MCMC traces
#'
#' Trace plots make it possible to see whether saved MCMC draws move around a
#' stable region. They are a first diagnostic, not a proof of convergence.
#'
#' @param fit A btsbm_fit object.
#' @param quantities Optional character vector of trace names to show.
#'
#' @return A ggplot object.
#' @export
plot_mcmc_traces <- function(fit, quantities = NULL) {
  .require_ggplot2()
  if (!inherits(fit, "btsbm_fit")) {
    stop("'fit' must be a btsbm_fit object.", call. = FALSE)
  }
  traces <- .fit_diagnostic_traces(fit)
  if (!is.null(quantities)) {
    missing_quantities <- setdiff(quantities, names(traces))
    if (length(missing_quantities)) {
      stop("No trace named: ", paste(missing_quantities, collapse = ", "), ".", call. = FALSE)
    }
    traces <- traces[quantities]
  }
  if (!length(traces)) {
    stop("This fit has no scalar traces to plot.", call. = FALSE)
  }
  trace_data <- do.call(rbind, lapply(names(traces), function(quantity) {
    data.frame(
      saved_draw = seq_along(traces[[quantity]]),
      value = as.numeric(traces[[quantity]]),
      quantity = quantity,
      stringsAsFactors = FALSE
    )
  }))
  trace_data$quantity <- factor(trace_data$quantity, levels = names(traces))

  ggplot2::ggplot(trace_data, ggplot2::aes(x = .data$saved_draw, y = .data$value)) +
    ggplot2::geom_line(linewidth = 0.35, colour = "#0072B2") +
    ggplot2::facet_wrap(~quantity, scales = "free_y", ncol = 1) +
    ggplot2::labs(
      x = "Saved MCMC draw",
      y = NULL,
      title = "MCMC traces",
      caption = "Look for sustained movement without long-term drift or a stuck trace."
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

#' Plot how often items appear at each ranking position
#'
#' This descriptive heatmap is useful before fitting a PL model. It shows which
#' items commonly appear near the top of observed rankings. It is not adjusted
#' for uncertainty or latent grouping.
#'
#' @param rankings A btsbm_ranking_data object created by as_rankings().
#'
#' @return A ggplot object.
#' @export
plot_ranking_positions <- function(rankings) {
  .require_ggplot2()
  if (!inherits(rankings, "btsbm_ranking_data")) {
    stop("'rankings' must be created by as_rankings().", call. = FALSE)
  }
  ranking_matrix <- rankings$rankings
  count_matrix <- matrix(
    0L, nrow = rankings$n_items, ncol = ncol(ranking_matrix),
    dimnames = list(rankings$item_labels, paste0("Position ", seq_len(ncol(ranking_matrix))))
  )
  for (position in seq_len(ncol(ranking_matrix))) {
    count_matrix[, position] <- tabulate(
      ranking_matrix[, position], nbins = rankings$n_items
    )
  }
  plot_data <- expand.grid(
    item = rankings$item_labels,
    position = colnames(count_matrix),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  plot_data$count <- as.vector(count_matrix)
  item_order <- names(sort(rowSums(count_matrix * rev(seq_len(ncol(count_matrix)))),
                           decreasing = TRUE))
  plot_data$item <- factor(plot_data$item, levels = rev(item_order))
  plot_data$position <- factor(plot_data$position, levels = colnames(count_matrix))

  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$position, y = .data$item)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$count), colour = "white", linewidth = 0.15) +
    ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#238B45") +
    ggplot2::labs(
      x = "Reported rank position",
      y = NULL,
      fill = "Rankings",
      title = "Observed ranking positions",
      subtitle = "Darker cells indicate items that appeared more often in that position"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

#' Plot observed pairwise win counts
#'
#' This descriptive heatmap is useful before fitting a Bradley--Terry model.
#' Rows are winners and columns are losers; it displays the number of observed
#' wins, not a fitted probability.
#'
#' @param pairwise A btsbm_pairwise_data object created by as_bt_data().
#'
#' @return A ggplot object.
#' @export
plot_pairwise_outcomes <- function(pairwise) {
  .require_ggplot2()
  if (!inherits(pairwise, "btsbm_pairwise_data")) {
    stop("'pairwise' must be created by as_bt_data().", call. = FALSE)
  }
  labels <- pairwise$item_labels
  plot_data <- expand.grid(
    winner = labels, loser = labels,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  plot_data$wins <- as.vector(pairwise$wins)
  plot_data$winner <- factor(plot_data$winner, levels = rev(labels))
  plot_data$loser <- factor(plot_data$loser, levels = labels)

  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$loser, y = .data$winner)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$wins), colour = "white", linewidth = 0.15) +
    ggplot2::scale_fill_gradient(low = "#FFF7EC", high = "#D7301F") +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      x = "Loser",
      y = "Winner",
      fill = "Wins",
      title = "Observed pairwise outcomes",
      subtitle = "Rows record wins over the column item"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )
}
