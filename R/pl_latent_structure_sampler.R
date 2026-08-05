# Plackett--Luce latent-structure samplers.
#
# These follow the exponential-race augmentation used in PLuce, but use the
# package's prepared data objects and descriptive state names. The allocation
# sweeps remain readable R; the shared augmentation/exposure kernel is a
# registered Rcpp implementation checked against the R reference below.

.pl_win_matrix <- function(rankings, n_items) {
  wins <- matrix(0, nrow = nrow(rankings), ncol = n_items)
  for (ranking_index in seq_len(nrow(rankings))) {
    wins[ranking_index, ] <- tabulate(rankings[ranking_index, ], nbins = n_items)
  }
  wins
}

.pl_draw_grouped_augmentation_r <- function(rankings, ranking_cluster, latent_strength) {
  rankings <- as.matrix(rankings)
  latent_strength <- as.matrix(latent_strength)
  n_rankings <- nrow(rankings)
  n_items <- ncol(latent_strength)
  rank_length <- ncol(rankings)
  if (length(ranking_cluster) != n_rankings || any(ranking_cluster < 1L) ||
      any(ranking_cluster > nrow(latent_strength))) {
    stop("Ranking clusters and latent strengths are incompatible.", call. = FALSE)
  }

  augmentation <- matrix(0, nrow = n_rankings, ncol = rank_length)
  exposure <- matrix(0, nrow = n_rankings, ncol = n_items)
  for (ranking_index in seq_len(n_rankings)) {
    item_strength <- latent_strength[ranking_cluster[ranking_index], ]
    total_strength <- sum(item_strength)
    chosen_strength <- 0
    elapsed <- 0
    for (position in seq_len(rank_length)) {
      item <- rankings[ranking_index, position]
      risk_set_strength <- total_strength - chosen_strength
      if (!is.finite(risk_set_strength) || risk_set_strength <= 0) {
        stop("PL risk-set strength became non-positive.", call. = FALSE)
      }
      waiting_time <- stats::rexp(1L, rate = risk_set_strength)
      augmentation[ranking_index, position] <- waiting_time
      elapsed <- elapsed + waiting_time
      exposure[ranking_index, item] <- exposure[ranking_index, item] + elapsed
      chosen_strength <- chosen_strength + item_strength[item]
    }
    unranked_items <- setdiff(seq_len(n_items), rankings[ranking_index, ])
    if (length(unranked_items)) {
      exposure[ranking_index, unranked_items] <-
        exposure[ranking_index, unranked_items] + elapsed
    }
  }
  list(augmentation = augmentation, exposure = exposure)
}

.pl_draw_grouped_augmentation <- function(rankings, ranking_cluster, latent_strength) {
  pl_grouped_augmentation_cpp(
    rankings = matrix(as.integer(rankings), nrow = nrow(rankings)),
    ranking_cluster = as.integer(ranking_cluster),
    latent_strength = as.matrix(latent_strength)
  )
}

.initial_latent_partition <- function(n_entities, partition_prior, initial = NULL,
                                      initial_name = "initial cluster") {
  if (!is.null(initial)) {
    if (length(initial) != n_entities || any(!is.finite(initial)) ||
        any(initial < 1L) || any(initial != floor(initial))) {
      stop("`", initial_name, "` must contain one positive integer label per entity.",
           call. = FALSE)
    }
    labels <- .compact_partition(as.integer(initial))
  } else {
    initial_clusters <- if (partition_prior$type == "dm") {
      min(2L, partition_prior$K, n_entities)
    } else {
      min(2L, n_entities)
    }
    labels <- .compact_partition(sample.int(initial_clusters, n_entities, replace = TRUE))
  }
  if (partition_prior$type == "dm" && length(unique(labels)) > partition_prior$K) {
    stop("The initial partition exceeds `max_clusters`.", call. = FALSE)
  }
  labels
}

.collapsed_strength_increment <- function(base_wins, base_exposure, added_wins,
                                          added_exposure, shape, rate) {
  sum(.collapsed_pl_predictive(
    base_wins, base_exposure, added_wins, added_exposure, shape, rate
  ))
}

.reassign_partition <- function(labels, observation_wins, observation_exposure,
                                partition_prior, shape, rate) {
  n_entities <- length(labels)
  for (entity_index in seq_len(n_entities)) {
    other_indices <- (seq_len(n_entities))[-entity_index]
    other_labels <- labels[-entity_index]
    occupied <- sort(unique(other_labels))
    n_occupied <- length(occupied)
    cluster_sizes <- tabulate(match(other_labels, occupied), nbins = n_occupied)
    allocation_weights <- .partition_predictive_weights(cluster_sizes, partition_prior)
    log_probability <- rep(-Inf, n_occupied + 1L)

    if (n_occupied > 0L) {
      for (cluster_index in seq_len(n_occupied)) {
        members <- other_indices[other_labels == occupied[cluster_index]]
        log_probability[cluster_index] <- log(allocation_weights[cluster_index]) +
          .collapsed_strength_increment(
            colSums(observation_wins[members, , drop = FALSE]),
            colSums(observation_exposure[members, , drop = FALSE]),
            observation_wins[entity_index, ],
            observation_exposure[entity_index, ],
            shape, rate
          )
      }
    }
    if (allocation_weights[n_occupied + 1L] > 0) {
      log_probability[n_occupied + 1L] <- log(allocation_weights[n_occupied + 1L]) +
        .collapsed_strength_increment(
          rep(0, ncol(observation_wins)), rep(0, ncol(observation_exposure)),
          observation_wins[entity_index, ], observation_exposure[entity_index, ],
          shape, rate
        )
    }
    draw <- .draw_log_probabilities(log_probability)
    labels[entity_index] <- if (draw <= n_occupied) {
      occupied[draw]
    } else if (n_occupied == 0L) {
      1L
    } else {
      max(occupied) + 1L
    }
    labels <- .compact_partition(labels)
  }
  labels
}

.sample_group_strength <- function(labels, observation_wins, observation_exposure,
                                   shape, rate) {
  n_groups <- length(unique(labels))
  n_features <- ncol(observation_wins)
  out <- matrix(NA_real_, nrow = n_groups, ncol = n_features)
  for (group_index in seq_len(n_groups)) {
    members <- which(labels == group_index)
    out[group_index, ] <- stats::rgamma(
      n_features,
      shape = shape + colSums(observation_wins[members, , drop = FALSE]),
      rate = rate + colSums(observation_exposure[members, , drop = FALSE])
    )
  }
  out
}

.normalise_strength_rows <- function(latent_strength, method) {
  t(vapply(
    seq_len(nrow(latent_strength)),
    function(group_index) .normalise_strength(latent_strength[group_index, ], method),
    numeric(ncol(latent_strength))
  ))
}

.pl_mixture_log_likelihood <- function(rankings, ranking_cluster, latent_strength) {
  out <- numeric(nrow(rankings))
  for (ranking_index in seq_len(nrow(rankings))) {
    out[ranking_index] <- .pl_log_lik_one(
      rankings[ranking_index, , drop = FALSE],
      latent_strength[ranking_cluster[ranking_index], ]
    )
  }
  out
}

.gibbs_pl_mixture <- function(data, model, control) {
  rankings <- data$rankings
  prior <- model$latent_strength
  ranking_prior <- model$ranking_clustering
  n_rankings <- nrow(rankings)
  n_items <- data$n_items
  n_saved <- control$iter - control$warmup
  ranking_wins <- .pl_win_matrix(rankings, n_items)
  ranking_cluster <- .initial_latent_partition(
    n_rankings, ranking_prior, control$init_ranking_cluster, "init_ranking_cluster"
  )
  latent_strength <- matrix(
    stats::rgamma(length(unique(ranking_cluster)) * n_items,
                  shape = prior$shape + 1, rate = prior$rate + 1),
    ncol = n_items
  )

  ranking_cluster_draws <- matrix(NA_integer_, nrow = n_saved, ncol = n_rankings)
  cluster_count_trace <- integer(n_saved)
  latent_strength_draws <- vector("list", n_saved)
  reported_strength_draws <- vector("list", n_saved)
  log_likelihood_trace <- numeric(n_saved)
  augmentation_draws <- if (identical(control$store, "augmented")) {
    array(NA_real_, dim = c(n_saved, n_rankings, ncol(rankings)))
  } else NULL

  saved_index <- 0L
  for (iteration in seq_len(control$iter)) {
    augmented <- .pl_draw_grouped_augmentation(
      rankings, ranking_cluster, latent_strength
    )
    ranking_cluster <- .reassign_partition(
      ranking_cluster, ranking_wins, augmented$exposure, ranking_prior,
      prior$shape, prior$rate
    )
    latent_strength <- .sample_group_strength(
      ranking_cluster, ranking_wins, augmented$exposure, prior$shape, prior$rate
    )

    if (iteration > control$warmup) {
      saved_index <- saved_index + 1L
      ranking_cluster_draws[saved_index, ] <- ranking_cluster
      cluster_count_trace[saved_index] <- nrow(latent_strength)
      latent_strength_draws[[saved_index]] <- latent_strength
      reported_strength_draws[[saved_index]] <- .normalise_strength_rows(
        latent_strength, model$identifiability$method
      )
      log_likelihood_trace[saved_index] <- sum(.pl_mixture_log_likelihood(
        rankings, ranking_cluster, latent_strength
      ))
      if (!is.null(augmentation_draws)) {
        augmentation_draws[saved_index, , ] <- augmented$augmentation
      }
    }
    if (control$verbose && iteration %% control$progress_every == 0L) {
      cat("PL mixture | iter", iteration, "occupied ranking groups =",
          length(unique(ranking_cluster)), "\n")
    }
  }

  colnames(ranking_cluster_draws) <- data$ranking_labels
  .new_fit(
    draws = list(
      ranking_cluster = ranking_cluster_draws,
      n_ranking_clusters = cluster_count_trace,
      latent_strength = latent_strength_draws,
      component_strength = reported_strength_draws,
      augmented = augmentation_draws
    ),
    diagnostics = list(
      identifiability = model$identifiability$method,
      traces = list(
        n_ranking_clusters = cluster_count_trace,
        log_likelihood = log_likelihood_trace
      )
    )
  )
}

.row_block_statistics <- function(ranking_wins, ranking_exposure, item_cluster,
                                  n_item_clusters) {
  n_rankings <- nrow(ranking_wins)
  wins <- matrix(0, nrow = n_rankings, ncol = n_item_clusters)
  exposure <- matrix(0, nrow = n_rankings, ncol = n_item_clusters)
  for (item_cluster_index in seq_len(n_item_clusters)) {
    items <- which(item_cluster == item_cluster_index)
    wins[, item_cluster_index] <- rowSums(ranking_wins[, items, drop = FALSE])
    exposure[, item_cluster_index] <- rowSums(ranking_exposure[, items, drop = FALSE])
  }
  list(wins = wins, exposure = exposure)
}

.item_block_statistics <- function(ranking_wins, ranking_exposure, ranking_cluster,
                                   n_ranking_clusters) {
  n_items <- ncol(ranking_wins)
  wins <- matrix(0, nrow = n_items, ncol = n_ranking_clusters)
  exposure <- matrix(0, nrow = n_items, ncol = n_ranking_clusters)
  for (ranking_cluster_index in seq_len(n_ranking_clusters)) {
    rankings_in_cluster <- which(ranking_cluster == ranking_cluster_index)
    wins[, ranking_cluster_index] <- colSums(
      ranking_wins[rankings_in_cluster, , drop = FALSE]
    )
    exposure[, ranking_cluster_index] <- colSums(
      ranking_exposure[rankings_in_cluster, , drop = FALSE]
    )
  }
  list(wins = wins, exposure = exposure)
}

.sample_lbm_strength <- function(ranking_cluster, item_cluster, ranking_wins,
                                 ranking_exposure, shape, rate) {
  n_ranking_clusters <- length(unique(ranking_cluster))
  n_item_clusters <- length(unique(item_cluster))
  latent_strength <- matrix(NA_real_, nrow = n_ranking_clusters,
                            ncol = n_item_clusters)
  for (ranking_cluster_index in seq_len(n_ranking_clusters)) {
    ranking_indices <- which(ranking_cluster == ranking_cluster_index)
    for (item_cluster_index in seq_len(n_item_clusters)) {
      item_indices <- which(item_cluster == item_cluster_index)
      latent_strength[ranking_cluster_index, item_cluster_index] <- stats::rgamma(
        1L,
        shape = shape + sum(ranking_wins[ranking_indices, item_indices, drop = FALSE]),
        rate = rate + sum(ranking_exposure[ranking_indices, item_indices, drop = FALSE])
      )
    }
  }
  latent_strength
}

.pl_lbm_log_likelihood <- function(rankings, ranking_cluster, item_cluster,
                                   latent_strength) {
  out <- numeric(nrow(rankings))
  for (ranking_index in seq_len(nrow(rankings))) {
    item_strength <- latent_strength[
      ranking_cluster[ranking_index], item_cluster
    ]
    out[ranking_index] <- .pl_log_lik_one(
      rankings[ranking_index, , drop = FALSE], item_strength
    )
  }
  out
}

.gibbs_pl_lbm <- function(data, model, control) {
  rankings <- data$rankings
  prior <- model$latent_strength
  n_rankings <- nrow(rankings)
  n_items <- data$n_items
  n_saved <- control$iter - control$warmup
  ranking_wins <- .pl_win_matrix(rankings, n_items)
  ranking_cluster <- .initial_latent_partition(
    n_rankings, model$ranking_clustering, control$init_ranking_cluster,
    "init_ranking_cluster"
  )
  item_cluster <- .initial_latent_partition(
    n_items, model$item_clustering, control$init_item_cluster, "init_item_cluster"
  )
  latent_strength <- matrix(
    stats::rgamma(length(unique(ranking_cluster)) * length(unique(item_cluster)),
                  shape = prior$shape + 1, rate = prior$rate + 1),
    nrow = length(unique(ranking_cluster)), ncol = length(unique(item_cluster))
  )

  ranking_cluster_draws <- matrix(NA_integer_, nrow = n_saved, ncol = n_rankings)
  item_cluster_draws <- matrix(NA_integer_, nrow = n_saved, ncol = n_items)
  ranking_cluster_count_trace <- integer(n_saved)
  item_cluster_count_trace <- integer(n_saved)
  latent_strength_draws <- vector("list", n_saved)
  reported_strength_draws <- vector("list", n_saved)
  log_likelihood_trace <- numeric(n_saved)
  augmentation_draws <- if (identical(control$store, "augmented")) {
    array(NA_real_, dim = c(n_saved, n_rankings, ncol(rankings)))
  } else NULL

  saved_index <- 0L
  for (iteration in seq_len(control$iter)) {
    augmented <- .pl_draw_grouped_augmentation(
      rankings, ranking_cluster, latent_strength[, item_cluster, drop = FALSE]
    )

    row_statistics <- .row_block_statistics(
      ranking_wins, augmented$exposure, item_cluster, length(unique(item_cluster))
    )
    ranking_cluster <- .reassign_partition(
      ranking_cluster, row_statistics$wins, row_statistics$exposure,
      model$ranking_clustering, prior$shape, prior$rate
    )

    item_statistics <- .item_block_statistics(
      ranking_wins, augmented$exposure, ranking_cluster,
      length(unique(ranking_cluster))
    )
    item_cluster <- .reassign_partition(
      item_cluster, item_statistics$wins, item_statistics$exposure,
      model$item_clustering, prior$shape, prior$rate
    )
    latent_strength <- .sample_lbm_strength(
      ranking_cluster, item_cluster, ranking_wins, augmented$exposure,
      prior$shape, prior$rate
    )

    if (iteration > control$warmup) {
      saved_index <- saved_index + 1L
      ranking_cluster_draws[saved_index, ] <- ranking_cluster
      item_cluster_draws[saved_index, ] <- item_cluster
      ranking_cluster_count_trace[saved_index] <- length(unique(ranking_cluster))
      item_cluster_count_trace[saved_index] <- length(unique(item_cluster))
      latent_strength_draws[[saved_index]] <- latent_strength
      reported_strength_draws[[saved_index]] <- .normalise_strength_rows(
        latent_strength, model$identifiability$method
      )
      log_likelihood_trace[saved_index] <- sum(.pl_lbm_log_likelihood(
        rankings, ranking_cluster, item_cluster, latent_strength
      ))
      if (!is.null(augmentation_draws)) {
        augmentation_draws[saved_index, , ] <- augmented$augmentation
      }
    }
    if (control$verbose && iteration %% control$progress_every == 0L) {
      cat("PL LBM | iter", iteration, "ranking groups =",
          length(unique(ranking_cluster)), "item groups =",
          length(unique(item_cluster)), "\n")
    }
  }

  colnames(ranking_cluster_draws) <- data$ranking_labels
  colnames(item_cluster_draws) <- data$item_labels
  .new_fit(
    draws = list(
      ranking_cluster = ranking_cluster_draws,
      item_cluster = item_cluster_draws,
      n_ranking_clusters = ranking_cluster_count_trace,
      n_item_clusters = item_cluster_count_trace,
      latent_strength = latent_strength_draws,
      cell_strength = reported_strength_draws,
      augmented = augmentation_draws
    ),
    diagnostics = list(
      identifiability = model$identifiability$method,
      traces = list(
        n_ranking_clusters = ranking_cluster_count_trace,
        n_item_clusters = item_cluster_count_trace,
        log_likelihood = log_likelihood_trace
      )
    )
  )
}
