# Internal Plackett--Luce helpers.  These functions use the exponential-race
# augmentation of Caron and Doucet (2012).  The first implementation is kept in
# R for auditability and small examples; the matching PLuce C++ kernels can
# replace the inner loops without changing this interface.

.fit_pl_model <- function(data, model, control) {
  if (model$structure == "none") return(.gibbs_pl(data, model, control))
  if (model$structure == "item_sbm") return(.gibbs_pl_sbm(data, model, control))
  if (model$structure == "ranking_mixture") return(.gibbs_pl_mixture(data, model, control))
  if (model$structure == "lbm") return(.gibbs_pl_lbm(data, model, control))
  stop("Unknown PL model structure.", call. = FALSE)
}

.pl_draw_augmentation <- function(rankings, item_strength) {
  L <- nrow(rankings)
  m <- ncol(rankings)
  n <- length(item_strength)
  Z <- matrix(0, nrow = L, ncol = m)
  S_i <- numeric(n)
  total_rate <- sum(item_strength)

  for (ell in seq_len(L)) {
    cumulative_weight <- 0
    cumulative_time <- 0
    for (q in seq_len(m)) {
      item <- rankings[ell, q]
      rate <- total_rate - cumulative_weight
      if (!is.finite(rate) || rate <= 0) {
        stop("PL risk-set rate became non-positive.", call. = FALSE)
      }
      z <- stats::rexp(1L, rate = rate)
      Z[ell, q] <- z
      cumulative_time <- cumulative_time + z
      S_i[item] <- S_i[item] + cumulative_time
      cumulative_weight <- cumulative_weight + item_strength[item]
    }
    unranked <- setdiff(seq_len(n), rankings[ell, ])
    if (length(unranked)) S_i[unranked] <- S_i[unranked] + cumulative_time
  }
  list(Z = Z, exposure = S_i)
}

.pl_win_counts <- function(rankings, n_items) tabulate(as.integer(rankings), nbins = n_items)

.normalise_strength <- function(strength, method) {
  strength <- pmax(as.numeric(strength), .Machine$double.eps)
  if (identical(method, "none")) return(strength)
  strength / exp(mean(log(strength)))
}

.pl_log_lik_one <- function(rankings, item_strength) {
  L <- nrow(rankings)
  m <- ncol(rankings)
  out <- numeric(L)
  for (ell in seq_len(L)) {
    remaining <- sum(item_strength)
    for (q in seq_len(m)) {
      item <- rankings[ell, q]
      out[ell] <- out[ell] + log(item_strength[item]) - log(remaining)
      remaining <- remaining - item_strength[item]
    }
  }
  out
}

.gibbs_pl <- function(data, model, control) {
  rankings <- data$rankings
  n <- data$n_items
  prior <- model$latent_strength
  S <- control$iter - control$warmup
  wins <- .pl_win_counts(rankings, n)

  strength <- stats::rgamma(n, shape = prior$shape + 1, rate = prior$rate + 1)
  latent_strength_draws <- matrix(NA_real_, nrow = S, ncol = n)
  item_strength_draws <- matrix(NA_real_, nrow = S, ncol = n)
  log_likelihood_trace <- numeric(S)
  z_store <- if (identical(control$store, "augmented")) {
    array(NA_real_, dim = c(S, nrow(rankings), ncol(rankings)))
  } else NULL

  save_i <- 0L
  for (iter in seq_len(control$iter)) {
    augmented <- .pl_draw_augmentation(rankings, strength)
    strength <- stats::rgamma(n, shape = prior$shape + wins, rate = prior$rate + augmented$exposure)

    if (iter > control$warmup) {
      save_i <- save_i + 1L
      latent_strength_draws[save_i, ] <- strength
      item_strength_draws[save_i, ] <- .normalise_strength(strength, model$identifiability$method)
      log_likelihood_trace[save_i] <- sum(.pl_log_lik_one(rankings, strength))
      if (!is.null(z_store)) z_store[save_i, , ] <- augmented$Z
    }
    if (control$verbose && iter %% control$progress_every == 0L) {
      cat("PL | iter", iter, "\n")
    }
  }
  colnames(latent_strength_draws) <- data$item_labels
  colnames(item_strength_draws) <- data$item_labels
  .new_fit(
    draws = list(
      latent_strength = latent_strength_draws,
      item_strength = item_strength_draws,
      ability = latent_strength_draws,
      item_ability = item_strength_draws,
      augmented = z_store
    ),
    diagnostics = list(
      identifiability = model$identifiability$method,
      traces = list(log_likelihood = log_likelihood_trace)
    )
  )
}

.compact_partition <- function(x) match(x, sort(unique(x)))

.initial_item_partition <- function(n_items, prior, init = NULL) {
  if (!is.null(init)) {
    if (length(init) != n_items || any(!is.finite(init)) || any(init < 1) || any(init != floor(init))) {
      stop("`init_item_cluster` must contain one positive integer label per item.", call. = FALSE)
    }
    x <- .compact_partition(as.integer(init))
  } else {
    initial_K <- if (prior$type == "dm") min(2L, prior$K, n_items) else min(2L, n_items)
    x <- sample.int(initial_K, n_items, replace = TRUE)
    x <- .compact_partition(x)
  }
  if (prior$type == "dm" && length(unique(x)) > prior$K) {
    stop("The initial item partition exceeds the finite prior's `K`.", call. = FALSE)
  }
  x
}

.partition_predictive_weights <- function(sizes, prior) {
  switch(
    prior$type,
    gnedin = urn_GN(sizes, prior$gamma),
    dp = urn_DP(sizes, prior$alpha),
    py = urn_PY(sizes, prior$alpha, prior$sigma),
    dm = urn_DM(sizes, prior$beta, prior$K),
    stop("Unsupported partition prior.", call. = FALSE)
  )
}

.collapsed_pl_predictive <- function(w_base, s_base, w_add, s_add, shape, rate) {
  alpha <- shape + w_base
  beta <- rate + s_base
  lgamma(alpha + w_add) - lgamma(alpha) + alpha * log(beta) -
    (alpha + w_add) * log(beta + s_add)
}

.draw_log_probabilities <- function(log_weights) {
  if (all(!is.finite(log_weights))) {
    stop("All PL--SBM allocation probabilities are zero.", call. = FALSE)
  }
  max_log <- max(log_weights)
  probability <- exp(log_weights - max_log)
  probability[!is.finite(probability)] <- 0
  probability <- probability / sum(probability)
  sample.int(length(probability), size = 1L, prob = probability)
}

.gibbs_pl_sbm <- function(data, model, control) {
  rankings <- data$rankings
  n <- data$n_items
  prior <- model$latent_strength
  partition <- model$item_clustering
  S <- control$iter - control$warmup
  wins <- .pl_win_counts(rankings, n)
  x <- .initial_item_partition(n, partition, control$init_item_cluster)
  K_curr <- length(unique(x))
  strength <- stats::rgamma(K_curr, shape = prior$shape + 1, rate = prior$rate + 1)

  item_cluster <- matrix(NA_integer_, nrow = S, ncol = n)
  K_trace <- integer(S)
  latent_strength_draws <- vector("list", S)
  item_strength_draws <- matrix(NA_real_, nrow = S, ncol = n)
  log_likelihood_trace <- numeric(S)
  z_store <- if (identical(control$store, "augmented")) {
    array(NA_real_, dim = c(S, nrow(rankings), ncol(rankings)))
  } else NULL

  save_i <- 0L
  for (iter in seq_len(control$iter)) {
    augmented <- .pl_draw_augmentation(rankings, strength[x])
    exposure <- augmented$exposure

    for (i in seq_len(n)) {
      x_without_i <- x[-i]
      occupied <- sort(unique(x_without_i))
      sizes <- tabulate(match(x_without_i, occupied), nbins = length(occupied))
      weights <- .partition_predictive_weights(sizes, partition)
      H <- length(occupied)
      logp <- rep(-Inf, H + 1L)

      if (H > 0L) {
        for (h in seq_len(H)) {
          members <- which(x_without_i == occupied[h])
          source_index <- (seq_len(n))[-i][members]
          contribution <- .collapsed_pl_predictive(
            sum(wins[source_index]), sum(exposure[source_index]), wins[i], exposure[i],
            prior$shape, prior$rate
          )
          logp[h] <- log(weights[h]) + contribution
        }
      }
      if (weights[H + 1L] > 0) {
        logp[H + 1L] <- log(weights[H + 1L]) +
          .collapsed_pl_predictive(0, 0, wins[i], exposure[i], prior$shape, prior$rate)
      }

      draw <- .draw_log_probabilities(logp)
      x[i] <- if (draw <= H) occupied[draw] else if (H == 0L) 1L else max(occupied) + 1L
      x <- .compact_partition(x)
    }

    K_curr <- length(unique(x))
    strength <- numeric(K_curr)
    for (k in seq_len(K_curr)) {
      members <- which(x == k)
      strength[k] <- stats::rgamma(
        1L,
        shape = prior$shape + sum(wins[members]),
        rate = prior$rate + sum(exposure[members])
      )
    }

    if (iter > control$warmup) {
      save_i <- save_i + 1L
      item_cluster[save_i, ] <- x
      K_trace[save_i] <- K_curr
      latent_strength_draws[[save_i]] <- strength
      item_strength_draws[save_i, ] <- .normalise_strength(strength[x], model$identifiability$method)
      log_likelihood_trace[save_i] <- sum(.pl_log_lik_one(rankings, strength[x]))
      if (!is.null(z_store)) z_store[save_i, , ] <- augmented$Z
    }
    if (control$verbose && iter %% control$progress_every == 0L) {
      cat("PL--SBM | iter", iter, "occupied item blocks =", K_curr, "\n")
    }
  }
  colnames(item_cluster) <- data$item_labels
  colnames(item_strength_draws) <- data$item_labels
  .new_fit(
    draws = list(
      item_cluster = item_cluster,
      n_item_clusters = K_trace,
      latent_strength = latent_strength_draws,
      item_strength = item_strength_draws,
      ability = latent_strength_draws,
      item_ability = item_strength_draws,
      augmented = z_store
    ),
    diagnostics = list(
      identifiability = model$identifiability$method,
      traces = list(
        n_item_clusters = K_trace,
        log_likelihood = log_likelihood_trace
      )
    )
  )
}

.pl_log_lik_fit <- function(fit) {
  rankings <- fit$data$rankings
  S <- .fit_n_draws(fit)
  out <- matrix(NA_real_, S, nrow(rankings))
  for (s in seq_len(S)) {
    item_strength <- if (fit$model$structure == "none") {
      fit$draws$latent_strength[s, ]
    } else if (fit$model$structure == "item_sbm") {
      fit$draws$latent_strength[[s]][fit$draws$item_cluster[s, ]]
    }
    if (fit$model$structure %in% c("none", "item_sbm")) {
      out[s, ] <- .pl_log_lik_one(rankings, item_strength)
    } else if (fit$model$structure == "ranking_mixture") {
      out[s, ] <- .pl_mixture_log_likelihood(
        rankings, fit$draws$ranking_cluster[s, ], fit$draws$latent_strength[[s]]
      )
    } else if (fit$model$structure == "lbm") {
      out[s, ] <- .pl_lbm_log_likelihood(
        rankings, fit$draws$ranking_cluster[s, ], fit$draws$item_cluster[s, ],
        fit$draws$latent_strength[[s]]
      )
    } else {
      stop("PL likelihood extraction is not implemented for this structure.", call. = FALSE)
    }
  }
  attr(out, "unit_index") <- fit$data$ranking_id
  out
}
