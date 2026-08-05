#' Create validated Bradley--Terry pairwise data
#'
#' Converts a square directed win-count matrix into the package's validated
#' pairwise-data representation.  This is the boundary between the legacy
#' matrix API and the common fitting interface.
#'
#' @param wins Numeric square matrix.  Entry `wins[i, j]` is the number of
#'   wins for item `i` over item `j`.
#' @param item_labels Optional character vector of item names.  By default,
#'   matrix row names are used when available.
#'
#' @return An object of class `btsbm_pairwise_data`.
#' @export
as_bt_data <- function(wins, item_labels = NULL) {
  if (!is.matrix(wins) || nrow(wins) != ncol(wins) || nrow(wins) < 2L) {
    stop("`wins` must be a square matrix with at least two items.", call. = FALSE)
  }
  if (!is.numeric(wins) || any(!is.finite(wins))) {
    stop("`wins` must contain finite numeric counts.", call. = FALSE)
  }
  if (any(wins < 0) || any(diag(wins) != 0)) {
    stop("`wins` must be non-negative with a zero diagonal.", call. = FALSE)
  }

  n_items <- nrow(wins)
  if (is.null(item_labels)) item_labels <- rownames(wins)
  if (is.null(item_labels)) item_labels <- paste0("Item_", seq_len(n_items))
  if (!is.character(item_labels) || length(item_labels) != n_items ||
      anyNA(item_labels) || anyDuplicated(item_labels)) {
    stop("`item_labels` must be a unique, non-missing character vector of length n_items.", call. = FALSE)
  }

  dimnames(wins) <- list(item_labels, item_labels)
  structure(
    list(wins = wins, n_items = n_items, item_labels = item_labels),
    class = "btsbm_pairwise_data"
  )
}

#' Create validated Plackett--Luce ranking data
#'
#' Rows of `rankings` are preference orderings and columns are successive
#' positions, from most to least preferred. Each row may be a complete ranking
#' or a top-`m` prefix from one common universe of items. Internally this is
#' represented in the compact integer format used by the PLuce reference code,
#' but users never need to work with that internal name.
#'
#' @param rankings Integer-like matrix. Row `l` lists item ids from most to
#'   least preferred in ranking `l`.
#' @param item_count Total number of available items. If `NULL`, it is inferred
#'   as `max(rankings)`; supply it when a top-`m` dataset has never-ranked
#'   items.
#' @param ranking_labels Optional unique character vector identifying ranking rows.
#' @param item_labels Optional unique character vector of length `item_count`.
#'
#' @return An object of class `btsbm_ranking_data`.
#' @references Caron, F. and Doucet, A. (2012). Efficient Bayesian inference
#'   for generalized Bradley--Terry models. *Journal of Computational and
#'   Graphical Statistics*, 21(1), 174--196.
#'   \doi{10.1080/10618600.2012.638220}
#' @export
as_rankings <- function(rankings, item_count = NULL, ranking_labels = NULL,
                        item_labels = NULL) {
  if (!is.matrix(rankings) || nrow(rankings) < 1L || ncol(rankings) < 1L) {
    stop("`rankings` must be a non-empty matrix with one row per ranking.", call. = FALSE)
  }
  if (!is.numeric(rankings) || any(!is.finite(rankings)) || any(rankings != floor(rankings))) {
    stop("`rankings` must contain finite integer item ids.", call. = FALSE)
  }
  rankings <- matrix(as.integer(rankings), nrow = nrow(rankings), ncol = ncol(rankings),
                     dimnames = dimnames(rankings))

  if (is.null(item_count)) item_count <- max(rankings)
  if (length(item_count) != 1L || !is.finite(item_count) || item_count != floor(item_count) || item_count < 2L) {
    stop("`item_count` must be one integer of at least two.", call. = FALSE)
  }
  n_items <- as.integer(item_count)
  if (any(rankings < 1L | rankings > n_items)) {
    stop("Every `rankings` entry must lie in 1:item_count.", call. = FALSE)
  }
  if (any(apply(rankings, 1L, anyDuplicated))) {
    stop("Each ranking row in `rankings` must contain distinct item ids.", call. = FALSE)
  }

  n_rankings <- nrow(rankings)
  if (is.null(ranking_labels)) ranking_labels <- rownames(rankings)
  if (is.null(ranking_labels)) ranking_labels <- paste0("Ranking_", seq_len(n_rankings))
  if (!is.character(ranking_labels) || length(ranking_labels) != n_rankings ||
      anyNA(ranking_labels) || anyDuplicated(ranking_labels)) {
    stop("`ranking_labels` must be a unique, non-missing character vector of length nrow(rankings).", call. = FALSE)
  }

  if (is.null(item_labels)) item_labels <- paste0("Item_", seq_len(n_items))
  if (!is.character(item_labels) || length(item_labels) != n_items ||
      anyNA(item_labels) || anyDuplicated(item_labels)) {
    stop("`item_labels` must be a unique, non-missing character vector of length n_items.", call. = FALSE)
  }

  rownames(rankings) <- ranking_labels
  colnames(rankings) <- paste0("Rank_", seq_len(ncol(rankings)))
  structure(
    list(
      rankings = rankings,
      rho = rankings,
      n_items = n_items,
      n_rankings = n_rankings,
      rank_length = ncol(rankings),
      ranking_labels = ranking_labels,
      ranking_id = ranking_labels,
      item_labels = item_labels
    ),
    class = "btsbm_ranking_data"
  )
}

#' @rdname as_rankings
#' @param rho Deprecated name for `rankings`.
#' @param n_items Deprecated name for `item_count`.
#' @param ranking_id Deprecated name for `ranking_labels`.
#' @export
as_pl_data <- function(rho, n_items = NULL, ranking_id = NULL, item_labels = NULL) {
  as_rankings(
    rankings = rho,
    item_count = n_items,
    ranking_labels = ranking_id,
    item_labels = item_labels
  )
}

#' @export
print.btsbm_pairwise_data <- function(x, ...) {
  cat("<btsbm_pairwise_data>", x$n_items, "items and", sum(x$wins), "wins\n")
  invisible(x)
}

#' @export
print.btsbm_ranking_data <- function(x, ...) {
  cat("<btsbm_ranking_data>", x$n_rankings, "rankings of length", x$rank_length,
      "from", x$n_items, "items\n")
  invisible(x)
}

#' @export
print.btsbm_strength_prior <- function(x, ...) {
  cat("<latent_strength_prior> Gamma(shape =", x$shape, ", rate =", x$rate, ")\n")
  invisible(x)
}

#' @export
print.btsbm_partition_prior <- function(x, ...) {
  description <- switch(
    x$type,
    gnedin = paste0("Gnedin(gnedin_hyperparameter = ", x$gnedin_hyperparameter, ")"),
    dp = paste0("Dirichlet process(concentration = ", x$concentration, ")"),
    py = paste0("Pitman--Yor(concentration = ", x$concentration,
                ", discount = ", x$discount, ")"),
    dm = paste0("finite Dirichlet(max_clusters = ", x$max_clusters,
                ", concentration = ", x$concentration, ")")
  )
  cat("<clustering_prior>", description, "\n")
  invisible(x)
}

#' Define a Gamma prior for positive latent strengths
#'
#' A Plackett--Luce likelihood only identifies relative strengths. This prior
#' is defined on the positive, unnormalised strengths used by the sampler;
#' `pl_identifiability()` controls the scale used when strengths are reported.
#'
#' @param shape Positive Gamma shape for each latent strength.
#' @param rate Positive Gamma rate for each latent strength.
#'
#' @return An object of class `btsbm_strength_prior`.
#' @export
latent_strength <- function(shape = 1, rate = 1) {
  if (length(shape) != 1L || !is.finite(shape) || shape <= 0 ||
      length(rate) != 1L || !is.finite(rate) || rate <= 0) {
    stop("`shape` and `rate` must be finite positive scalars.", call. = FALSE)
  }
  structure(list(family = "gamma", shape = as.numeric(shape), rate = as.numeric(rate)),
            class = c("btsbm_strength_prior", "btsbm_ability_prior"))
}

#' @rdname latent_strength
#' @export
gamma_ability <- function(shape = 1, rate = 1) {
  latent_strength(shape = shape, rate = rate)
}

#' Choose a prior for a latent partition
#'
#' This is the common prior selector for item blocks in PL--SBM and ranking
#' groups in a PL mixture. `family = "gnedin"` is a flexible finite-cluster
#' prior and is the default. A Dirichlet process and Pitman--Yor process allow
#' a potentially unbounded number of groups. A finite Dirichlet prior fixes an
#' upper bound through `max_clusters`.
#'
#' @param family One of `"gnedin"`, `"dirichlet_process"`, `"pitman_yor"`,
#'   or `"finite_dirichlet"`.
#' @param gnedin_hyperparameter Gnedin hyperparameter in `(0, 1)`. Larger
#'   values favour more occupied groups a priori.
#' @param concentration Concentration for the Dirichlet-process, Pitman--Yor,
#'   or finite-Dirichlet prior.
#' @param discount Pitman--Yor discount in `[0, 1)`.
#' @param max_clusters Required maximum number of groups for
#'   `family = "finite_dirichlet"`.
#'
#' @return An object of class `btsbm_partition_prior`.
#' @export
clustering_prior <- function(
    family = c("gnedin", "dirichlet_process", "pitman_yor", "finite_dirichlet"),
    gnedin_hyperparameter = 0.5, concentration = 1, discount = 0.2,
    max_clusters = NULL) {
  family <- match.arg(family)
  type <- switch(
    family,
    gnedin = "gnedin",
    dirichlet_process = "dp",
    pitman_yor = "py",
    finite_dirichlet = "dm"
  )
  if (type == "gnedin") {
    if (!is.finite(gnedin_hyperparameter) || gnedin_hyperparameter <= 0 || gnedin_hyperparameter >= 1) {
      stop("`gnedin_hyperparameter` must lie in (0, 1).", call. = FALSE)
    }
    values <- list(
      type = type, family = family,
      gnedin_hyperparameter = as.numeric(gnedin_hyperparameter),
      gamma = as.numeric(gnedin_hyperparameter)
    )
  } else if (type == "dp") {
    if (!is.finite(concentration) || concentration <= 0) {
      stop("`concentration` must be positive for a Dirichlet-process prior.", call. = FALSE)
    }
    values <- list(type = type, family = family, concentration = as.numeric(concentration),
                   alpha = as.numeric(concentration), discount = 0, sigma = 0)
  } else if (type == "py") {
    if (!is.finite(discount) || discount < 0 || discount >= 1 ||
        !is.finite(concentration) || concentration <= -discount) {
      stop("Pitman--Yor requires `discount` in [0, 1) and `concentration > -discount`.", call. = FALSE)
    }
    values <- list(type = type, family = family, concentration = as.numeric(concentration),
                   alpha = as.numeric(concentration), discount = as.numeric(discount),
                   sigma = as.numeric(discount))
  } else {
    if (!is.finite(concentration) || concentration <= 0 || length(max_clusters) != 1L ||
        !is.finite(max_clusters) || max_clusters < 1 || max_clusters != floor(max_clusters)) {
      stop("A finite-Dirichlet prior requires positive `concentration` and integer `max_clusters >= 1`.", call. = FALSE)
    }
    values <- list(type = type, family = family, concentration = as.numeric(concentration),
                   beta = as.numeric(concentration), max_clusters = as.integer(max_clusters),
                   K = as.integer(max_clusters))
  }
  structure(values, class = "btsbm_partition_prior")
}

#' Gnedin prior for a latent partition
#'
#' @param gnedin_hyperparameter Gnedin hyperparameter in `(0, 1)`.
#' @return An object of class `btsbm_partition_prior`.
#' @export
gnedin_prior <- function(gnedin_hyperparameter = 0.5) {
  clustering_prior("gnedin", gnedin_hyperparameter = gnedin_hyperparameter)
}

#' Dirichlet-process prior for a latent partition
#'
#' @param concentration Positive concentration parameter.
#' @return An object of class `btsbm_partition_prior`.
#' @export
dirichlet_process_prior <- function(concentration = 1) {
  clustering_prior("dirichlet_process", concentration = concentration)
}

#' Pitman--Yor prior for a latent partition
#'
#' @param concentration Concentration parameter.
#' @param discount Discount parameter in `[0, 1)`.
#' @return An object of class `btsbm_partition_prior`.
#' @export
pitman_yor_prior <- function(concentration = 1, discount = 0.2) {
  clustering_prior("pitman_yor", concentration = concentration, discount = discount)
}

#' Define a finite-Dirichlet prior for a latent partition
#'
#' @param max_clusters Positive integer upper bound on the number of groups.
#' @param concentration Positive symmetric Dirichlet concentration.
#' @return An object of class `btsbm_partition_prior`.
#' @export
finite_partition <- function(max_clusters, concentration = 1) {
  clustering_prior("finite_dirichlet", concentration = concentration,
                   max_clusters = max_clusters)
}

#' Legacy partition-prior interface
#'
#' @param type One of `"gnedin"`, `"dp"`, `"py"`, or `"dm"`.
#' @param gamma Gnedin hyperparameter.
#' @param alpha DP/Pitman--Yor concentration.
#' @param sigma Pitman--Yor discount.
#' @param beta Finite-Dirichlet concentration.
#' @param K Maximum finite-Dirichlet groups.
#' @return An object of class `btsbm_partition_prior`.
#' @export
partition_prior <- function(type = c("gnedin", "dp", "py", "dm"),
                            gamma = 0.5, alpha = 1, sigma = 0.2,
                            beta = 1, K = NULL) {
  type <- match.arg(type)
  switch(
    type,
    gnedin = clustering_prior("gnedin", gnedin_hyperparameter = gamma),
    dp = clustering_prior("dirichlet_process", concentration = alpha),
    py = clustering_prior("pitman_yor", concentration = alpha, discount = sigma),
    dm = clustering_prior("finite_dirichlet", concentration = beta, max_clusters = K)
  )
}

#' Set the identifiable representation for Plackett--Luce abilities
#'
#' @param method Normalisation method.  `"logmean0"` gives each reported
#'   ability vector geometric mean one.
#'
#' @return An object of class `btsbm_pl_identifiability`.
#' @export
pl_identifiability <- function(method = c("logmean0", "none")) {
  structure(list(method = match.arg(method)), class = "btsbm_pl_identifiability")
}

#' Set MCMC controls for a BTSBM fit
#'
#' @param iter Total MCMC iterations.
#' @param warmup Number of warmup iterations.
#' @param seed Optional integer seed.
#' @param store Whether to store the minimal state or latent augmentations.
#' @param init_item_cluster Optional initial item partition.
#' @param init_ranking_cluster Optional initial ranking partition for PL mixture
#'   and PL--LBM models.
#' @param verbose Whether to report progress.
#' @param progress_every Progress frequency when `verbose` is `TRUE`.
#' @param partition_moves Optional result of [partition_moves()].
#'
#' @return An object of class `btsbm_mcmc_control`.
#' @export
mcmc_control <- function(iter = 2000, warmup = floor(iter / 2), seed = NULL,
                         store = c("minimal", "augmented"),
                         init_item_cluster = NULL, init_ranking_cluster = NULL,
                         verbose = FALSE,
                         progress_every = 1000L, partition_moves = NULL) {
  store <- match.arg(store)
  if (length(iter) != 1L || iter < 2 || iter != floor(iter) ||
      length(warmup) != 1L || warmup < 0 || warmup != floor(warmup) || warmup >= iter) {
    stop("Require integer `iter >= 2` and `0 <= warmup < iter`.", call. = FALSE)
  }
  if (!is.null(seed) && (length(seed) != 1L || !is.finite(seed) || seed != floor(seed))) {
    stop("`seed` must be NULL or one integer.", call. = FALSE)
  }
  if (length(progress_every) != 1L || !is.finite(progress_every) ||
      progress_every < 1L || progress_every != floor(progress_every)) {
    stop("`progress_every` must be a positive integer.", call. = FALSE)
  }
  structure(
    list(
      iter = as.integer(iter), warmup = as.integer(warmup), seed = seed,
      store = store, init_item_cluster = init_item_cluster,
      init_ranking_cluster = init_ranking_cluster,
      verbose = isTRUE(verbose), progress_every = as.integer(progress_every),
      partition_moves = partition_moves
    ),
    class = "btsbm_mcmc_control"
  )
}

#' Set partition-move controls
#'
#' @param split_merge Number of generic split--merge moves per iteration.
#' @param ranking_split_merge,item_split_merge PL--LBM moves for the ranking
#'   and item partitions respectively.
#' @param pair_proposal Pair-selection rule for a future split--merge kernel.
#' @param restricted_scans Number of restricted Gibbs scans.
#'
#' @return An object of class `btsbm_partition_moves`.
#' @export
partition_moves <- function(split_merge = 0L, ranking_split_merge = 0L,
                            item_split_merge = 0L,
                            pair_proposal = c("uniform", "balanced", "informed"),
                            restricted_scans = 1L) {
  values <- c(split_merge, ranking_split_merge, item_split_merge, restricted_scans)
  if (any(!is.finite(values)) || any(values < 0) || any(values != floor(values)) || restricted_scans < 1L) {
    stop("Move counts must be non-negative integers and `restricted_scans >= 1`.", call. = FALSE)
  }
  structure(
    list(
      split_merge = as.integer(split_merge),
      ranking_split_merge = as.integer(ranking_split_merge),
      item_split_merge = as.integer(item_split_merge),
      pair_proposal = match.arg(pair_proposal),
      restricted_scans = as.integer(restricted_scans)
    ),
    class = "btsbm_partition_moves"
  )
}

.default_latent_strength <- function(shape = 1, rate = 1) {
  latent_strength(shape = shape, rate = rate)
}

.new_btsbm_model <- function(likelihood, structure, latent_strength, item_clustering = NULL,
                              ranking_clustering = NULL, identifiability = NULL) {
  if (!inherits(latent_strength, "btsbm_strength_prior")) {
    stop("`latent_strength` must be created by `latent_strength()`.", call. = FALSE)
  }
  if (!is.null(item_clustering) && !inherits(item_clustering, "btsbm_partition_prior")) {
    stop("`item_clustering` must be created by `clustering_prior()` or a named prior helper.", call. = FALSE)
  }
  if (!is.null(ranking_clustering) && !inherits(ranking_clustering, "btsbm_partition_prior")) {
    stop("`ranking_clustering` must be created by `clustering_prior()` or a named prior helper.", call. = FALSE)
  }
  structure(
    list(
      likelihood = likelihood, structure = structure,
      latent_strength = latent_strength,
      item_clustering = item_clustering,
      ranking_clustering = ranking_clustering,
      identifiability = identifiability,
      # Compatibility fields for code written before the public API rename.
      ability_prior = latent_strength, item_prior = item_clustering,
      ranking_prior = ranking_clustering
    ),
    class = "btsbm_model"
  )
}

#' Define a simple Bradley--Terry model
#' @param latent_strength Gamma prior for positive item strengths.
#' @param ability_prior Deprecated alias for `latent_strength`.
#' @return A `btsbm_model`.
#' @export
bt_model <- function(latent_strength = .default_latent_strength(shape = 0.01, rate = 0.1),
                     ability_prior = NULL) {
  if (!is.null(ability_prior)) latent_strength <- ability_prior
  .new_btsbm_model("bt", "none", latent_strength)
}

#' Define a Bradley--Terry stochastic block model
#' @param latent_strength Gamma prior for positive block strengths.
#' @param item_clustering Prior on latent item blocks.
#' @param ability_prior Deprecated alias for `latent_strength`.
#' @param item_prior Deprecated alias for `item_clustering`.
#' @return A `btsbm_model`.
#' @export
bt_sbm_model <- function(latent_strength = .default_latent_strength(shape = 4, rate = exp(digamma(4))),
                         item_clustering = gnedin_prior(gnedin_hyperparameter = 0.5),
                         ability_prior = NULL, item_prior = NULL) {
  if (!is.null(ability_prior)) latent_strength <- ability_prior
  if (!is.null(item_prior)) item_clustering <- item_prior
  .new_btsbm_model("bt", "item_sbm", latent_strength,
                    item_clustering = item_clustering)
}

#' Define a simple Plackett--Luce model
#' @param latent_strength Gamma prior for positive item strengths.
#' @param identifiability PL reporting constraint.
#' @param ability_prior Deprecated alias for `latent_strength`.
#' @return A `btsbm_model`.
#' @export
pl_model <- function(latent_strength = .default_latent_strength(),
                     identifiability = pl_identifiability(), ability_prior = NULL) {
  if (!is.null(ability_prior)) latent_strength <- ability_prior
  .new_btsbm_model("pl", "none", latent_strength, identifiability = identifiability)
}

#' Define a Plackett--Luce stochastic block model
#' @param latent_strength Gamma prior for positive block strengths.
#' @param item_clustering Prior on latent item blocks. Use [gnedin_prior()],
#'   [dirichlet_process_prior()], [pitman_yor_prior()], or [finite_partition()].
#' @param identifiability PL reporting constraint.
#' @param ability_prior Deprecated alias for `latent_strength`.
#' @param item_prior Deprecated alias for `item_clustering`.
#' @return A `btsbm_model`.
#' @export
pl_sbm_model <- function(latent_strength = .default_latent_strength(),
                         item_clustering = gnedin_prior(gnedin_hyperparameter = 0.5),
                         identifiability = pl_identifiability(), ability_prior = NULL,
                         item_prior = NULL) {
  if (!is.null(ability_prior)) latent_strength <- ability_prior
  if (!is.null(item_prior)) item_clustering <- item_prior
  .new_btsbm_model("pl", "item_sbm", latent_strength, item_clustering = item_clustering,
                    identifiability = identifiability)
}

#' Define a Plackett--Luce ranking-mixture model
#'
#' Ranking rows are assigned to latent preference groups. Each group has a
#' separate vector of item strengths.
#'
#' @param latent_strength Gamma prior for positive component strengths.
#' @param ranking_clustering Prior on latent ranking groups.
#' @param identifiability PL reporting constraint.
#' @param ability_prior Deprecated alias for `latent_strength`.
#' @param ranking_prior Deprecated alias for `ranking_clustering`.
#' @return A `btsbm_model`.
#' @export
pl_mixture_model <- function(latent_strength = .default_latent_strength(),
                             ranking_clustering = gnedin_prior(gnedin_hyperparameter = 0.5),
                             identifiability = pl_identifiability(), ability_prior = NULL,
                             ranking_prior = NULL) {
  if (!is.null(ability_prior)) latent_strength <- ability_prior
  if (!is.null(ranking_prior)) ranking_clustering <- ranking_prior
  .new_btsbm_model("pl", "ranking_mixture", latent_strength,
                    ranking_clustering = ranking_clustering,
                    identifiability = identifiability)
}

#' Define a Plackett--Luce latent block model
#'
#' The PL--LBM jointly partitions ranking rows and items. Each pair of latent
#' groups has its own positive strength.
#'
#' @param latent_strength Gamma prior for positive cell strengths.
#' @param ranking_clustering Prior on latent ranking groups.
#' @param item_clustering Prior on latent item blocks.
#' @param identifiability PL reporting constraint.
#' @param ability_prior Deprecated alias for `latent_strength`.
#' @param ranking_prior Deprecated alias for `ranking_clustering`.
#' @param item_prior Deprecated alias for `item_clustering`.
#' @return A `btsbm_model`.
#' @export
pl_lbm_model <- function(latent_strength = .default_latent_strength(),
                         ranking_clustering = gnedin_prior(gnedin_hyperparameter = 0.5),
                         item_clustering = gnedin_prior(gnedin_hyperparameter = 0.5),
                         identifiability = pl_identifiability(), ability_prior = NULL,
                         ranking_prior = NULL, item_prior = NULL) {
  if (!is.null(ability_prior)) latent_strength <- ability_prior
  if (!is.null(ranking_prior)) ranking_clustering <- ranking_prior
  if (!is.null(item_prior)) item_clustering <- item_prior
  .new_btsbm_model("pl", "lbm", latent_strength,
                    item_clustering = item_clustering,
                    ranking_clustering = ranking_clustering,
                    identifiability = identifiability)
}

#' @export
print.btsbm_model <- function(x, ...) {
  likelihood <- switch(x$likelihood, bt = "Bradley--Terry", pl = "Plackett--Luce")
  structure <- switch(
    x$structure,
    none = "item-specific strengths",
    item_sbm = "latent item clustering",
    ranking_mixture = "latent ranking clustering",
    lbm = "joint ranking and item clustering",
    x$structure
  )
  cat("<btsbm_model>", likelihood, "with", structure, "\n")
  cat("  latent strength: Gamma(shape =", x$latent_strength$shape,
      ", rate =", x$latent_strength$rate, ")\n")
  if (!is.null(x$item_clustering)) {
    cat("  item clustering:")
    print(x$item_clustering)
  }
  if (!is.null(x$ranking_clustering)) {
    cat("  ranking clustering:")
    print(x$ranking_clustering)
  }
  invisible(x)
}
