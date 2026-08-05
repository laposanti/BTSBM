#' Simulate strict Plackett--Luce rankings from a weight matrix
#' @param strength Matrix with one positive latent-strength vector per ranking row.
#' @param rank_length Number of retained ranking positions.
#' @return Integer ranking matrix.
#' @keywords internal
.simulate_pl_rankings <- function(strength, rank_length = ncol(strength)) {
  strength <- as.matrix(strength)
  if (any(!is.finite(strength)) || any(strength <= 0)) {
    stop("All ranking strengths must be finite and positive.", call. = FALSE)
  }
  rank_length <- as.integer(rank_length)
  if (rank_length < 1L || rank_length > ncol(strength)) {
    stop("`rank_length` must lie in 1:ncol(strength).", call. = FALSE)
  }
  rankings <- matrix(NA_integer_, nrow(strength), rank_length)
  for (ell in seq_len(nrow(strength))) {
    race <- stats::rexp(ncol(strength), rate = strength[ell, ])
    rankings[ell, ] <- order(race)[seq_len(rank_length)]
  }
  rankings
}

#' Sushi-inspired toy ranking data
#'
#' A reproducible, synthetic teaching dataset with two preference profiles over
#' ten sushi types. It is designed to demonstrate `as_rankings()`, `pl_model()`,
#' PL mixtures, and PL--LBM data shapes without redistributing respondent-level
#' data from the Sushi Preference Dataset.
#'
#' @param n_rankings Number of synthetic rankings.
#' @param rank_length Number of reported positions per ranking.
#' @param seed Integer seed.
#'
#' @return A `btsbm_ranking_data` object.  Its `"truth"` attribute contains
#'   synthetic ranking-cluster labels and the generating ability matrix.
#' @references Kamishima, T. (2003). Nantonac collaborative filtering:
#'   Recommendation based on order responses. *Proceedings of the Ninth ACM
#'   SIGKDD International Conference on Knowledge Discovery and Data Mining*.
#'   The original Sushi Preference Dataset must be obtained from its authors;
#'   its licence does not allow redistribution.
#' @export
sushi_toy <- function(n_rankings = 40L, rank_length = 5L, seed = 2026L) {
  if (length(n_rankings) != 1L || n_rankings < 2L || n_rankings != floor(n_rankings)) {
    stop("`n_rankings` must be an integer of at least two.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  items <- c(
    "egg", "shrimp", "tuna", "squid", "sea_urchin",
    "salmon_roe", "fatty_tuna", "eel", "cucumber_roll", "inari"
  )
  profiles <- rbind(
    light = c(4.5, 4.0, 3.5, 2.8, 1.6, 2.2, 1.3, 1.1, 3.3, 3.0),
    rich = c(1.5, 2.0, 4.0, 3.2, 4.5, 4.0, 5.2, 4.8, 1.4, 1.2)
  )
  ranking_cluster <- rep(seq_len(nrow(profiles)), length.out = n_rankings)
  ranking_cluster <- sample(ranking_cluster)
  strength <- profiles[ranking_cluster, , drop = FALSE]
  data <- as_rankings(
    .simulate_pl_rankings(strength, rank_length = rank_length),
    item_count = length(items),
    item_labels = items
  )
  attr(data, "truth") <- list(
    ranking_cluster = ranking_cluster,
    latent_strength = profiles,
    ability = profiles,
    source = "synthetic sushi-inspired teaching data; not the Sushi Preference Dataset"
  )
  data
}

#' Eurovision-inspired toy ranking data
#'
#' A synthetic preference-ranking dataset for a song-contest setting. Two
#' latent voting profiles produce different rankings of eight fictional acts.
#' It is designed to demonstrate [pl_mixture_model()] without redistributing
#' data from a real contest.
#'
#' @param n_rankings Number of synthetic jury or viewer rankings.
#' @param rank_length Number of reported positions per ranking.
#' @param seed Integer seed.
#'
#' @return A btsbm_ranking_data object. Its truth attribute contains the
#'   synthetic voter-group labels and generating latent strengths.
#' @export
eurovision_toy <- function(n_rankings = 40L, rank_length = 5L, seed = 2030L) {
  if (length(n_rankings) != 1L || n_rankings < 2L || n_rankings != floor(n_rankings)) {
    stop("'n_rankings' must be an integer of at least two.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  acts <- c(
    "Aurora", "Baltica", "Caspia", "Danubia",
    "Etruria", "Fjordland", "Galicia", "Helvetia"
  )
  latent_strength <- rbind(
    televote = c(5.0, 1.3, 3.8, 1.2, 4.4, 2.4, 2.0, 1.1),
    jury = c(2.0, 4.8, 1.4, 4.3, 1.5, 3.6, 2.7, 3.0)
  )
  ranking_cluster <- sample(rep(seq_len(nrow(latent_strength)), length.out = n_rankings))
  rankings <- .simulate_pl_rankings(
    latent_strength[ranking_cluster, , drop = FALSE], rank_length = rank_length
  )
  data <- as_rankings(rankings, item_count = length(acts), item_labels = acts)
  attr(data, "truth") <- list(
    ranking_cluster = ranking_cluster,
    latent_strength = latent_strength,
    source = "synthetic Eurovision-inspired teaching data; not contest voting data"
  )
  data
}

#' Tennis-league toy paired-comparison data
#'
#' A synthetic, tiered round-robin tournament for demonstrating the
#' Bradley--Terry and BT--SBM interfaces. It contains no real player results.
#'
#' @param matches_per_pair Positive number of matches played by each pair.
#' @param seed Integer seed.
#'
#' @return A btsbm_pairwise_data object. Its truth attribute contains the
#'   generating player strengths and latent tiers.
#' @export
tennis_toy <- function(matches_per_pair = 8L, seed = 2031L) {
  if (length(matches_per_pair) != 1L || matches_per_pair < 1L ||
      matches_per_pair != floor(matches_per_pair)) {
    stop("'matches_per_pair' must be a positive integer.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  player_labels <- c("Ace", "Baseline", "Crosscourt", "Dropshot", "Elite", "Forehand")
  latent_strength <- c(5.5, 4.7, 2.5, 2.1, 1.2, 0.9)
  player_tier <- c(1L, 1L, 2L, 2L, 3L, 3L)
  wins <- matrix(0, nrow = length(player_labels), ncol = length(player_labels),
                 dimnames = list(player_labels, player_labels))
  for (first_player in seq_len(nrow(wins) - 1L)) {
    for (second_player in (first_player + 1L):nrow(wins)) {
      win_probability <- latent_strength[first_player] /
        (latent_strength[first_player] + latent_strength[second_player])
      first_wins <- stats::rbinom(1L, size = matches_per_pair, prob = win_probability)
      wins[first_player, second_player] <- first_wins
      wins[second_player, first_player] <- matches_per_pair - first_wins
    }
  }
  data <- as_bt_data(wins, item_labels = player_labels)
  attr(data, "truth") <- list(
    latent_strength = stats::setNames(latent_strength, player_labels),
    item_cluster = player_tier,
    source = "synthetic tennis-league teaching data; not professional match data"
  )
  data
}

#' Structured Plackett--Luce latent-block-model toy data
#'
#' Generates rankings from a two-way latent block model with known ranking-row
#' clusters, item blocks, and a `C × K` latent-strength matrix. It is a compact test
#' and vignette fixture for the future PL--LBM sampler.
#'
#' @param n_rankings Number of ranking rows.
#' @param rank_length Number of reported positions.
#' @param seed Integer seed.
#'
#' @return A `btsbm_ranking_data` object with a `"truth"` attribute containing
#'   the generating ranking and item partitions.
#' @export
pl_lbm_toy <- function(n_rankings = 30L, rank_length = 4L, seed = 2027L) {
  if (length(n_rankings) != 1L || n_rankings < 2L || n_rankings != floor(n_rankings)) {
    stop("`n_rankings` must be an integer of at least two.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  item_labels <- paste0("Item_", seq_len(6L))
  item_cluster <- c(1L, 1L, 2L, 2L, 3L, 3L)
  latent_strength <- rbind(
    c(5.0, 2.0, 0.8),
    c(0.8, 5.0, 2.0)
  )
  ranking_cluster <- sample(rep(1:2, length.out = n_rankings))
  strength <- latent_strength[ranking_cluster, item_cluster, drop = FALSE]
  data <- as_rankings(
    .simulate_pl_rankings(strength, rank_length = rank_length),
    item_count = length(item_labels), item_labels = item_labels
  )
  attr(data, "truth") <- list(
    ranking_cluster = ranking_cluster,
    item_cluster = item_cluster,
    latent_strength = latent_strength,
    ability = latent_strength
  )
  data
}
