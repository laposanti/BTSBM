#' ATP Paired-Comparison Panels, 2000–2022
#'
#' A yearly panel of head-to-head counts for the top 105 ATP players, suitable
#' for Bradley–Terry/Stochastic Block Model analyses. For each season, the data
#' include (i) the number of wins by player \eqn{i} against player \eqn{j}
#' and (ii) the number of matches played between \eqn{i} and \eqn{j}, along
#' with a per-season player metadata tibble.
#'
#' @format A named list of length 23, with elements `"2000"`, `"2001"`, …, `"2022"`.
#'   Each yearly element is a list of length 3:
#'   \describe{
#'     \item{`Y_ij`}{A numeric matrix \eqn{105 \times 105}. Entry \code{Y_ij[i, j]}
#'       is the count of matches in which player \eqn{i} defeated player \eqn{j}
#'       in that calendar year (nonnegative integer; diagonal is zero).}
#'     \item{`N_ij`}{A numeric matrix \eqn{105 \times 105}. Entry \code{N_ij[i, j]}
#'       is the total number of matches between players \eqn{i} and \eqn{j} that year
#'       (nonnegative integer; symmetric by construction; diagonal is zero).}
#'     \item{`players_df`}{A tibble/data frame with 105 rows and 7 columns
#'       describing the player index used in the matrices for that year:
#'       \describe{
#'         \item{`player`}{Integer player identifier (row/column index used in
#'           `Y_ij` and `N_ij`).}
#'         \item{`worst_rank`}{Worst (numerically largest) ATP ranking attained
#'           by the player during the year.}
#'         \item{`median_rank`}{Median ATP ranking across the player's ranking
#'           snapshots in that year.}
#'         \item{`last_rank`}{ATP ranking at the last snapshot available in the
#'           year (e.g., year-end ranking).}
#'         \item{`age_year`}{Approximate age (in years) for that season (e.g.,
#'           at mid-season).}
#'         \item{`ht_year`}{Player height in centimeters (season-level value).}
#'         \item{`player_slug`}{Character identifier (URL-safe or underscored
#'           name) for the player.}
#'       }}
#'   }
#'
#' @details
#' The player ordering in `players_df` defines the row/column indexing of
#' `Y_ij` and `N_ij` for the corresponding year. The diagonal entries of both
#' matrices are zero by definition. In typical usage for Bradley–Terry-type
#' models, one can treat \code{Y_ij[i, j]} as the number of “successes” for
#' \eqn{i} vs. \eqn{j}, with the binomial denominator \code{N_ij[i, j]}.
#'
#' @note
#' The matrices may be sparse for many player pairs. Ensure any model code
#' guards against divisions by zero when \code{N_ij[i, j] = 0}.
#'
#' @source
#' Aggregated by the package author from public ATP results (e.g., the
#' \emph{tennis\_atp} datasets by Jeff Sackmann) and internal preprocessing.
#' See the package vignette for provenance and cleaning steps.
#'
#' @examples
#' data(ATP_2000_2022)
#' names(ATP_2000_2022)
#' year <- "2000"
#' str(ATP_2000_2022[[year]])
#'
#' # Player i's total wins that year:
#' i <- 1
#' sum(ATP_2000_2022[[year]]$Y_ij[i, ], na.rm = TRUE)
#'
#' # Total matches between i and j:
#' j <- 2
#' ATP_2000_2022[[year]]$N_ij[i, j]
#'
#' # Join player metadata to indices used in the matrices:
#' head(ATP_2000_2022[[year]]$players_df)
#'
"ATP_2000_2022"
