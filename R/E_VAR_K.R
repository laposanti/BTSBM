#' Expected number of clusters under the Gnedin prior
#'
#' Compute the closed-form mean of \eqn{K_n} (number of clusters) under the
#' Gnedin finite-type model (\eqn{\sigma=-1}).
#'
#' For sample size \eqn{n \ge 1} and parameter \eqn{\gamma \in (0,1)}, the mean is
#' \deqn{
#'   \mathbb{E}[K_n] \;=\; \frac{\Gamma(n+1)\,\Gamma(1+\gamma)}{\Gamma(n+\gamma)}.
#' }
#' The implementation uses \code{lgamma} for numerical stability and is vectorized
#' over \code{n} and \code{gamma} (with standard R recycling rules).
#'
#' @param n Integer vector of sample sizes (each \eqn{\ge 1}).
#' @param gamma Numeric vector of Gnedin parameters, each in \eqn{(0,1)}.
#'
#' @return A numeric vector with \eqn{\mathbb{E}[K_n]} for each pair \code{(n, gamma)}.
#'
#' @details
#' The formula follows directly from standard Gibbs–type manipulations using
#' factorial moments and the Chu–Vandermonde identity specialized to \eqn{\sigma=-1}.
#'
#' @references
#' Gnedin, A. (2010). A species sampling model with finitely many types.
#' \emph{Electronic Communications in Probability}, 15, 79–88.
#'
#' Favaro, S., Lijoi, A., & Prünster, I. (2013).
#' Extending the class of Gibbs-type priors: theoretical properties and new examples.
#' \emph{Annals of Applied Probability}, 23(4), 1729–1754.
#'
#' Pitman, J. (2006). \emph{Combinatorial Stochastic Processes}. Springer.
#'
#' @seealso [gnedin_K_var()]
#'
#' @examples
#' # Scalar inputs
#' gnedin_K_mean(105, 0.8)
#'
#' # Vectorized over gamma
#' gnedin_K_mean(50, c(0.3, 0.5, 0.8))
#'
#' # Vectorized over n and gamma (recycling rules)
#' gnedin_K_mean(c(20, 50, 100), 0.5)
#'
#' @export
gnedin_K_mean <- function(n, gamma) {
  .gn_check_inputs(n, gamma)
  # E[K_n] = Gamma(n+1) * Gamma(1+gamma) / Gamma(n+gamma)
  exp(lgamma(n + 1) + lgamma(1 + gamma) - lgamma(n + gamma))
}

#' Variance of the number of clusters under the Gnedin prior
#'
#' Compute the closed-form variance of \eqn{K_n} under the Gnedin finite-type model
#' (\eqn{\sigma=-1}).
#'
#' Using the factorial–ordinary moment relation
#' \eqn{\mathbb{E}[K_n^2]=\mathbb{E}[K_n(K_n-1)]+\mathbb{E}[K_n]}, the variance is
#' \deqn{
#'   \mathrm{Var}(K_n) \;=\; \mathbb{E}[K_n(K_n-1)] + \mathbb{E}[K_n] - \{\mathbb{E}[K_n]\}^2,
#' }
#' where
#' \deqn{
#'   \mathbb{E}[K_n(K_n-1)]
#'   \;=\; n(n-1)(1-\gamma)\,\frac{\Gamma(n)\,\Gamma(1+\gamma)}{\Gamma(n+\gamma)}.
#' }
#' The implementation uses \code{lgamma} for numerical stability and is vectorized
#' over \code{n} and \code{gamma} (with standard R recycling rules).
#'
#' @param n Integer vector of sample sizes (each \eqn{\ge 1}).
#' @param gamma Numeric vector of Gnedin parameters, each in \eqn{(0,1)}.
#'
#' @return A numeric vector with \eqn{\mathrm{Var}(K_n)} for each pair \code{(n, gamma)}.
#'
#' @references
#' Gnedin, A. (2010). A species sampling model with finitely many types.
#' \emph{Electronic Communications in Probability}, 15, 79–88.
#'
#' Favaro, S., Lijoi, A., & Prünster, I. (2013).
#' Extending the class of Gibbs-type priors: theoretical properties and new examples.
#' \emph{Annals of Applied Probability}, 23(4), 1729–1754.
#'
#' Pitman, J. (2006). \emph{Combinatorial Stochastic Processes}. Springer.
#'
#' @seealso [gnedin_K_mean()]
#'
#' @examples
#' # Scalar inputs
#' gnedin_K_var(105, 0.8)
#'
#' # Consistency check: variance is nonnegative
#' all(gnedin_K_var(20, c(0.3, 0.5, 0.8)) >= 0)
#'
#' @export
gnedin_K_var <- function(n, gamma) {
  .gn_check_inputs(n, gamma)
  EK  <- gnedin_K_mean(n, gamma)
  EK2f <- n * (n - 1) * (1 - gamma) * exp(lgamma(n) + lgamma(1 + gamma) - lgamma(n + gamma))
  EK2f + EK - EK^2
}

# ---- internal helpers --------------------------------------------------------

# Validate inputs for Gnedin functions
.gn_check_inputs <- function(n, gamma) {
  if (length(n) == 0L || length(gamma) == 0L) {
    stop("'n' and 'gamma' must have positive length.", call. = FALSE)
  }
  if (!is.numeric(n) || any(!is.finite(n))) {
    stop("'n' must be numeric and finite.", call. = FALSE)
  }
  if (any(n < 1) || any(abs(n - round(n)) > .Machine$double.eps^0.5)) {
    stop("'n' must contain integers >= 1.", call. = FALSE)
  }
  if (!is.numeric(gamma) || any(!is.finite(gamma))) {
    stop("'gamma' must be numeric and finite.", call. = FALSE)
  }
  if (any(gamma <= 0 | gamma >= 1)) {
    stop("'gamma' must lie in (0, 1) for the Gnedin model.", call. = FALSE)
  }
  invisible(TRUE)
}
