#' Gibbs sampler for the Bradley–Terry Stochastic Block Model (BT–SBM)
#'
#' Runs a Gibbs sampler for the BT–SBM with optional DP/PY/DM/GN priors on
#' the partition. Returns raw draws plus minimal summaries; `n_ij` is
#' computed internally as `w_ij + t(w_ij)`.
#'
#' Row names of \code{w_ij} (if present) are propagated to item-level outputs;
#' otherwise items are named \code{Item_1, ..., Item_n}.
#'
#' @param w_ij Integer or numeric square matrix \eqn{n \times n} of directed wins
#'   (i over j). Must be nonnegative with zero diagonal. The function builds
#'   \eqn{n_{ij} = w_{ij} + w_{ji}} internally.
#' @param a Positive shape parameter for the Gamma prior
#'   \eqn{\lambda_k \sim \mathrm{Gamma}(a,b)}. The algorithm uses
#'   \eqn{b = \exp(\psi(a))} so that \eqn{\mathbb{E}[\log \lambda_k] = 0} a priori.
#' @param prior Character scalar, one of \code{"DP"}, \code{"PY"}, \code{"DM"}, \code{"GN"}.
#' @param alpha_PY,sigma_PY Hyperparameters for Pitman–Yor / Dirichlet Process.
#'   For \code{prior="DP"} use \code{alpha_PY} (with \code{sigma_PY} ignored).
#'   For \code{prior="PY"} use both \code{alpha_PY} and \code{sigma_PY \in (0,1)}.
#' @param beta_DM,K_DM Hyperparameters for the finite Dirichlet–Multinomial prior.
#'   \code{K_DM} is the maximum number of allowed clusters.
#' @param gamma_GN Hyperparameter for the Gnedin prior.
#' @param T_iter,T_burn Integers: total iterations and burn-in. Must satisfy \code{T_burn < T_iter}.
#' @param init_x Optional integer vector of length \code{n} with initial labels (1-based).
#' @param store_z Logical; if \code{TRUE}, store latent \code{Z} draws (heavy).
#' @param verbose Logical; if \code{TRUE}, prints progress every 1000 iterations.
#'
#' @return A \code{list} with:
#' \itemize{
#'   \item \code{x_samples}: integer matrix \eqn{S \times n} of raw labels (\eqn{S = T_{\mathrm{iter}}-T_{\mathrm{burn}}}).
#'   \item \code{lambda_samples}: list of length \eqn{S}; each element is a numeric vector
#'         of length \eqn{L_{\mathrm{cap}}} for that draw, with \code{NA} at empty labels.
#'   \item \code{K_per_iter}: integer vector length \eqn{S} (occupied cluster count per saved draw).
#'   \item \code{L_cap_per_iter}: integer vector length \eqn{S} (label capacity trace).
#'   \item \code{z_samples}: if \code{store_z=TRUE}, a numeric array \eqn{S \times n \times n}; otherwise \code{NULL}.
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 6L
#' w <- matrix(0L, n, n)
#' w[lower.tri(w)] <- rpois(sum(lower.tri(w)), 2)
#' diag(w) <- 0
#' rownames(w) <- colnames(w) <- paste0("P", seq_len(n))
#' fit <- gibbs_bt_sbm(
#'   w_ij = w, prior = "GN", gamma_GN = 0.5,
#'   T_iter = 200, T_burn = 100, verbose = FALSE
#' )
#' }
#' @importFrom mcclust comp.psm
#' @importFrom mcclust.ext minVI
#' @export
gibbs_bt_sbm <- function(
    w_ij,
    a = 4,
    prior = c("DP", "PY", "DM", "GN"),
    alpha_PY = NA_real_,
    sigma_PY = NA_real_,
    beta_DM  = NA_real_,
    K_DM    = NA_integer_,
    gamma_GN = NA_real_,
    T_iter = 2000, T_burn = 1000,
    init_x = NULL,
    store_z = FALSE,
    verbose = TRUE
) {
  ## ---- input checks & derived quantities ---------------------------------
  stopifnot(is.matrix(w_ij))
  n <- nrow(w_ij)
  stopifnot(n == ncol(w_ij))
  if (any(w_ij < 0)) stop("w_ij: negative entries not allowed.")
  if (any(diag(w_ij) != 0)) stop("w_ij: diagonal must be zero.")
  n_ij <- w_ij + t(w_ij)
  if (!isTRUE(all.equal(n_ij, t(n_ij)))) stop("n_ij must be symmetric.")
  if (T_burn >= T_iter) stop("Require T_burn < T_iter.")
  if (a <= 0) stop("Gamma(a,b): need a > 0.")

  prior <- match.arg(prior)
  if (prior == "DP") {
    if (!is.finite(alpha_PY)) stop("DP requires alpha_PY.")
  } else if (prior == "PY") {
    if (!is.finite(alpha_PY) || !is.finite(sigma_PY) || sigma_PY <= 0 || sigma_PY >= 1)
      stop("PY requires alpha_PY and sigma_PY in (0,1).")
  } else if (prior == "DM") {
    if (!is.finite(beta_DM) || !is.finite(K_DM) || K_DM < 1)
      stop("DM requires beta_DM > 0 and K_DM >= 1.")
  } else if (prior == "GN") {
    if (!is.finite(gamma_GN)) stop("GN requires gamma_GN.")
  }

  urn_fun <- switch(
    prior,
    DP = function(v) urn_DP(v, alpha_PY),
    PY = function(v) urn_PY(v, alpha_PY, sigma_PY),
    DM = function(v) urn_DM(v, beta_DM, K_DM),
    GN = function(v) urn_GN(v, gamma_GN)
  )

  ## precompute row wins
  w_i <- rowSums(w_ij)

  ## ---- init labels and rates ---------------------------------------------
  L_cap <- if (is.finite(K_DM)) min(n, K_DM) else n  # initial capacity
  if (is.null(init_x)) {
    x_curr <- sample.int(L_cap, n, replace = TRUE)
  } else {
    stopifnot(length(init_x) == n)
    x_curr <- as.integer(init_x)
    L_cap <- max(L_cap, max(x_curr))
  }
  a_curr <- a
  b_eff  <- as.numeric(exp(digamma(a_curr)))  # b := exp(psi(a))

  lambda_curr <- rep(NA_real_, L_cap)
  csize0 <- tabulate(x_curr, nbins = L_cap)
  occ0   <- which(csize0 > 0L)
  if (length(occ0)) {
    lambda_curr[occ0] <- rgamma(length(occ0), shape = a_curr, rate = b_eff)
  }

  ## latent Z
  Z_curr <- matrix(0, n, n)
  S <- T_iter - T_burn
  z_store <- if (store_z) array(0, dim = c(S, n, n)) else NULL

  ## storage
  x_samples    <- matrix(NA_integer_, S, n)
  lambda_list  <- vector("list", S)   # ragged storage for lambdas
  L_cap_trace  <- integer(S)
  K_trace      <- integer(S)

  ## ---- main loop ----------------------------------------------------------
  save_i <- 0L
  for (iter in seq_len(T_iter)) {
    ## 1) Update Z and row sums (C++ helper expects length(lambda) >= max(x))
    outZ <- updateZ_and_rowsums(
      n_ij  = structure(n_ij, .Dim = dim(n_ij)),
      x     = as.integer(x_curr),
      lambda= lambda_curr
    )
    Z_curr    <- outZ$Z
    Z_row_sum <- outZ$rowSumsZ

    ## 2) Single-site updates for x_i
    for (i in seq_len(n)) {
      csize_minus <- tabulate(x_curr[-i], nbins = L_cap)
      occupied <- which(csize_minus > 0L)
      H <- length(occupied)
      v_minus <- csize_minus[occupied]

      # prior masses for existing and new
      prior_vec <- if (H > 0L) urn_fun(v_minus) else c(0, urn_fun(integer(0))[1L])
      new_weight <- prior_vec[H + 1L]
      existing_weights <- if (H > 0L) prior_vec[seq_len(H)] else numeric(0)

      Zi <- Z_row_sum[i]
      wi <- w_i[i]

      logp <- rep(-Inf, H + 1L)

      if (H > 0L) {
        lam_h  <- lambda_curr[occupied]
        llh    <- wi * .safe_log(lam_h) - lam_h * Zi
        lprior <- .safe_log(existing_weights)
        logp[seq_len(H)] <- lprior + llh
      }

      # "new" cluster option (anchored integral)
      if (new_weight > 0) {
        logp[H + 1L] <- .safe_log(new_weight) +
          .new_cluster_integral_log_anchored(wi, Zi, a_curr, b_eff)
      }

      choice <- .sample_from_logweights(logp)

      if (choice <= H) {
        x_curr[i] <- occupied[choice]
      } else {
        # create / reuse empty label id
        csize_all <- tabulate(x_curr[-i], nbins = L_cap)
        empties <- which(csize_all == 0L)

        if (length(empties) == 0L) {
          # grow capacity if allowed
          if (is.finite(K_DM) && L_cap >= K_DM) {
            # cannot create new; fall back to a random existing cluster
            x_curr[i] <- occupied[ if (H > 0L) sample.int(H, 1L) else 1L ]
          } else {
            L_cap <- L_cap + 1L
            lambda_curr <- c(lambda_curr, NA_real_)
            x_curr[i] <- L_cap
          }
        } else {
          x_curr[i] <- empties[1L]
        }

        # initialize lambda for the (newly) used label
        lbl <- x_curr[i]
        shape_k <- a_curr + wi
        rate_k  <- b_eff + Zi
        lambda_curr[lbl] <- rgamma(1L, shape = shape_k, rate = rate_k)
      }
    }

    ## 3) Update lambda_k for occupied clusters (conjugate)
    csize <- tabulate(x_curr, nbins = L_cap)
    occ <- which(csize > 0L)
    if (length(occ)) {
      for (k in occ) {
        members <- which(x_curr == k)
        shape_k <- a_curr + sum(w_i[members])
        rate_k  <- b_eff + sum(Z_row_sum[members])
        lambda_curr[k] <- rgamma(1L, shape = shape_k, rate = rate_k)
      }
    }
    lambda_curr[csize == 0L] <- NA_real_

    ## Numerically stable global rescaling so that E(log lambda) ~ 0
    if (length(occ) > 0L) {
      g_mean <- exp(mean(log(lambda_curr[occ])))
      lambda_curr[occ] <- lambda_curr[occ] / g_mean
    }

    ## 4) Store
    if (iter > T_burn) {
      save_i <- save_i + 1L
      x_samples[save_i, ]      <- x_curr
      lambda_list[[save_i]]    <- lambda_curr
      K_trace[save_i]          <- length(occ)
      if (store_z) z_store[save_i, , ] <- Z_curr
    }

    if (verbose && iter %% 1000L == 0L) {
      cat("iter", iter, "occupied =", length(unique(x_curr)), "\n")
    }
  }

  list(
    x_samples        = x_samples,
    lambda_samples   = lambda_list,      # ragged storage
    K_per_iter       = K_trace,
    z_samples        = if (store_z) z_store else NULL
  )
}
