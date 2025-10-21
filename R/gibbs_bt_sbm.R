#' Gibbs sampler for a Bradley–Terry Stochastic Block Model (BT–SBM)
#'
#' Runs the Gibbs sampler for the BT–SBM.
#' It then returns relabeled samples and computes minimal posterior summaries.
#' Row names of \code{w_ij} (if present) are propagated to all item-level outputs;
#' otherwise items are named \code{Item_1, ..., Item_N}.
#'
#' @param w_ij Integer or numeric square matrix \eqn{n \times n} of directed wins
#'   (e.g., \code{w_ij[i,j]} is wins of \eqn{i} over \eqn{j}). Must have nonnegative
#'   entries and \code{n_ij <= w_ij + t(w_ij)} elementwise.
#' @param a positive parameter for the Gamma prior \eqn{\lambda_k \sim \mathrm{Gamma}(a,b)}.
#'   The algorithm uses \eqn{b=\exp(\psi(a))} to align the prior with the likelihood expectation on log \eqn{\lambda_k}
#' @param prior Character scalar, one of \code{"DP"}, \code{"PY"}, \code{"DM"}, \code{"GN"}.
#' @param alpha_PY,sigma_PY Hyperparameters for Pitman–Yor/Dirichlet Process.
#'   For \code{prior="DP"} use \code{alpha_PY}; for \code{prior="PY"} use both
#'   \code{alpha_PY} and \code{sigma_PY in (0,1)}.
#' @param beta_DM,H_DM Hyperparameters for the Dirichlet–Multinomial prior
#'   where \code{H_DM} gives the maximum number of allowed clusters.
#' @param gamma_GN Hyperparameter for the Gnedin prior.
#' @param n_iter,burnin Integers, total iterations and burn-in. Must satisfy \code{burnin < n_iter}.
#' @param init_x Optional integer vector of length \code{n} with initial labels (1-based).
#' @param store_z Logical, store latent \code{Z} draws (heavy; defaults \code{FALSE}).
#' @param verbose Logical, progress every 1000 iterations.
#'
#' @return A plain \code{list} with the following elements:
#' list with \code{x} (matrix \code{S x n} of raw labels),
#'         \code{K} (vector length \code{S})
#'
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 6L
#' w <- matrix(0L, n, n)
#' w[lower.tri(w)] <- rpois(sum(lower.tri(w)), 2)
#' rownames(w) <- colnames(w) <- paste0("P", seq_len(n))
#' fit <- gibbs_bt_sbm(w, prior="GN", gamma_GN=0.5, n_iter=200, burnin=100)
#' }
#' @importFrom mcclust comp.psm
#' @importFrom mcclust.ext minVI
#' @export
#'

# ----------------------- Main Gibbs sampler ----------------------------------
gibbs_bt_sbm <- function(
    w_ij,
    a = 4,
    prior = c("DP", "PY", "DM", "GN"),
    alpha_PY = NA_real_,
    sigma_PY = NA_real_,
    beta_DM  = NA_real_,
    H_DM     = NA_integer_,
    gamma_GN = NA_real_,
    n_iter = 2000, burnin = 1000,
    init_x = NULL,
    store_z = FALSE,
    verbose = TRUE
) {
  stopifnot(is.matrix(w_ij))
  N <- nrow(w_ij)
  stopifnot(N == ncol(w_ij))
  if (any(w_ij < 0)) stop("Negative entries not allowed.")
  n_ij <- w_ij + t(w_ij)
  if (any(w_ij > n_ij)) stop("w_ij must be <= n_ij elementwise.")
  if (burnin >= n_iter) stop("burnin < n_iter required.")
  if (a <= 0) stop("Gamma(a,b) needs a>0")

  prior <- match.arg(prior)
  if (prior == "DP" && !is.finite(alpha_PY)) stop("DP requires alpha_PY.")
  if (prior == "PY" && (!is.finite(alpha_PY) || !is.finite(sigma_PY) || sigma_PY <= 0 || sigma_PY >= 1))
    stop("PY requires alpha_PY and sigma_PY in (0,1).")
  if (prior == "DM" && (!is.finite(beta_DM) || !is.finite(H_DM) || H_DM < 1))
    stop("DM requires beta_DM>0 and H_DM>=1.")
  if (prior == "GN" && !is.finite(gamma_GN)) stop("GN requires gamma_GN.")

  urn_fun <- switch(
    prior,
    DP = function(v) urn_DP(v, alpha_PY),
    PY = function(v) urn_PY(v, alpha_PY, sigma_PY),
    DM = function(v) urn_DM(v, beta_DM, H_DM),
    GN = function(v) urn_GN(v, gamma_GN)
  )

  # precompute row wins
  w_i <- rowSums(w_ij)

  # --- init labels and rates ------------------------------------------------------
  L <- N                            # label capacity (can grow)
  if (is.null(init_x)) x_curr <- sample.int(L, N, replace = TRUE) else {
    stopifnot(length(init_x) == N)
    x_curr <- as.integer(init_x)
    L <- max(L, max(x_curr))
  }
  a_curr <- a;
  b_eff <- as.numeric(exp(digamma(a_curr)))   # b := exp(psi(a))

  lambda_curr <- rep(NA_real_, L)
  csize0 <- tabulate(x_curr, nbins = L)
  occ0 <- which(csize0 > 0L)
  lambda_curr[occ0] <- rgamma(length(occ0), shape = a_curr, rate = b_eff)

  # latent Z
  Z_curr <- matrix(0, N, N)
  if (store_z) z_store <- array(0, dim = c(n_iter - burnin, N, N)) else z_store <- NULL

  # storage
  S <- n_iter - burnin
  x_samples <- matrix(NA_integer_, S, N)
  lambda_list <- vector("list", S)         # ragged storage for lambdas
  L_hist <- integer(S)
  K_occ <- integer(S)


  save_i <- 0L
  for (iter in seq_len(n_iter)) {
    # ---- 1) Update Z and row sums in C++ (expects lambda length >= max(x))
    outZ <- updateZ_and_rowsums(n_ij = structure(n_ij, .Dim = dim(n_ij)),
                                x = as.integer(x_curr),
                                lambda = lambda_curr)
    Z_curr <- outZ$Z
    Z_row_sum <- outZ$rowSumsZ

    # ---- 2) Single-site updates for x_i
    for (i in seq_len(N)) {
      csize_minus <- tabulate(x_curr[-i], nbins = L)
      occupied <- which(csize_minus > 0L)
      H <- length(occupied)
      v_minus <- csize_minus[occupied]

      prior_vec <- if (H > 0L) urn_fun(v_minus) else c(0, urn_fun(integer(0))[1L])
      new_weight <- prior_vec[H + 1L]
      existing_weights <- if (H > 0L) prior_vec[seq_len(H)] else numeric(0)

      Zi <- Z_row_sum[i]
      wi <- w_i[i]

      logp <- rep(-Inf, H + 1L)

      if (H > 0L) {
        lam_h <- lambda_curr[occupied]
        llh <- wi * .safe_log(lam_h) - lam_h * Zi
        lprior <- .safe_log(existing_weights)
        logp[seq_len(H)] <- lprior + llh
      }
      # new (anchored):
      if (new_weight > 0) {
        logp[H + 1L] <- .safe_log(new_weight) +
          .new_cluster_integral_log_anchored(wi, Zi, a_curr, b_eff)
      }


      choice <- .sample_from_logweights(logp)

      if (choice <= H) {
        x_curr[i] <- occupied[choice]
      } else {
        # create / reuse empty label id
        csize_all <- tabulate(x_curr[-i], nbins = L)
        empties <- which(csize_all == 0L)
        if (length(empties) == 0L) {
          L <- L + 1L
          lambda_curr <- c(lambda_curr, NA_real_)
          new_label <- L
        } else new_label <- empties[1L]
        x_curr[i] <- new_label

        shape_k <- a_curr + wi
        rate_k  <- b_eff + Zi
        lambda_curr[new_label] <- rgamma(1L, shape = shape_k, rate = rate_k)
      }
    }

    # ---- 3) Update lambda_k for occupied clusters (conjugate)
    csize <- tabulate(x_curr, nbins = L)
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


    #numerically stable global rescaling so that E(log(lambda))=0
    occupied <- which(csize > 0L)
    if (length(occupied) > 0L) {
      g_mean <- exp(mean(log(lambda_curr[occupied])))
      lambda_curr[occupied] <- lambda_curr[occupied] / g_mean
    }

    # ---- 7) Store
    if (iter > burnin) {
      save_i <- save_i + 1L
      x_samples[save_i, ] <- x_curr
      lambda_list[[save_i]] <- lambda_curr
      L_hist[save_i] <- L
      K_occ[save_i] <- length(occ)
      if (store_z) z_store[save_i, , ] <- Z_curr
    }

    if (verbose && iter %% 1000L == 0L) {
      cat("iter", iter, "occupied =", length(unique(x_curr)), "\n")
    }
  }

  list(
    x_samples = x_samples,
    lambda_samples = lambda_list,  # ragged storage
    L_hist = L_hist,
    K_occ = K_occ,
    z_samples = if (store_z) z_store else NULL
    )
}




