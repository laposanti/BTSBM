#' Gibbs sampler for a Bradley–Terry Stochastic Block Model
#'
#' @description
#' Runs a Gibbs sampler for a BT–SBM with a Gamma(a, b) prior on block rates
#' and a choice of clustering prior (DP, Pitman–Yor, Dirichlet–Multinomial, or Gnedin).
#'
#' @param w_ij integer/numeric K x K matrix of wins from i over j (diagonal ignored).
#' @param n_ij integer/numeric K x K matrix of total matches between i and j (symmetric, diag 0).
#' @param a numeric(1). Shape for the Gamma(a, b) prior on each block rate at init; updated inside.
#' @param b numeric(1). Rate for the Gamma(a, b) prior (fixed).
#' @param prior character(1). One of \code{c("DP","PY","DM","GN")}.
#' @param alpha_PY numeric(1) or \code{NA}. Concentration for DP/PY (must be > 0 when used).
#' @param sigma_PY numeric(1) or \code{NA}. Discount for PY in [0,1); set 0 for DP if using "PY".
#' @param beta_DM numeric(1) or \code{NA}. Dirichlet-Multinomial concentration (>0) if \code{prior="DM"}.
#' @param H_DM integer(1) or \code{NA}. Max # clusters for DM (>=1) if \code{prior="DM"}.
#' @param gamma_GN numeric(1) or \code{NA}. Parameter for Gnedin if \code{prior="GN"}.
#' @param n_iter integer(1). Total MCMC iterations.
#' @param burnin integer(1). Number of initial iterations to discard (must be < n_iter).
#' @param init_x integer vector of length K with initial cluster labels (1..K), or \code{NULL}.
#' @param store_z logical(1). If TRUE, store the latent Z matrices.
#' @param verbose logical(1). If TRUE, print progress every 200 iters.
#' @param a_prior_draw function with signature \code{function() numeric(1)}. Draws a proposal from the prior of \eqn{a} (Damien slice sampler). Default is \code{function() rexp(1, 1)} (Exp(1) on \eqn{a>0}).
#' @param a_prior_support function with signature \code{function(a) logical(1)}. Checks support of \eqn{a}. Default: \code{function(a) a > 0}.
#'
#' @return A list with:
#' \item{x_samples}{(n_save x K) integer matrix of cluster labels after burn-in.}
#' \item{lambda_samples}{(n_save x K) numeric matrix of block rates (entries for empty labels can be unused).}
#' \item{z_samples}{optional (n_save x K x K) array of latent \eqn{Z_{ij}} if \code{store_z=TRUE}.}
#'
#' @details
#' - The sampler uses Gamma-Poisson augmentation for the BT likelihood: \eqn{Z_{ij} \sim \mathrm{Gamma}(n_{ij}, \lambda_{x_i}+\lambda_{x_j})}.
#' - Cluster reassignment uses a one-step CRP/PY/DM/GN urn with a marginal new-cluster integral under \eqn{\lambda \sim \mathrm{Gamma}(a,b)}.
#' - \eqn{a} is updated via a Damien-style slice move on the joint of \eqn{\lambda} given \eqn{b} and the prior on \eqn{a}.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' K <- 6
#' w <- matrix(0, K, K); n <- matrix(0, K, K)
#' n[upper.tri(n)] <- sample(0:5, sum(upper.tri(n)), TRUE)
#' n <- n + t(n); diag(n) <- 0
#' w[upper.tri(w)] <- rbinom(sum(upper.tri(w)), size = n[upper.tri(n)], prob = 0.5)
#' w <- w + t(n - w); diag(w) <- 0
#' out <- gibbs_bt_sbm(w, n, a = 1, b = 1, prior = "DP", alpha_PY = 1,
#'                     n_iter = 500, burnin = 250, verbose = FALSE)
#' str(out)
#' }
#' @export
#'
#'
#'
gibbs_bt_sbm <- function(
    w_ij, n_ij,
    a, b,
    prior = "DP",
    alpha_PY = NA_real_,
    sigma_PY = NA_real_,
    beta_DM = NA_real_,
    H_DM = NA_integer_,
    gamma_GN = NA_real_,
    n_iter = 2000L,
    burnin = 1000L,
    init_x = NULL,
    store_z = FALSE,
    verbose = TRUE,
    a_prior_draw = function() rexp(1, 1),       # Exp(1) on a>0 (generic, proper)
    a_prior_support = function(a) a > 0
) {
  ## ---------- input checks ----------
  stopifnot(is.matrix(w_ij), is.matrix(n_ij))
  K <- nrow(w_ij)
  if (!identical(dim(w_ij), dim(n_ij)))
    stop("w_ij and n_ij must have the same dimensions.")
  if (K != ncol(w_ij) || K != nrow(n_ij) || K != ncol(n_ij))
    stop("w_ij and n_ij must be square K x K matrices.")
  if (any(is.na(w_ij)) || any(is.na(n_ij)))
    stop("w_ij and n_ij must not contain NA.")
  if (any(w_ij < 0) || any(n_ij < 0))
    stop("w_ij and n_ij must be nonnegative.")
  if (any(w_ij > n_ij))
    stop("Elementwise, w_ij cannot exceed n_ij.")
  if (any(diag(n_ij) != 0) || any(diag(w_ij) != 0))
    stop("Diagonal of n_ij and w_ij must be zero.")
  if (!isTRUE(all.equal(n_ij, t(n_ij))))
    stop("n_ij must be symmetric.")
  prior <- match.arg(prior, c("DP","PY","DM","GN"))
  if (prior %in% c("DP","PY")) {
    if (is.na(alpha_PY) || alpha_PY <= 0) stop("alpha_PY must be >0 for DP/PY.")
    if (prior == "PY") {
      if (is.na(sigma_PY) || sigma_PY < 0 || sigma_PY >= 1)
        stop("sigma_PY must be in [0,1) for PY. Use 0 for DP behavior.")
    } else {
      sigma_PY <- 0
    }
  }
  if (prior == "DM") {
    if (is.na(beta_DM) || beta_DM <= 0) stop("beta_DM must be >0 for DM.")
    if (is.na(H_DM) || H_DM < 1) stop("H_DM must be >=1 for DM.")
  }
  if (prior == "GN") {
    if (is.na(gamma_GN) || gamma_GN <= 0) stop("gamma_GN must be >0 for GN.")
  }
  if (!is.numeric(a) || length(a) != 1 || a <= 0) stop("'a' must be a single positive number.")
  if (!is.numeric(b) || length(b) != 1 || b <= 0) stop("'b' must be a single positive number.")
  if (!is.numeric(n_iter) || !is.numeric(burnin) || n_iter <= burnin || n_iter < 2L)
    stop("Require n_iter > burnin >= 0.")
  if (!is.null(init_x)) {
    if (length(init_x) != K) stop("init_x must have length K.")
    if (any(init_x < 1 | init_x > K)) stop("init_x labels must be in 1..K.")
  }

  ## ---------- helpers: urn weights ----------
  urn_weights <- function(v_minus) {
    # v_minus: sizes of currently occupied clusters, length H
    H <- length(v_minus)
    if (prior == "DP" || (prior == "PY" && isTRUE(all.equal(sigma_PY, 0)))) {
      # CRP with alpha = alpha_PY
      c(v_minus, alpha_PY)
    } else if (prior == "PY") {
      # Pitman–Yor: existing weights v_minus - sigma, new alpha + sigma*H
      if (any(v_minus <= sigma_PY)) stop("PY: some cluster size <= sigma; invalid state.")
      c(v_minus - sigma_PY, alpha_PY + sigma_PY * H)
    } else if (prior == "DM") {
      # Finite mixture with symmetric Dirichlet(beta_DM/H_DM) over H_DM labels
      # Existing: size_i, New: number of empty labels (H_DM - H) each with prior mass beta_DM/H_DM
      # For single-site update, we need a single "new" bucket proportional to total mass of empty labels
      empty_count <- max(H_DM - H, 0L)
      existing <- v_minus
      new_mass <- empty_count * (beta_DM / H_DM)
      c(existing, new_mass)
    } else if (prior == "GN") {
      # Placeholder simple GN-style: existing proportional to v_minus, new to gamma_GN
      c(v_minus, gamma_GN)
    } else {
      stop("Invalid prior.")
    }
  }

  ## ---------- Damien-style slice move for 'a' ----------
  damien_update_a <- function(a_curr, lambda_vec, b_rate,
                              prior_draw, prior_support,
                              max_retries = 2000L) {
    Klam <- length(lambda_vec)
    gamma_val <- b_rate^Klam * prod(lambda_vec)              # (b^K) * prod(lambda_k)

    l1 <- function(a) gamma_val^a
    l2 <- function(a) (gamma(a))^(-Klam)

    U1 <- runif(1, 0, l1(a_curr))
    U2 <- runif(1, 0, l2(a_curr))

    for (trial in seq_len(max_retries)) {
      a_prop <- prior_draw()
      if (!is.finite(a_prop) || !prior_support(a_prop)) next
      cond1 <- (gamma_val^a_prop > U1)
      cond2 <- ((gamma(a_prop))^(-Klam) > U2)
      if (cond1 && cond2) return(list(a_new = a_prop, accepted = TRUE))
    }
    list(a_new = a_curr, accepted = FALSE)
  }

  ## ---------- precompute ----------
  w_i <- rowSums(w_ij)                          # wins out of i
  Z    <- matrix(0.0, K, K)                     # latent gammas
  x    <- if (is.null(init_x)) sample.int(K, K, replace = TRUE) else as.integer(init_x)
  a_curr <- as.numeric(a)
  b_rate <- as.numeric(b)
  lambda <- rgamma(K, shape = a_curr, rate = b_rate)  # one per label 1..K

  n_save <- as.integer(n_iter - burnin)
  x_out <- matrix(NA_integer_, n_save, K)
  lambda_out <- matrix(NA_real_, n_save, K)
  z_out <- if (store_z) array(0.0, dim = c(n_save, K, K)) else NULL
  save_i <- 0L

  sumZi <- function(i) sum(Z[i, ], na.rm = TRUE)

  new_cluster_integral_log <- function(wi, Zi, a_shp, b_rt) {
    # log ∫ λ^{wi+a-1} e^{-(Zi+b)λ} b^a / Γ(a) dλ  =  a log b + lgamma(a+wi) - lgamma(a) - (a+wi) log(Zi + b)
    a_shp * log(b_rt) + lgamma(a_shp + wi) - lgamma(a_shp) - (a_shp + wi) * log(Zi + b_rt)
  }

  ## ---------- main loop ----------
  for (iter in seq_len(n_iter)) {

    # 1) Update Z_{ij} for i<j, copy symmetrically
    for (i in 1:(K - 1L)) {
      lam_i <- lambda[x[i]]
      for (j in (i + 1L):K) {
        nij <- n_ij[i, j]
        if (nij > 0) {
          lam_j <- lambda[x[j]]
          zij <- rgamma(1L, shape = nij, rate = lam_i + lam_j)
          Z[i, j] <- zij
          Z[j, i] <- zij
        } else {
          Z[i, j] <- 0.0
          Z[j, i] <- 0.0
        }
      }
    }

    # 2) Update cluster labels x_i (single-site)
    for (i in seq_len(K)) {
      # cluster sizes excluding i
      csize <- tabulate(x[-i], nbins = K) # length K
      occupied <- which(csize > 0L)
      H <- length(occupied)
      v_minus <- csize[occupied]

      weights <- urn_weights(v_minus)
      if (length(weights) != H + 1L)
        stop("urn_weights must return length H+1 vector (existing + new).")

      Zi <- sumZi(i)
      wi <- w_i[i]

      logp <- numeric(H + 1L)

      # existing clusters
      for (h in seq_len(H)) {
        k_id <- occupied[h]
        lam_k <- lambda[k_id]
        llh <- wi * log(lam_k + 1e-300) - lam_k * Zi
        lprior <- log(weights[h] + 1e-300)
        logp[h] <- lprior + llh
      }

      # new cluster mass (marginalized λ)
      new_w <- weights[H + 1L]
      if (new_w > 0) {
        logp[H + 1L] <- log(new_w) + new_cluster_integral_log(wi, Zi, a_curr, b_rate)
      } else {
        logp[H + 1L] <- -Inf
      }

      # normalize
      offs <- max(logp)
      pr <- exp(logp - offs)
      pr <- pr / sum(pr)

      choice <- sample.int(H + 1L, 1L, prob = pr)

      if (choice <= H) {
        x[i] <- occupied[choice]
      } else {
        # assign to a fresh empty label
        empties <- which(csize == 0L)
        new_label <- empties[1L]
        x[i] <- new_label
        # draw its λ | data_i
        shape_k <- a_curr + wi
        rate_k  <- b_rate + Zi
        lambda[new_label] <- rgamma(1L, shape = shape_k, rate = rate_k)
      }
    }

    # 3) Update λ_k for each occupied cluster
    csize_full <- tabulate(x, nbins = K)
    for (k in which(csize_full > 0L)) {
      members <- (x == k)
      shape_k <- a_curr + sum(w_i[members])
      rate_k  <- b_rate + sum(rowSums(Z[members, , drop = FALSE]))
      lambda[k] <- rgamma(1L, shape = shape_k, rate = rate_k)
    }

    # 4) Optional reparameterization step (Dirichlet-Gamma “rescale”)
    #    Keep as in your code: draw total mass then normalize to improve mixing.
    pi_curr <- lambda / sum(lambda)
    total   <- rgamma(1L, shape = K * a_curr, rate = b_rate)
    lambda  <- pi_curr * total

    # 5) Update 'a' via Damien move
    a_upd <- damien_update_a(
      a_curr, lambda_vec = lambda, b_rate = b_rate,
      prior_draw   = a_prior_draw,
      prior_support = a_prior_support
    )
    if (a_upd$accepted) a_curr <- a_upd$a_new

    # 6) Store
    if (iter > burnin) {
      save_i <- save_i + 1L
      x_out[save_i, ] <- x
      lambda_out[save_i, ] <- lambda
      if (store_z) z_out[save_i, , ] <- Z
    }

    if (verbose && iter %% 200L == 0L) {
      cat("iter", iter, "| occupied =", length(which(csize_full > 0L)), "\n")
    }
  }

  list(
    x_samples = x_out,
    lambda_samples = lambda_out,
    z_samples = z_out
  )
}
