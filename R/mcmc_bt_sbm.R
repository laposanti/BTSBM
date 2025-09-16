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
    w_ij, n_ij,              # K x K pairwise data
    a, b,                    # Gamma(a,b) prior for each block rate
    prior = "DP",            # which prior among c("DP","PY","DM","GN")
    alpha_PY = NA,           # DP/PY parameter
    sigma_PY = NA,           # discount for PY
    beta_DM = NA,            # param for DM
    H_DM = NA,               # max # clusters for DM
    gamma_GN = NA,           # param for Gnedin
    n_iter = 2000,
    burnin = 1000,
    init_x = NULL,           # optional initialization
    store_z = FALSE,         # store the latent Z?
    verbose = TRUE
) {
  K <- nrow(w_ij)
  stopifnot(K == ncol(w_ij), K == nrow(n_ij), K == ncol(n_ij))

  # Precompute w_i
  w_i <- rowSums(w_ij)

  # Init cluster assignment
  if (is.null(init_x)) {
    x_curr <- sample.int(K, K, replace = TRUE)
  } else {
    x_curr <- init_x
  }
  a_current = a
  b=b
  # We keep a block-rate vector of length K
  lambda_curr <- rgamma(K, shape = a_current, rate = b)

  # Latent Z
  Z_curr <- matrix(0, K, K)
  z_store <- if (store_z) array(0, dim = c(n_iter-burnin, K, K)) else NULL

  x_samples <- matrix(NA, n_iter-burnin, K)
  lambda_samples <- matrix(NA, n_iter-burnin, K)

  # ------------------------------------------------------------------
  # "urn_fun" picks the correct prior function to produce c(...)
  # existing cluster weights, new cluster weight
  # ------------------------------------------------------------------
  urn_fun <- function(v_minus) {
    if (prior == "DP") {
      urn_DP(v_minus, alpha_PY)
    } else if (prior == "PY") {
      urn_PY(v_minus, alpha_PY, sigma_PY)
    } else if (prior == "DM") {
      urn_DM(v_minus, beta_DM, H_DM)
    } else if (prior == "GN") {
      urn_GN(v_minus, gamma_GN)
    } else {
      stop("Invalid prior type chosen.")
    }
  }

  #damien sampler for the hyperparameters

  damien_sampler_a <- function(a,                # current a
                               lambda_vec,       # vector of lambdas (length K)
                               b,                # known scale hyperparam
                               prior_draw,       # function() that draws from p(a)
                               prior_support,    # function(a) that returns TRUE if in supp
                               max_retries = 1000) {
    # 1) Compute l1(a_current) and l2(a_current) for the U_1, U_2 draws
    K <- length(lambda_vec)
    gamma_val <- b^K * prod(lambda_vec)         # "b^K * prod(lambda_i)"

    l1 <- function(a) { gamma_val^a }
    l2 <- function(a) { (gamma( a ))^(-K) }

    # 2) Sample U1, U2 from Uniform(0, l1(a)), Uniform(0, l2(a)) respectively
    U1 <- runif(1, min=0, max=l1(a_current))
    U2 <- runif(1, min=0, max=l2(a_current))

    # 3) We now want to sample 'a_new' from p(a) restricted to
    #       { l1(a) > U1 } intersect { l2(a) > U2 } = { gamma_val^a > U1 } intersect { gamma(a)^(-K) > U2 }.
    #    We'll do a simple rejection loop.  If your prior is well-behaved, this should be efficient enough.

    for (trial in seq_len(max_retries)) {
      a_prop <- prior_draw()             # draw from p(a)

      # Check if in prior support
      if (!prior_support(a_prop)) next

      # Check the slice constraints:
      cond1 <- ( gamma_val^a_prop > U1 )
      cond2 <- ( (gamma(a_prop))^(-K) > U2 )  # i.e. gamma(a_prop) < U2^(-1/K)

      if (cond1 && cond2) {
        # success!
        return(list(a_new = a_prop, u1 = U1, u2 = U2, accepted=TRUE))
      }
    }

    # If we got here, we never found a suitable a in 'max_retries' tries
    return(list(a_new = a, u1 = U1, u2 = U2, accepted=FALSE))
  }
  # "prior_draw" draws from that Gamma
  prior_draw <- function() runif(1,0,1)

  # "prior_support" just requires a>0
  prior_support <- function(a) (a>0)

  # sumZi
  sumZi <- function(i) sum(Z_curr[i, ])

  # integral over a new lambda from Gamma(a,b)
  # for node i alone in a new block:

  # Log version
  new_cluster_integral_log <- function(w_i_val, Z_i_val,a_current) {
    # log( top / bot ) = log(top) - log(bot)
    log_top <- (a_current*log(b)) + lgamma(a_current + w_i_val)
    log_bot <- lgamma(a_current) + (a_current + w_i_val)*log(b + Z_i_val)
    log_top - log_bot
  }
  save_i = 0
  for (iter in seq_len(n_iter)) {

    # 1) Update Z_{ij}
    for (i in 1:(K-1)) {
      lam_i <- lambda_curr[x_curr[i]]
      for (j in (i+1):K) {
        nij <- n_ij[i,j]
        if (nij > 0) {
          lam_j <- lambda_curr[x_curr[j]]
          Z_val <- rgamma(1, shape = nij, rate = lam_i + lam_j)
          Z_curr[i,j] <- Z_val
          Z_curr[j,i] <- Z_val
        } else {
          Z_curr[i,j] <- 0
          Z_curr[j,i] <- 0
        }
      }
    }


    # 3) Single-site update for x_i
    for (i in seq_len(K)) {
      old_k <- x_curr[i]

      # cluster sizes ignoring i
      csize <- integer(K)
      for (k in seq_len(K)) {
        csize[k] <- sum(x_curr[-i] == k)
      }

      occupied <- which(csize > 0)
      H <- length(occupied)
      v_minus <- csize[occupied]

      # get prior weights from urn function
      # => vector of length H+1
      prior_vec <- urn_fun(v_minus)
      new_weight <- prior_vec[H+1]
      existing_weights <- prior_vec[1:H]

      existing_ids <- occupied
      Zi_val <- sumZi(i)
      wi_val <- w_i[i]

      # We'll build a vector log_probs of length H+1
      log_probs <- numeric(H+1)

      # existing
      for (h in seq_len(H)) {
        k_id <- existing_ids[h]
        lam_k <- lambda_curr[k_id]

        # log-likelihood factor = w_i[i]*log(lam_k) - lam_k * Zi_val
        # log-prior factor = log( existing_weights[h] )
        # so total is:
        # log_probs[h] = log( existing_weights[h] ) + [ w_i[i]*log(lam_k) - lam_k*Zi_val ]


        llh_ex <- wi_val*log(lam_k) - lam_k*Zi_val
        lprior_ex <- log(existing_weights[h] + 1e-300)  # +1e-300 in case it's 0
        log_probs[h] <- lprior_ex + llh_ex
      }


      # new cluster
      if (new_weight > 0) {
        # log of new_weight + log of new_cluster_integral
        lprior_new <- log(new_weight)
        lint_new <- new_cluster_integral_log(wi_val, Zi_val,a_current)
        log_probs[H+1] <- lprior_new + lint_new
      } else {
        log_probs[H+1] <- -Inf
      }

      # Numeric stability:
      offset <- max(log_probs)
      shifted <- log_probs - offset
      unnorm <- exp(shifted)
      denom <- sum(unnorm)
      probs <- unnorm / denom

      choice <- sample.int(H+1, size=1, prob=probs)

      if (choice <= H) {
        x_curr[i] <- existing_ids[choice]
      } else {
        # new cluster => pick an empty label
        csize <- integer(K)
        for (k in seq_len(K)) {
          csize[k] <- sum(x_curr[-i] == k)
        }

        empties   <- which(csize == 0)
        new_label <- empties[1]
        x_curr[i] <- new_label
        shape_k   <- a_current + w_i[i]
        rate_k    <- b
        rate_k    <- rate_k + sumZi(i)

        lambda_curr[new_label] <- rgamma(1, shape = shape_k, rate = rate_k)
        #2.1) Update a for Gamma
        csize <- integer(K)
        for (k in seq_len(K)) {
          csize[k] <- sum(x_curr == k)
        }


      }
    }



    # 2) Update lambda_k for each occupied cluster
    for (k in seq_len(K)) {
      members_k <- which(x_curr == k)
      if (length(members_k) > 0) {
        shape_k <- a_current + sum(w_i[members_k])
        rate_k  <- b
        for (i in members_k) {
          rate_k <- rate_k + sumZi(i)
        }
        lambda_curr[k] <- rgamma(1, shape = shape_k, rate = rate_k)
      }
    }

    # cluster sizes ignoring i
    csize <- integer(K)
    for (k in seq_len(K)) {
      csize[k] <- sum(x_curr == k)
    }

    pi_curr = lambda_curr/sum(lambda_curr)

    lambda_curr = pi_curr*rgamma(1,K*a_current,b)

    #2.1) Update a for Gamma prior
    # out <- damien_sampler_a(a=a_current,
    #                         lambda_vec=lambda_curr,
    #                         b=b,
    #                         prior_draw=prior_draw,
    #                         prior_support=prior_support)
    #
    # #
    # if (out$accepted) {
    #   a_current <- out$a_new
    # } else {
    #   # fallback: keep old value if we run out of tries
    # }

    #
    # store
    if(iter>burnin)
      save_i = save_i+1
    x_samples[save_i, ] <- x_curr
    lambda_samples[save_i, ] <- lambda_curr
    if (store_z) {
      z_store[save_i,,] <- Z_curr
    }

    if (verbose && iter %% 200 == 0) {
      cat("Iteration:", iter,
          "- #occupied =", length(unique(x_curr)), "\n")
    }
  }

  list(
    x_samples      = x_samples,
    lambda_samples = lambda_samples,
    z_samples      = z_store
  )
}
