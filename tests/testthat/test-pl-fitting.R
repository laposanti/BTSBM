test_that("simple PL fitting returns a common fit object and ranking-level likelihood", {
  data <- sushi_toy(n_rankings = 16, rank_length = 4, seed = 1)
  fit <- fit_btsbm(
    data,
    pl_model(),
    mcmc_control(iter = 40, warmup = 20, seed = 2)
  )

  expect_s3_class(fit, "btsbm_fit")
  expect_equal(dim(fit$draws$latent_strength), c(20, 10))
  expect_equal(dim(fit$draws$item_strength), c(20, 10))
  expect_equal(dim(implied_ability(fit)), c(20, 10))
  expect_equal(length(posterior_strength(fit, summary = "mean")), 10)

  ll <- log_lik(fit)
  expect_equal(dim(ll), c(20, 16))
  expect_true(all(is.finite(ll)))
  expect_equal(attr(ll, "unit_index"), data$ranking_id)
  expect_equal(summary(fit)$n_draws, 20L)
})

test_that("PL--SBM fitting records item partition and compatible likelihood", {
  data <- sushi_toy(n_rankings = 12, rank_length = 4, seed = 3)
  fit <- fit_btsbm(
    data,
    pl_sbm_model(item_clustering = gnedin_prior(gnedin_hyperparameter = 0.6)),
    mcmc_control(iter = 35, warmup = 15, seed = 4)
  )

  expect_s3_class(fit, "btsbm_fit")
  expect_equal(dim(fit$draws$item_cluster), c(20, 10))
  expect_equal(length(fit$draws$n_item_clusters), 20)
  expect_true(all(fit$draws$n_item_clusters >= 1L))
  expect_equal(dim(implied_ability(fit)), c(20, 10))
  expect_true(all(vapply(fit$draws$latent_strength, length, integer(1)) == fit$draws$n_item_clusters))
  expect_true(all(is.finite(log_lik(fit))))

  similarity <- posterior_similarity(fit)
  expect_equal(dim(similarity), c(10, 10))
  expect_equal(unname(diag(similarity)), rep(1, 10))
  expect_equal(similarity, t(similarity))
})

test_that("PL--SBM accepts every exposed clustering-prior family", {
  data <- sushi_toy(n_rankings = 10, rank_length = 3, seed = 12)
  priors <- list(
    gnedin_prior(gnedin_hyperparameter = 0.6),
    dirichlet_process_prior(concentration = 1),
    pitman_yor_prior(concentration = 1, discount = 0.2),
    finite_partition(max_clusters = 3, concentration = 1)
  )

  fits <- lapply(priors, function(item_clustering) {
    fit_btsbm(
      data,
      pl_sbm_model(item_clustering = item_clustering),
      mcmc_control(iter = 14, warmup = 7, seed = 13)
    )
  })

  expect_true(all(vapply(fits, function(fit) {
    all(is.finite(posterior_strength(fit, summary = "mean")))
  }, logical(1))))
})

test_that("PL mixture stores ranking groups, likelihoods, and diagnostics", {
  data <- sushi_toy(n_rankings = 10, rank_length = 3, seed = 9)
  control <- mcmc_control(iter = 16, warmup = 8, seed = 10)
  model <- pl_mixture_model(
    ranking_clustering = dirichlet_process_prior(concentration = 1)
  )
  expect_no_warning(fit <- fit_btsbm(data, model, control))

  expect_equal(dim(fit$draws$ranking_cluster), c(8, 10))
  expect_equal(length(fit$draws$n_ranking_clusters), 8)
  expect_true(all(fit$draws$n_ranking_clusters >= 1L))
  expect_true(all(vapply(fit$draws$latent_strength, is.matrix, logical(1))))
  expect_true(all(is.finite(log_lik(fit))))
  expect_equal(dim(posterior_similarity(fit, target = "ranking")), c(10, 10))

  diagnostic <- mcmc_diagnostics(fit)
  expect_setequal(diagnostic$quantity, c("n_ranking_clusters", "log_likelihood"))
  expect_true(all(is.finite(diagnostic$ess)))
})

test_that("PL--LBM stores both partitions and has deterministic seeded draws", {
  data <- pl_lbm_toy(n_rankings = 10, rank_length = 3, seed = 11)
  model <- pl_lbm_model(
    ranking_clustering = finite_partition(max_clusters = 3),
    item_clustering = finite_partition(max_clusters = 3)
  )
  control <- mcmc_control(iter = 16, warmup = 8, seed = 12)
  expect_no_warning(first <- fit_btsbm(data, model, control))
  expect_no_warning(second <- fit_btsbm(data, model, control))

  expect_equal(dim(first$draws$ranking_cluster), c(8, 10))
  expect_equal(dim(first$draws$item_cluster), c(8, 6))
  expect_equal(first$draws$ranking_cluster, second$draws$ranking_cluster)
  expect_equal(first$draws$item_cluster, second$draws$item_cluster)
  expect_true(all(is.finite(log_lik(first))))
  expect_equal(dim(posterior_similarity(first, target = "item")), c(6, 6))
  expect_equal(dim(posterior_similarity(first, target = "ranking")), c(10, 10))

  diagnostic <- mcmc_diagnostics(first)
  expect_setequal(
    diagnostic$quantity,
    c("n_ranking_clusters", "n_item_clusters", "log_likelihood")
  )
})
