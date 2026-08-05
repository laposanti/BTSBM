test_that("common fitting interface wraps the simple BT sampler", {
  wins <- matrix(c(
    0, 3, 2,
    1, 0, 4,
    0, 1, 0
  ), nrow = 3, byrow = TRUE)
  data <- as_bt_data(wins, item_labels = c("A", "B", "C"))
  fit <- fit_btsbm(
    data,
    bt_model(gamma_ability(shape = 1, rate = 1)),
    mcmc_control(iter = 30, warmup = 10, seed = 12)
  )

  expect_s3_class(fit, "btsbm_fit")
  expect_equal(dim(fit$draws$ability), c(20, 3))
  expect_equal(dim(log_lik(fit)), c(20, 3))
  expect_true(all(is.finite(log_lik(fit))))
})

test_that("common fitting interface respects the BT--SBM Gamma rate and capacity trace", {
  wins <- matrix(c(
    0, 3, 2,
    1, 0, 4,
    0, 1, 0
  ), nrow = 3, byrow = TRUE)
  data <- as_bt_data(wins)
  fit <- fit_btsbm(
    data,
    bt_sbm_model(
      ability_prior = gamma_ability(shape = 2, rate = 3),
      item_prior = partition_prior("gnedin", gamma = 0.5)
    ),
    mcmc_control(iter = 24, warmup = 12, seed = 13)
  )

  expect_s3_class(fit, "btsbm_fit")
  expect_equal(dim(fit$draws$item_cluster), c(12, 3))
  expect_equal(length(fit$draws$n_item_clusters), 12)
  expect_true(all(is.finite(log_lik(fit))))

  raw <- gibbs_bt_sbm(
    wins, a = 2, b = 3, prior = "GN", gamma_GN = 0.5,
    T_iter = 12, T_burn = 6, verbose = FALSE
  )
  expect_equal(length(raw$L_cap_per_iter), 6)
  expect_true(all(raw$L_cap_per_iter >= raw$K_per_iter))
})

test_that("BT--SBM simulation honours its supplied seed", {
  first <- sample_from_BTSBM(n_players = 6, K = 2, seed = 100)
  second <- sample_from_BTSBM(n_players = 6, K = 2, seed = 100)
  third <- sample_from_BTSBM(n_players = 6, K = 2, seed = 101)

  expect_equal(first$w, second$w)
  expect_false(identical(first$w, third$w))
})
