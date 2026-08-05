test_that("PL latent-structure reference workloads stay bounded", {
  skip_on_cran()

  mixture_data <- sushi_toy(n_rankings = 30, rank_length = 5, seed = 21)
  mixture_time <- system.time({
    mixture_fit <- fit_btsbm(
      mixture_data,
      pl_mixture_model(),
      mcmc_control(iter = 24, warmup = 12, seed = 22)
    )
  })[["elapsed"]]

  lbm_data <- pl_lbm_toy(n_rankings = 24, rank_length = 4, seed = 23)
  lbm_time <- system.time({
    lbm_fit <- fit_btsbm(
      lbm_data,
      pl_lbm_model(),
      mcmc_control(iter = 20, warmup = 10, seed = 24)
    )
  })[["elapsed"]]

  # These deliberately broad limits catch accidental quadratic state growth or
  # an unbounded allocation loop without making normal CI timing brittle.
  expect_lt(mixture_time, 20)
  expect_lt(lbm_time, 20)
  expect_lt(as.numeric(object.size(mixture_fit)), 2e6)
  expect_lt(as.numeric(object.size(lbm_fit)), 2e6)
})

test_that("registered augmentation is no slower than the R reference workload", {
  skip_on_cran()
  set.seed(31)
  rankings <- t(replicate(400L, sample.int(12L, 6L)))
  ranking_cluster <- rep(1:4, length.out = nrow(rankings))
  latent_strength <- matrix(
    rep(seq(0.8, 2.4, length.out = 12L), 4L), nrow = 4L, byrow = TRUE
  )

  set.seed(32)
  r_elapsed <- system.time(
    BTSBM:::.pl_draw_grouped_augmentation_r(rankings, ranking_cluster, latent_strength)
  )[["elapsed"]]
  set.seed(32)
  cpp_elapsed <- system.time(
    BTSBM:::.pl_draw_grouped_augmentation(rankings, ranking_cluster, latent_strength)
  )[["elapsed"]]

  # A generous multiplier accommodates noisy shared CI while detecting a
  # regression that accidentally routes the production kernel through R.
  expect_lte(cpp_elapsed, 2 * r_elapsed + 0.05)
})
