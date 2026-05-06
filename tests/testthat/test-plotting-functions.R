test_that("block_winprob_table_from_lambdas returns expected block summaries", {
  lambda_draws <- rbind(
    c(2, 1, 4, 3),
    c(2, 2, 5, 3)
  )
  z <- c(1, 1, 2, 2)
  N_ij <- matrix(c(
    0, 1, 2, 1,
    1, 0, 1, 2,
    2, 1, 0, 1,
    1, 2, 1, 0
  ), nrow = 4, byrow = TRUE)

  out <- block_winprob_table_from_lambdas(
    lambda_item = lambda_draws,
    z = z,
    N_ij = N_ij,
    include_diag = TRUE
  )

  expect_equal(dim(out$P_block), c(2, 2))
  expect_equal(out$block_sizes, c(Block_1 = 2L, Block_2 = 2L))
  expect_true(all(out$lambda_point > 0))
  expect_true(is.finite(out$P_block[1, 2]))
  expect_true(is.finite(out$P_block[2, 1]))
})

test_that("compute_expected_wins_rank_posterior summarises posterior ranks", {
  lambda_draws <- rbind(
    c(5, 3, 1),
    c(4, 2, 1),
    c(6, 3, 2)
  )

  out <- compute_expected_wins_rank_posterior(
    lambda_item = lambda_draws,
    player_names = c("A", "B", "C")
  )

  expect_equal(colnames(out$rank_samples), c("A", "B", "C"))
  expect_equal(out$summary$player, c("A", "B", "C"))
  expect_true(out$summary$rank_mean[1] < out$summary$rank_mean[2])
  expect_equal(nrow(out$rank_prob), 3)
})

test_that("exploratory_adjacency returns ggplot objects without preview", {
  skip_if_not_installed("ggplot2")
  Y_ij <- matrix(c(
    0, 1, 2,
    0, 0, 1,
    1, 1, 0
  ), nrow = 3, byrow = TRUE)
  rownames(Y_ij) <- colnames(Y_ij) <- c("alpha_player", "beta_player", "gamma_player")

  p <- exploratory_adjacency(
    Y_ij = Y_ij,
    order_by = "none",
    show_margins = FALSE,
    bw_preview = FALSE
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot_assignment_probabilities returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  fit <- list(
    x_samples = rbind(c(1, 1, 2), c(1, 2, 2), c(1, 1, 2)),
    lambda_samples = list(c(4, 2), c(4, 2), c(5, 2)),
    item_cluster_assignment_probs = data.frame(
      Cluster_1 = c(1, 2 / 3, 0),
      Cluster_2 = c(0, 1 / 3, 1),
      Cluster_3 = c(0, 0, 0)
    ),
    minVI_partition = c(1, 1, 2)
  )
  w_ij <- matrix(c(
    0, 2, 1,
    1, 0, 0,
    1, 2, 0
  ), nrow = 3, byrow = TRUE)
  rownames(w_ij) <- colnames(w_ij) <- c("alpha", "beta", "gamma")

  p <- plot_assignment_probabilities(
    fit = fit,
    w_ij = w_ij,
    x_hat = fit$minVI_partition
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot_lambda_uncertainty and plot_rank_intervals return ggplot objects", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("coda")

  fit <- list(
    lambda_samples_relabel = rbind(
      c(4, 3, 1),
      c(5, 2.5, 1.2),
      c(4.5, 2.7, 1.1)
    ),
    x_samples_relabel = rbind(
      c(1, 1, 2),
      c(1, 1, 2),
      c(1, 1, 2)
    ),
    minVI_partition = c(1, 1, 2)
  )
  w_ij <- matrix(c(
    0, 2, 3,
    1, 0, 2,
    0, 1, 0
  ), nrow = 3, byrow = TRUE)
  rownames(w_ij) <- colnames(w_ij) <- c("alpha", "beta", "gamma")

  p_lambda <- plot_lambda_uncertainty(
    fit = fit,
    w_ij = w_ij,
    x_hat = fit$minVI_partition
  )
  p_rank <- plot_rank_intervals(fit)

  expect_s3_class(p_lambda, "ggplot")
  expect_s3_class(p_rank, "ggplot")
  expect_true(is.data.frame(attr(p_lambda, "player_summary")))
})

test_that("plot_block_adjacency returns ggplot when optional packages are available", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("reshape2")
  skip_if_not_installed("ggside")

  fit <- list(minVI_partition = c(1, 1, 2))
  w_ij <- matrix(c(
    0, 2, 3,
    1, 0, 1,
    0, 2, 0
  ), nrow = 3, byrow = TRUE)
  rownames(w_ij) <- colnames(w_ij) <- c("alpha", "beta", "gamma")

  p <- plot_block_adjacency(
    fit = fit,
    w_ij = w_ij,
    bw_preview = FALSE
  )

  expect_s3_class(p, "ggplot")
})
