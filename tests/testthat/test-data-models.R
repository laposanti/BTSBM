test_that("as_bt_data validates a paired-comparison matrix", {
  wins <- matrix(c(
    0, 3, 1,
    2, 0, 4,
    1, 0, 0
  ), nrow = 3, byrow = TRUE)
  rownames(wins) <- colnames(wins) <- c("A", "B", "C")

  data <- as_bt_data(wins)
  expect_s3_class(data, "btsbm_pairwise_data")
  expect_equal(data$item_labels, c("A", "B", "C"))
  expect_error(as_bt_data(matrix(0, 1, 1)), "at least two")
  wins[1, 1] <- 1
  expect_error(as_bt_data(wins), "zero diagonal")
})

test_that("as_rankings validates complete and top-m ranking matrices", {
  rankings <- rbind(
    c(1, 2, 3),
    c(3, 1, 2),
    c(2, 3, 1)
  )
  data <- as_rankings(rankings, item_count = 4, item_labels = paste0("S", 1:4))

  expect_s3_class(data, "btsbm_ranking_data")
  expect_equal(data$n_items, 4L)
  expect_equal(data$rank_length, 3L)
  expect_equal(unname(data$rankings), rankings)
  expect_equal(data$item_labels, paste0("S", 1:4))
  expect_error(as_rankings(rbind(c(1, 1), c(2, 3))), "distinct")
  expect_error(as_rankings(rbind(c(1, 4), c(2, 3)), item_count = 3), "1:item_count")
  expect_error(as_rankings(matrix(c(1, 2.5), nrow = 1)), "integer")

  legacy <- as_pl_data(rankings, n_items = 4, item_labels = paste0("S", 1:4))
  expect_equal(legacy$rankings, data$rankings)
})

test_that("PL constructors use descriptive strength and clustering names", {
  prior <- gnedin_prior(gnedin_hyperparameter = 0.6)
  strength <- latent_strength(shape = 2, rate = 3)
  model <- pl_sbm_model(latent_strength = strength, item_clustering = prior)
  control <- mcmc_control(
    iter = 50, warmup = 20, seed = 8,
    partition_moves = partition_moves(split_merge = 1, pair_proposal = "informed")
  )

  expect_s3_class(model, "btsbm_model")
  expect_equal(model$likelihood, "pl")
  expect_equal(model$structure, "item_sbm")
  expect_identical(model$latent_strength, strength)
  expect_identical(model$item_clustering, prior)
  expect_output(print(strength), "latent_strength_prior")
  expect_output(print(prior), "Gnedin")
  expect_output(print(model), "Plackett--Luce")
  expect_s3_class(control, "btsbm_mcmc_control")
  expect_equal(control$iter, 50L)
  expect_error(gnedin_prior(gnedin_hyperparameter = 1), "hyperparameter")
  expect_error(mcmc_control(iter = 10, warmup = 10), "warmup")
})

test_that("all clustering-prior choices have intuitive parameter names", {
  gnedin <- gnedin_prior(gnedin_hyperparameter = 0.6)
  dp <- dirichlet_process_prior(concentration = 2)
  py <- pitman_yor_prior(concentration = 2, discount = 0.25)
  finite <- finite_partition(max_clusters = 3, concentration = 0.75)

  expect_equal(gnedin$family, "gnedin")
  expect_equal(gnedin$gnedin_hyperparameter, 0.6)
  expect_equal(dp$family, "dirichlet_process")
  expect_equal(dp$concentration, 2)
  expect_equal(py$discount, 0.25)
  expect_equal(finite$max_clusters, 3L)
  expect_error(finite_partition(max_clusters = 0), "max_clusters")
  expect_error(pitman_yor_prior(concentration = 1, discount = 1), "discount")
})

test_that("toy ranking datasets are self-contained and carry truth", {
  sushi <- sushi_toy(n_rankings = 12, rank_length = 4, seed = 4)
  eurovision <- eurovision_toy(n_rankings = 12, rank_length = 4, seed = 5)
  lbm <- pl_lbm_toy(n_rankings = 12, rank_length = 4, seed = 5)
  tennis <- tennis_toy(matches_per_pair = 4, seed = 6)

  expect_s3_class(sushi, "btsbm_ranking_data")
  expect_equal(dim(sushi$rankings), c(12, 4))
  expect_equal(length(attr(sushi, "truth")$ranking_cluster), 12)
  expect_s3_class(eurovision, "btsbm_ranking_data")
  expect_equal(dim(eurovision$rankings), c(12, 4))
  expect_equal(dim(attr(eurovision, "truth")$latent_strength), c(2, 8))
  expect_s3_class(lbm, "btsbm_ranking_data")
  expect_equal(dim(attr(lbm, "truth")$latent_strength), c(2, 3))
  expect_equal(length(attr(lbm, "truth")$item_cluster), 6)
  expect_s3_class(tennis, "btsbm_pairwise_data")
  expect_equal(dim(tennis$wins), c(6, 6))
  expect_equal(unname(tennis$wins + t(tennis$wins)), matrix(4, 6, 6) - diag(4, 6))
})
