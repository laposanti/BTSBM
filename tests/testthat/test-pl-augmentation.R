test_that("C++ grouped augmentation matches the audited R reference", {
  rankings <- rbind(
    c(1L, 2L, 3L),
    c(3L, 1L, 2L),
    c(2L, 3L, 1L)
  )
  ranking_cluster <- c(1L, 2L, 1L)
  latent_strength <- rbind(
    c(2.0, 1.0, 3.0),
    c(1.5, 3.0, 0.75)
  )

  set.seed(101)
  reference <- BTSBM:::.pl_draw_grouped_augmentation_r(
    rankings, ranking_cluster, latent_strength
  )
  set.seed(101)
  accelerated <- BTSBM:::.pl_draw_grouped_augmentation(
    rankings, ranking_cluster, latent_strength
  )

  expect_equal(accelerated, reference, tolerance = 0)
})

test_that("grouped augmentation rejects invalid strength state", {
  rankings <- matrix(c(1L, 2L), nrow = 1L)
  expect_error(
    BTSBM:::pl_grouped_augmentation_cpp(rankings, 1L, matrix(c(1, 0), nrow = 1L)),
    "strictly positive"
  )
})
