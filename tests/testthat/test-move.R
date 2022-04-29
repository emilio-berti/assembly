test_that("Metropolis-Hastings", {
  expect_true(metropolis.hastings(1, 1, 0))
  expect_false(metropolis.hastings(1, 2, 0))
  expect_type(metropolis.hastings(1, 1, 1), "logical")
})

test_that("One move and t = 0", {
  set.seed(1234)
  n <- sample(seq(10, 50), 1)
  sp <- draw_random_species(n, colnames(adirondack))
  sp <- resource_filtering(sp, adirondack)
  sp <- assembly:::.move(sp, adirondack)
  expect_length(sp, n)
  # force isolated species and expect error
  fw <- adirondack
  fw[sp[1], ] <- 0
  fw[, sp[1]] <- 0
  expect_error(similarity_filtering(sp, t = 0.1, fw))
  expect_identical(sp, assembly:::.move(sp, fw))
})
