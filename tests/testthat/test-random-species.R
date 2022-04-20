test_that("Random species correct length", {
  set.seed(1234)
  n <- sample(seq(10, 50), 1)
  sp <- draw_random_species(n, colnames(adirondack))
  expect_length(sp, n)
})
