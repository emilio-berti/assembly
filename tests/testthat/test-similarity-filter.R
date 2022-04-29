test_that("Limiting similarity", {
  set.seed(12) #seed changed to throw a first error
  n <- sample(seq(10, 20), 1)
  sp <- draw_random_species(n, colnames(adirondack))
  sp <- resource_filtering(sp, adirondack)
  # and then returning a valid output
  sp <- similarity_filtering(sp, adirondack, t = .1, max.iter = 100)
  expect_length(sp, n)
})
