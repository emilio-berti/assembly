test_that("Resource filtering: no isolated species", {
  set.seed(1234)
  n <- sample(seq(10, 50), 1)
  sp <- draw_random_species(n, colnames(adirondack))
  sp <- resource_filtering(sp, adirondack, keep.n.basal = FALSE)
  fw <- adirondack[sp, sp]
  basals <- assembly:::.basals(adirondack)
  no_prey <- setdiff(colnames(fw)[colSums(fw) == 0], basals)
  consumers <- assembly:::.consumers(adirondack)
  no_cons <- setdiff(colnames(fw)[rowSums(fw) == 0], consumers)
  expect_length(no_prey, 0)
  expect_length(no_cons, 0)
})

test_that("Resource filtering: keep number of basal species", {
  set.seed(1234)
  n <- sample(seq(10, 50), 1)
  sp <- draw_random_species(n, colnames(adirondack))
  start_basals <- intersect(assembly:::.basals(adirondack), sp)
  sp <- resource_filtering(sp, adirondack, keep.n.basal = TRUE)
  end_basals <- intersect(assembly:::.basals(adirondack), sp)
  expect_equal(length(start_basals), length(end_basals))
})
