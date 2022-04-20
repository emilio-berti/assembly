test_that("Hidden functions", {
  set.seed(1234)
  fw <- matrix(as.numeric(runif(1e2) > .9), 1e1, 1e1)
  fw <- name_metaweb(fw)
  expect_identical(rownames(fw)[colSums(fw) == 0],
                   assembly:::.basals(fw))
  expect_identical(rownames(fw)[colSums(fw) > 0],
                   assembly:::.consumers(fw))
  expect_identical(rownames(fw)[rowSums(fw) == 0 & colSums(fw) > 0],
                   assembly:::.top(fw))
  expect_identical(rownames(fw)[rowSums(fw) == 0 & colSums(fw) == 0],
                   assembly:::.find_isolated(colnames(fw), fw))
})
