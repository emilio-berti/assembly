test_that("Test wJaccard Cpp", {
  A <- matrix(runif(100), 10, 10)
  sim <- matrix(0, 10, 10)
  for (i in 1:nrow(sim)) {
    for (j in 1:ncol(sim)) {
      Min <- apply(A[, c(i, j)], MARGIN = 1, min)
      Max <- apply(A[, c(i, j)], MARGIN = 1, max)
      sim[i, j] <- sum(Min) / sum(Max)
      sim[j, i] <- sim[i, j]
    }
  }
  diag(sim) <- 1
  w <- wJaccard(A)$Jaccard
  expect_identical(sim, w)
})
