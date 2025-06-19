test_that("Predictions: Basic setup",{
  preds <- NULL
  expect_no_error(preds <- predictTfBinding(modTest, fmTest, sparsify=FALSE))
  expect_true(is(preds, "RangedSummarizedExperiment"))

  rangesObs <- rowRanges(preds)
  rangesObs <- rangesObs[order(rangesObs)]

  rangesExp <- rowRanges(fmTest)
  rangesExp <- rangesExp[rangesExp$context==levels(rangesExp$context)[[1]]]
  rangesExp <- rangesExp[order(rangesExp)]

  expect_equal(rangesObs, rangesExp)
  expect_contains(names(assays(preds)),paste(PREDPREFIX, MODELNAMES, sep="_"))
  expect_equal(colnames(preds), levels(rangesExp$context))

  fmTest2 <- fmTest
  metadata(fmTest2)[[TFNAMECOL]] <- "YY1"
  expect_error(predictTfBinding(modTest, fmTest2, sparsify=FALSE))
})

test_that("Predictions: sparsification",{
  preds <- NULL
  expect_no_error(preds <- predictTfBinding(modTest, fmTest, sparsify=TRUE))
})


test_that("Predictions: chunking",{
  preds <- NULL
  expect_no_error(preds <- predictTfBinding(modTest, fmTest, chunk=10))
})
