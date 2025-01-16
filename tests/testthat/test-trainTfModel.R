test_that("Arguments check: Error if mismatch between transcription factor specified and feature matrix provided.", {
  fm <- Matrix(runif(100),ncol=10, nrow=10)
  attributes(fm)$transcription_factor <- "JUN"
  expect_error(trainBagged("YY1", fm))
})
