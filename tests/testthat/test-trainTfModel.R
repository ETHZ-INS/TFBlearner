test_that("Arguments check: Error if mismatch between transcription factor specified and feature matrix provided.", {
  fm <- Matrix(runif(100),ncol=10, nrow=10)
  attributes(fm)$transcription_factor <- "JUN"
  expect_error(trainBagged("YY1", fm, evalRounds=2, BPPARAM=SerialParam()))
})

test_that("Arguments check: Basic training setup",{

  fm <- Matrix(runif(1e7), ncol=10)
  motifScore <- Matrix(runif(1e6, 0.5, 10), ncol=1)
  labelCol <- Matrix(sample(c(1,0), 1e6, prob=c(0.05,0.95), replace=TRUE), ncol=1)
  cellTypeCol <- Matrix(sample(c(1,2), 1e6, replace=TRUE), ncol=1)
  countCol <- Matrix(sample(1:100, 1e6, replace=TRUE), ncol=1)
  widthCol <- Matrix(sample(50:200, 1e6, replace=TRUE), ncol=1)
  fm <- cbind(fm, motifScore, labelCol, cellTypeCol, countCol, widthCol)
  colnames(fm) <- c(paste("feature", 1:10, sep="_"),
                    "motif_YY1",
                    "contextFeat_label",
                    "context",
                    "total_overlaps",
                    "siteFeat_width")
  fm <- as(fm, "CsparseMatrix")
  attributes(fm)$transcription_factor <- "YY1"
  expect_no_error(trainBagged("YY1", fm, evalRounds=2, BPPARAM=SerialParam()))
})

#test_that("Arguments check: Basic training setup - small feature matrix",{})
