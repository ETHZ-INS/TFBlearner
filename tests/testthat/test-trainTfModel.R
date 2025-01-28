test_that("Arguments check: Error if mismatch between transcription factor specified and feature matrix provided.", {
  fm <- Matrix(runif(100),ncol=10, nrow=10)
  attributes(fm)$transcription_factor <- "JUN"
  expect_error(trainBagged("YY1", fm, evalRounds=2, BPPARAM=SerialParam()))
})

test_that("Arguments check: Basic training setup",{

  fm <- Matrix(runif(1e5), ncol=10)
  motifScore <- Matrix(runif(1e4, 0.5, 10), ncol=1)
  labelCol <- Matrix(sample(c(1,0), 1e4, prob=c(0.05,0.95), replace=TRUE), ncol=1)
  cellTypeCol <- Matrix(sample(c(1,2), 1e4, replace=TRUE), ncol=1)
  countCol <- Matrix(sample(1:100, 1e4, replace=TRUE), ncol=1)
  widthCol <- Matrix(sample(50:200, 1e4, replace=TRUE), ncol=1)
  fm <- cbind(fm, motifScore, labelCol, cellTypeCol, countCol, widthCol)
  colnames(fm) <- c(paste("feature", 1:10, sep="_"),
                    "motif_YY1",
                    "contextTfFeat_label",
                    "context",
                    "total_overlaps",
                    "siteFeat_width")
  fm <- as(fm, "CsparseMatrix")
  attributes(fm)$transcription_factor <- "YY1"
  expect_no_error(trainBagged("YY1", fm, evalRounds=2, BPPARAM=SerialParam()))
})

test_that("Correct assignment of positive and negative fractions",{
  fm <- Matrix(runif(1e7), ncol=10)
  motifScore <- Matrix(runif(1e6, 0.5, 10), ncol=1)
  labelCol <- Matrix(sample(c(1,0), 1e6, prob=c(0.05,0.95), replace=TRUE), ncol=1)
  cellTypeCol <- Matrix(sample(c(1,2), 1e6, replace=TRUE), ncol=1)
  countCol <- Matrix(sample(1:100, 1e6, replace=TRUE), ncol=1)
  widthCol <- Matrix(sample(50:200, 1e6, replace=TRUE), ncol=1)
  fm <- cbind(fm, motifScore, labelCol, cellTypeCol, countCol, widthCol)
  colnames(fm) <- c(paste("feature", 1:10, sep="_"),
                    "motif_YY1",
                    "contextTfFeat_label",
                    "context",
                    "total_overlaps",
                    "siteFeat_width")
  fm <- as(fm, "CsparseMatrix")
  attributes(fm)$transcription_factor <- "YY1"

  mod <- trainBagged("YY1", fm, evalRounds=1, tuneHyperparams=FALSE,
                     BPPARAM=SerialParam())

  nPos <- mod$top_weighted_pos$params$n_pos_train
  nTot <- mod$top_weighted_pos$params$n_neg_train+nPos
  expect_equal(nPos/nTot, 0.25, tolerance=0.01)

  nPos <- mod$med_weighted_pos$params$n_pos_train
  nTot <- mod$med_weighted_pos$params$n_neg_train+nPos
  expect_equal(nPos/nTot, 0.25, tolerance=0.01)

  nPos <- mod$all_weighted_pos$params$n_pos_train
  nTot <- mod$all_weighted_pos$params$n_neg_train+nPos
  expect_equal(nPos/nTot, 0.25, tolerance=0.01)

  nPos <- mod$all_pos$params$n_pos_train
  nTot <- mod$all_pos$params$n_neg_train+nPos
  expect_equal(nPos/nTot, 0.25, tolerance=0.01)

  nPos <- mod$top_weighted_pos$params$n_pos_val1
  nTot <- mod$top_weighted_pos$params$n_neg_val1+nPos
  expect_equal(nPos/nTot, 0.01, tolerance=0.001)

  nPos <- mod$med_weighted_pos$params$n_pos_val1
  nTot <- mod$med_weighted_pos$params$n_neg_val1+nPos
  expect_equal(nPos/nTot, 0.01, tolerance=0.001)

  nPos <- mod$all_weighted_pos$params$n_pos_val1
  nTot <- mod$all_weighted_pos$params$n_neg_val1+nPos
  expect_equal(nPos/nTot, 0.01, tolerance=0.001)

  nPos <- mod$all_pos$params$n_pos_val1
  nTot <- mod$all_pos$params$n_neg_val1+nPos
  expect_equal(nPos/nTot, 0.01, tolerance=0.001)

  expect_equal(mod$top_weighted_pos$params$n_overlap_train_val1, integer(0))
  expect_equal(mod$top_weighted_pos$params$n_overlap_train_val2, integer(0))

  expect_equal(mod$med_weighted_pos$params$n_overlap_train_val1, integer(0))
  expect_equal(mod$med_weighted_pos$params$n_overlap_train_val2, integer(0))

  expect_equal(mod$all_weighted_pos$params$n_overlap_train_val1, integer(0))
  expect_equal(mod$all_weighted_pos$params$n_overlap_train_val2, integer(0))

  expect_equal(mod$all_pos$params$n_overlap_train_val1, integer(0))
  expect_equal(mod$all_pos$params$n_overlap_train_val2, integer(0))

  expect_true(mod$top_weighted_pos$params$n_pos_train<=mod$med_weighted_pos$params$n_pos_train)
  expect_true(mod$med_weighted_pos$params$n_pos_train<=mod$all_weighted_pos$params$n_pos_train)
  expect_true(mod$all_weighted_pos$params$n_pos_train<=mod$all_pos$params$n_pos_train)
})

test_that("Arguments check: Basic training - too small matrix",{
  expect_error(trainBagged("CTCF", fmTest, evalRounds=2, BPPARAM=SerialParam()))
})
