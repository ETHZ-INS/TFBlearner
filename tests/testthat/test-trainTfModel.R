test_that("Arguments check: Error if mismatch between transcription factor specified and feature matrix provided.", {
  fm <- Matrix(runif(100),ncol=10, nrow=10)
  fm <- SummarizedExperiment(assays=list(features=fm))
  metadata(fm)$tf_name <- "JUN"
  expect_error(trainBagged("YY1", fm, evalRounds=2, BPPARAM=SerialParam()))
})

test_that("Arguments check: Basic training setup",{

  fm <- Matrix(runif(1e5), ncol=10)
  motifScore <- Matrix(runif(1e4, 0.5, 10), ncol=1)
  labelCol <- Matrix(sample(c(1,0), 1e4, prob=c(0.05,0.95), replace=TRUE), ncol=1)
  cellTypeCol <- Matrix(c(rep(1, 5e3), rep(2, 5e3)), ncol=1)
  countCol <- Matrix(sample(1:100, 1e4, replace=TRUE), ncol=1)
  widthCol <- Matrix(sample(50:200, 1e4, replace=TRUE), ncol=1)
  fm <- cbind(fm, motifScore, labelCol, cellTypeCol, countCol, widthCol)
  tfName <- "YY1"
  annoCol <- "context"
  colnames(fm) <- c(paste("feature", 1:10, sep="_"),
                    paste(assocMotifPrefix, tfName, sep="_"),
                    paste(contextTfFeat, labelFeatName, sep="_"),
                    annoCol,
                    totalOverlapsName,
                    paste(siteFeat, widthFeatName, sep="_"))
  fm <- as(fm, "CsparseMatrix")
  coords <- refCoords[sample(1:length(refCoords),
                  nrow(fm), replace=TRUE)]
  coords@elementMetadata[[annoCol]] <- cellTypeCol
  fm <- SummarizedExperiment(assays=list(features=fm),
                             rowRanges=coords)
  metadata(fm)[[annoCol]] <- unique(cellTypeCol[,1,drop=TRUE])
  metadata(fm)$tf_name <- tfName
  mods <- NULL
  expect_no_error(mods <- trainBagged(tfName, fm, evalRounds=2, BPPARAM=SerialParam()))
  expect_no_error(trainStacked(fm, mods, stackingStrat="last", BPPARAM=SerialParam()))
  expect_no_error(trainStacked(fm, mods, stackingStrat="wLast", BPPARAM=SerialParam()))
  expect_no_error(trainStacked(fm, mods, stackingStrat="wMean",
                               subSample=100, BPPARAM=SerialParam()))
})

test_that("Correct assignment of positive and negative fractions during training",{
  fm <- Matrix(runif(1e7), ncol=10)
  motifScore <- Matrix(runif(1e6, 0.5, 10), ncol=1)
  labelCol <- Matrix(sample(c(1,0), 1e6, prob=c(0.05,0.95), replace=TRUE), ncol=1)
  cellTypeCol <- Matrix(c(rep(1, 5e5), rep(2, 5e5)), ncol=1)
  countCol <- Matrix(sample(1:100, 1e6, replace=TRUE), ncol=1)
  widthCol <- Matrix(sample(50:200, 1e6, replace=TRUE), ncol=1)
  fm <- cbind(fm, motifScore, labelCol, cellTypeCol, countCol, widthCol)
  tfName <- "YY1"
  annoCol <- "context"
  colnames(fm) <- c(paste("feature", 1:10, sep="_"),
                    paste(assocMotifPrefix, tfName, sep="_"),
                    paste(contextTfFeat, labelFeatName, sep="_"),
                    annoCol,
                    totalOverlapsName,
                    paste(siteFeat, widthFeatName, sep="_"))
  fm <- as(fm, "CsparseMatrix")
  coords <- refCoords[sample(1:length(refCoords),
                             nrow(fm), replace=TRUE)]
  coords@elementMetadata[[annoCol]] <- cellTypeCol
  fm <- SummarizedExperiment(assays=list(features=fm),
                             rowRanges=coords)
  metadata(fm)[[annoCol]] <- unique(cellTypeCol[,1,drop=TRUE])
  metadata(fm)$tf_name <- tfName

  mod <- trainBagged(tfName, fm, evalRounds=1, tuneHyperparams=FALSE,
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

  expect_true(mod$top_weighted_pos$params$n_pos_train<=mod$med_weighted_pos$params$n_pos_train)
  expect_true(mod$med_weighted_pos$params$n_pos_train<=mod$all_weighted_pos$params$n_pos_train)
  expect_true(mod$all_weighted_pos$params$n_pos_train<=mod$all_pos$params$n_pos_train)
})
