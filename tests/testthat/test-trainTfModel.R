test_that("Arguments check: Error if mismatch between transcription factor specified and feature matrix provided.", {
  fm <- Matrix(runif(100),ncol=10, nrow=10)
  fm <- SummarizedExperiment(assays=list(features=fm))
  metadata(fm)$tf_name <- "JUN"
  expect_error(trainTfModel("YY1", fm, evalRounds=2, stackingStrat="last",
                            BPPARAM=SerialParam()))
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
                    MOTIFFEATCOLNAME,
                    LABELCOLNAME,
                    annoCol,
                    COUNTCOLNAME,
                    paste(SITEFEAT, WIDTHFEATNAME, sep="_"))
  fm <- as(fm, "CsparseMatrix")
  coords <- refCoords[sample(1:length(refCoords),
                      floor(nrow(fm)/2), replace=TRUE)]
  coords <- as.data.table(coords)
  coords[,seqnames:=sample(paste0("chr", 1:20), .N, replace=TRUE)]
  coords <- makeGRangesFromDataFrame(as.data.frame(coords),
                                     keep.extra.columns=TRUE)
  coords <- c(coords, coords)
  coords@elementMetadata[[annoCol]] <- cellTypeCol[,1,drop=TRUE]
  fm <- SummarizedExperiment(assays=list(features=fm),
                             rowRanges=coords)
  metadata(fm)[[annoCol]] <- unique(cellTypeCol[,1,drop=TRUE])
  metadata(fm)[[TFNAMECOL]] <- tfName
  mods <- NULL
  modsBaggedNames <- c(MODELTOPWEIGHTNAME,
                       MODELMEDWEIGHTNAME,
                       MODELALLWEIGHTNAME,
                       MODELALLNAME)

  expect_no_error(mods <- trainTfModel(tfName, fm, evalRounds=2,
                                       stackingStrat=c("wMean"),
                                       subSample=100,
                                       valChrs=c("chr9", "chr4"),
                                       BPPARAM=SerialParam()))
  expect_contains(names(mods), c(modsBaggedNames,
                                 paste(MODELSTACKEDSUFFIX, "wMean", sep="_"),
                                 "stacking_strategy"))
  expect_no_error(.trainStacked(fm, mods[modsBaggedNames],
                                stackingStrat="wLast"))
  expect_no_error(.trainStacked(fm, mods[modsBaggedNames],
                                stackingStrat="last"))

  # test saving functionality
  outDir <- tempdir()
  modFilePath <- file.path(outDir, "testModels.txt")
  saveModels(mods, filePath=modFilePath)
  expect_true(file.exists(modFilePath))
  modLoad <- loadModels(modFilePath)
  expect_contains(names(modLoad), c(modsBaggedNames, "stacking_strategy"))
  expect_equal(mods[[MODELMEDWEIGHTNAME]]$sparse_thr,
               modLoad[[MODELMEDWEIGHTNAME]]$sparse_thr)
  expect_equal(mods[[MODELMEDWEIGHTNAME]]$tf,
               modLoad[[MODELMEDWEIGHTNAME]]$tf)
  expect_equal(mods[[MODELMEDWEIGHTNAME]]$stacking_weights,
               modLoad[[MODELMEDWEIGHTNAME]]$stacking_weights)

  # test saved thresholds
  expect_equal(mods[[DICHOTTHRESH]], modLoad[[DICHOTTHRESH]])
  expect_equal(mods[[RECENTRY]], modLoad[[RECENTRY]])
  expect_equal(mods[[PRENTRY]], modLoad[[PRENTRY]])
  expect_equal(mods[[AUCPRENTRY]], modLoad[[AUCPRENTRY]])

  file.remove(modFilePath)
})


test_that("Saving and loading motif information",{

  # check that contained in models
  expect_contains(names(modTest[[MODELALLNAME]]$params),
                  c(PRESELMOTIFCOL, PRESELACTCOL))
  tfMotifName <- paste(TFMOTIFPREFIX, 1, sep="_")
  expect_equal(modTest[[MODELALLNAME]]$params[[PRESELMOTIFCOL]][[tfMotifName]],
               modTest[[MODELALLNAME]]$params$tf)
  expect_contains(names(modTest[[MODELTOPWEIGHTNAME]]$params),
                  c(PRESELMOTIFCOL, PRESELACTCOL))
  expect_equal(modTest[[MODELTOPWEIGHTNAME]]$params[[PRESELMOTIFCOL]][[tfMotifName]],
               modTest[[MODELTOPWEIGHTNAME]]$params$tf)

  # check saving
  outDir <- tempdir()
  modFilePath <- file.path(outDir, "testModels.txt")
  saveModels(modTest, filePath=modFilePath)
  expect_true(file.exists(modFilePath))

  # check correct loading
  modLoad <- loadModels(filePath=modFilePath)
  expect_contains(names(modLoad[[MODELALLNAME]]$params),
                  c(PRESELMOTIFCOL, PRESELACTCOL))
  tfMotifName <- paste(TFMOTIFPREFIX, 1, sep="_")
  expect_equal(modLoad[[MODELALLNAME]]$params[[PRESELMOTIFCOL]][[tfMotifName]],
               modLoad[[MODELALLNAME]]$params$tf)
  expect_contains(names(modLoad[[MODELTOPWEIGHTNAME]]$params),
                  c(PRESELMOTIFCOL, PRESELACTCOL))
  expect_equal(modLoad[[MODELTOPWEIGHTNAME]]$params[[PRESELMOTIFCOL]][[tfMotifName]],
               modLoad[[MODELTOPWEIGHTNAME]]$params$tf)

  file.remove(modFilePath)
})

test_that("Saving and loading maintains functioning model",{
  outDir <- tempdir()
  modFilePath <- file.path(outDir, "testModels.txt")
  saveModels(modTest, filePath=modFilePath)

  # check correct loading
  modLoad <- loadModels(filePath=modFilePath)
  expect_no_error(preds <- predictTfBinding(modLoad, fmTest, sparsify=FALSE))
  file.remove(modFilePath)
})

test_that("Saving and loading package version",{
  expect_contains(names(modTest[[MODELALLNAME]]$params), "package_version")
  expVersion <- .getPackageVersion()
  expect_equal(modTest[[MODELALLNAME]]$params$package_version, expVersion)

  # check saving
  outDir <- tempdir()
  modFilePath <- file.path(outDir, "testModels.txt")
  saveModels(modTest, filePath=modFilePath)

  # check correct loading
  modLoad <- loadModels(filePath=modFilePath)
  expect_contains(names(modLoad[[MODELALLNAME]]$params), "package_version")
  expect_equal(modLoad[[MODELALLNAME]]$params$package_version, expVersion)
  expect_contains(names(modLoad[[MODELTOPWEIGHTNAME]]$params),"package_version")
  expect_contains(modLoad[[MODELTOPWEIGHTNAME]]$params, expVersion)

  file.remove(modFilePath)
})

test_that("Correct assignment of positive and negative fractions during training",{
  fm <- Matrix(runif(1e6), ncol=10)
  motifScore <- Matrix(runif(1e5, 0.5, 10), ncol=1)
  labelCol <- Matrix(sample(c(1,0), 1e5, prob=c(0.05,0.95), replace=TRUE), ncol=1)
  cellTypeCol <- Matrix(c(rep(1, 5e4), rep(2, 5e4)), ncol=1)
  countCol <- Matrix(sample(1:100, 1e5, replace=TRUE), ncol=1)
  widthCol <- Matrix(sample(50:200, 1e5, replace=TRUE), ncol=1)
  fm <- cbind(fm, motifScore, labelCol, cellTypeCol, countCol, widthCol)
  tfName <- "YY1"
  annoCol <- "context"
  colnames(fm) <- c(paste("feature", 1:10, sep="_"),
                    MOTIFFEATCOLNAME,
                    LABELCOLNAME,
                    annoCol,
                    COUNTCOLNAME,
                    paste(SITEFEAT, WIDTHFEATNAME, sep="_"))
  fm <- as(fm, "CsparseMatrix")
  coords <- refCoords[sample(1:length(refCoords),
                             floor(nrow(fm)/2), replace=TRUE)]
  coords <- c(coords, coords)
  coords@elementMetadata[[annoCol]] <- cellTypeCol[,1,drop=TRUE]
  fm <- SummarizedExperiment(assays=list(features=fm),
                             rowRanges=coords)
  metadata(fm)[[annoCol]] <- unique(cellTypeCol[,1,drop=TRUE])
  metadata(fm)[[TFNAMECOL]] <- tfName

  mod <- trainTfModel(tfName, fm, evalRounds=1,
                      tuneHyperparams=FALSE, stackingStrat="last",
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
  expect_equal(nPos/nTot, 0.01, tolerance=0.01)

  nPos <- mod$med_weighted_pos$params$n_pos_val1
  nTot <- mod$med_weighted_pos$params$n_neg_val1+nPos
  expect_equal(nPos/nTot, 0.01, tolerance=0.01)

  nPos <- mod$all_weighted_pos$params$n_pos_val1
  nTot <- mod$all_weighted_pos$params$n_neg_val1+nPos
  expect_equal(nPos/nTot, 0.01, tolerance=0.01)

  nPos <- mod$all_pos$params$n_pos_val1
  nTot <- mod$all_pos$params$n_neg_val1+nPos
  expect_equal(nPos/nTot, 0.01, tolerance=0.01)

  expect_true(mod$top_weighted_pos$params$n_pos_train<=mod$med_weighted_pos$params$n_pos_train)
  expect_true(mod$med_weighted_pos$params$n_pos_train<=mod$all_weighted_pos$params$n_pos_train)
  expect_true(mod$all_weighted_pos$params$n_pos_train<=mod$all_pos$params$n_pos_train)
})


test_that("Sampling of additional points for training (when no available)",{
  posFracExp <- 0.01
  n <- 1e5
  set <- 1:n
  weights <- sample(runif(5,0.001,1), n, replace=TRUE)
  weights <- weights[order(weights, decreasing=TRUE)]
  labels <- sample(c(1,0), n, replace=TRUE, prob=c(0.1,0.9))

  subSet <- set[1:floor(n/10)]
  valSet <-  sample(subSet, floor(length(subSet)*0.15))
  trainSet <- setdiff(subSet, valSet)

  valAdd <- .addInst(weights, labels, subSet, setdiff(set, c(valSet, trainSet)),
                     posFrac=posFracExp, nPos=100)
  trainAdd <- .addInst(weights, labels, subSet, setdiff(set, c(valSet,
                                                               valAdd,
                                                               trainSet)),
                    posFrac=posFracExp, nPos=100)
  trainSet <- c(trainSet, trainAdd)
  trainSub <- .ensureFrac(labels[trainSet], weights[trainSet],
                          posFrac=posFracExp)
  trainSet <- trainSet[trainSub]

  valSet <- c(valSet, valAdd)
  valSub <- .ensureFrac(labels[valSet], weights[valSet],
                        posFrac=posFracExp)
  valSet <- valSet[valSub]

  expect_length(intersect(trainSet, valSet), 0)
  expect_equal(sum(labels[valSet]==1)/length(valSet),
               posFracExp, tolerance=0.01)
  expect_equal(sum(labels[trainSet]==1)/length(trainSet),
               posFracExp, tolerance=0.01)
})


test_that("Sampling of additional points for training when available",{
  posFracExp <- 0.01
  n <- 1e5
  set <- 1:n
  weights <- sample(seq(0.5,1,length.out=5), n, replace=TRUE)
  weights <- weights[order(weights, decreasing=TRUE)]
  labels <- sample(c(1,0), n, replace=TRUE, prob=c(0.1,0.9))

  subSet <- set[1:floor(n/10)]
  valSet <-  sample(subSet, floor(length(subSet)*0.15))
  trainSet <- setdiff(subSet, valSet)

  nPosAdd <- 100
  valAdd <- .addInst(weights, labels, subSet, setdiff(set, c(valSet, trainSet)),
                     posFrac=posFracExp, nPos=nPosAdd)
  trainAdd <- .addInst(weights, labels, subSet, setdiff(set, c(valSet,
                                                               valAdd,
                                                               trainSet)),
                       posFrac=posFracExp, nPos=nPosAdd)
  trainSetOrig <- trainSet
  trainSet <- c(trainSetOrig, trainAdd)
  trainSub <- .ensureFrac(labels[trainSet], weights[trainSet],
                          posFrac=posFracExp)
  trainSet <- trainSet[trainSub]

  valSetOrig <- valSet
  valSet <- c(valSet, valAdd)
  valSub <- .ensureFrac(labels[valSet], weights[valSet],
                        posFrac=posFracExp)
  valSet <- valSet[valSub]

  expect_length(intersect(trainSet, valSet), 0)
  expect_true(length(valSetOrig)<length(valSet))
  expect_true(length(trainSetOrig)<length(trainSet))
  expect_true(sum(labels[valAdd]==1)<=nPosAdd)
  expect_true(sum(labels[trainAdd]==1)<=nPosAdd)
  expect_equal(sum(labels[valAdd]==1)/length(valAdd),
               posFracExp, tolerance=0.001)
  expect_equal(sum(labels[trainAdd]==1)/length(trainAdd),
               posFracExp, tolerance=0.01)
  expect_equal(sum(labels[valSet]==1)/length(valSet),
               posFracExp, tolerance=0.01)
  expect_equal(sum(labels[trainSet]==1)/length(trainSet),
               posFracExp, tolerance=0.01)
})


# do the same weighted
test_that("Sampling of additional points for training unweighted",{
  posFracExp <- 0.01
  n <- 1e5
  set <- 1:n
  weights <- rep(1,n)
  labels <- sample(c(1,0), n, replace=TRUE, prob=c(0.1,0.9))

  subSet <- set[1:floor(n/10)]
  valSet <-  sample(subSet, floor(length(subSet)*0.15))
  trainSet <- setdiff(subSet, valSet)

  valAdd <- .addInst(weights, labels, subSet, setdiff(set, c(valSet, trainSet)),
                     posFrac=posFracExp, nPos=100)
  trainAdd <- .addInst(weights, labels, subSet, setdiff(set, c(valSet,
                                                               valAdd,
                                                               trainSet)),
                       posFrac=posFracExp, nPos=100)
  trainSetOrig <- trainSet
  trainSet <- c(trainSetOrig, trainAdd)
  trainSub <- .ensureFrac(labels[trainSet], weights[trainSet],
                          posFrac=posFracExp)
  trainSet <- trainSet[trainSub]

  valSetOrig <- valSet
  valSet <- c(valSet, valAdd)
  valSub <- .ensureFrac(labels[valSet], weights[valSet],
                        posFrac=posFracExp)
  valSet <- valSet[valSub]

  expect_length(intersect(trainSet, valSet), 0)
  expect_true(length(valSetOrig)<length(valSet))
  expect_true(length(trainSetOrig)<length(trainSet))
  expect_equal(sum(labels[valAdd]==1)/length(valAdd),
               posFracExp, tolerance=0.001)
  expect_equal(sum(labels[trainAdd]==1)/length(trainAdd),
               posFracExp, tolerance=0.01)
  expect_equal(sum(labels[valSet]==1)/length(valSet),
               posFracExp, tolerance=0.01)
  expect_equal(sum(labels[trainSet]==1)/length(trainSet),
               posFracExp, tolerance=0.01)
})
