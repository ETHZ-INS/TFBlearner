test_that("Feature Matrix: Basic functionality", {
  fm <- getFeatureMatrix(maeTest, tfName="CTCF",
                         addLabels=FALSE,
                         saveHdf5=FALSE)
  contexts <- getContexts(maeTest, tfName="CTCF", which="ATAC")

  expect_s4_class(assays(fm)$features, "CsparseMatrix")
  expect_equal(nrow(fm), length(contexts)*length(example_coords))
})

test_that("Feature Matrix: Basic functionality - HDF5", {
  fm <- getFeatureMatrix(maeTestHdf5, tfName="CTCF",
                         addLabels=FALSE,
                         saveHdf5=FALSE)
  contexts <- getContexts(maeTest, tfName="CTCF", which="ATAC")

  expect_s4_class(assays(fm)$features, "CsparseMatrix")
  expect_equal(nrow(fm), length(contexts)*length(example_coords))
})

test_that("Feature Matrix: Basic functionality - Saving to HDF5", {
  outDir <- tempdir()
  prefix <- "test"
  fm <- getFeatureMatrix(maeTest, tfName="CTCF",
                                  addLabels=FALSE,
                                  saveHdf5=TRUE,
                                  outDir=outDir,
                                  prefix=prefix)

  fmFilePath <- file.path(outDir,
                          paste(prefix, "feature_matrix", "CTCF", sep="_"))
  fmFilePath <- paste0(fmFilePath, ".h5")
  expect_s4_class(assays(fm)$features, "DelayedMatrix")
  expect_true(file.exists(fmFilePath))
  file.remove(fmFilePath)
})

test_that("Feature Matrix: Correct context selection - only for training contexts", {
  fm <- getFeatureMatrix(maeTest, tfName="CTCF",
                         whichCol="OnlyTrain",
                         addLabels=FALSE,
                         saveHdf5=FALSE)

  testContext <- subset(colData(maeTest), is_testing)$context
  contexts <- getContexts(maeTest, tfName="CTCF", which="ATAC")
  trainContexts <- setdiff(contexts, testContext)

  expect_equal(nrow(fm), length(trainContexts)*length(example_coords))
  expect_equal(rowRanges(fm),
               rep(rowRanges(maeTest[["ATAC"]]), length(trainContexts)))
  expect_equal(metadata(fm)$cellular_contexts, trainContexts)
})

test_that("Feature Matrix: Correct context selection - only for specified context", {
  fm <- getFeatureMatrix(maeTest, tfName="CTCF",
                         whichCol="Col",
                         colSel="A549",
                         addLabels=FALSE,
                         saveHdf5=FALSE)

  expect_equal(nrow(fm), length(example_coords))
  expect_equal(rowRanges(fm), rowRanges(maeTest[["ATAC"]]))
  expect_equal(metadata(fm)$cellular_contexts, "A549")
  expect_equal(metadata(fm)$cofactors, "JUN")
  expect_contains(names(metadata(fm)), "associated_motifs")
})

test_that("Feature Matrix: Column names corresponding to R conventions", {
  fm <- getFeatureMatrix(maeTest, tfName="CTCF",
                         whichCol="Col",
                         colSel="A549",
                         addLabels=FALSE,
                         saveHdf5=FALSE)

  # should pass anyways given current setup script, for later
  expect_equal(colnames(fm), make.names(colnames(fm), unique=TRUE))
})

test_that("Feature Matrix: Correct context selection - no context provided", {
  expect_error(getFeatureMatrix(maeTest, tfName="CTCF",
                                whichCol="Col",
                                addLabels=FALSE,
                                saveHdf5=FALSE))
})

test_that("Feature Matrix: warning if features have not been computed for provided cellular contexts", {
  experiments(maeTest)$contextTfFeat <- NULL
  maeTest <- contextTfFeatures(maeTest, tfName="CTCF",
                               whichCol="Col",
                               colSel="K562",
                               features=c("Inserts"))

  expect_warning(getFeatureMatrix(maeTest, tfName="CTCF",
                                  whichCol="Col",
                                  colSel=c("K562", "A549"),
                                  addLabels=FALSE,
                                  saveHdf5=FALSE),
                 "Not all the cellular-contexts requested have contextTfFeats compiled.\n
                   Missing are: A549")
})
