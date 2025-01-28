test_that("Feature Matrix: Basic functionality", {
  fm <- getFeatureMatrix(maeTest, tfName="CTCF",
                         addLabels=FALSE,
                         saveHdf5=FALSE,
                         BPPARAM=SerialParam())
  contexts <- getContexts(maeTest, tfName="CTCF")

  expect_s4_class(fm, "CsparseMatrix")
  expect_equal(nrow(fm), length(contexts)*length(example_coords))
})

test_that("Feature Matrix: Basic functionality - HDF5", {
  fm <- getFeatureMatrix(maeTestHdf5, tfName="CTCF",
                         addLabels=FALSE,
                         saveHdf5=FALSE,
                         BPPARAM=SerialParam())
  contexts <- getContexts(maeTest, tfName="CTCF")

  expect_s4_class(fm, "CsparseMatrix")
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
  expect_s4_class(fm, "DelayedMatrix")
  expect_true(file.exists(fmFilePath))
  file.remove(fmFilePath)
})

test_that("Feature Matrix: Basic functionality - Saving to HDF5", {
  fm <- getFeatureMatrix(maeTest, tfName="CTCF",
                         which="OnlyTrain",
                         addLabels=FALSE,
                         saveHdf5=FALSE)

  testContext <- subset(colData(maeTest), is_testing)$context
  contexts <- getContexts(maeTest, tfName="CTCF")
  trainContexts <- setdiff(contexts, testContext)

  expect_equal(nrow(fm), length(trainContexts)*length(example_coords))
  expect_equal(attributes(fm)$cellular_contexts, trainContexts)
})
