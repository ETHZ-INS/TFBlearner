test_that("Context-TF-features: Basic functionality", {
  experiments(maeTest2)[[contextTfFeat]] <- NULL
  maeTest2 <- contextTfFeatures(maeTest2, tfName="CTCF",
                                features=c("Inserts", "Weighted_Inserts",
                                           "ChromVAR_Scores"))

  expect_s4_class(maeTest2, "MultiAssayExperiment")
  expect_contains(names(experiments(maeTest2)), contextTfFeat)
})

test_that("Context-TF-features: Basic functionality - HDF5", {
  experiments(maeTest)[[contextTfFeat]] <- NULL
  maeTestHdf5 <- contextTfFeatures(maeTestHdf5, tfName="CTCF",
                        features=c("Inserts", "Weighted_Inserts"))

  expect_s4_class(maeTestHdf5, "MultiAssayExperiment")
  expect_contains(names(experiments(maeTestHdf5)), contextTfFeat)
})

test_that("Context-TF-features: Correct training context selection", {
  experiments(maeTest)[[contextTfFeat]] <- NULL
  maeTest <- contextTfFeatures(maeTest, tfName="CTCF",
                               whichCol="OnlyTrain",
                               features=c("Inserts", "Weighted_Inserts"))

  expect_s4_class(maeTest, "MultiAssayExperiment")
  expect_equal(rownames(colData(maeTest[[contextTfFeat]])), "K562_CTCF")
})

test_that("Assays are preserved when computing for new TF", {
  assayNamesOrig <- names(assays(maeTest[[contextTfFeat]]))
  maeTest <- tfFeatures(maeTest, tfName="JUN",
                        features=c("CTCF", "MAX"))
  maeTest <- contextTfFeatures(maeTest, tfName="JUN",
                               features=c("Inserts", "Weighted_Inserts"),
                               addLabels=TRUE)
  assayNamesNew <- names(assays(maeTest[[contextTfFeat]]))
  expect_equal(assayNamesNew, assayNamesOrig)
})

test_that("Error if features have not been computed for provided TF", {
  tfName="JUN"
  expect_error(contextTfFeatures(maeTest, tfName=tfName))
})
