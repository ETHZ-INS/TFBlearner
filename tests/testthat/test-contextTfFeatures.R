test_that("Context-TF-features: Basic functionality", {
  experiments(maeTest)$contextTfFeat <- NULL
  maeTest <- contextTfFeatures(maeTest, tfName="CTCF",
                        features=c("Inserts", "Weighted_Inserts",
                                   "Cofactor_Inserts"))

  expect_s4_class(maeTest, "MultiAssayExperiment")
  expect_contains(names(experiments(maeTest)), "contextTfFeat")
})

test_that("Context-TF-features: Basic functionality - HDF5", {
  experiments(maeTest)$contextTfFeat <- NULL
  maeTestHdf5 <- contextTfFeatures(maeTestHdf5, tfName="CTCF",
                        features=c("Inserts", "Weighted_Inserts",
                                   "Cofactor_Inserts"))

  expect_s4_class(maeTestHdf5, "MultiAssayExperiment")
  expect_contains(names(experiments(maeTestHdf5)), "contextTfFeat")
})

test_that("Context-TF-features: Correct training context selection", {
  experiments(maeTest)$contextTfFeat <- NULL
  maeTest <- contextTfFeatures(maeTest, tfName="CTCF",
                               whichCol="OnlyTrain",
                               features=c("Inserts", "Weighted_Inserts",
                                          "Cofactor_Inserts"))

  expect_s4_class(maeTest, "MultiAssayExperiment")
  expect_equal(rownames(colData(experiments(maeTest)$contextTfFeat)),
               "K562_CTCF")
})

test_that("Context-TF-features: save pre-computed ChromVAR parameters in colData of siteFeat", {
  experiments(maeTest)$contextTfFeat <- NULL
  colData(experiments(maeTest)$siteFeat)$ChromVAR_expectations <- NULL
  colData(experiments(maeTest)$siteFeat)$ChromVAR_sub_ind <- NULL
  colData(experiments(maeTest)$siteFeat)$ChromVAR_background_peaks <- NULL
  maeTest <- contextTfFeatures(maeTest, tfName="CTCF",
                               whichCol="OnlyTrain",
                               features=c("Inserts", "ChromVAR_Scores"))

  expect_contains(colnames(colData(experiments(maeTest)$siteFeat)),
                  c("ChromVAR_sub_ind", "ChromVAR_expectations",
                    "ChromVAR_background_peaks"))
})

test_that("Context-TF-features: message that pre-computed ChromVAR parameters are used", {
  experiments(maeTest)$contextTfFeat <- NULL
  colData(experiments(maeTest)$siteFeat)$ChromVAR_expectations <- NULL
  colData(experiments(maeTest)$siteFeat)$ChromVAR_sub_ind <- NULL
  colData(experiments(maeTest)$siteFeat)$ChromVAR_background_peaks <- NULL
  maeTest <- contextTfFeatures(maeTest, tfName="CTCF",
                               whichCol="OnlyTrain",
                               features=c("Inserts", "ChromVAR_Scores"))

  msgs <- list()
  withCallingHandlers(contextTfFeatures(maeTest, tfName="CTCF",
                                        whichCol="OnlyTrain",
                                        features=c("Inserts",
                                                   "ChromVAR_Scores")),
                      message = function(msg) {msgs <<- append(msgs,
                                                               msg$message)})
  msgs <- unlist(msgs)

  expect_true("ChromVAR features have been pre-computed\n" %in% msgs)
  expect_true("ChromVAR-Activity ATAC associations have been pre-computed\n"
               %in% msgs)
})

test_that("Context-TF-features: save pre-computed ChromVAR-Activity and ATAC association in rowData of tfFeat", {
  experiments(maeTest)$contextTfFeat <- NULL
  tfName <- "CTCF"
  maeTest <- contextTfFeatures(maeTest, tfName=tfName,
                               whichCol="OnlyTrain",
                               features=c("Inserts", "ChromVAR_Scores"))

  # expected column names
  expCols <- paste("ChromVAR_ATAC",
                   gsub('_[0-9]+',"", c("Pearson","Cohen_Kappa")),
                   tfName, sep="_")
  expect_contains(colnames(rowData(experiments(maeTest)$tfFeat)), expCols)
})


# add dimensionality checks and mapping checks
