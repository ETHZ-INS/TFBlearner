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

test_that("Context-TF-features: message that pre-computed parameters are used", {
  experiments(maeTest)$contextTfFeat <- NULL
  colData(experiments(maeTest)$siteFeat)$ChromVAR_expectations <- NULL
  colData(experiments(maeTest)$siteFeat)$ChromVAR_sub_ind <- NULL
  colData(experiments(maeTest)$siteFeat)$ChromVAR_background_peaks <- NULL
  maeTest <- contextTfFeatures(maeTest, tfName="CTCF",
                               whichCol="All",
                               addLabels=FALSE,
                               features=c("Inserts", "ChromVAR_Scores",
                                          "MDS_Context"))

  msgs <- list()
  withCallingHandlers(contextTfFeatures(maeTest, tfName="CTCF",
                                        whichCol="All",
                                        addLabels=FALSE,
                                        features=c("Inserts",
                                                   "ChromVAR_Scores",
                                                   "MDS_Context")),
                      message = function(msg) {msgs <<- append(msgs,
                                                               msg$message)})
  msgs <- unlist(msgs)

  expect_true("ChromVAR features have been pre-computed\n" %in% msgs)
  expect_true("ChromVAR-Activity ATAC associations have been pre-computed\n"
               %in% msgs)
  expect_true("Using pre-computed MDS-dimensions\n" %in% msgs)
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


test_that("Context-TF-features: save pre-computed MDS-dimensions in colData of ATAC", {
  experiments(maeTest)$contextTfFeat <- NULL
  tfName <- "CTCF"
  maeTest <- contextTfFeatures(maeTest, tfName=tfName,
                               whichCol="All",
                               addLabels=FALSE,
                               features=c("Inserts",
                                          "ChromVAR_Scores",
                                          "MDS_Context"))

  # expected column names
  expect_contains(colnames(colData(experiments(maeTest)$ATAC)),
                  c("MDS_Context_1", "MDS_Context_2"))
})

# add dimensionality checks and mapping checks
