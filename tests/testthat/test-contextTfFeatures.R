# to muffle the warning "In cor(...): the standard deviation is zero"
# when computing ChromVAR-ATAC association.
# This warning is supposed to appear in the given test setup but
# is not informative for the test cases.
suppressSdWarning <- function(fun, args){

  msg <- "the standard deviation is zero"
  withCallingHandlers(
          res <- do.call(fun, args),
          warning=function(w){
          if(grepl(msg, conditionMessage(w))){
            invokeRestart("muffleWarning")}
  })
  return(res)
}

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
  colData(experiments(maeTest)$siteFeat)$top_var_sites <- NULL
  maeTest <- suppressSdWarning(contextTfFeatures, list(mae=maeTest,
                                                       tfName="CTCF",
                                                       whichCol="OnlyTrain",
                                                       features=c("Inserts",
                                                            "ChromVAR_Scores")))

  expect_contains(colnames(colData(experiments(maeTest)$siteFeat)),
                  c("ChromVAR_expectations",
                    "ChromVAR_background_peaks",
                    "top_var_sites"))
})

test_that("Context-TF-features: message that pre-computed parameters/features are used", {
  experiments(maeTest2)$contextTfFeat <- NULL
  colData(experiments(maeTest2)$siteFeat)$ChromVAR_expectations <- NULL
  colData(experiments(maeTest2)$siteFeat)$top_var_sites <- NULL
  colData(experiments(maeTest2)$siteFeat)$ChromVAR_background_peaks <- NULL

  maeTest2 <- suppressSdWarning(contextTfFeatures, list(mae=maeTest2,
                                                        tfName="CTCF",
                                                        whichCol="All",
                                                        addLabels=FALSE,
                                                        features=c("Inserts",
                                                            "ChromVAR_Scores",
                                                            "MDS_Context",
                                                            "Max_ATAC_Signal")))

  msgs <- list()
  withCallingHandlers(contextTfFeatures(maeTest2, tfName="CTCF",
                                        whichCol="All",
                                        addLabels=FALSE,
                                        features=c("Inserts",
                                                   "ChromVAR_Scores",
                                                   "MDS_Context",
                                                   "Max_ATAC_Signal")),
                      message = function(msg) {msgs <<- append(msgs,
                                                               msg$message)})
  msgs <- unlist(msgs)

  expect_true("ChromVAR features have been pre-computed\n" %in% msgs)
  expect_true("ChromVAR-Activity ATAC associations have been pre-computed\n"
               %in% msgs)
  expect_true("Using pre-computed MDS-dimensions\n" %in% msgs)
  expect_true("Using pre-computed maximal ATAC signals\n" %in% msgs)
})

test_that("Context-TF-features: save pre-computed ChromVAR-Activity and ATAC association in rowData of tfFeat", {
  experiments(maeTest)$contextTfFeat <- NULL
  tfName <- "CTCF"

  maeTest <- suppressSdWarning(contextTfFeatures, list(mae=maeTest,
                                                       tfName=tfName,
                                                       whichCol="OnlyTrain",
                                                       features=c("Inserts",
                                                            "ChromVAR_Scores")))

  # expected column names
  expCols <- paste("ChromVAR_ATAC",
                   gsub('_[0-9]+',"", c("Pearson","Cohen_Kappa")),
                   tfName, sep="_")
  expect_contains(colnames(rowData(experiments(maeTest)$tfFeat)), expCols)
})


test_that("Context-TF-features: save pre-computed MDS-dimensions in colData of ATAC", {
  experiments(maeTest2)$contextTfFeat <- NULL
  tfName <- "CTCF"
  maeTest2 <- contextTfFeatures(maeTest2, tfName=tfName,
                                whichCol="All",
                                addLabels=FALSE,
                                features=c("Inserts",
                                           "MDS_Context"))

  # expected column names
  expect_contains(colnames(colData(experiments(maeTest2)$ATAC)),
                  c("MDS_Context_1", "MDS_Context_2"))
})

test_that("Context-TF-features: save pre-computed maximal ATAC-signal in rowData of ATAC", {
  experiments(maeTest)$contextTfFeat <- NULL
  tfName <- "CTCF"
  maeTest <- contextTfFeatures(maeTest, tfName=tfName,
                               whichCol="All",
                               addLabels=FALSE,
                               features=c("Inserts",
                                          "Max_ATAC_Signal"))

  # expected column names
  expect_contains(colnames(rowData(experiments(maeTest)$ATAC)),
                  c("Max_ATAC_Signal"))
})

test_that("Context-TF-features: save pre-computed ATAC-signal variance in rowData of ATAC", {
  experiments(maeTest)$contextTfFeat <- NULL
  tfName <- "CTCF"
  maeTest <- contextTfFeatures(maeTest, tfName=tfName,
                               whichCol="All",
                               addLabels=FALSE,
                               features=c("Inserts",
                                          "ATAC_Variance"))

  # expected column names
  expect_contains(colnames(rowData(experiments(maeTest)$ATAC)),
                  c("ATAC_Variance"))
})

test_that("Assays are preserved when computing for new TF", {
  assayNamesOrig <- names(assays(maeTest[["contextTfFeat"]]))
  maeTest <- tfFeatures(maeTest, tfName="JUN",
                        features=c("CTCF", "MAX"))
  maeTest <- contextTfFeatures(maeTest, tfName="JUN",
                               features=c("Inserts", "Weighted_Inserts"),
                               addLabels=TRUE)
  assayNamesNew <- names(assays(maeTest[["contextTfFeat"]]))
  expect_equal(assayNamesNew, assayNamesOrig)
})

# add dimensionality checks and mapping checks
