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
  experiments(maeTest)[[contextTfFeat]] <- NULL
  maeTest <- contextTfFeatures(maeTest, tfName="CTCF",
                        features=c("Inserts", "Weighted_Inserts",
                                   "Cofactor_Inserts"))

  expect_s4_class(maeTest, "MultiAssayExperiment")
  expect_contains(names(experiments(maeTest)), contextTfFeat)
})

test_that("Context-TF-features: Basic functionality - HDF5", {
  experiments(maeTest)[[contextTfFeat]] <- NULL
  maeTestHdf5 <- contextTfFeatures(maeTestHdf5, tfName="CTCF",
                        features=c("Inserts", "Weighted_Inserts",
                                   "Cofactor_Inserts"))

  expect_s4_class(maeTestHdf5, "MultiAssayExperiment")
  expect_contains(names(experiments(maeTestHdf5)), contextTfFeat)
})

test_that("Context-TF-features: Correct training context selection", {
  experiments(maeTest)[[contextTfFeat]] <- NULL
  maeTest <- contextTfFeatures(maeTest, tfName="CTCF",
                               whichCol="OnlyTrain",
                               features=c("Inserts", "Weighted_Inserts",
                                          "Cofactor_Inserts"))

  expect_s4_class(maeTest, "MultiAssayExperiment")
  expect_equal(rownames(colData(maeTest[[contextTfFeat]])), "K562_CTCF")
})

test_that("Context-TF-features: save pre-computed ChromVAR parameters in colData of siteFeat", {
  experiments(maeTest)[[contextTfFeat]] <- NULL
  colData(experiments(maeTest)[[siteFeat]])[[chromVarExpName]] <- NULL
  colData(experiments(maeTest)[[siteFeat]])[[chromVarBgName]] <- NULL
  colData(experiments(maeTest)[[siteFeat]])[[topVarSitesName]] <- NULL
  maeTest <- suppressSdWarning(contextTfFeatures, list(mae=maeTest,
                                                       tfName="CTCF",
                                                       whichCol="OnlyTrain",
                                                       features=c("Inserts",
                                                            "ChromVAR_Scores")))

  expect_contains(colnames(colData(experiments(maeTest)[[siteFeat]])),
                  c(chromVarExpName,
                    chromVarBgName,
                    topVarSitesName))
})

test_that("Context-TF-features: message that pre-computed parameters/features are used", {
  experiments(maeTest2)[[contextTfFeat]] <- NULL
  colData(experiments(maeTest2)[[siteFeat]])[[chromVarExpName]] <- NULL
  colData(experiments(maeTest2)[[siteFeat]])[[topVarSitesName]] <- NULL
  colData(experiments(maeTest2)[[siteFeat]])[[chromVarBgName]] <- NULL

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
  experiments(maeTest)[[contextTfFeat]] <- NULL
  tfName <- "CTCF"

  tfMotifNames <- unlist(colData(maeTest[["tfFeat"]])$preselected_motifs)
  tfMotifNames <- names(tfMotifNames)[grepl("tf_motif", names(tfMotifNames))]
  maeTest <- suppressSdWarning(contextTfFeatures, list(mae=maeTest,
                                                       tfName=tfName,
                                                       whichCol="OnlyTrain",
                                                       features=c("Inserts",
                                                            "ChromVAR_Scores")))

  # expected column names
  expCols <- paste(chromVarAssocSuffix,
                   c(assocPearsonPrefix, assocCohenPrefix), tfName, sep="_")
  expect_contains(colnames(rowData(experiments(maeTest)[[tfFeat]])), expCols)
})


test_that("Context-TF-features: save pre-computed MDS-dimensions in colData of ATAC", {
  experiments(maeTest2)[[contextTfFeat]] <- NULL
  tfName <- "CTCF"
  maeTest2 <- contextTfFeatures(maeTest2, tfName=tfName,
                                whichCol="All",
                                addLabels=FALSE,
                                features=c("Inserts",
                                           "MDS_Context"))

  # expected column names
  expect_contains(colnames(colData(experiments(maeTest2)[[atacExp]])),
                  paste(mdsDimFeatName, 1:2, sep="_"))
})

test_that("Context-TF-features: save pre-computed maximal ATAC-signal in rowData of ATAC", {
  experiments(maeTest)[[contextTfFeat]] <- NULL
  tfName <- "CTCF"
  maeTest <- contextTfFeatures(maeTest, tfName=tfName,
                               whichCol="All",
                               addLabels=FALSE,
                               features=c("Inserts",
                                          "Max_ATAC_Signal"))

  # expected column names
  expect_contains(colnames(rowData(experiments(maeTest)[[atacExp]])),
                  maxAtacFeatName)
})

test_that("Context-TF-features: save pre-computed ATAC-signal variance in rowData of ATAC", {
  experiments(maeTest)[[contextTfFeat]] <- NULL
  tfName <- "CTCF"
  maeTest <- contextTfFeatures(maeTest, tfName=tfName,
                               whichCol="All",
                               addLabels=FALSE,
                               features=c("Inserts",
                                          "ATAC_Variance"))

  # expected column names
  expect_contains(colnames(rowData(experiments(maeTest)[[atacExp]])),
                  atacVarFeatName)
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

# add dimensionality checks and mapping checks
