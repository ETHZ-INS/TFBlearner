test_that("Context-features: chromVAR-Activity scores precomputation functionality", {
  experiments(maeTest2)[[actExp]] <- NULL
  experiments(maeTest2)[[assocExp]] <- NULL
  colData(maeTest2[[atacExp]])[[chromVarExpCol]] <- NULL
  metadata(maeTest2[[atacExp]])[[chromVarBgCol]] <- NULL
  maeTest2 <- suppressSdWarning(panContextFeatures, list(mae=maeTest2))

  expect_s4_class(maeTest2, "MultiAssayExperiment")
  expect_contains(names(experiments(maeTest2)), actExp)
  expect_contains(names(experiments(maeTest2)), assocExp)
  expect_contains(names(experiments(maeTest2)), assocExp)
  expect_contains(colnames(rowData(maeTest2[[atacExp]])),
                  c(chromVarExpCol))
  expect_contains(names(metadata(maeTest2[[atacExp]])),
                  chromVarBgCol)
})

test_that("Context-features: pan-context feature computations", {
  experiments(maeTest2)[[actExp]] <- NULL
  experiments(maeTest2)[[assocExp]] <- NULL
  colData(maeTest2[[atacExp]])[[chromVarExpCol]] <- NULL
  metadata(maeTest2[[atacExp]])[[chromVarBgCol]] <- NULL
  maeTest2 <- suppressSdWarning(panContextFeatures, list(mae=maeTest2))

  expect_contains(colnames(rowData(maeTest2[[atacExp]])), atacVarFeatName)
  expect_contains(colnames(colData(maeTest2[[atacExp]])), paste(mdsDimFeatName, 1:2, sep="_"))
  expect_contains(colnames(rowData(maeTest2[[atacExp]])), maxAtacFeatName)
})

# add intest_that("Context-features: Message re-using precomputed chromVAR-Activity scores", {
test_that("Context-features: Message re-using precomputed chromVAR-Activity scores", {
  maeTest2 <- suppressSdWarning(panContextFeatures, list(mae=maeTest2))
  msgs <- list()
  withCallingHandlers(suppressSdWarning(panContextFeatures, list(mae=maeTest2)),
                      message = function(msg) {msgs <<- append(msgs,
                                                               msg$message)})
  msgs <- unlist(msgs)
  expect_true(sum(grepl("have been pre-computed and will reused", msgs))>0)
})

# add message in case features have been precomputed
test_that("Context-features: Save assocations as .h5-file", {
  experiments(maeTest2)[[actExp]] <- NULL
  experiments(maeTest2)[[assocExp]] <- NULL
  colData(maeTest2[[atacExp]])[[chromVarExpCol]] <- NULL
  metadata(maeTest2[[atacExp]])[[chromVarBgCol]] <- NULL

  outDir <-  tempdir()
  fileName  <- paste0(paste(assocExp, "mapped", sep="_"), ".h5")
  filePath <- file.path(outDir, fileName)
  maeTest2 <- suppressSdWarning(panContextFeatures, list(mae=maeTest2,
                                                         outDir=outDir,
                                                         saveHdf5=TRUE))
  expect_true(file.exists(filePath))
  file.remove(filePath)
})
