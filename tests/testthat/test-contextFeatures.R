test_that("Context-features: chromVAR-Activity scores precomputation functionality", {
  experiments(maeTest2)[[ACTEXP]] <- NULL
  experiments(maeTest2)[[ASSOCEXP]] <- NULL
  colData(maeTest2[[ATACEXP]])[[CHROMVAREXPCOL]] <- NULL
  metadata(maeTest2[[ATACEXP]])[[CHROMVARBGCOL]] <- NULL
  maeTest2 <- suppressSdWarning(panContextFeatures, list(mae=maeTest2))

  expect_s4_class(maeTest2, "MultiAssayExperiment")
  expect_contains(names(experiments(maeTest2)), ACTEXP)
  expect_contains(names(experiments(maeTest2)), ASSOCEXP)
  expect_contains(names(experiments(maeTest2)), ASSOCEXP)
  expect_contains(colnames(rowData(maeTest2[[ATACEXP]])),
                  c(CHROMVAREXPCOL))
  expect_contains(names(metadata(maeTest2[[ATACEXP]])),
                  CHROMVARBGCOL)
})

test_that("Context-features: pan-context feature computations", {
  experiments(maeTest2)[[ACTEXP]] <- NULL
  experiments(maeTest2)[[ASSOCEXP]] <- NULL
  colData(maeTest2[[ATACEXP]])[[CHROMVAREXPCOL]] <- NULL
  metadata(maeTest2[[ATACEXP]])[[CHROMVARBGCOL]] <- NULL
  maeTest2 <- suppressSdWarning(panContextFeatures, list(mae=maeTest2))

  expect_contains(colnames(rowData(maeTest2[[ATACEXP]])), ATACVARFEATNAME)
  expect_contains(colnames(colData(maeTest2[[ATACEXP]])), paste(MDSDIMFEATNAME, 1:2, sep="_"))
  expect_contains(colnames(rowData(maeTest2[[ATACEXP]])), MAXATACFEATNAME)
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
  experiments(maeTest2)[[ACTEXP]] <- NULL
  experiments(maeTest2)[[ASSOCEXP]] <- NULL
  colData(maeTest2[[ATACEXP]])[[CHROMVAREXPCOL]] <- NULL
  metadata(maeTest2[[ATACEXP]])[[CHROMVARBGCOL]] <- NULL

  outDir <-  tempdir()
  fileName  <- paste0(paste(ASSOCEXP, "mapped", sep="_"), ".h5")
  filePath <- file.path(outDir, fileName)
  maeTest2 <- suppressSdWarning(panContextFeatures, list(mae=maeTest2,
                                                         outDir=outDir,
                                                         saveHdf5=TRUE))
  expect_true(file.exists(filePath))
  file.remove(filePath)
})
