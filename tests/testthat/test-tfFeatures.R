test_that("TF-features: Basic functionality", {
  experiments(maeTest)$tfFeat <- NULL
  maeTest <- tfFeatures(maeTest, tfName="CTCF",
                        features=c("Binding_Patterns",
                                   "Promoter_Association", "C_Score",
                                   "Cooccuring_Motifs"))

  expect_s4_class(maeTest, "MultiAssayExperiment")
  expect_contains(names(experiments(maeTest)), "tfFeat")
})

test_that("TF-features: Basic functionality - HDF5", {
  experiments(maeTestHdf5)$tfFeat <- NULL
  maeTestHdf5 <- tfFeatures(maeTestHdf5, tfName="CTCF",
                        features=c("Binding_Patterns",
                                   "Promoter_Association", "C_Score",
                                   "Cooccuring_Motifs"))

  expect_s4_class(maeTestHdf5, "MultiAssayExperiment")
  expect_contains(names(experiments(maeTestHdf5)), "tfFeat")
})

test_that("Assays are preserved when computing for new TF", {
  assayNamesOrig <- names(assays(maeTest[["tfFeat"]]))
  maeTest <- tfFeatures(maeTest, tfName="JUN",
                        features=c("Binding_Patterns",
                                   "Promoter_Association", "C_Score",
                                   "Cooccuring_Motifs"))
  assayNamesNew <- names(assays(maeTest[["tfFeat"]]))
  expect_equal(assayNamesNew, assayNamesOrig)
})


# add dimensionality checks and mapping checks
