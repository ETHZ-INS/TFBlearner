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

# add dimensionality checks and mapping checks
