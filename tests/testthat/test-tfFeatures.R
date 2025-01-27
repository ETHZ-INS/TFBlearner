test_that("Site features: Basic functionality", {
  maeFeat <- siteFeatures(maeTest)
  maeFeat <- tfFeatures(maeFeat, tfName="CTCF",
                        features=c("Binding_Patterns",
                                   "Promoter_Association", "C_Score",
                                   "Cooccuring_Motifs"))

  expect_s4_class(maeFeat, "MultiAssayExperiment")
  expect_contains(names(experiments(maeFeat)), "tfFeat")
})

test_that("Site features: Basic functionality", {
  maeFeat <- siteFeatures(maeTestHdf5)
  maeFeat <- tfFeatures(maeFeat, tfName="CTCF",
                        features=c("Binding_Patterns",
                                   "Promoter_Association", "C_Score",
                                   "Cooccuring_Motifs"))

  expect_s4_class(maeFeat, "MultiAssayExperiment")
  expect_contains(names(experiments(maeFeat)), "tfFeat")
})

# add dimensionality checks and mapping checks
