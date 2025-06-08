test_that("TF-features: Basic functionality", {
  experiments(maeTest)[[tfFeat]] <- NULL
  maeTest <- tfFeatures(maeTest, tfName="CTCF",
                        features=c("Binding_Patterns",
                                   "Promoter_Association", "C_Score",
                                   "Cooccuring_Motifs"))

  expect_s4_class(maeTest, "MultiAssayExperiment")
  expect_contains(names(experiments(maeTest)), tfFeat)
})

test_that("TF-features: Basic functionality - HDF5", {
  experiments(maeTestHdf5)[[tfFeat]] <- NULL
  maeTestHdf5 <- tfFeatures(maeTestHdf5, tfName="CTCF",
                        features=c("Binding_Patterns",
                                   "Promoter_Association", "C_Score",
                                   "Cooccuring_Motifs"))

  expect_s4_class(maeTestHdf5, "MultiAssayExperiment")
  expect_contains(names(experiments(maeTestHdf5)), tfFeat)
})

test_that("Assays are preserved when computing for new TF", {
  assayNamesOrig <- names(assays(maeTest[[tfFeat]]))
  maeTest <- tfFeatures(maeTest, tfName="JUN",
                        features=c("Binding_Patterns",
                                   "Promoter_Association", "C_Score",
                                   "Cooccuring_Motifs"))
  assayNamesNew <- names(assays(maeTest[[tfFeat]]))
  expect_equal(assayNamesNew, assayNamesOrig)
})

test_that("Preselected motifs are saved in colData", {
  experiments(maeTest2)[[tfFeat]] <- NULL
  maeTest2 <- tfFeatures(maeTest2, tfName="JUN",
                         features=c("Binding_Patterns",
                                    "Promoter_Association", "C_Score",
                                    "Cooccuring_Motifs",
                                    "Associated_Motifs",
                                    "Associated_Motif_Activity"))
  expect_contains(colnames(colData(maeTest2[[tfFeat]])),
                  c(preSelMotifCol, preSelActCol))

  preSelMotifs <- colData(maeTest2[[tfFeat]])[[preSelMotifCol]][[1]]
  expect_true(is.vector(preSelMotifs))
  expect_equal(preSelMotifs[[paste(tfMotifPrefix, 1, sep="_")]][[1]], "JUN")

  preSelActMotifs <- colData(maeTest2[[tfFeat]])[[preSelActCol]][[1]]
  expect_true(is.vector(preSelActMotifs))
  expect_equal(preSelActMotifs[[paste(tfMotifPrefix, 1, sep="_")]], "JUN")
})
