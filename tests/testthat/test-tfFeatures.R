test_that("TF-features: Site & (pan)context-features are there", {
  experiments(maeTest)[[ACTEXP]] <- NULL
  experiments(maeTest)[[ASSOCEXP]] <- NULL
  expect_error(tfFeatures(maeTest, tfName="CTCF"),
               regexp="Context-specific")

  experiments(maeTest)[[SITEFEAT]] <- NULL
  expect_error(tfFeatures(maeTest, tfName="CTCF"),
               regexp="Site-specific features")
})

test_that("TF-features: Basic functionality", {
  experiments(maeTest)[[TFFEAT]] <- NULL
  maeTest <- tfFeatures(maeTest, tfName="CTCF",
                        features=c("Binding_Patterns",
                                   "Promoter_Association", "C_Score",
                                   "Cooccuring_Motifs"))

  expect_s4_class(maeTest, "MultiAssayExperiment")
  expect_contains(names(experiments(maeTest)), TFFEAT)
})

test_that("TF-features: Basic functionality - HDF5", {
  experiments(maeTestHdf5)[[TFFEAT]] <- NULL
  maeTestHdf5 <- tfFeatures(maeTestHdf5, tfName="CTCF",
                        features=c("Binding_Patterns",
                                   "Promoter_Association", "C_Score",
                                   "Cooccuring_Motifs"))

  expect_s4_class(maeTestHdf5, "MultiAssayExperiment")
  expect_contains(names(experiments(maeTestHdf5)), TFFEAT)
})

test_that("Assays are preserved when computing for new TF", {
  assayNamesOrig <- names(assays(maeTest[[TFFEAT]]))
  maeTest <- tfFeatures(maeTest, tfName="JUN",
                        features=c("Binding_Patterns",
                                   "Promoter_Association", "C_Score",
                                   "Cooccuring_Motifs"))
  assayNamesNew <- names(assays(maeTest[[TFFEAT]]))
  expect_equal(assayNamesNew, assayNamesOrig)
})

test_that("Preselected motifs are saved in colData", {
  experiments(maeTest2)[[TFFEAT]] <- NULL
  maeTest2 <- tfFeatures(maeTest2, tfName="JUN",
                         features=c("Binding_Patterns",
                                    "Promoter_Association", "C_Score",
                                    "Cooccuring_Motifs",
                                    "Associated_Motifs",
                                    "Associated_Motif_Activity"))
  expect_contains(colnames(colData(maeTest2[[TFFEAT]])),
                  c(PRESELMOTIFCOL, PRESELACTCOL))

  preSelMotifs <- colData(maeTest2[[TFFEAT]])[[PRESELMOTIFCOL]][[1]]
  expect_true(is.vector(preSelMotifs))
  expect_equal(preSelMotifs[[paste(TFMOTIFPREFIX, 1, sep="_")]][[1]], "JUN")

  preSelActMotifs <- colData(maeTest2[[TFFEAT]])[[PRESELACTCOL]][[1]]
  expect_true(is.vector(preSelActMotifs))
  expect_equal(preSelActMotifs[[paste(TFMOTIFPREFIX, 1, sep="_")]], "JUN")
})
