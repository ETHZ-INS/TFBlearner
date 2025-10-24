test_that("Context-TF-features: Basic functionality", {
  experiments(maeTest2)[[CONTEXTTFFEAT]] <- NULL
  maeTest2 <- contextTfFeatures(maeTest2, tfName="CTCF",
                                features=c("Inserts", "Weighted_Inserts",
                                           "ChromVAR_Scores"))

  expect_s4_class(maeTest2, "MultiAssayExperiment")
  expect_contains(names(experiments(maeTest2)), CONTEXTTFFEAT)
})

test_that("Context-TF-features: Basic functionality - HDF5", {
  experiments(maeTest)[[CONTEXTTFFEAT]] <- NULL
  maeTestHdf5 <- contextTfFeatures(maeTestHdf5, tfName="CTCF",
                        features=c("Inserts", "Weighted_Inserts"))

  expect_s4_class(maeTestHdf5, "MultiAssayExperiment")
  expect_contains(names(experiments(maeTestHdf5)), CONTEXTTFFEAT)
})

test_that("Context-TF-features: Correct training context selection", {
  experiments(maeTest)[[CONTEXTTFFEAT]] <- NULL
  maeTest <- contextTfFeatures(maeTest, tfName="CTCF",
                               whichCol="OnlyTrain",
                               features=c("Inserts", "Weighted_Inserts"))

  expect_s4_class(maeTest, "MultiAssayExperiment")
  expect_equal(rownames(colData(maeTest[[CONTEXTTFFEAT]])), "K562_CTCF")
})

test_that("Context-TF-features: Correct labelling",{
  maeTest <- tfFeatures(maeTest, tfName="JUN", tfCofactors="CTCF",
                        features="Binding_Patterns")
  maeTest <- contextTfFeatures(maeTest, tfName="JUN",
                               whichCol="OnlyTrain",
                               features=c("Inserts", "Weighted_Inserts"))
  expect_equal(assays(maeTest[[CHIPEXP]])[[PEAKASSAY]][,"K562_JUN",
                                                       drop=TRUE],
               assays(maeTest[[CONTEXTTFFEAT]])[[LABELCOLNAME]][,"K562_JUN",
                                                                drop=TRUE])
})

test_that("Assays are preserved when computing for new TF", {
  assayNamesOrig <- names(assays(maeTest[[CONTEXTTFFEAT]]))
  maeTest <- tfFeatures(maeTest, tfName="JUN",
                        features=c("CTCF", "MAX"))
  maeTest <- contextTfFeatures(maeTest, tfName="JUN",
                               features=c("Inserts", "Weighted_Inserts"),
                               addLabels=TRUE)
  assayNamesNew <- names(assays(maeTest[[CONTEXTTFFEAT]]))
  expect_equal(assayNamesNew, assayNamesOrig)
})

test_that("Error if features have not been computed for provided TF", {
  tfName="JUN"
  expect_error(contextTfFeatures(maeTest, tfName=tfName))
})

test_that("Using precomputed profile", {
  experiments(maeTest)[[CONTEXTTFFEAT]] <- NULL
  profile <- data.table(rel_pos=-200:200)
  profile[,w:=1/nrow(profile)]
  profile <- list("CTCF"=profile)

  expect_message(contextTfFeatures(maeTest, tfName="CTCF",
                                   insertionProfile=profile),
                 regexp="Using pre-computed insertion-profiles")
  expect_message(maeTest <- contextTfFeatures(maeTest, tfName="CTCF",
                                   insertionProfile=profile),
                 regexp="Skipped insertion-profiles computation. Using provided pre-computed ones")
  expect_equal(sum(is.na(assays(maeTest[["contextTfFeat"]])$contextTfFeat_weightedInserts.margin_tfMotif_1)), 0)
  expect_equal(sum(is.na(assays(maeTest[["contextTfFeat"]])$contextTfFeat_weightedInserts.within_tfMotif_1)), 0)

  expect_no_warning(contextTfFeatures(maeTest, tfName="CTCF",
                                      insertionProfile=profile))
  expect_no_message(contextTfFeatures(maeTest, tfName="CTCF",
                                      insertionProfile=profile),
                    message="Computing insertion-profiles")
})

test_that("Warning when using precomputed profile - not maching the motifRanges by name", {
  experiments(maeTest)[[CONTEXTTFFEAT]] <- NULL
  profile <- data.table(rel_pos=-200:200)
  profile[,w:=1/nrow(profile)]
  profile <- list("ATF2"=profile)

  expect_warning(contextTfFeatures(maeTest, tfName="CTCF",
                                   whichCol="OnlyTrain",
                                   insertionProfile=profile),
                 regexp="*Not all motif-ranges*")
  expect_message(suppressWarnings(contextTfFeatures(maeTest, tfName="CTCF",
                                  insertionProfile=profile)),
                 regexp="Computing insertion-profiles")
})
