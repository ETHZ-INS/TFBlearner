test_that("Site features: Basic functionality", {
  experiments(maeTest)[[SITEFEAT]] <- NULL
  maeFeat <- siteFeatures(maeTest)

  expect_s4_class(maeFeat, "MultiAssayExperiment")
  expect_contains(names(experiments(maeFeat)), SITEFEAT)
})

test_that("Site features: Basic functionality - HDF5", {
  experiments(maeTestHdf5)[[SITEFEAT]] <- NULL
  maeFeat <- siteFeatures(maeTestHdf5)

  expect_s4_class(maeFeat, "MultiAssayExperiment")
  expect_contains(names(experiments(maeFeat)), SITEFEAT)
})

# add dimensionality checks and mapping checks
