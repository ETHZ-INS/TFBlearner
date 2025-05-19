test_that("Site features: Basic functionality", {
  experiments(maeTest)[[siteFeat]] <- NULL
  maeFeat <- siteFeatures(maeTest)

  expect_s4_class(maeFeat, "MultiAssayExperiment")
  expect_contains(names(experiments(maeFeat)), siteFeat)
})

test_that("Site features: Basic functionality - HDF5", {
  experiments(maeTestHdf5)[[siteFeat]] <- NULL
  maeFeat <- siteFeatures(maeTestHdf5)

  expect_s4_class(maeFeat, "MultiAssayExperiment")
  expect_contains(names(experiments(maeFeat)), siteFeat)
})

# add dimensionality checks and mapping checks
