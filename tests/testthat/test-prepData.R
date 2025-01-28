test_that("Shifting checks", {
  #TODO: move this test to the corresponding file
  assayTableSimple1$strand <- c("+", "-", "+", "-")
  atacFragShifted <- .processData(assayTableSimple1, shift=TRUE)
  expect_equal(atacFragShifted$start, assayTableSimple1$start+c(4,0,4,0))
  expect_equal(atacFragShifted$end, assayTableSimple1$end+c(0,-5,0,-5))
})

test_that("Object construction: Basic functionality", {
  expect_no_error(suppressMessages({prepData(example_coords,
                           motifData=exampleMotif,
                           atacData=exampleATAC,
                           chIPData=exampleChIP)}))
})

test_that("Object construction: dimensionality check", {
  mae <- suppressMessages({prepData(example_coords,
                                    motifData=exampleMotif,
                                    atacData=exampleATAC,
                                    chIPData=exampleChIP)})

  expect_s4_class(mae, "MultiAssayExperiment")
  expect_equal(dim(experiments(mae)$Motifs), c(length(example_coords),
                                               length(exampleMotif)))
  expect_equal(dim(experiments(mae)$ChIP), c(length(example_coords),
                                             length(exampleChIP)))
  expect_equal(dim(experiments(mae)$ATAC), c(length(example_coords),
                                             length(exampleATAC)))
})

test_that("Object construction: Column naming", {
  mae <- suppressMessages({prepData(example_coords,
                                    motifData=exampleMotif,
                                    atacData=exampleATAC,
                                    chIPData=exampleChIP)})

  expect_equal(colnames(experiments(mae)$Motifs), names(exampleMotif))
  expect_equal(colnames(experiments(mae)$ChIP), names(exampleChIP))
  expect_equal(colnames(experiments(mae)$ATAC), names(exampleATAC))
})

test_that("Object construction: Train test assignment", {
  mae <- suppressMessages({prepData(example_coords,
                                    motifData=exampleMotif,
                                    atacData=exampleATAC,
                                    chIPData=exampleChIP,
                                    testSet="A549")})

  expect_equal(unique(subset(colData(mae),is_testing)$context), "A549")
  expect_equal(unique(subset(colData(mae),is_training)$context), "K562")
})

test_that("Object construction: Saving as hdf5", {
  outDir <-  tempdir()
  mae <- suppressWarnings({prepData(example_coords,
                                    motifData=exampleMotif,
                                    atacData=exampleATAC,
                                    chIPData=exampleChIP,
                                    saveHdf5=TRUE,
                                    outDir=outDir)})
  expect_true(file.exists(file.path(outDir, "ATAC_mapped")))
  expect_true(file.exists(file.path(outDir, "ChIP_mapped")))
  expect_true(file.exists(file.path(outDir, "Motifs_mapped")))

  expect_equal(dim(experiments(mae)$Motifs), c(length(example_coords),
                                               length(exampleMotif)))
  expect_equal(dim(experiments(mae)$ChIP), c(length(example_coords),
                                             length(exampleChIP)))
  expect_equal(dim(experiments(mae)$ATAC), c(length(example_coords),
                                             length(exampleATAC)))

  file.remove(file.path(outDir, "ATAC_mapped"))
  file.remove(file.path(outDir, "ChIP_mapped"))
  file.remove(file.path(outDir, "Motifs_mapped"))
})


test_that("Mappings check - addATACData", {
  # check that all mappings get retained
})

test_that("Coldata check - addATACData", {
  # check that is_testing, is_training is set correctly
})

test_that("Features check - addATACData", {
  #TODO: move this test to the corresponding file
  # check that features get computed for correct tf
})

test_that("Arguments check - addATACData", {
  # check if correct specification are provided
})


