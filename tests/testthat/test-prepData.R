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

test_that("Object construction: Motif scores in column data", {
  mae <- suppressMessages({prepData(example_coords,
                                    motifData=exampleMotif,
                                    atacData=exampleATAC,
                                    chIPData=exampleChIP)})
  expect_contains(colnames(colData(experiments(mae)$Motifs)), "max_score")

  motMat <- assays(experiments(mae)$Motifs)$match_scores
  exp <- lapply(colnames(motMat), function(col) max(motMat[,col,drop=TRUE]))
  exp <- unlist(exp)
  names(exp) <- names(colnames(motMat))

  cd <-  colData(experiments(mae)$Motifs)
  cd <- cd[order(match(cd$motif,names(exp))),,drop=FALSE]
  obs <- cd$max_score

  expect_equivalent(obs, exp)
})

test_that("Object construction: Saving as hdf5", {
  outDir <-  tempdir()
  mae <- suppressMessages({prepData(example_coords,
                                    motifData=exampleMotif,
                                    atacData=exampleATAC,
                                    chIPData=exampleChIP,
                                    saveHdf5=TRUE,
                                    outDir=outDir)})
  expect_true(file.exists(file.path(outDir, "ATAC_mapped.h5")))
  expect_true(file.exists(file.path(outDir, "ChIP_mapped.h5")))
  expect_true(file.exists(file.path(outDir, "Motif_mapped.h5")))

  expect_equal(dim(experiments(mae)$Motifs), c(length(example_coords),
                                               length(exampleMotif)))
  expect_equal(dim(experiments(mae)$ChIP), c(length(example_coords),
                                             length(exampleChIP)))
  expect_equal(dim(experiments(mae)$ATAC), c(length(example_coords),
                                             length(exampleATAC)))

  file.remove(file.path(outDir, "ATAC_mapped.h5"))
  file.remove(file.path(outDir, "ChIP_mapped.h5"))
  file.remove(file.path(outDir, "Motif_mapped.h5"))
})


test_that("Mappings check - addATACData", {
  atacData <- exampleATAC
  names(atacData) <- c("MCF7", "JURKAT")

  maeAdd <- addATACData(maeTest, atacData, shift=FALSE)
  expect_s4_class(maeTest, "MultiAssayExperiment")

  expect_equal(c(colnames(experiments(maeTest)$ATAC),  c("MCF7", "JURKAT")),
               colnames(experiments(maeAdd)$ATAC))
  motifs <- subset(sampleMap(maeAdd), primary=="MCF7" & assay=="Motifs")$colname
  motifs <- motifs[order(motifs)]
  motifsAll <- subset(sampleMap(maeTest), assay=="Motifs")$colname
  motifsAll <- motifsAll[order(motifsAll)]
  expect_equal(unique(motifs), unique(motifsAll))
})

test_that("Mappings check HDF5 - addATACData", {
  atacData <- exampleATAC
  names(atacData) <- c("MCF7", "JURKAT")

  maeAddHdf5 <- addATACData(maeTestHdf5, atacData, shift=FALSE)
  expect_s4_class(maeAddHdf5, "MultiAssayExperiment")

  motifs <- subset(sampleMap(maeAddHdf5), primary=="MCF7" & assay=="Motifs")$colname
  motifs <- motifs[order(motifs)]
  motifsAll <- subset(sampleMap(maeAddHdf5), assay=="Motifs")$colname
  motifsAll <- motifsAll[order(motifsAll)]
  expect_equal(unique(motifs), unique(motifsAll))
})

test_that("Coldata check - addATACData", {
  atacData <- exampleATAC
  names(atacData) <- c("MCF7", "JURKAT")

  maeAdd <- addATACData(maeTest, atacData, shift=FALSE, testSet="MCF7")
  testCons <- subset(colData(maeAdd), is_testing)$context
  trainCons <- subset(colData(maeAdd), is_training)$context

  expect_contains(testCons, "MCF7")
  expect_true(!("MCF7" %in% trainCons))
})

test_that("Features check - addATACData", {
  atacData <- exampleATAC
  names(atacData) <- c("MCF7", "JURKAT")

  maeAdd <- addATACData(maeTest, atacData, shift=FALSE,
                        tfName="CTCF",
                        computeFeatures=TRUE,
                        features=c("Inserts"))

  expect_contains(colnames(experiments(maeAdd)$contextTfFeat),
                  c("MCF7_CTCF", "JURKAT_CTCF"))
})

test_that("Features check  HDF5- addATACData", {
  atacData <- exampleATAC
  names(atacData) <- c("MCF7", "JURKAT")

  maeAdd <- addATACData(maeTestHdf5, atacData, shift=FALSE,
                        tfName="CTCF",
                        computeFeatures=TRUE,
                        features=c("Inserts"))

  expect_contains(colnames(experiments(maeAdd)$contextTfFeat),
                  c("MCF7_CTCF", "JURKAT_CTCF"))
})
