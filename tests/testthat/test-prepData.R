test_that("Shifting checks", {
  #TODO: move this test to the corresponding file
  assayTableSimple1$strand <- c("+", "-", "+", "-")
  atacFragShifted <- .processData(assayTableSimple1, shift=TRUE)
  expect_equal(atacFragShifted$start, assayTableSimple1$start+c(4L,4L,4L,4L))
  expect_equal(atacFragShifted$end, assayTableSimple1$end+c(-4L,-4L,-4L,-4L))
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
  expect_equal(dim(mae[[motifExp]]), c(length(example_coords),
                                       length(exampleMotif)))
  expect_equal(dim(mae[[chIPExp]]), c(length(example_coords),
                                      length(exampleChIP)))
  expect_equal(dim(mae[[atacExp]]), c(length(example_coords),
                                      length(exampleATAC)))
})

test_that("Object construction: Column naming", {
  mae <- suppressMessages({prepData(example_coords,
                                    motifData=exampleMotif,
                                    atacData=exampleATAC,
                                    chIPData=exampleChIP)})

  expect_equal(colnames(mae[[motifExp]]), names(exampleMotif))
  expect_equal(colnames(mae[[chIPExp]]), names(exampleChIP))
  expect_equal(colnames(mae[[atacExp]]), names(exampleATAC))
})

test_that("Object construction: Train test assignment", {
  mae <- suppressMessages({prepData(example_coords,
                                    motifData=exampleMotif,
                                    atacData=exampleATAC,
                                    chIPData=exampleChIP,
                                    testSet="A549")})

  expect_true(file.exists(subset(colData(mae[[motifExp]]),
                          get(motifNameCol)==names(exampleMotif)[[1]])$origin))
  expect_true(file.exists(subset(colData(mae[[motifExp]]),
                          get(motifNameCol)==names(exampleMotif)[[2]])$origin))
})

test_that("Object construction: Train test assignment", {
  mae <- suppressMessages({prepData(example_coords,
                                    motifData=exampleMotif,
                                    atacData=exampleATAC,
                                    chIPData=exampleChIP,
                                    testSet="A549")})

  expect_equal(unique(subset(colData(mae), get(isTestCol))$context), "A549")
  expect_equal(unique(subset(colData(mae), get(isTrainCol))$context), "K562")
})

test_that("Object construction: Motif scores in column data", {
  mae <- suppressMessages({prepData(example_coords,
                                    motifData=exampleMotif,
                                    atacData=exampleATAC,
                                    chIPData=exampleChIP)})
  expect_contains(colnames(colData(mae[[motifExp]])), maxScoreCol)

  motMat <- assays(mae[[motifExp]])[[matchAssay]]
  exp <- lapply(colnames(motMat), function(col) max(motMat[,col,drop=TRUE]))
  exp <- unlist(exp)
  names(exp) <- names(colnames(motMat))

  cd <-  colData(mae[[motifExp]])
  cd <- cd[order(match(cd[[motifNameCol]],names(exp))),,drop=FALSE]
  obs <- cd[[maxScoreCol]]

  expect_equal(obs, exp, ignore_attr = TRUE)
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

  expect_equal(dim(mae[[motifExp]]), c(length(example_coords),
                                       length(exampleMotif)))
  expect_equal(dim(mae[[chIPExp]]), c(length(example_coords),
                                      length(exampleChIP)))
  expect_equal(dim(mae[[atacExp]]), c(length(example_coords),
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

  expect_equal(c(colnames(maeTest[[atacExp]]),  c("MCF7", "JURKAT")),
                 colnames(maeAdd[[atacExp]]))
  motifs <- subset(sampleMap(maeAdd), primary=="MCF7" & assay==motifExp)$colname
  motifs <- motifs[order(motifs)]
  motifsAll <- subset(sampleMap(maeTest), assay==motifExp)$colname
  motifsAll <- motifsAll[order(motifsAll)]
  expect_equal(unique(motifs), unique(motifsAll))
})

test_that("Mappings check HDF5 - addATACData", {
  atacData <- exampleATAC
  names(atacData) <- c("MCF7", "JURKAT")

  maeAddHdf5 <- addATACData(maeTestHdf5, atacData, shift=FALSE)
  expect_s4_class(maeAddHdf5, "MultiAssayExperiment")

  motifs <- subset(sampleMap(maeAddHdf5), primary=="MCF7" & assay==motifExp)$colname
  motifs <- motifs[order(motifs)]
  motifsAll <- subset(sampleMap(maeAddHdf5), assay==motifExp)$colname
  motifsAll <- motifsAll[order(motifsAll)]
  expect_equal(unique(motifs), unique(motifsAll))
})

test_that("Coldata check - addATACData", {
  atacData <- exampleATAC
  names(atacData) <- c("MCF7", "JURKAT")

  maeAdd <- addATACData(maeTest, atacData, shift=FALSE, testSet="MCF7")
  testCons <- subset(colData(maeAdd), get(isTestCol))$context
  trainCons <- subset(colData(maeAdd), get(isTrainCol))$context

  expect_contains(testCons, "MCF7")
  expect_true(!("MCF7" %in% trainCons))
})

test_that("Features simple check - addATACData", {
  atacData <- exampleATAC
  names(atacData) <- c("MCF7", "JURKAT")

  maeAdd <- suppressSdWarning(addATACData,
                              list(mae=maeTest, atacData=atacData,
                                   shift=FALSE,
                                   computeFeatures="simple"))

  # check that col / row data preserved
  expect_true(all(c(atacVarFeatName, maxAtacFeatName)
                  %in% colnames(rowData(maeAdd[[atacExp]]))))
  expect_contains(colnames(maeAdd[[contextTfFeat]]),
                  c("MCF7_CTCF", "JURKAT_CTCF"))
})

test_that("Features simple check including MDS-projection - addATACData", {
  atacData <- exampleATAC
  names(atacData) <- c("MCF7", "JURKAT")

  maeAdd <- suppressSdWarning(addATACData,
                              list(mae=maeTest2, atacData=atacData,
                                   shift=FALSE,
                                   computeFeatures="simple"))

  # check that col / row data preserved
  expect_true(all(c(atacVarFeatName, maxAtacFeatName)
                  %in% colnames(rowData(maeAdd[[atacExp]]))))
  expect_contains(colnames(maeAdd[[contextTfFeat]]),
                  c("MCF7_CTCF", "JURKAT_CTCF"))
  mdsDims <- paste(mdsDimFeatName, 1:2,sep="_")
  expect_contains(colnames(colData(maeAdd[[atacExp]])), mdsDims)
  expect_equal(nrow(colData(maeAdd[[atacExp]])), 5)
  expect_equal(sum(is.na(colData(maeAdd[[atacExp]])[,mdsDims])), 0)
})

test_that("Features from scratch check - addATACData", {
  atacData <- exampleATAC
  names(atacData) <- c("MCF7", "JURKAT")

  maeAdd <- suppressSdWarning(addATACData,
                              list(mae=maeTest, atacData=atacData,
                                   shift=FALSE,
                                   computeFeatures="scratch"))

  # check that col / row data preserved
  expect_true(all(c(atacVarFeatName, maxAtacFeatName)
                  %in% colnames(rowData(maeAdd[[atacExp]]))))
  expect_contains(colnames(maeAdd[[contextTfFeat]]),
                  c("MCF7_CTCF", "JURKAT_CTCF"))
})

test_that("Features none check - addATACData", {
  atacData <- exampleATAC
  names(atacData) <- c("MCF7", "JURKAT")

  maeAdd <- suppressSdWarning(addATACData,
                              list(mae=maeTest, atacData=atacData,
                                   shift=FALSE,
                                   computeFeatures="none"))

  expect_true(all(c(atacVarFeatName, maxAtacFeatName)
                  %in% colnames(rowData(maeAdd[[atacExp]]))))
})

test_that("Features check  HDF5- addATACData", {
  atacData <- exampleATAC
  names(atacData) <- c("MCF7", "JURKAT")

  maeAdd <- addATACData(maeTestHdf5, atacData, shift=FALSE,
                        computeFeatures="simple")

  expect_contains(colnames(maeAdd[[contextTfFeat]]),
                  c("MCF7_CTCF", "JURKAT_CTCF"))
})
