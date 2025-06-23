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
  expect_equal(dim(mae[[MOTIFEXP]]), c(length(example_coords),
                                       length(exampleMotif)))
  expect_equal(dim(mae[[CHIPEXP]]), c(length(example_coords),
                                      length(exampleChIP)))
  expect_equal(dim(mae[[ATACEXP]]), c(length(example_coords),
                                      length(exampleATAC)))
})

test_that("Object construction: Column naming", {
  mae <- suppressMessages({prepData(example_coords,
                                    motifData=exampleMotif,
                                    atacData=exampleATAC,
                                    chIPData=exampleChIP)})

  expect_equal(colnames(mae[[MOTIFEXP]]), names(exampleMotif))
  expect_equal(colnames(mae[[CHIPEXP]]), names(exampleChIP))
  expect_equal(colnames(mae[[ATACEXP]]), names(exampleATAC))
})

test_that("Object construction: Train test assignment", {
  mae <- suppressMessages({prepData(example_coords,
                                    motifData=exampleMotif,
                                    atacData=exampleATAC,
                                    chIPData=exampleChIP,
                                    testSet="A549")})

  baseDir <- metadata(colData(mae[[MOTIFEXP]]))[[BASEDIRCOL]]
  mmPath <- subset(colData(mae[[MOTIFEXP]]),
                   get(MOTIFNAMECOL)==names(exampleMotif)[[1]])$origin
  expect_true(file.exists(file.path(baseDir, mmPath)))

  mmPath <- subset(colData(mae[[MOTIFEXP]]),
                   get(MOTIFNAMECOL)==names(exampleMotif)[[2]])$origin
  expect_true(file.exists(file.path(baseDir, mmPath)))
})

test_that("Object construction: Train test assignment", {
  mae <- suppressMessages({prepData(example_coords,
                                    motifData=exampleMotif,
                                    atacData=exampleATAC,
                                    chIPData=exampleChIP,
                                    testSet="A549")})

  expect_equal(unique(subset(colData(mae), get(ISTESTCOL))$context), "A549")
  expect_equal(unique(subset(colData(mae), get(ISTRAINCOL))$context), "K562")
})

test_that("Object construction: ATAC-frag path assignment", {
  mae <- suppressMessages({prepData(example_coords,
                                    motifData=exampleMotif,
                                    atacData=exampleATAC,
                                    chIPData=exampleChIP,
                                    testSet="A549")})
  baseDirs <- lapply(exampleATAC, dirname)
  baseDir <- unique(unlist(baseDirs))

  obsBaseDir <- metadata(colData(mae[[ATACEXP]]))[[BASEDIRCOL]]
  expect_equal(obsBaseDir, baseDir)

  expFileNames <- unlist(lapply(colData(mae[[ATACEXP]])$origin, names))
  expect_equal(colData(mae[[ATACEXP]])[["context"]], expFileNames)
  expect_true(all(file.exists(file.path(obsBaseDir,
                                        colData(mae[[ATACEXP]])$origin))))
})

test_that("Object construction: Motif scores in column data", {
  mae <- suppressMessages({prepData(example_coords,
                                    motifData=exampleMotif,
                                    atacData=exampleATAC,
                                    chIPData=exampleChIP)})
  expect_contains(colnames(colData(mae[[MOTIFEXP]])), MAXSCORECOL)

  motMat <- assays(mae[[MOTIFEXP]])[[MATCHASSAY]]
  exp <- lapply(colnames(motMat), function(col) max(motMat[,col,drop=TRUE]))
  exp <- unlist(exp)
  names(exp) <- names(colnames(motMat))

  cd <-  colData(mae[[MOTIFEXP]])
  cd <- cd[order(match(cd[[MOTIFNAMECOL]],names(exp))),,drop=FALSE]
  obs <- cd[[MAXSCORECOL]]

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

  expect_equal(dim(mae[[MOTIFEXP]]), c(length(example_coords),
                                       length(exampleMotif)))
  expect_equal(dim(mae[[CHIPEXP]]), c(length(example_coords),
                                      length(exampleChIP)))
  expect_equal(dim(mae[[ATACEXP]]), c(length(example_coords),
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

  expect_equal(c(colnames(maeTest[[ATACEXP]]),  c("MCF7", "JURKAT")),
                 colnames(maeAdd[[ATACEXP]]))
  motifs <- subset(sampleMap(maeAdd), primary=="MCF7" & assay==MOTIFEXP)$colname
  motifs <- motifs[order(motifs)]
  motifsAll <- subset(sampleMap(maeTest), assay==MOTIFEXP)$colname
  motifsAll <- motifsAll[order(motifsAll)]
  expect_equal(unique(motifs), unique(motifsAll))
})

test_that("Mappings check HDF5 - addATACData", {
  atacData <- exampleATAC
  names(atacData) <- c("MCF7", "JURKAT")

  maeAddHdf5 <- addATACData(maeTestHdf5, atacData, shift=FALSE)
  expect_s4_class(maeAddHdf5, "MultiAssayExperiment")

  motifs <- subset(sampleMap(maeAddHdf5), primary=="MCF7" & assay==MOTIFEXP)$colname
  motifs <- motifs[order(motifs)]
  motifsAll <- subset(sampleMap(maeAddHdf5), assay==MOTIFEXP)$colname
  motifsAll <- motifsAll[order(motifsAll)]
  expect_equal(unique(motifs), unique(motifsAll))
})

test_that("Coldata check - addATACData", {
  atacData <- exampleATAC
  names(atacData) <- c("MCF7", "JURKAT")

  maeAdd <- addATACData(maeTest, atacData, shift=FALSE, testSet="MCF7")
  testCons <- subset(colData(maeAdd), get(ISTESTCOL))$context
  trainCons <- subset(colData(maeAdd), get(ISTRAINCOL))$context

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
  expect_true(all(c(ATACVARFEATNAME, MAXATACFEATNAME)
                  %in% colnames(rowData(maeAdd[[ATACEXP]]))))
  expect_contains(colnames(maeAdd[[CONTEXTTFFEAT]]),
                  c("MCF7_CTCF", "JURKAT_CTCF"))
})

test_that("Features simple check ATAC-frag paths - addATACData",{

  atacData <- exampleATAC
  names(atacData) <- c("MCF7", "JURKAT")

  maeAdd <- suppressSdWarning(addATACData,
                              list(mae=maeTest, atacData=atacData,
                                   shift=FALSE,
                                   computeFeatures="simple"))

  partialPaths <-  colData(maeTest[[ATACEXP]])$origin
  fullPaths <- file.path(metadata(colData(maeTest[[ATACEXP]]))[[BASEDIRCOL]],
                         partialPaths)
  names(fullPaths) <- unlist(lapply(partialPaths, names))

  expBaseDir <- .commonDir(c(atacData, fullPaths))
  obsBaseDir <- metadata(colData(maeAdd[[ATACEXP]]))[[BASEDIRCOL]]
  expect_equal(obsBaseDir, expBaseDir)

  expFileNames <- unlist(lapply(colData(maeAdd[[ATACEXP]])$origin, names))
  expect_equal(colData(maeAdd[[ATACEXP]])[["context"]], expFileNames)
  expect_true(all(file.exists(file.path(obsBaseDir,
                                        colData(maeAdd[[ATACEXP]])$origin))))
})

test_that("Features simple check including MDS-projection - addATACData", {
  atacData <- exampleATAC
  names(atacData) <- c("MCF7", "JURKAT")

  maeAdd <- suppressSdWarning(addATACData,
                              list(mae=maeTest2, atacData=atacData,
                                   shift=FALSE,
                                   computeFeatures="simple"))

  # check that col / row data preserved
  expect_true(all(c(ATACVARFEATNAME, MAXATACFEATNAME)
                  %in% colnames(rowData(maeAdd[[ATACEXP]]))))
  expect_contains(colnames(maeAdd[[CONTEXTTFFEAT]]),
                  c("MCF7_CTCF", "JURKAT_CTCF"))
  mdsDims <- paste(MDSDIMFEATNAME, 1:2,sep="_")
  expect_contains(colnames(colData(maeAdd[[ATACEXP]])), mdsDims)
  expect_equal(nrow(colData(maeAdd[[ATACEXP]])), 5)
  expect_equal(sum(is.na(colData(maeAdd[[ATACEXP]])[,mdsDims])), 0)
})

test_that("Features from scratch check - addATACData", {
  atacData <- exampleATAC
  names(atacData) <- c("MCF7", "JURKAT")

  maeAdd <- suppressSdWarning(addATACData,
                              list(mae=maeTest, atacData=atacData,
                                   shift=FALSE,
                                   computeFeatures="scratch"))

  # check that col / row data preserved
  expect_true(all(c(ATACVARFEATNAME, MAXATACFEATNAME)
                  %in% colnames(rowData(maeAdd[[ATACEXP]]))))
  expect_contains(colnames(maeAdd[[CONTEXTTFFEAT]]),
                  c("MCF7_CTCF", "JURKAT_CTCF"))
})

test_that("Features from scratch check ATAC-frag paths - addATACData",{

  atacData <- exampleATAC
  names(atacData) <- c("MCF7", "JURKAT")

  maeAdd <- suppressSdWarning(addATACData,
                              list(mae=maeTest, atacData=atacData,
                                   shift=FALSE,
                                   computeFeatures="scratch"))

  partialPaths <-  colData(maeTest[[ATACEXP]])$origin
  fullPaths <- file.path(metadata(colData(maeTest[[ATACEXP]]))[[BASEDIRCOL]],
                         partialPaths)
  names(fullPaths) <- unlist(lapply(partialPaths, names))

  expBaseDir <- .commonDir(c(atacData, fullPaths))
  obsBaseDir <- metadata(colData(maeAdd[[ATACEXP]]))[[BASEDIRCOL]]
  expect_equal(obsBaseDir, expBaseDir)

  expFileNames <- unlist(lapply(colData(maeAdd[[ATACEXP]])$origin, names))
  expect_equal(colData(maeAdd[[ATACEXP]])[["context"]], expFileNames)
  expect_true(all(file.exists(file.path(obsBaseDir,
                                        colData(maeAdd[[ATACEXP]])$origin))))
})

test_that("Features none check - addATACData", {
  atacData <- exampleATAC
  names(atacData) <- c("MCF7", "JURKAT")

  maeAdd <- suppressSdWarning(addATACData,
                              list(mae=maeTest, atacData=atacData,
                                   shift=FALSE,
                                   computeFeatures="none"))

  expect_true(all(c(ATACVARFEATNAME, MAXATACFEATNAME)
                  %in% colnames(rowData(maeAdd[[ATACEXP]]))))
})

test_that("Features check  HDF5- addATACData", {
  atacData <- exampleATAC
  names(atacData) <- c("MCF7", "JURKAT")

  maeAdd <- addATACData(maeTestHdf5, atacData, shift=FALSE,
                        computeFeatures="simple")

  expect_contains(colnames(maeAdd[[CONTEXTTFFEAT]]),
                  c("MCF7_CTCF", "JURKAT_CTCF"))
})
