test_that("Feature Matrix: Basic functionality", {
  fm <- getFeatureMatrix(maeTest, tfName="CTCF",
                         addLabels=TRUE,
                         saveHdf5=FALSE)
  contexts <- getContexts(maeTest, tfName="CTCF", which="Both")

  expect_s4_class(assays(fm)$features, "CsparseMatrix")
  expect_equal(nrow(fm), length(contexts)*length(example_coords))
})

test_that("Feature Matrix: Basic functionality without labels", {
  fm <- getFeatureMatrix(maeTest, tfName="CTCF",
                         addLabels=FALSE,
                         saveHdf5=FALSE)
  contexts <- getContexts(maeTest, tfName="CTCF", which="ATAC")

  expect_s4_class(assays(fm)$features, "CsparseMatrix")
  expect_equal(nrow(fm), length(contexts)*length(example_coords))
})

test_that("Feature Matrix: Basic functionality - HDF5", {
  fm <- getFeatureMatrix(maeTestHdf5, tfName="CTCF",
                         addLabels=FALSE,
                         saveHdf5=FALSE)
  contexts <- getContexts(maeTest, tfName="CTCF", which="ATAC")

  expect_s4_class(assays(fm)$features, "CsparseMatrix")
  expect_equal(nrow(fm), length(contexts)*length(example_coords))
})

test_that("Feature Matrix: Basic functionality - Saving to HDF5", {
  outDir <- tempdir()
  prefix <- "test"
  fm <- getFeatureMatrix(maeTest, tfName="CTCF",
                                  addLabels=FALSE,
                                  saveHdf5=TRUE,
                                  outDir=outDir,
                                  prefix=prefix)

  fmFilePath <- file.path(outDir,
                          paste(prefix, "feature_matrix", "CTCF", sep="_"))
  fmFilePath <- paste0(fmFilePath, ".h5")
  expect_s4_class(assays(fm)$features, "DelayedMatrix")
  expect_true(file.exists(fmFilePath))
  file.remove(fmFilePath)
})

test_that("Feature Matrix: Correct context selection - only for training contexts", {
  annoCol <- "context"
  fm <- getFeatureMatrix(maeTest, tfName="CTCF",
                         whichCol="OnlyTrain",
                         addLabels=FALSE,
                         annoCol=annoCol,
                         saveHdf5=FALSE)

  testContext <- subset(colData(maeTest), get(ISTESTCOL))$context
  contexts <- getContexts(maeTest, tfName="CTCF", which="ATAC")
  trainContexts <- setdiff(contexts, testContext)

  obsRanges <-  rowRanges(fm)
  mcols(obsRanges) <- NULL
  names(obsRanges) <- NULL

  expRanges <- rep(rowRanges(maeTest[[ATACEXP]]), length(trainContexts))
  mcols(expRanges) <- NULL

  expect_equal(as.character(levels(rowRanges(fm)@elementMetadata[[annoCol]])),
               trainContexts)
  expect_equal(nrow(fm), length(trainContexts)*length(example_coords))
  expect_equal(obsRanges, expRanges, ignore_attr=TRUE)
  expect_equal(metadata(fm)[[annoCol]], trainContexts)
})

test_that("Feature Matrix: Correct context selection - only for specified context", {
  annoCol <- "context"
  tfName <- "CTCF"
  fm <- getFeatureMatrix(maeTest, tfName=tfName,
                         whichCol="Col",
                         colSel="A549",
                         addLabels=FALSE,
                         annoCol=annoCol,
                         saveHdf5=FALSE)

  obsRanges <- rowRanges(fm)
  mcols(obsRanges) <- NULL
  names(obsRanges) <- NULL
  expRanges <- rowRanges(maeTest[[ATACEXP]])
  mcols(expRanges) <- NULL
  names(expRanges) <- NULL

  expect_equal(nrow(fm), length(example_coords))
  expect_equal(as.character(levels(rowRanges(fm)@elementMetadata[[annoCol]])),
               "A549")
  expect_equal(obsRanges, expRanges)
  expect_equal(metadata(fm)[[annoCol]], "A549")
  expect_equal(metadata(fm)[[TFCOFACTORSCOL]], "JUN")
})

test_that("Feature Matrix: Correct metadata assignment", {

  annoCol <- "context"
  tfName <- "CTCF"
  fm <- getFeatureMatrix(maeTest, tfName=tfName,
                         whichCol="Col",
                         colSel="A549",
                         addLabels=FALSE,
                         annoCol=annoCol,
                         saveHdf5=FALSE)

  expect_contains(names(metadata(fm)), PRESELMOTIFCOL)
  preSelMotifs <- metadata(fm)[[PRESELMOTIFCOL]]
  expect_equal(preSelMotifs[[paste(TFMOTIFPREFIX, 1, sep="_")]], tfName)

  expect_contains(names(metadata(fm)), PRESELACTCOL)
  preSelActMotifs <- metadata(fm)[[PRESELACTCOL]]
  expect_equal(preSelActMotifs[[paste(TFMOTIFPREFIX, 1, sep="_")]], tfName)
})

test_that("Feature Matrix: Column names corresponding to R conventions", {
  fm <- getFeatureMatrix(maeTest, tfName="CTCF",
                         whichCol="Col",
                         colSel="A549",
                         addLabels=FALSE,
                         saveHdf5=FALSE)

  # should pass anyways given current setup script, for later
  expect_equal(colnames(fm), make.names(colnames(fm), unique=TRUE))
})

test_that("Feature Matrix: Correct context selection - no context provided", {
  expect_error(getFeatureMatrix(maeTest, tfName="CTCF",
                                whichCol="Col",
                                addLabels=FALSE,
                                saveHdf5=FALSE))
})

test_that("Feature Matrix: warning if features have not been computed for provided cellular contexts", {
  experiments(maeTest)[[CONTEXTTFFEAT]] <- NULL
  maeTest <- contextTfFeatures(maeTest, tfName="CTCF",
                               whichCol="Col",
                               colSel="K562",
                               features=c("Inserts"))

  expect_warning(getFeatureMatrix(maeTest, tfName="CTCF",
                                  whichCol="Col",
                                  colSel=c("K562", "A549"),
                                  addLabels=FALSE,
                                  saveHdf5=FALSE),
                 "Not all the cellular-contexts requested have contextTfFeats compiled.\n
                   Missing are: A549")
})

test_that("Feature Matrix: error if features have not been computed for provided TF", {
  tfName="JUN"
  expect_error(getFeatureMatrix(maeTest, tfName=tfName,
                                  whichCol="Col",
                                  colSel=c("K562", "A549"),
                                  addLabels=FALSE,
                                  saveHdf5=FALSE))
})

test_that("Feature Matrix: columns are named correctly", {
  annoCol <- "context"
  fm <- getFeatureMatrix(maeTest, tfName="CTCF",
                         addLabels=TRUE,
                         annoCol=annoCol,
                         saveHdf5=FALSE)

  lDt <- listFeatures()
  enumFeats <- unlist(tstrsplit(colnames(fm), split="_"))
  enumFeats <- as.integer(enumFeats[grepl("^[-+]?[0-9]+$", enumFeats)])
  maxEnum <- max(enumFeats)
  colNamesAll <- lapply(1:maxEnum, gsub, pattern="<i>",
                        unlist(lDt$feature_matrix_column_names))
  colNamesAll <- unlist(colNamesAll)
  colNamesAll <- c(colNamesAll,
                   paste(colNamesAll, NORMEDAFFIX, sep="_"),
                   paste(colNamesAll, NORMEDMAXAFFIX, sep="_"),
                   paste(colNamesAll, GCNORMEDAFFIX, sep="_"))

  expect_setequal(setdiff(colnames(fm), colNamesAll),
                  c(annoCol, LABELCOLNAME))
})


test_that("Feature Matrix: correct features are normalized", {
  annoCol <- "context"
  fm <- getFeatureMatrix(maeTest, tfName="CTCF",
                         addLabels=TRUE,
                         annoCol=annoCol,
                         saveHdf5=FALSE)

  lDt <- listFeatures()
  lDt <- subset(lDt, !context_normed)
  enumFeats <- unlist(tstrsplit(colnames(fm), split="_"))
  enumFeats <- as.integer(enumFeats[grepl("^[-+]?[0-9]+$", enumFeats)])
  maxEnum <- max(enumFeats)
  colNamesAll <- lapply(1:maxEnum, gsub, pattern="<i>",
                        unlist(lDt$feature_matrix_column_names))
  colNamesAll <- unlist(colNamesAll)
  colNamesAll <- c(colNamesAll,
                   paste(colNamesAll, NORMEDAFFIX, sep="_"),
                   paste(colNamesAll, NORMEDMAXAFFIX, sep="_"),
                   paste(colNamesAll, GCNORMEDAFFIX, sep="_"))

  normedCols <- colnames(fm)[grepl(NORMEDAFFIX,colnames(fm))]
  expect_length(setdiff(normedCols, colNamesAll), 0)
})
