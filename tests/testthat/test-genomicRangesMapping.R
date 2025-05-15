# the assayTable_test used for testing has been created like this:
# TODO: Put this in a setup-script !!!

# Arguments checks -------------------------------------------------------------

test_that("Arguments check: Error if no byCols provided", {
  expect_error(genomicRangesMapping(refCoords, assayTableTest))
})

test_that("Arguments check: Error if byCols not in assayTable are provided", {
  expect_error(genomicRangesMapping(refCoords, assayTableTest,
                                    byCols="condition"))
})

# Dimensionality checks --------------------------------------------------------

test_that("Dimensionality checks: Non-overlapping table", {
  nonOvAssayTable <- data.table(chr=c("chr10", "chr2"),
                                start=c(1,1), end=c(100,100),
                                group=c("TRT", "CTRL"))
  mapRes <- genomicRangesMapping(refCoords, nonOvAssayTable,
                                 byCols=c("group"))
  expect_identical(nrow(mapRes), length(refCoords))
})

test_that("Dimensionality checks: Number of ranges", {
  mapRes <- genomicRangesMapping(refCoords, assayTableTest,
                                 byCols=c("group"))
  expect_identical(nrow(mapRes), length(refCoords))
})

test_that("Dimensionality checks: Number of columns", {
  mapRes <- genomicRangesMapping(refCoords, assayTableTest,
                                 byCols=c("group"))
  expect_identical(ncol(mapRes), length(unique(assayTableTest$group)))
})

test_that("Dimensionality checks: Number of list elements", {
  mapRes <- genomicRangesMapping(refCoords, assayTableTest,
                                 byCols=c("tf","group"))
  expect_identical(length(mapRes), length(unique(assayTableTest$tf)))

  mapRes <- genomicRangesMapping(refCoords, assayTableTest,
                                 byCols=c("group","tf"))
  expect_identical(length(mapRes), length(unique(assayTableTest$group)))
})

# Correct aggregation ----------------------------------------------------------

test_that("Aggregation checks: Simple case mean aggregation - one byCol",{
  mapRes <- genomicRangesMapping(refCoords, assayTableSimple1,
                                 byCols=c("group"),
                                 scoreCol="score",
                                 aggregationFun=mean)
  expect_equal(mapRes@x, c(2.5,6.0))
  expect_equal(mapRes[1,1][[1]], 2.5)
  expect_equal(mapRes[10,2][[1]], 6.0)
})

test_that("Aggregation checks: Simple case mean aggregation - two byCols",{
  assayTableSimple <- rbind(assayTableSimple1, assayTableSimple2)
  mapRes <- genomicRangesMapping(refCoords, assayTableSimple,
                                 byCols=c("tf","group"),
                                 scoreCol="score",
                                 aggregationFun=mean)

  expect_equal(mapRes$YY1[1,1][[1]], 2.5)
  expect_equal(mapRes$YY1[10,2][[1]], 6.0)
  expect_equal(mapRes$GR[1,1][[1]], 5.0)
  expect_equal(mapRes$GR[10,2][[1]], 4.0)
})

test_that("Aggregation checks: Simple case sum aggregation - one byCol",{
  mapRes <- genomicRangesMapping(refCoords, assayTableSimple1,
                                 byCols=c("group"),
                                 scoreCol="score",
                                 aggregationFun=sum)

  expect_equal(mapRes[1,1][[1]], 5.0)
  expect_equal(mapRes[10,2][[1]], 12.0)
})

test_that("Aggregation checks: Simple case sum aggregation - two byCols",{
  assayTableSimple <- rbind(assayTableSimple1, assayTableSimple2)
  mapRes <- genomicRangesMapping(refCoords, assayTableSimple,
                                 byCols=c("tf","group"),
                                 scoreCol="score",
                                 aggregationFun=sum)

  expect_equal(mapRes$YY1[1,1][[1]], 5.0)
  expect_equal(mapRes$YY1[10,2][[1]], 12.0)
  expect_equal(mapRes$GR[1,1][[1]], 10.0)
  expect_equal(mapRes$GR[10,2][[1]], 8.0)
})

test_that("Aggregation checks: Simple case no aggregation function provided - one byCol",{
  mapRes <- genomicRangesMapping(refCoords, assayTableSimple1,
                                 byCols=c("group"),
                                 scoreCol="score")

  expect_equal(mapRes[1,1][[1]], 2.0)
  expect_equal(mapRes[10,2][[1]], 2.0)
})

test_that("Aggregation checks: Simple case no aggregation function provided - two byCols",{
  assayTableSimple <- rbind(assayTableSimple1, assayTableSimple2)
  mapRes <- genomicRangesMapping(refCoords, assayTableSimple,
                                 byCols=c("tf","group"),
                                 scoreCol="score")

  expect_equal(mapRes$YY1[1,1][[1]], 2.0)
  expect_equal(mapRes$YY1[10,2][[1]], 2.0)
  expect_equal(mapRes$GR[1,1][[1]], 2.0)
  expect_equal(mapRes$GR[10,2][[1]], 2.0)
})


test_that("Aggregation checks: Simple case no score column provided - one byCol",{
  mapRes <- genomicRangesMapping(refCoords, assayTableSimple1,
                                 byCols=c("group"),
                                 scoreCol="score")

  expect_equal(mapRes[1,1][[1]], 2.0)
  expect_equal(mapRes[10,2][[1]], 2.0)
})

test_that("Aggregation checks: Simple case no score column provided - two byCols",{
  assayTableSimple <- rbind(assayTableSimple1, assayTableSimple2)
  mapRes <- genomicRangesMapping(refCoords, assayTableSimple,
                                 byCols=c("tf","group"))

  expect_equal(mapRes$YY1[1,1][[1]], 2.0)
  expect_equal(mapRes$YY1[10,2][[1]], 2.0)
  expect_equal(mapRes$GR[1,1][[1]], 2.0)
  expect_equal(mapRes$GR[10,2][[1]], 2.0)
})

test_that("Aggregation checks: Reference implementation mean aggregation - one byCol",{

  mapRes <- genomicRangesMapping(refCoords, assayTableTest,
                                 byCols=c("group"),
                                 scoreCol="score",
                                 aggregationFun=mean)

  assayTableTestRanges <- makeGRangesFromDataFrame(as.data.frame(assayTableTest))
  ov <- findOverlaps(refCoords, assayTableTestRanges)
  subAssayTableTest <- assayTableTest[subjectHits(ov), ]
  subAssayTableTest$ref <- queryHits(ov)
  subAssayTableTest <- subAssayTableTest[,.(mean(score)), by=.(ref, group)]
  subAssayTableTest$group <- factor(subAssayTableTest$group,
                                    levels=c("CTRL", "TRT"), ordered=TRUE)
  refRes <- sparseMatrix(i=subAssayTableTest$ref,
                         j=subAssayTableTest$group,
                         x=subAssayTableTest$V1,
                         dims=c(length(refCoords),
                                length(unique(assayTableTest$group))))
  colnames(refRes) <- levels(subAssayTableTest$group)

  expect_equal(as.matrix(mapRes), as.matrix(refRes))
})

test_that("Aggregation checks: Reference implementation mean aggregation - two byCols",{

  mapRes <- genomicRangesMapping(refCoords, assayTableTest,
                                 byCols=c("tf","group"),
                                 scoreCol="score",
                                 aggregationFun=mean)

  tfs <- unique(assayTableTest$tf)
  refResMats <- lapply(tfs, function(tf_i){
    subAssayTableTest <- subset(assayTableTest, tf==tf_i)
    assayTableTestRanges <- makeGRangesFromDataFrame(as.data.frame(subAssayTableTest))
    ov <- findOverlaps(refCoords, assayTableTestRanges)
    subAssayTableTest <- subAssayTableTest[subjectHits(ov), ]
    subAssayTableTest$ref <- queryHits(ov)
    subAssayTableTest <- subAssayTableTest[,.(mean(score)), by=.(ref, group)]
    subAssayTableTest$group <- factor(subAssayTableTest$group,
                                      levels=c("CTRL", "TRT"), ordered=TRUE)
    refRes <- sparseMatrix(i=subAssayTableTest$ref,
                           j=subAssayTableTest$group,
                           x=subAssayTableTest$V1,
                           dims=c(length(refCoords),
                                  length(unique(assayTableTest$group))))
    colnames(refRes) <- c("CTRL", "TRT")
    refRes
  })
  names(refResMats) <- tfs

  expect_equal(as.matrix(mapRes$JUN), as.matrix(refResMats$JUN), ignore_attr=TRUE)
  expect_equal(as.matrix(mapRes$YY1), as.matrix(refResMats$YY1), ignore_attr=TRUE)
  expect_equal(as.matrix(mapRes$GR), as.matrix(refResMats$GR), ignore_attr=TRUE)
})

test_that("Aggregation checks: Reference implementation sum aggregation - one byCol",{
  mapRes <- genomicRangesMapping(refCoords, assayTableTest,
                                 byCols=c("group"),
                                 scoreCol="score",
                                 aggregationFun=sum)

  assayTableTestRanges <- makeGRangesFromDataFrame(as.data.frame(assayTableTest))
  ov <- findOverlaps(refCoords, assayTableTestRanges)
  subAssayTableTest <- assayTableTest[subjectHits(ov), ]
  subAssayTableTest$ref <- queryHits(ov)
  subAssayTableTest <- subAssayTableTest[,.(sum(score)), by=.(ref, group)]
  subAssayTableTest$group <- factor(subAssayTableTest$group,
                                    levels=c("CTRL", "TRT"), ordered=TRUE)
  refRes <- sparseMatrix(i=subAssayTableTest$ref,
                         j=subAssayTableTest$group,
                         x=subAssayTableTest$V1,
                         dims=c(length(refCoords),
                                length(unique(assayTableTest$group))))
  colnames(refRes) <- levels(subAssayTableTest$group)

  expect_equal(as.matrix(mapRes), as.matrix(refRes))
})

test_that("Aggregation checks: Reference implementation sum aggregation - two byCols",{

  mapRes <- genomicRangesMapping(refCoords, assayTableTest,
                                 byCols=c("tf","group"),
                                 scoreCol="score",
                                 aggregationFun=sum)

  tfs <- unique(assayTableTest$tf)
  refResMats <- lapply(tfs, function(tf_i){
    subAssayTableTest <- subset(assayTableTest, tf==tf_i)
    assayTableTestRanges <- makeGRangesFromDataFrame(as.data.frame(subAssayTableTest))
    ov <- findOverlaps(refCoords, assayTableTestRanges)
    subAssayTableTest <- subAssayTableTest[subjectHits(ov), ]
    subAssayTableTest$ref <- queryHits(ov)
    subAssayTableTest <- subAssayTableTest[,.(sum(score)), by=.(ref, group)]
    subAssayTableTest$group <- factor(subAssayTableTest$group,
                                      levels=c("CTRL", "TRT"), ordered=TRUE)
    refRes <- sparseMatrix(i=subAssayTableTest$ref,
                           j=subAssayTableTest$group,
                           x=subAssayTableTest$V1,
                           dims=c(length(refCoords),
                                  length(unique(assayTableTest$group))))
    colnames(refRes) <- c("CTRL", "TRT")
    refRes
  })
  names(refResMats) <- tfs

  expect_equal(as.matrix(mapRes$JUN), as.matrix(refResMats$JUN), ignore_attr=TRUE)
  expect_equal(as.matrix(mapRes$YY1), as.matrix(refResMats$YY1), ignore_attr=TRUE)
  expect_equal(as.matrix(mapRes$GR), as.matrix(refResMats$GR), ignore_attr=TRUE)
})

test_that("Aggregation checks: Reference implementation no aggregation function provided - one byCol",{
  mapRes <- genomicRangesMapping(refCoords, assayTableTest,
                                 byCols=c("group"),
                                 scoreCol="score")

  assayTableTestRanges <- makeGRangesFromDataFrame(as.data.frame(assayTableTest))
  ov <- findOverlaps(refCoords, assayTableTestRanges)
  subAssayTableTest <- assayTableTest[subjectHits(ov), ]
  subAssayTableTest$ref <- queryHits(ov)
  subAssayTableTest <- subAssayTableTest[,.(.N), by=.(ref, group)]
  subAssayTableTest$group <- factor(subAssayTableTest$group,
                                    levels=c("CTRL", "TRT"), ordered=TRUE)
  refRes <- sparseMatrix(i=subAssayTableTest$ref,
                         j=subAssayTableTest$group,
                         x=subAssayTableTest$N,
                         dims=c(length(refCoords),
                                length(unique(assayTableTest$group))))
  colnames(refRes) <- levels(subAssayTableTest$group)

  expect_equal(as.matrix(mapRes), as.matrix(refRes))
})

test_that("Aggregation checks: Reference implementation no aggregation function provided - two byCols",{

  mapRes <- genomicRangesMapping(refCoords, assayTableTest,
                                 byCols=c("tf","group"),
                                 scoreCol="score")

  tfs <- unique(assayTableTest$tf)
  refResMats <- lapply(tfs, function(tf_i){
    subAssayTableTest <- subset(assayTableTest, tf==tf_i)
    assayTableTestRanges <- makeGRangesFromDataFrame(as.data.frame(subAssayTableTest))
    ov <- findOverlaps(refCoords, assayTableTestRanges)
    subAssayTableTest <- subAssayTableTest[subjectHits(ov), ]
    subAssayTableTest$ref <- queryHits(ov)
    subAssayTableTest <- subAssayTableTest[,.(.N), by=.(ref, group)]
    subAssayTableTest$group <- factor(subAssayTableTest$group,
                                      levels=c("CTRL", "TRT"), ordered=TRUE)
    refRes <- sparseMatrix(i=subAssayTableTest$ref,
                           j=subAssayTableTest$group,
                           x=subAssayTableTest$N,
                           dims=c(length(refCoords),
                                  length(unique(assayTableTest$group))))
    colnames(refRes) <- c("CTRL", "TRT")
    refRes
  })
  names(refResMats) <- tfs

  expect_equal(as.matrix(mapRes$JUN), as.matrix(refResMats$JUN), ignore_attr=TRUE)
  expect_equal(as.matrix(mapRes$YY1), as.matrix(refResMats$YY1), ignore_attr=TRUE)
  expect_equal(as.matrix(mapRes$GR), as.matrix(refResMats$GR), ignore_attr=TRUE)
})


test_that("Aggregation checks: Reference implementation no aggregation function provided & chunked - one byCol",{
  mapRes <- genomicRangesMapping(refCoords, assayTableTest,
                                 byCols=c("group"),
                                 scoreCol="score",
                                 chunk=TRUE)

  assayTableTestRanges <- makeGRangesFromDataFrame(as.data.frame(assayTableTest))
  ov <- findOverlaps(refCoords, assayTableTestRanges)
  subAssayTableTest <- assayTableTest[subjectHits(ov), ]
  subAssayTableTest$ref <- queryHits(ov)
  subAssayTableTest <- subAssayTableTest[,.(.N), by=.(ref, group)]
  subAssayTableTest$group <- factor(subAssayTableTest$group,
                                    levels=c("CTRL", "TRT"), ordered=TRUE)
  refRes <- sparseMatrix(i=subAssayTableTest$ref,
                         j=subAssayTableTest$group,
                         x=subAssayTableTest$N,
                         dims=c(length(refCoords),
                                length(unique(assayTableTest$group))))
  colnames(refRes) <- levels(subAssayTableTest$group)

  expect_equal(as.matrix(mapRes[,colnames(refRes)]), as.matrix(refRes))
})

test_that("Aggregation checks: Reference implementation no aggregation function provided & chunked - two byCols",{

  mapRes <- genomicRangesMapping(refCoords, assayTableTest,
                                 byCols=c("tf","group"),
                                 scoreCol="score",
                                 chunk=TRUE)

  tfs <- unique(assayTableTest$tf)
  refResMats <- lapply(tfs, function(tf_i){
    subAssayTableTest <- subset(assayTableTest, tf==tf_i)
    assayTableTestRanges <- makeGRangesFromDataFrame(as.data.frame(subAssayTableTest))
    ov <- findOverlaps(refCoords, assayTableTestRanges)
    subAssayTableTest <- subAssayTableTest[subjectHits(ov), ]
    subAssayTableTest$ref <- queryHits(ov)
    subAssayTableTest <- subAssayTableTest[,.(.N), by=.(ref, group)]
    subAssayTableTest$group <- factor(subAssayTableTest$group,
                                      levels=c("CTRL", "TRT"), ordered=TRUE)
    refRes <- sparseMatrix(i=subAssayTableTest$ref,
                           j=subAssayTableTest$group,
                           x=subAssayTableTest$N,
                           dims=c(length(refCoords),
                                  length(unique(assayTableTest$group))))
    colnames(refRes) <- c("CTRL", "TRT")
    refRes
  })
  names(refResMats) <- tfs

  expect_equal(as.matrix(mapRes$JUN), as.matrix(refResMats$JUN), ignore_attr=TRUE)
  expect_equal(as.matrix(mapRes$YY1), as.matrix(refResMats$YY1), ignore_attr=TRUE)
  expect_equal(as.matrix(mapRes$GR), as.matrix(refResMats$GR), ignore_attr=TRUE)
})


test_that("Aggregation checks: Change of overlap type to within",{
  mapRes <- genomicRangesMapping(refCoords, assayTableTest,
                                 type="within",
                                 byCols=c("group"),
                                 scoreCol="score")

  assayTableTestRanges <- makeGRangesFromDataFrame(as.data.frame(assayTableTest))
  ov <- findOverlaps(refCoords, assayTableTestRanges, type="within")
  subAssayTableTest <- assayTableTest[subjectHits(ov), ]
  subAssayTableTest$ref <- queryHits(ov)
  subAssayTableTest <- subAssayTableTest[,.(.N), by=.(ref, group)]
  subAssayTableTest$group <- factor(subAssayTableTest$group,
                                    levels=c("CTRL", "TRT"), ordered=TRUE)
  refRes <- sparseMatrix(i=subAssayTableTest$ref,
                         j=subAssayTableTest$group,
                         x=subAssayTableTest$N,
                         dims=c(length(refCoords),
                                length(unique(assayTableTest$group))))
  colnames(refRes) <- levels(subAssayTableTest$group)

  expect_equal(as.matrix(mapRes), as.matrix(refRes))
})
