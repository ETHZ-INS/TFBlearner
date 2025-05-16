
test_that("Check correct ranges: Insertion counts", {
  insRes <- getInsertionProfiles(assayTableTest, motifCoords)
  insRanges <- unique(makeGRangesFromDataFrame(as.data.frame(insRes$motifScores)))
  isWithin <- length(subsetByOverlaps(insRanges,motifCoords, type="equal"))
  expect_identical(isWithin, length(insRanges))
})

test_that("Dimensionality checks: Insertion profile", {
  margin <- 100
  insRes <- getInsertionProfiles(assayTableTest, motifCoords, margin=margin)

  profPos <- unique(insRes$profile$rel_pos)
  profPos <- profPos[order(profPos)]

  expect_equal(profPos, seq(-margin, margin))
})

test_that("Correct insertion counting checks: With insertion profile calculation", {
  insRes <- getInsertionProfiles(assayTableSimple1, motifCoords, calcProfile=TRUE)
  insCounts <- insRes$motifScores
  setorder(insCounts, motif_match_id)

  expect_identical(insCounts$insert_counts, c(4L,4L))
})

test_that("Correct insertion counting checks: Without insertion profile calculation", {
  insRes <- getInsertionProfiles(assayTableSimple1, motifCoords, calcProfile=FALSE)
  insCounts <- insRes$motifScores
  setorder(insCounts, motif_match_id)

  expect_identical(insCounts$insert_counts, c(4L,4L))
})

test_that("Correct relative insertion counts: Insertion profile calculation", {
  margin <- 100
  insRes <- getInsertionProfiles(assayTableSimple1, motifCoords,
                                 margin=margin, calcProfile=TRUE)
  profile <- insRes$profile
  motifMarginRanges <- as.data.table(GenomicRanges::resize(motifCoords,
                                                           width=2*margin,
                                                           fix="center"))
  motifData <- as.data.table(motifCoords)
  motifData[,motif_center:=floor((end-start)/2)+start]
  motifData <- cbind(motifData, motifMarginRanges)

  relPos <- c(assayTableSimple1$start, assayTableSimple1$end)-motifData$motif_center
  relPos <- unique(relPos[which(abs(relPos)<=margin)])

  expect_equal(sum(profile$pos_count_global), 4L)
  expect_equal(subset(profile, rel_pos %in% relPos)$pos_count_global, c(2L,2L))
  expect_equal(insRes$motifScores$chi2, c(0L,0L))
})

test_that("Symmetric profile calculation check", {
  insRes <- getInsertionProfiles(assayTableTestLarge, motifCoords,
                                 calcProfile=TRUE, symmetric=TRUE,
                                 margin=10)

  wNeg <- subset(insRes$profile, rel_pos<0)
  setorder(wNeg, -rel_pos)
  wPos <- subset(insRes$profile, rel_pos>0)
  setorder(wPos, rel_pos)

  expect_identical(wNeg$w, wPos$w)
})

test_that("Simplified output format", {
  setnames(assayTableTestLarge, "group", "sample")
  insRes <- getInsertionProfiles(assayTableTestLarge, motifCoords,
                                 calcProfile=TRUE,
                                 simplified=TRUE)

  expect_s4_class(insRes, "RangedSummarizedExperiment")
  expect_identical(colnames(insRes),
                   unique(assayTableTestLarge$sample)[order(colnames(insRes))])
  expect_equal(names(assays(insRes)), c("insert_counts",
                                        "weighted_insert_counts",
                                        "chi2"))
  expect_identical(motifCoords, rowRanges(insRes))
})
