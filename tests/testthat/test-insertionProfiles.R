
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
