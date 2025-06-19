
test_that("Check correct ranges: Insertion counts", {
  insRes <- getInsertionProfiles(assayTableTest, motifCoords)
  insRanges <- unique(makeGRangesFromDataFrame(as.data.frame(insRes[[RETSCORESNAME]])))
  isWithin <- length(subsetByOverlaps(insRanges,motifCoords, type="equal"))
  expect_identical(isWithin, length(insRanges))
})

test_that("Dimensionality checks: Insertion profile", {
  margin <- 100
  insRes <- getInsertionProfiles(assayTableTest, motifCoords, margin=margin)

  profPos <- unique(insRes[[REPROFILENAME]]$rel_pos)
  profPos <- profPos[order(profPos)]

  expect_equal(profPos, seq(-margin, margin))
})

test_that("Correct insertion counting checks: With insertion profile calculation", {
  insRes <- getInsertionProfiles(assayTableSimple1, motifCoords, calcProfile=TRUE)
  insCounts <- insRes[[RETSCORESNAME]]
  setorder(insCounts, motif_match_id)

  expect_identical(insCounts[[INSERTFEATNAME]], c(4L,4L))
})

test_that("Correct insertion counting checks: Without insertion profile calculation", {
  insRes <- getInsertionProfiles(assayTableSimple1, motifCoords, calcProfile=FALSE)
  insCounts <- insRes[[RETSCORESNAME]]
  setorder(insCounts, motif_match_id)

  expect_identical(insCounts[[INSERTFEATNAME]], c(4L,4L))
})

test_that("Correct relative insertion counts: Insertion profile calculation", {
  margin <- 100
  insRes <- getInsertionProfiles(assayTableSimple1, motifCoords,
                                 margin=margin, calcProfile=TRUE)
  profile <- insRes[[REPROFILENAME]]
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
  expect_equal(insRes[[RETSCORESNAME]][[DEVFEATNAME]], c(0L,0L))
})

test_that("Symmetric profile calculation check", {
  insRes <- getInsertionProfiles(assayTableTestLarge, motifCoords,
                                 calcProfile=TRUE, symmetric=TRUE,
                                 margin=10)

  wNeg <- subset(insRes[[REPROFILENAME]], rel_pos<0)
  setorder(wNeg, -rel_pos)
  wPos <- subset(insRes[[REPROFILENAME]], rel_pos>0)
  setorder(wPos, rel_pos)

  expect_identical(wNeg$w, wPos$w)
})

test_that("Using precomputed profiles", {

  profile <- data.table(rel_pos=-200:200,w=1)
  profile[,w:=w/sum(w)]
  profiles <- list("JUN"=profile, "YY1"=profile, "USF1"=profile)

  motifCoords1 <- motifCoords
  motifCoords1$motif_id <- "YY1"
  motifCoords2 <- motifCoords
  motifCoords2$motif_id <- "JUN"
  motifCoords <- c(motifCoords1, motifCoords2)

  insRes <- NULL
  expect_message(insRes <- getInsertionProfiles(assayTableTestLarge, motifCoords,
                                                calcProfile=FALSE,
                                                symmetric=TRUE,
                                                profiles=profiles,
                                                margin=10),
                 regexp="Using precomputed profiles")
  profiles <- rbindlist(profiles, idcol="motif_id")
  expect_equal(insRes$insertProfiles, profiles)
})

test_that("No matching inserts", {
  margin <- 10
  atacFrag <- data.table(chr="chr1", start=100, end=200, sample="sample1")
  insRes <- NULL
  expect_no_error(insRes <- getInsertionProfiles(atacFrag,
                                                 motifCoords,
                                                 calcProfile=TRUE,
                                                 symmetric=TRUE,
                                                 margin=margin))
  nPos <- 2*margin+1
  expect_equal(insRes[[REPROFILENAME]]$w, rep(1/nPos, nPos))

  atacFrag <- data.table(chr="chr5", start=100, end=200, sample="sample1")
  expect_error(getInsertionProfiles(atacFrag, motifCoords,
                                    calcProfile=TRUE, symmetric=TRUE,
                                    margin=margin),
               regexp="No common chromosomes")
})

test_that("Simplified output format", {
  setnames(assayTableTestLarge, "group", "sample")
  insRes <- getInsertionProfiles(assayTableTestLarge, motifCoords,
                                 calcProfile=TRUE,
                                 simplified=TRUE)

  expect_s4_class(insRes, "RangedSummarizedExperiment")
  expect_identical(colnames(insRes),
                   unique(assayTableTestLarge$sample)[order(colnames(insRes))])
  expect_equal(names(assays(insRes)), c(INSERTFEATNAME,
                                        WINSERTSFEATNAME,
                                        DEVFEATNAME))
  expect_identical(motifCoords, rowRanges(insRes))
})
