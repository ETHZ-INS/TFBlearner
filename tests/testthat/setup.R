library(data.table)
library(GenomicRanges)
library(BiocParallel)

refCoords <- readRDS(test_path("./test_data/refCoords_test.rds"))

# Setup for genomicRangesMapping tests------------------------------------------

n <- runif(1,90,100)
tfs <- c("YY1", "GR", "JUN")
groups <- c("CTRL", "TRT")

assayTableTest <- data.table(chr=paste0("chr", sample(1:20, n, replace=TRUE)),
                             group=c(groups, sample(groups, n-2, replace=TRUE)),
                             tf=c(tfs, sample(tfs, n-3, replace=TRUE)),
                             score=runif(n))
assayTableTest[, start:=round(sample(start(refCoords),n,replace=TRUE)+runif(n,-100,100))]
assayTableTest[, end:=round(start+runif(n, 1,100))]

# contains overlap to only the first and last range
assayTableSimple1 <- data.table(chr=c("chr1", "chr1", "chr1", "chr1"),
                                start=c(10070,10070,79415,79415),
                                group=c("CTRL", "CTRL", "TRT", "TRT"),
                                tf=rep("YY1", 4),
                                score=c(1,4,5,7))
assayTableSimple1[,end:=start+10]
assayTableSimple2 <- data.table(chr=c("chr1", "chr1", "chr1", "chr1"),
                                start=c(10070,10070,79415, 79415),
                                group=c("CTRL", "CTRL", "TRT", "TRT"),
                                tf=rep("GR", 4),
                                score=c(2,8,2,6))
assayTableSimple2[,end:=start+10]

#countOverlaps(GRanges(seqnames="chr1", ranges=IRanges(start=2, end=100)),
#              GRanges(seqnames="chr1", ranges=IRanges(start=50, end=70)))
#countOverlaps(GRanges(seqnames="chr1", ranges=IRanges(start=2, end=100)),
#              GRanges(seqnames="chr1", ranges=IRanges(start=1, end=110)))

# Setup for insertionProfiles tests---------------------------------------------
nLarge <- 1e4

motifCoords <- data.table(chr=c("chr1", "chr1"),
                          start=c(10069, 79414),
                          motif_id="YY1")
motifCoords[,end:=start+100]
motifCoords <- makeGRangesFromDataFrame(as.data.frame(motifCoords))

assayTableTestLarge <- data.table(chr=paste0("chr", sample(1:20, nLarge, replace=TRUE)),
                             group=c(groups, sample(groups, nLarge-2, replace=TRUE)),
                             tf=c(tfs, sample(tfs, nLarge-3, replace=TRUE)),
                             score=runif(nLarge))
assayTableTestLarge[, start:=round(sample(start(motifCoords),nLarge,replace=TRUE)+runif(nLarge,-100,100))]
assayTableTestLarge[, end:=round(start+runif(nLarge, 1,100))]

