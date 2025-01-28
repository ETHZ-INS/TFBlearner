library(data.table)
library(GenomicRanges)
library(BiocParallel)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(phastCons100way.UCSC.hg38)

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

# example data for mae for feature functions -----------------------------------

load(test_path("../../data/example_coords.rda"))
exampleMotif <- list(CTCF=test_path("../../inst/extdata/ctcf_motif.tsv"),
                     JUN=test_path("../../inst/extdata/jun_motif.tsv"))
exampleMotif <- list(CTCF="../../inst/extdata/ctcf_motif.tsv",
                     JUN="../../inst/extdata/jun_motif.tsv")
exampleATAC <-  list(A549="../../inst/extdata/example_atac_A549.bed",
                     K562="../../inst/extdata/example_atac_K562.bed")
exampleChIP <-  list(K562_CTCF="../../inst/extdata/example_chIP_K562_ctcf.tsv",
                     A549_CTCF="../../inst/extdata/example_chIP_A549_ctcf.tsv",
                     K562_JUN="../../inst/extdata/example_chIP_K562_jun.tsv")

# maeTest <- suppressMessages({prepData(example_coords,
#                                        motifData=exampleMotif,
#                                        atacData=exampleATAC,
#                                        chIPData=exampleChIP,
#                                        testSet="A549")})
# maeTest <- siteFeatures(maeTest)
# maeTest <- tfFeatures(maeTest, tfName="CTCF", tfCofactors="JUN")
# maeTest <- contextTfFeatures(maeTest, tfName="CTCF", subSample=20,
#                               features=c("Inserts", "Weighted_Inserts",
#                                          "Cofactor_Inserts"),
#                               addLabels=TRUE,
#                               BPPARAM=SerialParam())
# saveRDS(maeTest, "./test_data/maeTest.rds")

maeTest <- readRDS(test_path("./test_data/maeTest.rds"))

# maeTestHdf5 <- suppressMessages({prepData(example_coords,
#                                        motifData=exampleMotif,
#                                        atacData=exampleATAC,
#                                        chIPData=exampleChIP,
#                                        testSet="A549",
#                                        saveHdf5=TRUE,
#                                        outDir="./test_data")})
# maeTestHdf5 <- siteFeatures(maeTestHdf5)
# maeTestHdf5 <- tfFeatures(maeTestHdf5, tfName="CTCF", tfCofactors="JUN")
# maeTestHdf5 <- contextTfFeatures(maeTestHdf5, tfName="CTCF", subSample=20,
#                                  features=c("Inserts", "Weighted_Inserts",
#                                              "Cofactor_Inserts"),
#                                  addLabels=TRUE,
#                                  BPPARAM=SerialParam())
# saveRDS(maeTestHdf5, "./test_data/maeTest_hdf5.rds")

maeTestHdf5 <- readRDS(test_path("./test_data/maeTest_hdf5.rds"))

# Training & Prediction testing ------------------------------------------------

# fm <- getFeatureMatrix(maeTest, tfName="CTCF",
#                        addLabels=TRUE,
#                        saveHdf5=FALSE)
# saveRDS(fm, "./test_data/fmTest.rds")

fmTest <- readRDS(test_path("./test_data/fmTest.rds"))

