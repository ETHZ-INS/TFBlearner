library(data.table)
library(GenomicRanges)
library(BiocParallel)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(phastCons100way.UCSC.hg38)

refCoords <- readRDS(test_path("./test_data/refCoords_test.rds"))

# to muffle the warning "In cor(...): the standard deviation is zero"
# when computing ChromVAR-ATAC association.
# This warning is supposed to appear in the given test setup but
# is not informative for the test cases.
suppressSdWarning <- function(fun, args){

  msg <- "the standard deviation is zero"
  withCallingHandlers(
    res <- do.call(fun, args),
    warning=function(w){
      if(grepl(msg, conditionMessage(w))){
        invokeRestart("muffleWarning")}
    })
  return(res)
}

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

load(file.path(system.file("data", package="TFBlearner"), "example_coords.rda"))
load(file.path(system.file("data", package="TFBlearner"), "example_pwms.rda"))

exampleATAC <- list(A549=system.file("extdata", "example_atac_A549.bed", package = "TFBlearner"),
                    K562=system.file("extdata", "example_atac_K562.bed", package = "TFBlearner"))
exampleChIP <- list(K562_CTCF=system.file("extdata", "example_chIP_K562_ctcf.tsv", package = "TFBlearner"),
                    A549_CTCF=system.file("extdata", "example_chIP_A549_ctcf.tsv", package = "TFBlearner"),
                    K562_JUN=system.file("extdata", "example_chIP_K562_jun.tsv", package = "TFBlearner"))

exampleMotif <- prepMotifs(example_coords, example_pwms,
                           genome = BSgenome.Hsapiens.UCSC.hg38,
                           outputFolder="./test_data")

maeTest <- suppressMessages({prepData(example_coords,
                                       motifData=exampleMotif,
                                       atacData=exampleATAC,
                                       chIPData=exampleChIP,
                                       testSet="A549")})
maeTest <- siteFeatures(maeTest)
maeTest <- tfFeatures(maeTest, tfName="CTCF", tfCofactors="JUN")
maeTest <- contextTfFeatures(maeTest, tfName="CTCF", subSample=20,
                              features=c("Inserts", "Weighted_Inserts"),
                              addLabels=TRUE,
                              BPPARAM=SerialParam())

exampleATAC2 <- list(A549=system.file("extdata", "example_atac_A549.bed", package = "TFBlearner"),
                     K562=system.file("extdata", "example_atac_K562.bed", package = "TFBlearner"),
                     HepG2=system.file("extdata", "example_atac_K562.bed", package = "TFBlearner")) # dummy type for some tests

maeTest2 <- suppressMessages({prepData(example_coords,
                                      motifData=exampleMotif,
                                      atacData=exampleATAC2,
                                      chIPData=exampleChIP,
                                      testSet="A549")})
maeTest2 <- siteFeatures(maeTest2)
maeTest2 <- tfFeatures(maeTest2, tfName="CTCF", tfCofactors="JUN")

maeTestHdf5 <- suppressMessages({prepData(example_coords,
                                       motifData=exampleMotif,
                                       atacData=exampleATAC,
                                       chIPData=exampleChIP,
                                       testSet="A549",
                                       saveHdf5=TRUE,
                                       outDir="./test_data")})
maeTestHdf5 <- siteFeatures(maeTestHdf5)
maeTestHdf5 <- tfFeatures(maeTestHdf5, tfName="CTCF", tfCofactors="JUN")
maeTestHdf5 <- contextTfFeatures(maeTestHdf5, tfName="CTCF", subSample=20,
                                 features=c("Inserts", "Weighted_Inserts",
                                             "Cofactor_Inserts"),
                                 addLabels=TRUE,
                                 BPPARAM=SerialParam())

# Training & Prediction testing ------------------------------------------------

fmSe <- getFeatureMatrix(maeTest, tfName="CTCF",
                         addLabels=TRUE,
                         saveHdf5=FALSE)
fmTest <- assays(fmSe)$features
