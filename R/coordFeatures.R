.getSeqFeatures <- function(refRanges,
                            phast=phastCons100way.UCSC.hg38,
                            BSgenome=BSgenome.Hsapiens.UCSC.hg38,
                            ...)
{
  # CpG content
  seqlevelsStyle(refRanges) <- "UCSC"
  cpgDens <- Repitools::cpgDensityCalc(refRanges, BSgenome)
  gcCont <- Repitools::gcContentCalc(refRanges, BSgenome)

  # PhastCon scores
  phastScores <- gscores(phast, refRanges, scores.only=TRUE)
  phastScores <- phastScores$default
  phastScores <- fifelse(is.na(phastScores), 0, phastScores)

  # sequence features
  # seqFeat <- data.table(phast_scores=phastScores$default,
  #                       cpg_density=cpgDens,
  #                       gc_cont=gcCont)
  # seqFeat <- as.matrix(seqFeat)
  #
  # colnames(seqFeat) <- c("conservation_scores", "cpg_density", "gc_content")
  # seqFeat <- as(seqFeat, "sparseMatrix")

  seqFeat <- list(phast_scores=Matrix::Matrix(phastScores, ncol=1),
                  cpg_density=Matrix::Matrix(cpgDens, ncol=1),
                  gc_cont=Matrix::Matrix(gcCont, ncol=1))

  names(seqFeat) <- c("conservation_scores", "cpg_density", "gc_content")

  return(seqFeat)
}

#' Adds an experiment with features specific for each site, such as CpG-content, conservation scores etc.
#' @name coordFeatures
#' @param mae MultiAssayExperiment as construced by [prepData()] containing ATAC- and ChIP-seq data.
#' @param aggFun Aggregation function for aggregating genomic scores overlapping the rowRanges of the provided MultiAssayExperiment object.
#' @param annoData Further data to be aggregated across the rowRanges of the provided MultiAssayExperiment object.
#'  List of paths pointing to .bed / .bam files or data.frames / data.tables or GRanges containing chr, start and end columns and score columns.
#' @param scoreCols Column names of the scores of annoData to be aggregated.
#' @param features Names of features to be added. Can be all or some of "Sequence", "Width", "Annot".
#' @import MultiAssayExperiment
#' @import Repitools
#' @import GenomicScores
#' @export
coordFeatures <- function(mae,
                          aggFun=max,
                          annoData=NULL,
                          scoreCols=NULL,
                          features=c("Sequence", "Width", "Annot"),
                          ...){

  features <- match.arg(features,
                        choices=c("Sequence", "Width", "Annot"),
                        several.ok=TRUE)
  featMats <- list()

  #  object validator
  .checkObject(mae)

  # Retrieve the coordinates
  coords <- rowRanges(experiments(mae)$Motifs)

  # seqFeatures: GC content etc
  if("Sequence" %in% features){
    seqFeat <- .getSeqFeatures(coords, ...)
    featMats <- append(featMats, seqFeat)}

  # reference coordinates width
  if("Width" %in% features){
    wi <- list(Matrix(width(coords), ncol=1))
    names(wi) <- "width"
    featMats <- append(featMats, wi)
  }

  if(!is.null(annoData) & "Annot" %in% features){

    if(length(scoreCols) < length(annoData)){
      message("Fewer colnames (scoreCols) to be aggregated provided, than
             annoData elements. Assuming all columns of elements of annoData
             are named the same")
      scoreCols <- rep(scoreCols, length(annoData))
    }

    annotFeats <- mapply(function(ad, an, ...){
      annoDt <- .processData(data, ...)
      annoDt <- ad
      annoDt$type <- an
      annoFeat <- genomicRangesMapping(coords,
                                        annoDt,
                                        byCols=c("type"),
                                        scoreCol=scoreCol,
                                        aggregationFun=aggFun,
                                        ...)},
      annoData, names(annoData), scoreCols, ...)
    #annotFeats <- Reduce(cbind, annotFeats[-1], annotFeats[[1]])
    names(annotFeats) <- names(annot)
    featMats <- append(featMats, annotFeats)
  }

  if(length(featMats)>0){
    #featMat <- Reduce("cbind", featMats[-1], featMats[[1]])
    mae <- .addFeatures(mae, featMats, mapTo="All", prefix="coordFeat", ...)}

  return(mae)
}
