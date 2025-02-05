.getSeqFeatures <- function(refRanges,
                            phast=phastCons100way.UCSC.hg38,
                            genome=BSgenome.Hsapiens.UCSC.hg38)
{
  # CpG content
  seqlevelsStyle(refRanges) <- "UCSC"
  cpgDens <- Repitools::cpgDensityCalc(refRanges, genome)
  gcCont <- Repitools::gcContentCalc(refRanges, genome)

  # PhastCon scores
  phastScores <- gscores(phast, refRanges, scores.only=TRUE)
  phastScores <- phastScores$default
  phastScores <- fifelse(is.na(phastScores), 0, phastScores)

  seqFeat <- list(phast_scores=Matrix::Matrix(phastScores, ncol=1),
                  cpg_density=Matrix::Matrix(cpgDens, ncol=1),
                  gc_cont=Matrix::Matrix(gcCont, ncol=1))

  names(seqFeat) <- c("conservation_scores", "cpg_density", "gc_content")

  return(seqFeat)
}

#' Site-specific features
#'
#' Adds an experiment with features specific for each site, such as CpG-content and conservation scores, to the provided MultiAssayExperiment object.
#' @name siteFeatures
#' @param mae [MultiAssayExperiment::MultiAssayExperiment-class] as construced by [TFBlearner::prepData()] containing Motif, ATAC- and ChIP-seq data.
#' @param aggregationFun Aggregation function for aggregating genomic scores overlapping the rowRanges of the provided MultiAssayExperiment object.
#' @param annoData Further data to be aggregated across the rowRanges of the provided MultiAssayExperiment object.
#'  Named list of paths pointing to .bed / .bam files or data.frames / [data.table::data.table] or [GenomicRanges::GRanges-class] containing chr, start and end columns and score columns.
#' @param scoreCols Column names of the scores of annoData to be aggregated. Needs to be same length as annoData.
#' If of length one, it is assumed that score columns of all list elements have the same name.
#' @param features Names of features to be added. Can be all or some of "Sequence", "Width", "Annot".
#' Features are stored in the assays of the added experiment. See [TFBlearner::listFeatures] for an overview of the features.
#' @param phast phastCons conservation scores to be used.
#' @param genome [BSgenome::BSgenome-class] to be used.
#' @param ... Arguments to be passed to [TFBlearner::genomicRangesMapping] for aggregation of annoData.
#' @return [MultiAssayExperiment::MultiAssayExperiment-class] object with an experiment containing site-specific features added to [MultiAssayExperiment::experiments].
#' @import MultiAssayExperiment
#' @importFrom Repitools cpgDensityCalc gcContentCalc
#' @importFrom GenomicScores gscores
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @export
siteFeatures <- function(mae,
                         aggregationFun=max,
                         annoData=NULL,
                         scoreCols=NULL,
                         features=c("Sequence", "Width", "Annot"),
                         phast=phastCons100way.UCSC.hg38,
                         genome=BSgenome.Hsapiens.UCSC.hg38,
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
    seqFeat <- .getSeqFeatures(coords, genome=genome, phast=phast)
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

    annotFeats <- mapply(function(ad, an, scoreCol, ...){
      annoDt <- .processData(ad, readAll=TRUE)
      annoDt$type <- an
      annoFeat <- genomicRangesMapping(coords,
                                       annoDt,
                                       byCols=c("type"),
                                       scoreCol=scoreCol,
                                       aggregationFun=aggregationFun,
                                       ...)},
      annoData, names(annoData), scoreCols, ...)

    names(annotFeats) <- names(annoData)
    featMats <- append(featMats, annotFeats)
  }

  names(featMats) <- paste("siteFeat", names(featMats), sep="_")
  seSiteFeat <- SummarizedExperiment(assays=featMats, rowRanges=coords)
  colnames(seSiteFeat) <- "all"
  colData(seSiteFeat)$feature_type <- "siteFeature"
  colsToMap <- unique(sampleMap(mae)$primary)
  mae <- .addFeatures(mae, seSiteFeat, colsToMap=colsToMap, prefix="siteFeat")

  return(mae)
}
