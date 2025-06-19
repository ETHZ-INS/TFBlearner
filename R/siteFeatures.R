.getSeqFeatures <- function(refRanges,
                            phast=phastCons100way.UCSC.hg38,
                            genome=BSgenome.Hsapiens.UCSC.hg38)
{
  # CpG content
  seqlevelsStyle(refRanges) <- "UCSC"

  seqs <- Biostrings::getSeq(genome, refRanges)
  diFreq <- Biostrings::dinucleotideFrequency(x=seqs, as.prob=TRUE)
  cpgDens <- diFreq[,"CG"]+diFreq[,"GC"]
  gcCont <- Biostrings::letterFrequency(x=seqs, letters="GC", as.prob=TRUE)

  # PhastCon scores
  phastScores <- gscores(phast, refRanges, scores.only=TRUE)
  phastScores <- phastScores$default
  phastScores <- fifelse(is.na(phastScores), 0, phastScores)

  seqFeat <- list(phast_scores=Matrix::Matrix(phastScores, ncol=1),
                  cpg_density=Matrix::Matrix(cpgDens, ncol=1),
                  gc_cont=Matrix::Matrix(gcCont))

  names(seqFeat) <- c(CONSSCOREFEATNAME, CPGDENSEATNAME, GCCONTFEATNAME)

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
#' @importFrom Biostrings getSeq dinucleotideFrequency letterFrequency
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
    names(wi) <- WIDTHFEATNAME
    featMats <- append(featMats, wi)
  }

  if(!is.null(annoData) & "Annot" %in% features){

    if(length(scoreCols) < length(annoData)){
      message("Fewer colnames (scoreCols) to be aggregated provided, than
             annoData elements. Assuming all columns of elements of annoData
             are named the same")
      scoreCols <- rep(scoreCols, length(annoData))
    }

    annotFeats <- mapply(function(ad, an, SCORECOL, ...){
      annoDt <- .processData(ad, readAll=TRUE)
      annoDt$type <- an
      annoFeat <- genomicRangesMapping(coords,
                                       annoDt,
                                       byCols=c("type"),
                                       SCORECOL=SCORECOL,
                                       aggregationFun=aggregationFun,
                                       ...)},
      annoData, names(annoData), scoreCols, ...)

    names(annotFeats) <- names(annoData)
    featMats <- append(featMats, annotFeats)
  }

  # get gc content of promoter coordinates
  if(ATACPROMEXP %in% names(experiments(mae))){
    # prune to standard chromosomes
    promCoords <- rowRanges(mae[[ATACPROMEXP]])
    mae[[ATACPROMEXP]] <- keepStandardChromosomes(mae[[ATACPROMEXP]],
                                                   pruning.mode="coarse")

    promSeqs <- Biostrings::getSeq(genome, promCoords)
    gcProm <- Biostrings::letterFrequency(x=promSeqs, letters="GC", as.prob=TRUE)
    rowData(mae[[ATACPROMEXP]])[[GCCONTFEATNAME]] <- gcProm[,1]
  }

  names(featMats) <- paste(SITEFEAT, names(featMats), sep="_")
  seSiteFeat <- SummarizedExperiment(assays=featMats, rowRanges=coords)
  colnames(seSiteFeat) <- "all"
  colData(seSiteFeat)[[FEATTYPECOL]] <- SITEFEAT
  colsToMap <- unique(sampleMap(mae)$primary)
  mae <- .addFeatures(mae, seSiteFeat, colsToMap=colsToMap, prefix=SITEFEAT)

  return(mae)
}
