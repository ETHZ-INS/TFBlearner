#' prepMotifs
#'
#' Prepare match coordinates files for a list of motifs.
#'
#' @param refCoords Coordinates across which to identify matches. This should be
#'   the same as will be inputted to \code{\link{prepData}}.
#' @param motifs A named list of motifs in `universalmotif` format, or an object
#'   of class `PWMatrixList`. The motif names should be gene names (or whatever
#'   was used to prepare and annotate the ChIP data).
#' @param outputFolder A folder where the coordinates of the motif matches will
#'   be stored (one file per motif will be created).
#' @param genome The genome to use, either a BSgenome or FaFile object.
#' @param resume Logical, whether to only prepare files that aren't already
#'   there.
#' @param topPlusCoOnly Logical, whether to report only the top match per
#'   `refCoords`, along with the number of other matches. If FALSE, will report
#'   the coordinates of all matches.
#' @param BPPARAM BiocParallel BPPARAM object for parallelization.
#'
#' @details
#' This will save a Rds file for each motif, containing a GRanges of its matches
#' within the regions defined by `refCoords`. If `topPlusCoOnly=TRUE`, only the
#' top match per `refCoords` is saved, along with the number of other
#' significant matches in the regions. If `topPlusCoOnly=TRUE`, all matches
#' passing threshold are saved, and the top per region is flagged using an
#' `isTop` column. Note that to identify the top match per region, a very low
#' threshold is used in order for most regions to end up with a match.
#'
#' @author Pierre-Luc Germain
#' @returns A named vector of paths to the files created, which can be used for
#'   input to \code{\link{prepData}}.
#' @importFrom motifmatchr matchMotifs
#' @importFrom GenomicRanges findOverlaps countOverlaps
#' @importFrom IRanges overlapsAny
#' @importFrom GenomicRanges score score<-
#' @export
prepMotifs <- function(refCoords, motifs, outputFolder, genome, resume=FALSE,
                       topPlusCoOnly=TRUE, BPPARAM=SerialParam(progress=TRUE)){

  # check inputs
  stopifnot(is(refCoords, "GRanges"))
  stopifnot(is(genome, "BSgenome") || is(genome, "FaFile"))
  if(!dir.exists(outputFolder))
    stop("`outputFolder` does not seem to exist.")

  stopifnot(!is.null(names(motifs)))
  if(is.list(motifs) && all(sapply(motifs, is, class2="universalmotif"))){
    motifs <- lapply(motifs, universalmotif::convert_motifs,
                     class="TFBSTools-PWMatrix")
    motifs <- do.call(TFBSTools::PWMatrixList,motifs)
  } else if(!is(motifs, "PWMatrixList")){
    stop("`motifs` should be a named list of universalmotif objects or a ",
         "PWMatrixList object.")
  }

  # run scan and save results for each motif
  res <- bplapply(names(motifs), BPPARAM=BPPARAM, FUN=function(xn){
    x <- motifs[[xn]]
    fn <- file.path(outputFolder,paste0(make.names(xn),".rds"))
    if(!resume || !file.exists(fn)){
      a <- .doOneScan(x, refCoords, genome=genome, topPlusCoOnly=topPlusCoOnly)
      saveRDS(a, fn)
    }
    fn
  })
  setNames(unlist(res), names(motifs))
}

.doOneScan <- function(mo, dhs, genome, lowp=1e-02, p=5e-05, minOL=3L,
                       topPlusCoOnly=TRUE){
  # low-confidence matches
  po <- matchMotifs(mo, subject=dhs, genome=genome, out="positions",
                    p.cutoff=lowp)[[1]]
  # trace back the original DHS behind each
  po$orig <- to(findOverlaps(po, dhs, minoverlap=minOL))
  # keep only the strongest hit per DHS
  po <- po[order(seqnames(po), -po$score)]
  po <- po[!duplicated(po$orig)]
  if(!topPlusCoOnly) po$isTop <- TRUE

  # find all below the default cutoff
  po2 <- matchMotifs(mo, subject=dhs, genome=genome, out="positions",
                     p.cutoff=p)[[1]]
  names(po2@elementMetadata) <- SCORECOL

  if(!topPlusCoOnly) po2$isTop <- FALSE
  # remove overlapping motifs
  po2 <- po2[!overlapsAny(po2, po, minoverlap=minOL, ignore.strand=TRUE)]
  if(topPlusCoOnly){
    po@elementMetadata[[COOCCURRENCECOL]] <- countOverlaps(dhs[po$orig], po2,
                                                          minoverlap=minOL,
                                                          ignore.strand=TRUE)
    po$orig <- NULL
  }else{
    po$orig <- NULL
    po <- c(po,po2)
  }
  score(po) <- round(score(po),2)
  sort(po)
}
