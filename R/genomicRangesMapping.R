#' Mapping & aggregation of modalities with genomic coordinates to a set of reference
#' coordinates.
#'
#' Convencience function for mapping different modality scores with
#' cell type and further labels such as tfs to reference coordinates. The resulting
#' table will have dimension ref coords x byCols (or ref coord x byCols`[1]` x byCols`[2]`).
#' ByCols can be for instance cell type labels and/or transcription factor names.
#'
#' @name genomicRangesMapping
#' @param refRanges GRanges object with reference coordinates
#' @param assayTable  List of [GenomicRanges::GRanges-class], [data.table::data.table], data.frames or paths to .bed /. bam files
#' containing assay data (e.g. fragments, motif scores, peaks) to be aggregated across the reference ranges provided.
#' Need to contain genomic coordinates (e.g. a chr/seqnames, start and end column).
#' @param byCols Variables across which the assays will be aggregated.
#' Will be the columns of the resulting [Matrix::Matrix-class]. If byCols is a vector with two elements,
#' the first one will constitute the list elements, the second the columns of the matrices being the list elements.
#' @param scoreCol name of the score column (e.g. motif matching scores, atac fragment counts) to be aggregated.
#' If it is NULL, the number of assay ranges overlapping the reference ranges will be counted.
#' @param aggregationFun function (e.g. mean, median, sum) used to aggregate.
#' If it is NULL, the number of assay ranges overlapping the reference ranges will be counted.
#' @param minoverlap Minimal overlap between refRanges and the assay ranges
#' Passed to [GenomicRanges::findOverlaps()]
#' @param shift If Tn5 insertion bias should be considered (only if strand column is provided).
#' @param BPPARAM Parallel back-end to be used. Passed to [BiocParallel::bplapply()].
#' @return [Matrix::Matrix-class] or list of Matrices with rows corresponding to the reference ranges and columns (and list elements) to byCols.
#' @import data.table
#' @import Matrix
#' @importFrom GenomicRanges findOverlaps GRanges GRangesList
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom BiocParallel bplapply MulticoreParam SerialParam SnowParam
#' @importFrom GenomicAlignments readGAlignmentPairs start end strand
#' @importFrom Rsamtools ScanBamParam
#' @export
genomicRangesMapping <- function(refRanges,
                                 assayTable,
                                 byCols=c("tf_name",
                                          "cellular_context"),
                                 scoreCol=NULL,
                                 aggregationFun=NULL,
                                 minoverlap=1,
                                 shift=FALSE,
                                 BPPARAM=SerialParam()){

  # TODO: - add warning for integer overflows - data.table size
  assayTable <- .processData(copy(assayTable), readAll=TRUE, shift=shift,
                             seqLevelStyle=seqlevelsStyle(refRanges))

  seqNamesCol <- "chr"
  startCol <- "start"
  endCol <- "end"

  if(sum(!(byCols %in% colnames(assayTable)))>0 | is.null(byCols)){
    stop("byCols needed to be provided and column names of assayTable.")
  }

  # attribute generic names to dimensionalities
  if(length(byCols)==2)
  {
    setnames(assayTable, byCols, c("col_depth", "col_width"))
    byCols <- c("col_depth", "col_width")
    multiTf <- TRUE
  }
  else
  {
    setnames(assayTable, byCols, c("col_width"))
    byCols <- c("col_width")
    multiTf <- FALSE
  }

  # get dimensions of tables
  nRefs <- length(refRanges)
  colsWidth <- unique(assayTable$col_width)
  nColsWidth <- length(colsWidth)
  colsDepth <- unique(assayTable$col_depth)

  # convert to integer for speed-up
  levels <- unique(assayTable$col_width)
  assayTable[,col_width:=as.integer(factor(assayTable$col_width, levels=levels))]

  # convert to GRanges for faster overlap finding
  suppressWarnings(assayTable$width <- NULL)
  suppressWarnings(assayTable$strand <- NULL)
  assayRanges <- makeGRangesFromDataFrame(as.data.frame(assayTable),
                                          keep.extra.columns=TRUE,
                                          seqnames.field=seqNamesCol,
                                          start.field=startCol,
                                          end.field=endCol,
                                          ignore.strand=TRUE)

  # find overlaps with ref. coordinates
  overlapTable <- as.data.table(GenomicRanges::findOverlaps(refRanges,
                                                            assayRanges,
                                                            type="any",
                                                            minoverlap=minoverlap,
                                                            ignore.strand=TRUE))
  rm(refRanges, assayRanges)

  # retrieve tf and cell type ids
  overlapTable <- cbind(overlapTable$queryHits,
                        assayTable[overlapTable$subjectHits,
                                   c(byCols, scoreCol),
                                   with=FALSE])

  threads <- floor(getDTthreads())/BPPARAM$workers

  if(multiTf)
  {
    setkey(overlapTable, V1, col_width)
    if(!is.null(scoreCol)) setnames(overlapTable, scoreCol, "scoreCol")
    overlapTable <- split(overlapTable, by=c("col_depth"))

    overlapTable <- BiocParallel::bplapply(overlapTable, function(table,
                                                                  scoreCol,
                                                                  aggregationFun,
                                                                  nRefs,
                                                                  nColsWidth,
                                                                  threads){

      data.table::setDTthreads(threads)

      if(is.null(scoreCol) | is.null(aggregationFun)){
        table <- table[,.(value=.N),
                       by=c("V1", "col_width")]}
      else{
        table <- table[,.(value=aggregationFun(scoreCol)),
                       by=c("V1", "col_width")]}

      # convert to sparse matrix
      table <- sparseMatrix(i=table$V1,
                            j=as.integer(table$col_width),
                            dims=c(nRefs, nColsWidth),
                            x=table$value)
      colnames(table) <- levels
      return(table)}, scoreCol=scoreCol, aggregationFun=aggregationFun,
                      nRefs=nRefs, nColsWidth=nColsWidth, threads=threads,
      BPPARAM=BPPARAM)
  }
  else
  {
    # setkeys for speed-up
    overlapTable[,V1:=as.integer(V1)]
    setkey(overlapTable, col_width, V1)

    # overlap with ref. coordinates
    if(is.null(scoreCol) | is.null(aggregationFun)){
      overlapTable <- overlapTable[,.(scoreCol=.N),
                                   by=c("col_width", "V1")]}
    else{
      setnames(overlapTable, scoreCol, "scoreCol")
      overlapTable <- overlapTable[,.(scoreCol=aggregationFun(scoreCol)),
                                   by=c("col_width", "V1")]}

    # convert to sparse matrix
    overlapTable <- Matrix::sparseMatrix(i=overlapTable$V1,
                                         j=overlapTable$col_width,
                                         dims=c(nRefs, nColsWidth),
                                         x=overlapTable$scoreCol)

    colnames(overlapTable) <- levels
  }

  # add combinations with zero overlaps
  missingDepthCols <- setdiff(colsDepth, names(overlapTable))
  if(length(missingDepthCols)>0){
    missingMat <- Matrix(0,nrow=nRefs, ncol=nColsWidth, doDiag=FALSE)
    colnames(missingMat) <- colsWidth
    missingTables <- replicate(length(missingDepthCols),
                               missingMat)
    names(missingTables) <- missingDepthCols
    overlapTable <- c(overlapTable, missingTables)
  }

  gc()
  return(overlapTable)
}
