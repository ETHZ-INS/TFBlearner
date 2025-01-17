# TODO: could define generic to also accept data.table / data.frame / GRanges as list elements of input
# Could be useful for ChIP-data. / motif-data

# simplified import function => define generic
.processData <- function(data, readAll=FALSE, shift=FALSE, seqLevelStyle="UCSC"){
  if(is.character(data)){
    if(grepl(".bam", basename(data), fixed=TRUE))
    {
      param <- Rsamtools::ScanBamParam(what=c('pos', 'qwidth', 'isize'))
      readPairs <- GenomicAlignments::readGAlignmentPairs(data, param=param)

      # get fragment coordinates from read pairs
      seqDat <- GRanges(seqnames(readPairs@first),
                        IRanges(start=pmin(GenomicAlignments::start(readPairs@first),
                                           GenomicAlignments::start(readPairs@last)),
                                end=pmax(GenomicAlignments::end(readPairs@first),
                                         GenomicAlignments::end(readPairs@last))),
                        strand=GenomicAlignments::strand(readPairs))
      seqDat <- granges(seqDat, use.mcols=TRUE)
      seqDat <- as.data.table(seqDat)
      setnames(seqDat, c("seqnames"), c("chr"))
    }
    else if(grepl(".bed", basename(data), fixed=TRUE)){
      if(readAll) seqDat <- fread(data)
      else{
        seqDat <- fread(data, select=c(1:3,6),
                        col.names=c("chr", "start", "end", "strand"),
                        stringsAsFactors=TRUE)}
    }
    else if(grepl(".tsv", basename(data), fixed=TRUE)){
      if(readAll) seqDat <- fread(data)
      else{
        seqDat <- fread(data, select=c(1:3),
                        col.names=c("chr", "start", "end"),
                        stringsAsFactors=TRUE)}
    }
    else if(grepl(".rds", basename(data), fixed=TRUE)){
      seqDat <- as.data.table(readRDS(data))
      if("seqnames" %in% colnames(seqDat)) setnames(seqDat, "seqnames", "chr")
    }
  }
  else{
    seqDat <- as.data.table(data)
    if("seqnames" %in% colnames(seqDat)) setnames(seqDat, "seqnames", "chr")
  }

  # Match seqlevelstyle to reference
  if(sum(grepl("chr", seqDat$chr))==0 & seqLevelStyle=="UCSC" |
     sum(grepl("chr", seqDat$chr))>0 & seqLevelStyle=="NCBI"){
    seqDat <- makeGRangesFromDataFrame(as.data.frame(seqDat),
                                       keep.extra.columns=TRUE)
    seqlevelsStyle(seqDat) <- seqLevelStyle
    seqDat <- as.data.table(seqDat)
    setnames(seqDat, "seqnames", "chr")}

  # Insert ATAC shift
  if(shift & "strand" %in% colnames(seqDat))
  {
    seqDat[, start:=fifelse(strand=="+", start+4, start)]
    seqDat[, end:=fifelse(strand=="-", end-5, end)]
  }
  else if(shift){
    warning("Did not shift as no column named strand was not found")
  }

  seqDat[, start:=as.integer(start)]
  seqDat[, end:=as.integer(end)]

  return(seqDat)
}

.getType <- function(atacFrag, cuts=c("nucleosome_free"=0,
                                      "mononucleosome"=120,
                                      "dinucleosome"=300,
                                      "multinucleosome"=500),
                     label=FALSE) {
  atacFrag[,width:=end-start]
  if(label){
    labels <- names(cuts)
  }
  else{
    labels <- FALSE
  }
  atacFrag[,frag_type:=cut(width, breaks=c(cuts, Inf), labels=labels,
                      right=TRUE, include.lowest=TRUE)]
  return(atacFrag)
}

.mapSeqData <- function(data, refCoords, type=c("ATAC", "ChIP", "Motif"),
                        annoCol="context",
                        scoreCol="score",
                        weightCol=NULL,
                        isUncertainCol=NULL,
                        aggregationFun=max,
                        shift=FALSE,
                        BPPARAM=SerialParam()){
  type <- match.arg(type, choices=c("ATAC", "ChIP", "Motif"))

  if(is.data.table(data) | is.matrix(data) |
     is(data, 'sparseMatrix') | is.data.frame(data) |
     is(data, "DelayedMatrix")){

    if(nrow(data)==length(refCoords))
    {
      message("Reference coordinates and the number of rows of the provided
               data are the same - assuming data is already mapped")

      dataColData <- data.table(names=colnames(data))
      if(type=="Motif"){
        setnames(dataColData, "names", "motif")
      }
      else{
        setnames(dataColData, "names", annoCol)
      }
      se <- SummarizedExperiment(assays=data,
                                 rowRanges=refCoords,
                                 colData=dataColData)

    }
    else{
      data <- list(data)
      names(data) <- "unnamed"
      message("Not tested yet")
    }
  }

  if(type=="ATAC"){
    mappedSe <- .mapAtacData(data, refCoords,
                             annoCol=annoCol, shift=shift, BPPARAM=BPPARAM)
  }
  else if(type=="ChIP"){
    mappedSe <- .mapChIPData(data, refCoords,
                             weightCol=weightCol,
                             isUncertainCol=isUncertainCol,
                             aggregationFun=aggregationFun,
                             BPPARAM=BPPARAM)
  }
  else if(type=="Motif"){
    mappedSe <- .mapMotifData(data, refCoords, scoreCol=scoreCol,
                              aggregationFun=aggregationFun,
                              BPPARAM=BPPARAM)
  }

  return(mappedSe)
}

.mapMotifData <- function(data,
                          refCoords,
                          scoreCol="score",
                          aggregationFun=max,
                          BPPARAM=SerialParam()){
  threads <- floor(getDTthreads())/BPPARAM$workers

  motifScores <- BiocParallel::bpmapply(function(d, name, refCoords, scoreCol,
                                                 aggregationFun, threads){
    motifScore <- .processData(d, readAll=TRUE, shift=FALSE,
                               seqLevelStyle=seqlevelsStyle(refCoords))
    motifScore$motif_name <- name

    motifScore <- genomicRangesMapping(refCoords,
                                       motifScore,
                                       scoreCol=scoreCol,
                                       byCols="motif_name",
                                       aggregationFun=aggregationFun,
                                       BPPARAM=SerialParam())
    return(motifScore)
  }, data, names(data),
  MoreArgs=list(refCoords=refCoords,
                scoreCol=scoreCol,
                aggregationFun=aggregationFun,
                thread=threads),
  BPPARAM=BPPARAM)

  motifScores <- Reduce("cbind", motifScores[-1], motifScores[[1]])

  motifColData <- data.table(motif=colnames(motifScores))
  motifSe <- SummarizedExperiment(assays=list(match_scores=motifScores),
                                  rowRanges=refCoords,
                                  colData=motifColData)
}

.mapAtacData <- function(data,
                         refCoords,
                         annoCol="context",
                         shift=TRUE,
                         BPPARAM=SerialParam())
{
  threads <- floor(getDTthreads())/BPPARAM$workers

  # looped processing of ATAC data
  #typeCountNames <- c("total_ov_counts", "type_ov_counts",
  #                    "total_ins_counts","type_ins_counts")
  atacCounts <- BiocParallel::bplapply(data, function(d, refCoords,
                                                      shift, threads){

    setDTthreads(threads)
    atacFrag <- .processData(d, shift=shift,
                             seqLevelStyle=seqlevelsStyle(refCoords))
    atacFrag <- .getType(atacFrag, label=TRUE)

    # inserts counts
    atacFrag <- split(atacFrag, by="frag_type")
    atacIns <- suppressMessages(lapply(atacFrag, getInsertionProfiles,
                                       refCoords, margin=0, calcProfile=FALSE,
                                       shift=FALSE))

    atacIns <- lapply(atacIns, `[[`, "motifScores")
    atacFrag <- rbindlist(atacFrag)
    atacIns <- rbindlist(atacIns, idcol="frag_type")

    atacTypeOvs <- genomicRangesMapping(refCoords,
                                         atacFrag[,c("chr", "start",
                                                     "end", "frag_type"),
                                                  with=FALSE],
                                         byCols="frag_type",
                                         BPPARAM=SerialParam())

    atacTotalOvs <- Matrix::Matrix(Matrix::rowSums(atacTypeOvs), ncol=1)
    colnames(atacTotalOvs) <- "total_overlaps"
    typeNames <- colnames(atacTypeOvs)
    atacTypeOvs <- lapply(typeNames,
                           function(col) atacTypeOvs[,col, drop=FALSE])
    names(atacTypeOvs) <- paste(typeNames, "overlaps", sep="_")

    atacTypeIns <- genomicRangesMapping(refCoords,
                                         atacIns[,c("chr", "start",
                                                    "end", "frag_type"),
                                                 with=FALSE],
                                         byCols="frag_type",
                                         BPPARAM=SerialParam())
    atacTotalIns <- Matrix::Matrix(Matrix::rowSums(atacTypeIns), ncol=1)
    colnames(atacTotalIns) <- "total_inserts"
    atacTypeIns <- lapply(typeNames,
                          function(col) atacTypeIns[,col, drop=FALSE])
    names(atacTypeIns) <- paste(typeNames, "inserts", sep="_")

    atacAssays <- c(list("total_overlaps"=atacTotalOvs),
                    atacTypeOvs,
                    list("total_inserts"=atacTotalIns),
                    atacTypeIns)
    return(atacAssays)},
    refCoords, shift, threads, BPPARAM=BPPARAM)

  typeCountNames <- names(atacCounts[[1]])

  # merge replicates
  annoNames <- unique(names(atacCounts))
  atacCounts <- lapply(annoNames, function(anno){
    typeCounts <- unlist(atacCounts[names(atacCounts)==anno], recursive=FALSE)
    names(typeCounts) <- unlist(tstrsplit(names(typeCounts), split=".",
                                          fixed=TRUE, keep=2))

    typeCounts <- lapply(typeCountNames, function(type){
      counts <- typeCounts[names(typeCounts)==type]
      counts <- Reduce("cbind", counts[-1], counts[[1]])
      #column normalize?: counts <- (counts/(colSums(counts)+1e-04))
      counts <- Matrix::Matrix(Matrix::rowSums(counts), ncol=1)
    })
    names(typeCounts) <- typeCountNames
    typeCounts
  })
  names(atacCounts) <- annoNames

  # convert to assay matrices
  atacCounts <- unlist(atacCounts, recursive=TRUE)
  atacAssays <- lapply(typeCountNames, function(type){
    assayMat <- atacCounts[grepl(type, names(atacCounts))]
    cn <- unlist(tstrsplit(names(assayMat), split=".", keep=1, fixed=TRUE))
    assayMat <- Reduce("cbind", assayMat[-1], assayMat[[1]])
    colnames(assayMat) <- cn
    assayMat
  })
  names(atacAssays) <- typeCountNames

  atacColData <- data.table(annoNames)
  colnames(atacColData) <- annoCol
  dataOrigin <- lapply(data, function(d) {
    if(is.character(d))
    {
      return(d)
    }
    else return(NULL)})
  names(dataOrigin) <- names(data)
  atacColData[,origin:=lapply(get(annoCol), function(x) data[names(data)==x])]

  atacSe <- SummarizedExperiment(assays=atacAssays,
                                 rowRanges=refCoords,
                                 colData=atacColData)
  return(atacSe)
}

# data named list with context_tfName
.mapChIPData <- function(data,
                         refCoords,
                         aggregationFun=max,
                         annoCol="context",
                         weightCol=NULL,
                         isUncertainCol=NULL,
                         BPPARAM=SerialParam()){

  threads <- floor(getDTthreads())/BPPARAM$workers

  chIPPeaks <- BiocParallel::bpmapply(function(d, comb, refCoords, threads){

    chIPPeaks <- .processData(d, readAll=TRUE, shift=FALSE,
                              seqLevelStyle=seqlevelsStyle(refCoords))

    if(!is.null(isUncertainCol)){
      #chIPPeaks[,is_uncertain:=fifelse(is.na(is_uncertain), TRUE, is_uncertain)]
      byCols <- isUncertainCol}
    else{
      chIPPeaks$combination <- comb
      byCols <- "combination"
    }

    if(is.null(weightCol)){
      isPeak <- function(x){
        if(sum(x>0)>0){
          return(1)}
        else{
         return(0)}
      }
      aggregationFun <- isPeak
    }

    # map chIPPeaks
    chIPPeaks <- genomicRangesMapping(refCoords, chIPPeaks,
                                       aggregationFun=aggregationFun,
                                       scoreCol=weightCol,
                                       byCols=byCols,
                                       BPPARAM=SerialParam())
    # mark uncertains
    if(!is.null(isUncertainCol)){
      chIPPeaks <- as(chIPPeaks, "TsparseMatrix")[,c("FALSE", "TRUE")]
      chIPPeaks <- Matrix::Matrix(fifelse(chIPPeaks[,"TRUE", drop=TRUE]>0,
                                          -chIPPeaks[,"TRUE", drop=TRUE],
                                          chIPPeaks[,"FALSE", drop=TRUE]),
                                  ncol=1)}

    colnames(chIPPeaks) <- comb

    return(chIPPeaks)
  }, data, names(data), MoreArgs=list(refCoords=refCoords,
                                      threads=threads),
  BPPARAM=BPPARAM)

  # merge replicates:
  # Max score across replicates and uncertainty flag removed if exists in one
  repMergeFun <- function(x){
    if(sum(x>0, na.rm=TRUE)>1){ # replicated in more than one
      w <- 1
    }
    else if(sum(x>0, na.rm=TRUE)==1){ #replicated in one
      w <- max(x)
    }
    else if(sum(x<0, na.rm=TRUE)>0){ # uncertain in one
      w <- -max(abs(x))
    }
    else{ # found in none
      w <- 0
    }
    return(w)
  }

  annoNames <- unique(names(chIPPeaks))
  chIPPeaks <- lapply(annoNames, function(anno){
    repPeaks <- chIPPeaks[names(chIPPeaks)==anno]
    repPeaks <- as(Reduce("cbind", repPeaks[-1], repPeaks[[1]]), "TsparseMatrix")
    repPeaksDt <- data.table(i=repPeaks@i, j=repPeaks@j, x=repPeaks@x)
    repPeaksDt <- repPeaksDt[,.(x=repMergeFun(x)), by=i]
    repPeaksDt <- subset(repPeaksDt, x!=0)
    aggPeaksMat <- sparseMatrix(i=repPeaksDt$i+1, j=rep(1,nrow(repPeaksDt)),
                                x=repPeaksDt$x,
                                dims=c(length(refCoords), 1))
  })

  chIPPeaks <- Reduce("cbind", chIPPeaks[-1], chIPPeaks[[1]])
  colnames(chIPPeaks) <- annoNames

  chIPColData <- data.table(combination=annoNames)
  chIPColData[,c(annoCol, "tf_name"):=tstrsplit(combination, split="_")]
  chIPColData[,origin:=lapply(combination, function(x) data[names(data)==x])]

  chIPSe <- SummarizedExperiment(assays=list(peaks=chIPPeaks),
                                 rowRanges=refCoords,
                                 colData=chIPColData)

  return(chIPSe)
}
