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
  if((sum(grepl("chr", seqDat$chr))==0 & seqLevelStyle=="UCSC") |
     (sum(grepl("chr", seqDat$chr))>0 & seqLevelStyle=="NCBI")){
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
                        saveHdf5=FALSE,
                        outDir=NULL,
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
      assayList <- list(data)
      names(assayList) <- type
      se <- SummarizedExperiment(assays=assayList,
                                 rowRanges=refCoords,
                                 colData=dataColData)
      return(se)
    }
    else{
      data <- list(data)
      names(data) <- "unnamed"
      message("Not tested yet")
    }
  }

  if(type=="ATAC"){
    mappedSe <- .mapAtacData(data, refCoords,
                             annoCol=annoCol, shift=shift,
                             saveHdf5=saveHdf5, outDir=outDir,
                             BPPARAM=BPPARAM)
  }
  else if(type=="ChIP"){
    mappedSe <- .mapChIPData(data, refCoords,
                             weightCol=weightCol,
                             isUncertainCol=isUncertainCol,
                             aggregationFun=aggregationFun,
                             saveHdf5=saveHdf5, outDir=outDir,
                             BPPARAM=BPPARAM)
  }
  else if(type=="Motif"){
    mappedSe <- .mapMotifData(data, refCoords, scoreCol=scoreCol,
                              aggregationFun=aggregationFun,
                              saveHdf5=saveHdf5, outDir=outDir,
                              BPPARAM=BPPARAM)
  }

  return(mappedSe)
}

.mapMotifData <- function(data,
                          refCoords,
                          scoreCol="score",
                          aggregationFun=max,
                          saveHdf5=FALSE,
                          outDir=NULL,
                          BPPARAM=SerialParam()){
  threads <- floor(getDTthreads())/BPPARAM$workers

  # reshape list
  colNames <- unique(names(data))
  data <- lapply(colNames, function(name) data[names(data)==name])
  names(data) <- colNames

  if(saveHdf5)
  {
    if(is.null(outDir)) outDir <- getwd() # .
    fileName <- "Motifs_mapped"
    hdf5FileName <- file.path(outDir, fileName)

    if(file.exists(hdf5FileName)){
      h5delete(file=hdf5FileName, name="match_scores")
      file.remove(hdf5FileName)
    }

    h5createFile(hdf5FileName)
    h5createDataset(file=hdf5FileName, dataset="match_scores",
                    dims=c(length(refCoords), length(colNames)),
                    storage.mode="double",
                    chunk=c(1e5, length(colNames)))
  }
  else{
    hdf5FileName <- NULL
  }

  motifScores <- BiocParallel::bplapply(data,
                                        function(d, refCoords, scoreCol,
                                                 aggregationFun, threads,
                                                 saveHdf5, hdf5FileName,
                                                 colNames){
    data.table::setDTthreads(threads)
    name <- unique(names(d))
    motifScore <- lapply(d, .processData, readAll=TRUE, shift=FALSE,
                         seqLevelStyle=seqlevelsStyle(refCoords))
    motifScore <- rbindlist(motifScore)
    motifScore$motif_name <- name

    motifScore <- genomicRangesMapping(refCoords,
                                       motifScore,
                                       scoreCol=scoreCol,
                                       byCols="motif_name",
                                       aggregationFun=aggregationFun,
                                       BPPARAM=SerialParam())
    maxScore <- max(motifScore@x)

    if(saveHdf5){
      i <- which(colNames==name)
      h5write(as.matrix(motifScore),
              file=hdf5FileName,
              name="match_scores",
              createnewfile=FALSE,
              index=list(1:length(refCoords), i))
      motifScore <- head(motifScore, n=1)
    }

    return(list(motifScore, maxScore))
  }, refCoords=refCoords, scoreCol=scoreCol,
     aggregationFun=aggregationFun, thread=threads,
     saveHdf5=saveHdf5, hdf5FileName=hdf5FileName, colNames=colNames,
     BPPARAM=BPPARAM)

  maxScores <- unlist(lapply(motifScores, `[[`, 2))
  motifScores <- lapply(motifScores, `[[`, 1)

  if(saveHdf5){
    H5close()
    motifScores <- HDF5Array(filepath=hdf5FileName, name="match_scores")
  }
  else{
    motifScores <- Reduce("cbind", motifScores[-1], motifScores[[1]])
  }

  colnames(motifScores) <- colNames
  motifColData <- data.table(motif=colnames(motifScores),
                             max_score=maxScores)
  motifSe <- SummarizedExperiment(assays=list(match_scores=motifScores),
                                  rowRanges=refCoords,
                                  colData=motifColData)
}

.mapAtacData <- function(data,
                         refCoords,
                         annoCol="context",
                         shift=TRUE,
                         saveHdf5=FALSE,
                         outDir=NULL,
                         BPPARAM=SerialParam())
{
  threads <- floor(getDTthreads())/BPPARAM$workers

  # reshape list
  colNames <- unique(names(data))
  data <- lapply(colNames, function(name) data[names(data)==name])
  names(data) <- colNames

  if(saveHdf5)
  {
    if(is.null(outDir)) outDir <- getwd() # .
    fileName <- "ATAC_mapped"
    hdf5FileName <- file.path(outDir, fileName)

    datasets <- c("total_overlaps", "nucleosome_free_overlaps",
                  "dinucleosome_overlaps", "mononucleosome_overlaps",
                  "multinucleosome_overlaps", "total_inserts",
                  "nucleosome_free_inserts", "dinucleosome_inserts",
                  "mononucleosome_inserts", "multinucleosome_inserts")

    if(file.exists(hdf5FileName)){
      lapply(datasets, h5delete, file=hdf5FileName)
      file.remove(hdf5FileName)
    }
    h5createFile(hdf5FileName)
    lapply(datasets,h5createDataset, file=hdf5FileName,
           dims=c(length(refCoords), length(colNames)), storage.mode="integer",
           chunk=c(1e5, length(colNames)))
  }
  else{
    hdf5FileName <- NULL
  }

  # looped processing of ATAC data
  atacCounts <- BiocParallel::bplapply(data, function(d, refCoords,
                                                      shift, threads,
                                                      saveHdf5, hdf5FileName,
                                                      colNames){

    data.table::setDTthreads(threads)
    atacFrag <- lapply(d, .processData, shift=shift,
                       seqLevelStyle=seqlevelsStyle(refCoords))
    atacFrag <- rbindlist(atacFrag)
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
    if(saveHdf5){
      i <- which(colNames==unique(names(d)))
      lapply(names(atacAssays), function(assay){
        h5write(as.matrix(atacAssays[[assay]]),
                file=hdf5FileName,
                name=assay,
                createnewfile=FALSE,
                index=list(1:length(refCoords), i))

      })
      atacAssays <- lapply(atacAssays, head, n=1)
    }

    return(atacAssays)},
    refCoords, shift, threads, saveHdf5, hdf5FileName=hdf5FileName,
    colNames=names(data), BPPARAM=BPPARAM)

  typeCountNames <- names(atacCounts[[1]])

  if(saveHdf5){
    H5close()
    atacAssays <- lapply(typeCountNames, HDF5Array, filepath=hdf5FileName,
                         as.sparse=TRUE)
    atacAssays <- lapply(atacAssays, `colnames<-`, colNames)
    names(atacAssays) <- typeCountNames
  }
  else{
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
  }

  atacColData <- data.table(colNames)
  colnames(atacColData) <- annoCol
  dataOrigin <- lapply(data, function(d) {
    if(is.character(d))
    {
      return(d)
    }
    else return(NULL)})
  names(dataOrigin) <- colNames


  atacColData[,origin:=lapply(get(annoCol), function(x){
    ds <- unlist(data[names(data)==x])
    names(ds) <- unlist(tstrsplit(names(ds), split=".", keep=2, fixed=TRUE))
    ds <- lapply(ds, function(d) fifelse(is.character(d), d, NA))
    return(ds)})]

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
                         saveHdf5=FALSE,
                         outDir=NULL,
                         BPPARAM=SerialParam()){

  threads <- floor(getDTthreads())/BPPARAM$workers

  # reshape list
  colNames <- unique(names(data))
  data <- lapply(colNames, function(name) data[names(data)==name])
  names(data) <- colNames

  if(saveHdf5)
  {
    if(is.null(outDir)) outDir <- getwd() # .
    fileName <- "ChIP_mapped"
    hdf5FileName <- file.path(outDir, fileName)

    if(file.exists(hdf5FileName)){
      h5delete(file=hdf5FileName, name="peaks")
      file.remove(hdf5FileName)
    }

    h5createFile(hdf5FileName)
    h5createDataset(file=hdf5FileName, dataset="peaks",
                    dims=c(length(refCoords), length(colNames)),
                    storage.mode="double",
                    chunk=c(1e5, length(colNames)))
  }
  else{
    hdf5FileName <- NULL
  }

  chIPPeaks <- BiocParallel::bplapply(data, function(d, refCoords, threads,
                                                     saveHdf5, hdf5FileName,
                                                     colNames){

    data.table::setDTthreads(threads)

    comb <- unique(names(d))
    chIPPeaks <- lapply(d, .processData, readAll=TRUE, shift=FALSE,
                        seqLevelStyle=seqlevelsStyle(refCoords))
    chIPPeaks <- rbindlist(chIPPeaks)


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
      missCols <- setdiff(c("FALSE", "TRUE"), colnames(chIPPeaks))
      if(length(missCols)>0){
        missMat <- Matrix(0, nrow=length(refCoords), ncol=length(missCols))
        colnames(missMat) <- missCols
        chIPPeaks <- cbind(chIPPeaks, missMat)
      }

      chIPPeaks <- as(chIPPeaks, "TsparseMatrix")[,c("FALSE", "TRUE")]
      chIPPeaks <- Matrix::Matrix(fifelse(chIPPeaks[,"TRUE", drop=TRUE]>0,
                                          -chIPPeaks[,"TRUE", drop=TRUE],
                                          chIPPeaks[,"FALSE", drop=TRUE]),
                                  ncol=1)}

    colnames(chIPPeaks) <- comb
    if(saveHdf5){
      i <- which(colNames==comb)
      h5write(as.matrix(chIPPeaks),
                file=hdf5FileName,
                name="peaks",
                createnewfile=FALSE,
                index=list(1:length(refCoords), i))
      chIPPeaks <- head(chIPPeaks, n=1)
    }

    return(chIPPeaks)
  }, refCoords=refCoords, threads=threads, saveHdf5=saveHdf5,
     hdf5FileName=hdf5FileName, colNames=names(data), BPPARAM=BPPARAM)

  if(saveHdf5){
    H5close()
    chIPPeaks <- HDF5Array(filepath=hdf5FileName, name="peaks")
  }
  else{
    chIPPeaks <- Reduce("cbind", chIPPeaks[-1], chIPPeaks[[1]])
  }

  colnames(chIPPeaks) <- colNames
  chIPColData <- data.table(combination=colNames)
  chIPColData[,c(annoCol, "tf_name"):=tstrsplit(combination, split="_")]
  chIPColData[,origin:=lapply(combination, function(x){
    ds <- unlist(data[names(data)==x])
    names(ds) <- unlist(tstrsplit(names(ds), split=".", keep=2, fixed=TRUE))
    ds <- lapply(ds, function(d) fifelse(is.character(d), d, NA))
    return(ds)})]

  chIPSe <- SummarizedExperiment(assays=list(peaks=chIPPeaks),
                                 rowRanges=refCoords,
                                 colData=chIPColData)

  return(chIPSe)
}
