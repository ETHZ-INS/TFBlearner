.processData <- function(data, readAll=FALSE, shift=FALSE,
                         subSample=NULL, seqLevelStyle="UCSC"){
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

        readBed <- function(data){
          tryCatch(
            {
              seqDat <- fread(data, select=c(1:3,6),
                              col.names=c("chr", "start", "end", "strand"),
                              stringsAsFactors=TRUE)
              return(seqDat)},
            error = function(cond){
              seqDat <- fread(data, select=c(1:3),
                              col.names=c("chr", "start", "end"),
                              stringsAsFactors=TRUE)
              return(seqDat)
            })}
        seqDat <- readBed(data)
      }
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

  if(!is.null(subSample) & is.numeric(subSample)){
    message("Subsampling file")
    subSample <- as.integer(subSample)
    seqDat <- seqDat[sample(1:nrow(seqDat), min(nrow(seqDat), subSample)),]
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
  if(shift)
  {
    seqDat[, start:=start-4]
    seqDat[, end:=end-4]
  }
  else if(shift){
    warning("Did not shift as no column named strand was not found")
  }

  seqDat[, start:=as.integer(start)]
  seqDat[, end:=as.integer(end)]
  if("width" %in% colnames(seqDat)) seqDat$width <- NULL

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
  atacFrag$width <- NULL
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
                        fileName=NULL,
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

  if(is.null(fileName)) fileName <- paste(type, "mapped", sep="_")

  if(type=="ATAC"){
    mappedSe <- .mapAtacData(data, refCoords,
                             annoCol=annoCol, shift=shift,
                             saveHdf5=saveHdf5,
                             fileName=fileName,
                             outDir=outDir,
                             BPPARAM=BPPARAM)
  }
  else if(type=="ChIP"){
    mappedSe <- .mapChIPData(data, refCoords,
                             weightCol=weightCol,
                             isUncertainCol=isUncertainCol,
                             aggregationFun=aggregationFun,
                             saveHdf5=saveHdf5,
                             fileName=fileName,
                             outDir=outDir,
                             BPPARAM=BPPARAM)
  }
  else if(type=="Motif"){
    mappedSe <- .mapMotifData(data, refCoords, scoreCol=scoreCol,
                              aggregationFun=aggregationFun,
                              saveHdf5=saveHdf5,
                              fileName=fileName,
                              outDir=outDir,
                              BPPARAM=BPPARAM)
  }

  return(mappedSe)
}

.mapMotifData <- function(data,
                          refCoords,
                          scoreCol="score",
                          aggregationFun=max,
                          saveHdf5=FALSE,
                          fileName=NULL,
                          outDir=NULL,
                          BPPARAM=SerialParam()){
  threads <- floor(getDTthreads())/BPPARAM$workers

  # reshape list
  colNames <- unique(names(data))
  data <- lapply(colNames, function(name) data[names(data)==name])
  names(data) <- colNames

  motifScores <- BiocParallel::bplapply(data,
                                        function(d, refCoords, scoreCol,
                                                 aggregationFun, threads,
                                                 saveHdf5, outDir){
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
    motifScore <- .roundingCompression(motifScore, factor=1e4)
    maxScore <- max(motifScore@x)

    if(saveHdf5){
      dataList <- list(motifScore)
      names(dataList) <- matchAssayName
      .writeToHdf5(dataList,
                   paste0(file.path(outDir, name), ".h5"),
                   storage="integer")
      motifScore <- head(motifScore, n=1)
    }

    return(list(motifScore, maxScore))
  }, refCoords=refCoords, scoreCol=scoreCol,
     aggregationFun=aggregationFun, thread=threads,
     saveHdf5=saveHdf5, outDir=outDir,
     BPPARAM=BPPARAM)

  maxScores <- unlist(lapply(motifScores, `[[`, 2))
  motifScores <- lapply(motifScores, `[[`, 1)

  if(saveHdf5){
    hdf5FilesCollect <- unlist(lapply(colNames, function(f) paste0(file.path(outDir,f), ".h5")))
    motifScores <- .collectHdf5Files(hdf5FilesCollect,
                                   hdf5FileName=paste0(file.path(outDir, fileName), ".h5"),
                                   storage="integer", asSparse=FALSE)
    motifScores <- motifScores[[1]]
  }
  else{
    motifScores <- Reduce("cbind", motifScores[-1], motifScores[[1]])
  }

  colnames(motifScores) <- colNames
  motifColData <- data.table(colnames(motifScores), maxScores)
  colnames(motifColData) <- c(motifNameCol, maxScoreCol)

  assayList <- list(motifScores)
  names(assayList) <- matchAssayName
  motifSe <- SummarizedExperiment(assays=assayList,
                                  rowRanges=refCoords,
                                  colData=motifColData)

  return(motifSe)

}

.writeToHdf5 <- function(datasets, hdf5FileName, storage="integer"){
  fid <- H5Fcreate(hdf5FileName)
  on.exit(H5Fclose(fid))

  mapply(function(d, datasetName, storage, fid){
    h5createDataset(file=fid, dataset=datasetName,
                    dims=c(nrow(d), 1), storage.mode=storage,
                    chunk=c(min(1e6, nrow(d)),1),
                    level=0)
    h5writeDataset(as.matrix(d), h5loc=fid, name=datasetName,
                   level=0)},
    datasets, names(datasets),
    MoreArgs=list(storage=storage,
                  fid=fid))
  H5garbage_collect()
  return(TRUE)
}

# see: https://www.bioconductor.org/packages/devel/bioc/vignettes/rhdf5/inst/doc/practical_tips.html#3_Writing_in_parallel
.collectHdf5Files <- function(hdf5Files, hdf5FileName,
                              storage="integer", asSparse=TRUE){

  #TODO:  check requirement that all files have same dimension and datasets

  checkFile <- H5Fopen(hdf5Files[1]) # open first file
  info <- rhdf5::h5ls(checkFile)
  datasetNames <- info$name
  dims <- as.numeric(strsplit(info$dim, " x ")[[1]])
  H5Fclose(checkFile)

  if(file.exists(hdf5FileName)){
    file.remove(hdf5FileName)
  }

  # create combined object
  fid <- H5Fcreate(name=hdf5FileName)
  lapply(datasetNames, h5createDataset, file=fid,
         dims=c(dims[1], length(hdf5Files)), storage.mode=storage,
         level=6, chunk=c(min(1e6, dims[1]), 1))
  toFile <- H5Fopen(hdf5FileName)

  for(i in seq_len(length(hdf5Files))){
     file <- hdf5Files[i]
     origFile <- H5Fopen(file)
     datasets <- rhdf5::h5ls(origFile)$name
     H5Fclose(origFile)
     for(d in datasets){
       mat <- h5read(file, name=d)
       h5writeDataset(mat, h5loc=toFile, name=d,
                      index=list(1:dims[1], i), level=6)
     }
     file.remove(file)
  }
  H5Fclose(toFile)
  h5closeAll()
  H5garbage_collect()

  assays <- lapply(datasetNames, HDF5Array, filepath=hdf5FileName, as.sparse=asSparse)
  names(assays) <- datasetNames

  return(assays)
}

.mapAtacData <- function(data,
                         refCoords,
                         annoCol="context",
                         shift=TRUE,
                         saveHdf5=FALSE,
                         fileName=NULL,
                         outDir=NULL,
                         BPPARAM=SerialParam())
{
  threads <- floor(getDTthreads())/BPPARAM$workers

  # reshape list
  colNames <- unique(names(data))
  data <- lapply(colNames, function(name) data[names(data)==name])
  names(data) <- colNames

  # looped processing of ATAC data
  atacCounts <- BiocParallel::bplapply(data, function(d, refCoords,
                                                      shift, threads,
                                                      saveHdf5, outDir){

    data.table::setDTthreads(threads)

    atacFrag <- lapply(d, .processData, shift=shift,
                       seqLevelStyle=seqlevelsStyle(refCoords))
    atacFrag <- rbindlist(atacFrag)
    atacFrag <- .getType(atacFrag, label=TRUE)

    # inserts counts
    atacFrag <- split(atacFrag, by="frag_type")
    atacIns <- mapply(function(atacFragType, type){
      if(nrow(atacFragType)>0){
        ins <- suppressMessages({
          getInsertionProfiles(atacFragType, refCoords, margin=0,
                               calcProfile=FALSE, shift=FALSE)})
        ins <- ins[[retScoresName]]
        ins <- ins[,c("chr", "start", "end", insertFeatName),with=FALSE]}
      else{
        ins <- data.table(chr=character(), start=numeric(), end=numeric(),
                          insert_counts=numeric())
      }
      ins$frag_type <- factor(type)
      ins
    }, atacFrag, names(atacFrag), SIMPLIFY=FALSE)
    atacFrag <- rbindlist(atacFrag)
    atacIns <- rbindlist(atacIns)

    atacFrag$strand <- NULL
    atacTypeOvs <- genomicRangesMapping(refCoords,
                                        atacFrag,
                                        byCols="frag_type",
                                        BPPARAM=SerialParam())
    atacTotalOvs <- Matrix::Matrix(Matrix::rowSums(atacTypeOvs), ncol=1)
    colnames(atacTotalOvs) <- totalOverlapsName
    typeNames <- colnames(atacTypeOvs)
    atacTypeOvs <- lapply(typeNames,
                          function(col) atacTypeOvs[,col, drop=FALSE])
    names(atacTypeOvs) <- paste(typeNames, typeOverlapSuffix, sep="_")

    atacTypeIns <- genomicRangesMapping(refCoords,
                                        atacIns,
                                        scoreCol="insert_counts",
                                        byCols="frag_type",
                                        aggregationFun=sum,
                                        BPPARAM=SerialParam())
    atacTotalIns <- Matrix::Matrix(Matrix::rowSums(atacTypeIns), ncol=1)
    colnames(atacTotalIns) <- totalInsertsName
    atacTypeIns <- lapply(typeNames,
                          function(col) atacTypeIns[,col, drop=FALSE])
    names(atacTypeIns) <- paste(typeNames, typeInsertsSuffix, sep="_")

    atacTotalOvs <- list(atacTotalOvs)
    names(atacTotalOvs) <- totalOverlapsName
    atacTotalIns <- list(atacTotalIns)
    names(atacTotalIns) <- totalInsertsName

    atacAssays <- c(atacTotalOvs,
                    atacTypeOvs,
                    atacTotalIns,
                    atacTypeIns)
    if(saveHdf5){
      .writeToHdf5(atacAssays, paste0(file.path(outDir, unique(names(d))), ".h5"),
                   storage="integer")
      atacAssays <- lapply(atacAssays, head, n=1)
    }

    return(atacAssays)},
    refCoords, shift, threads, saveHdf5, outDir=outDir,
    BPPARAM=BPPARAM)

  typeCountNames <- names(atacCounts[[1]])

  if(saveHdf5){
    hdf5FilesCollect <- unlist(lapply(colNames, function(f) paste0(file.path(outDir,f), ".h5")))
    atacAssays <- .collectHdf5Files(hdf5FilesCollect,
                                    hdf5FileName=paste0(file.path(outDir, fileName), ".h5"),
                                    storage="integer", asSparse=TRUE)
    atacAssays <- lapply(atacAssays, `colnames<-`, colNames)
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
                         fileName=NULL,
                         outDir=NULL,
                         BPPARAM=SerialParam()){

  threads <- floor(getDTthreads())/BPPARAM$workers

  # reshape list
  colNames <- unique(names(data))
  data <- lapply(colNames, function(name) data[names(data)==name])
  names(data) <- colNames

  chIPPeaks <- BiocParallel::bplapply(data, function(d, refCoords, threads,
                                                     saveHdf5, outDir){

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
      dataList <- list(chIPPeaks)
      names(dataList) <- peakAssayName
      .writeToHdf5(dataList,
                   paste0(file.path(outDir, unique(names(d))), ".h5"),
                   storage="double")
      chIPPeaks <- head(chIPPeaks, n=1)
    }

    return(chIPPeaks)
  }, refCoords=refCoords, threads=threads, saveHdf5=saveHdf5,
     outDir=outDir, BPPARAM=BPPARAM)

  if(saveHdf5){
    hdf5FilesCollect <- unlist(lapply(colNames, function(f) paste0(file.path(outDir,f), ".h5")))
    chIPPeaks <- .collectHdf5Files(hdf5FilesCollect,
                                   hdf5FileName=paste0(file.path(outDir, fileName), ".h5"),
                                   storage="double", asSparse=TRUE)
    chIPPeaks <- chIPPeaks[[1]]
  }
  else{
    chIPPeaks <- Reduce("cbind", chIPPeaks[-1], chIPPeaks[[1]])
  }

  colnames(chIPPeaks) <- colNames
  chIPColData <- data.table(combination=colNames)
  chIPColData[,c(annoCol, tfNameCol):=tstrsplit(combination, split="_")]
  chIPColData[,origin:=lapply(combination, function(x){
    ds <- unlist(data[names(data)==x])
    names(ds) <- unlist(tstrsplit(names(ds), split=".", keep=2, fixed=TRUE))
    ds <- lapply(ds, function(d) fifelse(is.character(d), d, NA))
    return(ds)})]

  assayList <- list(chIPPeaks)
  names(assayList) <- peakAssayName
  chIPSe <- SummarizedExperiment(assays=assayList,
                                 rowRanges=refCoords,
                                 colData=chIPColData)

  return(chIPSe)
}
