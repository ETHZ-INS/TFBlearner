.roundingCompression <- function(mat, factor=1e7){
  colNames <- colnames(mat)
  mat <- Matrix::Matrix(as.integer(ceiling(factor*mat)),
                        ncol=ncol(mat),
                        nrow=nrow(mat))
  colnames(mat) <- colNames
  return(mat)
}

getfeatureMatrix <- function(mae,
                             tfName,
                             addLabels=TRUE,
                             whichCol=c("All", "OnlyTrain", "Col"),
                             colNorm=TRUE,
                             saveHdf5=TRUE,
                             convertInteger=FALSE,
                             outDir=NULL,
                             prefix=NULL,
                             annoCol="context",
                             BPPARAM=SerialParam()){

  .checkObject(mae, checkFor=c("Coord", "TF", "Context"))

  whichCol <- match.arg(whichCol, choices=c("All", "OnlyTrain", "Col"))
  if(whichCol=="OnlyTrain"){
    maeSub <- mae[,colData(mae)$is_training]
  }
  else if(whichCol=="Col"){
    if(is.null(colSel)) stop("If features should be computed only for some columns (e.g. cellular contexts,
                              (whichCol=Col), please do provide the names via colSel.")
    maeSub <- mae[,colSel,]
  }
  else{
    maeSub <- mae
  }

  sampleMapDt <- as.data.table(sampleMap(mae))

  message("Attaching Coordinate Features")
  coordFeatMat <- Reduce("cbind", assays(experiments(mae)$coordFeat)[-1],
                         assays(experiments(mae)$coordFeat)[[1]])
  colnames(coordFeatMat) <- names(assays(experiments(mae)$coordFeat))

  message("Attaching TF Features")
  seTf <- experiments(mae)$tfFeat
  seTf <- seTf[,colData(seTf)$tf_name==tfName]

  tfFeatMat <- Reduce("cbind", assays(seTf)[-1], assays(seTf)[[1]])
  colnames(tfFeatMat) <- names(assays(seTf))

  selMotifs <- subset(colData(experiments(mae)$ChIP),
                      tf_name==tfName)$preselected_motifs
  selMotifs <- unique(unlist(selMotifs))
  motifMat <- assays(experiments(mae)$Motifs)$match_scores[,selMotifs, drop=FALSE]
  colnames(motifMat) <- paste("motif", colnames(motifMat), sep="_")

  nonContextFeat <- list(coordFeatMat, tfFeatMat, motifMat)
  nonContextFeat <- Reduce("cbind", nonContextFeat[-1], nonContextFeat[[1]])

  message("Attaching cellular context-specific features")
  contexts <- getContexts(mae, tfName)
  seTfContext <- experiments(mae)$contextFeat[,paste("contextFeat", tfName, contexts, sep="_")]
  seAtac <- experiments(mae)$ATAC[,contexts]

  # get the number of features
  nFeats <- length(assays(experiments(mae)$coordFeat))+
    length(assays(experiments(mae)$tfFeat))+
    length(assays(experiments(mae)$contextFeat))+
    length(assays(experiments(mae)$ATAC))+
    length(selMotifs)+1 # for contextCol
  if(!addLabels) nFeats <- nFeats-1

  #TODO: check better heueristic for this
  if(saveHdf5)
  {
    saveChunk <- fifelse(length(contexts)>3, TRUE, FALSE)
    if(is.null(outDir)) outDir <- getwd() # .
    if(!is.null(prefix)){
      fileName <- paste(prefix, "feature_matrix", tfName, sep="_")}
    else{
      fileName <- paste("feature_matrix", tfName, sep="_")}

    hdf5FileName <- file.path(outDir, fileName)

    if(file.exists(hdf5FileName)){
      h5delete(hdf5FileName, name="feature_matrix")
      file.remove(hdf5FileName)
    }

    store <- fifelse(convertInteger, "integer", "double")

    h5createFile(hdf5FileName)
    h5createDataset(hdf5FileName, "feature_matrix",
                    dims=c(nrow(seAtac)*length(contexts),nFeats),
                    storage.mode=store, # or integer ??
                    chunk=c(1e5,nFeats))
  }
  else{
    saveChunk <- FALSE
    hdf5FileName <- NULL
  }

  if(convertInteger){
    nonContextFeat <- .roundingCompression(nonContextFeat, factor=1e5)}

  featMats <- BiocParallel::bplapply(contexts, function(context, seAtac,
                                                        seTfContext,
                                                        otherFeatMat,
                                                        colNorm,
                                                        saveChunk,
                                                        hdf5FileName,
                                                        addLabels,
                                                        convertInteger){

    # get context & TF-specific features
    featsTfContext <- lapply(assays(seTfContext), function(assayMat){
      assayMat[,paste("contextFeat", tfName, context, sep="_"), drop=FALSE]})
    if(!addLabels) featsTfContext[["contextFeat_label"]] <- NULL

    featsTfContext <- Reduce("cbind", featsTfContext[-1], featsTfContext[[1]])
    colnames(featsTfContext) <- names(assays(seTfContext))

    # get TF-specific features
    featsContext <- lapply(assays(seAtac), function(assayMat){
      assayMat[,context, drop=FALSE]})
    featsContext <- Reduce("cbind", featsContext[-1], featsContext[[1]])
    colnames(featsContext) <- names(assays(seAtac))

    featsContextMat <- cbind(featsTfContext, featsContext)
    labelCol <- featsContextMat[,"contextFeat_label",drop=FALSE]
    featsContextMat <- featsContextMat[,setdiff(colnames(featsContextMat),
                                                "contextFeat_label")]
    if(colNorm){
      # column normalization
      featsContextMat <- Matrix::t(Matrix::t(featsContextMat)/colSums(featsContextMat))
    }

    featsContextMat <- cbind(featsContextMat, labelCol)

    if(convertInteger){
      featsContextMat <- .roundingCompression(featsContextMat, factor=1e10)
    }

    featsMat <- cbind(featsContextMat, otherFeatMat)

    if(saveChunk){
      i <- which(colnames(seAtac)==context)
      h5write(as.matrix(featsMat), file=hdf5FileName,
              name="feature_matrix",
              createnewfile=FALSE,
              index=list(((i-1)*nrow(featsMat)+1):(nrow(featsMat)*i),
                         1:ncol(featsMat)))
      #TODO: Dirty
      return(head(featsMat,1))
    }
    else{
      return(featsMat)
    }
  }, seAtac, seTfContext, nonContextFeat,
  colNorm, saveChunk, hdf5FileName,
  addLabels, convertInteger, BPPARAM=BPPARAM)

  contextFact <- as.integer(factor(contexts, levels=contexts))
  contextCol <- Matrix::Matrix(rep(contextFact, each=nrow(featMats[[1]])), ncol=1)
  colnames(contextCol) <- annoCol
  featMats <- Reduce("rbind", featMats[-1], featMats[[1]])
  featMats <- cbind(featMats, contextCol)
  #featMats <- Matrix(featMats)

  if(saveHdf5){
    featNames <- colnames(featMats)
    if(!saveChunk){
      # save the full chunk
      h5write(as.matrix(featMats),
              file=hdf5FileName,
              name="feature_matrix",
              createnewfile=FALSE)
    }
    H5close()
    featMats <- HDF5Array(hdf5FileName, "feature_matrix", as.sparse=TRUE)
    colnames(featMats) <- featNames
  }

  return(featMats)
}
