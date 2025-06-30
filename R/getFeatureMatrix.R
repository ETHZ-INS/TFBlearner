.roundingCompression <- function(mat, factor=NULL){
  colNames <- colnames(mat)
  if(is.null(factor)){
    maxVal <- max(mat, na.rm=TRUE)
    factor <- floor(.Machine$integer.max / maxVal)
    message(paste("Using", factor, "as a factor for integer conversion"))
  }

  mat <- Matrix::Matrix(as.integer(ceiling(factor*mat)),
                        ncol=ncol(mat),
                        nrow=nrow(mat))
  colnames(mat) <- colNames
  return(mat)
}

.robustNormalization <- function(mat){
   qs <- .marginQuant(mat, probs=c(0.25,0.5,0.75), margin="col")
   Matrix::Matrix(t(t(sweep(mat, 2, qs[2,], "-"))/max((qs[3,]-qs[1,]),1e-5)))
}

.minMaxNormalization <- function(mat, useMax=FALSE){
  if(useMax){
    qs <- .marginQuant(mat, probs=c(0.0,1.0), margin="col")
  }
  else{
    qs <- .marginQuant(mat, probs=c(0.0,0.9), margin="col")
  }
  Matrix::Matrix(t(t(mat)/max((qs[2,]-qs[1,]),1e-5)))
}

.contextNormalization <- function(mat, method=c("robust", "min-max",
                                                "column", "none")){

  method <- match.arg(method, choices=c("robust", "min-max",
                                        "column", "none"))
  if(method=="column"){
    normMat <- Matrix::t(Matrix::t(mat)/pmax(colSums(mat), rep(1e-5, nrow(mat))))
  }
  else if(method=="min-max"){
    normMat <- .minMaxNormalization(mat)
  }
  else if(method=="robust"){
    normMat <- .robustNormalization(mat)
  }
  else if(method=="none")
  {
    normMat <- mat
  }

  return(normMat)
}

.cbindFeat <- function(mats, isDelayed=FALSE){
  colNames <- mapply(function(mat,name){
    if(ncol(mat)>1){
      return(colnames(mat))
    }else{
      return(name)}}, mats, names(mats))
  colNames <- unlist(colNames)

  if(isDelayed) mats <- lapply(mats, .convertToMatrix)

  # homogenize class
  mats <- lapply(mats, drop0)
  xs <- lapply(mats, function(m){
      n <- ncol(m)*nrow(m)
      (length(m@x)/n)})
  nCols <- unlist(lapply(mats, ncol))
  sp <- weighted.mean(unlist(xs), nCols/sum(nCols))

  mats <- lapply(mats, as, "TsparseMatrix")
  js <- cumsum(nCols)
  indDt <- lapply(1:length(mats), function(i){
      j <- js[[i]]
      indDt <- data.table(x=mats[[i]]@x,
                          i=mats[[i]]@i+1,
                          j=rep(j, length(mats[[i]]@x)))})
  indDt <- rbindlist(indDt)

  # Get matrix
  # 2x faster than Reduce("cbind") or do.call("cbind") on the cost of a slight memory overhead
  mat <- sparseMatrix(i=indDt$i, j=indDt$j, x=indDt$x,
                      dims=c(nrow(mats[[1]]), sum(nCols)))
  if(sp>0.5) mat <- suppressWarnings({Matrix::Matrix(mat)})
  colnames(mat) <- colNames

  return(mat)
}

#' Feature matrix construction
#'
#' Compiles feature matrix based on pre-computed features stored in experiments of the provided [MultiAssayExperiment::MultiAssayExperiment-class] object.
#'
#' Columns in the feature matrix returned are named the following:
#' `<featureType>_<featureName>` depending of the feature followed by (`<motifSuffix>, <i>, <normedSuffix>`).
#'
#' @name getFeatureMatrix
#' @param mae [MultiAssayExperiment::MultiAssayExperiment-class] as construced by [TFBlearner::prepData()] containing Motif, ATAC-, ChIP-seq,
#' site-specific features as obtained by [TFBlearner::siteFeatures()], transcription factor-specific features as obtained by [TFBlearner::tfFeatures()]
#' and transcription-factor and cellular-context specific features as obtained by [TFBlearner::contextTfFeatures()].
#' @param tfName Name of transcription factor to compute features for.
#' @param addLabels Should ChIP-seq peak labels be added to the feature matrix.
#' @param whichCol Should feature matrix be calculated for all cellular contexts (`"All"`), only the training data (`"OnlyTrain"`)
#' or only for some specific cellular contexts (`"Col"`) specified in `colSel`.
#' @param colSel If `whichCol="colSel"`, name of the cellular context to compute the feature matrix for.
#' @param norm Normalization strategy to be used for features of different cellular contexts,
#' can be one of"robust" - substraction of median followed by division by IQR, "min-max" - divison by range between minimum and 0.9-quantile, "column" - division by column sums or "none".
#' Default is "robust".
#' @param convertInteger If feature matrix should be converted to integer (to lower memory footprint).
#' Not recommend to be used, might lower predictive performance due to loss of information.
#' @param saveHdf5 If feature matrix should be saved as HDF5 file.
#' @param outDir Directory to save HDF5 file to.
#' @param prefix Prefix added to filename of feature matrix in case saved as HDF5 file.
#' @param annoCol Name of column indicating cellular contexts in colData.
#' @return [SummarizedExperiment::RangedSummarizedExperiment-class] object with the feature matrix stored in the assays.
#' Assay is either saved as [Matrix::Matrix-class] or as [HDF5Array::HDF5Array-class] (if `saveHdf5=TRUE`).
#' @import Matrix
#' @importFrom rhdf5 h5createFile h5createDataset h5delete h5write H5close H5garbage_collect
#' @importFrom MatrixGenerics colQuantiles rowQuantiles
#' @importFrom IRanges ranges
#' @importClassesFrom HDF5Array HDF5Array
#' @author Emanuel Sonder
#' @export
getFeatureMatrix <- function(mae,
                             tfName,
                             addLabels=TRUE,
                             whichCol=c("All", "OnlyTrain", "Col"),
                             colSel=NULL,
                             norm=c("robust", "min-max",
                                    "column", "none"),
                             convertInteger=FALSE,
                             saveHdf5=TRUE,
                             outDir=NULL,
                             prefix=NULL,
                             annoCol="context"){

  .checkObject(mae, checkFor=c("site", "context", "tf", "tf-context"),
               tfName=tfName)

  norm <- match.arg(norm, choices=c("robust", "min-max",
                                    "column", "none"))

  # retrieve specified contexts
  whichCol <- match.arg(whichCol, choices=c("All", "OnlyTrain", "Col"))

  if(whichCol=="OnlyTrain"){
    colSel <- subset(sampleMap(mae), get(ISTRAINCOL) & assay==ATACEXP)$primary
  }
  else if(whichCol=="All"){
    colSel <- subset(sampleMap(mae), assay==ATACEXP)$primary
  }

  colSel <- unique(colSel)
  whichContexts <- fifelse(addLabels, "Both", "ATAC")
  contexts <- getContexts(mae, tfName, which=whichContexts)
  contexts <- intersect(contexts, colSel)

  featContexts <- unlist(tstrsplit(colnames(mae[[CONTEXTTFFEAT]]),
                                   split="_", keep=1))
  if(!all(contexts %in% featContexts)){
    missing <- setdiff(contexts, featContexts)
    warning(paste("Not all the cellular-contexts requested have contextTfFeats compiled.\n
                   Missing are:", paste(missing, collapse=",")))
    contexts <- intersect(contexts, featContexts)
  }

  # get the cofactors
  tfCofactors <- unique(unlist(subset(colData(mae[[TFFEAT]]),
                                      get(TFNAMECOL)==tfName)[[TFCOFACTORSCOL]]))

  message("Attaching Site & TF-Features")
  selMotifs <- subset(colData(mae[[TFFEAT]]),
                      get(TFNAMECOL)==tfName)[[PRESELMOTIFCOL]]
  selMotifs <- unlist(selMotifs)
  motifMat <- assays(mae[[MOTIFEXP]])[[MATCHASSAY]][,selMotifs,drop=FALSE]
  colnames(motifMat) <- paste(TFFEAT, MOTIFFEATNAME, names(selMotifs), sep="_")

  selActMotifs <- subset(colData(mae[[TFFEAT]]),
                         get(TFNAMECOL)==tfName)[[PRESELACTCOL]]
  selActMotifs <- unlist(selActMotifs)
  actMat <- assays(mae[[ASSOCEXP]])[[ASSOCASSAY]][,selActMotifs,drop=FALSE]
  colnames(actMat) <- paste(TFFEAT, ACTASSOCFEATNAME, names(selActMotifs), sep="_") #TODO: should be saved with actual motif name or name like names(selActMotifs)

  # remove sole NA columns
  seTf <- mae[[TFFEAT]]
  seTf <- seTf[,colData(seTf)[[TFNAMECOL]]==tfName]
  isCovered <- lapply(assays(seTf), function(mat) sum(!is.na(mat))>0)
  isCovered <- unlist(isCovered)
  tfAssaysToKeep <- assays(seTf)[isCovered]

  # add the maxATAC / variance here
  seAtac <- mae[[ATACEXP]][,contexts]
  colsStats <- intersect(c(ATACVARFEATNAME, MAXATACFEATNAME),
                         colnames(rowData(seAtac)))
  atacStats <- Matrix::Matrix(as.matrix(rowData(mae[[ATACEXP]])[,colsStats]))
  colnames(atacStats) <- paste(PANCONTEXTFEAT,
                               colsStats, sep="_")

  nonContextTfFeat <- c(assays(mae[[SITEFEAT]]), tfAssaysToKeep, list(motifMat),
                        list(actMat), list(atacStats))
  nonContextTfFeat <- .cbindFeat(nonContextTfFeat, isDelayed=TRUE)
  names(colnames(nonContextTfFeat)) <- NULL

  message("Attaching cellular context-specific features")
  seTfContext <- mae[[CONTEXTTFFEAT]][, paste(contexts, tfName, sep="_")]
  coords <- rowRanges(seAtac)

  # remove sole NA columns
  isCovered <- lapply(assays(seTfContext), function(mat) sum(!is.na(mat))>0)
  isCovered <- unlist(isCovered)
  assaysToKeep <- assays(seTfContext)[isCovered]
  seTfContext <- SummarizedExperiment(assays=assaysToKeep)

  # get the number of features
  nFeats <- ncol(nonContextTfFeat)+
            length(assays(seTfContext))+
            length(assays(seAtac))+
            length(intersect(colnames(colData(seAtac)),
                   paste(MDSDIMFEATNAME, 1:2, sep="_")))

  if(MAXATACCOLNAME %in% colnames(nonContextTfFeat)){
   nFeats <- nFeats+sum(grepl(INSERTFEATNAME, names(assays(seTfContext))))+
                    sum(TOTALOVERLAPSFEATNAME %in% names(assays(seAtac)))
  }
  if(addLabels) nFeats <- nFeats+1 # for context-label column

  if(saveHdf5)
  {
    saveChunk <- fifelse(length(contexts)>3, TRUE, FALSE)
    if(is.null(outDir)) outDir <- getwd() # .
    if(!is.null(prefix)){
      fileName <- paste(prefix, "feature_matrix", tfName, sep="_")}
    else{
      fileName <- paste("feature_matrix", tfName, sep="_")}

    hdf5FileName <- file.path(outDir, paste0(fileName, ".h5"))

    if(file.exists(hdf5FileName)){
      h5delete(hdf5FileName, name="feature_matrix")
      file.remove(hdf5FileName)
    }

    store <- fifelse(convertInteger, "integer", "double")

    h5createFile(hdf5FileName)
    h5createDataset(hdf5FileName, "feature_matrix",
                    dims=c(nrow(seAtac)*length(contexts),nFeats),
                    storage.mode=store,
                    level=1, chunk=c(min(1e6, nrow(seAtac)), nFeats))
  }
  else{
    saveChunk <- FALSE
    hdf5FileName <- NULL
  }

  if(convertInteger){
    nonContextTfFeat<- .roundingCompression(nonContextTfFeat)}

  gc()
  featMats <- lapply(contexts, function(context, seAtac,
                                        seTfContext,
                                        tfName, tfCofactors,
                                        otherFeatMat,
                                        norm, saveChunk,
                                        hdf5FileName,
                                        annoCol,
                                        addLabels,
                                        convertInteger){

    # get context- & TF-specific features
    assayNames <- names(assays(seTfContext))
    if(!addLabels) assayNames <- setdiff(assayNames, LABELCOLNAME)
    featsTfContext <- lapply(assayNames, function(assayName){
                             assayMat <- assays(seTfContext)[[assayName]]
                             assayMat[,paste(context, tfName, sep="_"),drop=FALSE]})
    names(featsTfContext) <- featNames <- assayNames
    featsTfContext <- .cbindFeat(featsTfContext)

    # get context-specific features
    atacAssays <- lapply(assays(seAtac), function(assayMat) assayMat[,context,drop=FALSE])
    featsContext <- .cbindFeat(atacAssays, isDelayed=TRUE)
    mdsCols <- intersect(paste(MDSDIMFEATNAME, 1:2, sep="_"),
                         colnames(colData(seAtac)))
    mdsContext <- colData(seAtac)[context,mdsCols]
    mdsContext <- mdsContext[rep(1, nrow(featsContext)),]
    featsContext <- cbind(featsContext, Matrix::Matrix(as.matrix(mdsContext)))
    colnames(featsContext) <- paste(CONTEXTFEAT, colnames(featsContext), sep="_")

    featsContextMat <- cbind(featsTfContext, featsContext)

    if(addLabels){
      if(!(LABELCOLNAME %in% colnames(featsContextMat))){
        stop("Labels are not available for this context.
              Make sure ChIP-seq labels exist for this TF and cellular context (`getContexts()`).
              Rerun `contextTfFeatures()` with `addLabels=TRUE`")
      }
      labelCol <- featsContextMat[,LABELCOLNAME,drop=FALSE]
      featsContextMat <- featsContextMat[,setdiff(colnames(featsContextMat),
                                                  LABELCOLNAME)]
    }
    else{
      labelCol <- NULL
    }

    # determine sub-mat to be normalized
    featsNormed <- unlist(subset(listFeatures(),
                                 (feature_type=="context-tf-Feature" |
                                  feature_type=="context-Feature") &
                                 context_normed)$feature_matrix_column_names)
    enumFeats <- unlist(tstrsplit(colnames(featsContextMat), split="_"))
    enumFeats <- as.integer(enumFeats[grepl("^[-+]?[0-9]+$", enumFeats)])
    enumFeats <- enumFeats[is.finite(enumFeats) & !is.na(enumFeats)]
    if(length(enumFeats)>0){maxEnum <- max(max(enumFeats),1)}
    else{maxEnum <- 1}

    featsNormed <- lapply(1:maxEnum, gsub, pattern="<i>", featsNormed)
    featsNormed <- unique(unlist(featsNormed))

    featsNormed <- intersect(colnames(featsContextMat),
                             c(featsNormed, paste(CONTEXTFEAT,
                                                  TOTALOVERLAPSFEATNAME,
                                                  GCNORMEDAFFIX, sep="_")))
    featsNormedMat <- featsContextMat[,featsNormed, drop=FALSE]
    featsContextMat <- featsContextMat[, setdiff(colnames(featsContextMat),
                                                 colnames(featsNormedMat))]

    # normalize by maximum ATAC-signal
    if(MAXATACCOLNAME %in% colnames(nonContextTfFeat)){
      whichCol <- grepl(paste(CONTEXTTFFEAT, INSERTFEATNAME, sep="_"),
                        colnames(featsContextMat))
      countCols <- c(colnames(featsContextMat)[whichCol],
                   paste(CONTEXTFEAT, TOTALOVERLAPSFEATNAME, sep="_"))
      scaledSig <- .minMaxNormalization(featsContextMat[,countCols, drop=FALSE])
      maxSig <- nonContextTfFeat[,MAXATACCOLNAME, drop=TRUE]
      maxScaledMat <- scaledSig / pmax(maxSig, 1e-4)
      colnames(maxScaledMat) <- paste(colnames(maxScaledMat),
                                      NORMEDMAXAFFIX, sep="_")
    }
    else{
      maxScaledMat <- NULL
    }

    # normalize context-specific features
    featsContextMat <- .contextNormalization(featsContextMat, method=norm)
    colnames(featsContextMat) <- paste(colnames(featsContextMat),
                                       NORMEDAFFIX, sep="_")

    featsContextMat <- cbind(featsContextMat, featsNormedMat)
    featsContextMat <- cbind(featsContextMat, maxScaledMat)
    featsContextMat <- cbind(featsContextMat, labelCol)
    gc()

    if(convertInteger){
      featsContextMat <- .roundingCompression(featsContextMat)
    }

    matClass <- fifelse(class(otherFeatMat)=="dgCMatrix", "CsparseMatrix", "dgeMatrix")
    featsContextMat <- as(featsContextMat, matClass)
    featsMat <- cbind(featsContextMat, otherFeatMat)

    i <- which(colnames(seAtac)==context)
    contextCol <-  Matrix::Matrix(i, nrow=nrow(featsMat), ncol=1)
    contextCol <- as(contextCol, matClass)
    colnames(contextCol) <- annoCol
    featsMat <- cbind(featsMat, contextCol)

    if(saveChunk){
      toFile <- H5Fopen(hdf5FileName)
      h5writeDataset(as.matrix(featsMat), h5loc=toFile,
                     name="feature_matrix",
                     index=list(((i-1)*nrow(featsMat)+1):(nrow(featsMat)*i),
                                1:ncol(featsMat)))
      H5Fclose(toFile)
      H5garbage_collect()
      featsMat <- head(featsMat,1)
      gc()
      return(featsMat)
    }
    else{
      return(featsMat)
    }
  }, seAtac, seTfContext,
     tfName, tfCofactors,
     nonContextTfFeat, norm,
     saveChunk, hdf5FileName,
     annoCol,
     addLabels, convertInteger)

  featMats <- Reduce("rbind", featMats[-1], featMats[[1]])
  featMats <- suppressWarnings({Matrix::Matrix(featMats)})
  colnames(featMats) <- make.names(colnames(featMats), unique=TRUE)

  if(saveHdf5){
    featNames <- colnames(featMats)
    if(!saveChunk){
      toFile <- H5Fopen(hdf5FileName)
      h5writeDataset(as.matrix(featMats),
                     h5loc=toFile,
                     name="feature_matrix")
      h5closeAll()
    }
    asSparse<- fifelse(is(featMats, "CsparseMatrix"), TRUE, FALSE)
    featMats <- HDF5Array(hdf5FileName, "feature_matrix", as.sparse=asSparse)
    colnames(featMats) <- featNames
    gc()
  }

  # add colData and metadata
  coords <- rep(coords, length(contexts))
  coords@elementMetadata[[annoCol]] <- factor(contexts[featMats[,annoCol]],
                                              levels=contexts)
  fmSe <- SummarizedExperiment(assays=list("features"=featMats),
                               rowRanges=coords)
  metadata(fmSe)[[TFNAMECOL]] <- tfName
  metadata(fmSe)[[TFCOFACTORSCOL]] <- tfCofactors
  metadata(fmSe)[[annoCol]] <- contexts
  metadata(fmSe)[[PRESELMOTIFCOL]] <- selMotifs
  metadata(fmSe)[[PRESELACTCOL]] <- selActMotifs

  return(fmSe)
}
