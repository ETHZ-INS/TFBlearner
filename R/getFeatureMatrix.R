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
   qs <- apply(mat,2, quantile, c(0.25,0.5,0.75))
   Matrix::Matrix(t(t(sweep(mat, 2, qs[2,], "-"))/(qs[3,]-qs[1,])))
}

.minMaxNormalization <- function(mat, useMax=FALSE){
  if(useMax){
    qs <- apply(mat,2, quantile, c(0.0,1.0))
  }
  else{
    qs <- apply(mat,2, quantile, c(0.0,0.9))
  }
  Matrix::Matrix(t(t(mat)/(qs[2,]-qs[1,])))
}

.contextNormalization <- function(mat, method=c("robust", "min-max",
                                                "column", "none")){

  method <- match.arg(method, choices=c("robust", "min-max",
                                        "column", "none"))
  if(method=="column"){
    normMat <- Matrix::t(Matrix::t(mat)/colSums(mat))
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

#' Feature matrix construction
#'
#' Compiles feature matrix based on pre-computed features stored in experiments of the provided [MultiAssayExperiment::MultiAssayExperiment-class] object.
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
#' @return Feature Matrix either as [Matrix::Matrix-class] or as [HDF5Array::HDF5Array-class] (if `saveHdf5=TRUE`).
#' @import Matrix
#' @importFrom rhdf5 h5createFile h5createDataset h5delete h5write H5close H5garbage_collect
#' @importClassesFrom HDF5Array HDF5Array
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

  .checkObject(mae, checkFor=c("Site", "TF", "Context"))

  norm <- match.arg(norm, choices=c("robust", "min-max",
                                    "column", "none"))


  # retrieve specified contexts
  whichCol <- match.arg(whichCol, choices=c("All", "OnlyTrain", "Col"))

  if(whichCol=="OnlyTrain"){
    colSel <- subset(sampleMap(mae), is_training & assay=="ATAC")$primary
  }
  else if(whichCol=="All"){
    colSel <- subset(sampleMap(mae), assay=="ATAC")$primary
  }

  colSel <- unique(colSel)
  whichContexts <- fifelse(addLabels, "Both", "ATAC")
  contexts <- getContexts(mae, tfName, which=whichContexts)
  contexts <- intersect(contexts, colSel)

  featContexts <- unlist(tstrsplit(colnames(mae[["contextTfFeat"]]),
                                   split="_", keep=1))
  if(!all(contexts %in% featContexts)){
    missing <- setdiff(contexts, featContexts)
    warning(paste("Not all the cellular-contexts requested have contextTfFeats compiled.\n
                   Missing are:", paste(missing, collapse=",")))
    contexts <- intersect(contexts, featContexts)
  }

  # get the cofactors
  tfCofactors <- unique(unlist(subset(colData(mae[["tfFeat"]]),
                                      tf_name==tfName)$tf_cofactors))

  message("Attaching Coordinate Features")
  siteFeatMat <- Reduce("cbind", assays(mae[["siteFeat"]])[-1],
                                 assays(mae[["siteFeat"]])[[1]])
  colnames(siteFeatMat) <- names(assays(mae[["siteFeat"]]))

  message("Attaching TF Features")
  seTf <- mae[["tfFeat"]]
  seTf <- seTf[,colData(seTf)$tf_name==tfName]

  tfFeatMat <- Reduce("cbind", assays(seTf)[-1], assays(seTf)[[1]])
  colnames(tfFeatMat) <- names(assays(seTf))

  selMotifs <- subset(colData(mae[["tfFeat"]]),
                      tf_name==tfName)$preselected_motifs
  selMotifs <- unique(unlist(selMotifs))
  motifMat <- as(assays(mae[["Motifs"]])$match_scores[,selMotifs,drop=FALSE],
                 "CsparseMatrix")
  colnames(motifMat) <- paste("motif", colnames(motifMat), sep="_")

  nonContextTfFeat<- list(siteFeatMat, tfFeatMat, motifMat)
  nonContextTfFeat<- Reduce("cbind", nonContextTfFeat[-1], nonContextTfFeat[[1]])

  message("Attaching cellular context-specific features")
  seTfContext <- mae[["contextTfFeat"]][, paste(contexts, tfName, sep="_")]
  seAtac <- mae[["ATAC"]][,contexts]

  # get the number of features
  contextAssayNames <- c(names(assays(mae[["contextTfFeat"]])),
                         names(assays(mae[["ATAC"]])))
  nFeats <- length(assays(mae[["siteFeat"]]))+
            length(assays(mae[["tfFeat"]]))+
            length(assays(mae[["contextTfFeat"]]))+
            length(assays(mae[["ATAC"]]))+
            length(selMotifs)

  if("contextTfFeat_Max_ATAC_Signal" %in% contextAssayNames){
   nFeats <- nFeats+sum(grepl("insert_", contextAssayNames))+
                    sum("total_overlaps" %in% contextAssayNames)
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
                    storage.mode=store, # or integer ??
                    level=3, chunk=c(min(1e6, nrow(seAtac)),nFeats))
  }
  else{
    saveChunk <- FALSE
    hdf5FileName <- NULL
  }

  if(convertInteger){
    nonContextTfFeat<- .roundingCompression(nonContextTfFeat)}

  featMats <- lapply(contexts, function(context, seAtac,
                                        seTfContext,
                                        tfName, tfCofactors,
                                        selMotifs,
                                        otherFeatMat,
                                        norm, saveChunk,
                                        hdf5FileName,
                                        annoCol,
                                        addLabels,
                                        convertInteger){

    # get context- & TF-specific features
    assayNames <- names(assays(mae[["contextTfFeat"]]))
    if(!addLabels) assayNames <- setdiff(assayNames, "contextTfFeat_label")
    featsTfContext <- lapply(assayNames, function(assayName){
                        assayMat <- assays(mae[["contextTfFeat"]])[[assayName]]
                        assayMat[,paste(context, tfName, sep="_"),drop=FALSE]})
    names(featsTfContext) <- featNames <- assayNames
    featsTfContext <- Reduce("cbind", featsTfContext[-1], featsTfContext[[1]])
    colnames(featsTfContext) <- featNames

    # get context-specific features
    featsContext <- lapply(assays(seAtac), function(assayMat){
                           as(assayMat[,context, drop=FALSE], "CsparseMatrix")})
    featsContext <- Reduce("cbind", featsContext[-1], featsContext[[1]])
    colnames(featsContext) <- names(assays(seAtac))

    featsContextMat <- cbind(featsTfContext, featsContext)
    if(addLabels){
      labelCol <- featsContextMat[,"contextTfFeat_label",drop=FALSE]
      featsContextMat <- featsContextMat[,setdiff(colnames(featsContextMat),
                                                  "contextTfFeat_label")]
    }
    else{
      labelCol <- NULL
    }

    # determine sub-mat to be normalized
    featsNormed <- unlist(subset(listFeatures(),
                                 feature_type=="context-tf-Feature" &
                                 context_normed)$feature_matrix_column_names)
    featsNormed <- lapply(c(tfName, tfCofactors, selMotifs), gsub,
                          pattern="<tf_name>|<cofactor_name>", featsNormed)
    featsNormed <- unique(unlist(featsNormed))
    featsNormed <- intersect(colnames(featsContextMat),
                             c(featsNormed, "norm_total_overlaps"))
    featsNormedMat <- featsContextMat[,featsNormed, drop=FALSE]
    featsContextMat <- featsContextMat[, setdiff(colnames(featsContextMat),
                                                 colnames(featsNormedMat))]

    # normalize by maximum ATAC-signal
    if("contextTfFeat_Max_ATAC_Signal" %in% featsNormed){
    countCols <- c(colnames(featsContextMat)[grepl("insert_",
                                                 colnames(featsContextMat))],
                   "total_overlaps")
    scaledSig <- .minMaxNormalization(featsContextMat[,countCols, drop=FALSE])
    maxSig <- featsNormedMat[,"contextTfFeat_Max_ATAC_Signal", drop=TRUE]
    maxScaledMat <- scaledSig / pmax(maxSig, 1e-4)
    colnames(maxScaledMat) <- paste("maxATACscaled",
                                    colnames(maxScaledMat), sep="_")
    }
    else{
      maxScaledMat <- NULL
    }

    # normalize context-specific features
    featsContextMat <- .contextNormalization(featsContextMat, method=norm)

    featsContextMat <- cbind(featsContextMat, featsNormedMat)
    featsContextMat <- cbind(featsContextMat, maxScaledMat)
    featsContextMat <- cbind(featsContextMat, labelCol)

    if(convertInteger){
      featsContextMat <- .roundingCompression(featsContextMat)
    }

    featsMat <- cbind(featsContextMat, otherFeatMat)

    i <- which(colnames(seAtac)==context)
    contextCol <-  Matrix(i, nrow=nrow(featsMat), ncol=1)
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
     selMotifs,
     nonContextTfFeat, norm,
     saveChunk, hdf5FileName,
     annoCol,
     addLabels, convertInteger)

  featMats <- Reduce("rbind", featMats[-1], featMats[[1]])
  featMats <- Matrix::Matrix(featMats)
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
  }

  # add attributes of feature matrix
  attr(featMats, "transcription_factor") <- tfName
  attr(featMats, "cellular_contexts") <- contexts
  attr(featMats, "cofactors") <- tfCofactors

  return(featMats)
}
