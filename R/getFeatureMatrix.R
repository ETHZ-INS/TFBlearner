.roundingCompression <- function(mat, factor=1e7){
  colNames <- colnames(mat)
  mat <- Matrix::Matrix(as.integer(ceiling(factor*mat)),
                        ncol=ncol(mat),
                        nrow=nrow(mat))
  colnames(mat) <- colNames
  return(mat)
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
#' @param colNorm If cellular-context specific features should be column-normalized.
#' @param convertInteger If feature matrix should be converted to integer (to lower memory footprint).
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
                             colNorm=TRUE,
                             convertInteger=FALSE,
                             saveHdf5=TRUE,
                             outDir=NULL,
                             prefix=NULL,
                             annoCol="context"){

  .checkObject(mae, checkFor=c("Site", "TF", "Context"))

  whichCol <- match.arg(whichCol, choices=c("All", "OnlyTrain", "Col"))
  if(whichCol=="OnlyTrain"){
    cols <- lapply(experiments(mae),
                   function(n){colnames(n)[colnames(n) %in% unique(subset(sampleMap(mae),
                                                                          is_training)$colname)]})
    mae <- subsetByColumn(mae, cols)
  }
  else if(whichCol=="Col"){
    if(is.null(colSel)) stop("If features should be computed only for some columns (e.g. cellular contexts,
                              (whichCol=Col), please do provide the names via colSel.")
    mae <- mae[,colSel,]
    contexts <- colSel
  }

  if(whichCol!="Col"){
    whichContexts <- fifelse(addLabels, "Both", "ATAC")
    contexts <- getContexts(mae, tfName, which=whichContexts)
  }


  sampleMapDt <- as.data.table(sampleMap(mae))

  message("Attaching Coordinate Features")
  siteFeatMat <- Reduce("cbind", assays(experiments(mae)$siteFeat)[-1],
                         assays(experiments(mae)$siteFeat)[[1]])
  colnames(siteFeatMat) <- names(assays(experiments(mae)$siteFeat))

  message("Attaching TF Features")
  seTf <- experiments(mae)$tfFeat
  seTf <- seTf[,colData(seTf)$tf_name==tfName]

  tfFeatMat <- Reduce("cbind", assays(seTf)[-1], assays(seTf)[[1]])
  colnames(tfFeatMat) <- names(assays(seTf))

  selMotifs <- subset(colData(experiments(mae)$tfFeat),
                      tf_name==tfName)$preselected_motifs
  selMotifs <- unique(unlist(selMotifs))
  motifMat <- as(assays(experiments(mae)$Motifs)$match_scores[,selMotifs,drop=FALSE], "CsparseMatrix")
  colnames(motifMat) <- paste("motif", colnames(motifMat), sep="_")

  noncontextTfFeat <- list(siteFeatMat, tfFeatMat, motifMat)
  noncontextTfFeat <- Reduce("cbind", noncontextTfFeat[-1], noncontextTfFeat[[1]])

  message("Attaching cellular context-specific features")
  seTfContext <- experiments(mae)$contextTfFeat[,paste(contexts, tfName, sep="_")]
  seAtac <- experiments(mae)$ATAC[,contexts]

  # get the number of features
  nFeats <- length(assays(experiments(mae)$siteFeat))+
    length(assays(experiments(mae)$tfFeat))+
    length(assays(experiments(mae)$contextTfFeat))+
    length(assays(experiments(mae)$ATAC))+
    length(selMotifs)+1 # for contextCol

  if(!addLabels &
     "contextTfFeat_label" %in% names(assays(experiments(mae)$contextTfFeat))){
    nFeats <- nFeats-1}

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
                    level=3, chunk=c(min(1e6, nrow(seAtac)),nFeats))
  }
  else{
    saveChunk <- FALSE
    hdf5FileName <- NULL
  }

  if(convertInteger){
    noncontextTfFeat <- .roundingCompression(noncontextTfFeat, factor=1e5)}

  featMats <- lapply(contexts, function(context, seAtac,
                                        seTfContext, otherFeatMat,
                                        colNorm, saveChunk,
                                        hdf5FileName,
                                        annoCol,
                                        addLabels,
                                        convertInteger){

    # get context & TF-specific features
    featsTfContext <- lapply(assays(seTfContext), function(assayMat){
      assayMat[,paste(context, tfName, sep="_"), drop=FALSE]})
    names(featsTfContext) <- featNames <- names(assays(seTfContext))
    if(!addLabels){
      featsTfContext[["contextTfFeat_label"]] <- NULL
      featNames <- setdiff(featNames, "contextTfFeat_label")
    }

    featsTfContext <- Reduce("cbind", featsTfContext[-1], featsTfContext[[1]])
    colnames(featsTfContext) <- featNames

    # get TF-specific features
    featsContext <- lapply(assays(seAtac), function(assayMat){
      as(assayMat[,context, drop=FALSE], "CsparseMatrix")})
    featsContext <- Reduce("cbind", featsContext[-1], featsContext[[1]])
    colnames(featsContext) <- names(assays(seAtac))

    featsContextMat <- cbind(featsTfContext, featsContext)
    if(addLabels){
      labelCol <- featsContextMat[,"contextTfFeat_label",drop=FALSE]
    }
    else{
      labelCol <- NULL
    }

    featsContextMat <- featsContextMat[,setdiff(colnames(featsContextMat),
                                                "contextTfFeat_label")]

    if(colNorm){
      # column normalization
      featsContextMat <- Matrix::t(Matrix::t(featsContextMat)/colSums(featsContextMat))
    }

    featsContextMat <- cbind(featsContextMat, labelCol)

    if(convertInteger){
      featsContextMat <- .roundingCompression(featsContextMat, factor=1e10)
    }

    featsMat <- cbind(featsContextMat, otherFeatMat)

    i <- which(colnames(seAtac)==context)
    contextCol <-  Matrix(i, nrow=nrow(featsMat), ncol=1)
    colnames(contextCol) <- annoCol
    featsMat <- cbind(featsMat, contextCol)

    if(saveChunk){
      h5write(as.matrix(featsMat), file=hdf5FileName,
              name="feature_matrix",
              createnewfile=FALSE,
              index=list(((i-1)*nrow(featsMat)+1):(nrow(featsMat)*i),
                         1:ncol(featsMat)))
      H5garbage_collect()

      return(head(featsMat,1))
    }
    else{
      return(featsMat)
    }
  }, seAtac, seTfContext, noncontextTfFeat, colNorm, saveChunk, hdf5FileName,
     annoCol,
     addLabels, convertInteger)

  featMats <- Reduce("rbind", featMats[-1], featMats[[1]])

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

  # add attributes of feature matrix
  tfCofactors <- unique(unlist(subset(colData(experiments(mae)$ChIP),
                                      tf_name==tfName)$tf_cofactors))

  attr(featMats, "transcription_factor") <- tfName
  attr(featMats, "cellular_contexts") <- contexts
  attr(featMats, "cofactors") <- tfCofactors

  return(featMats)
}
