#' Predicts transcription-factor binding
#'
#' Predicts bindings with each model in the bag (as obtained by [TFBlearner::trainBagged]) seperately.
#'
#' @name predictTfBinding
#' @param models List of tree-based gradient boosting [lightgbm::lightgbm] models as obtained by [TFBlearner::trainBagged]
#' @param fm [SummarizedExperiment::RangedSummarizedExperiment-class] object containing features as obtained by [TFBlearner::getFeatureMatrix].
#' @param annoCol Name of column indicating cellular contexts.
#' @param simplified If predictions should be returned as a [SummarizedExperiment::RangedSummarizedExperiment-class] object.
#' @param chunk If predictions should be performed on chunks of the data to lower the memory footprint (chunk-size: 1e5)
#' @param sparsify Should predictions be sparsified. Very small binding probabilities will be rounded to zero.
#' @param numThreads Total number of threads to be used. In case [BiocParallel::MulticoreParam] or [BiocParallel::SnowParam] with several workers are
#' are specified as parallel back-ends, `floor(numThreads/nWorker)` threads are used per worker.
#' @param BPPARAM Parallel back-end to be used. Passed to [BiocParallel::bplapply()].
#' @return [SummarizedExperiment::RangedSummarizedExperiment-class] object with predictions of the models saved in the assays (if `simplified=TRUE`)
#' and columns corresponding to cellular contexts.
#' If `simplified=FALSE` a [Matrix::Matrix] is returned with columns corresponding to predictions the models.
#' @import Matrix
#' @importFrom BiocParallel bplapply SerialParam MulticoreParam SnowParam
#' @export
predictTfBinding <- function(models,
                             fm,
                             annoCol="context",
                             simplified=TRUE,
                             chunk=NULL,
                             sparsify=FALSE,
                             numThreads=40,
                             BPPARAM=SerialParam()){

  fmTfName <- metadata(fm)$tf_name
  modTfName <- models[[1]]$params$tf
  if(fmTfName!=modTfName){
    stop(paste("Feature matrix has been computed for", fmTfName, "and model trained for", modTfName))
  }
  featMat <- assays(fm)$features

  data <- featMat
  labelCol <- paste(contextTfFeat, labelFeatName, sep="_")
  contextCol <- annoCol
  cScoreCol <- paste(tfFeat, cScoreFeatName, sep="_")

  npDt <- data.table(ind=1:nrow(data))
  if(!is.null(chunk) & is.numeric(chunk)){
    chunkSize <- min(chunk, nrow(data))
    npDt[,chunk:=cut(ind, ceiling(nrow(npDt)/chunkSize))]
    npDt <- split(npDt, by="chunk")
  }

  if(stackingStratEntry %in% names(models)){
    stackingStrat <- models[[stackingStratEntry]]
    stackedModelName <- paste(modelStackedSuffix, stackingStrat, sep="_")
    modelNamesBag <- setdiff(names(models), c(stackedModelName,
                                              stackingStratEntry))
  }else{
    modelNamesBag <- setdiff(names(models), stackingStratEntry)
  }

  modelsBag <- models[modelNamesBag]
  nWorker <- BPPARAM$workers
  preds <- bpmapply(function(model, name, data,
                             npDt, chunk,
                             labelCol, contextCol,
                             cScoreCol){
    data.table::setDTthreads(floor(numThreads/nWorker))

    allFeats <- listFeatures()
    colsToRemove <- unlist(subset(allFeats,
                                  !included_in_training)$feature_matrix_column_names)
    colsToRemove <- c(colsToRemove, labelCol, contextCol)
    if(name!=modelAllName){
      colsToRemove <- setdiff(colsToRemove, cScoreCol)
    }

    if(!is.null(chunk) & is.numeric(chunk)){
    predDts <- lapply(npDt, function(indDt){
      pred  <- predict(model, as.matrix(data[indDt$ind, !(colnames(data) %in% colsToRemove)]))
      predDt <- data.table(pred=pred)

      if(sparsify) predDt[,pred:=fifelse(pred*scalFactPred>model$params$sparse_thr,pred,0L)]
      predDt
    })
      predDt <- rbindlist(predDts)
      predMat <- Matrix::Matrix(as.matrix(predDt))
    }
    else{
      preds  <- predict(model, as.matrix(data[,!(colnames(data) %in% colsToRemove)]))
      predDt <- data.table(pred=preds)
      if(sparsify) predDt[,pred:=fifelse(pred*scalFactPred>model$params$sparse_thr,pred,0L)]
      predMat <- Matrix::Matrix(as.matrix(predDt))
    }
    gc()
    predMat
  },
  modelsBag,
  names(modelsBag),
  MoreArgs=list(
    data=data,
    npDt=npDt,
    chunk=chunk,
    labelCol=labelCol,
    contextCol=contextCol,
    cScoreCol=cScoreCol),
  BPPARAM=BPPARAM)

  preds <- Reduce(cbind, preds[-1], preds[[1]])
  colnames(preds) <- paste(predPrefix, names(modelsBag), sep="_")

  if(labelCol %in% colnames(data)){
    message("Adding labels")
    label <- data[,labelCol]
    labelBin <- fifelse(label>0,1L,0L)
    labelBin <- fifelse(label<0,-1L,labelBin)
    labelBin <- Matrix(labelBin, ncol=1)
    colnames(labelBin) <- binLabelName

    outPreds <- list(preds, labelBin, data[,labelCol], data[,contextCol])
    cnPreds <- c(colnames(preds), binLabelName, labelCol, contextCol)
  }
  else{
    outPreds <- list(preds, data[,contextCol])
    cnPreds <- c(colnames(preds), contextCol)
  }

  preds <- Matrix::Matrix(Reduce(cbind, outPreds[-1], outPreds[[1]]))
  colnames(preds) <- cnPreds

  if(stackingStratEntry %in% names(models)){
    modelStacked <- models[[stackedModelName]]
    predsStacked <- .predictTfBindingStacked(models, fm,
                                             predsBagged=preds,
                                             annoCol=annoCol)
    predsStackedCol <- paste(predPrefix, stackedSuffix, sep="_")
    preds <- cbind(preds, predsStacked[,predsStackedCol, drop=FALSE])
  }

  if(simplified){
    contextNames <- metadata(fm)[[annoCol]]
    preds <- lapply(setdiff(cnPreds, annoCol), function(col){
      predCol <- as.data.table(as.matrix(preds[,c(col, annoCol),drop=FALSE]))
      predCol <- split(predCol, by=annoCol)
      predCol <- lapply(predCol, function(p) as.numeric(p[[col]]))
      predMat <- Matrix::Matrix(unlist(predCol), nrow=length(predCol[[1]]))
      colnames(predMat) <- contextNames
      predMat
    })
    coords <- rowRanges(fm)[rowRanges(fm)@elementMetadata[[annoCol]]==1]
    names(preds) <- setdiff(cnPreds, annoCol)
    preds <- SummarizedExperiment(assays=preds,
                                  rowRanges=coords)
  }

  return(preds)
}

#' Predicts transcription-factor binding
#'
#' Predicts bindings for the stacked model as obtained by [TFBlearner::trainStacked]
#'
#' @name .predictTfBindingStacked
#' @param modelStacked Stacked model as obtained by [TFBlearner::trainStacked].
#' @param fm [SummarizedExperiment::RangedSummarizedExperiment-class] object containing features as obtained by [TFBlearner::getFeatureMatrix].
#' @param predsBagged Predictions of the single models in the bag. Only needed if the stacked model has been obtained with the weighted mean strategy.
#' @param annoCol Name of column indicating cellular contexts.
#' @param ... Additional arguments passed to [TFBlearner::predictTfBinding].
#' @return Matrix with predicted binding probabilities.
#' @import data.table
#' @import Matrix
.predictTfBindingStacked <- function(models, fm,
                                     predsBagged=NULL, annoCol=NULL, ...){

  tfName <- models[[1]]$params$tf
  stackingStrat <- models[[stackingStratEntry]]
  labelCol <- paste(contextTfFeat, labelFeatName, sep="_")
  cScoreCol <- paste(tfFeat, cScoreFeatName, sep="_")
  predCol <- paste(predPrefix, stackedSuffix, sep="_")
  stackedModelName <- paste(modelStackedSuffix, stackingStrat, sep="_")
  modelStacked <- models[[stackedModelName]]

  if(stackingStrat %in% c("last", "wLast")){
    models <- list(modelStacked)
    names(models) <- stackingStrat
    preds <- predictTfBinding(models, fm, simplified=FALSE, ...)
    preds <- preds[,paste(predPrefix, stackingStrat, sep="_"), drop=FALSE]
  }
  else if(stackingStrat=="wMean"){
    if(is.null(predsBagged)) stop("Provide bagged predictions")
    modelNamesBag <- setdiff(names(models), c(stackedModelName,
                                              stackingStratEntry))
    preds <- lapply(modelNamesBag, function(modelName){
      modelWeight <- models[[modelName]]$params$stacking_weights
      predsBagged[,paste(predPrefix, modelName, sep="_"),drop=FALSE]*modelWeight
    })
    preds <- Matrix(rowSums(Reduce("cbind", preds[-1], preds[[1]])), ncol=1)
  }
  else if(stackingStrat=="boostTree"){
    if(is.null(predsBagged)) stop("Provide bagged predictions")

    # get feature for stacked model
    colSel <- c(paste(siteFeat, widthFeatName, sep="_"),
                annoCol, totalOverlapsName,
                paste(siteFeat, gcContFeatName),
                labelCol,
                paste(assocMotifPrefix, tfName, sep="_"))
    colSel <- intersect(colSel,colnames(fm))
    assaysPred <- lapply(colSel,function(col){
      assays(fm)$features[,col,drop=FALSE]})

    predBagCols <- paste(predPrefix, c(modelTopWeightName, modelMedWeightName,
                         modelAllWeigthName, modelAllName), sep="_")
    predsBagged <- lapply(predBagCols, function(col){
      predsBagged[,col,drop=FALSE]})
    assaysPred <- append(assaysPred, predsBagged)
    assaysPred <- list(Reduce("cbind",assaysPred[-1], assaysPred[[1]]))
    names(assaysPred) <- "features"
    fmStack <- SummarizedExperiment(assays=assaysPred,
                                    rowRanges=rowRanges(fm))
    metadata(fmStack)$tf_name <- tfName

    models <- list(modelStacked)
    names(models) <- stackingStrat
    preds <- predictTfBinding(models, fmStack, simplified=FALSE, ...)
    preds <- preds[,paste(predPrefix, stackingStrat, sep="_"), drop=FALSE]
  }
  colnames(preds) <- predCol

  return(preds)
}
