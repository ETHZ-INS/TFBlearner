.predictTfBinding <- function(model,
                             data){

  pred  <- predict(model, as.matrix(data[, colnames(data)!="label" & colnames(data)!="context"]))
  predDt <- data.table(pred=pred,
                       label=data[,"label"],
                       state=data[,"context"])

  return(predDt)
}

#' Predicts transcription-factor binding
#'
#' Predicts bindings with each model in the bag (as obtained by [TFBlearner::trainBagged]) seperately.
#'
#' @name predictTfBindingBagged
#' @param models List of tree-based gradient boosting [lightgbm::lightgbm] models as obtained by [TFBlearner::trainBagged]
#' @param featMat Feature matrix to predict on, as obtained by [TFBlearner::getFeatureMatrix] for prediction.
#' @param chunk If predictions should be performed on chunks of the data to lower the memory footprint (chunk-size: 1e5)
#' @param sparsify Should predictions be sparsified. Very small binding probabilities will be rounded to zero.
#' @param numThreads Total number of threads to be used. In case [BiocParallel::MulticoreParam] or [BiocParallel::SnowParam] with several workers are
#' are specified as parallel back-ends, `floor(numThreads/nWorker)` threads are used per worker.
#' @param BPPARAM Parallel back-end to be used. Passed to [BiocParallel::bplapply()].
#' @return Matrix with predicted binding probabilities.
#' @import Matrix
#' @importFrom BiocParallel bplapply SerialParam MulticoreParam SnowParam
#' @export
predictTfBindingBagged <- function(models,
                                   featMat,
                                   chunk=FALSE,
                                   sparsify=TRUE,
                                   numThreads=40,
                                   BPPARAM=SerialParam()){

  fmTfName <- attributes(featMat)$transcription_factor
  modTfName <- mod$all_pos$params$tf
  if(fmTfName!=modTfName){
    stop(paste("Feature matrix has been computed for", fmTfName, "and model trained for", modTfName))
  }

  data <- featMat
  labelCol <- "contextTfFeat_label"
  contextCol <- "context"

  npDt <- data.table(ind=1:nrow(data))
  npDt[,chunk:=cut(ind, ceiling(nrow(npDt)/1e5))]
  npDt <- split(npDt, by="chunk")
  nWorker <- BPPARAM$workers

  preds <- bplapply(models, function(model, data, npDt){
    data.table::setDTthreads(floor(numThreads/nWorker))

    allFeats <- listFeatures()
    colsToRemove <- unlist(subset(allFeats,
                                  !included_in_training)$feature_matrix_column_names)
    colsToRemove <- c(colsToRemove, labelCol, contextCol)

    if(chunk){
    predDts <- lapply(npDt, function(indDt){
      pred  <- predict(model, as.matrix(data[indDt$ind, !(colnames(data) %in% colsToRemove)]))
      predDt <- data.table(pred=pred)

      if(sparsify) predDt[,pred:=fifelse(pred>model$params$sparse_thr,pred,0L)]
      predDt
    })
      predDt <- rbindlist(predDts)
      predMat <- Matrix::Matrix(as.matrix(predDt))
    }
    else{
      preds  <- predict(model, as.matrix(data[,!(colnames(data) %in% colsToRemove)]))
      predDt <- data.table(pred=preds)
      if(sparsify) predDt[,pred:=fifelse(pred>model$params$sparse_thr,pred,0L)]
      predMat <- Matrix::Matrix(as.matrix(predDt))
    }
    gc()
    predMat
  },
  data=data,
  npDt=npDt,
  BPPARAM=BPPARAM)

  preds <- Reduce(cbind, preds[-1], preds[[1]])
  colnames(preds) <- names(models)

  if(labelCol %in% colnames(data)){
    message("add labels")
    label <- data[,labelCol]
    labelBin <- fifelse(label>0,1L,0L)
    labelBin <- fifelse(label<0,-1L,labelBin)
    labelBin <- Matrix(labelBin, ncol=1)
    colnames(labelBin) <- "label_bin"

    outPreds <- list(preds, labelBin, data[,labelCol], data[,contextCol], data[,"tfFeat_n_chIP_peaks"])
    cnPreds <- c(colnames(preds), "label_bin", labelCol, contextCol, "tfFeat_n_chIP_peaks")
  }
  else{
    outPreds <- list(preds, data[,contextCol], data[,"tfFeat_n_chIP_peaks"])
    cnPreds <- c(colnames(preds), contextCol, "tfFeat_n_chIP_peaks")
  }

  preds <- Reduce(cbind, outPreds[-1], outPreds[[1]])
  colnames(preds) <- cnPreds

  return(preds)
}

#' Predicts transcription-factor binding
#'
#' Predicts bindings for the stacked model as obtained by [TFBlearner::trainStacked]
#'
#' @name predictTfBindingStacked
#' @param modelStacked Stacked model as obtained by [TFBlearner::trainStacked].
#' @param featMat Feature matrix to predict on, as obtained by [TFBlearner::getFeatureMatrix] for prediction.
#' @param predsBagged Predictions of the single models in the bag. Only needed if the stacked model has been obtained with the weighted mean strategy.
#' @param ... Additional arguments passed to [TFBlearner::predictTfBindingBagged].
#' @return Matrix with predicted binding probabilities.
#' @import data.table
#' @import Matrix
#' @export
predictTfBindingStacked <- function(modelStacked, featMat, predsBagged=NULL, ...){
  stackingStrat <- attributes(modelStacked)$stacking_strategy
  labelCol <- "contextTfFeat_label"

  if(stackingStrat %in% c("last", "weighted_last")){
    preds <- predictTfBindingBagged(list("preds_stacked"=modelStacked), featMat, ...)
  }
  else if(stackingStrat=="weighted_mean"){
    if(is.null(predsBagged)) stop("Provide bagged predictions")
    preds <- lapply(names(modelStacked), function(modelName){
      predsBagged[,modelName,drop=FALSE]*modelStacked[[modelName]]
    })
    preds <- Matrix(rowSums(Reduce("cbind", preds[-1], preds[[1]])), ncol=1)
    colnames(preds) <- "preds_stacked"

    if(labelCol %in% colnames(featMat)){
      message("add labels")
      label <- featMat[,labelCol]
      labelBin <- fifelse(label>0,1L,0L)
      labelBin <- fifelse(label<0,-1L,labelBin)
      labelBin <- Matrix(labelBin, ncol=1)
      colnames(labelBin) <- "label_bin"

      preds <- list(preds, labelBin, featMat[,labelCol],
                    featMat[,contextCol], featMat[,"tfFeat_n_chIP_peaks"])
      cnPreds <- c("preds_stacked", "label_bin", labelCol, contextCol, "tfFeat_n_chIP_peaks")
    }
    else{
      preds <- list(preds, featMat[,contextCol], featMat[,"tfFeat_n_chIP_peaks"])
      cnPreds <- c("preds_stacked", contextCol, "tfFeat_n_chIP_peaks")
    }
    preds <- Reduce(cbind, preds[-1], preds[[1]])
    colnames(preds) <- cnPreds
  }
  else if(stackingStrat=="boosted_tree"){
    if(is.null(predsBagged)) stop("Provide bagged predictions")
    featMat <- cbind(featMat, predsBagged[,c("top_weighted_pos", "med_weighted_pos",
                                             "all_weighted_pos", "all_pos"),
                                          with=FALSE])
    preds <- predictTfBindingBagged(list("preds_stacked"=modelStacked), featMat, ...)
  }

  return(preds)
}
