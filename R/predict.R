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
#' Not done yet
#'
#' @name predictTfBindingBagged
#' @param models List of tree-based gradient boosting [lightgbm::lightgbm] models as obtained by [TFBlearner::trainBagged]
#' @param featMat Feature matrix as obtained by [TFBlearner::getFeatureMatrix] for prediction.
#' @param isLabelled Does the feature matrix contain a label column
#' @param chunk If predictions should be performed on chunks of the data to lower the memory footprint (chunk-size: 1e5)
#' @param sparsify Should predictions be sparsified. Very small binding probabilities will be rounded to zero.
#' @param stackingStrat Stacking strategy to use, either "none", the last model trained (e.g. the one seeing the most data points),
#' the last models using weights ("wLast"), weighted mean of the models based on performance on hold-out data ("wMean"),
# training on ligthGBM data on the predictions of the single trees and hold-out data.
#' @param dataStack Hold-out data to be used for training stacked models (e.g. if `stackingStrat` is either "wMean" or "boostTree").
#' @param numThreads Total number of threads to be used. In case [BiocParallel::MulticoreParam] or [BiocParallel::SnowParam] with several workers are
#' are specified as parallel back-ends, `floor(numThreads/nWorker)` threads are used per worker.
#' @param BPPARAM Parallel back-end to be used. Passed to [BiocParallel::bplapply()].
#' @import Matrix
#' @importFrom BiocParallel bplapply SerialParam MulticoreParam SnowParam
#' @export
predictTfBindingBagged <- function(models,
                                   featMat,
                                   isLabelled=FALSE,
                                   chunk=FALSE,
                                   sparsify=TRUE,
                                   stackingStrat=c("none",
                                                   "last", "wlast",
                                                   "wMean", "boostTree"),
                                   dataStack=NULL,
                                   numThreads=40,
                                   BPPARAM=SerialParam()){
  data <- featMat
  labelCol <- "contextFeat_label"
  contextCol <- "context"


  stackingStrat <- match.arg(stackingStrat,
                             choices=c("none", "last", "wlast",
                                       "wMean", "boostTree"))

  npDt <- data.table(ind=1:nrow(data))
  npDt[,chunk:=cut(ind, ceiling(nrow(npDt)/1e5))]
  npDt <- split(npDt, by="chunk")
  nWorker <- BPPARAM$workers

  preds <- bplapply(models, function(model, data, npDt){
    setDTthreads(floor(numThreads/nWorker))

    if(chunk){
    predDts <- lapply(npDt, function(indDt){
      pred  <- predict(model, as.matrix(data[indDt$ind, colnames(data)!=labelCol &
                                                        colnames(data)!=contextCol &
                                                        colnames(data)!="width" &
                                                        colnames(data)!="weight" &
                                                        colnames(data)!="n_chIP_peaks" &
                                                        colnames(data)!="density_n_chIP_peaks" &
                                                        colnames(data)!="row_id"]))
      predDt <- data.table(pred=pred)

      if(sparsify) predDt[,pred:=fifelse(pred>model$params$sparse_thr,pred,0L)]
      predDt
    })
      predDt <- rbindlist(predDts)
      predMat <- Matrix(as.matrix(predDt))
    }
    else{
      preds  <- predict(model, as.matrix(data[, colnames(data)!=labelCol &
                                                colnames(data)!=contextCol &
                                                colnames(data)!="width" &
                                                colnames(data)!="weight" &
                                                colnames(data)!="n_chIP_peaks" &
                                                colnames(data)!="density_n_chIP_peaks" &
                                                colnames(data)!="row_id"]))
      predDt <- data.table(pred=preds)
      if(sparsify) predDt[,pred:=fifelse(pred>model$params$sparse_thr,pred,0L)]
      predMat <- Matrix(as.matrix(predDt))
    }
    gc()
    predMat
  },
  data=data,
  npDt=npDt,
  BPPARAM=BPPARAM)

  preds <- Reduce(cbind, preds[-1], preds[[1]])
  colnames(preds) <- names(models)

  if(isLabelled){
    label <- data[,labelCol]
    labelBin <- fifelse(label>0,1L,0L)
    labelBin <- fifelse(label<0,-1L,labelBin)

    outPreds <- list(preds, labelBin, data[,contextCol], data[,"n_chIP_peaks"])
    cnPreds <- c(colnames(preds), "label", "context", "n_chIP_peaks")
  }
  else{
    outPreds <- list(preds, data[,contextCol], data[,"n_chIP_peaks"])
    cnPreds <- c(colnames(preds), "context", "n_chIP_peaks")
  }

  preds <- Reduce(cbind, outPreds[-1], outPreds[[1]])
  colnames(preds) <- cnPreds

  # stack
  if(stackingStrat!="none" & isLabelled & !is.null(dataStack)){
    # predict on stacked data
    predsStack <- predictTfBindingBagged(models, dataStack,
                                         labelCol=labelCol,
                                         contextCol=contextCol,
                                         isLabelled=isLabelled,
                                         chunk=chunk,
                                         sparsify=sparsify,
                                         stackingStrat="none",
                                         dataStack=NULL,
                                         BPPARAM=BPPARAM)
    tfName <- models[[1]]$params$tf
    preds <- .predictStacked(tfName, preds, predsStack, dataStack,
                             numThreads=numThreads)
  }
  else if(stackingStrat!="none"){
    stop("Please provide data for training a stacked model")
  }

  return(preds)
}

.predictStacked <- function(tfName,
                            predsPred,
                            dataPred,
                            predsStack,
                            dataStack,
                            numThreads=numThreads,
                            labelCol="chIP_label",
                            stackingStrat=c("last", "wlast",
                                             "wMean", "boostTree")){
  # implement stacking strategies in function!!!
  stackingStrat <- match.arg(stackingStrat,
                             choices=c("last", "wlast", "wMean", "boostTree"))

  if(stackingStrat=="last"){
    predStack <- predsPred[,3]
  }
  else if(stackingStrat=="wLast"){
    predStack <- predsPred[,4]
  }
  else if(stackingStrat=="wMean"){
    # determine weights on leave-out stacking data
    predsStack <- as.data.table(as.matrix(predsStack))
    predsStack <- subset(predsStack, label_bin>=0)
    predsStackLong <- melt(predsStack,
                           id.vars=c("label_bin", "state", "c_score"))

    predsPred$row_id <- 1:nrow(predsPred)
    predsPredLong <- melt(predsPred,
                          id.vars=c("row_id","label_bin", "state", "c_score"))

    # get weights based on performance on hold-out data
    auprDt <- getRocs(predsStackLong,
                      labels="label_bin", scores="value",
                      models="variable", posClass=1, negClass=0)
    auprDt <- auprDt[,.(auc=first(auc_pr_mod)), by=variable]
    auprDt[,w:=auc/sum(auc)]

    predsPredLong <- merge(predsPredLong,
                           auprDt[,w:=auc/sum(auc)],
                           by.x=c("variable"), by.y=c("variable"))

    predsPredLong <- predsPredLong[,.(w_mean_pred=weighted.mean(value, w),
                                      state=first(state),
                                      row_id=first(row_id),
                                      c_score=first(c_score),
                                      label=first(label_bin)), by=row_id]

    # predict on prediction data
    predStack <- predsPredLong$w_mean_pred
    rm(predsPredLong)
  }
  else{
    # boostedTree learning
    predStack <- .trainStackedBoostedTree(tfName=tfName,
                             predsPred=predsPred,
                             dataPred=dataPred,
                             predsStack=predsStack,
                             dataStack=dataStack,
                             labelCol=labelCol,
                             numThreads=numThreads)
  }

  pred <- cbind(predsPred, predStack)

  return(pred)
}

.trainStackedBoostedTree <- function(tfName,
                                     predsPred,
                                     dataPred,
                                     predsStack,
                                     dataStack,
                                     labelCol, numThreads=20){

  # train a model on the stacked data
  colSel <- c("width", "context",
              "norm_counts", "gc_content",
              paste(tfName, "motif", sep="_"))
  fmStack <- dataStack[,colSel]
  fmStack <- Matrix(fmStack, nrow=nrow(fmStack),
                    ncol=ncol(fmStack))
  colnames(fmStack) <- colSel

  # predictions of models in bag
  fmStack <- cbind(fmStack,
                   predsStack[,1:4])
  labelsStack <- predsStack[,"label"]
  nonFlank <- labelsStack>=0 & !is.na(labelsStack)
  fmStack <- fmStack[nonFlank,] # remove flanking regions for training
  labelsBinStack <- labelsStack[nonFlank]

  # select hyperparameters
  set <- 1:nrow(fmStack)
  hp <- .getTunersBagged(set,
                         isWeighted=FALSE,
                         data=as.matrix(fmStack), # is that conversion needed?
                         earlyStoppingRounds=20,
                         evalRounds=100,
                         numThreads=numThreads,
                         measure=msr("classif.aupcr"),
                         contexts=dataStack[nonFlank ,contextCol],
                         labels=labelsBinStack)

  colSel <- c("norm_counts", "gc_content",
              paste(tfName, "motif", sep="_"))
  modStacked <- .fitModel(hp, set,
                          isWeighted=FALSE,
                          featMat=fmStack[,colSel],
                          labels=labelsBinStack,
                          earlyStoppingRounds=20,
                          numThreads=numThreads,
                          posFrac=0.25)

  # get data for prediction
  fmPred <- dataPred[,colSel]
  fmPred <- Matrix(fmPred, nrow=nrow(fmPred),
                           ncol=ncol(fmPred))
  fmPred <- cbind(fmPred,
                  predsPred[,1:4])

  # fit stacked model
   predsStackedMod <- predict(modStacked, fmPred)

   return(predsStackedMod)
}
