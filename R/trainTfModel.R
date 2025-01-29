 MeasureAupr = R6::R6Class("MeasureAucPr",
                          inherit = mlr3::MeasureClassif,
                          public = list(
                            initialize = function(){
                              super$initialize(
                                id = "classif.aucPr",
                                packages = c("data.table", "PRROC"),
                                properties = character(),
                                predict_type = "prob",
                                range = c(0, 1),
                                minimize = FALSE)
                            }
                          ),

                          private = list(
                            # define score as private method
                            .score = function(prediction, ...) {
                              # define loss
                              aupcr=function(prediction){
                                truth1 <- ifelse(prediction$truth==levels(prediction$truth)[1], 1, 0)
                                auc <-  PRROC::pr.curve(scores.class0=prediction$prob[, 1],
                                                        weights.class0=truth1)$auc.integral
                                return(auc)
                              }

                              # call loss function
                              aupcr(prediction)
                            }))

.getTunersBagged <- function(set,
                             isWeighted,
                             data,
                             contexts,
                             measure,
                             earlyStoppingRounds,
                             evalRounds,
                             numThreads,
                             seed=42,
                             tuneHyperparams=TRUE,
                             loContext=FALSE,
                             labels=NULL,
                             weights=NULL){

   set.seed(seed)
   data.table::setDTthreads(numThreads)
   mlr3::set_threads(numThreads)
   subIds <- set

   # all training data
   if(is(data,"sparseMatrix")){
     trainTable <- data.table::as.data.table(as.matrix(data[subIds, ]))}
   else{
     trainTable <- data.table::as.data.table(data[subIds, ])}

   trainTable$label <- labels[subIds]
   if(isWeighted){
    trainTable$weights <- weights[subIds]
   }
   else{
     trainTable$weights <- rep(1,length(subIds))
   }

   trainTable$label <- factor(trainTable$label, levels=c(0,1))

   sumWeightsPos <- sum(subset(trainTable, label==1)$weights)
   sumWeightsNeg <- sum(subset(trainTable, label==0)$weights)
   weightsFracLow <- sumWeightsPos/sumWeightsNeg
   scalePosMin <- (1/weightsFracLow)

   print("Size of training table")
   print(dim(trainTable))

   if(tuneHyperparams){
   task <- mlr3::TaskClassif$new(id="tfb_tree_param",
                                 backend=trainTable,
                                 target="label",
                                 positive="1")

   task$set_col_roles("weights", add_to="weight", remove_from="feature")
   #task$set_col_roles("siteFeat_width", remove_from="feature")

   if(loContext){
     task$set_col_roles("context", add_to="group", remove_from="feature")
     rs <- mlr3::rsmp("loo")
   }
   else{
     task$set_col_roles("context", remove_from="feature")
     rs <- mlr3::rsmp("cv", folds=5)
   }
   rs$instantiate(task)

   learner <- mlr3::lrn("classif.lightgbm",
                        predict_type="prob",
                        objective="binary",
                        verbose=-1,
                        max_bin=63,
                        min_data_in_leaf=20,
                        first_metric_only=TRUE,
                        num_iterations=50,
                        bagging_fraction=0.7,
                        learning_rate=1,
                        seed=seed,
                        num_threads=numThreads)

   if(nrow(trainTable)<775){
    stop("Too few datapoints for training. A minimum of 775 datapoints is required")
   }

   paramSet <- paradox::ParamSet$new(list(
     num_leaves=paradox::p_int(lower=31,upper=min(floor((nrow(trainTable)*0.8)/20), 8192)),
     lambda_l1=paradox::p_dbl(lower=0.001,upper=1000),
     lambda_l2=paradox::p_dbl(lower=0.001,upper=1000),
     feature_fraction=paradox::p_dbl(lower=0.5,upper=1.0),
     bagging_freq=paradox::p_int(lower=1,upper=10),
     scale_pos_weight=paradox::p_dbl(lower=scalePosMin,upper=scalePosMin*1000)))

   terminator <- mlr3tuning::trm("combo",
                     list(mlr3tuning::trm("stagnation", iters=10, threshold=0.001),
                          mlr3tuning::trm("evals", n_evals=evalRounds)), any=FALSE)

   lgbmInst <- mlr3tuning::ti(task=task,
                              learner=learner,
                              resampling=rs,
                              measure=measure,
                              search_space=paramSet,
                              terminator=terminator)

   # Define the tuner using mlr3mbo
   # lgr::get_logger("mlr3")$set_threshold("warn")
   tuner <- mlr3mbo::TunerMbo$new()
   suppressMessages(tuner$optimize(lgbmInst))
   hp <- data.table::as.data.table(lgbmInst$result)}
   else{
     hp <- data.table(num_leaves=min(floor((nrow(trainTable)*0.8)/20), 8192),
                      scale_pos_weight=scalePosMin,
                      lambda_l1=1,
                      lambda_l2=1,
                      feature_fraction=0.75,
                      bagging_freq=5,
                      scale_pos_weight=scalePosMin)
   }

   return(hp)
 }

 .fitModel <- function(hyperparams,
                       set,
                       isWeighted,
                       featMat,
                       labels,
                       numThreads,
                       earlyStoppingRounds,
                       posFrac,
                       weights=NULL,
                       seed=42){
   data.table::setDTthreads(numThreads)
   set.seed(seed)

   if(!isWeighted){
     weights <- rep(1, length(labels))
   }

   posFracMed <- posFrac

   # add validation datapoints keeping similar weight distributions
   valSet <-  sample(set, floor(length(set)*0.15))
   trainSet <- setdiff(set, valSet)

   valSet <- .addUnused(weights, labels, set, valSet,
                        posFrac=posFracMed, nPos=4000,
                        weighted=isWeighted, seed=seed)

   # add training datapoints keeping similar weight distributions
   trainSet <- .addUnused(weights, labels, unique(c(set, valSet)),
                          trainSet, posFrac=0.25, weighted=isWeighted,
                          seed=seed)

   allFeats <- TFBlearner::listFeatures()
   colsToRemove <- unlist(subset(allFeats,
                                 !included_in_training)$feature_matrix_column_names)
   colsSel <- !(colnames(featMat) %in% c("context", "contextFeat_label",
                                         colsToRemove))

   trainData <- lgb.Dataset(data=as.matrix(featMat[trainSet,colsSel,drop=FALSE]),
                            label=labels[trainSet],
                            weight=weights[trainSet],
                            params=list(max_bin=63))
   gc()

   validData1 <- lgb.Dataset(data=as.matrix(featMat[valSet,colsSel,drop=FALSE]),
                             label=labels[valSet],
                             weight=weights[valSet],
                             params=list(max_bin=63))
   gc()

   # get second validation set with observed rate of weights (just recorded)
   topPosFrac <- sum(labels[set]==1)/sum(labels==1)
   nPosSub <- floor((1-topPosFrac)*sum(labels==1)*0.15)

   # sample positives
   indPos <- setdiff(which(labels==1), c(trainSet, valSet))
   indSubPos <- sample(indPos, min(nPosSub, length(indPos)))

   # sample negatives
   indNeg <- setdiff(which(labels==0), trainSet)
   indSubNeg <- sample(indNeg, min(2e5, length(indNeg)))

   indSubPosAll <- c(indSubPos, valSet[labels[valSet]==1])
   nPos <- floor(length(indSubNeg)*((posFracMed)/(1-posFracMed)))
   indSubPos <- sample(indSubPosAll, min(nPos, length(indSubPosAll)))
   if(length(indSubPos)<nPos){
     nNeg <- ((1-posFracMed)/posFracMed)*length(indSubPos)
     indSubNeg <- sample(indNeg, nNeg)
   }

   valSet2 <- c(indSubNeg, indSubPos)
   validData2 <- lgb.Dataset(data=as.matrix(featMat[valSet2,colsSel,drop=FALSE]),
                             label=labels[valSet2],
                             weight=weights[valSet2],
                             params=list(max_bin=63))
   gc()

   # best model is chosen best on the first set
   valids= list("valid_sim_pos"=validData1,
                "valid_all_pos"=validData2)

   # early stopping / best iteration selection metric
   aupr <- function(pred, dval){
     labels <- dval$get_field("label")
     auc <- PRROC::pr.curve(scores.class0=pred,
                            weights.class0=as.integer(labels),
                            curve=FALSE)$auc.integral

     return(list(name="aucpr", value=auc, higher_better=TRUE))
   }

   mod <- suppressMessages(lightgbm(trainData,
                                    params=list(objective="binary",
                                                learning_rate=0.1,
                                                num_leaves=hyperparams$num_leaves,
                                                lambda_l1=hyperparams$lambda_l1,
                                                lambda_l2=hyperparams$lambda_l2,
                                                scale_pos_weight=hyperparams$scale_pos_weight,
                                                feature_fraction=hyperparams$feature_fraction,
                                                bagging_fraction=0.7,
                                                #bagging_fraction=hyperparams$bagging_fraction,
                                                bagging_freq=hyperparams$bagging_freq, # * 10
                                                first_metric_only=TRUE,
                                                metric="custom",#c("average_precision","cross_entropy"),
                                                seed=seed,
                                                num_threads=numThreads),
                                    early_stopping_rounds=earlyStoppingRounds, # early_stopping_rounds !!
                                    #early_stopping_min_delta=
                                    #eval="average_precision",
                                    eval=list("aupr"=aupr),
                                    eval_freq=1,
                                    valids=valids,
                                    verbose=-1,
                                    nrounds=2500))

   # determine sparsification threshold
   predTrain <- predict(mod, as.matrix(featMat[trainSet,colsSel,drop=FALSE]))
   predsTrainDt <- data.table(pred=predTrain*1e3,
                              label=labels[trainSet])
   th <- .findSparsificationThreshold(predsTrainDt, retAll=FALSE)
   rm(predsTrainDt)
   gc()

   mod$params$sparse_thr <- th # sparsification threshold parameter
   mod$params$n_pos_train <- sum(trainData$get_field("label")==1)
   mod$params$n_neg_train <- sum(trainData$get_field("label")==0)
   mod$params$n_pos_val1 <-sum(validData1$get_field("label")==1)
   mod$params$n_pos_val2 <-sum(validData2$get_field("label")==1)
   mod$params$n_neg_val1 <-sum(validData1$get_field("label")==0)
   mod$params$n_neg_val2 <- sum(validData2$get_field("label")==0)
   mod$params$n_overlap_train_val1 <- intersect(trainSet, valSet)
   mod$params$n_overlap_train_val2 <- intersect(trainSet, valSet2)

   return(mod)
 }

 .addUnused <- function(weights,
                        labels,
                        usedSet,
                        resizeSet,
                        posFrac,
                        nPos=NULL,
                        weighted=FALSE,
                        seed){

   set.seed(seed)

   # weights of selected sets
   usedPos <- usedSet[which(labels[usedSet]==1)]
   weightsRefDt <- data.table(bin=cut(weights[usedPos], breaks=10),
                              w=weights[usedPos],
                              ordered=TRUE)
   binRefDt <- weightsRefDt[,.(n=.N,
                               frac=.N/nrow(weightsRefDt),
                               lb=min(w),
                               ub=max(w)), by=bin]
   binRefDt$bin <- factor(binRefDt$bin, ordered=TRUE)
   setkey(binRefDt, lb, ub)

   # weights of unused positives
   notUsed <- setdiff(which(labels==1), usedPos)
   weightsSampDt <- data.table(w_s=weights[notUsed],
                               w_e=weights[notUsed],
                               i=notUsed)
   weightsSampDt <- foverlaps(weightsSampDt,
                              binRefDt,
                              by.x=c("w_s", "w_e"),
                              by.y=c("lb", "ub"),
                              type="within",
                              mult="all",
                              nomatch=NULL)
   weightsSampDt$bin <- factor(weightsSampDt$bin, ordered=TRUE)

   # positives in set to resize
   resizePos <- resizeSet[which(labels[resizeSet]==1)]

   if(nrow(weightsSampDt>0)){
       matchLevels <- head(levels(binRefDt$bin), n=1)==head(levels(weightsSampDt$bin), n=1) &
                      tail(levels(binRefDt$bin), n=1)==tail(levels(weightsSampDt$bin), n=1)}
   else{
       matchLevels <- FALSE
   }

   if(matchLevels){

     if(is.null(nPos))
     {
       nPosAvail <- nrow(weightsSampDt)+length(resizePos)
       nPos <- fifelse(nPosAvail>5e4, 5e4, nPosAvail) # 5e4
     }
     nPosAdd <- nPos-length(resizePos)

     # sample the additional positives
     if(!weighted & nPosAdd>0){
       resizePosAdd <- sample(weightsSampDt$i, nPosAdd)
     }
     else if(nPosAdd>0){
       weightsSampDt[,n_samp:=floor(frac*nPosAdd)]
       weightsSampDt[,n_avail:=.N, by=bin]

       underPopBinsDt <- weightsSampDt[n_samp>n_avail,]
       if(nrow(underPopBinsDt)>0){
         underPopBinsDt[,missing_frac:=n_samp/n_avail]

         underPopBinDt <- subset(underPopBinsDt,
                                 missing_frac==max(underPopBinsDt$missing_frac))
         nMinAvail <- data.table::first(underPopBinDt$n_avail)
         underPopFrac <- data.table::first(underPopBinDt$frac)

         nPosAdd <- floor(nMinAvail/underPopFrac)
         weightsSampDt[,n_samp:=floor(frac*nPosAdd)]
       }

       addPos <- weightsSampDt[,.SD[sample(.N, unique(n_samp))], by=bin]
       resizePosAdd <- addPos$i
     }
     else resizePosAdd <- c()

     # add positives used already
     resizePos <- c(resizePosAdd, resizePos)
     nPosTot <- length(resizePos)
     if(nPosTot>nPos) resizePos <- sample(resizePos, nPos)
   }
   else{
     # get requested number of positives
     if(!is.null(nPos)){
        nPosTot <- length(resizePos)
      if(nPos<nPosTot) resizePos <- sample(resizePos, nPos)
     }
   }

   nPosSamp <- length(resizePos)

   # sample negatives to obtain the desired rate of positives
   nNeg <- ((1-posFrac)/posFrac)*nPosSamp
   resizeNeg <-  resizeSet[which(labels[resizeSet]==0)]
   nNegSamp <- nNeg-length(resizeNeg)

   if(nNegSamp>0){
     allNeg <- setdiff(which(labels==0), c(usedSet, resizeNeg))
     resizeNegAdd <- sample(allNeg, nNegSamp)
     resizeNeg <- c(resizeNeg, resizeNegAdd)
   }
   else{
     resizeNeg <- sample(resizeNeg, nNeg)
   }

   resizeSet <- c(resizePos, resizeNeg)
   return(resizeSet)
 }

 .chooseBags <- function(featMat,
                         weights,
                         labels,
                         nModels,
                         motifName="motif_1",
                         countCol="total_overlaps",
                         posFrac=0.25,
                         seed){

   set.seed(seed)
   nTotPos <- sum(labels)
   posIds <-  which(labels==1)

   if(nModels>1){
    nSampPos <- floor(seq(from=sqrt(min(4000, nTotPos)),
                          to=min(200,sqrt(nTotPos)),
                          length.out=nModels)^2)
    randSamp <- FALSE
   }
   else{
     nSampPos <- min(4e4, nTotPos)
     randSamp <- TRUE
   }

   # sample positive instances -------------------------------------------------
   posSamps <- lapply(1:nModels, function(i){
     if(randSamp)
     {
       ps <- sample(posIds, min(4e4, nTotPos))
     }
     else{
       subWeights <- weights[posIds]
       rand <- sample(1:length(subWeights), length(subWeights))
       samp <- order(subWeights, rand, decreasing=TRUE)[1:nSampPos[i]]
       ps <- posIds[samp]
     }
      ps
   })

   # sample negative instances -------------------------------------------------

   negFact <- (1-posFrac)/posFrac
   negIds <- which(labels==0)
   motifScores <- featMat[,motifName]
   atacFrags <- featMat[,countCol]
   ids <- 1:nrow(featMat)

   .chooseNeg <- function(motifScores, atacFrags, ids,
                          labels, nNeg){

     nEasyNeg <- floor(nNeg*0.3)
     nHardNeg <- nNeg-nEasyNeg

     isNegative <- labels==0
     isEasyNegative <- isNegative &
                       (motifScores < quantile(motifScores, 0.4, na.rm=TRUE) |
                       motifScores < 1 |
                       atacFrags <= quantile(atacFrags, 0.3, na.rm=TRUE))
     isHardNegative <- !isEasyNegative & isNegative

     if(nHardNeg>sum(isHardNegative)){
       diff <- nHardNeg-  sum(isHardNegative)
       nHardNeg <- sum(isHardNegative)
       nEasyNeg <- nEasyNeg+diff
     }

     if(nEasyNeg>sum(isEasyNegative)){
       diff <- nEasyNeg-sum(isEasyNegative)
       nEasyNeg <- sum(isEasyNegative)
       nHardNeg <- nHardNeg+diff
     }

     easyNegIds <- sample(ids[isEasyNegative], nEasyNeg)
     hardNegIds <- sample(ids[isHardNegative], nHardNeg)
     sampNeg <- c(easyNegIds, hardNegIds)

     return(sampNeg)
   }

   negSamps <- lapply(nSampPos, function(p){
     nNeg <- min(p*negFact, length(negIds))

     negSamp <- .chooseNeg(motifScores, atacFrags, ids,
                           labels, nNeg)
     negSamp
   })

   # compose sets for hyperparameter selection
   idSets <- lapply(1:nModels, function(i){
     ids <- c(posSamps[[i]], negSamps[[i]])
   })
   names(idSets) <- paste("nPos", nSampPos, 1:nModels, sep="_")

   return(idSets)
 }

 #' Training transcription factor-specific tree-based gradient boosting Models
 #'
 #' Trains a bag of four tree-based gradient boosting models of the [lightgbm::lightgbm] library.
 #' Hyperparameter selection is performed for each model seperately using model-based optimization by deploying the [mlr3tuning] library.
 #' The lightgbm classification learner used for the hyperparameter selection has been copied from the GitHub repository [https://github.com/mlr-org/mlr3extralearners](https://github.com/mlr-org/mlr3extralearners), whichs
 #' contains the remote package [mlr3extralearners] developed by Raphael Sonabend and Patrick Schratz and Sebastian Fischer.
 #'
 #' @name trainBagged
 #' @param tfName Name of transcription factor to train model for.
 #' @param featMat Feature matrix as constructed by [TFBlearner::getFeatureMatrix].
 #' @param measureName Measure used for hyperparameter selection.
 #' Either area under the precision-recall curve computed using [PRROC::pr.curve] ("classif.aucpr") or
 #' logloss as implemented by [mlr3measures::logloss].
 #' @param evalRounds Number of evaluation rounds for the hyperparameter selection rounds.
 #' @param earlyStoppingRounds Number of early stopping rounds for the hyperparameter selection and training of the [lightgbm::lightgbm] model.
 #' @param posFrac Fraction of positives to use for the training of the model.
 #' Negatives will be supsampled to achieve the specified fraction.
 #' @param loContext Should cellular-contexts be used for leave-one-context-out rounds during hyperparameter selection.
 #' Only works if more than one cellular-context is contained within the feature matrix.
 #' @param tuneHyperparams If hyperparameters should be tuned with [mlr3mbo::TunerMbo]. Recommend to have this turned on (`tuneHyperparams=TRUE`).
 #' Otherwise (hopefully) sensible defaults are used.
 #' @param seed Integer value for setting the seed for random number generation with [base::set.seed].
 #' @param numThreads Total number of threads to be used. In case [BiocParallel::MulticoreParam] or [BiocParallel::SnowParam] with several workers are
 #' are specified as parallel back-ends, `floor(numThreads/nWorker)` threads are used per worker.
 #' @param BPPARAM Parallel back-end to be used. Passed to [BiocParallel::bpmapply()].
 #' @return A list of four [lightgbm::lightgbm] models trained on different strata of the data.
 #' @import mlr3
 #' @import data.table
 #' @import Matrix
 #' @importFrom R6 R6Class
 #' @importFrom mlr3tuning trm ti
 #' @importFrom mlr3mbo TunerMbo
 #' @importFrom mlr3measures logloss
 #' @importFrom mlr3misc crate named_list
 #' @importFrom paradox CondAnyOf CondEqual ParamSet ps p_int p_dbl p_int p_lgl p_uty p_fct
 #' @importFrom BiocParallel bpmapply SerialParam MulticoreParam SnowParam register
 #' @importFrom lightgbm lgb.Dataset lightgbm
 #' @importFrom PRROC pr.curve
 #' @importFrom MatrixGenerics colMaxs
 #' @export
trainBagged <- function(tfName,
                        featMat,
                        measureName=c("classif.aucpr",
                                      "classif.logloss"),
                        evalRounds=100,
                        earlyStoppingRounds=100,
                        posFrac=0.25,
                        loContext=FALSE,
                        tuneHyperparams=TRUE,
                        seed=42,
                        numThreads=10,
                        BPPARAM=SerialParam){
  set.seed(seed)

  fmTfName <- attributes(featMat)$transcription_factor
  if(fmTfName!=tfName){
    stop(paste("Feature matrix has been computed for", fmTfName, "and not for", tfName))
  }

  measureName <- match.arg(measureName, choices=c("classif.aucpr",
                                                  "classif.logloss"))
  if(measureName=="classif.aucpr"){
    mlr3::mlr_measures$add("classif.aucpr", MeasureAupr)}

  labelCol <- "contextTfFeat_label"
  cellTypeCol <- "context"
  countCol <- "total_overlaps"

  measure <- msr(measureName)

  data.table::setDTthreads(numThreads)
  nWorker <- BPPARAM$workers

  motifName <- paste("motif", tfName, sep="_")

  # remove rows with missing values in essential columns -----------------------
  message("Preparing training data")

  featMat <- featMat[!is.na(featMat[,labelCol]) &
                     !is.na(featMat[,countCol]) &
                     !is.na(featMat[,motifName]),]

  # remove peak-flanks
  featMat <- featMat[!featMat[,labelCol]<0, ]

  # weighting of sites ---------------------------------------------------------

  contexts <- unique(featMat[,cellTypeCol])
  scalFact <- max(featMat[,labelCol])
  weights <- featMat[,labelCol]/scalFact

  if(!loContext | length(contexts)==1){
    # scale to weight positives of contexts overall the same
    weightsDt <- data.table(w=weights, con=featMat[,cellTypeCol])
    weightsDt[,sum_w:=sum(w), by=.(con)]
    weightsDt[,scaled_w:=(w/sum(w))*1e3, by=.(con)]
    weights <- weightsDt$scaled_w
  }
  else if(loContext)
  {
    weightsDt <- data.table(w=weights, con=featMat[,cellTypeCol])
    weightsDt$ind <- 1:nrow(weightsDt)
    weightsDt <- subset(weightsDt, w>0)
    weightsDt[,n_pos:=.N, by=con]

    # quantile normalize the top positives
    weightsDt[,w_rank:=frank(-w, ties.method="random"), by=con]
    minNPos <- min(weightsDt$n_pos)
    weightsQuantDt <- lapply(unique(weightsDt$con),
                             function(i){
                               subset(weightsDt, con==i & w_rank<=minNPos)})
    weightsQuantDt <- rbindlist(weightsQuantDt)
    quant <- quantile(weightsQuantDt$w,
                      probs=seq(0,1, length.out=nrow(weightsQuantDt)))
    weightsQuantDt[,scaled_w:=rev(quant)[w_rank]]

    # add low (weight) positives
    weightsLowDt <- subset(weightsDt, w_rank>minNPos)
    weightsLowDt[,w_rank:=frank(w, ties.method="random"), by=con]
    weightsLowDt[,n_pos:=.N, by=con]
    quant2 <- quantile(weightsLowDt$w,
                       probs=seq(0,1, length.out=max(weightsLowDt$n_pos)))
    weightsLowDt$scaled_w <- quant2[pmin(weightsLowDt$w_rank, length(quant))]

    # bind all positives together
    weightsDt <- rbind(weightsQuantDt, weightsLowDt, use.names=TRUE)
    weights <- rep(0, length(weights))
    weights[weightsDt$ind] <- weightsDt$scaled_w
    }

  # label and weight the instances
  weights <- fifelse(weights==0,1,weights)
  labelInd <- as(as.matrix(featMat[,labelCol, drop=FALSE]), "TsparseMatrix")
  labels <- rep(0, nrow(featMat))
  labels[labelInd@i+1] <- 1

  allFeats <- TFBlearner::listFeatures()
  colsToRemove <- unlist(subset(allFeats, !included_in_training)$feature_matrix_column_names)
  colsToRemove <- c(colsToRemove, labelCol)
  featMat <- featMat[,!(colnames(featMat) %in% colsToRemove)]

  setsWeighted <- .chooseBags(featMat,
                              weights,
                              labels,
                              nModels=3,
                              countCol=countCol,
                              motifName=motifName,
                              posFrac=posFrac,
                              seed=seed)
  setsUnweighted <- .chooseBags(featMat,
                                rep(1, length(weights)),
                                labels,
                                nModels=1,
                                countCol=countCol,
                                motifName=motifName,
                                posFrac=posFrac,
                                seed=seed)
  sets <- c(setsWeighted, setsUnweighted)
  isWeighted <- c(TRUE,TRUE,TRUE,FALSE)

  # hyperparameter selection
  message("Hyperparameter selection")
  ptm <- proc.time()
  res <- BiocParallel::bpmapply(.getTunersBagged, sets, isWeighted,
                                MoreArgs=list(data=featMat,
                                              labels=labels,
                                              contexts=contexts,
                                              measure=measure,
                                              weights=weights,
                                              earlyStoppingRounds=earlyStoppingRounds,
                                              evalRounds=evalRounds,
                                              loContext=loContext,
                                              tuneHyperparams=tuneHyperparams,
                                              numThreads=floor(numThreads/nWorker),
                                              seed=seed),
                                  BPPARAM=BPPARAM,
                                  SIMPLIFY=FALSE)
  message(paste("Time elapsed for hyperparameter optimization:", round((proc.time()-ptm)[2],1), "\n"))
  gc()

  message("Training Model")
  ptm <- proc.time()
  fits <- BiocParallel::bpmapply(.fitModel,
                                 res,
                                 sets,
                                 isWeighted,
                 MoreArgs=list(featMat=featMat,
                               labels=labels,
                               weights=weights,
                               numThreads=floor(numThreads/nWorker),
                               earlyStoppingRounds=earlyStoppingRounds,
                               posFrac=0.01,
                               seed=seed),
                 SIMPLIFY=FALSE,
                 BPPARAM=BPPARAM)

  message(paste("Time elapsed for training the model:", round((proc.time()-ptm)[2],1), "\n"))

  fits <- lapply(fits, function(mod){mod$params$tf <- tfName
                                     mod})

  names(fits) <- c("top_weighted_pos",
                   "med_weighted_pos",
                   "all_weighted_pos",
                   "all_pos")

  return(fits)
}

.findSparsificationThreshold <- function(x, AUPRC.cost.factor=5, retAll=TRUE,
                                         triedQuantiles=sqrt((1:29)/30)){
  names(th) <- th <- unique(pmax(1L,quantile(x$pred,triedQuantiles)))
  w <- which(x$label>0L) # treated as true
  res <- t(sapply(th, FUN=function(h){
    p <- x$pred
    p[which(p<h)] <- 0L # set to zero values below the tested threshold
    if(!is.null(x$model)){
      # calculate AUC separately for the different sets
      auc <- sapply(split(seq_along(p), x$model), FUN=function(i){
        pr.curve(p[intersect(w,i)], p[setdiff(i,w)], dg.compute=FALSE)$auc.integral
      })
      return(c(thres=h, AUC.min=min(auc), AUC.median=median(auc), AUC.mean=mean(auc),
               AUC=pr.curve(p[w], p[-w], dg.compute=FALSE)$auc.integral,
               nonZero=sum(p>0L)/length(p)))
    }else{
      # single AUC
      return(c(thres=h, AUC=pr.curve(p[w], p[-w], dg.compute=FALSE)$auc.integral,
               nonZero=sum(p>0L)/length(p)))
    }
  }))
  # if there are more than one AUC per threshold, we do the mean of the minimum,
  # median and mean AUC:
  aucs <- res[,grep("AUC",colnames(res)),drop=FALSE]
  res <- as.data.frame(res)
  res$AUC.cost <- AUPRC.cost.factor*Matrix::rowMeans(aucs)-MatrixGenerics::colMaxs(as.matrix(aucs))
  res$nonZero.cost <- res$nonZero-min(res$nonZero)
  res$cost <- res$nonZero.cost - res$AUC.cost
  res <- res[order(res$cost),]
  if(retAll) return(res)
  res$thres[1]
}


#' Training transcription factor-specific tree-based gradient boosting Models
#'
#' Trains a stacked model provided a bag of four tree-based gradient boosting models as obtained by [TFBlearner::trainBagged].
#' For different stacking strategies can be used.
#' The lightgbm classification learner used for the hyperparameter selection has been copied from the GitHub repository [https://github.com/mlr-org/mlr3extralearners](https://github.com/mlr-org/mlr3extralearners), whichs
#' contains the remote package [mlr3extralearners] developed by Raphael Sonabend and Patrick Schratz and Sebastian Fischer.
#'
#' @name trainStacked
#' @param featMat Labelled feature matrix as obtained with [TFBlearner::getFeatureMatrix]. Ideally not used for training the bagged models.
#' @param modsBagged Bag of models trained on different stratas of the data, as obtained by [TFBlearner::trainBagged].
#' @param stackingStrat Stacking strategy to use. `last`, chooses the the model which has been trained on all (most) positives not using
#' observational weights for the ChIP-seq peaks. `wLast` using the last model which has seen most positives and has been trained with observational weights.
#' `wMean` weighted mean all models based on performance on the feature matrix provided.
#' `boostTree` Trains a lightgbm model on the predictions of the models in the bag, together with some additional features (e.g. gc_content, total_overlaps).
#' @param subSample Number of rows of featMat whoich should be used for computing performance estimates. Only used if `stackingStrat="wMean"`-
#' @param evalRounds Number of evaluation rounds for the hyperparameter selection rounds. Only used if `stackingStrat="wBoost"`.
#' @param earlyStoppingRounds Number of early stopping rounds for the hyperparameter selection and training of the [lightgbm::lightgbm] model. Only used if `stackingStrat="wBoost"`.
#' @param seed Integer value for setting the seed for random number generation with [base::set.seed].
#' @param numThreads Total number of threads to be used. In case [BiocParallel::MulticoreParam] or [BiocParallel::SnowParam] with several workers are
#' are specified as parallel back-ends, `floor(numThreads/nWorker)` threads are used per worker.
#' @param BPPARAM Parallel back-end to be used. Passed to [BiocParallel::bpmapply()].
#' @return Stacked model. Depending on the strategy either a [lightgbm:lightgbm] model (`last`, `wLast`, `boostTree`)
#' or a vector with weights for the models in the provided bag (`wMean`).
#' @import mlr3
#' @import data.table
#' @import Matrix
#' @importFrom R6 R6Class
#' @importFrom mlr3tuning trm ti
#' @importFrom mlr3mbo TunerMbo
#' @importFrom mlr3measures logloss
#' @importFrom mlr3misc crate named_list
#' @importFrom paradox CondAnyOf CondEqual ParamSet ps p_int p_dbl p_int p_lgl p_uty p_fct
#' @importFrom BiocParallel bpmapply SerialParam MulticoreParam SnowParam register
#' @importFrom lightgbm lgb.Dataset lightgbm
#' @importFrom PRROC pr.curve
#' @importFrom MatrixGenerics colMaxs
#' @export
trainStacked <- function(featMat, modsBagged,
                         stackingStrat=c("last", "wLast",
                                         "wMean", "boostTree"),
                         subSample=1e5,
                         evalRounds=100,
                         earlyStoppingRounds=10,
                         seed=42,
                         numThreads=10,
                         BPPARAM=SerialParam()){

  stackingStrat <- match.arg(stackingStrat, choices=c("last", "wLast",
                                                      "wMean", "boostTree"))
  preds <- predictTfBindingBagged(modsBagged, featMat, BPPARAM=BPPARAM)

  if(stackingStrat=="last"){
    stackMod <- modsBagged[[3]]
    attr(stackMod, "stacking_strategy") <- "last"
  }
  else if(stackingStrat=="wLast"){
    stackMod <- modsBagged[[4]]
    attr(stackMod, "stacking_strategy") <- "weighted_last"
  }
  else if(stackingStrat=="wMean"){

    stackMod <- .trainWeightedMean(preds, subSample=subSample, seed=seed)
    attr(stackMod, "stacking_strategy") <- "weighted_mean"
  }
  else{
    stackMod <- .trainStackedBoostedTree(featMat, preds, evalRounds=evalRounds,
                             earlyStoppingRounds=earlyStoppingRounds,
                             numThreads=numThreads, seed=seed)
    attr(stackMod, "stacking_strategy") <- "boosted_tree"
  }

  return(stackMod)
}

.trainWeightedMean <- function(preds, subSample=1e5, seed=42){

  set.seed(seed)
  contextCol <- "context"
  labelCol <- "contextTfFeat_label"

  preds <- preds[preds[,"label_bin"]>=0,]
  preds <- preds[,setdiff(colnames(preds), labelCol)]
  subRows <- sample(nrow(preds))
  predsDt <- as.data.table(as.matrix(preds[subRows,]))

  #predsDt <- melt(predsDt, id.vars=c("label_bin", contextCol,
  #                                   "tfFeat_n_chIP_peaks"))

  predsDt$row_id <- 1:nrow(predsDt)
  predsDt <- melt(predsDt, id.vars=c("row_id","label_bin", contextCol,
                                     "tfFeat_n_chIP_peaks"))

  # get weights based on auc-pr
  auprDt <- .getRocs(predsDt, labels="label_bin", scores="value",
                     models="variable", posClass=1, negClass=0, seed=seed)
  auprDt <- auprDt[,.(auc=data.table::first(auc_pr_mod)), by=variable]
  auprDt[,w:=auc/sum(auc)]

  modelWeights <-  as.list(auprDt$w)
  names(modelWeights) <- auprDt$variable

  return(modelWeights)
}

.trainStackedBoostedTree <- function(featMat,
                                     preds,
                                     evalRounds=100,
                                     earlyStoppingRounds=10,
                                     seed=42,
                                     numThreads=10){
  set.seed(seed)
  message("Training stacked model")
  mlr3::mlr_measures$add("classif.aucpr", MeasureAupr)

  labelCol <- "contextTfFeat_label"
  contextCol <- "context"

  colSel <- c("siteFeat_width", contextCol,
              "total_overlaps", "siteFeat_gc_content",
              labelCol,
              paste("motif", tfName, sep="_"))
  fmStack <- featMat[,colSel]

  # predictions of models in bag
  fmStack <- cbind(fmStack,
                   preds[,c("top_weighted_pos", "med_weighted_pos",
                            "all_weighted_pos", "all_pos")])

  labelsStack <- fmStack[,labelCol]
  nonFlank <- labelsStack>=0 & !is.na(labelsStack)
  fmStack <- fmStack[nonFlank,] # remove flanking regions for training
  labelsBinStack <- labelsStack[nonFlank]
  labelsBinStack <- fifelse(labelsBinStack>0,1,0)
  colSel <- setdiff(colSel, labelCol)

  # select hyperparameters
  set <- 1:nrow(fmStack)
  hp <- .getTunersBagged(set,
                         isWeighted=FALSE,
                         data=as.matrix(fmStack[,colSel]), # is that conversion needed?
                         earlyStoppingRounds=earlyStoppingRounds,
                         evalRounds=evalRounds,
                         numThreads=numThreads,
                         measure=msr("classif.aucpr"),
                         contexts=fmStack[,contextCol],
                         labels=labelsBinStack,
                         seed=seed)

  stackedMod <- .fitModel(hp, set,
                          isWeighted=FALSE,
                          featMat=fmStack[,colSel],
                          labels=labelsBinStack,
                          earlyStoppingRounds=earlyStoppingRounds,
                          numThreads=numThreads,
                          posFrac=0.25,
                          seed=seed)

  return(stackedMod)
}
