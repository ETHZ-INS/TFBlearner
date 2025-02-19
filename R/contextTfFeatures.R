.getChromVARScores <- function(atacMat,
                               matchScores,
                               gcContent,
                               subInd=NULL,
                               expectations=NULL,
                               backgroundPeaks=NULL,
                               seed=42,
                               BPPARAM=SerialParam(),
                               ...){

    set.seed(seed)
    register(BPPARAM) # for usage with chromVAR
    threads <- floor(getDTthreads())/BPPARAM$workers
    data.table::setDTthreads(threads)
    nPeaks <- nrow(atacMat)
    atacFullMat <- atacMat

    # Association with motif score
    if(is.null(subInd))
    {
      nonZero <- which(Matrix::rowSums(atacMat)>0)

      nSub <- min(5e5, length(nonZero))
      subInd <- sample(nonZero, nSub)
    }

    atacMat <- atacMat[subInd,,drop=FALSE]
    matchScores <- matchScores[subInd,,drop=FALSE]
    gcContent <- gcContent[subInd]

    addArgs <- list(...)
    addArgs <- addArgs[names(addArgs) %in% c("niterations", "w", "bs")]
    args <- c(list(object=atacMat, bias=gcContent), addArgs)

    if(is.null(expectations) | is.null(backgroundPeaks)){
      # compute chromVAR deviations
      expectations <- chromVAR::computeExpectations(atacMat)
      backgroundPeaks <- do.call(chromVAR::getBackgroundPeaks,args)
    }

    # compute deviatons
    motExp <- (t(matchScores) %*% expectations)@x
    y <- ((t(matchScores) %*% atacMat) - motExp)/motExp
    yp <- sapply(1:ncol(backgroundPeaks), function(i){
        B <- sparseMatrix(i=1:nrow(backgroundPeaks), j=backgroundPeaks[,i], x=1,
                          dims=c(nrow(backgroundPeaks), nrow(backgroundPeaks)))
        yp <- ((t(matchScores) %*% B) %*% atacMat - ((t(matchScores) %*% B) %*% expectations)@x)/motExp
        yp@x
      }, simplify=TRUE)
    scoreMat <- as((y-mean(yp))/sd(yp), "CsparseMatrix")

    # get association between chromVAR-scores and ATAC-signal
    atacFullMat <- as(atacFullMat, "CsparseMatrix")
    assocMat <- .getAssociation(atacFullMat, scoreMat)
    colnames(assocMat) <- paste("ChromVAR_ATAC",
                                gsub('_[0-9]+',"", colnames(assocMat)),
                                rownames(scoreMat), sep="_")

    # reformat Activity-scores
    colNamesScores <- paste(rownames(scoreMat), rep(colnames(scoreMat),
                                                each=nrow(scoreMat)), sep="_")
    colNamesScores <- paste("ChromVAR_Score", colNamesScores, sep="_")
    activityMats <- lapply(colnames(scoreMat), function(col){
      activityMat <- matrix(scoreMat[,col], nrow=nPeaks,
                                            ncol=nrow(scoreMat),
                                byrow=TRUE)
      colnames(activityMat) <- paste("ChromVAR_score", rownames(scoreMat), sep="_")
      activityMat <- cbind(assocMat, activityMat)
      featNames <- colnames(activityMat)
      activityMat <- lapply(colnames(activityMat), function(featCol){
                                              activityMat[,featCol,drop=FALSE]})
      names(activityMat) <- featNames

      activityMat
    })
    names(activityMats) <- colnames(scoreMat)

    return(list("activity_matrix"=activityMats,
                "sub_ind"=subInd,
                "expectations"=expectations,
                "background_peaks"=backgroundPeaks))
}

#' Cellular context and transcription factor-specific features
#'
#' Adds an experiment with features specific for the specified transcription factor and the cellular-contexts it is covered for,
#' such as Tn5 insertion around motif matches and footprint profiles.
#'
#' @name contextTfFeatures
#' @param mae [MultiAssayExperiment::MultiAssayExperiment-class] as construced by [TFBlearner::prepData()] containing Motif, ATAC-, ChIP-seq,
#' site-specific features as obtained by [TFBlearner::siteFeatures()] and transcription factor-specific features as obtained by [TFBlearner::tfFeatures()].
#' @param tfName Name of transcription factor to compute features for.
#' @param addLabels Should ChIP-seq peak labels be added to the features
#' @param whichCol Should features be calculated for all cellular contexts (`"All"`), only the training data (`"OnlyTrain"`)
#' or only for some specific cellular contexts (`"Col"`) specified in `colSel`.
#' @param colSel If `whichCol="colSel"`, name of the cellular context to compute the features for.
#' @param features Names of features to be added. Can be all or some of "Inserts", "Weighted_Inserts", "Cofactor_Inserts".
#' See [TFBlearner::listFeatures] for an overview of the features.
#' @param annoCol Name of column indicating cellular contexts in colData.
#' @param insertionProfile Pre-computed insertion footprint profile for the specified transcription factor.
#' Needs to contain coordinate (chr/seqnames, start, end) columns and weight column (termed "w").
#' @param subSample If "ChromVAR_Scores" amongst features, the number of sites to subsample for computing the scores. If `NULL` no subsampling performed.
#' @param aggregationFun function (e.g. mean, median, sum) used to aggregate features across the rowRanges of experiments of the
#' provided [MultiAssayExperiment::MultiAssayExperiment-class] object.
#' @param seed Integer value for setting the seed for random number generation with [base::set.seed].
#' @param BPPARAM Parallel back-end to be used. Passed to [BiocParallel::bpmapply()] & [BiocParallel::bplapply()].
#' @param ... Arguments passed to [TFBlearner::getInsertionProfiles] and [chromVAR::getBackgroundPeaks].
#' @return [MultiAssayExperiment::MultiAssayExperiment-class] object with an experiment containing transcription factor- and cellular context-specific features added to [MultiAssayExperiment::experiments].
#' If already an experiment of this feature group exists, columns are added to it.
#' @import MultiAssayExperiment
#' @importFrom BiocParallel bpmapply bplapply SerialParam MulticoreParam SnowParam register
#' @importFrom chromVAR getBackgroundPeaks computeExpectations
#' @export
contextTfFeatures <- function(mae,
                              tfName,
                              addLabels=FALSE,
                              whichCol=c("All", "OnlyTrain", "Col"), # evt. rename these arguments
                              colSel=NULL,
                              features=c("Inserts", "Weighted_Inserts",
                                         "Cofactor_Inserts", "ChromVAR_Scores",
                                         "Cofactor_ChromVAR_Scores"),
                              annoCol="context",
                              insertionProfile=NULL,
                              subSample=1e4,
                              aggregationFun=sum,
                              seed=42,
                              BPPARAM=SerialParam(),
                              ...){

  .checkObject(mae, checkFor=c("Site", "TF"))

  whichCol <- match.arg(whichCol, choices=c("All", "OnlyTrain", "Col"))
  whichContexts <- fifelse(addLabels, "Both", "ATAC")
  if(whichCol=="OnlyTrain"){
    cols <- lapply(experiments(mae),
                   function(n){colnames(n)[colnames(n) %in% unique(subset(sampleMap(mae),
                                                                          is_training)$colname)]})
    maeSub <- subsetByColumn(mae, cols)
    # get all cellular contexts covered for that TF
    contexts <- getContexts(maeSub, tfName, which=whichContexts)
  }
  else if(whichCol=="Col"){
    if(is.null(colSel)) stop("If features should be computed only for some columns (e.g. cellular contexts,
                              (whichCol=Col), please do provide the names via colSel.")
    maeSub <- mae[,colSel,]
    contexts <- colSel
  }
  else{
    maeSub <- mae
    # get all cellular contexts covered for that TF
    contexts <- getContexts(maeSub, tfName, which=whichContexts)
  }

  coords <- rowRanges(experiments(maeSub)$Motifs)

  if(!is.null(insertionProfile)){
    insertionProfile <- .processData(insertionProfile, readAll=TRUE,
                                     shift=FALSE,
                                     seqLevelStyle=seqlevelsStyle(coords))
  }

  features <- match.arg(features, choices=c("Inserts", "Weighted_Inserts",
                                            "Cofactor_Inserts",
                                            "ChromVAR_Scores",
                                            "Cofactor_ChromVAR_Scores"),
                        several.ok=TRUE)

  tfCofactors <- unique(unlist(subset(colData(experiments(maeSub)$tfFeat),
                                      tf_name==tfName)$tf_cofactors))
  if(("Cofactor_Inserts" %in% features |
      "Cofactor_ChromVAR_Scores" %in% features) & is.null(tfCofactors)){
    msg <- c("No cofactors have been specified when computing transcription ",
             "factor-specific features for ", tfName, ". ", "\n",
             "Cofactor_Inserts will not be computed. ", "\n",
             "Re-compute transcription factor-specific features (tfFeatures()) ",
             "for ", tfName, " with argument `tfCofactors` specified ", "\n",
             "if Cofactor_Insert features are required. ")

    msg <- paste0(msg)
    warning(msg)
  }

  atacFragPaths <- unlist(subset(colData(experiments(maeSub)$ATAC),
                                 get(annoCol) %in% contexts)$origin)

  # get list of motif ranges
  motifRanges <-  names(experiments(maeSub)[grepl("match_ranges",
                                                  names(experiments(maeSub)))])
  motifRanges <- lapply(experiments(maeSub)[motifRanges], rowRanges)
  names(motifRanges) <- unlist(tstrsplit(names(motifRanges), split="_", keep=3))
  motifRanges <- motifRanges[c(tfName, tfCofactors)]

  if(addLabels){
    colDataChIP <- colData(experiments(maeSub)$ChIP)
    colDataChIP <- subset(colDataChIP, tf_name==tfName)
    labelCols <- colDataChIP$combination
    names(labelCols) <- colDataChIP[[annoCol]]
    labels <- lapply(labelCols, function(col){
      as(assays(experiments(maeSub)$ChIP)$peaks[,col,drop=TRUE], "CsparseMatrix")})
  }
  else{
    labels <- vector(mode = "list", length = 2)
    names(labels) <- contexts
  }

  # loop over contexts to get the features
  labels <- labels[contexts] # ensure ordering
  threads <- floor(getDTthreads())/BPPARAM$workers
  feats <- BiocParallel::bpmapply(function(context,
                                           labels,
                                           coords,
                                           atacFrag,
                                           motifRanges,
                                           features,
                                           profile,
                                           tfName,
                                           aggregationFun,
                                           threads, ...){
    data.table::setDTthreads(threads)

    calcProfile <- FALSE
    if("Weighted_Inserts" %in% features & is.null(profile)){
      calcProfile <- TRUE
    }

    atacFrag <- atacFrag[names(atacFrag)==context]

    addArgs <- list(...)
    addArgs <- addArgs[names(addArgs) %in% c("margin", "shift", "symmetric", "stranded")]
    args <- c(list(atacData=atacFrag, motifRanges=motifRanges,
                   profiles=profile, calcProfile=calcProfile), addArgs)
    insRes <- suppressMessages(do.call(getInsertionProfiles, args))

    scoreCols <- c()
    if("Weighted_Inserts" %in% features) scoreCols <- c(scoreCols,
                                                        "weighted_insert_counts")
    if("Inserts" %in% features) scoreCols <- c(scoreCols, "insert_counts")

    # aggregate features across ranges of interest
    insFeats <- lapply(scoreCols, function(scoreCol){
      feats <- genomicRangesMapping(coords,
                                    insRes$motifScores,
                                    scoreCol=scoreCol,
                                    byCols="type",
                                    aggregationFun=aggregationFun)
      colnames(feats) <- fifelse(colnames(feats)=="0",
                                 paste(tfName, scoreCol, "margin", sep="_"),
                                 paste(tfName, scoreCol, "within", sep="_"))
      feats
    })
    insFeats <- Reduce("cbind", insFeats[-1], insFeats[[1]])
    namesFeats <- colnames(insFeats)
    insFeats <- lapply(namesFeats, function(col) insFeats[,col,drop=FALSE])
    names(insFeats) <- namesFeats

    if(!is.null(labels)){
      insFeats <- append(insFeats, list("label"=Matrix::Matrix(labels, ncol=1)))
    }

    return(insFeats)
  }, contexts, labels,
     MoreArgs=list(coords=coords, atacFrag=atacFragPaths,
                   motifRanges=motifRanges[[tfName]], features=features,
                   profile=insertionProfile[[tfName]], tfName=tfName,
                   threads=threads, aggregationFun=aggregationFun, ...),
     SIMPLIFY=FALSE,
     BPPARAM=BPPARAM)

  tfCofactorsSub <- intersect(tfCofactors, names(motifRanges))
  if("Cofactor_Inserts" %in% features & length(tfCofactorsSub)>0){

    coFeats <- BiocParallel::bplapply(contexts,
                                      function(context,
                                               coords,
                                               atacFrag,
                                               tfCofactors,
                                               motifRanges,
                                               features,
                                               profile,
                                               aggregationFun,
                                               threads, ...){
      data.table::setDTthreads(threads)

      calcProfile <- fifelse("Weighted_Inserts" %in% features, TRUE, FALSE)
      scoreCols <- c()
      if("Weighted_Inserts" %in% features){
        scoreCols <- c(scoreCols, "weighted_insert_counts")}
      if("Inserts" %in% features) scoreCols <- c(scoreCols, "insert_counts")

      atacFrag <- atacFrag[names(atacFrag)==context]

      insCoFeats <- lapply(1:length(tfCofactors), function(i){
          coFact <- tfCofactors[[i]]
          cofactRanges <- motifRanges[[coFact]]
          profile <- insertionProfile[[coFact]]
          if(!is.null(profile)) calcProfile <- FALSE

          addArgs <- list(...)
          addArgs <- addArgs[names(addArgs) %in% c("margin", "shift", "symmetric", "stranded")]
          args <- c(list(atacData=atacFrag, motifRanges=cofactRanges,
                         profiles=profile, calcProfile=calcProfile), addArgs)
          insResCo <- suppressMessages(do.call(getInsertionProfiles, args))

          insFeats <- lapply(scoreCols, function(scoreCol){
            feats <- genomicRangesMapping(coords,
                                           insResCo$motifScores,
                                           scoreCol=scoreCol,
                                           byCols="type",
                                           aggregationFun=aggregationFun)
            colnames(feats) <- fifelse(colnames(feats)=="0",
                                       paste(coFact, i, scoreCol, "margin", sep="_"),
                                       paste(coFact, i, scoreCol, "within", sep="_"))
            feats
          })
         insFeats <- Reduce("cbind", insFeats[-1], insFeats[[1]])
      })

      insCoFeats <- Reduce("cbind", insCoFeats[-1], insCoFeats[[1]])
      namesFeats <- colnames(insCoFeats)
      insCoFeats <- lapply(namesFeats,
                           function(col) insCoFeats[,col,drop=FALSE])
      names(insCoFeats) <- namesFeats

      return(insCoFeats)
    }, coords=coords, atacFrag=atacFragPaths, tfCofactors=tfCofactorsSub,
       motifRanges=motifRanges[tfCofactorsSub], features=features,
       profile=insertionProfile[tfCofactorsSub],
       aggregationFun=aggregationFun, threads=threads, ...,
       BPPARAM=BPPARAM)
    names(coFeats) <- contexts

    feats <- lapply(contexts, function(context){
      c(coFeats[[context]], feats[[context]])
    })
    names(feats) <- contexts
  }

  if("ChromVAR_Scores" %in% features){
    message("Get chromVAR features")

    if(!is.null(colData(experiments(maeSub)$siteFeat)$ChromVAR_sub_ind) &
       !is.null(colData(experiments(maeSub)$siteFeat)$ChromVAR_expectations) &
       !is.null(colData(experiments(maeSub)$siteFeat)$ChromVAR_background_peaks)){

      subInd <- unlist(colData(experiments(maeSub)$siteFeat)$ChromVAR_sub_ind)
      expectations <- unlist(colData(experiments(maeSub)$siteFeat)$ChromVAR_expectations)
      backgroundPeaks <- colData(experiments(maeSub)$siteFeat)$ChromVAR_background_peaks[[1]]

      atacMat <- as(assays(experiments(maeSub)$ATAC)$total_overlaps[,contexts,drop=FALSE],
                    "CsparseMatrix")
    }
    else{
      if(whichCol=="Col"){
        cols <- lapply(experiments(mae),
                     function(n){colnames(n)[colnames(n) %in% unique(subset(sampleMap(mae),
                                                                            is_training)$colname)]})
        maeTrain <- subsetByColumn(mae, cols)}
      else{
        maeTrain <- maeSub
      }
      atacMat <- as(assays(experiments(maeTrain)$ATAC)$total_overlaps,"CsparseMatrix")
      rm(maeTrain)

      subInd <- expectations <- backgroundPeaks <- NULL
    }


    if("Cofactor_ChromVAR_Scores" %in% features){
      tfCols <- c(tfName, tfCofactors)}
    else{
      tfCols <- tfName}

    cols <- grep(paste(tfCols,collapse="|"),
                 colnames(experiments(maeSub)$Motifs), value=TRUE)
    matchScores <- as(as(assays(experiments(maeSub)$Motif)$match_scores[,cols,drop=FALSE],
                      "CsparseMatrix"), "TsparseMatrix")
    colDataMotifs <- subset(colData(experiments(maeSub)$Motifs), motif %in% cols)
    colDataMotifs <- colDataMotifs[order(match(colDataMotifs$motif, cols)),,
                                   drop=FALSE]
    thr <- colDataMotifs$max_score/2

    matchScores@x[matchScores@x<thr[matchScores@j + 1] & matchScores@x<4e4] <- 0
    gcContent <-  assays(experiments(maeSub)$siteFeat)$siteFeat_gc_content[,,drop=TRUE]
    res <- .getChromVARScores(atacMat, matchScores, gcContent,
                              subInd=subInd, expectations=expectations,
                              backgroundPeaks=backgroundPeaks,
                              seed=seed, BPPARAM=BPPARAM, ...)
    actFeatMats <- res$activity_matrix

    feats <- lapply(contexts, function(context){
      c(feats[[context]], actFeatMats[[context]])
    })
  }

  names(feats) <- contexts

  # add features back to the full object
  for(context in contexts){
    featMats <- lapply(feats[[context]], `colnames<-`, NULL)
    names(featMats) <- paste("contextTfFeat", names(featMats), sep="_")

    seTfFeat <- SummarizedExperiment(assays=featMats,
                                     rowRanges=coords)
    colnames(seTfFeat) <- paste(context, tfName, sep="_")
    colData(seTfFeat)$feature_type <- "contextTfFeat"
    colData(seTfFeat)$tf_name <- tfName
    mae <- .addFeatures(mae, seTfFeat, colsToMap=context,
                        prefix="contextTfFeat")
  }

  if("ChromVAR_Scores" %in% features |
     "Cofactor_ChromVAR_Scores" %in% features){
    # add chromVAR parameters to object
    colData(experiments(mae)$siteFeat)$ChromVAR_sub_ind <- list(res$sub_ind)
    colData(experiments(mae)$siteFeat)$ChromVAR_expectations <- list(res$expectations)
    colData(experiments(mae)$siteFeat)$ChromVAR_background_peaks <- list(res$background_peaks)}

  return(mae)
}
