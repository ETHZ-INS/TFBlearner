.getChromVARScores <- function(atacMat,
                               matchScores,
                               gcContent,
                               subInd=NULL,
                               expectations=NULL,
                               backgroundPeaks=NULL,
                               assocMat=NULL,
                               seed=42,
                               BPPARAM=SerialParam(),
                               ...){

  set.seed(seed)
  register(SerialParam()) # for usage with chromVAR
  threads <- floor(getDTthreads())/BPPARAM$workers
  data.table::setDTthreads(threads)
  nPeaks <- nrow(atacMat)
  atacFullMat <- atacMat

  # Association with motif score
  if(is.null(subInd))
  {
    nonZero <- which(Matrix::rowSums(atacMat)>0)

    nSub <- min(1e4, length(nonZero))
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
  motExp <- (t(matchScores) %*% expectations)
  y <- sweep((t(matchScores) %*% atacMat),1,motExp, FUN="-")
  y <- sweep(y,1,motExp, FUN="/")
  yps <- sapply(1:ncol(backgroundPeaks), function(i){
    B <- sparseMatrix(i=1:nrow(backgroundPeaks), j=backgroundPeaks[,i], x=1,
                      dims=c(nrow(backgroundPeaks), nrow(backgroundPeaks)))
    yp <- ((t(matchScores) %*% B) %*% atacMat - ((t(matchScores) %*% B) %*% expectations)@x)
    yp <- as.matrix(sweep(yp,1,motExp, FUN="/"))
  }, simplify=FALSE)

  yp <- array(unlist(yps), dim = c(ncol(matchScores),
                                   ncol(atacMat), length(yps)))
  ypMean <- apply(yp, c(1, 2), mean)
  ypSd <- apply(yp, c(1, 2), sd)
  scoreMat <- as((y-ypMean)/ypSd, "CsparseMatrix")

  return(list("activity_matrix"=scoreMat,
              "sub_ind"=subInd,
              "expectations"=expectations,
              "background_peaks"=backgroundPeaks))
}

.getContextProjection <- function(atacMat, subSample=1e4){
  if(!is.null(subSample)){
    atacMat <- atacMat[sample(1:nrow(atacMat), min(subSample, nrow(atacMat))), ]
  }

  atacMat <- t(atacMat)

  # Taken from stackexchange:
  # https://stats.stackexchange.com/questions/31565/compute-a-cosine-dissimilarity-matrix-in-r
  sim <- atacMat / sqrt(rowSums(atacMat * atacMat))
  sim <- sim %*% t(sim)
  dSim <- as.dist(1-sim)
  mdsRes <- stats::cmdscale(dSim)
  colnames(mdsRes) <- c("MDS_Context_1", "MDS_Context_2") # add to imports
  rownames(mdsRes) <- rownames(atacMat)

  return(mdsRes)
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
#' @param addLabels Should ChIP-seq peak labels be added to the features.
#' If TRUE features will only be computed for cellular contexts with both, ATAC- and ChIP-seq data.
#' @param whichCol Should features be calculated for all cellular contexts (`"All"`), only the training data (`"OnlyTrain"`)
#' or only for some specific cellular contexts (`"Col"`) specified in `colSel`.
#' @param colSel If `whichCol="colSel"`, name of the cellular context to compute the features for.
#' @param features Names of features to be added. Can be all or some of "Inserts", "Weighted_Inserts", "ChromVAR_Scores", "Cofactor_ChromVAR_Scores", "MDS_Context", "Max_ATAC_Signal".
#' "Insert" features will always be computed.
#' See [TFBlearner::listFeatures] for an overview of the features.
#' @param annoCol Name of column indicating cellular contexts in colData.
#' @param insertionProfile Pre-computed insertion footprint profile for the specified transcription factor.
#' Needs to contain coordinate (chr/seqnames, start, end) columns and weight column (termed "w").
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
#' @importFrom stats cmdscale
#' @export
contextTfFeatures <- function(mae,
                              tfName,
                              addLabels=FALSE,
                              whichCol=c("All", "OnlyTrain", "Col"), # evt. rename these arguments
                              colSel=NULL,
                              features=c("Inserts", "Weighted_Inserts",
                                         "ChromVAR_Scores",
                                         "Cofactor_ChromVAR_Scores",
                                         "MDS_Context",
                                         "Max_ATAC_Signal"),
                              annoCol="context",
                              insertionProfile=NULL,
                              aggregationFun=sum,
                              seed=42,
                              BPPARAM=SerialParam(),
                              ...){

  .checkObject(mae, checkFor=c("Site", "TF"))

  whichCol <- match.arg(whichCol, choices=c("All", "OnlyTrain", "Col"))
  whichContexts <- fifelse(addLabels, "Both", "ATAC")
  if(whichCol=="OnlyTrain"){
    trainCols <- unique(subset(sampleMap(mae), is_training)$colname)
    cols <- lapply(experiments(mae),
                   function(n){colnames(n)[colnames(n) %in% trainCols]})
    maeSub <- subsetByColumn(mae, cols)
    contexts <- getContexts(maeSub, tfName, which=whichContexts)
  }
  else if(whichCol=="Col"){
    if(is.null(colSel)) stop("If features should be computed only for some columns (e.g. cellular contexts,
                              (whichCol=Col), please do provide the names via colSel.")
    maeSub <- mae
    if(addLabels){
      contexts <- getContexts(maeSub, tfName, which=whichContexts)
      contexts <- intersect(colSel, contexts)
    }
    else{
      contexts <- intersect(colSel, colnames(mae[["ATAC"]]))
    }
  }
  else{
    maeSub <- mae
    contexts <- getContexts(maeSub, tfName, which=whichContexts)
  }

  coords <- rowRanges(experiments(maeSub)$Motifs)
  features <- match.arg(features, choices=c("Inserts", "Weighted_Inserts",
                                            "ChromVAR_Scores",
                                            "Cofactor_ChromVAR_Scores",
                                            "MDS_Context",
                                            "Max_ATAC_Signal"),
                        several.ok=TRUE)
  features <- unique(c(features, "Inserts"))

  tfCofactors <- unique(unlist(subset(colData(maeSub[["tfFeat"]]),
                                      tf_name==tfName)$tf_cofactors))

  if(("Cofactor_ChromVAR_Scores" %in% features) & is.null(tfCofactors)){
    msg <- c("No cofactors have been specified when computing transcription ",
             "factor-specific features for ", tfName, ". ")
    msg <- paste0(msg)
    warning(msg)
  }

  atacFragPaths <- unlist(subset(colData(maeSub[["ATAC"]]),
                                 get(annoCol) %in% contexts)$origin)

  # get list of motif ranges
  motifRanges <-  names(experiments(maeSub)[grepl("match_ranges",
                                                  names(experiments(maeSub)))])
  motifRanges <- lapply(experiments(maeSub)[motifRanges], rowRanges)
  names(motifRanges) <- unlist(tstrsplit(names(motifRanges), split="_", keep=3))
  motifRanges <- motifRanges[c(tfName, tfCofactors)]

  if(addLabels){
    colDataChIP <- colData(mae[["ChIP"]])
    colDataChIP <- subset(colDataChIP, tf_name==tfName)
    labelCols <- colDataChIP$combination
    names(labelCols) <- colDataChIP[[annoCol]]
    labels <- lapply(labelCols, function(col){
      as(assays(mae[["ChIP"]])$peaks[,col,drop=TRUE], "CsparseMatrix")})
  }
  else{
    labels <- vector(mode = "list", length = length(contexts))
    names(labels) <- contexts
  }

  # loop over contexts to get the features
  message("Get insert features")
  labels <- labels[contexts] # ensure ordering
  threads <- floor(getDTthreads())/BPPARAM$workers
  feats <- mapply(function(context, labels, coords, atacFrag, motifRanges,
                           features, profile, tfName, aggregationFun,
                           threads, BPPARAM, ...){
    data.table::setDTthreads(threads)

    calcProfile <- FALSE
    if("Weighted_Inserts" %in% features & is.null(profile)){
      calcProfile <- TRUE
    }
    else if("Weighted_Inserts" %in% features & !is.null(profile)){
      message("Using pre-computed insertion-profiles")
    }

    atacFrag <- atacFrag[names(atacFrag)==context]

    addArgs <- list(...)
    addArgs <- addArgs[names(addArgs) %in% c("margin", "shift",
                                             "symmetric", "stranded")]
    args <- c(list(atacData=atacFrag, motifRanges=motifRanges,
                   profiles=profile, calcProfile=calcProfile,
                   subSample=1e8, BPPARAM=BPPARAM),
              addArgs)
    insRes <- suppressWarnings(suppressMessages(do.call(getInsertionProfiles,
                                                        args)))

    scoreCols <- c()
    if("Weighted_Inserts" %in% features) scoreCols <- c(scoreCols,
                                                        "weighted_insert_counts",
                                                        "chi2")
    if("Inserts" %in% features) scoreCols <- c(scoreCols, "insert_counts")

    # aggregate features across ranges of interest
    insFeats <- lapply(scoreCols, function(scoreCol){
      feats <- genomicRangesMapping(coords,
                                    insRes$motifScores,
                                    scoreCol=scoreCol,
                                    byCols="type",
                                    aggregationFun=aggregationFun,
                                    BPPARAM=BPPARAM)
      colnames(feats) <- fifelse(colnames(feats)=="0",
                                 paste(tfName, scoreCol, "margin", sep="_"),
                                 paste(tfName, scoreCol, "within", sep="_"))
      if(scoreCol=="chi2"){
        feats <- list(Matrix::Matrix(rowSums(feats), ncol=1))
        names(feats) <- "chi2_dev_profile"
      }else{
        namesFeats <- colnames(feats)
        feats <- lapply(colnames(feats), function(col) feats[,col,drop=FALSE])
        names(feats) <- namesFeats
      }
      feats
    })
    insFeats <- do.call(c, insFeats)

    if(!is.null(labels)){
      insFeats <- append(insFeats, list("label"=Matrix::Matrix(labels)))
    }

    return(insFeats)
  }, contexts, labels,
     MoreArgs=list(coords=coords, atacFrag=atacFragPaths,
                   motifRanges=motifRanges[[tfName]], features=features,
                   profile=insertionProfile[[tfName]], tfName=tfName,
                   threads=threads, aggregationFun=aggregationFun,
                   BPPARAM=BPPARAM, ...),
     SIMPLIFY=FALSE)
  names(feats) <- contexts

  if("ChromVAR_Scores" %in% features){
    message("Get chromVAR features")

    subInd <- unlist(colData(maeSub[["siteFeat"]])$ChromVAR_sub_ind)
    expectations <- unlist(colData(maeSub[["siteFeat"]])$ChromVAR_expectations)
    bgPeaks <- colData(maeSub[["siteFeat"]])$ChromVAR_background_peaks[[1]]

    if(!is.null(subInd) & !is.null(expectations) & !is.null(bgPeaks)){
     message("ChromVAR features have been pre-computed")
     atacMat <- .convertToMatrix(assays(mae[["ATAC"]])$total_overlaps[,contexts])
    }
    else{
      # full matrix needed for expectation & background calculation
      atacMat <- .convertToMatrix(assays(mae[["ATAC"]])$total_overlaps)
    }

    # get relevant motif columns
    if("Cofactor_ChromVAR_Scores" %in% features){
      tfCols <- c(tfName, tfCofactors)
    }else{
      tfCols <- tfName
    }

    tfReg <- paste0("(^|[^A-Za-z])([.:/_]?(",
                    paste(tfCols, collapse = "|"), "))([^A-Za-z]|$)")
    tfCols <- grep(tfReg, colnames(maeSub[["Motifs"]]), value=TRUE)
    tfCols <- unique(tfCols)

    matchScores <- assays(maeSub[["Motifs"]])$match_scores[,tfCols,drop=FALSE]
    matchScores <- as(as(matchScores, "CsparseMatrix"), "TsparseMatrix")

    colDataMotifs <- subset(colData(maeSub[["Motifs"]]), motif %in% tfCols)
    colDataMotifs <- colDataMotifs[order(match(colDataMotifs$motif, tfCols)),]

    # subset motif matches
    thr <- colDataMotifs$max_score/2
    matchScores@x[matchScores@x<thr[matchScores@j + 1] & matchScores@x<4e4] <- 0

    # retrieve GC-content
    gcContent <- c(assays(experiments(maeSub)$siteFeat)$siteFeat_gc_content[,1])

    res <- .getChromVARScores(atacMat, matchScores, gcContent,
                              colsSel=contexts,
                              subInd=subInd,
                              expectations=expectations,
                              backgroundPeaks=bgPeaks,
                              seed=seed, BPPARAM=BPPARAM, ...)

    chromDevMat <- res$activity_matrix

    assocFeatNames <- paste("ChromVAR_ATAC",
                            gsub('_[0-9]+',"",c("Pearson","Cohen_Kappa")),
                            tfName, sep="_")
    # compute or retrieve association between accessibility and activation across contexts
    if(!all(assocFeatNames %in% colnames(rowData(mae[["tfFeat"]])))){

      # full atac matrix required
      if(ncol(atacMat)!=ncol(experiments(mae)$ATAC)){
        atacMat <- .convertToMatrix(assays(mae[["ATAC"]])$total_overlaps)
      }

      # TODO: Fix this when defining naming conventions
      # motifName <- paste(tfName, "motif", sep="_")
      assocMat <- .getAssociation(atacMat,
                                  chromDevMat[rownames(chromDevMat)==tfName,,
                                              drop=FALSE])
      colnames(assocMat) <- assocFeatNames}
    else{
      message("ChromVAR-Activity ATAC associations have been pre-computed")

      # association features have been pre-computed
      assocMat <- as.matrix(rowData(mae[["tfFeat"]])[,assocFeatNames])
      assocMat <- Matrix::Matrix(assocMat)
      colnames(assocMat) <- assocFeatNames
    }

      actFeats <- lapply(contexts,function(col){
        scoreMat <- chromDevMat[,col, drop=FALSE]
        activityMat <- matrix(scoreMat, nrow=nrow(atacMat), ncol=nrow(scoreMat),
                              byrow=TRUE)
        colnames(activityMat) <- paste("ChromVAR_score",
                                       rownames(scoreMat), sep="_")
        activityMat <- cbind(activityMat, assocMat)
        featNames <- colnames(activityMat)
        activityMat <- lapply(colnames(activityMat), function(featCol){
          activityMat[,featCol,drop=FALSE]})
        names(activityMat) <- featNames
        activityMat
      })
      names(actFeats) <- contexts

      feats <- lapply(contexts, function(context){
        c(feats[[context]], actFeats[[context]])
      })
      names(feats) <- contexts
  }

    if("MDS_Context" %in% features){
      message("Get MDS-Dimensions")
      # check if is present of colData of ATAC
      if(!all(c("MDS_Context_1", "MDS_Context_2") %in%
              colnames(colData(mae[["ATAC"]])))){
        # if not check if full atacMAT is around (as for associ)
        if(!exists("atacMat") || ncol(atacMat)!=ncol(experiments(mae)$ATAC)){
          atacMat <- .convertToMatrix(assays(mae[["ATAC"]])$total_overlaps)
        }

        mdsDim <- .getContextProjection(atacMat)
        mdsDimSub <- as.matrix(mdsDim[contexts,])
      }
      else{
        message("Using pre-computed MDS-dimensions")
        colAtac <- subset(colData(mae[["ATAC"]]), get(annoCol) %in% contexts)
        mdsDim <- colData(mae[["ATAC"]])[,c("MDS_Context_1", "MDS_Context_2")]
        mdsDimSub <- as.matrix(colAtac[match(contexts, colAtac[[annoCol]]),
                                  c("MDS_Context_1", "MDS_Context_2")])
      }
      feats <- lapply(contexts, function(context){
        mdsFeat <- Matrix::Matrix(mdsDimSub[context,],
                                  nrow=length(coords),
                                  ncol=2, byrow=TRUE)
        colnames(mdsFeat) <- colnames(mdsDimSub)
        mdsFeatMat <- lapply(colnames(mdsFeat), function(featCol){
          mdsFeat[,featCol,drop=FALSE]})
        names(mdsFeatMat) <- colnames(mdsFeat)
        c(feats[[context]], mdsFeatMat)
      })
      names(feats) <- contexts
    }

    if("Max_ATAC_Signal" %in% features){
      message("Get maximal ATAC-signal per site")
      if(!("Max_ATAC_Signal" %in% colnames(rowData(mae[["ATAC"]])))){
        if(!exists("atacMat") || ncol(atacMat)!=ncol(experiments(mae)$ATAC)){
          atacMat <-.convertToMatrix(assays(mae[["ATAC"]])$total_overlaps)
        }
        atacMat <- .minMaxNormalization(atacMat, BPPARAM=BPPARAM)
        maxFeat <- Matrix::Matrix(apply(atacMat,1, max), ncol=1)
        colnames(maxFeat) <- "Max_ATAC_Signal"
        maxFeat <- list(maxFeat)
      }
      else{
        message("Using pre-computed maximal ATAC signals")
        maxFeat <- Matrix::Matrix(rowData(mae[["ATAC"]])[,"Max_ATAC_Signal"],
                                  ncol=1)
        colnames(maxFeat) <- "Max_ATAC_Signal"
        maxFeat <- list(maxFeat)
      }
      names(maxFeat) <- "Max_ATAC_Signal"

      feats <- lapply(contexts, function(context){
        c(feats[[context]], maxFeat)
      })
      names(feats) <- contexts
    }

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
      colData(mae[["siteFeat"]])$ChromVAR_sub_ind <- list(res$sub_ind)
      colData(mae[["siteFeat"]])$ChromVAR_expectations <- list(res$expectations)
      colData(mae[["siteFeat"]])$ChromVAR_background_peaks <- list(res$background_peaks)

      rowData(mae[["tfFeat"]]) <- cbind(rowData(mae[["tfFeat"]]),
                                        as.data.frame(as.matrix(assocMat)))
    }

    if("MDS_Context" %in% features){
      colAtac <- colData(mae[["ATAC"]])
      colAtac <- colAtac[,setdiff(colnames(colAtac), c("MDS_Context_1",
                                                       "MDS_Context_2"))]
      colData(mae[["ATAC"]]) <- cbind(colAtac,
                                      mdsDim[match(colAtac[[annoCol]],
                                                   rownames(mdsDim)),])
    }

    if("Max_ATAC_Signal" %in% features){
      rowAtac <- rowData(mae[["ATAC"]])
      rowAtac <- rowAtac[,setdiff(colnames(rowAtac), "Max_ATAC_Signal")]
      maxFeat <- as.data.frame(as.matrix(maxFeat[[1]]))

      if(!is.null(rowAtac)){
        rowAtac <- cbind(rowAtac, maxFeat)}
      else{
        rowAtac <- maxFeat
      }
      rowData(mae[["ATAC"]]) <- rowAtac
    }

  return(mae)
}
