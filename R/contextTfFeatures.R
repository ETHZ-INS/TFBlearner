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
  # if(is.null(subInd))
  # {
  #   nonZero <- which(Matrix::rowSums(atacMat)>0)
  #   nSub <- min(1e4, length(nonZero))
  #   subInd <- sample(nonZero, nSub)
  # }

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
  colnames(mdsRes) <- paste(mdsDimFeatName, 1:2, sep="_")
  rownames(mdsRes) <- rownames(atacMat)

  return(mdsRes)
}

.getVariableSites <- function(atacMat){
  atacNormMat <- .minMaxNormalization(atacMat, useMax=TRUE)
  # entropy <- function(x){
  #   b <- cut(x, breaks=10, labels=FALSE)
  #   p <- table(b)/length(b)
  #   p <- -sum(log(p)*p)
  # }
  #ent <- apply(atacNormMat, 1, entropy)

  rowVars <- apply(atacNormMat, 1, var)
  varMat <- Matrix::Matrix(rowVars, ncol=1)
  colnames(varMat) <- atacVarFeatName
  return(varMat)
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
#' @param features Names of features to be added. Can be all or some of "Inserts", "Weighted_Inserts", "ChromVAR_Scores", "Cofactor_ChromVAR_Scores", "MDS_Context", "Max_ATAC_Signal", "ATAC_Variance".
#' "Insert" features will always be computed.
#' See [TFBlearner::listFeatures] for an overview of the features.
#' @param annoCol Name of column indicating cellular contexts in colData.
#' @param insertionProfile Pre-computed insertion footprint profile for the specified transcription factor.
#' Needs to contain coordinate (chr/seqnames, start, end) columns and weight column (termed "w").
#' @param aggregationFun function (e.g. mean, median, sum) used to aggregate features across the rowRanges of experiments of the
#' provided [MultiAssayExperiment::MultiAssayExperiment-class] object.
#' @param nVarSites Number of sites with highest ATAC-signal variance to include for ChromVAR deviation scores and MDS projections.
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
                                         "Max_ATAC_Signal",
                                         "ATAC_Variance"),
                              annoCol="context",
                              insertionProfile=NULL,
                              aggregationFun=sum,
                              nVarSites=1e5,
                              seed=42,
                              BPPARAM=SerialParam(),
                              ...){

  .checkObject(mae, checkFor=c("Site", "TF"), tfName=tfName)

  whichCol <- match.arg(whichCol, choices=c("All", "OnlyTrain", "Col"))
  whichContexts <- fifelse(addLabels, "Both", "ATAC")
  if(whichCol=="OnlyTrain"){
    trainCols <- unique(subset(sampleMap(mae), get(isTrainCol))$colname)
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
      contexts <- intersect(colSel, colnames(mae[[atacExp]]))
    }
  }
  else{
    maeSub <- mae
    contexts <- getContexts(maeSub, tfName, which=whichContexts)
  }

  coords <- rowRanges(maeSub[[motifExp]])
  nVarSites <- min(nVarSites, length(coords))

  features <- match.arg(features, choices=c("Inserts", "Weighted_Inserts",
                                            "ChromVAR_Scores",
                                            "Cofactor_ChromVAR_Scores",
                                            "MDS_Context",
                                            "Max_ATAC_Signal",
                                            "ATAC_Variance"),
                        several.ok=TRUE)
  features <- unique(c(features, "Inserts"))

  tfCofactors <- unique(unlist(subset(colData(maeSub[[tfFeat]]),
                                    get(tfNameCol)==tfName)[[tfCofactorsCol]]))

  if(("Cofactor_ChromVAR_Scores" %in% features) & is.null(tfCofactors)){
    msg <- c("No cofactors have been specified when computing transcription ",
             "factor-specific features for ", tfName, ". ")
    msg <- paste0(msg)
    warning(msg)
  }

  atacFragPaths <- unlist(subset(colData(maeSub[[atacExp]]),
                                 get(annoCol) %in% contexts)$origin)

  # get list of motif ranges, this will eventually be refactored anyways
  motifRanges <-  names(experiments(maeSub)[grepl("match_ranges",
                                                  names(experiments(maeSub)))])
  motifRanges <- lapply(experiments(maeSub)[motifRanges], rowRanges)
  names(motifRanges) <- unlist(tstrsplit(names(motifRanges), split="_", keep=3))
  motifRanges <- motifRanges[c(tfName, tfCofactors)]

  if(addLabels){
    colDataChIP <- colData(mae[[chIPExp]])
    colDataChIP <- subset(colDataChIP, get(tfNameCol)==tfName)
    labelCols <- colDataChIP$combination
    names(labelCols) <- colDataChIP[[annoCol]]
    labels <- lapply(labelCols, function(col){
     as(assays(mae[[chIPExp]])[[peakAssayName]][,col,drop=TRUE],
        "CsparseMatrix")})
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
                                                        wInsertsFeatName,
                                                        devFeatName)
    if("Inserts" %in% features) scoreCols <- c(scoreCols, insertFeatName)

    # aggregate features across ranges of interest
    insFeats <- lapply(scoreCols, function(scoreCol){
      feats <- genomicRangesMapping(coords,
                                    insRes$motifScores,
                                    scoreCol=scoreCol,
                                    byCols="type",
                                    aggregationFun=aggregationFun,
                                    BPPARAM=BPPARAM)
      colnames(feats) <- fifelse(colnames(feats)=="0",
                                 paste("tf", scoreCol, "margin", sep="_"),
                                 paste("tf", scoreCol, "within", sep="_"))
      if(scoreCol==devFeatName){
        feats <- list(Matrix::Matrix(rowSums(feats), ncol=1))
        names(feats) <- devFeatName
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

  if("ATAC_Variance" %in% features | "ChromVAR_Scores" %in% features |
     "MDS_Context" %in% features){
    message("Get site-specific variance in ATAC-signal")
    if(!("ATAC_Variance" %in% colnames(rowData(mae[[atacExp]])))){
     atacMat <- .convertToMatrix(assays(mae[[atacExp]])$total_overlaps)
     varFeat <- .getVariableSites(atacMat)
     rs <- Matrix::rowSums(atacMat)

     topVarSites <- setdiff(order(-varFeat[,1])[1:nVarSites], which(rs==0))
     colData(mae[[siteFeat]])[[topVarSitesName]] <- list(topVarSites)
     varFeat <- list(varFeat)
    }
    else{
     message("Using pre-computed variances")
     varFeat <- Matrix::Matrix(rowData(mae[[atacExp]])[,atacVarFeatName],ncol=1)
     topVarSites <- unlist(colData(mae[[siteFeat]])[[topVarSitesName]])
     colnames(varFeat) <- atacVarFeatName
     varFeat <- list(varFeat)
    }
    names(varFeat) <- atacVarFeatName
  }

  if("ATAC_Variance" %in% features){
    feats <- lapply(contexts, function(context){
      c(feats[[context]], varFeat)
    })
    names(feats) <- contexts
  }

  if("ChromVAR_Scores" %in% features){
    message("Get chromVAR features")

    subInd <- topVarSites
    expectations <- colData(maeSub[[siteFeat]])[[chromVarExpName]][[1]]
    bgPeaks <- colData(maeSub[[siteFeat]])[[chromVarBgName]][[1]]

    if(!is.null(subInd) & !is.null(expectations) & !is.null(bgPeaks)){
     message("ChromVAR features have been pre-computed")
     atacMat <- assays(mae[[atacExp]])[[totalOverlapsName]][,contexts]
     atacMat <- .convertToMatrix(atacMat)
    }
    else if(!exists("atacMat") || ncol(atacMat)!=ncol(mae[[atacExp]])){
      atacMat <-.convertToMatrix(assays(mae[[atacExp]])[[totalOverlapsName]])
    }

    # get relevant motif columns
    if("Cofactor_ChromVAR_Scores" %in% features){
      tfCols <- unlist(subset(colData(mae[[tfFeat]]),
                       get(tfNameCol)==tfName)[[preSelMotifCol]])
    }else{
      selMotifs <- unlist(subset(colData(mae[[tfFeat]]),
                                 get(tfNameCol)==tfName)[[preSelMotifCol]])
      tfCols <- selMotifs[grep(tfMotifPrefix, names(selMotifs))]
    }

    matchScores <- assays(mae[[motifExp]])[[matchAssayName]][,tfCols,drop=FALSE]
    matchScores <- as(as(matchScores, "CsparseMatrix"), "TsparseMatrix")

    colDataMotifs <- subset(colData(maeSub[[motifExp]]),
                            get(motifNameCol) %in% tfCols)
    colDataMotifs <- colDataMotifs[order(match(colDataMotifs[[motifNameCol]],
                                               tfCols)),]

    # subset motif matches
    thr <- colDataMotifs[[maxScoreCol]]/2
    matchScores@x[matchScores@x<thr[matchScores@j + 1] & matchScores@x<4e4] <- 0

    # retrieve GC-content
    gcContent <- assays(maeSub[[siteFeat]])[[paste(siteFeat,
                                                   gcContFeatName,
                                                   sep="_")]][,,drop=TRUE]
    colnames(matchScores) <- names(tfCols)

    res <- .getChromVARScores(atacMat, matchScores, gcContent,
                              colsSel=contexts,
                              subInd=subInd,
                              expectations=expectations,
                              backgroundPeaks=bgPeaks,
                              seed=seed, BPPARAM=BPPARAM, ...)

    chromDevMat <- res$activity_matrix

    # TODO: Fix motif naming convention
    isTfCol <- which(tfCols==tfName)
    assocFeatNames <- paste(chromVarAssocSuffix,
                            c(assocPearsonPrefix, assocCohenPrefix),
                            tfName, sep="_")
    # compute or retrieve association between accessibility and activation across contexts
    if(!all(assocFeatNames %in% colnames(rowData(mae[[tfFeat]])))){

      # full atac matrix required
      if(ncol(atacMat)!=ncol(mae[[atacExp]])){
        atacMat <- .convertToMatrix(assays(mae[[ataxExp]])[[totalOverlapsName]])
      }

      # TODO: Fix this when defining naming conventions
      # motifName <- paste(tfName, "motif", sep="_")
      assocMat <- .getAssociation(atacMat, chromDevMat[isTfCol,,drop=FALSE])}
    else{
      message("ChromVAR-Activity ATAC associations have been pre-computed")

      # association features have been pre-computed
      assocMat <- as.matrix(rowData(mae[[tfFeat]])[,assocFeatNames])
      assocMat <- Matrix::Matrix(assocMat)}

      colnames(assocMat) <- paste(colnames(assocMat),
                                  chromVarAssocSuffix, "tf", sep="_")

      actFeats <- lapply(contexts,function(col){
        scoreMat <- chromDevMat[,col, drop=FALSE]
        activityMat <- matrix(scoreMat, nrow=nrow(atacMat), ncol=nrow(scoreMat),
                              byrow=TRUE)
        colnames(activityMat) <- paste(chromVarScoreSuffix,
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
      mdsDimNames <- paste(mdsDimFeatName, 1:2, sep="_")
      if(!all(mdsDimNames %in% colnames(colData(mae[[atacExp]])))){
        # if not check if full atacMAT is around (as for associ)
        if(!exists("atacMat") || ncol(atacMat)!=ncol(mae[[atacExp]])){
          atacMat <- assays(mae[[atacExp]])[[totalOverlapsName]]
          atacMat <- .convertToMatrix(atacMat)
        }

        mdsDim <- .getContextProjection(atacMat[topVarSites,], subSample=NULL)
        mdsDimSub <- as.matrix(mdsDim[contexts,])
      }
      else{
        message("Using pre-computed MDS-dimensions")
        colAtac <- subset(colData(mae[[atacExp]]), get(annoCol) %in% contexts)
        mdsDim <- colData(mae[[atacExp]])[,mdsDimNames]
        mdsDimSub <- as.matrix(colAtac[match(contexts, colAtac[[annoCol]]),
                                       mdsDimNames])
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
      if(!("Max_ATAC_Signal" %in% colnames(rowData(mae[[atacExp]])))){
        if(!exists("atacMat") || ncol(atacMat)!=ncol(mae[[atacExp]])){
          atacMat <- assays(mae[[atacExp]])[[totalOverlapsName]]
          atacMat <-.convertToMatrix(atacMat)
        }
        atacMat <- .minMaxNormalization(atacMat)
        maxFeat <- Matrix::Matrix(apply(atacMat,1, max), ncol=1)
        colnames(maxFeat) <- maxAtacFeatName
        maxFeat <- list(maxFeat)
      }
      else{
        message("Using pre-computed maximal ATAC signals")
        maxFeat <- Matrix::Matrix(rowData(mae[[atacExp]])[,maxAtacFeatName],
                                  ncol=1)
        colnames(maxFeat) <- maxAtacFeatName
        maxFeat <- list(maxFeat)
      }
      names(maxFeat) <- maxAtacFeatName

      feats <- lapply(contexts, function(context){
        c(feats[[context]], maxFeat)
      })
      names(feats) <- contexts
    }

    # add features back to the full object
    for(context in contexts){
      featMats <- lapply(feats[[context]], `colnames<-`, NULL)
      names(featMats) <- paste(contextTfFeat, names(featMats), sep="_")

      seTfFeat <- SummarizedExperiment(assays=featMats,
                                       rowRanges=coords)
      colnames(seTfFeat) <- paste(context, tfName, sep="_")
      colData(seTfFeat)[[featTypeCol]] <- contextTfFeat
      colData(seTfFeat)[[tfNameCol]] <- tfName
      mae <- .addFeatures(mae, seTfFeat, colsToMap=context,
                          prefix=contextTfFeat)
    }

    if("ChromVAR_Scores" %in% features |
       "Cofactor_ChromVAR_Scores" %in% features){

      # add chromVAR parameters to object
      colData(mae[[siteFeat]])[[chromVarExpName]] <- list(res$expectations)
      colData(mae[[siteFeat]])[[chromVarBgName]] <- list(res$background_peaks)

      assocFeat <- as.data.frame(as.matrix(assocMat))
      colnames(assocFeat) <- assocFeatNames
      rowData(mae[[tfFeat]]) <- cbind(rowData(mae[[tfFeat]]), assocFeat)
    }

    if("MDS_Context" %in% features){
      colAtac <- colData(mae[[atacExp]])
      colAtac <- colAtac[,setdiff(colnames(colAtac), mdsDimNames)]
      colData(mae[[atacExp]]) <- cbind(colAtac,
                                      mdsDim[match(colAtac[[annoCol]],
                                                   rownames(mdsDim)),])
    }

    if("Max_ATAC_Signal" %in% features){
      rowAtac <- rowData(mae[[atacExp]])
      rowAtac <- rowAtac[,setdiff(colnames(rowAtac), maxAtacFeatName)]
      maxFeat <- as.data.frame(as.matrix(maxFeat[[1]]))

      if(!is.null(rowAtac)){
        rowAtac <- cbind(rowAtac, maxFeat)}
      else{
        rowAtac <- maxFeat
      }
      rowData(mae[[atacExp]]) <- rowAtac
    }

  if("ATAC_Variance" %in% features){
    rowAtac <- rowData(mae[[atacExp]])
    rowAtac <- rowAtac[,setdiff(colnames(rowAtac), atacVarFeatName)]
    varFeat <- as.data.frame(as.matrix(varFeat[[1]]))
    colnames(varFeat) <- atacVarFeatName

    if(!is.null(rowAtac)){
      rowAtac <- cbind(rowAtac, varFeat)}
    else{
      rowAtac <- varFeat
    }
    rowData(mae[[atacExp]]) <- rowAtac
  }

  return(mae)
}
