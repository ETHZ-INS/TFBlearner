contextTfFeatures <- function(mae,
                              tfName,
                              tfCofactors=NULL,
                              addLabels=FALSE,
                              whichCol=c("All", "OnlyTrain", "Col"), # evt. rename these arguments
                              colSel=NULL,
                              features=c("Inserts", "Weighted_Inserts",
                                         "Cofactor_Inserts"),
                              annoCol="context",
                              insertionProfile=NULL,
                              aggregationFun=sum,
                              BPPARAM=SerialParam(),
                              ...){

  .checkObject(mae, checkFor=c("Coord", "TF"))

  whichCol <- match.arg(whichCol, choices=c("All", "OnlyTrain", "Col"))
  if(whichCol=="OnlyTrain"){
    maeSub <- mae[,colData(mae)$is_training]
  }
  else if(whichCol=="Col"){
    if(is.null(colSel)) stop("If features should be computed only for some columns (e.g. cellular contexts,
                              (whichCol=Col), please do provide the names via colSel.")
    maeSub <- mae[,colSel,]
  }
  else{
    maeSub <- mae
  }

  features <- match.arg(features, choices=c("Inserts", "Weighted_Inserts",
                                            "Cofactor_Inserts"),
                        several.ok=TRUE)

  # get all cellular contexts covered for that TF
  contexts <- getContexts(maeSub, tfName)
  atacFragPaths <- unlist(subset(colData(experiments(maeSub)$ATAC),
                                 get(annoCol) %in% contexts)$origin)

  coords <- rowRanges(experiments(maeSub)$Motifs)

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
      assays(experiments(maeSub)$ChIP)$peaks[,col,drop=TRUE]})
  }
  else{
    labels <- rep(NULL, length(contexts))
    names(labels) <- contexts
  }

  # loop over contexts to get the features
  labels <- labels[contexts] # ensure ordering
  feats <- BiocParallel::bpmapply(function(context,
                                           labels,
                                           coords,
                                           atacFrag,
                                           motifRanges,
                                           features,
                                           profile,
                                           tfName=tfName,
                                           aggregationFun, ...){
    calcProfile <- FALSE
    if("Weighted_Inserts" %in% features & is.null(profile)){
      calcProfile <- TRUE
    }

    atacFrag <- atacFrag[names(atacFrag)==context]
    insRes <- suppressMessages(getInsertionProfiles(atacFrag, motifRanges,
                                                    profiles=profile,
                                                    calcProfile=calcProfile,
                                                    ...))

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
      colnames(feats) <- paste(tfName, scoreCol, c("margin", "within"), sep="_")
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
                   aggregationFun=aggregationFun, ...),
     SIMPLIFY=FALSE,
     BPPARAM=BPPARAM)

  if("Cofactor_Inserts" %in% features){
    coFeats <- BiocParallel::bplapply(contexts,
                                      function(context,
                                               coords,
                                               atacFrag,
                                               tfCofactors,
                                               motifRanges,
                                               features,
                                               profile,
                                               aggregationFun, ...){

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

          insResCo <- getInsertionProfiles(atacFrag, cofactRanges,
                                           profiles=profile,
                                           calcProfile=calcProfile, ...)
          insFeats <- lapply(scoreCols, function(scoreCol){
            feats <- genomicRangesMapping(coords,
                                           insResCo$motifScores,
                                           scoreCol=scoreCol,
                                           byCols="type",
                                           aggregationFun=aggregationFun)
            colnames(feats) <- paste(coFact, i, scoreCol, c("margin", "within"), sep="_")
            feats
          })
         insFeats <- Reduce("cbind", insFeats[-1], insFeats[[1]])
      })

      insCoFeats <- Reduce(insCoFeats, insCoFeats[-1], insCoFeats[[1]])
      namesFeats <- colnames(insCoFeats)
      insCoFeats <- lapply(namesFeats,
                           function(col) insCoFeats[,col,drop=FALSE])
      names(insCoFeats) <- namesFeats

      return(insCoFeats)
    }, coords=coords, atacFrag=atacFragPaths, tfCofactors=tfCofactors,
       motifRanges=motifRanges[tfCofactors], features=features,
       profile=insertionProfile[tfCofactors],
       aggregationFun=aggregationFun, ...,
       BPPARAM=BPPARAM)
    names(coFeats) <- contexts

    feats <- lapply(contexts, function(context){
      c(coFeats[[context]], feats[[context]])
    })
  }

  names(feats) <- contexts

  # add features back to the full object
  for(context in contexts){
    mae <- .addFeatures(mae, feats[[context]],
                        mapTo="Col", prefix="contextFeat",
                        tfName=tfName,
                        colsToMap=context)}

  return(mae)
}
