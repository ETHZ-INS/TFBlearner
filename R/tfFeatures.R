.getMotifMatches <- function(coords,
                             tfName,
                             tfCofactors,
                             genome=BSgenome.Hsapiens.UCSC.hg38,
                             seqStyle="UCSC"){

  # retrieve motifs
  seqlevelsStyle(coords) <- "UCSC"
  if(genome@metadata$organism=="Homo sapiens"){
    species <- "Hsapiens"
  }
  else if(genome@metadata$organism=="Mus musculus"){
    species <- "Mmusculus"
  }
  else{
    stop("Currently BSgenome needs to be either from Homo sapiens or Mus musculus")
  }

  motifs <- .getNonRedundantMotifs(c(tfName, tfCofactors), species=species)
  motifNames <- names(motifs)
  motifNames[grepl("YY1", motifNames)] <- "YY1"
  names(motifs) <- motifNames

  # get motif matching scores & positions
  matchRanges <- motifmatchr::matchMotifs(motifs,
                                          subject=coords,
                                          out="positions",
                                          genome=genome)

  matchScores <- motifmatchr::matchMotifs(motifs,
                                          subject=coords,
                                          out="scores",
                                          genome=genome)

  # match seqLevelStyle
  motifCoords <- matchRanges[[tfName]]
  tfCofactors <- intersect(names(matchRanges), tfCofactors)
  cofactCoords <- matchRanges[tfCofactors]

  # match seqLevelStyle
  if(!is.null(motifCoords)) seqlevelsStyle(motifCoords) <- seqStyle
  if(!is.null(cofactCoords)) seqlevelsStyle(cofactCoords) <- seqStyle
  if(!is.null(matchScores)) seqlevelsStyle(matchScores) <- seqStyle

  return(list(motifCoords=motifCoords,
              cofactCoords=cofactCoords,
              matchScores=matchScores))
}

.binMat <- function(chIPMat, threshold=NULL){
  cn <- colnames(chIPMat)
  chIPMat <- as(chIPMat, "TsparseMatrix")

  if(!is.null(threshold)){
    #maxV <- max(chIPMat@x)
    #trueL <- which(chIPMat@x==maxV)

    trueL <- which(chIPMat@x>threshold)
    i <- chIPMat@i[trueL]
    j <- chIPMat@i[trueL]

    chIPMat <- sparseMatrix(i=chIPMat@i[trueL]+1,
                            j=chIPMat@j[trueL]+1,
                            x=rep(1, length(trueL)),
                            dims=c(nrow(chIPMat),
                                   ncol(chIPMat)))}
  else{
    chIPMat <- sparseMatrix(i=chIPMat@i+1,
                            j=chIPMat@j+1,
                            x=rep(1, length(chIPMat@x)),
                            dims=c(nrow(chIPMat),
                                   ncol(chIPMat)))
  }
  colnames(chIPMat) <- cn

  return(chIPMat)
}

.aggregate <- function(chIPMat, aggVar=c("tf", "context")){

  aggVar <- match.arg(aggVar, choices=c("tf", "context"))

  cn <- colnames(chIPMat)
  chIPMat <- as(chIPMat, "TsparseMatrix")

  cmDt <- data.table(i=chIPMat@i,
                     j=chIPMat@j,
                     x=chIPMat@x,
                     name=cn[chIPMat@j+1])
  cmDt[,c("tf", "context"):=tstrsplit(name, split="_", keep=1:2)]

  aggDt <- cmDt[,.(value=sum(x>0)),by=c("i", aggVar)]
  cols <- unique(aggDt[[aggVar]])
  cols <- cols[order(cols)]
  aggDt[,j:=as.integer(factor(get(aggVar), levels=cols, ordered=TRUE))]
  aggMat <- sparseMatrix(i=aggDt$i+1,
                         j=aggDt$j,
                         x=aggDt$value,
                         dims=c(nrow(chIPMat),
                                length(cols)))
  aggMat
}

.getBindingPatterns <- function(chIPMat,
                                tfName,
                                tfCofactors,
                                nPatterns=10,
                                aggregate=TRUE,
                                L1=c(0.5,0.5)){

    # Take out TF of interest
    cn <- colnames(chIPMat)
    chIPMat <- chIPMat[, tstrsplit(cn, split="_", keep=2)[[1]]!=tfName, drop=FALSE]
    cn <- colnames(chIPMat)

    cm <- .binMat(chIPMat, threshold=0)

    if(aggregate){
      aggTf <- .aggregate(cm, aggVar="tf")
      nmfTfRes <- suppressMessages(RcppML::nmf(aggTf, k=nPatterns, L1=L1))

      aggCon <- .aggregate(cm, aggVar="context")
      nmfConRes <- suppressMessages(RcppML::nmf(cm, k=nPatterns, L1=L1))

      wTf <- nmfTfRes$w
      colnames(wTf) <- paste("pattern", "tf", 1:ncol(wTf), sep="_")
      wCon <- nmfConRes$w
      colnames(wCon) <- paste("pattern", "context", 1:ncol(wCon), sep="_")
      weights <- list(wTf, wCon)
      #sim <- .jaccard(wTf, aggCon)
      #wTf <- wTf[,]
      #sim <- .jaccard(wCon, aggTf)

      #eventually add coefficients H
      fm <- Reduce(cbind, weights[-1], weights[[1]])

    }
    else{
      nmfRes <- suppressMessages(RcppML::nmf(cm, k=nPatterns, L1=L1))
      fm <- nmfRes$w
    }

    fm <- Matrix::Matrix(fm)

    return(fm)
}

.jaccard <- function(set1, set2)
{
  set1 <- as(set1, "TsparseMatrix")
  set2 <- as(set2, "TsparseMatrix")

  set1Dt <- data.table(i=set1@i,
                       set1_col=set1@j,
                       set1_label=set1@x) #w@x>0
  # get the number of single matches
  set1Dt[,n_set1:=.N, by=set1_col]

  set2Dt <- data.table(i=set2@i,
                       set2_col=set2@j,
                       set2_label=set2@x) # ms@x!=0 !=0
  # get the number of single matches
  set2Dt[,n_set2:=.N, by=set2_col]

  # get the matching
  indDt <- merge(set1Dt, set2Dt, by=c("i"), allow.cartesian=TRUE) # => still needs to be all=TRUE, in case of no matches => then its just zero (alt. all=TRUE, indDt[is.na(indDt)] <- FALSE)
  # indDt[,is_match:=topic_load==tf_match] # it will be always a binary case

  # jaccard index
  cont <- indDt[,.(cont=.N/(data.table::first(n_set1)+data.table::first(n_set2)-.N)), by=c("set1_col", "set2_col")]
  #cont <- indDt[,.(cont=sum(is_match)/(first(n_tf)+first(n_topic))), by=c("tf", "topic")]

  # get names
  cont[,set2_col:=colnames(set2)[set2_col+1]]
  cont[,set1_col:=colnames(set1)[set1_col+1]]

  return(cont)
}


.getNBindings <- function(chIPMat, tfName){

  chIPMat <- chIPMat[,unlist(tstrsplit(colnames(chIPMat), split="_", keep=1))!=tfName]
  chIPMat <- .binMat(chIPMat)

  tfs <- tstrsplit(colnames(chIPMat), split="_", keep=2)[[1]]
  colnames(chIPMat) <- tfs

  # get per TF means
  rm <- sapply(tfs, \(tf){
    Matrix::rowMeans(chIPMat[,colnames(chIPMat)==tf, drop=FALSE])})
  nb <- Matrix::rowSums(rm)
  nb <- Matrix::Matrix(nb, ncol=1)
  colnames(nb) <- "n_chIP_peaks"

  return(nb)
}

.selectMotifs <- function(matchScores, labels, nMotifs=10, subSample=10000)
{
  labels <- .binMat(labels, threshold=0)
  labels <- rowMaxs(labels)

  # per motif threshold: checkout H5sparseMarix
  #thr <- vapply(colnames(matchScores), function(col){
  #  thr <- quantile(matchScores[,col], prob=0.9)
  #})
  thr <- DelayedMatrixStats::colQuantiles(matchScores, probs=0.9)
  #ind <- which(matchScores != 0, arr.ind = TRUE)

  subRows <- sample(1:nrow(matchScores), subSample)
  matchSubScores <- matchScores[subRows,,drop=FALSE]
  matchSubScores <- as(as.matrix(matchSubScores), "TsparseMatrix")
  labels <- labels[subRows]

  # top motif scores
  matchCoScores <- matchSubScores
  matchCoScores@x[matchCoScores@x < thr[matchCoScores@i + 1] ] <- 0
  matchCoScores@x[matchCoScores@x >= thr[matchCoScores@i + 1] ] <- 1

  # get mutually exclusive motif scores
  zeroInd <- which(matchSubScores==0, arr.ind = TRUE)
  matchExScores <- sparseMatrix(i=zeroInd[,1],
                                j=zeroInd[,2],
                                x=rep(1, nrow(zeroInd)),
                                dims=c(nrow(matchSubScores),
                                       ncol(matchSubScores)))

  # jaccard index of mutually exclusive and top co-occuring motifs
  labels <- matrix(labels, nrow=length(labels), ncol=1)
  matchCo <- .jaccard(matchCoScores, labels)
  matchCo[,motif_id:=1:.N]
  setorder(matchCo, -cont)
  topCoMotif <- matchCo$motif_id[1:nMotifs]

  matchEx <- .jaccard(matchExScores, labels)
  matchEx[,motif_id:=1:.N]
  setorder(matchEx, -cont)
  topExMotif <- matchEx$motif_id[1:nMotifs]

  # select top ones
  selectedMotifs <- colnames(matchScores)[unique(c(topCoMotif, topExMotif))]
  # can happen if the motif-matches matrix has less columns than motifs to select
  selectedMotifs <- selectedMotifs[!is.na(selectedMotifs)]

  return(selectedMotifs)
}


.GCSmoothQuantile <- function(gc,
                              counts,
                              g=20,
                              round=FALSE) {
  gcg <- Hmisc::cut2(gc, g=g)
  countsGCSQ <- .gcqn_qsmooth(counts, gcg, round=round)

  countsGCSQ <- Matrix::Matrix(countsGCSQ, sparse = TRUE)

  return(countsGCSQ)
}

.gcqn_qsmooth <- function(counts, gcg, round=FALSE){

  gcn <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts),
                dimnames=list(rownames(counts), colnames(counts)))

  for(ii in 1:nlevels(gcg)){
    id <- which(gcg==levels(gcg)[ii])
    countBin <- counts[id,,drop=FALSE]
    # auantile norm
    normCountBin <- preprocessCore::normalize.quantiles(as.matrix(countBin), copy=FALSE)

    if(round) normCountBin <- round(normCountBin)
    normCountBin[normCountBin<0L] <- 0L
    gcn[id,] <- as.matrix(normCountBin)
  }

  return(gcn)
}

.matrixKappa <- function(mat1, mat2, adjust=TRUE) {
    stopifnot(ncol(mat1)==ncol(mat2))
    # Compute pairwise agreements
    tp <- tcrossprod(mat1,mat2) # 1-1 matches
    fn <- tcrossprod((1-mat1),mat2) # 0-1 mismatches
    fp <- tcrossprod(mat1,(1-mat2)) # 1-0 mismatches
    tn <- ncol(mat1)-(tp+fn+fp) # 0-0 matches
    #tn <- tcrossprod(abs(1-mat1),abs(1-mat2))
    if(!adjust) return((tp+tn)/ncol(mat1))
    Po <- (tp+tn)/ncol(mat1)
    Pe <- (((tp+fn)/ncol(mat1)) * ((tp+fp)/ncol(mat1))) +
      (((tn+fp)/ncol(mat1)) * ((tn+fn)/ncol(mat1)))
    return((Po-Pe)/(1-Pe))
}

.getCoOccuringMotifs <- function(motifCoords, refCoords){
  coCounts <- lapply(names(motifCoords), function(tf){
    mDt <- as.data.table(motifCoords[[tf]])
    mDt$tf <- tf
    coCount <- genomicRangesMapping(refCoords,
                                    mDt[,c("seqnames", "start", "end", "tf")],
                                    byCols=c("tf"))
    colnames(coCount) <- paste("n_cooccuring_motifs", tf, sep="_")
    coCount
  })
  coCounts <- Reduce(cbind, coCounts[-1], coCounts[[1]])

  return(coCounts)
}

.getAssociation <- function(atacMat1, atacMat2){

  pearCor <- cor(Matrix::t(as.matrix(atacMat1)), Matrix::t(as.matrix(atacMat2)))
  colnames(pearCor) <- paste("Pearson", 1:ncol(pearCor),
                             sep="_")
  pearCor <- Matrix::Matrix(pearCor)

  # max-scaling
  atacScalMat1 <- Matrix::t(Matrix::t(atacMat1)/MatrixGenerics::colMaxs(atacMat1))
  atacScalMat2 <- Matrix::t(Matrix::t(atacMat2)/MatrixGenerics::colMaxs(atacMat2))

  # binarizing
  atacBinMat1 <- .binMat(atacScalMat1, threshold=0.3)
  atacBinMat2 <- .binMat(atacScalMat2, threshold=0.3)

  cohKappa <- .matrixKappa(atacBinMat1, atacBinMat2)
  colnames(cohKappa) <- paste("Cohen_Kappa", 1:ncol(cohKappa),
                              sep="_")

  promMat <- cbind(pearCor, cohKappa)
  promMat
}

#' Transcription factor-specific features
#'
#' Adds an experiment with features specific for the specified transcription factor, such as co-occuring motifs and binding patterns, to the provided MultiAssayExperiment object.
#'
#' @name tfFeatures
#' @param mae [MultiAssayExperiment::MultiAssayExperiment-class] as construced by [TFBlearner::prepData()] containing Motif, ATAC-, ChIP-seq
#' and site-specific features as obtained by [TFBlearner::siteFeatures()].
#' @param tfName Name of transcription factor to compute features for.
#' Provided MultiAssayExperiment object needs to contain ChIP-seq data for the specified transcription factor (check `colData(experiments(mae)$ChIP`).
#' @param tfCofactors Names of cofactors (other transcription factors) of the specified transcription factor.
#' @param features Names of features to be added. Can be all or some of "Binding_Patterns", "Promoter_Association",
#' "C_Score", "Cooccuring_Motifs", "Associated_Motifs". Features are stored in the assays of the added experiment.
#' See [TFBlearner::listFeatures] for an overview of the features.
#' @param nPatterns Number of non-negative matrix factorization (NMF) components to consider for the decomposition of the ChIP-seq peaks matrix.
#' Passed to [RcppML::nmf].
#' @param L1 LASSO penalties for lower rank matrices w and h resulting from NMF. Vector with two elements for w and h.
#' Passed to [RcppML::nmf]
#' @param genome [BSgenome::BSgenome-class] to be used.
#' @return [MultiAssayExperiment::MultiAssayExperiment-class] object with an experiment containing transcription factor-specific features added to [MultiAssayExperiment::experiments].
#' If already an experiment with transcription factor-specific features exists in the object, columns for the specified transcription factor are added to it.
#' @param nMotifs Number of associated motifs to select. Chooses `nMotifs` co-occuring and `nMotifs` mutually exclusive motifs.
#' @import MultiAssayExperiment
#' @importFrom motifmatchr matchMotifs
#' @importFrom universalmotif convert_motifs
#' @importFrom TFBSTools PFMatrixList PWMatrixList
#' @importFrom MotifDb MotifDb
#' @importFrom RcppML nmf
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom Hmisc cut2
#' @importFrom DelayedMatrixStats colQuantiles
#' @importFrom MatrixGenerics colMaxs
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @export
tfFeatures <- function(mae,
                       tfName,
                       tfCofactors=NULL,
                       features=c("Binding_Patterns", "Promoter_Association",
                                  "C_Score", "Cooccuring_Motifs",
                                  "Associated_Motifs"),
                       nPatterns=10,
                       L1=c(0.5,0.5),
                       nMotifs=10,
                       genome=BSgenome.Hsapiens.UCSC.hg38){

  features <- match.arg(features, choices=c("Binding_Patterns",
                                            "Promoter_Association",
                                            "C_Score",
                                            "Cooccuring_Motifs",
                                            "Associated_Motifs"),
                        several.ok=TRUE)
  featMats <- list()

  # object validator
  .checkObject(mae, checkFor="Site")

  # reference coordinates
  coords <- rowRanges(experiments(mae)$Motifs)

  # subset to training data
  maeTrain <- mae[,colData(mae)$is_training]

  # get assays: ChIP & ATAC
  atacMat <- assays(experiments(maeTrain)$ATAC)$total_overlaps
  chIPMat <- assays(experiments(maeTrain)$ChIP)$peaks

  # Normalize ATAC-seq data
  message("GC Normalization") # relevant for the computation of other features
  siteFeatName <- "siteFeat" #names(experiments(mae))[[length(experiments(mae))]] # get the experiment added the last
  gcFeatName <- paste(siteFeatName, "gc_content", sep="_")
  gc <- assays(experiments(maeTrain)[[siteFeatName]])[[gcFeatName]][,,drop=TRUE]

  atacNormMat <- .GCSmoothQuantile(gc, atacMat, g=20, round=TRUE)
  assays(experiments(maeTrain)$ATAC)$norm_total_overlaps <- atacNormMat

  if("ATAC_promoters" %in% names(experiments(maeTrain))){
    # prune to standard chromosomes
    promCoords <- rowRanges(experiments(maeTrain)$ATAC_promoters)
    experiments(maeTrain)$ATAC_promoters <-
       keepStandardChromosomes(experiments(maeTrain)$ATAC_promoters,
                               pruning.mode="coarse")

    gc <- Repitools::gcContentCalc(promCoords, genome)
    atacPromMat <- assays(experiments(maeTrain)$ATAC_promoters)$total_overlaps
    atacNormPromMat <- .GCSmoothQuantile(gc, atacPromMat, g=20, round=TRUE)
    assays(experiments(maeTrain)$ATAC_promoters)$norm_total_overlaps <- atacNormPromMat
  }

  message("Get motif match coordinates")
  matchData <- .getMotifMatches(coords,
                                tfName,
                                tfCofactors,
                                seqStyle=seqlevelsStyle(coords),
                                genome=genome)
  motifCoords <- append(GRangesList(matchData$motifCoords),
                        matchData$cofactCoords, after=1)
  names(motifCoords) <- c(tfName, names(matchData$cofactCoords))

  # If not present yet add motif matching ranges to object
  for(tf in names(motifCoords)){
    matchCoords <- motifCoords[[tf]]
    scoreMat <- Matrix::Matrix(matchCoords$score, ncol=1)
    colnames(scoreMat) <-  paste("match_score", tf, sep="_")
    matchFeats <- list(scoreMat)
    prefix <- paste("match_ranges", tf, sep="_")
    names(matchFeats) <- prefix
    if(!(prefix) %in% names(experiments(mae))){
      mae <- .addFeatures(mae, matchFeats, mapTo="Tf",
                          prefix=paste("match_ranges", tf, sep="_"),
                          tfName=tf, coords=matchCoords)}
  }

  # NMF decomposition Binding Patterns
  if("Binding_Patterns" %in% features){
    message("Binding pattern Features")
    bindPattern <- .getBindingPatterns(chIPMat, tfName=tfName,
                                       nPatterns=nPatterns, L1=L1)
    namesBindPattern <- colnames(bindPattern)
    bindPattern <- lapply(namesBindPattern,
                          function(col) bindPattern[,col,drop=FALSE])
    names(bindPattern) <- namesBindPattern
    featMats <- append(featMats, bindPattern)
  }

  if("Promoter_Association" %in% features &
     "ATAC_promoters" %in% names(experiments(maeTrain))){
    message("Promoter association Features")

    isProm <- which(rowData(experiments(mae)$ATAC_promoters)$tf_name==tfName)
    atacSubPromMat <- atacNormPromMat[isProm,,drop=FALSE]

    promAssoc <- .getAssociation(atacNormMat, atacSubPromMat)
    colnames(promAssoc) <- paste("Promoter", colnames(promAssoc), sep="_")
    namesPromAssoc <- colnames(promAssoc)

    promAssoc <- lapply(namesPromAssoc,
                        function(col) promAssoc[,col,drop=FALSE])
    names(promAssoc) <- namesPromAssoc
    featMats <- append(featMats, promAssoc)
  }

  if("C_Score" %in% features){
    message("Crowdedness Scores")
    cScore <- list(.getNBindings(chIPMat, tfName=tfName))
    names(cScore) <- colnames(cScore[[1]])
    featMats <- append(featMats, cScore)
  }

  if("Cooccuring_Motifs" %in% features){
    message("Co-occuring motifs counts")

    coCounts <- .getCoOccuringMotifs(motifCoords, coords)
    if(ncol(coCounts)>1){
      namesCoCounts <- c("n_cooccuring_tf_motifs",
                         paste("n_cooccuring_cofactor_motifs",
                               1:max(1, ncol(coCounts)-1), sep="_"))
    }
    else{
      namesCoCounts <- "n_cooccuring_tf_motifs"
    }

    colnames(coCounts) <- namesCoCounts

    coCounts <- lapply(namesCoCounts,
                       function(col) coCounts[,col,drop=FALSE])
    names(coCounts) <- namesCoCounts
    featMats <- append(featMats, coCounts)
  }

  if("Associated_Motifs" %in% features){
    message("Select motifs")

    isTf <- colData(experiments(maeTrain)$ChIP)$tf_name==tfName
    chIPLabels <- chIPMat[,isTf, drop=FALSE]

    matchScores <- assays(experiments(maeTrain)$Motifs)$match_scores
    tfCols <- grep(paste(c(tfName, tfCofactors),collapse="|"),
                   colnames(matchScores), value=TRUE) # these will be included anyways

    # select motifs
    matchScoresSub <- matchScores[,!c(colnames(matchScores) %in% tfCols), drop=FALSE]
    selectedMotifs <- .selectMotifs(matchScoresSub, chIPLabels, nMotifs=nMotifs)
    selectedMotifs <- unique(c(selectedMotifs, tfCols))

    # important its added back to the full object!
    # might actually not be needed as long as its just done on the training matrix
    colDataChIP <- colData(experiments(mae)$ChIP)
    if(is.null(colDataChIP$preselected_motifs)){
      colDataChIP$preselected_motifs <- replicate(nrow(colDataChIP), list)}

    colData(experiments(mae)$ChIP)$preselected_motifs <-
      fifelse(colDataChIP$tf_name==tfName, list(selectedMotifs),
              colDataChIP$preselected_motifs)

    # add cofactors too
    if(is.null(colDataChIP$tf_cofactors)){
      colDataChIP$tf_cofactors <- replicate(nrow(colDataChIP), list)}

    colData(experiments(mae)$ChIP)$tf_cofactors <-
      fifelse(colDataChIP$tf_name==tfName, list(tfCofactors),
              colDataChIP$tf_cofactors)
  }

  # TODO: Normalize by width
  # normalize full mat by width => have the width normalized features too
  # how to do this => check how this was done earlier

  # add features full mae object: Features are only computed on training data
  if(length(featMats)>0){
    mae <- .addFeatures(mae, featMats, mapTo="Tf",
                        prefix="tfFeat", tfName=tfName)}

  # check if pre-selected motifs are kept: => it is
  # .checkObject(mae, checkFor="TF")

  return(mae)
}
