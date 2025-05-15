.doOneScan <- function(mo, coords, genome, lowp=1e-03){ #1e-02
  # low-confidence matches
  po <- matchMotifs(mo, subject=coords, genome=genome,
                    out="positions", p.cutoff=lowp)[[1]]
  # trace back the original DHS behind each
  po$orig <- to(findOverlaps(po, coords, minoverlap=2L))
  # keep only the strongest hit per DHS
  po <- po[order(seqnames(po), -po$score)]
  po <- po[!duplicated(po$orig)]
  po$orig <- NULL
  # find all below the default cutoff
  po2 <- matchMotifs(mo, subject=coords, genome=genome, out="positions")[[1]]
  po <- po[!overlapsAny(po,po2)]
  sort(c(po,po2))
}

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

  # get motif matching positions
  matchRanges <- lapply(motifs, .doOneScan, coords=coords, genome=genome)
  matchRanges <- as(matchRanges, "GRangesList")

  # match seqLevelStyle
  motifCoords <- matchRanges[[tfName]]
  tfCofactors <- intersect(names(matchRanges), tfCofactors)
  cofactCoords <- matchRanges[tfCofactors]

  # match seqLevelStyle
  if(length(motifCoords)>0) seqlevelsStyle(motifCoords) <- seqStyle
  if(length(cofactCoords)>0) seqlevelsStyle(cofactCoords) <- seqStyle

  return(list(motifCoords=motifCoords,
              cofactCoords=cofactCoords))
}

.binMat <- function(mat, threshold=NULL){
  cn <- colnames(mat)
  mat <- as(mat, "TsparseMatrix")

  if(is.null(threshold)) threshold <- 0
  ind <- which(mat@x>threshold)
  mat <- sparseMatrix(i=mat@i[ind]+1,
                      j=mat@j[ind]+1,
                      x=rep(1, length(ind)),
                      dims=c(nrow(mat), ncol(mat)))
  mat <- Matrix::Matrix(mat)
  colnames(mat) <- cn

  return(mat)
}

.toSparseMatrix <- function(mat, threshold=NULL, binarize=TRUE){
  if(is.null(threshold)) threshold <- 0
  inds <- DelayedArray::blockApply(mat, function(block, binarize, threshold){
    ind <- which(as.logical(block>threshold))
    if(binarize){
      x <- rep(1, length(ind))}
    else{
      x <- block[ind]}
    indDt <- data.table(i=ind, x=x)
  }, grid=colAutoGrid(mat, ncol=1),
  binarize=binarize,
  threshold=threshold)

  indDt <- rbindlist(inds, idcol="j")
  mat <- Matrix::sparseMatrix(i=indDt$i, j=indDt$j, x=indDt$x)

  gc()

  return(mat)
}

.aggregate <- function(chIPMat,
                       threshold=NULL,
                       aggFun=function(x){sum(x)/length(x)},
                       aggVar=c("tf", "context")){

   chIPMat <- as(chIPMat, "TsparseMatrix")
   cn <- as.factor(colnames(chIPMat))
   aggVar <- match.arg(aggVar, choices=c("tf", "context"))
   var <- factor(unlist(tstrsplit(cn, split="_",
                                  keep=fifelse(aggVar=="tf", 2, 1))))

   if(is.null(threshold)) threshold <- 0

   toKeep <- which(chIPMat@x>threshold)
   colsInd <- chIPMat@j[toKeep]+1
   chIPIndDt <- data.table(i=chIPMat@i[toKeep]+1,
                           j=colsInd,
                           x=chIPMat@x[toKeep],
                           aggVar=var[colsInd])

   aggDt <- chIPIndDt[,.(value=aggFun(x)), by=.(i, aggVar)]
   rm(chIPIndDt)

   aggCols <- unique(aggDt$aggVar)
   aggDt[,j:=as.integer(factor(aggVar, levels=aggCols, ordered=TRUE))]
   aggMat <- sparseMatrix(i=aggDt$i,
                          j=aggDt$j,
                          x=aggDt$value,
                          dims=c(nrow(chIPMat),
                                 length(aggCols)))
   aggMat <- Matrix::Matrix(aggMat)
   colnames(aggMat) <- aggCols
   aggMat
}

.getBindingPatterns <- function(chIPMat,
                                nPatterns=50,
                                aggFun=max,
                                binarize=FALSE,
                                L1=c(0.5,0.5),
                                seed=42){
    set.seed(seed)

    if(binarize){
      chIPMat <- .binMat(chIPMat, threshold=0)
    }

    if(is.null(aggFun)){
      nmfRes <- suppressMessages(RcppML::nmf(chIPMat, k=nPatterns,
                                             L1=L1, seed=seed))
      fm <- nmfRes$w
      colnames(fm) <- paste("pattern", "tf", 1:ncol(fm), sep="_")
    }
    else{
      aggMatTf <- .aggregate(chIPMat, aggVar="tf", aggFun=aggFun)
      nmfTfRes <- suppressMessages(RcppML::nmf(aggMatTf, k=nPatterns, L1=L1,
                                               seed=seed))
      gc()

      aggMatCon <- .aggregate(chIPMat, aggVar="context", aggFun=aggFun)
      nmfConRes <- suppressMessages(RcppML::nmf(aggMatCon, k=nPatterns, L1=L1,
                                                seed=seed))
      gc()

      wTf <- nmfTfRes$w
      colnames(wTf) <- paste("pattern", "tf", 1:ncol(wTf), sep="_")
      wCon <- nmfConRes$w
      colnames(wCon) <- paste("pattern", "context", 1:ncol(wCon), sep="_")
      weights <- list(wTf, wCon)

      #eventually add coefficients H
      fm <- Reduce(cbind, weights[-1], weights[[1]])
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
                       set1_label=set1@x)

  # get the number of single matches
  set1Dt[,n_set1:=.N, by=set1_col]

  set2Dt <- data.table(i=set2@i,
                       set2_col=set2@j,
                       set2_label=set2@x)

  # get the number of single matches
  set2Dt[,n_set2:=.N, by=set2_col]

  # get the matching
  indDt <- merge(set1Dt, set2Dt, by=c("i"), allow.cartesian=TRUE)

  # jaccard index
  cont <- indDt[,.(cont=.N/(data.table::first(n_set1)+
                            data.table::first(n_set2)-.N)),
                by=c("set1_col", "set2_col")]
  # get names
  cont[,set2_col:=colnames(set2)[set2_col+1]]
  cont[,set1_col:=colnames(set1)[set1_col+1]]

  return(cont)
}

# pass mae non-test in here => in all the functions it makes it much more readible
.getCrowdedness <- function(chIPMat){
  aggMat <- .aggregate(chIPMat, aggVar="tf")
  cScoreMat <- Matrix::Matrix(Matrix::rowSums(aggMat), ncol=1)
  colnames(cScoreMat) <- "c_score"

  return(cScoreMat)
}

.selectMotifs <- function(matchScores, maxScores, labels, nMotifs=10,
                          subSample=10000)
{
  labels <- .binMat(labels, threshold=0)
  labels <- as(labels, "CsparseMatrix")
  labels <- sparseMatrixStats::rowMaxs(labels)

  thr <- maxScores/2

  subRows <- sample(1:nrow(matchScores), min(subSample, nrow(matchScores)))
  matchSubScores <- matchScores[subRows,,drop=FALSE]
  matchSubScores <- as(as.matrix(matchSubScores), "TsparseMatrix")
  labels <- labels[subRows]

  # top motif scores
  matchCoScores <- matchSubScores
  matchCoScores@x[matchCoScores@x < thr[matchCoScores@j + 1] &
                    matchCoScores@x<4e4] <- 0
  matchCoScores@x[matchCoScores@x >= thr[matchCoScores@j + 1] &
                    matchCoScores@x>=4e4] <- 1

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

.getCofactorBindings <- function(chIPMat, tfCofactors){
  tfCols <- unlist(tstrsplit(colnames(chIPMat), split="_", keep=2))
  tfCofactors <- intersect(tfCols, tfCofactors)

  if(length(tfCofactors)>0){
    cofactBindings <- lapply(tfCofactors, function(tfCol){
      cofactBinding <- Matrix::Matrix(
        Matrix::rowMeans(chIPMat[,tfCols==tfCol, drop=FALSE]), ncol=1)
      colnames(cofactBinding) <- paste("cofactor_binding", tfCol, sep="_")
      cofactBinding})
    names(cofactBindings) <- paste("cofactor_binding", tfCofactors, sep="_")
    return(cofactBindings)}
  else{
    return(NULL)
  }
}

.GCSmoothQuantile <- function(gc, counts, nBins=20, round=FALSE) {
  gcBins <- cut(gc, breaks=nBins)
  counts <- as.matrix(counts)

  # loop over the bins
  for(b in 1:nlevels(gcBins)){
    binId <- which(gcBins==levels(gcBins)[b])
    countBin <- counts[binId,,drop=FALSE]

    # auantile normalize
    normCountBin <- preprocessCore::normalize.quantiles(countBin, copy=FALSE)
    if(round) normCountBin <- as.integer(round(normCountBin))
    normCountBin[normCountBin<0L] <- 0L
    counts[binId,] <- normCountBin
  }
  normCounts <- Matrix::Matrix(counts)
  return(normCounts)
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

  # subset to common contexts
  commonContexts <- intersect(colnames(atacMat1), colnames(atacMat2))
  atacMat1 <- atacMat1[,commonContexts,drop=FALSE]
  atacMat2 <- atacMat2[,commonContexts,drop=FALSE]

  atacMat1 <- as(atacMat1, "CsparseMatrix")
  atacMat2 <- as(atacMat2, "CsparseMatrix")

  pearCor <- cor(Matrix::t(as.matrix(atacMat1)), Matrix::t(as.matrix(atacMat2)))
  colnames(pearCor) <- paste("Pearson", 1:ncol(pearCor), sep="_")
  pearCor <- Matrix::Matrix(pearCor)

  # max-scaling
  atacScalMat1 <- Matrix::t(Matrix::t(atacMat1)/sparseMatrixStats::colMaxs(atacMat1))
  atacScalMat2 <- Matrix::t(Matrix::t(atacMat2)/sparseMatrixStats::colMaxs(atacMat2))

  # binarizing
  atacBinMat1 <- .binMat(atacScalMat1, threshold=0.3)
  atacBinMat2 <- .binMat(atacScalMat2, threshold=0.3)

  cohKappa <- .matrixKappa(atacBinMat1, atacBinMat2)
  colnames(cohKappa) <- paste("Cohen_Kappa", 1:ncol(cohKappa), sep="_")

  promMat <- cbind(pearCor, cohKappa)
  promMat <- Matrix::Matrix(promMat)
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
#' "C_Score", "Cofactor_Binding", "CTCF_Signal", "Cooccuring_Motifs", "Associated_Motifs". Features are stored in the assays of the added experiment.
#' See [TFBlearner::listFeatures] for an overview of the features.
#' @param nPatterns Number of non-negative matrix factorization (NMF) components to consider for the decomposition of the ChIP-seq peaks matrix.
#' Passed to [RcppML::nmf].
#' @param L1 LASSO penalties for lower rank matrices w and h resulting from the NMF decomposition. Vector with two elements for w and h.
#' Passed to [RcppML::nmf]
#' @param nMotifs Number of associated motifs to select. Chooses `nMotifs` co-occuring and `nMotifs` mutually exclusive motifs.
#' @param seed Integer value for setting the seed for random number generation with [base::set.seed].
#' @param genome [BSgenome::BSgenome-class] to be used.
#' @return [MultiAssayExperiment::MultiAssayExperiment-class] object with an experiment containing transcription factor-specific features added to [MultiAssayExperiment::experiments].
#' If already an experiment with transcription factor-specific features exists in the object, columns for the specified transcription factor are added to it.
#' @import MultiAssayExperiment
#' @importFrom motifmatchr matchMotifs
#' @importFrom universalmotif convert_motifs
#' @importFrom TFBSTools PFMatrixList PWMatrixList
#' @importFrom MotifDb MotifDb
#' @importFrom RcppML nmf
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom sparseMatrixStats colMaxs rowMaxs
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @export
tfFeatures <- function(mae,
                       tfName,
                       tfCofactors=NULL,
                       features=c("Binding_Patterns", "Promoter_Association",
                                  "C_Score", "Cooccuring_Motifs",
                                  "Cofactor_Binding", "CTCF_Signal",
                                  "Associated_Motifs"),
                       nPatterns=50,
                       L1=c(0.5,0.5),
                       nMotifs=10,
                       seed=42,
                       genome=BSgenome.Hsapiens.UCSC.hg38){
  set.seed(seed)

  allTfs <- unique(colData(mae[["ChIP"]])$tf_name)
  if(!(tfName %in% allTfs)){
    stop(paste("Transcription factor provided -", tfName,
    "- has no corresponding experiments in the provided MultiAssayExperiment
    object.\n Tfs are named in the following way:",
    paste(head(allTfs), collapse=","), "..."))
  }

  features <- match.arg(features, choices=c("Binding_Patterns",
                                            "Promoter_Association",
                                            "C_Score",
                                            "Cofactor_Binding",
                                            "CTCF_Signal",
                                            "Cooccuring_Motifs",
                                            "Associated_Motifs"),
                        several.ok=TRUE)
  featMats <- list()

  # validate MultiAssayExperiment object
  .checkObject(mae, checkFor="Site")

  # reference coordinates
  coords <- rowRanges(experiments(mae)$Motifs)

  # assay-matrices
  atacMat <- .convertToMatrix(assays(mae[["ATAC"]])$total_overlaps)
  colnames(atacMat) <- colnames(mae[["ATAC"]])
  whichCol <- which(mae[["ChIP"]]$tf_name!=tfName)
  chIPMat <- as(assays(mae[["ChIP"]])$peaks[,whichCol],"CsparseMatrix")
  colnames(chIPMat) <- paste(colData(mae[["ChIP"]])[whichCol,c("context")],
                             colData(mae[["ChIP"]])[whichCol,c("tf_name")],
                             sep="_")

  # Normalize ATAC-seq data
  message("GC Normalization")
  gc <- assays(mae[["siteFeat"]])$siteFeat_gc_content[,,drop=TRUE]

  atacMat <- .GCSmoothQuantile(gc, atacMat, nBins=20, round=TRUE)
  assays(mae[["ATAC"]], withDimnames=FALSE)$norm_total_overlaps <- atacMat

  if("ATAC_promoters" %in% names(experiments(mae))){
    # prune to standard chromosomes
    promCoords <- rowRanges(mae[["ATAC_promoters"]])
    mae[["ATAC_promoters"]] <- keepStandardChromosomes(mae[["ATAC_promoters"]],
                                                       pruning.mode="coarse")

    gc <- Repitools::gcContentCalc(promCoords, genome)
    atacPromMat <- assays(mae[["ATAC_promoters"]])$total_overlaps
    atacPromMat <- .GCSmoothQuantile(gc, atacPromMat, nBins=20, round=TRUE)
    assays(mae[["ATAC_promoters"]])$norm_total_overlaps <- atacPromMat
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

  # Add motif matching ranges to object
  for(tf in names(motifCoords)){
    matchCoords <- motifCoords[[tf]]
    scoreMat <- Matrix::Matrix(matchCoords$score, ncol=1)
    colnames(scoreMat) <-  "all"
    scoreMat <- list(scoreMat)
    names(scoreMat) <- paste("match_score", tf, sep="_")
    prefix <- paste("match_ranges", tf, sep="_")
    if(!(prefix) %in% names(experiments(mae))){
      seMatch <- SummarizedExperiment(assays=scoreMat, rowRanges=matchCoords)
      colData(seMatch)$feature_type <- "motif_matches"
      rownames(colData(seMatch)) <- "all"
      contextsTf <- getContexts(mae, tfName=tfName, which="ChIP")
      mae <- .addFeatures(mae, seMatch, colsToMap=contextsTf, prefix=prefix)
    }
  }

  # NMF decomposition Binding Patterns
  if("Binding_Patterns" %in% features){
    message("Binding pattern Features")
    bindPatterns <- .getBindingPatterns(chIPMat, nPatterns=nPatterns,
                                        L1=L1, seed=seed)
    colNamesPatterns <- colnames(bindPatterns)
    bindPatterns <- lapply(colNamesPatterns,
                           function(col){ bindPatterns[,col,drop=FALSE]})
    names(bindPatterns) <- colNamesPatterns
    featMats <- append(featMats, bindPatterns)
  }

  if("Promoter_Association" %in% features &
     "ATAC_promoters" %in% names(experiments(mae))){
    message("Promoter association Features")

    isProm <- which(rowData(mae[["ATAC_promoters"]])$tf_name==tfName)
    atacPromMat <- atacPromMat[isProm,,drop=FALSE]

    promAsc <- .getAssociation(atacMat, atacPromMat)
    colnames(promAsc) <- paste("promoter", colnames(promAsc), sep="_")
    colNamesPromAsc <- colnames(promAsc)

    promAsc <- lapply(colNamesPromAsc, function(col) promAsc[,col,drop=FALSE])
    names(promAsc) <- colNamesPromAsc
    featMats <- append(featMats, promAsc)
  }

  if("Cofactor_Binding" %in% features){
    message("Cofactor Bindings")
    if(is.null(tfCofactors)){
      stop("Please provide cofactor names (`tfCofactors`) if Cofactor_Bindings should be computed.")}
    cofactBindings <- .getCofactorBindings(chIPMat, tfCofactors)
    if(!is.null(cofactBindings)){
      featMats <- append(featMats, cofactBindings)
    }
  }

  if("C_Score" %in% features){
    message("Crowdedness Scores")
    cScore <- list(.getCrowdedness(chIPMat))

    names(cScore) <- colnames(cScore[[1]])
    featMats <- append(featMats, cScore)
  }

  if("Cooccuring_Motifs" %in% features){
    message("Co-occuring motifs counts")

    coCounts <- .getCoOccuringMotifs(motifCoords, coords)
    namesCoCounts <- colnames(coCounts)

    coCounts <- lapply(namesCoCounts, function(col) coCounts[,col,drop=FALSE])
    names(coCounts) <- namesCoCounts
    featMats <- append(featMats, coCounts)
  }

  if("Associated_Motifs" %in% features){
    message("Select motifs")

    # get (except from test experiments) ChIP-labels of TF
    chIPCols <- colnames(mae[["ChIP"]])
    tfCols <- chIPCols[colData(mae[["ChIP"]])$tf_name==tfName]
    isTesting <- subset(sampleMap(mae), assay=="ChIP" & is_testing)$colname
    tfCols <- intersect(chIPCols[!(chIPCols %in% isTesting)], tfCols)

    if(length(tfCols)>0){
      labels <- .convertToMatrix(assays(mae[["ChIP"]])$peaks[,tfCols,drop=FALSE])
    }
    else{
      # TODO: motif naming
      tfCol <- paste(tfName, "motif", sep="_")
      labels <- Matrix::Matrix(
        assays(mae[["Motifs"]])$match_scores[,tfCol,drop=FALSE], ncol=1)
    }

    matchScores <- assays(mae[["Motifs"]])$match_scores
    tfMotifCols <- grep(paste(c(tfName, tfCofactors),collapse="|"),
                        colnames(matchScores), value=TRUE)

    # select motifs co-occuring around ChIP-peaks or motif matches of TF of interest
    matchScores <- matchScores[,!c(colnames(matchScores) %in% tfMotifCols), drop=FALSE]
    colDataMotifs <- colData(mae[["Motifs"]])
    colDataMotifs <- subset(colDataMotifs, motif %in% colnames(matchScores))
    colDataMotifs <- colDataMotifs[order(match(colDataMotifs$motif,
                                               colnames(matchScores))),,
                                   drop=FALSE]

    maxScores <- colDataMotifs$max_score
    selMotifs <- .selectMotifs(matchScores, maxScores, labels, nMotifs=nMotifs)
    selMotifs <- unique(c(selMotifs, tfMotifCols))
  }
  else{
    matchScores <- assays(mae[["Motifs"]])$match_scores
    selMotifs <- unique(grep(paste(c(tfName, tfCofactors),collapse="|"),
                        colnames(matchScores), value=TRUE))
  }

  # Add CTCF-Features()
  if("CTCF_Signal" %in% features & tf!="CTCF"){
    message("CTCF Signal")

    # add the motif to selected motifs
    matchScores <- assays(mae[["Motifs"]])$match_scores
    motifCols <- grep("CTCF", colnames(matchScores), value=TRUE)
    selMotifs <- unique(selMotifs, motifCols)

    tfCols <- colnames(chIPMat)[grepl("CTCF", colnames(chIPMat))]
    sig <- sparseMatrixStats::rowMaxs(chIPMat[,tfCols, drop=FALSE])
    sig <- Matrix::Matrix(sig, ncol=1)

    colnames(sig) <- "CTCF_signal"
    ctcfFeat <- list(sig)
    names(ctcfFeat) <- "CTCF_signal"

    featMats <- append(featMats, ctcfFeat)
  }

  # add features full mae object
  featMats <- lapply(featMats, `colnames<-`, NULL)
  names(featMats) <- paste("tfFeat", names(featMats), sep="_")
  seTfFeat <- SummarizedExperiment(assays=featMats, rowRanges=coords)
  colnames(seTfFeat) <- tfName
  colData(seTfFeat)$feature_type <- "tfFeat"
  colData(seTfFeat)$tf_name <- tfName

  colsToMap <- getContexts(mae, tfName, which="ChIP")
  mae <- .addFeatures(mae, seTfFeat, colsToMap=colsToMap, prefix="tfFeat")

  # add cofactors for later use
  colDataTf <- colData(experiments(mae)$tfFeat)
  if(is.null(colDataTf$tf_cofactors)){
      colDataTf$tf_cofactors <- replicate(nrow(colDataTf), list)}

  colData(experiments(mae)$tfFeat)$tf_cofactors <-
    fifelse(colDataTf$tf_name==tfName, list(tfCofactors),
            colDataTf$tf_cofactors)

  # add associated motifs to colData
  if(is.null(colDataTf$preselected_motifs)){
      colDataTf$preselected_motifs <- replicate(nrow(colDataTf), list)}

  colData(experiments(mae)$tfFeat)$preselected_motifs <-
  fifelse(colDataTf$tf_name==tfName, list(selMotifs), colDataTf$preselected_motifs)

  return(mae)
}
