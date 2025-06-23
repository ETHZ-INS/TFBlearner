.binMat <- function(mat, threshold=NULL){
  cn <- colnames(mat)
  mat <- as(mat, "TsparseMatrix")

  if(is.null(threshold)) threshold <- 0L
  ind <- which(mat@x>threshold)
  mat <- sparseMatrix(i=mat@i[ind]+1L,
                      j=mat@j[ind]+1L,
                      x=rep(1L, length(ind)),
                      dims=c(nrow(mat), ncol(mat)))
  mat <- Matrix::Matrix(mat)
  colnames(mat) <- cn

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

   if(is.null(threshold)) threshold <- 0L

   toKeep <- which(chIPMat@x>threshold)
   colsInd <- chIPMat@j[toKeep]+1L
   chIPIndDt <- data.table(i=chIPMat@i[toKeep]+1L,
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
      chIPMat <- as(chIPMat, "CsparseMatrix")
      nmfRes <- suppressMessages(RcppML::nmf(chIPMat, k=nPatterns,
                                             L1=L1, seed=seed))
      fm <- nmfRes$w
      colnames(fm) <- paste(PATTFFEATNAME, 1:ncol(fm), sep="_")
    }
    else{
      aggMatTf <- .aggregate(chIPMat, aggVar="tf", aggFun=aggFun)
      aggMatTf <- as(aggMatTf, "CsparseMatrix")
      nmfTfRes <- suppressMessages(RcppML::nmf(aggMatTf, k=nPatterns, L1=L1,
                                               seed=seed))
      gc()

      aggMatCon <- .aggregate(chIPMat, aggVar="context", aggFun=aggFun)
      aggMatCon <- as(aggMatCon, "CsparseMatrix")
      nmfConRes <- suppressMessages(RcppML::nmf(aggMatCon, k=nPatterns, L1=L1,
                                                seed=seed))
      gc()

      wTf <- nmfTfRes$w
      colnames(wTf) <- paste(PATTFFEATNAME, 1:ncol(wTf), sep="_")
      wCon <- nmfConRes$w
      colnames(wCon) <- paste(PATCONTEXTFEATNAME, 1:ncol(wCon), sep="_")
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
  colnames(cScoreMat) <- CSCOREFEATNAME

  return(cScoreMat)
}

.selectMotifs <- function(matchScores,
                          maxScores,
                          labels,
                          addThr=4,
                          nMotifs=10,
                          subSample=10000)
{
  labels <- .binMat(labels, threshold=0L)
  labels <- .marginMax(labels, margin="row")
  thr <- maxScores/2

  subRows <- sample(1:nrow(matchScores), min(subSample, nrow(matchScores)))
  matchSubScores <- matchScores[subRows,,drop=FALSE]
  matchSubScores <- as(as.matrix(matchSubScores), "TsparseMatrix")
  labels <- labels[subRows]

  # top motif scores
  matchCoScores <- matchSubScores
  matchCoScores@x[matchCoScores@x < thr[matchCoScores@j + 1] &
                    matchCoScores@x<addThr*scalFactMotif] <- 0
  matchCoScores@x[matchCoScores@x >= thr[matchCoScores@j + 1] &
                    matchCoScores@x>=addThr*scalFactMotif] <- 1

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

  topExMotif <- intersect(topExMotif, colnames(matchScores))
  topCoMotif <- intersect(topCoMotif, colnames(matchScores))
  selectedMotifs <- c(topCoMotif, topExMotif)
  if(length(selectedMotifs)>0){
  names(selectedMotifs) <- c(paste(COMOTIFAFFIX, 1:length(topCoMotif), sep="_"),
                             paste(EXMOTIFAFFIX, 1:length(topExMotif), sep="_"))}

  # can happen if the motif-matches matrix has less columns than motifs to select
  selectedMotifs <- unique(selectedMotifs[!is.na(selectedMotifs)])

  return(selectedMotifs)
}

.getCofactorBindings <- function(chIPMat, tfCofactors){
  tfCols <- unlist(tstrsplit(colnames(chIPMat), split="_", keep=2))
  namesSub <- names(tfCofactors)[which(tfCofactors %in% tfCols)]
  tfCofactors <- intersect(tfCols, tfCofactors)
  names(tfCofactors) <- namesSub

  if(length(tfCofactors)>0){
    cofactBindings <- lapply(tfCofactors, function(tfCol){
      cofactBinding <- Matrix::Matrix(
        Matrix::rowMeans(chIPMat[,tfCols==tfCol, drop=FALSE]), ncol=1)
      colnames(cofactBinding) <- paste(COBINDFEATNAME, tfCol, sep=".")
      cofactBinding})
    names(cofactBindings) <- paste(COBINDFEATNAME,
                                   gsub(MOTIFAFFIX, "", namesSub), sep="_")
    return(cofactBindings)}
  else{
    return(NULL)
  }
}

# Adapted from: https://github.com/ETHZ-INS/DTFAB/blob/main/Scripts/GCnorm.R
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
    fn <- tcrossprod((1L-mat1),mat2) # 0-1 mismatches
    fp <- tcrossprod(mat1,(1L-mat2)) # 1-0 mismatches
    tn <- ncol(mat1)-(tp+fn+fp) # 0-0 matches
    #tn <- tcrossprod(abs(1-mat1),abs(1-mat2))
    if(!adjust) return((tp+tn)/ncol(mat1))
    Po <- (tp+tn)/ncol(mat1)
    Pe <- (((tp+fn)/ncol(mat1)) * ((tp+fp)/ncol(mat1))) +
      (((tn+fp)/ncol(mat1)) * ((tn+fn)/ncol(mat1)))
    return((Po-Pe)/(1-Pe))
}

.getAssociation <- function(atacMat1, atacMat2){

  # subset to common contexts
  commonContexts <- intersect(colnames(atacMat1), colnames(atacMat2))
  atacMat1 <- atacMat1[,commonContexts,drop=FALSE]
  atacMat2 <- atacMat2[,commonContexts,drop=FALSE]

  pearCor <- cor(Matrix::t(as.matrix(atacMat1)), Matrix::t(as.matrix(atacMat2)))
  colnames(pearCor) <- paste(PEARSONAFFIX, 1:ncol(pearCor), sep="_")
  pearCor <- Matrix::Matrix(pearCor)

  # max-scaling
  atacScalMat1 <- Matrix::t(Matrix::t(atacMat1)/.marginMax(atacMat1, margin="col"))
  atacScalMat2 <- Matrix::t(Matrix::t(atacMat2)/.marginMax(atacMat2, margin="col"))

  # binarizing
  atacBinMat1 <- .binMat(atacScalMat1, threshold=0.3)
  atacBinMat2 <- .binMat(atacScalMat2, threshold=0.3)

  cohKappa <- .matrixKappa(atacBinMat1, atacBinMat2)
  colnames(cohKappa) <- paste(COHENAFFIX, 1:ncol(cohKappa), sep="_")

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
#' "C_Score", "Cofactor_Binding", "CTCF_Signal", "Cooccuring_Motifs", "Associated_Motifs", "Associated_Motif_Activity".
#' Features are stored in the assays of the added experiment.
#' See [TFBlearner::listFeatures] for an overview of the features.
#' @param nPatterns Number of non-negative matrix factorization (NMF) components to consider for the decomposition of the ChIP-seq peaks matrix.
#' Passed to [RcppML::nmf].
#' @param L1 LASSO penalties for lower rank matrices w and h resulting from the NMF decomposition. Vector with two elements for w and h.
#' Passed to [RcppML::nmf]
#' @param nMotifs Number of associated motifs to select. Chooses `nMotifs` co-occuring and `nMotifs` mutually exclusive motifs.
#' @param annoCol Name of column indicating cellular contexts in colData.
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
#' @importFrom MatrixGenerics colMaxs rowMaxs
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @author Emanuel Sonder
#' @export
tfFeatures <- function(mae,
                       tfName,
                       tfCofactors=NULL,
                       features=c("Binding_Patterns", "Promoter_Association",
                                  "C_Score", "Cooccuring_Motifs",
                                  "Cofactor_Binding", "CTCF_Signal",
                                  "Associated_Motifs",
                                  "Associated_Motif_Activity"),
                       nPatterns=50,
                       L1=c(0.5,0.5),
                       nMotifs=10,
                       annoCol="context",
                       seed=42,
                       genome=BSgenome.Hsapiens.UCSC.hg38){
  set.seed(seed)

  allTfs <- unique(colData(mae[[CHIPEXP]])[[TFNAMECOL]])
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
                                            "Associated_Motifs",
                                            "Associated_Motif_Activity"),
                        several.ok=TRUE)
  featMats <- list()

  # validate MultiAssayExperiment object
  .checkObject(mae, checkFor=c("site", "context"))

  # reference coordinates
  coords <- rowRanges(mae[[MOTIFEXP]])
  if(length(tfCofactors)>0){
    names(tfCofactors) <- paste(TFCOFACTORMOTIFPREFIX,
                                1:length(tfCofactors),sep="_")}

  # assay-matrices
  atacMat <- .convertToMatrix(assays(mae[[ATACEXP]])[[TOTALOVERLAPSFEATNAME]])
  colnames(atacMat) <- colnames(mae[[ATACEXP]])
  whichCol <- which(mae[[CHIPEXP]][[TFNAMECOL]]!=tfName)
  chIPMat <- as(assays(mae[[CHIPEXP]])$peaks[,whichCol],"CsparseMatrix")
  colnames(chIPMat) <- paste(colData(mae[[CHIPEXP]])[whichCol,annoCol],
                             colData(mae[[CHIPEXP]])[whichCol,TFNAMECOL],
                             sep="_")

  # Normalize ATAC-seq data
  message("GC Normalization")
  gc <- assays(mae[[SITEFEAT]])[[paste(SITEFEAT,
                                       GCCONTFEATNAME, sep="_")]][,,drop=TRUE]

  atacMat <- .GCSmoothQuantile(gc, atacMat, nBins=20, round=TRUE)
  normedCountsName <- paste(TOTALOVERLAPSFEATNAME, GCNORMEDAFFIX, sep="_")
  assays(mae[[ATACEXP]], withDimnames=FALSE)[[normedCountsName]] <- atacMat

  if(ATACPROMEXP %in% names(experiments(mae))){
    # prune to standard chromosomes
    gc <- rowData(mae[[ATACPROMEXP]])[[GCCONTFEATNAME]]
    atacPromMat <- assays(mae[[ATACPROMEXP]])[[TOTALOVERLAPSFEATNAME]]
    atacPromMat <- .GCSmoothQuantile(gc, atacPromMat, nBins=20, round=TRUE)
    assays(mae[[ATACPROMEXP]])[[normedCountsName]] <- atacPromMat
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

    isProm <- which(rowData(mae[[ATACPROMEXP]])[[TFNAMECOL]]==tfName)
    atacPromMat <- atacPromMat[isProm,,drop=FALSE]

    promAsc <- .getAssociation(atacMat, atacPromMat)
    colnames(promAsc) <- paste(promoterPrefix, colnames(promAsc), sep="_")
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
    tfs <- c(tfName, tfCofactors)
    tfs <- intersect(colData(mae[[MOTIFEXP]])[[MOTIFNAMECOL]],tfs)

    coCounts <- lapply(tfs, function(tf){
      mmPath <- subset(colData(mae[[MOTIFEXP]]),
                       get(MOTIFNAMECOL)==tf)$origin
      baseDir <- metadata(colData(mae[[MOTIFEXP]]))[[BASEDIRCOL]]
      mmPath <- file.path(baseDir, mmPath)
      mm <- as.data.table(readRDS(mmPath))
      mm$motif <- tf

      if(COOCCURRENCECOL %in% colnames(mm)){
        aggregationFun <- sum
      }
      else{
        aggregationFun <- NULL
      }
      mmc <- genomicRangesMapping(coords, mm, byCols="motif",
                                  scoreCol=COOCCURRENCECOL,
                                  aggregationFun=aggregationFun)
    })

    # ensure correct naming
    namesCoCounts <- lapply(tfs, function(tf){
      if(tf==tfName){name <- paste(TFMOTIFPREFIX, 1:length(tf), sep="_")}
      else{name <- names(tfCofactors)[which(tfCofactors==tf)]}
    })
    names(coCounts) <- unlist(namesCoCounts)
    names(coCounts)  <- paste(COCOUNTFEATNAME, names(coCounts) , sep="_")
    featMats <- append(featMats, coCounts)
  }

  # Choose motif-matches by name
  motifNames <- colData(mae[[MOTIFEXP]])[[MOTIFNAMECOL]]
  tfSimMotifCols <- unique(grep(tfName, motifNames, value=TRUE))
  tfSimMotifCols <- setdiff(tfSimMotifCols, tfName)
  if(length(tfSimMotifCols)>0){
    names(tfSimMotifCols) <- paste(PRIORMOTIFPREFIX, 1:length(tfSimMotifCols),
                                sep="_")}

  tfCofactorCols <- unique(grep(paste(tfCofactors,collapse="|"),
                                motifNames, value=TRUE))
  if(length(tfCofactorCols)>0){
    names(tfCofactorCols) <- paste(TFCOFACTORMOTIFPREFIX,
                                   1:length(tfCofactorCols), sep="_")}
  tfMotifCols <- intersect(tfName, motifNames)
  if(length(tfMotifCols)>0){
    names(tfMotifCols) <- paste(TFMOTIFPREFIX, 1:length(tfMotifCols), sep="_")}
  priorMotifCols <- c(tfMotifCols, tfSimMotifCols, tfCofactorCols)

  # Choose motif-matches by association
  if("Associated_Motifs" %in% features){
    message("Select motifs")

    # get (except from test experiments) ChIP-labels of TF
    chIPCols <- colnames(mae[[CHIPEXP]])
    tfCols <- chIPCols[colData(mae[[CHIPEXP]])[[TFNAMECOL]]==tfName]
    isTesting <- subset(sampleMap(mae), assay==CHIPEXP & get(ISTESTCOL))$colname
    tfCols <- intersect(chIPCols[!(chIPCols %in% isTesting)], tfCols)

    if(length(tfCols)>0){
      labels <- assays(mae[[CHIPEXP]])[[PEAKASSAY]][,tfCols,drop=FALSE]
      labels <- .convertToMatrix(labels)
    }
    else{
      labels <- assays(mae[[MOTIFEXP]])[[MATCHASSAY]][,tfName,drop=FALSE]
      labels <- .convertToMatrix(labels)
    }

    # select motifs co-occuring around ChIP-peaks or motif matches of TF of interest
    matchScores <- assays(mae[[MOTIFEXP]])[[MATCHASSAY]]
    matchScores <- matchScores[,!c(colnames(matchScores) %in% priorMotifCols), drop=FALSE]
    colDataMotifs <- colData(mae[[MOTIFEXP]])
    colDataMotifs <- subset(colDataMotifs, get(MOTIFNAMECOL) %in% colnames(matchScores))
    colDataMotifs <- colDataMotifs[order(match(colDataMotifs[[MOTIFNAMECOL]],
                                               colnames(matchScores))),,
                                   drop=FALSE]

    maxScores <- colDataMotifs[[MAXSCORECOL]]
    selMotifs <- .selectMotifs(matchScores, maxScores, labels, nMotifs=nMotifs)
    if(length(selMotifs)>0){
      names(selMotifs) <- paste0(SELMOTIFPREFIX, names(selMotifs))}
  }
  else{
    selMotifs <- NULL
  }
  selMotifs <- c(selMotifs, priorMotifCols)

  # Choose activity-associated motifs by name
  actMotifNames <- colnames(mae[["Activity.Association"]])
  tfSimMotifCols <- unique(grep(tfName, actMotifNames, value=TRUE))
  tfSimMotifCols <- setdiff(tfSimMotifCols, tfName)
  if(length(tfSimMotifCols)>0){
    names(tfSimMotifCols) <- paste(PRIORMOTIFPREFIX, 1:length(tfSimMotifCols),
                                   sep="_")}

  tfCofactorCols <- unique(grep(paste(tfCofactors,collapse="|"),
                                actMotifNames, value=TRUE))
  if(length(tfCofactorCols)>0){
    names(tfCofactorCols) <- paste(TFCOFACTORMOTIFPREFIX,
                                   1:length(tfCofactorCols), sep="_")}
  tfMotifCols <- intersect(tfName, actMotifNames)
  if(length(tfMotifCols)>0){
    names(tfMotifCols) <- paste(TFMOTIFPREFIX, 1:length(tfMotifCols), sep="_")}
  priorMotifCols <- c(tfMotifCols, tfSimMotifCols, tfCofactorCols)

  # Choose activity-associated motifs by association
  if("Associated_Motif_Activity" %in% features){
    message("Select motifs with associated activity")

    # get (except from test experiments) ChIP-labels of TF
    chIPCols <- colnames(mae[[CHIPEXP]])
    tfCols <- chIPCols[colData(mae[[CHIPEXP]])[[TFNAMECOL]]==tfName]
    isTesting <- subset(sampleMap(mae), assay==CHIPEXP & get(ISTESTCOL))$colname
    tfCols <- intersect(chIPCols[!(chIPCols %in% isTesting)], tfCols)

    if(length(tfCols)>0){
      labels <- assays(mae[[CHIPEXP]])[[PEAKASSAY]][,tfCols,drop=FALSE]
      labels <- .convertToMatrix(labels)
    }
    else{
      labels <- assays(mae[[MOTIFEXP]])[[MATCHASSAY]][,tfName,drop=FALSE]
      labels <- .convertToMatrix(labels)
    }

    # select motifs co-occuring around ChIP-peaks or motif matches of TF of interest
    actAssoc <- assays(mae[[ASSOCEXP]])[[ASSOCASSAY]]
    actAssoc <- actAssoc[,!c(colnames(actAssoc) %in% priorMotifCols), drop=FALSE]
    selActMotifs <- .selectMotifs(actAssoc, rep(1, ncol(actAssoc)), labels,
                                  addThr=0, nMotifs=nMotifs)
    if(length(selActMotifs)>0){
      names(selActMotifs) <- paste0(SELMOTIFPREFIX, names(selActMotifs))}
  }
  else{
    selActMotifs <- NULL
  }
  selActMotifs <- c(selActMotifs, priorMotifCols)

  # Add CTCF-Features()
  if("CTCF_Signal" %in% features & tolower(tfName)!="ctcf"){
    message("CTCF Signal")

    # add the motif to selected motifs
    matchScores <- assays(mae[[MOTIFEXP]])[[MATCHASSAY]]
    motifAllCols <- grep("CTCF([:_]|$)", colnames(matchScores),
                      value=TRUE, ignore.case=TRUE)
    motifCols <- setdiff(motifAllCols, selMotifs)
    if(length(motifCols)>0){
      names(motifCols) <- paste(CTCFMOTIFPREFIX,1:length(motifCols), sep="_")
      selMotifs <- c(selMotifs, motifCols)}

    # add the motif to the selected activity motifs
    motifCols <- setdiff(motifAllCols, selActMotifs)
    if(length(motifCols)>0){
      names(motifCols) <- paste(CTCFMOTIFPREFIX, 1:length(motifCols), sep="_")
      selActMotifs <- c(selActMotifs, motifCols)}

    # add some checks if it is in data
    tfCols <- colnames(chIPMat)[grepl("_CTCF$", colnames(chIPMat),
                                      ignore.case=TRUE)]
    if(length(tfCols)>0){
      sig <- .marginMax(chIPMat[,tfCols, drop=FALSE], margin="row")
      sig <- Matrix::Matrix(sig, ncol=1)

      colnames(sig) <- CTCFFEATNAME
      ctcfFeat <- list(sig)
      names(ctcfFeat) <- CTCFFEATNAME
      featMats <- append(featMats, ctcfFeat)
    }
  }

  # add features full mae object
  featMats <- lapply(featMats, `colnames<-`, NULL)
  names(featMats) <- paste(TFFEAT, names(featMats), sep="_")
  seTfFeat <- SummarizedExperiment(assays=featMats, rowRanges=coords)
  colnames(seTfFeat) <- tfName
  colData(seTfFeat)[[FEATTYPECOL]] <- TFFEAT
  colData(seTfFeat)[[TFNAMECOL]] <- tfName

  colsToMap <- getContexts(mae, tfName, which="ChIP")
  mae <- .addFeatures(mae, seTfFeat, colsToMap=colsToMap, prefix=TFFEAT)

  # add cofactors for later use
  colDataTf <- colData(mae[[TFFEAT]])
  if(is.null(colDataTf[[TFCOFACTORSCOL]])){
    colDataTf[[TFCOFACTORSCOL]] <- rep(NA, nrow(colDataTf))}

  tfRow <- colDataTf[[TFNAMECOL]]==tfName
  colDataTf[tfRow,][[TFCOFACTORSCOL]] <- list(tfCofactors)

  # add associated motifs to colData
  if(is.null(colDataTf[[PRESELMOTIFCOL]])){
    colDataTf[[PRESELMOTIFCOL]] <- rep(NA, nrow(colDataTf))}

  if(!is.null(selMotifs)){
    colDataTf[tfRow,][[PRESELMOTIFCOL]] <- list(selMotifs)}

  # add motifs with associated activity to colData
  if(is.null(colDataTf[[PRESELACTCOL]])){
    colDataTf[[PRESELACTCOL]] <- rep(NA, nrow(colDataTf))}

  if(!is.null(selActMotifs)){
    colDataTf[tfRow,][[PRESELACTCOL]] <- list(selActMotifs)}

  colData(mae[[TFFEAT]]) <- colDataTf

  return(mae)
}
