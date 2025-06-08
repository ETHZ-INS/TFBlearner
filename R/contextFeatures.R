#' Wrapper around chromVAR
#'
#' Do not forget to register(MulticoreParam(whatever)) before!
#'
#' @param atacSe The path to the atac SummarizedExperiment.
#' If it already contains expectations and background peaks in the rowData / metadata, these will be reused for the computations.
#' @param motifSe The path to the H5 motif scores SummarizedExperiment
#' @param nSites The number of sites to subsample to
#' @param chunkSize The number of motifs to process in one CV call
#' @param seed The random seed
#' @param trimTop The number of top sites to discard
#' @param trimBottom The number of bottom sites to discard
#' @param threshold1 The first minimum score to try
#' @param minMatches The minimum number of matches below which we resort to the
#'   second threshold (half of the max score)
#' @param ... Arguments passed to [chromVAR::getBackgroundPeaks].
#'
#' @importFrom chromVAR getBackgroundPeaks computeDeviations computeExpectations
#' @returns A chromVARDeviations object
.CVwrapper <- function(atacSe, motifSe,
                       nSites=5e5L, chunkSize=20L, seed=1234,
                       trimTop=100L, trimBottom=2e5L, threshold1=1e5L, minMatches=2000L,
                       ...){

  register(SerialParam())
  se <- atacSe
  if(is.character(se)) se <- readRDS(se)

  # subsample
  atac <- as(assay(se, totalOverlapsFeatName), "sparseMatrix")
  set.seed(seed)

  exp <- rowData(atacSe)[[chromVarExpCol]]
  bg <-  metadata(atacSe)[[chromVarBgCol]]

  if(is.null(exp) | is.null(bg)){

    # remove extreme sites
    if(trimBottom+trimTop>nrow(atac)){
      trimBottom <- floor(0.1*nrow(atac))
      trimTop <- floor(0.1*nrow(atac))
    }
    toRemove <- c(1:trimTop, (nrow(atac)-trimBottom):nrow(atac))
    nSamp <- min(nSites, nrow(atacSe)-length(toRemove))
    idx <- sample(order(rowMeans(atac),decreasing=TRUE)[-toRemove], nSamp)
    atac <- SummarizedExperiment(list(counts=atac[idx,]), colData=colData(se),
                                 rowRanges=rowRanges(se)[idx])

    # remove empty rows
    w <- which(rowSums(assay(atac))>0L)
    idx <- idx[w]
    atac <- atac[w,]
    bg <- getBackgroundPeaks(atac, bias=rowData(atac)[[gcContFeatName]], ...)
    exp <- computeExpectations(atac)}
  else{
    # in case bg atac & subInd have been computed beforehand
    idx <- which(!is.na(exp))
    exp <- exp[idx]
    bg <- bg[idx]
    atac <- atac[idx,]
  }

  se <- motifSe
  chunkSize <- min(chunkSize, floor(ncol(se)/2))
  if(is.character(se)) se <- readRDS(se)
  se <- se[idx,]
  chunks <- split(seq_len(ncol(se)),
                  cut(seq_len(ncol(se)), breaks=floor(ncol(se)/chunkSize)))
  dev <- lapply(chunks, FUN=function(i){
    tmp <- as.matrix(assay(se)[,i])
    tmp <- as(sapply(seq_len(ncol(tmp)), \(j){
      x2 <- tmp[,j]>threshold1
      if(sum(x2)>=minMatches) return(x2)
      tmp[,j]>as.integer(round(se[[maxScoreCol]][i[j]]/2L))
    }), "sparseMatrix")
    colnames(tmp) <- colnames(se)[i]
    computeDeviations(atac, tmp, background_peaks=bg, expectation=exp)
  })
  dev <- do.call(rbind, dev)
  colData(dev) <- colData(atac)
  rowData(dev)[[maxScoreCol]] <- se[[maxScoreCol]]
  dev <- dev[rowSums(is.na(assay(dev)))<ncol(dev),]
  assay(dev,devAssay)[is.na(assay(dev, devAssay))] <- 0
  assay(dev, normDevAssay) <- scale(assay(dev, devAssay))
  return(list(dev=dev, exp=exp, bg=bg, idx=idx))
}


#' CV2localAssociation
#'
#' Computes the pearson correlation between global TF activity (chromVAR scores)
#' and local accessibility. Stored are pearson correlations times 1000, rounded.
#'
#' @param cvSe The ChromVAR deviations SE, containing a `r normDevAssay` assay.
#' @param atacSe The ATAC SE, containing a `r totalOverlapsFeatName` assay.
#' @param h5path The path to the output H5 file.
#' @param sparTh Absolute sparsification threshold, default 200L (-0.2 to 0.2
#'   pearson correlations are set to zero).
#'
#' @importFrom rhdf5 h5createFile h5createDataset h5write h5closeAll
#' @returns Nothing
.CV2localAssociation <- function(cvSe, atacSe, outDir, sparTh=200L,
                                 saveHdf5=TRUE,
                                 BPPARAM=SerialParam()){
  if(is.character(cvSe)) cvSe <- readRDS(cvSe)
  if(is.character(atacSe)) atacSe <- readRDS(atacSe)
  atacMat <- TFBlearner:::.GCSmoothQuantile(
    counts=as.matrix(assay(atacSe, totalOverlapsFeatName)),
    gc=rowData(atacSe)[[gcContFeatName]])
  atacMat <- t(as.matrix(atacMat))

  threads <- floor(getDTthreads()/BPPARAM$workers)
  ascStats <- BiocParallel::bplapply(rownames(cvSe),
                                     function(motif, cvSe, atacMat,
                                              sparTh, threads,
                                              saveHdf5, outDir){

                                       data.table::setDTthreads(threads)
                                       x <- as.integer(round(1000*cor(atacMat, assay(cvSe, normDevAssay)[motif,])))
                                       q <- as.integer(round(quantile(x, prob=c(0,0.1,0.2,0.8,0.9,1), na.rm=TRUE)))
                                       x[abs(x)<sparTh] <- 0L
                                       asc <- Matrix::Matrix(x, ncol=1)
                                       colnames(asc) <- motif

                                       if(saveHdf5){
                                         dataList <- list(asc)
                                         names(dataList) <- assocAssay
                                         .writeToHdf5(dataList,
                                                      paste0(file.path(outDir, motif), ".h5"),
                                                      storage="integer")
                                         asc <- head(asc, n=1)}
                                       return(asc)},
                                     cvSe=cvSe, atacMat=atacMat,
                                     sparTh=sparTh,
                                     thread=threads,
                                     saveHdf5=saveHdf5, outDir=outDir,
                                     BPPARAM=BPPARAM)
  rm(atacMat)

  if(saveHdf5){
    fileName <- paste(assocExp, "mapped", sep="_")
    hdf5FilesCollect <- unlist(lapply(rownames(cvSe), function(f) paste0(file.path(outDir,f), ".h5")))
    ascStats <- .collectHdf5Files(hdf5FilesCollect,
                                  hdf5FileName=paste0(file.path(outDir, fileName), ".h5"),
                                  storage="integer", asSparse=FALSE)
    ascStats <- ascStats[[1]]
  }
  else{
    ascStats <- Reduce("cbind", ascStats[-1], ascStats[[1]])
  }
  colnames(ascStats) <- rownames(cvSe)

  assaysAsc <- list(ascStats)
  names(assaysAsc) <- assocAssay
  ascSe <- SummarizedExperiment(assays=assaysAsc, rowRanges=rowRanges(atacSe))

  return(ascSe)
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

  rowVars <- apply(atacNormMat, 1, var)
  varMat <- Matrix::Matrix(rowVars, ncol=1)
  colnames(varMat) <- atacVarFeatName
  return(varMat)
}

#' Features computed across cellular contexts
#'
#' Either site or cellular-context-specific features which are computed based on ATAC-seq data of all contexts combined.
#' Examples are site-specific variance in the ATAC-signal across contexts or a MDS projection of the cellular contexts.
#'
#' @name panContextFeatures
#' @param mae [MultiAssayExperiment::MultiAssayExperiment-class] as construced by [TFBlearner::prepData()] containing Motif, ATAC-, ChIP-seq,
#' site-specific features as obtained by [TFBlearner::siteFeatures()].
#' @param nVarSites Number of sites with highest ATAC-signal variance to include for the MDS projections.
#' @param features Names of features to be added. Can be all or some of "MDS_Context", "Max_ATAC_Signal", "ATAC_Variance".
#' ChromVAR activity scores and associations to the ATAC-signal will be computed as well if not already present in the experiments of the object.
#' @param annoCol Name of column indicating cellular contexts in colData.
#' @param seed Integer value for setting the seed for random number generation with [base::set.seed].
#' @param saveHdf5 If chromVar activity scors and associations experiments should be saved as HDF5 files.
#' @param outDir Directory to save HDF5 file to.
#' @param BPPARAM Parallel back-end to be used. Passed to [BiocParallel::bplapply()].
#' @param ... Arguments passed to [chromVAR::getBackgroundPeaks].
#' @return [MultiAssayExperiment::MultiAssayExperiment-class] object with columns added to rowData of the `r atacExp` experiment ("ATAC_Variance", "Max_ATAC_Signal")
#' or the column data of the `r atacExp`  experiment ("MDS_Context"). If not present in the object experiments for the chromVAR activity scores (name: `r actExp`)
#' and its association to the ATAC-signal of a site across contexts (name: `r assocExp`) are added.
#' @import MultiAssayExperiment
#' @importFrom BiocParallel bplapply SerialParam MulticoreParam SnowParam register
#' @importFrom chromVAR getBackgroundPeaks computeDeviations computeExpectations
#' @importFrom rhdf5 H5Fcreate h5createDataset h5writeDataset h5delete h5write h5ls H5Fopen H5Fclose H5garbage_collect
#' @importClassesFrom HDF5Array HDF5Array
#' @importClassesFrom SummarizedExperiment SummarizedExperiment RangedSummarizedExperiment
#' @importFrom stats cmdscale
#' @export
panContextFeatures <- function(mae,
                               nVarSites=1e5,
                               features=c("MDS_Context",
                                          "Max_ATAC_Signal",
                                          "ATAC_Variance"),
                               annoCol="context",
                               seed=42,
                               saveHdf5=FALSE,
                               outDir=getwd(),
                               BPPARAM=SerialParam(), ...){

  .checkObject(mae, checkFor="Site")

  threads <- floor(getDTthreads())/BPPARAM$workers
  data.table::setDTthreads(threads)

  features <- match.arg(features, choices=c("MDS_Context",
                                            "Max_ATAC_Signal",
                                            "ATAC_Variance"),
                        several.ok=TRUE)

  refCoords <- rowRanges(mae[[atacExp]])
  gcContent <- as.numeric(assays(mae[[siteFeat]])[[paste(siteFeat, gcContFeatName, sep="_")]])
  rowData(mae[[atacExp]])[[gcContFeatName]] <- gcContent

  atacMat <- .convertToMatrix(assays(mae[[atacExp]])[[totalOverlapsFeatName]])
  nVarSites <- min(nVarSites, nrow(atacMat))

  if("ATAC_Variance" %in% features){
    message("Get site-specific variance of ATAC-signal")

    varFeat <- .getVariableSites(atacMat)
    rs <- .marginSum(atacMat, margin="row")
    topVarSites <- setdiff(order(-varFeat[,1])[1:nVarSites], which(rs==0))

    rowData(mae[[atacExp]])[[topVarSitesCol]] <- fifelse(1:nrow(atacMat) %in% topVarSites,
                                                          TRUE, FALSE)
    rowData(mae[[atacExp]])[[atacVarFeatName]] <- varFeat[,1,drop=TRUE]
    subInd <- topVarSites
  }
  else{
    subInd <- sample(1:length(refCoords), nVarSites)
  }

  atacMatSub <- atacMat[subInd,]

  if("MDS_Context" %in% features){
    message("Get MDS-Dimensions")
    if(ncol(atacMatSub)>2){
    mdsDim <- .getContextProjection(atacMatSub, subSample=NULL)
    mdsDim <- as.data.table(mdsDim, keep.rownames=TRUE)

    # add to colData of ATAC
    co <- match(mdsDim$rn, colData(mae[[atacExp]])[[annoCol]])
    colData(mae[[atacExp]]) <-  cbind(colData(mae[[atacExp]]),
                                      mdsDim[co,])
    }
    else{
      message("Too view cellular contexts to compute MDS-dimensions. Skipping...")
    }
  }

  if("Max_ATAC_Signal" %in% features){
    message("Get maximal ATAC-signal per site")

    atacMat <- .minMaxNormalization(atacMat)
    maxFeat <- as.matrix(.marginMax(atacMat, margin="row"), ncol=1)
    colnames(maxFeat) <- maxAtacFeatName

    # add to rowData of ATAC
    rowData(mae[[atacExp]]) <- cbind(rowData(mae[[atacExp]]),
                                     as.data.table(maxFeat))
  }

  if(!all(c(actExp, assocExp) %in% names(mae))){
    message("Get chromVAR activity estimates across cellular contexts")
    res <- .CVwrapper(atacSe=mae[[atacExp]], motifSe=mae[[motifExp]],
                      seed=seed, ...)
    dev <- res$dev
    gc()

    # add as experiment
    actSe <- SummarizedExperiment(assays=assays(dev))
    mae <- .addFeatures(mae, actSe, colsToMap=colnames(actSe), prefix=actExp)

    # add expectations to rowdata in case context is added
    rowData(mae[[atacExp]])[res$id, chromVarExpCol] <- res$exp
    metadata(mae[[atacExp]])[[chromVarBgCol]] <- res$bg

    message("Get association between site-specific ATAC-signal and chromVAR activity estimates")
    ascSe <- .CV2localAssociation(actSe, mae[[atacExp]],
                                  saveHdf5=saveHdf5, outDir=outDir,
                                  BPPARAM=BPPARAM)
    mae <- .addFeatures(mae, ascSe, colsToMap=colnames(actSe), prefix=assocExp)
  }
  else{
    message(paste0("chromVAR activity estimates and associations to site-specific
                   ATAC-signal have been pre-computed and will reused in the following.
                   In case these should be re-computed please set mae[[",actExp,"]] <- NULL and ",
                   "mae[[", assocExp, "]] <- NULL"))
  }

  return(mae)
}
