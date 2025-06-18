.checkObject <- function(mae,
                         checkFor=c("none", "site", "context",
                                    "tf", "tf-context"),
                         tfName=NULL){

  if("none" %in% checkFor) checkFor <- "none"
  checkFor <- match.arg(checkFor, choices=c("none", "site",
                                            "context", "tf", "tf-context"),
                        several.ok=TRUE)

  if(!(is(mae, "MultiAssayExperiment"))){
    stop("Object provided needs to be a MultiAssayExperiment Object.
          For obtaining a object in the correct form use prepData()")
  }

  expectedExperiments <- c(atacExp, chIPExp, motifExp)
  experimentNames <- names(experiments(mae))
  if(!all(expectedExperiments %in% experimentNames)){
    stop("Object provided needs to be a MultiAssayExperiment Object containing
          experiments: ATAC, ChIP and Motifs.
          For obtaining a object in the correct form use prepData()")
  }

  # check that motifs match to all contexts
  contexts <- unique(sampleMap(mae)$primary)
  sampleMapMotif <- as.data.table(subset(sampleMap(mae), assay==motifExp))
  contextMotifMap <- sampleMapMotif[,.(n_context=length(unique(primary))),
                                       by=colname]
  if(nrow(subset(contextMotifMap, n_context!=length(contexts)))>0){
    stop("Motifs need to map to all primaries (cellular context) labels.
          For obtaining a object in the correct form use prepData()")
  }

  if("site" %in% checkFor){
    if(!(siteFeat %in% experimentNames)){
      stop("Site-specific features as obtained by using siteFeatures() on the object
            are required for further use.")
    }
  }

  if("context" %in% checkFor){
    if(!all(c(actExp, assocExp) %in% experimentNames)){
      stop("Context-specific & pan-context features as obtained by using panContextFeatures() on the object
            are required for further use.")
    }
  }

  if("tf" %in% checkFor){
    if(is.null(tfName)){
      warning("Please provide the name of a TF for which to check if features are present")}
    if(!(tfFeat %in% experimentNames) || !(tfName %in% colnames(mae[[tfFeat]]))){
      stop(paste0("TF-features as obtained by using tfFeatures(tfName=",tfName,",...) are required for further use."))
    }
  }

  if("tf-context" %in% checkFor){
    if(is.null(tfName)){
      warning("Please provide the name of a TF for which to check if features are present")}
    if(!(contextTfFeat %in% experimentNames)){
      stop("Context-features as obtained by using contextTffeatures() on the object
            are required for further use.")
    }
    else{
      tfs <- colData(mae[[contextTfFeat]])[[tfNameCol]]
      if(!(tfName %in% tfs)){
        stop(paste0("Context-features as obtained by using contextTffeatures(tfName=",tfName,",...), are required for further use."))
      }

    }

  }

  #TODO: Add further checks, current form merely captures the idea.
  # - that origin path is available in colData of ATAC-seq !!
  #message("The provided MultiAssayExperiment object is of correct form")

  return(TRUE)
}

.seToMae <- function(mae, seFeat, prefix, colsToMap, annoCol="context"){

  colFeatData <- colData(seFeat)
  if(length(setdiff(rownames(colFeatData), colsToMap))!=0){
    mapFeat <- data.frame(primary=rep(colsToMap, nrow(colFeatData)),
                          colname=rep(rownames(colFeatData),
                                    each=length(colsToMap)),
                          stringsAsFactors=FALSE)}
  else{
    mapFeat <- data.frame(primary=colsToMap,
                          colname=colsToMap,
                          stringsAsFactors=FALSE)
  }
  #contextMap <- data.frame(row.names=colsToMap,
  #                         description=rep(annoCol, length(colsToMap)))

  if(prefix %in% names(experiments(mae))){
    seOrig <- experiments(mae)[[prefix]]
    allAssays <- unique(c(names(assays(seOrig)), names(assays(seFeat))))
    missingNew <- setdiff(allAssays, names(assays(seFeat)))
    missingOrig <- setdiff(allAssays, names(assays(seOrig)))

    for(missing in missingOrig){
      mat <- Matrix::Matrix(nrow=nrow(seOrig), ncol=ncol(seOrig))
      colnames(mat) <- colnames(seOrig)
      assays(seOrig)[[missing]] <- mat
    }

    for(missing in missingNew){
      mat <- Matrix::Matrix(nrow=nrow(seFeat), ncol=ncol(seFeat))
      colnames(mat) <- colnames(seFeat)
      assays(seFeat)[[missing]] <- mat
    }

    # fill missing colData columns
    for(col in setdiff(colnames(colData(seOrig)), colnames(colData(seFeat)))){
      colData(seFeat)[[col]] <- rep(NA, nrow(colData(seFeat)))}

    # Avoid duplicate columns - keep newly added columns
    isDelayed <- unlist(lapply(assays(seOrig), is, "DelayedMatrix"))
    isDelayed <- fifelse(sum(isDelayed)>0, TRUE, FALSE)
    seFeat <- combineCols(
      seOrig[,!(colnames(seOrig) %in% colnames(seFeat))], seFeat,
      use.names=FALSE, delayed=isDelayed)

    mapOrig <- as.data.table(subset(sampleMap(mae), assay==prefix)[,c("primary",
                                                                      "colname")])
    mapFeat <- unique(rbind(mapOrig, mapFeat, use.names=TRUE, fill=TRUE))
    mapFeat <- as.data.frame(mapFeat)
  }

  contextMap <- data.frame(row.names=unique(mapFeat$primary))
  contextMap[[annoCol]] <- unique(mapFeat$primary)
  #contextMap$description <- annoCol

  listMap <- list(mapFeat)
  names(listMap) <- prefix
  dfMap <- listToMap(listMap)
  objList <- list(seFeat)
  names(objList) <- prefix

  # construct MultiAssayExperimentObject containing features
  maeFeat <- MultiAssayExperiment(experiments=objList,
                                  colData=contextMap,
                                  sampleMap=dfMap)
  mapFeat <- as.data.table(sampleMap(maeFeat))
  mapMae <- subset(as.data.table(sampleMap(mae)), assay==atacExp)

  mapFeat <- merge(mapFeat, mapMae[,c("primary", isTestCol, isTrainCol), with=FALSE],
                   by.x=c("primary"),
                   by.y=c("primary"), all.x=TRUE, all.y=FALSE,
                   allow.cartesian=TRUE)
  mapFeat[,eval(isTestCol):=fifelse(is.na(get(isTestCol)),FALSE,get(isTestCol))]
  mapFeat[,eval(isTrainCol):=fifelse(is.na(get(isTrainCol)),FALSE,get(isTrainCol))]

  sampleMap(maeFeat) <- mapFeat

  return(maeFeat)
}

.retainMappings <- function(mae){
  experimentNames <- intersect(names(experiments(mae)),
                               c(motifExp, siteFeat, tfFeat))
  for(experimentName in experimentNames){
    mapComb <- as.data.table(sampleMap(mae))
    mapMotif <- subset(mapComb, assay==experimentName)
    allComb <- data.table(expand.grid(primary=unique(mapComb$primary),
                                      colname=unique(mapMotif$colname)))
    allComb <- allComb[,lapply(.SD, as.character), .SDcols=colnames(allComb)]
    allComb$assay <- experimentName
    allComb[[isTrainCol]] <- FALSE
    allComb[[isTestCol]] <- FALSE
    sampleMap(mae) <- data.frame(rbind(mapComb,
                                       allComb[!mapMotif, on=.(primary, colname)],
                                       use.names=TRUE))
  }
  return(mae)
}

.addFeatures <- function(mae, seFeat,
                         colsToMap,
                         prefix="features",
                         annoCol="context",
                         ...){
  .checkObject(mae)
  maeFeat <- .seToMae(mae, seFeat=seFeat, prefix=prefix,
                      colsToMap=colsToMap, annoCol=annoCol)
  experiments(mae)[[prefix]] <- NULL

  # add missing cols in old somewhere here (with NAs)

  # ensure that motifs still map to every cellular context
  mae <- c(mae, maeFeat)
  mae <- .retainMappings(mae)

  return(mae)
}

#' Get cellular contexts
#'
#' Provides the names of cellular contexts for data of interest, e.g. a TF or modality.
#'
#' @param mae [MultiAssayExperiment::MultiAssayExperiment-class] as construced by [TFBlearner::prepData()] containing Motif, ATAC- and ChIP-seq data.
#' @param tfName Name of transcription factor to get cellular contexts for.
#' If `tfName=NULL` all celullar contexts of the requested modalities (`which`), will be provided.
#' @param which Which modality to get the cellular contexts names from.
#' If `which=both` only cellular contexts which have both ATAC and ChIP-seq data (for the specified TF, in case `tfName` is not NULL) will be returned.
#' @export
getContexts <- function(mae, tfName=NULL, which=c("ChIP", "ATAC", "Both")){
  .checkObject(mae)

  which <- match.arg(which, choices=c("ChIP", "Both", "ATAC"))
  contextsAtac <- rownames(colData(mae[[atacExp]]))

  if(which=="ATAC"){
    contexts <- contextsAtac
  }
  else{
    colsChIP <- colnames(mae[[chIPExp]])
    colsChIP <- as.data.table(tstrsplit(colsChIP, split="_"))
    if(!is.null(tfName)){
      contextsChIP <- subset(colsChIP, V2==tfName)$V1}
    else{
      contextsChIP <- colsChIP$V1}

    if(which=="Both"){
      contexts <- intersect(contextsChIP,contextsAtac)}
    else{
      contexts <- contextsChIP
    }
  }

  return(contexts)
}

#' Example Data: Coordinates
#'
#' Example genomic coordinates to construct a custom MultiAssayExperiment object and to compute features for TF-Binding predictions for.
#' The sample data contains 100 genomic ranges stored in GRanges object.
#'
#' @name example_coords
#' @format [GenomicRanges::GRanges] object
#' @docType data
#' @keywords GRanges
NULL

#' Example Data: Motifs
#'
#' Example motifs of JUN and CTCF as PWMatrixList. Can be used as an example input for prepMotifs().
#'
#' @name example_pwms
#' @format [TFBSTools::PWMatrixList] object
#' @docType data
#' @keywords PWMatrixList
NULL

#' Add further ATAC-seq datasets to object
#'
#' Adds additional columns for provided ATAC-seq datasets in the `r atacExp` experiment and
#' computes features if requested.
#'
#' @name addATACData
#' @param mae [MultiAssayExperiment::MultiAssayExperiment-class] as construced by [TFBlearner::prepData()] containing Motif, ATAC-, ChIP-seq.
#' @param atacData Named list of [data.table::data.table]/data.frame/[GenomicRanges::GRanges]
#' or paths to .bam/.bed files containing ATAC-seq fragments, names being the cellular contexts labels.
#' Need to contain genomic coordinates (e.g. a chr/seqnames, start and end column).
#' @param testSet Vector of cellular context labels used for testing.
#' @param computeFeatures If and how features for the newly added ATAC-seq dataset should be computed.
#'
#' Option `simple`: Features for newly added ATAC-seq dataset will be computed based on pre-computed chromVAR expectations & backgrounds in `rowData` of the `r atacExp` experiment.
#'
#' Option `scratch`: features will be re-computed from scratch including the added ATAC-seq dataset.
#' ChromVAR expectations and backgrounds will be recomputed by re-running [TFBlearner::panContextFeatures()].
#'
#' Option `none`: No features will be recomputed for the added cellular context.
#'
#' For both `simple` & `scratch`, context-tf-features for the newly added context will be computed for the same TFs as for the cellular contexts already contained in the object.
#' @param annoCol Name of column indicating cellular contexts in colData.
#' @param shift Only for ATAC-seq data, if Tn5 insertion bias should be considered (only if a strand column is provided).
#' @param seed Integer value for setting the seed for random number generation with [base::set.seed].
#' @param BPPARAM Parallel back-end to be used. Passed to [BiocParallel::bplapply()].
#' @param ... Additional arguments passed to [TFBlearner::contextTfFeatures()] and [TFBlearner::panContextFeatures()] (in case `computeFeatures="scratch"`).
#' @return [MultiAssayExperiment::MultiAssayExperiment-class] with added ATAC-seq experiments (and cellular context and TF-specific features if `computeFeatures=TRUE`).
#' @import MultiAssayExperiment
#' @importFrom BiocParallel bplapply SerialParam MulticoreParam SnowParam
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges colData assays
#' @importFrom S4Vectors combineCols
#' @importFrom GenomeInfoDb seqlevelsStyle
addATACData <- function(mae, atacData,
                        testSet=NULL,
                        computeFeatures=c("simple", "scratch", "none"),
                        annoCol="context",
                        shift=FALSE,
                        seed=42,
                        BPPARAM=SerialParam(), ...){
  set.seed(seed)
  .checkObject(mae)
  computeFeatures <- match.arg(computeFeatures)

  # get precomputed stats
  mdsDimFeats <- paste(mdsDimFeatName, 1:2, sep="_")
  compMDS <- all(mdsDimFeats %in% colnames(colData(mae[[atacExp]])))
  if(compMDS){
    mdsDims <- colData(mae[[atacExp]])[,mdsDimFeats]
  }
  exp <- rowData(mae[[atacExp]])[[chromVarExpCol]]
  bg <- metadata(mae[[atacExp]])[[chromVarBgCol]]

  coords <- rowRanges(mae[[motifExp]])
  seAtac <- .mapSeqData(atacData, coords, type="ATAC", annoCol=annoCol,
                        shift=shift, BPPARAM=BPPARAM)

  mae <- .addFeatures(mae, seAtac, names(atacData), prefix=atacExp,
                      annoCol=annoCol)
  matchedContexts <- intersect(colData(mae[[atacExp]])[[annoCol]],
                               colData(mae[[chIPExp]])[[annoCol]])

  annoDt <- as.data.table(sampleMap(mae))
  annoDt[,eval(isTestCol):=fifelse(primary %in% testSet, TRUE, get(isTestCol))]
  annoDt[,eval(isTestCol):=fifelse(is.na(get(isTestCol)), FALSE, get(isTestCol))]
  annoDt[,eval(isTrainCol):=fifelse(primary %in% matchedContexts &
                                    !get(isTrainCol), TRUE, FALSE)]
  sampleMap(mae) <- annoDt

  colData(mae)[[isTestCol]] <- fifelse(colData(mae)[[annoCol]] %in% testSet,
                                     TRUE, colData(mae)[[isTestCol]])
  colData(mae)[[isTestCol]] <- fifelse(is.na(colData(mae)[[isTestCol]]), FALSE,
                                       colData(mae)[[isTestCol]])
  colData(mae)[[isTrainCol]] <- fifelse(colData(mae)[[annoCol]] %in% matchedContexts &
                                       !colData(mae)[[isTestCol]],TRUE, FALSE)

  if(computeFeatures!="none"){
    if(contextTfFeat %in% names(experiments(mae))){
      tfs <- unique(colData(mae[[contextTfFeat]])[[tfNameCol]])
    }

    # get feature to compute
    lDt <- listFeatures()
    lPanDt <- subset(lDt, feature_type=="pan-context-Feature")
    lPanFeatures <- unlist(tstrsplit(unlist(lPanDt$feature_matrix_column_names),
                                     split="_", keep=2))
    preFeats <- colnames(rowData(mae[[atacExp]]))
    panConFeats <- lPanDt[lPanFeatures %in% preFeats, ]$feature_name
    if(compMDS){panConFeats <- c(panConFeats, "MDS_Context")}

    panConFeats <- intersect(panConFeats,
                             eval(formals(TFBlearner::panContextFeatures)$features))

    if(computeFeatures=="scratch" | !(actExp %in% names(experiments(mae)))){
       experiments(mae)[[actExp]] <- NULL
       experiments(mae)[[assocExp]] <- NULL

       extraArgs <- list(...)
       extraArgs$features <- NULL
       mae <- do.call(TFBlearner::panContextFeatures,
                      c(list(mae=mae,
                             seed=seed,
                             features = panConFeats), extraArgs))
    }
    else{
      if(compMDS){
        atacMat <- .convertToMatrix(assays(mae[[atacExp]])[[totalOverlapsFeatName]])
        idx <- rowData(mae[[atacExp]])[[mdsSubRowCol]]
        atacMat <- atacMat[idx,]
        atacOrigMat <- atacMat[,setdiff(colnames(atacMat),
                                        names(atacData)), drop=FALSE]
        mdsOrig <- metadata(mae)[[mdsDimStatsEntry]]
        mdsDimsNew <- mdsDims
        for(context in names(atacData)){
          mdsRow <- .projectOnMDS(mdsOrig, t(atacOrigMat),
                                  atacMat[,context,drop=TRUE])
          mdsRow <- t(data.frame(mdsRow))
          colnames(mdsRow) <- mdsDimFeats
          rownames(mdsRow) <- context
          mdsDimsNew <- rbind(mdsDimsNew, mdsRow)
        }
        co <- match(rownames(mdsDimsNew), colData(mae[[atacExp]])[[annoCol]])
        colData(mae[[atacExp]])[,mdsDimFeats] <- mdsDimsNew[co,]
      }

      rowData(seAtac)[[chromVarExpCol]] <- exp
      metadata(seAtac)[[chromVarBgCol]] <- bg
      res <- .CVwrapper(atacSe=seAtac, motifSe=mae[[motifExp]],
                        seed=seed, ...)
      dev <- res$dev

      actSe <- SummarizedExperiment(assays=assays(dev))
      mae <- .addFeatures(mae, actSe, colsToMap=colnames(actSe), prefix=actExp)
    }

    if(contextTfFeat %in% names(experiments(mae))){

      # get which features to recompute
      lconTfDt <- subset(lDt, feature_type=="context-tf-Feature")
      lconTfDt[,feat_name_affix:=tstrsplit(feature_matrix_column_names,
                                           split="_", keep=2)]
      preFeats <- unlist(tstrsplit(names(assays(mae[[contextTfFeat]])),
                                            split="_", keep=2))
      conTfFeats <- subset(lconTfDt, feat_name_affix %in% preFeats)$feature_name
      conTfFeats <- intersect(conTfFeats,
                              eval(formals(TFBlearner::contextTfFeatures)$features))

      extraArgs <- list(...)
      extraArgs$features <- NULL
      for(context in names(atacData)){
        for(tf in tfs){
          labelledContexts <- getContexts(mae, tfName=tf, which="Both")
          addLabels <- fifelse(context %in% labelledContexts, TRUE, FALSE)
          mae <- do.call(TFBlearner::contextTfFeatures,
                         c(list(mae=mae,
                                tfName=tf,
                                whichCol="Col",
                                colSel=context,
                                features=conTfFeats,
                                seed=seed,
                                addLabels=addLabels), extraArgs))}
      }
    }
  }
  return(mae)
}

#' Custom MultiAssayExperiment construction.
#'
#' Creates [MultiAssayExperiment::MultiAssayExperiment-class ] object for downstream feature construction.
#'
#' @name prepData
#' @param refCoords Coordinates across which to aggregate ATAC-, ChIP-seq and motif data.
#' Constitutes the rowRanges of the RangedSummarizedExperiment object.
#' @param motifData  Named list of [data.table::data.table]/data.frame/[GenomicRanges::GRanges] containing motif matches.
#' Ideally this argument is obtained by calling [TFBlearner::prepMotifs].
#' Needs to contain genomic coordinates (e.g. a chr/seqnames, start and end column), a `r scoreCol` column (motif matching score).
#' Ideally motifs corresponding to TFs for which `chIPData` is present are named the same in both arguments
#' (e.g `motifData=list(JUN=...)`, `chIPData=list(K562_JUN=...)`.
#' @param atacData Named list of [data.table::data.table]/data.frame/[GenomicRanges::GRanges]
#' or paths to .bam/.bed files containing ATAC-seq fragments, names being the cellular contexts labels.
#' Need to contain genomic coordinates (e.g. a chr/seqnames, start and end column).
#' @param chIPData Named list of [data.table::data.table]/data.frame/[GenomicRanges::GRanges]
#' or paths to .bam/.bed files containing ChIP-seq peaks, names being the combination of cellular context and TF labels (e.g. `K562_YY1`).
#' Need to contain genomic coordinates (e.g. a chr/seqnames, start and end column) and optionally a column with observational weights (`weightCol`)
#' and/or an uncertainty label (`isUncertainCol`).
#' @param testSet Vector of cellular context labels used for testing.
#' @param promoterCoords Optionally coordinates of promoters of transcription factor genes. Needs to contain a column in the metadata named `r tfNameCol`.
#' @param aggregationFun Function to aggregate
#' @param scoreCol Optional, column name of motif-matching data containing matching scores.
#' @param weightCol Optional, column name of ChIP-seq data containing observational weights for peaks.
#' @param isUncertainCol Optional, column name of ChIP-seq data labelling uncertain peaks (TRUE/FALSE).
#' @param annoCol Name of column indicating cellular contexts in colData.
#' @param shift Only for ATAC-seq data, if Tn5 insertion bias should be considered (only if a strand column is provided).
#' @param saveHdf5 If assays should be saved as HDF5 files.
#' @param outDir Directory to save HDF5 file to.
#' @param BPPARAM Parallel back-end to be used. Passed to [BiocParallel::bplapply()].
#' @return [MultiAssayExperiment::MultiAssayExperiment-class] with Motif, ATAC- & ChIP-seq experiments.
#' @import MultiAssayExperiment
#' @importFrom rhdf5 H5Fcreate h5createDataset h5writeDataset h5delete h5write h5ls H5Fopen H5Fclose H5garbage_collect
#' @importFrom BiocParallel bplapply SerialParam MulticoreParam SnowParam
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges colData cbind assays
#' @importClassesFrom SummarizedExperiment SummarizedExperiment RangedSummarizedExperiment
#' @importClassesFrom HDF5Array HDF5Array
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @export
prepData <- function(refCoords,
                     motifData,
                     atacData,
                     chIPData,
                     testSet=NULL,
                     promoterCoords=NULL,
                     aggregationFun=max,
                     scoreCol="score",
                     weightCol=NULL,
                     isUncertainCol=NULL,
                     annoCol="context",
                     shift=FALSE,
                     saveHdf5=FALSE,
                     outDir=FALSE,
                     BPPARAM=SerialParam()) {

  #TODO: Add flag for training (matched and not in test list),
  # prediction (not matched and not in test), test ones in vector

  # convert reference coordinates to GRanges
  if(is.data.frame(refCoords) || is.data.table(refCoords))
  {
    refCoords <- makeGRangesFromDataFrame(as.data.frame(refCoords),
                                          ignore.strand=TRUE)
  }

  # Preparing ChIP-seq data -----------------------------------------------------
  message("Processing ChIP-seq data")
  chIPSe <- .mapSeqData(chIPData, refCoords, type="ChIP",
                        weightCol=weightCol,
                        aggregationFun=aggregationFun,
                        isUncertainCol=isUncertainCol,
                        saveHdf5=saveHdf5, outDir=outDir,
                        BPPARAM=BPPARAM)
  covTfs <- unique(colData(chIPSe)[[tfNameCol]])
  chIPMap <- data.frame(primary=colData(chIPSe)[[annoCol]],
                        colname=colData(chIPSe)[["combination"]],
                        stringsAsFactors=FALSE)

  # Preparing ATAC-seq data ----------------------------------------------------
  message("Processing ATAC-seq data")
  atacSe <- .mapSeqData(atacData, refCoords, type="ATAC", annoCol=annoCol,
                        shift=shift, saveHdf5=saveHdf5, outDir=outDir,
                        BPPARAM=BPPARAM)
  atacMap <- data.frame(primary=colData(atacSe)[[annoCol]],
                        colname=colData(atacSe)[[annoCol]],
                        stringsAsFactors=FALSE)

  allContexts <- union(colData(chIPSe)[[annoCol]], colData(atacSe)[[annoCol]])
  matchedContexts <- intersect(colData(chIPSe)[[annoCol]], colData(atacSe)[[annoCol]])

  # Preparing motif data --------------------------------------------------------
  message("Processing Motif matching scores")
  motifSe <- .mapSeqData(motifData, refCoords, type="Motif",
                         aggregationFun=aggregationFun,
                         saveHdf5=saveHdf5, outDir=outDir,
                         scoreCol=scoreCol, BPPARAM=BPPARAM)

  motifMap <- data.frame(primary=rep(allContexts, ncol(motifSe)),
                         colname=rep(colnames(motifSe),
                                     each=length(allContexts)),
                         stringsAsFactors=FALSE)

  # Preparing ATAC-seq promoter data --------------------------------------------
  if(!is.null(promoterCoords)){
  message("Processing ATAC-seq promoter data")
    atacPromSe <- .mapSeqData(atacData, promoterCoords, type="ATAC",
                              annoCol=annoCol,
                              saveHdf5=saveHdf5,
                              fileName="ATAC_promoters_mapped",
                              outDir=outDir, BPPARAM=BPPARAM)
    atacPromMap <- data.frame(primary=colData(atacPromSe)[[annoCol]],
                              colname=colData(atacPromSe)[[annoCol]],
                              stringsAsFactors=FALSE)
    listMap <- list(atacMap, chIPMap, motifMap, atacPromMap)
    expNames <- c(atacExp, chIPExp, motifExp, atacPromExp)
    names(listMap) <- expNames
    objList <- list(atacSe, chIPSe, motifSe, atacPromSe)
    names(objList) <- expNames
  }
  else{
    listMap <- list(atacMap, chIPMap, motifMap)
    expNames <- c(atacExp, chIPExp, motifExp)
    names(listMap) <- expNames
    objList <- list(atacSe, chIPSe, motifSe)
    names(objList) <- expNames
  }

  dfMap <- listToMap(listMap)

  # Create MultiAssayExperiment Object -----------------------------------------
  annoData <- data.frame(row.names=allContexts)
  annoData[[annoCol]] <- allContexts

  mae <- MultiAssayExperiment(experiments=objList,
                              colData=annoData,
                              sampleMap=dfMap)

  # Annotate MultiAssayExperiment object ---------------------------------------

  # actually this is not strictly required in the sampleMap
  annoDt <- as.data.table(sampleMap(mae))
  annoDt[,eval(isTestCol):=fifelse(primary %in% testSet, TRUE, FALSE)]
  annoDt[,eval(isTrainCol):=fifelse(!(primary %in% testSet) &
                                      primary %in% matchedContexts, TRUE, FALSE)]

  primary <- colData(mae)[[annoCol]]
  colData(mae)[[isTestCol]] <- fifelse(primary %in% testSet, TRUE, FALSE)
  colData(mae)[[isTrainCol]] <- fifelse(!(primary %in% testSet) &
                                      primary %in% matchedContexts, TRUE, FALSE)

  sampleMap(mae) <- annoDt

  return(mae)
}
