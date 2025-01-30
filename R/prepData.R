.checkObject <- function(mae,
                         checkFor=c("None", "Site", "TF", "Context")){

  if("None" %in% checkFor) checkFor <- "None"
  checkFor <- match.arg(checkFor, choices=c("None", "Site", "TF", "Context"),
                        several.ok=TRUE)

  if(!(is(mae, "MultiAssayExperiment"))){
    stop("Object provided needs to be a MultiAssayExperiment Object.
          For obtaining a object in the correct form use prepData()")
  }

  expectedExperiments <- c("ATAC", "ChIP", "Motifs")
  experimentNames <- names(experiments(mae))
  if(length(intersect(experimentNames, expectedExperiments))!=3){
    stop("Object provided needs to be a MultiAssayExperiment Object containing
          experiments: ATAC, ChIP and Motifs.
          For obtaining a object in the correct form use prepData()")
  }

  # check that motifs match to all contexts
  contexts <- unique(sampleMap(mae)$primary)
  sampleMapMotif <- as.data.table(subset(sampleMap(mae), assay=="Motifs"))
  contextMotifMap <- sampleMapMotif[,.(n_context=length(unique(primary))),
                                     by=colname]
  if(nrow(subset(contextMotifMap, n_context!=length(contexts)))>0){
    stop("Motifs need to map to all primaries (cellular context) labels.
          For obtaining a object in the correct form use prepData()")
  }

  if("Site" %in% checkFor){
    if(sum(grepl("siteFeat", experimentNames, ignore.case=TRUE))==0){
      stop("Site-specific features as obtained by using siteFeatures() on the object
            are required for further use.")
    }
  }

  if("TF" %in% checkFor){
    if(sum(grepl("tfFeat", experimentNames, ignore.case=TRUE))==0){
      stop("TF-features as obtained by using tfFeatures() on the object
            are required for further use.")
    }
  }

  if("Context" %in% checkFor){
    if(sum(grepl("contextTfFeat", experimentNames, ignore.case=TRUE))==0){
      stop("Context-features as obtained by using contextTffeatures() on the object
            are required for further use.")
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
    commonAssays <- intersect(names(assays(seOrig)), names(assays(seFeat)))
    assays(seOrig) <- assays(seOrig)[commonAssays]
    assays(seFeat) <- assays(seFeat)[commonAssays]

    # fill missing colData columns
    for(col in setdiff(colnames(colData(seOrig)), colnames(colData(seFeat)))){
      colData(seFeat)[[col]] <- rep(NA, nrow(colData(seFeat)))}

    # Avoid duplicate columns - keep newly added columns
    seFeat <- combineCols(
      seOrig[,!(colnames(seOrig) %in% colnames(seFeat))], seFeat,
      use.names=FALSE)

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
  mapMae <- subset(as.data.table(sampleMap(mae)), assay=="ATAC")
  mapFeat <- merge(mapFeat, mapMae[,c("primary", "is_testing", "is_training")],
                   by.x=c("primary"),
                   by.y=c("primary"), all.x=TRUE, all.y=FALSE,
                   allow.cartesian=TRUE)
  mapFeat[,is_testing:=fifelse(is.na(is_testing), FALSE, is_testing)]
  mapFeat[,is_training:=fifelse(is.na(is_training), FALSE, is_training)]

  sampleMap(maeFeat) <- mapFeat

  return(maeFeat)
}

.retainMappings <- function(mae){
  experimentNames <- intersect(names(experiments(mae)),
                               c("Motifs", "siteFeat", "tfFeat"))
  for(experimentName in experimentNames){
    mapComb <- as.data.table(sampleMap(mae))
    mapMotif <- subset(mapComb, assay==experimentName)
    allComb <- data.table(expand.grid(primary=unique(mapComb$primary),
                                      colname=unique(mapMotif$colname)))
    allComb <- allComb[,lapply(.SD, as.character), .SDcols=colnames(allComb)]
    allComb$assay <- experimentName
    allComb$is_training <- FALSE
    allComb$is_testing <- FALSE
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

  # ensure that motifs still map to every cellular context
  mae <- c(mae, maeFeat)
  mae <- .retainMappings(mae)

  return(mae)
}

# helper function to get the cellular contexts per TF
getContexts <- function(mae, tfName=NULL){
  .checkObject(mae)

  colsChIP <- colnames(experiments(mae)$ChIP)
  colsChIP <- as.data.table(tstrsplit(colsChIP, split="_"))
  if(!is.null(tfName)){
    contexts <- subset(colsChIP, V2==tfName)$V1}
  else{
    contexts <- unique(sampleMap(mae)$primary)}

  return(contexts)
}


#' Add further ATAC-seq datasets to object
#'
#' Adds additional columns for provided ATAC-seq datasets in the ATAC-seq experiments and
#' computes features if required.
#'
#' @name addATACData
#' @param mae [MultiAssayExperiment::MultiAssayExperiment-class] as construced by [TFBlearner::prepData()] containing Motif, ATAC-, ChIP-seq. If features should be computed
#' the object also needs to contain site-specific features as obtained by [TFBlearner::siteFeatures()] and transcription factor-specific features as obtained by [TFBlearner::tfFeatures()].
#' @param atacData Named list of [data.table::data.table]/data.frame/[GenomicRanges::GRanges]
#' or paths to .bam/.bed files containing ATAC-seq fragments, names being the cellular contexts labels.
#' Need to contain genomic coordinates (e.g. a chr/seqnames, start and end column).
#' @param testSet Vector of cellular context labels used for testing.
#' @param computeFeatures If features should be computed ([TFBlearner::contextTfFeatures()]).
#' @param tfName Name of transcription factor to compute context-specific features for.
#' @param annoCol Name of column indicating cellular contexts in colData.
#' @param shift Only for ATAC-seq data, if Tn5 insertion bias should be considered (only if a strand column is provided).
#' @param BPPARAM Parallel back-end to be used. Passed to [BiocParallel::bplapply()].
#' @param ... Additional arguments passed to [TFBlearner::contextTfFeatures()].
#' @return [MultiAssayExperiment::MultiAssayExperiment-class] with added ATAC-seq experiments (and cellular context and TF-specific features if `computeFeatures=TRUE`).
#' @import MultiAssayExperiment
#' @importFrom BiocParallel bplapply SerialParam MulticoreParam SnowParam
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges colData assays
#' @importFrom S4Vectors combineCols
#' @importFrom GenomeInfoDb seqlevelsStyle
addATACData <- function(mae, atacData,
                        testSet=NULL,
                        computeFeatures=FALSE,
                        tfName=NULL,
                        annoCol="context",
                        shift=FALSE,
                        BPPARAM=SerialParam(), ...){
  .checkObject(mae)

  coords <- rowRanges(experiments(mae)$Motifs)
  seAtac <- .mapSeqData(atacData, coords, type="ATAC", annoCol=annoCol,
                        shift=shift, BPPARAM=BPPARAM)

  mae <- .addFeatures(mae, seAtac, names(atacData), prefix="ATAC",
                      annoCol=annoCol)
  matchedContexts <- intersect(colData(experiments(mae)$ATAC)[[annoCol]],
                               colData(experiments(mae)$ChIP)[[annoCol]])

  annoDt <- as.data.table(sampleMap(mae))
  annoDt[,is_testing:=fifelse(primary %in% testSet, TRUE, is_testing)]
  annoDt[,is_testing:=fifelse(is.na(is_testing), FALSE, is_testing)]
  annoDt[,is_training:=fifelse(primary %in% matchedContexts & !is_testing,
                               TRUE, FALSE)]
  sampleMap(mae) <- annoDt

  colData(mae)$is_testing <- fifelse(colData(mae)[[annoCol]] %in% testSet,
                                     TRUE, colData(mae)$is_testing)
  colData(mae)$is_testing <- fifelse(is.na(colData(mae)$is_testing), FALSE,
                                     colData(mae)$is_testing)
  colData(mae)$is_training <- fifelse(colData(mae)[[annoCol]] %in% matchedContexts &
                                      !colData(mae)$is_testing,TRUE, FALSE)

  if(computeFeatures){
    if(is.null(tfName)) stop("Please provide a TF to compute the features for")
    mae <- contextTfFeatures(mae, tfName=tfName, whichCol="Col",
                             colSel=unique(names(atacData)), ...)
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
#' @param motifData  Named list of [data.table::data.table]/data.frame/[GenomicRanges::GRanges]
#' or paths to .bam/.bed files containing motif-matching scores, names being the motif names.
#' Need to contain genomic coordinates (e.g. a chr/seqnames, start and end column).
#' @param atacData Named list of [data.table::data.table]/data.frame/[GenomicRanges::GRanges]
#' or paths to .bam/.bed files containing ATAC-seq fragments, names being the cellular contexts labels.
#' Need to contain genomic coordinates (e.g. a chr/seqnames, start and end column).
#' @param chIPData Named list of [data.table::data.table]/data.frame/[GenomicRanges::GRanges]
#' or paths to .bam/.bed files containing ChIP-seq peaks, names being the combination of cellular context and TF labels (e.g. K562_YY1).
#' Need to contain genomic coordinates (e.g. a chr/seqnames, start and end column) and optionally a column with observational weights (weightCol)
#' and/or a uncertainty label (isUncertainCol).
#' @param testSet Vector of cellular context labels used for testing.
#' @param promoterCoords Optionally coordinates of promoters of transcription factor genes.
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
#' @import HDF5Array
#' @import MultiAssayExperiment
#' @importFrom rhdf5 h5createFile h5createDataset h5delete h5write H5close
#' @importFrom BiocParallel bplapply SerialParam MulticoreParam SnowParam
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges colData cbind assays
#' @importClassesFrom SummarizedExperiment SummarizedExperiment RangedSummarizedExperiment
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

  # Preparing ATAC-seq data ----------------------------------------------------
  message("Processing ATAC-seq data")
  atacSe <- .mapSeqData(atacData, refCoords, type="ATAC", annoCol=annoCol,
                        shift=shift, saveHdf5=saveHdf5, outDir=outDir,
                        BPPARAM=BPPARAM)
  atacMap <- data.frame(primary=colData(atacSe)[[annoCol]],
                        colname=colData(atacSe)[[annoCol]],
                        stringsAsFactors=FALSE)

  # Preparing ChIP-seq data -----------------------------------------------------
  message("Processing ChIP-seq data")
  chIPSe <- .mapSeqData(chIPData, refCoords, type="ChIP",
                       weightCol=weightCol,
                       aggregationFun=aggregationFun,
                       isUncertainCol=isUncertainCol,
                       saveHdf5=saveHdf5, outDir=outDir,
                       BPPARAM=BPPARAM)
  covTfs <- unique(colData(chIPSe)$tf_name)
  chIPMap <- data.frame(primary=colData(chIPSe)[[annoCol]],
                        colname=colData(chIPSe)[["combination"]],
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
                             annoCol=annoCol, BPPARAM=BPPARAM)
    atacPromMap <- data.frame(primary=colData(atacPromSe)[[annoCol]],
                              colname=colData(atacPromSe)[[annoCol]],
                              stringsAsFactors=FALSE)
    listMap <- list(atacMap, chIPMap, motifMap, atacPromMap)
    names(listMap) <- c("ATAC", "ChIP", "Motifs", "ATAC_promoters")
    objList <- list("ATAC"=atacSe,
                    "ChIP"=chIPSe,
                    "Motifs"=motifSe,
                    "ATAC_promoters"=atacPromSe)
  }
  else{
    listMap <- list(atacMap, chIPMap, motifMap)
    names(listMap) <- c("ATAC", "ChIP", "Motifs")
    objList <- list("ATAC"=atacSe,
                    "ChIP"=chIPSe,
                    "Motifs"=motifSe)
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
  annoDt[,is_testing:=fifelse(primary %in% testSet, TRUE, FALSE)]
  annoDt[,is_training:=fifelse(!(primary %in% testSet) &
                                 primary %in% matchedContexts, TRUE, FALSE)]

  primary <- colData(mae)[[annoCol]]
  colData(mae)$is_testing <- fifelse(primary %in% testSet, TRUE, FALSE)
  colData(mae)$is_training <- fifelse(!(primary %in% testSet) &
                                      primary %in% matchedContexts, TRUE, FALSE)

  sampleMap(mae) <- annoDt

  return(mae)
}
