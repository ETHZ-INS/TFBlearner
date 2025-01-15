# Costum MultiassayExperiment object constructor
# if is path, data will be mapped.

# write a validator! to check if objects are of the correct shape !!
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
    if(sum(grepl("contextFeat", experimentNames, ignore.case=TRUE))==0){
      stop("Context-features as obtained by using contextTffeatures() on the object
            are required for further use.")
    }
  }

  #TODO: Add further checks, current form merely captures the idea.
  # - that origin path is available in colData of ATAC-seq !!
  #message("The provided MultiAssayExperiment object is of correct form")

  return(TRUE)
}

.addExperiment <- function(mae,
                           fm,
                           coords,
                           mapTo=c("All", "Tf", "Col"), # all essentially the same as coord
                           prefix="features",
                           annoCol="context",
                           tfName=NULL,
                           colsToMap=NULL, ...){

  # get coldata of feature matrix
  names(fm) <- paste(prefix, names(fm), sep="_")
  contexts <- unique(sampleMap(mae)$primary) #TODO: Check how protect is the primary naming

  if(mapTo=="All"){
    colFeatData <- DataFrame(feature_type=prefix,
                             row.names=prefix)

    rseFeat <- SummarizedExperiment(assays=fm,
                                    rowRanges=coords,
                                    colData=colFeatData,
                                    checkDimnames=FALSE)

    mapFeat <- data.frame(primary=rep(contexts, nrow(colFeatData)),
                          colname=rep(rownames(colFeatData),
                                      each=length(contexts)),
                          stringsAsFactors=FALSE)
    contextMap <- data.frame(row.names=contexts,
                             description=rep(annoCol, length(contexts)))
  }
  else{
    if(is.null(tfName)) stop("Provide the name of the TF for which features are added")

    if(mapTo=="Tf"){
       colDataChIP <- colData(experiments(mae)$ChIP)
       colsToMap <- unique(subset(colDataChIP, tf_name==tfName)[[annoCol]])
       colFeatData <- DataFrame(feature_type=prefix,
                                tf_name=tfName,
                                row.names=paste(prefix, tfName, sep="_"))
    }

    if(mapTo=="Col"){
      if(is.null(colsToMap)) stop("Column to match to needs to be provided")
        colFeatData <- DataFrame(feature_type=prefix,
                                 tf_name=tfName,
                                 row.names=paste(prefix, tfName,
                                                 colsToMap, sep="_"))
        colFeatData[[annoCol]] <- colsToMap
    }

    rseFeat <- SummarizedExperiment(assays=fm,
                                    rowRanges=coords,
                                    colData=colFeatData,
                                    checkDimnames=FALSE)

    mapFeat <- data.frame(primary=rep(colsToMap, nrow(colFeatData)),
                          colname=rep(rownames(colFeatData),
                                      each=length(colsToMap)),
                          stringsAsFactors=FALSE)
    contextMap <- data.frame(row.names=colsToMap,
                             description=rep(annoCol, length(colsToMap)))
  }
  # create multiassayexperiment object for feature to be added
  maeFeat <- .createMaeFeat(mapFeat, rseFeat, contextMap, mae, prefix=prefix)

  # Merge with existing object
  mae <- c(mae, maeFeat)
  return(mae)
}

.createMaeFeat <- function(mapFeat, rseFeat, contextMap, mae, prefix){

  listMap <- list(mapFeat)
  names(listMap) <- prefix
  dfMap <- listToMap(listMap)
  objList <- list(rseFeat)
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
  sampleMap(maeFeat) <- mapFeat

  return(maeFeat)
}

#' @import MultiAssayExperiment
.addColumn <- function(mae,
                       fm,
                       coords,
                       mapTo=c("All", "Tf", "Col"), # all essentially the same as coord
                       prefix="features",
                       annoCol="context",
                       tfName=NULL,
                       colsToMap=colsToMap, ...){
  warning("Assays which are not shared are dropped")
  seOrig <- experiments(mae)[[prefix]]
  names(fm) <- paste(prefix, names(fm), sep="_")

  if(mapTo=="Tf" | mapTo=="Col"){
    if(is.null(tfName)) stop("Provide the name of the TF for which features are added")
    if(mapTo=="Tf"){
      colDataChIP <- colData(experiments(mae)$ChIP)
      colsToMap <- unique(subset(colDataChIP, tf_name==tfName)[[annoCol]])

      colFeatData <- DataFrame(feature_type=prefix,
                               tf_name=tfName,
                               row.names=paste(prefix, tfName, sep="_"))
    }

    if(mapTo=="Col"){
      if(is.null(colsToMap)) stop("Column to match to needs to be provided")
      colFeatData <- DataFrame(feature_type=prefix,
                               tf_name=tfName,
                               row.names=paste(prefix, tfName,
                                               colsToMap, sep="_"))
      colFeatData[[annoCol]] <- colsToMap
    }

    seFeat <- SummarizedExperiment(assays=fm,
                                   rowRanges=coords,
                                   colData=colFeatData,
                                   checkDimnames=FALSE)

    mapFeat <- data.table(primary=rep(colsToMap, nrow(colFeatData)),
                          colname=rep(rownames(colFeatData),
                                      each=length(colsToMap)),
                          stringsAsFactors=FALSE)
    contextMap <- data.frame(row.names=colsToMap,
                             description=rep(annoCol, length(colsToMap)))
  }
  else{
    stop("Other cases like adding an assay are not implemented (yet)")
  }

  commonAssays <- intersect(names(assays(seOrig)), names(assays(seFeat)))
  assays(seOrig) <- assays(seOrig)[commonAssays]
  assays(seFeat) <- assays(seFeat)[commonAssays]
  seComb <- SummarizedExperiment::cbind(seOrig, seFeat)

  mapOrig <- as.data.table(subset(sampleMap(mae), assay==prefix)[,c("primary",
                                                                    "colname")])
  mapComb <- rbind(mapOrig, mapFeat, use.names=TRUE, fill=TRUE)
  mapComb <- as.data.frame(mapComb)
  contextMap <- data.frame(row.names=unique(mapComb$primary))
  contextMap$description <- annoCol

  # add back combined assay to original mae
  # create multiassayexperiment object for feature to be added
  maeFeat <- .createMaeFeat(mapComb, seComb, contextMap, mae, prefix=prefix)

  #TODO:  is there a better way of doing this?
  # remove experiments from mae
  experiments(mae)[[prefix]] <- NULL

  # Merge with existing object
  mae <- c(mae, maeFeat)

  return(mae)
}

# Function to help adding Features !!
.addFeatures <- function(mae,
                         fm,
                         mapTo=c("All", "Tf", "Col"), # all essentially the same as coord
                         prefix="features",
                         annoCol="context",
                         tfName=NULL,
                         colsToMap=NULL,
                         coords=NULL, ...){
  # object validator
  .checkObject(mae)

  mapTo <- match.arg(mapTo, choices=c("All", "Tf", "Col"))
  if(mapTo=="Tf" & is.null(tfName)) stop("Provide name of TF to map features to")

  if(is.null(coords)) coords <- rowRanges(experiments(mae)$Motifs)

  #TODO: If prefix exists as an experiment only add Columns!!!!
  if(!(prefix %in% names(experiments(mae)))){
    mae <- .addExperiment(mae, fm, coords, mapTo=mapTo, # all essentially the same as coord
                          prefix=prefix, annoCol=annoCol, tfName=tfName,
                          colsToMap=colsToMap, ...)
  }
  else{
    mae <- .addColumn(mae, fm, coords, mapTo=mapTo, # all essentially the same as coord
                      prefix=prefix, annoCol=annoCol, tfName=tfName,
                      colsToMap=colsToMap, ...)
  }


  #else{
  #  stop("Other cases are not implemented yet, e.g. adding")
  #}

  return(mae)
}

# helper function to get the cellular contexts per TF
getContexts <- function(mae, tfName){
  .checkObject(mae)

  colDatChIP <- colData(experiments(mae)$ChIP)
  contexts <- subset(colDatChIP, tf_name==tfName)$context
  return(contexts)
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
#' @param BPPARAM Parallel back-end to be used. Passed to [BiocParallel::bplapply()].
#' @return [MultiAssayExperiment::MultiAssayExperiment-class] with Motif, ATAC- & ChIP-seq experiments.
#' @import MultiAssayExperiment
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
                       shift=shift, BPPARAM=BPPARAM)
  atacMap <- data.frame(primary=colData(atacSe)[[annoCol]],
                        colname=colData(atacSe)[[annoCol]],
                        stringsAsFactors=FALSE)

  # Preparing ChIP-seq data -----------------------------------------------------
  message("Processing ChIP-seq data")
  chIPSe <- .mapSeqData(chIPData, refCoords, type="ChIP",
                       weightCol=weightCol,
                       aggregationFun=aggregationFun,
                       isUncertainCol=isUncertainCol,
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
  show(mae)

  return(mae)
}
