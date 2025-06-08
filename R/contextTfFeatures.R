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
#' @param features Names of features to be added. Can be all or some of "Inserts", "Weighted_Inserts", "ChromVAR_Scores".
#' "Insert" features will always be computed.
#' See [TFBlearner::listFeatures] for an overview of the features.
#' @param annoCol Name of column indicating cellular contexts in colData.
#' @param insertionProfile Pre-computed insertion footprint profile for the specified transcription factor.
#' Needs to contain coordinate (chr/seqnames, start, end) columns and weight column (termed "w").
#' @param aggregationFun function (e.g. mean, median, sum) used to aggregate features across the rowRanges of experiments of the
#' provided [MultiAssayExperiment::MultiAssayExperiment-class] object.
#' @param seed Integer value for setting the seed for random number generation with [base::set.seed].
#' @param BPPARAM Parallel back-end to be used. Passed to [BiocParallel::bpmapply()] & [BiocParallel::bplapply()].
#' @param ... Arguments passed to [TFBlearner::getInsertionProfiles].
#' @return [MultiAssayExperiment::MultiAssayExperiment-class] object with an experiment containing transcription factor- and cellular context-specific features added to [MultiAssayExperiment::experiments].
#' If already an experiment of this feature group exists, columns are added to it.
#' @import MultiAssayExperiment
#' @importFrom BiocParallel bpmapply bplapply SerialParam MulticoreParam SnowParam register
#' @export
contextTfFeatures <- function(mae,
                              tfName,
                              addLabels=TRUE,
                              whichCol=c("All", "OnlyTrain", "Col"), # evt. rename these arguments
                              colSel=NULL,
                              features=c("Inserts",
                                         "Weighted_Inserts",
                                         "ChromVAR_Scores"),
                              annoCol="context",
                              insertionProfile=NULL,
                              aggregationFun=sum,
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
  features <- match.arg(features, choices=c("Inserts", "Weighted_Inserts",
                                            "ChromVAR_Scores"),
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
  motifPath <- subset(colData(maeSub[[motifExp]]),
                      get(motifNameCol)==tfName)$origin
  motifRanges <- readRDS(motifPath)

  if(addLabels){
    colDataChIP <- colData(mae[[chIPExp]])
    colDataChIP <- subset(colDataChIP, get(tfNameCol)==tfName)
    labelCols <- colDataChIP$combination
    names(labelCols) <- colDataChIP[[annoCol]]
    labels <- lapply(labelCols, function(col){
      as(assays(mae[[chIPExp]])[[peakAssay]][,col,drop=TRUE],
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
                                    insRes[[retScoresName]],
                                    scoreCol=scoreCol,
                                    byCols="type",
                                    aggregationFun=aggregationFun,
                                    BPPARAM=BPPARAM)
      colnames(feats) <- paste(paste(scoreCol, colnames(feats), sep="."),
                               tfMotifPrefix, 1, sep="_")
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
      labels <- list(Matrix::Matrix(labels))
      names(labels) <- labelName
      insFeats <- append(insFeats, labels)
    }

    return(insFeats)
  }, contexts, labels,
  MoreArgs=list(coords=coords, atacFrag=atacFragPaths,
                motifRanges=motifRanges, features=features,
                profile=insertionProfile, tfName=tfName,
                threads=threads, aggregationFun=aggregationFun,
                BPPARAM=BPPARAM, ...),
  SIMPLIFY=FALSE)
  names(feats) <- contexts

  if("ChromVAR_Scores" %in% features){
    message("Get chromVAR features")

    if(!(actExp %in% names(experiments(mae))) |
       !(preSelActCol %in% colnames(colData(mae[[tfFeat]])))){
      warning("ChromVAR activity estimates can not be added if tfFeatures() with `Associated_Motif_Activity` and panContextFeatures() have not been called before")}
    else{
      selActMotifs <- unlist(subset(colData(mae[[tfFeat]]),
                                    get(tfNameCol)==tfName)[[preSelActCol]])
      devMat <- t(assays(mae[[actExp]])[[normDevAssay]][selActMotifs, contexts])
      devMat <- as.matrix(devMat)

      devMats <- lapply(1:length(contexts), function(i){
        idx <- rep(i, each=length(coords))
        dm <- Matrix::Matrix(devMat[idx,,drop=FALSE])
        return(dm)})
      names(devMats) <- contexts

      # add to features
      actFeatNames <-  paste(chromVarFeatName, names(selActMotifs), sep="_")
      feats <- lapply(contexts, function(context){
        devMat <- devMats[[context]]
        dev2Mat <- lapply(colnames(devMat), function(col){
          devMat[,col,drop=FALSE]})
        names(dev2Mat) <- actFeatNames
        c(feats[[context]], dev2Mat)
      })
      names(feats) <- contexts
    }
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
    mae <- .addFeatures(mae, seTfFeat, colsToMap=context, prefix=contextTfFeat)
  }

  return(mae)
}
