#' rebaseMaeH5paths - Change the paths of the HDF5 files in a MultiAssayExperiment object.
#'
#' @param mae A \code{\link[MultiAssayExperiment]{MultiAssayExperiment-class}}
#'   object.
#' @param newBase The new directory in which the HDF5 files are located. This
#'   can either be a single path, or a named vector of paths corresponding to
#'   the different experiments contained in `mae`.
#'
#' @returns The `mae` object with updated HDF5 file links.
#' @export
rebaseMaeH5paths <- function(mae, newBase=vector("character")){
  stopifnot(is.character(newBase) & length(newBase)>0)
  delayedExp <- names(experiments(mae))[sapply(experiments(mae), \(x){
    any(unlist(lapply(assays(x), inherits, what="DelayedArray")))
  })]
  if(length(delayedExp)==0) return(mae)
  if(length(newBase)==1)
    newBase <- sapply(setNames(delayedExp,delayedExp), \(x) newBase)
  if(length(missingExp <- setdiff(delayedExp, names(newBase)))>0)
    stop("The following experiments do not have a new base path specified: ",
         paste(missingExp, collapse=", "))
  if(length(notused <- setdiff(names(newBase), delayedExp))>0)
    message("The following new base paths were ignored for lack of matching ",
            "experiments:", paste(notused, collapse=", "))
  for(e in delayedExp){
    for(a in assayNames(experiments(mae)[[e]])){
      if(inherits(assay(experiments(mae)[[e]], a), "DelayedArray")){
        b <- basename(assay(experiments(mae)[[e]], a)@seed@seed@filepath)
        newPath <- file.path(newBase[[e]], b)
        assay(experiments(mae)[[e]], a)@seed@seed@filepath <- newPath
      }
    }
  }
  newPaths <- unlist(getMaeH5paths(mae))
  newPaths <- newPaths[!sapply(newPaths, is.null)]
  newPaths <- newPaths[!file.exists(newPaths)]
  if(length(newPaths)>0)
    warning("The object has been modified, however some of the specified files",
            " do not exist:\n", paste(newPaths, collapse="\n"))
  mae
}

#' getMaeH5paths - extract HDF5 file paths from MAE
#'
#' @param x A \code{\link[MultiAssayExperiment]{MultiAssayExperiment-class}}
#'   object, or an object inheriting
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}.
#'
#' @returns A named list containing the path(s) to the HDF5 file(s) underlying
#'   each experiment.
#' @export
getMaeH5paths <- function(x){
  if(inherits(x, "MultiAssayExperiment"))
    return(lapply(experiments(x), getMaeH5paths))
  if(!inherits(x, "SummarizedExperiment"))
    stop("`x` should inherit either SummarizeExperiment or MultiAssayExperiment.")
  p <- lapply(assays(x), \(y){
    if(!inherits(y, "DelayedArray")) return(NULL)
    y@seed@seed@filepath
  })
  unique(unlist(p[!sapply(p, is.null)]))
}

#' loadMae - load a MultiAssayExperiment object
#'
#' @param filepath The path to the RDS file containing a
#'   \code{\link[MultiAssayExperiment]{MultiAssayExperiment-class}} object.
#' @param keepAbsPaths Logical; whether to keep the absolute paths to HDF5 files
#'   indicated in the object, or to try to update them to the object's location.
#'   By default, the original paths will be preserved if all the HDF5 files are
#'   accessible, and otherwise looked for in the directory of `filepath`.
#'
#' @returns A \code{\link[MultiAssayExperiment]{MultiAssayExperiment-class}}
#'   object.
#' @export
loadMae <- function(filepath, keepAbsPaths=NULL){
  mae <- readRDS(filepath)
  if(isTRUE(keepAbsPaths)) return(mae)
  if(is.null(keepAbsPaths)){ # not specified
    # check if old absolute paths are still accessible, if so assume keep
    h5paths <- unlist(getMaeH5paths(mae))
    h5paths <- h5paths[!sapply(h5paths, is.null)]
    if(all(unlist(lapply(h5paths, file.exists)))) return(mae)
    # check if the files exist in the same folder:
    newPaths <- lapply(basename(h5paths), \(x){
      file.exists(file.path(dirname(filepath), x))
    })
    if(!all(unlist(newPaths)))
      stop("The object you are trying to load has delayed (i.e. on-disk) ",
           "assays pointing to HDF5 files that are not accessible. Please ",
           "ensure that these are placed in the same folder as the rds file.")
  }
  rebaseMaeH5paths(mae, dirname(filepath))
}

#' saveMae - Save a MultiAssayExperiment
#'
#' @param mae A \code{\link[MultiAssayExperiment]{MultiAssayExperiment-class}}
#'   object.
#' @param filepath The path of the RDS file in which to save the object.
#' @param copyH5 Logical; whether to make copies of any underlying HDF5 files
#'   in the directory of `filepath` and link the saved `mae` to those.
#'
#' @returns NULL.
#' @export
saveMae <- function(mae, filepath, copyH5=FALSE){
  if(!copyH5) return(base::saveRDS(mae, filepath))
  stopifnot(dir.exists(dirname(filepath)))
  h5paths <- h5paths[!sapply(h5paths, is.null)]
  for(f in unlist(getMaeH5paths(mae))){
    if(!is.null(f)){
      newpath <- file.path(dirname(filepath), basename(f))
      if(f==newpath) stop("Source and destination h5 are the same.")
      file.copy(f, newpath)
    }
  }
  mae <- rebaseMaeH5paths(mae, dirname(filepath))
  base::saveRDS(mae, filepath)
}

