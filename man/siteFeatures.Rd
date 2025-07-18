% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/siteFeatures.R
\name{siteFeatures}
\alias{siteFeatures}
\title{Site-specific features}
\usage{
siteFeatures(
  mae,
  aggregationFun = max,
  annoData = NULL,
  scoreCols = NULL,
  features = c("Sequence", "Width", "Annot"),
  phast = phastCons100way.UCSC.hg38,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  ...
)
}
\arguments{
\item{mae}{\link[MultiAssayExperiment:MultiAssayExperiment-class]{MultiAssayExperiment::MultiAssayExperiment} as construced by \code{\link[=prepData]{prepData()}} containing Motif, ATAC- and ChIP-seq data.}

\item{aggregationFun}{Aggregation function for aggregating genomic scores overlapping the rowRanges of the provided MultiAssayExperiment object.}

\item{annoData}{Further data to be aggregated across the rowRanges of the provided MultiAssayExperiment object.
Named list of paths pointing to .bed / .bam files or data.frames / \link[data.table:data.table]{data.table::data.table} or \link[GenomicRanges:GRanges-class]{GenomicRanges::GRanges} containing chr, start and end columns and score columns.}

\item{scoreCols}{Column names of the scores of annoData to be aggregated. Needs to be same length as annoData.
If of length one, it is assumed that score columns of all list elements have the same name.}

\item{features}{Names of features to be added. Can be all or some of "Sequence", "Width", "Annot".
Features are stored in the assays of the added experiment. See \link{listFeatures} for an overview of the features.}

\item{phast}{phastCons conservation scores to be used.}

\item{genome}{\link[BSgenome:BSgenome-class]{BSgenome::BSgenome} to be used.}

\item{...}{Arguments to be passed to \link{genomicRangesMapping} for aggregation of annoData.}
}
\value{
\link[MultiAssayExperiment:MultiAssayExperiment-class]{MultiAssayExperiment::MultiAssayExperiment} object with an experiment containing site-specific features added to \link[MultiAssayExperiment:MultiAssayExperiment-methods]{MultiAssayExperiment::experiments}.
}
\description{
Adds an experiment with features specific for each site, such as CpG-content and conservation scores, to the provided MultiAssayExperiment object.
}
\author{
Emanuel Sonder
}
