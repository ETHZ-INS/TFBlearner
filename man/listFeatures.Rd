% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/listFeatures.R
\name{listFeatures}
\alias{listFeatures}
\title{Overview and description of features.}
\usage{
listFeatures()
}
\value{
\link[data.table:data.table]{data.table::data.table} containing an overview of the features.
}
\description{
Provides an overview of the features.
Features are grouped into five types, which are explained in the following :
\enumerate{
\item site-Features: Specific for the sequence of a site.
\item context-Features: Features specific for the ATAC-profile of a cellular context.
\item pan-context-Features: Features of a site computed across ATAC-profiles of different cellular contexts.
\item tf-Features: Features specific for a TF but not the cellular context.
\item context-tf-Features: Features specific for a TF and the cellular context.
}
}
\details{
Columns in the feature matrix (as returned by \link{getFeatureMatrix}) are named the following:
\verb{<featureType>_<featureName>} depending of the feature followed by (\verb{<motifSuffix>, <i>, <normedSuffix>}).

The column \code{parameters} contains the names of parameters for the respective feature, that can be passed to function computing the feature (listed in \code{feature_function_exported}).
The column \code{context_normed} indicates if the feature does not need to be normalized across cellular contexts when constructing the feature matrix.
Note feature matrix column names and assay names of the MultiAssayExperiment object are the same.
}
\author{
Emanuel Sonder
}
