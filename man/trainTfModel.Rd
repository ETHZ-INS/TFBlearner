% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trainTfModel.R
\name{trainTfModel}
\alias{trainTfModel}
\title{Training transcription factor-specific tree-based gradient boosting Models}
\usage{
trainTfModel(
  tfName,
  fm,
  measureName = c("classif.aucpr", "classif.logloss"),
  evalRounds = 100,
  earlyStoppingRounds = 10,
  posFrac = 0.25,
  loContext = FALSE,
  tuneHyperparams = TRUE,
  stackingStrat = c("last", "wLast", "wMean", "boostTree"),
  subSample = 1e+05,
  valChrs = NULL,
  annoCol = "context",
  seed = 42,
  numThreads = 10,
  BPPARAM = SerialParam()
)
}
\arguments{
\item{tfName}{Name of transcription factor to train model for.}

\item{fm}{\link[SummarizedExperiment:RangedSummarizedExperiment-class]{SummarizedExperiment::RangedSummarizedExperiment} object containing features & labels as obtained by \link{getFeatureMatrix}.}

\item{measureName}{Measure used for hyperparameter selection.
Either area under the precision-recall curve computed using \link[PRROC:pr.curve]{PRROC::pr.curve} ("classif.aucpr") or
logloss as implemented by \link[mlr3measures:logloss]{mlr3measures::logloss}.}

\item{evalRounds}{Number of evaluation rounds for the hyperparameter selection rounds.}

\item{earlyStoppingRounds}{Number of early stopping rounds for the hyperparameter selection and training of the \link[lightgbm:lightgbm]{lightgbm::lightgbm} model.}

\item{posFrac}{Fraction of positives to use for the training of the model.
Negatives will be supsampled to achieve the specified fraction.}

\item{loContext}{Should cellular-contexts be used for leave-one-context-out rounds during hyperparameter selection.
Only works if more than one cellular-context is contained within the feature matrix.}

\item{tuneHyperparams}{If hyperparameters should be tuned with \link[mlr3mbo:mlr_tuners_mbo]{mlr3mbo::TunerMbo}. Recommend to have this turned on (\code{tuneHyperparams=TRUE}).
Otherwise (hopefully) sensible defaults are used.}

\item{stackingStrat}{Stacking strategy to use. \code{last}, chooses the the model which has been trained on all (most) positives not using
observational weights for the ChIP-seq peaks. \code{wLast} using the last model which has seen most positives and has been trained with observational weights.
\code{wMean} weighted mean all models based on performance on the feature matrix provided.
\code{boostTree} Trains a lightgbm model on the predictions of the models in the bag, together with some additional features (e.g. gc_content, total_overlaps).}

\item{subSample}{Number of rows of featMat which should be used for computing performance estimates. Only used if \code{stackingStrat="wMean"}.}

\item{valChrs}{Optional, holdout chromosomes to compute performance metrics (precision (P), recall (R), area under the precision-recall curve (AUC-PR))
and estimate the threshold for dichotomization of the probabilites.}

\item{annoCol}{Name of column indicating cellular contexts.}

\item{seed}{Integer value for setting the seed for random number generation with \link[base:Random]{base::set.seed}.}

\item{numThreads}{Total number of threads to be used. In case \link[BiocParallel:MulticoreParam-class]{BiocParallel::MulticoreParam} or \link[BiocParallel:SnowParam-class]{BiocParallel::SnowParam} with several workers are
are specified as parallel back-ends, \code{floor(numThreads/nWorker)} threads are used per worker.}

\item{BPPARAM}{Parallel back-end to be used. Passed to \code{\link[BiocParallel:bpmapply]{BiocParallel::bpmapply()}}.}
}
\value{
A list of four \link[lightgbm:lightgbm]{lightgbm::lightgbm} models trained on different strata of the data and stacked model combining the predictions of the later.
}
\description{
Trains a bag of four tree-based gradient boosting models of the \link[lightgbm:lightgbm]{lightgbm::lightgbm} library and stacked model combining the predictions of the latter.
Hyperparameter selection is performed for each model seperately using model-based optimization by deploying the \href{https://mlr3tuning.mlr-org.com}{mlr3tuning} library.
The lightgbm classification learner used for the hyperparameter selection has been copied from the GitHub repository \url{https://github.com/mlr-org/mlr3extralearners}, whichs
contains the package mlr3extralearners developed by Raphael Sonabend and Patrick Schratz and Sebastian Fischer.
}
\author{
Emanuel Sonder
}
