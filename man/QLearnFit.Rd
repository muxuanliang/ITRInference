% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QLearnFit.R
\name{QLearnFit}
\alias{QLearnFit}
\title{This function implements the Q-learning estimation for individualized treatment rule and the inference procedure based on the de-correlated score (see reference).}
\usage{
QLearnFit(data, intercept = FALSE, standardize = TRUE)
}
\arguments{
\item{data}{A list - list(predictor = x, treatment = trt, outcome = y)), where x is the covariate matrix, trt is 0 or 1 (1 indicates treatment), y is the outcome.}

\item{intercept}{includes intercept or not}

\item{standardize}{whether standardize the input covariant matrix.}
}
\value{
A list
 \describe{
\item{fit}{Use \code{glmnet} to fit a nigh-dimensional Q-learning.}
\item{pseudoPredictor}{Augmented design matrix with covariate and treatment interations for Q-learning.}
\item{pseudoTreatment}{Treatment coded in 1 and -1.}
\item{pseudoOutcome}{Outcome used in the inferenc step.}
}
}
\description{
This function implements the Q-learning estimation for individualized treatment rule and the inference procedure based on the de-correlated score (see reference).
}
\references{
Muxuan Liang, Young-Geun Choi, Yang Ning, Maureen Smith, Yingqi Zhao (2020). Estimation and inference on high-dimensional individualized treatment rule in observational data using split-and-pooled de-correlated score.
}
\author{
Muxuan Liang <mliang@fredhutch.org>
}
