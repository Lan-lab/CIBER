% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findFeatures_func.R
\name{dsRank}
\alias{dsRank}
\title{Rank genes according to diffScores}
\usage{
dsRank(diffScores)
}
\arguments{
\item{diffScores}{list of previously calculated diffScores}
}
\value{
list of ranked gene names
}
\description{
Rank genes according to diffScores
}
\examples{
gl.diffScore <- dsRank(ds.diffScore)
}
