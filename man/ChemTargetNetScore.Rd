% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ChemTargetNetScore.R
\name{ChemTargetNetScore}
\alias{ChemTargetNetScore}
\title{ChemTargetNetScore}
\usage{
ChemTargetNetScore(ChemData, diseaseGeneID, Methodnet = "KATZ", IFcombine = F)
}
\arguments{
\item{ChemData}{a data.frame with colnames=c('cid','geneID')}

\item{diseaseGeneID}{a character vector}

\item{Methodnet}{'KATZ', 'RW' or 'DIS'}

\item{IFcombine}{whether to combine the score}
}
\value{
a list Object
}
\description{
Calculate Compound Target Score based on Network based relevance score
}
\author{
Yang Ming
}
