% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeHerb_ChemNet.R
\name{makeHerb_ChemNet}
\alias{makeHerb_ChemNet}
\title{makeHerb_ChemNet}
\usage{
makeHerb_ChemNet(NetworkData, NetworkType = "HerbChem", INppi = T)
}
\arguments{
\item{NetworkData}{a dataframe with colnames:"herb","cid","chemical_name","geneID","score","inchikey", "smiles","database","QED_DES"}

\item{NetworkType}{'HerbChem','ChemGene',or 'HerbChemGene'}

\item{INppi}{whether to seed in PPI network}

\item{randomtimes}{randomtimes for permutaion test}
}
\value{
a list with igraph Object
}
\description{
Network for Visualization
}
\author{
Yang Ming
}
