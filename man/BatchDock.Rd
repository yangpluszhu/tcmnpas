% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BatchDock.R
\name{BatchDock}
\alias{BatchDock}
\title{BatchDock}
\usage{
BatchDock(
  chemList,
  chemListtype = "smiles",
  PDBList,
  PDBListtype = "file",
  dockType = "cross",
  python2PTAH = "~/mgltools_x86_64Linux2_1.5.6/bin/pythonsh",
  prepare_ligand4PATH = "~/ADtools",
  PSOvina = F
)
}
\arguments{
\item{chemList}{chemList\link{named character} OR file-->chemList:*.csv\link{colname:id,smiles}}

\item{chemListtype}{"smiles" OR "file"}

\item{PDBList}{PDBList:*.csv\link{colname:PDB,x,y,z,sizex,sizey,sizez} OR dataFrame-->PDBList\link{colname:PDB,x,y,z,sizex,sizey,sizez}}

\item{PDBListtype}{"file" OR "dataFrame"}

\item{dockType}{"cross" OR "parallel"}

\item{python2PTAH}{pythonPath for mgltools}

\item{prepare_ligand4PATH}{prepare_ligandPath for prepare functions}
}
\value{
a list Object
}
\description{
Batch Molecular Docking.
Note: mgltools must be installed and configured correctly.
}
\author{
Yang Ming
}
