% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepareLigand.R
\name{prepareLigand}
\alias{prepareLigand}
\title{prepareLigand}
\usage{
prepareLigand(
  molfileTXT,
  moltype = "smiles",
  preparedMOLfile = "outLigand.pdbqt",
  python2PTAH = "c:/Progra~2/MGLTools-1.5.6rc3/python2",
  prepare_ligand4PATH = system.file("tools", package = "tcmnpas")
)
}
\arguments{
\item{molfileTXT}{Molecular txt character}

\item{moltype}{Molecular type: smiles inchikey mol2 or pdb}

\item{preparedMOLfile}{prepared file name for output: *.pdbqt}

\item{python2PTAH}{python2.exe PTAH}

\item{prepare_ligand4PATH}{"prepare_ligand4.py" PTAH.The default Vina compound ligand preparation program (nphs–>: true to merge charges and remove non-polar hydrogens; lps–> true to merge charges and remove lone pairs of atoms) is used.}
}
\value{
a prepared molecular file
}
\description{
PrepareLigand for Molecular Docking.
Note: vina, PSOvina,and mgltools must be installed and configured correctly.
}
\author{
Yang Ming
}
