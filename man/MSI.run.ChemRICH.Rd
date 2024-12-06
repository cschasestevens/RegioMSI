% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4.2.chemrich.R
\name{MSI.run.ChemRICH}
\alias{MSI.run.ChemRICH}
\title{ChemRICH Analysis Function}
\usage{
MSI.run.ChemRICH(l.cr.input, pval.no, fc.no, md, cl.var)
}
\arguments{
\item{l.cr.input}{Formatted ChemRICH input data frame.}

\item{pval.no}{Numeric position of p-value column.}

\item{fc.no}{Numeric position of fold change column.}

\item{md}{Number of metadata columns present in ChemRICH data frame.}

\item{cl.var}{Cluster variable.}
}
\value{
A ChemRICH results data frame containing cluster size, increased/decreased proportion, and information for ChemRICH dot plots.
}
\description{
Conducts ChemRICH analysis for selected comparison.This enrichment analysis method is adapted
from Chemical Similarity Enrichment Analysis (ChemRICH), initially described
in Barupal et. al. 2017 (DOI: 10.1038/s41598-017-15231-w). Statistically enriched classes are determined by a weighted Kolmogorov-Smirnov
Test with false discovery rate correction.
}
\examples{
chemrich.plot <- MSI.run.ChemRICH(
chemrich.input,
7,
16,
3,
"Class"
)

}
