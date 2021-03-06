% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_StudySeqFunctions.R
\name{condition_haplos_no_cSNV}
\alias{condition_haplos_no_cSNV}
\title{Reduce haplos to contain non-cSNV data}
\usage{
condition_haplos_no_cSNV(haplos, RV_pool_loc)
}
\arguments{
\item{haplos}{sparseMatrix.  The sparseMatrix of genomes returned by \code{read_slimOut}.}

\item{RV_pool_loc}{The column locations of each SNV (in haplos)  for SNVs in the pool of cRVs.}
}
\value{
The reduced haplotype matrix.
}
\description{
Reduce haplos to contain non-cSNV data
}
\keyword{internal}
