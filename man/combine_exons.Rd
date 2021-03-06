% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combine_exons.R
\name{combine_exons}
\alias{combine_exons}
\title{Combine overlapping exons}
\usage{
combine_exons(exon_data)
}
\arguments{
\item{exon_data}{data.frame. This data frame must include named variables: \code{chrom}, a chromosome identifer; \code{exonStart}, the first position of the exon in base pairs; and \code{exonEnd}, the last position of the exon in base pairs.}
}
\value{
A data frame of combined exon segments.  This data frame includes the variables: \code{chrom}, a chromosome identifier; \code{exonStart}, the first position of the combined exon segment in base pairs; and \code{exonEnd}, the last position of the combined exon segment in base pairs.
}
\description{
Combine overlapping exons into a single observation
}
\examples{
# create a data frame that contains the
# the variables: chrom, exonStart, and exonEnd
exDat <- data.frame(chrom     = c(1, 1, 1, 2, 2, 2),
                    exonStart = c(1, 2, 5, 1, 3, 3),
                    exonEnd   = c(3, 4, 7, 4, 5, 6))

exDat

# supply exDat to combine_exons to combine
# overlapping exon segments
combine_exons(exDat)

}
