% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SlimFunctions.R
\name{create_slimMap}
\alias{create_slimMap}
\title{Create recombination map}
\usage{
create_slimMap(exon_df, mutation_rate = 1e-08, recomb_rate = 1e-08)
}
\arguments{
\item{exon_df}{Data frame. A data frame that contains the positions of each exon to simulate.  This data frame must contain the variables \code{chrom}, \code{exonStart}, and \code{exonEnd}.  See details.}

\item{mutation_rate}{Numeric.  The per-site per-generation mutation rate, assumed to be constant across the genome. By default, \code{mutation_rate= 1E-8}, as in Harris and Nielson (2016).}

\item{recomb_rate}{Numeric.  The per-site per-generation mutation rate, assumed to be constant across the genome. By default, \code{mutation_rate= 1E-8}, as in Harris and Nielson (2016)}
}
\value{
A recombination map that may be used in conjunction with SLiM (Haller and Messer 2017).  See details and example.
}
\description{
Create a recombination map that can be used with SLiM (Haller and Messer 2017).
}
\details{
The Eidos program SLiM (Haller and Messer 2017) is a versatile forwards-in-time evolutionary simulator.  SLiM simulates recombination hotspots by way of a user-specified recombination map.  This recombination map may be utilized to simulate mutations over unlinked regions (i.e. in different chromosomes) or in linked but non-contiguous regions (i.e in exon-only data).  The \code{create_slimMap} function may be used to generate the recombination map required by SLiM to simulate exon-only SNV data.

We expect that \code{exon_df} does not contain any overlapping segments.  Prior to supplying the exon data to \code{create_slimMap} users must combine overlapping exons into a single observation.  The \code{\link{combine_exons}} function may be used to accomplish this task.

The argument \code{exon_df} must contain the following variables:
\tabular{lll}{
\strong{name} \tab \strong{type} \tab \strong{description} \cr
\code{chrom} \tab numeric \tab chromosome identification number\cr
\code{exonStart} \tab numeric \tab the position of the first base pair in the exon\cr
\code{exonStop} \tab numeric \tab the position of the last base pair in the exon\cr
}

The data frame returned by \code{create_slimMap} contains variables required by SLiM to simulate exon-only data.  Additionally, the returned data frame also includes variables that are required to re-map mutations to their correct positions when importing SLiM data to \code{R}.  The variables contained in the returned data frame are described as follows.
\describe{
\item{\code{chrom}}{The chromosome number.}
\item{\code{segLength}}{The length of the segment in base pairs.  We assume that segments contain the positions listed in \code{exonStart} and \code{exonEnd}.  Therefore, for a combined exon segment, \code{segLength} is calculated as \code{exonEnd - exonStart + 1}.}
\item{\code{recRate}}{The per-site per-generation recombination rate.  Following Harris and Nielson (2016), segments between exons on the same chromosome are simulated as a single base pair with \code{rec_rate} equal to recombination rate multiplied by the number of base pairs in the segment.  For each chromosome, a single site is created between the last exon on the previous chromosome and the first exon of the current chromosome.  This site will have recombination rate 0.5 to accommodate unlinked chromosomes.}
\item{\code{mutRate}}{The per-site per-generation mutation rate.  Since we are interested in exon-only data, the mutation rate outside exons is set to zero.}
\item{\code{exon}}{A logical variable that is \code{TRUE} if the segment is an exon and \code{FALSE} otherwise.}
\item{\code{simDist}}{The simulated exon length, in base pairs. When \code{exon = TRUE}, \code{simDist = segLength}; however, when \code{exon = FALSE}, \code{simDist = 1} since segments between exons on the same chromosome are simulated as a single base pair.}
\item{\code{endPos}}{The simulated end position, in base pairs, of the segment.}
}

Only three of the variables returned by \code{create_slimMap} are required by SLiM to simulate exon-only data: \code{recRate}, \code{mutRate}, and \code{endPos}.  The other variables seen in the output above are used by the \code{\link{read_slim}} function to re-map mutations to their correct positions when importing SLiM data to \code{R}.

Please note: SLiM is written in a scripting language called Eidos. Unlike an \code{R} array, the first position in an Eidos array is 0.  Therefore, users must shift the variable \code{endPos} forward 1 unit before supplying this variable to SLiM. See example.
}
\examples{
#load hg_exons data
data(hg_exons)

#since the exons in hg_exons have already been combined into
#overlapping exons, we supply hg_exons to create_slimMap
slimMap <- create_slimMap(hg_exons)
head(slimMap)

# restrict output to the variables required by SLiM
slimMap <- slimMap[, c("recRate", "mutRate", "endPos")]

# shift endPos up by one unit
slimMap$endPos <- slimMap$endPos - 1

# print first four rows of slimMap
head(slimMap, n = 4)

}
\references{
Benjamin Haller and Phillip W. Messer (2017). \emph{Slim 2: Flexible, interactive forward genetic simulations}. Molecular Biology and Evolution; 34(1), pp. 230-240.

Kelly Harris and Rasmus Nielsen (2016). \emph{The genetic cost of neanderthal introgression}. Genetics, 203(2): pp. 881-891.
}
\seealso{
\code{\link{combine_exons}}
}
