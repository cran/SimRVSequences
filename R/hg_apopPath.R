#' Apoptosis sub-pathway dataset
#'
#' This data set catalogs combined exon segments from the 25 genes that have the highest interaction with the \emph{TNFSF10} gene, a known member of the human apoptosis pathway.
#'
#' The \code{hg_apopPath} data set catalogs the positions of exons residing in the 25 genes that have the highest interaction with the \emph{TNFSF10} gene.  The data contained in the \code{hg_apopPath} data set was collected from the hg 38 reference genome with the UCSC Genome Browser.  The 25 genes that have the highest interaction with the \emph{TNFSF10} gene were identified by the UCSC Genome Browser's Gene Interaction Tool.  In \code{hg_apopPath} overlapping exons have been combined into a single observation. When exons from genes with different NCBI accession numbers have been combined the variable \code{NCBIref} will contain multiple accession numbers separated by commas.  We note that different accession numbers may exist for transcript variants of the same gene.
#'
#' @docType data
#'
#' @references Karolchik, D., Hinrichs, A. S., Furey, T. S., Roskin, K. M., Sugnet, C. W., Haussler, D., and Ken, W. J. (2004). The UCSC Table Browser data retrieval tool. \emph{Nucleic Acids Res}. Accessed on 20 February 2018.
#' @references Kent, W. J., Sugnet, C. W., Furey, T. S., Roskin, K. M., Pringle, T. H., Zahler, A. M., and Haussler, D. (2002). The human genome browser at UCSC. \emph{Genome Res}, 12(6):996-1006.
#' @references Poon, H., Quirk, C., DeZiel, C., and Heckerman, D. (2014). Literome: Pubmed-scale genomic knowledge base in the cloud. \emph{Bioinformatics}, 30:2840-2842. Accessed on 20 February 2018.
#'
#' @format A data set with 253 rows and 5 variables:
#' \describe{
#'   \item{chrom}{Numeric. The chromosome number.}
#'   \item{exonStart}{Numeric. The position of the first base pair.}
#'   \item{exonStop}{Numeric. The position of the last base pair.}
#'   \item{NCBIref}{Character. The NCBI reference sequence accession number of the gene(s) in which the exon(s) reside.}
#'   \item{gene}{Character. The name of the gene.}
#' }
"hg_apopPath"
