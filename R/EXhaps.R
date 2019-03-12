#' Example Haplotypes dataset
#'
#' This data set contains 20,000 haplotypes (i.e. rows) spanning 500 single-nucleotide variants (SNVs).  This dataset is intended to accompany the \code{EXmuts} dataset; each row of \code{EXmuts} describes a column (i.e. SNV) in \code{EXhaps}.
#'
#' This dataset is intended to accompany the \code{EXmuts} dataset.  Together, the \code{EXmuts} and \code{EXhaps} datasets represent example output of the \code{read_slim} function.  The \code{EXhaps} data set represents the sparse matrix \code{Haplotypes} returned by \code{read_slim}, and the \code{EXmuts} data set represents the \code{Mutations} data frame returned by \code{read_slim}.  This toy data set, used primarily for demonstration, contains 50 SNVs which were randomly sampled from genes in the apoptosis sub-pathway, and 450 SNVs sampled from outside the pathway.
#'
#'
#' @docType data
#'
#' @seealso \code{\link{EXmuts}}, \code{\link{read_slim}}
#'
#' @format A sparseMatrix of class dgCMatrix with 20000 rows and 500 variables.  Each row represents an observed haplotype and each column represents an SNV locus.
"EXhaps"
