#' Constructor function for an object of class SNVdata
#'
#' @param Haplotypes sparseMatrix. A sparse matrix of haplotype data, which contains the haplotypes for unrelated individuals representing the founder population.  Rows are assumed to be haplotypes, while columns represent SNVs.
#' @param Mutations Data frame. A data frame that catalogs the SNVs in \code{Haplotypes}.
#' @param Samples An optional dataframe or matrix describing the individuals in \code{Haplotypes}.
#'
#' @return an object of class \code{SNVdata}.
#' @export
SNVdata <- function(Haplotypes, Mutations, Samples = NULL) {

  #check SNV_map for possible issues
  check_SNV_map(Mutations)

  if (!"marker" %in% colnames(Mutations)) {
    Mutations$marker <- make.unique(paste0(Mutations$chrom, sep = "_", Mutations$position))
  }

  if (nrow(Mutations) != ncol(Haplotypes)) {
    stop("\n nrow(Mutations) != ncol(Haplotypes). \n Mutations must catalog every SNV in Haplotypes.")
  }



  #create list containing all relavant of SNVdata information
  SNV_data = list(Haplotypes = Haplotypes,
                  Mutations = Mutations,
                  Samples = Samples)

  class(SNV_data) <- c("SNVdata", class(SNV_data))
  return(SNV_data)
}

#' Check to see if object is of class ped
#'
#' @param x An R object.
#'
#' @return Logical. Indicates if \code{x} is of class \code{ped}.
#'
#' @keywords internal
is.SNVdata <- function(x) {
  return(inherits(x, "SNVdata"))
}

#' Load pre-formatted 1000 Genomes Project exon data
#'
#' Load pre-formatted 1000 Genomes Project exon data
#'
#' The \code{load_1KG} is used to load pre-formatted, exon-only SNV data from any of the 22 human autosomes.  The original data was obtained from:
#'
#' http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/.
#'
#' The data was reduced to remove any related indiviuals, to accopmlish this we randomly sampled one relative from each set of related individuals.  This resulted in the removal of 22 individuals.  Additional information regaring the formatting of the 1000 Genomes Project data may be found at https://github.com/simrvprojects/1000-Genomes-Exon-Data/ in the pdf file entitled "Documentation for Creating Exon Data_090319.pdf".
#'
#' We expect that \code{pathwayDF} does not contain any overlapping segments.  Users may combine overlapping exons into a single observation with the \code{\link{combine_exons}} function.
#'
#' @param chrom Numeric.  The chromosome number(s).  A numeric list of chromosome numbers representing the 1000 Genomes Project exon-data to load.
#' @param pathway_df Data frame. (Optional) A data frame that contains the positions for each exon in a pathway of interest.  This data frame must contain the variables \code{chrom}, \code{exonStart}, and \code{exonEnd}. See Details.
#'
#' @return An object of class \code{SNVdata} containing the imported exon data.
#' @export
#' @importFrom utils read.csv
#' @seealso \code{\link{combine_exons}}
#'
#' @references 1000 Genomes Project (2010). \emph{A Map of Human Genome Variation from Population-Scale Sequencing}. Nature; 467:1061-1073.
#'
#' @examples
#' exdata = load_1KG(21:22)
#' unique(exdata$Mutations$chrom)
#'
#' head(exdata$Mutations)
#' exdata$Haplotypes[1:20, 1:10]
#' head(exdata$Samples)
load_1KG <- function(chrom, pathway_df = NULL){

  if (any(!chrom %in% seq(1:22))){
    stop("\n We expect 'chrom' to be a numeric list of automosomes.
        Note: We do not provide sex chromosomes (i.e. chromosomes X and Y)")
  }

  #store the object names for the imported data by chromosome number
  object_names <- paste0("SNVdata_chrom", chrom)


  #import the formatted date from github
  for (i in 1:length(chrom)){
    load(url(paste0("https://github.com/simrvprojects/1000-Genomes-Exon-Data/raw/master/Formatted-SNVdata/SNVdata_chrom",
                    chrom[i], ".rda", sep = "")))
  }

  #store each to SNVdata object to a list
  chrom_dat = lapply(object_names, function(x){get(x)})

  if (length(chrom) > 1){
    #combine the data from the different chromosomes, for ease of use
    Haplotypes <- do.call(cbind, lapply(chrom_dat, `[[`, 1))
    Mutations <- do.call(rbind, lapply(chrom_dat, `[[`, 2))
    #reformat colID
    Mutations$colID <- c(1:dim(Haplotypes)[2])
  } else {
    Haplotypes <- chrom_dat[[1]]$Haplotypes
    Mutations <- chrom_dat[[1]]$Mutations
  }

  #import sample data from GitHub
  SampleData <- read.csv(url("https://github.com/simrvprojects/1000-Genomes-Exon-Data/raw/master/Formatted-SNVdata/SampleData.csv"),
                         stringsAsFactors = FALSE)

  #----------------------#
  # Identify Pathway RVs #
  #----------------------#
  #if pathway data has been supplied, identify pathway SNVs
  if (!is.null(pathway_df)) {
    Mutations <- identify_pathwaySNVs(markerDF = Mutations,
                                      pathwayDF = pathway_df)
  }

  return(SNVdata(Haplotypes = Haplotypes, Mutations = Mutations, Samples = SampleData))
}
