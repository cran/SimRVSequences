## ------------------------------------------------------------------------
# load the SimRVSequences library
library(SimRVSequences)

# load the hg_exons dataset
data("hg_exons")

# print the first four rows of hg_exons
head(hg_exons, n = 4)

## ------------------------------------------------------------------------
# create recombination map for exon-only data using the hg_exons dataset 
s_map <- create_slimMap(exon_df = hg_exons)

# print first four rows of s_map 
head(s_map, n = 4)

## ------------------------------------------------------------------------
# restrict output to the variables required by SLiM
slimMap <- s_map[, c("recRate", "mutRate", "endPos")]

# shift endPos up by one unit
slimMap$endPos <- slimMap$endPos - 1

# print first four rows of slimMap 
head(slimMap, n = 4)

## ------------------------------------------------------------------------
# load the hg_apopPath data
data("hg_apopPath")

# View the first 4 observations of hg_apopPath
head(hg_apopPath, n = 4)

## ---- eval =  FALSE, echo = TRUE-----------------------------------------
#  # Let's suppose the output is saved in the
#  # current working directory and is named "slimOut.txt".
#  # We import the data using the read_slim function.
#  s_out <- read_slim(file_path  = "slimOut.txt",
#                     recomb_map = create_slimMap(hg_exons),
#                     pathway_df = hg_apopPath)

## ------------------------------------------------------------------------
# import the EXmuts dataset
data(EXmuts)

# view the first 4 observations of EXmuts
head(EXmuts, n = 4)

## ------------------------------------------------------------------------
# import the EXhaps dataset
data(EXhaps)

# dimensions of EXhaps
dim(EXhaps)

## ------------------------------------------------------------------------
# number of rows in EXmuts
nrow(EXmuts)

## ------------------------------------------------------------------------
# View the first 30 mutations of the first 15 haplotypes in EXhaps
EXhaps[1:15, 1:30]

## ---- echo = TRUE, eval = FALSE------------------------------------------
#  # load the vcfR package
#   library(vcfR)

## ---- echo = FALSE-------------------------------------------------------
load(url("https://github.com/simrvprojects/1000-Genomes-Exon-Data/raw/master/Vignette%20Data/vcf_chrom21.rda"))

## ---- eval = FALSE-------------------------------------------------------
#  # specify the file path
#  vcf_file_path <- "C:/Data/exons_chr21.vcf.gz"
#  
#  # import vcf file
#  vcf_chrom21 <- read.vcfR(vcf_file_path)

## ------------------------------------------------------------------------
# View the first five rows and columns of the genotype data
vcf_chrom21@gt[1:5, 1:5]

## ------------------------------------------------------------------------
# View first 4 mutations for the individual with ID "HG00096"
vcf_chrom21@gt[1:4, "HG00096"]

## ------------------------------------------------------------------------
# Remove the variable named "FORMAT" and store the resulting data as genos 
genos <- vcf_chrom21@gt[, -1]

# View the first 5 mutations (i.e. rows) for the first 3 individuals (i.e. columns)
genos[1:5, 1:3]

## ------------------------------------------------------------------------
# load the SimRVSequences package
library(SimRVSequences)

# Convert to sparseMatrix
haplotypes <- genos2sparseMatrix(genotypes = genos)

## ------------------------------------------------------------------------
# View haplotype data for the first 3 diploid individuals and first 5 SNVs
haplotypes[1:6, 1:5]

## ---- eval = FALSE-------------------------------------------------------
#  # View the first two rows of the fixed data
#  vcf_chrom21@fix[1:2, ]

## ---- echo = FALSE-------------------------------------------------------
# View the first two rows of the fixed data
vcf_chrom21@fix[1:2, ]

## ---- eval = FALSE-------------------------------------------------------
#  # View the meta data
#  vcf_chrom21@meta

## ---- echo = FALSE-------------------------------------------------------
# View the meta data
vcf_chrom21@meta

## ---- eval = FALSE-------------------------------------------------------
#  # View the eight item in the meta data
#  vcf_chrom21@meta[[8]]

## ---- echo = FALSE-------------------------------------------------------
# View the eight item in the meta data
vcf_chrom21@meta[[8]]

## ---- echo=FALSE, eval = TRUE--------------------------------------------
load(url("https://github.com/simrvprojects/1000-Genomes-Exon-Data/raw/master/Vignette%20Data/muts.rda"))

## ---- echo = TRUE, eval = FALSE------------------------------------------
#  # Extract and store the mutation data using the vcfR function vcfR2tidy
#  #
#  # NOTE: setting info_only = TRUE since we do not need to re-process the genotype data
#  #
#  # To include INFO variables, we supply a list of variable names to "info_fields".  The
#  # names of these variables must corresponds with the INFO variable names defined in the
#  # meta data. We specify each variable type, by variable name, using the "info_types"
#  # argument.  Each INFO variable type is available in the meta data. Note that the info
#  # types for float varibles are set to numeric, i.e. "n".
#  
#  muts <- vcfR2tidy(vcf_chrom21,
#                    info_only = TRUE,
#                    info_fields = c("AF", "AC", "NS", "AN",
#                                    "EAS_AF", "EUR_AF", "AFR_AF",
#                                    "AMR_AF", "SAS_AF", "DP"),
#                    info_types = c(AF = "n", AC = "i", NS = "i", AN = "i",
#                                   EAS_AF = "n", EUR_AF = "n", AFR_AF = "n",
#                                   AMR_AF = "n", SAS_AF = "n", DP = "i"))

## ------------------------------------------------------------------------
# View a summary of the output
summary(muts)

## ------------------------------------------------------------------------
# View fix (contained in output named muts)
muts$fix

## ------------------------------------------------------------------------
# View the first three rows of the fixed data in muts
head(muts$fix, n = 3)

## ------------------------------------------------------------------------
# Recalculate the alternate allele frequency
muts$fix$AF <- muts$fix$AC/muts$fix$AN

# View the first three rows of the fixed data in muts
head(muts$fix, n = 3)

## ------------------------------------------------------------------------
# store the fixed item in muts data as dataframe 
mutations <- as.data.frame(muts$fix)

## ------------------------------------------------------------------------
# Create the variable colID to identify the column position of the mutation.
# NOTE: Since mutations are already in the correct order, we accomplish this task 
# using the seq function.  Since the mutations are stored in the columns of
# the haplotypes matrix we determine the length of our sequence as the number of 
# columns in haplotypes.
mutations$colID <- seq(1:dim(haplotypes)[2])

# View the first three rows of mutations
head(mutations, n = 3)

## ------------------------------------------------------------------------
# Rename columns for consistency with expected format. 
# "AF" should be renamed "afreq",
# "CHROM" should be renamed "chrom", 
# and "POS" should be renamed "position".
colnames(mutations)[c(1, 2, 8)] = c("chrom", "position", "afreq")

# View the first three rows of mutations
head(mutations, n = 3)

## ---- tidy = TRUE--------------------------------------------------------
file_path <- 'https://raw.githubusercontent.com/simrvprojects/1000-Genomes-Exon-Data/master/Vignette%20Data/SampleInfo1.csv'
SampleInfo <- read.csv(file_path)

## ------------------------------------------------------------------------
# View the first 6 lines of the sample data
head(SampleInfo)

## ------------------------------------------------------------------------
file_path <- 'https://raw.githubusercontent.com/simrvprojects/1000-Genomes-Exon-Data/master/Formatted-SNVdata/SampleData.csv'
SampleData <- read.csv(file_path)

# view the first 6 lines of SampleData
head(SampleData)

## ------------------------------------------------------------------------
# View dimensions of haplotypes
dim(haplotypes)

## ------------------------------------------------------------------------
# Recall that the row names in the haplotypes matrix are the sample IDs
row.names(haplotypes)[1:5]

# Reduce the haplotype data to contain the 
# unrelated individuals described in SampleData 
haplotypes <- haplotypes[row.names(haplotypes) %in% SampleData$Sample, ]

# View the new dimensions of haplotypes
dim(haplotypes)

## ------------------------------------------------------------------------
# create SNVdata object for chromosome 21 without sample or meta data
SNVdata_chrom21 <- SNVdata(Haplotypes = haplotypes, 
                           Mutations = mutations,
                           Samples = SampleData)
str(SNVdata_chrom21)

## ------------------------------------------------------------------------
library(pryr)

# Determine size of vcfR object for chromosome 21
object_size(vcf_chrom21)

# Determine size of SNVdata object for chromosome 21
object_size(SNVdata_chrom21)

## ------------------------------------------------------------------------
# load the hg_apopPath dataset
data("hg_apopPath")

# import SNV data for chromosomes 21 and 22 and identify SNVs located in the 
# pathway defined by hg_apopPath 
EXdata = load_1KG(chrom = 21:22, 
                  pathway_df = hg_apopPath)

# determine the structure of EXdata
str(EXdata)

## ------------------------------------------------------------------------
# View the first 6 SNVs contained in the pathway of interest
head(EXdata$Mutations[EXdata$Mutations$pathwaySNV, ])

## ------------------------------------------------------------------------
# Recall that the variable pathwaySNV, which was created in section 3.2, is TRUE for 
# any SNVs in the pathway of interest. Here we tabulate the alternate allele frequencies 
# of the SNVs located in our pathway.
table(EXmuts$afreq[EXmuts$pathwaySNV == TRUE])

## ------------------------------------------------------------------------
# Create the variable 'is_CRV', which is TRUE for SNVs in our pathway with alternate 
# allele frequency 5e-05 or 1e-04, and FALSE otherwise
EXmuts$is_CRV <- EXmuts$pathwaySNV & EXmuts$afreq %in% c(5e-5, 1e-04)

# verify that sum of the alternate allele frequencies of causal SNVs is 0.001
sum(EXmuts$afreq[EXmuts$is_CRV])

# determine the number of variants in our pool of causal variants
sum(EXmuts$is_CRV)

## ------------------------------------------------------------------------
# view first 4 observations of EXmuts
head(EXmuts, n = 4)

## ---- eval = FALSE-------------------------------------------------------
#  # specify FamID variable
#  FamID = c(314, 314, 314, 314, 314, 314)

## ---- eval = FALSE-------------------------------------------------------
#  # specify ID variable
#  ID = c(1, 2, 3, 4, 5, 6)

## ---- eval = FALSE-------------------------------------------------------
#  # specify sex variable
#  sex = c(0, 1, 0, 1, 0, 1)

## ---- eval = FALSE-------------------------------------------------------
#  # specify dadID variable
#  dadID = c(NA, NA, 1, 1, NA, 5)
#  
#  # specify momID variable
#  momID = c(NA, NA, 2, 2, NA, 4)

## ---- eval = FALSE-------------------------------------------------------
#  # specify affected variable
#  affected = c(FALSE, FALSE, TRUE, TRUE, FALSE, TRUE)

## ---- eval = FALSE-------------------------------------------------------
#  # specify dadID variable
#  DA1 = c(0, 0, 0, 1, 0, 1)
#  
#  # specify momID variable
#  DA2 = c(1, 0, 0, 0, 1, 1)

## ---- eval = FALSE-------------------------------------------------------
#  # combine the required information for family 314 into a dataframe
#  fam314 <- data.frame(FamID = c(314, 314, 314, 314, 314, 314),
#                       ID = c(1, 2, 3, 4, 5, 6),
#                       sex = c(0, 1, 0, 1, 0, 1),
#                       dadID = c(NA, NA, 1, 1, NA, 5),
#                       momID = c(NA, NA, 2, 2, NA, 4),
#                       affected = c(FALSE, FALSE, TRUE, TRUE, FALSE, TRUE),
#                       DA1 = c(0, 0, 0, 1, 0, 1),
#                       DA2 = c(1, 0, 0, 0, 1, 1))

## ---- echo = TRUE, eval = FALSE------------------------------------------
#  # load the SimRVPedigree library
#  library(SimRVPedigree)
#  
#  # Create hazard object from AgeSpecific_Hazards data
#  data(AgeSpecific_Hazards)
#  my_HR = hazard(AgeSpecific_Hazards)
#  
#  # load libraries needed to simulate pedigrees in parallel.
#  library(doParallel)
#  library(doRNG)
#  
#  npeds <- 5    #set the number of pedigrees to generate
#  
#  cl <- makeCluster(2)   # create cluster
#  registerDoParallel(cl) # register cluster
#  
#  
#  # simulate a sample of five pedigrees using foreach
#  study_peds = foreach(i = seq(npeds), .combine = rbind,
#                    .packages = c("SimRVPedigree"),
#                    .options.RNG = 844090518
#  ) %dorng% {
#    # Simulate pedigrees ascertained for at least three disease-affected individuals,
#    # according to the age-specific hazard rates in the `AgeSpecific_Hazards` data
#    # set, ascertained from 1980 to 2010, with start year spanning
#    # from 1900 to 1920, stop year set to 2018, and with genetic relative-risk 50.
#    sim_RVped(hazard_rates = my_HR,
#              GRR = 50, FamID = i,
#              RVfounder = TRUE,
#              founder_byears = c(1900, 1920),
#              ascertain_span = c(1980, 2010),
#              stop_year = 2018,
#              recall_probs = c(1, 0.5, 0),
#              num_affected = 3)[[2]]}
#  
#  stopCluster(cl) #shut down cluster

## ---- fig.height = 8, fig.width = 8.6, fig.align = 'center'--------------
# import study peds
data("study_peds")

# Plot the pedigree with FamID 3
plot(study_peds[study_peds$FamID == 3, ],
     ref_year = 2018)

## ------------------------------------------------------------------------
# View the first 4 rows of study_peds
head(study_peds, n = 4)

## ------------------------------------------------------------------------
# construct an object of class SNVdata
ex_data <- SNVdata(Haplotypes = EXhaps,
                   Mutations = EXmuts)

## ------------------------------------------------------------------------
set.seed(11956)

# simulate SNV data with sim_RVstudy
study_seq <- sim_RVstudy(ped_files = study_peds,
                         SNV_data = ex_data)


# determine the class of study seq
class(study_seq)

## ---- fig.height = 4.5, fig.width = 5.5, fig.align = 'center'------------
# plot the pedigree returned by sim_RVseq for family 3.
# Since we used the affected_only option the pedigree has
# been reduced to contain only disease-affected relatives.
plot(study_seq$ped_files[study_seq$ped_files$FamID == 3, ],
     location = "bottomleft")

## ------------------------------------------------------------------------
# View the first 15 haplotypes in ped_haplos
study_seq$ped_haplos[1:15, ]

## ------------------------------------------------------------------------
# Determine the dimensions of ped_haplos
dim(study_seq$ped_haplos)

# Determine the dimensions of EXhaps
dim(EXhaps)

## ------------------------------------------------------------------------
# View the first 4 observations in SNV_map
head(study_seq$SNV_map, n = 4)

## ------------------------------------------------------------------------
# view the first observations in haplo_map
head(study_seq$haplo_map)

## ------------------------------------------------------------------------
# to quickly view the causal rare variant for 
# each family we supply the appropriate columns of 
# haplo_map to unique
unique(study_seq$haplo_map[, c("FamID", "FamCRV")])

## ------------------------------------------------------------------------
# supply the famStudy object, returned by 
# sim_RVstudy, to the summary function
study_summary <- summary(study_seq)

# determine the class of study_summary
class(study_summary)

# view the names of the items in study_summary
names(study_summary)

## ------------------------------------------------------------------------
# view fam_allele_count
study_summary$fam_allele_count

## ------------------------------------------------------------------------
# view pathway_count
study_summary$pathway_count

## ------------------------------------------------------------------------
# create an example data frame that contains the 
# the variables: chrom, exonStart, and exonEnd
exDat <- data.frame(chrom     = c(1, 1, 1, 2, 2, 2),
                    exonStart = c(1, 2, 5, 1, 3, 3),
                    exonEnd   = c(3, 4, 7, 4, 5, 6))

# View exDat data set
exDat

## ------------------------------------------------------------------------
# supply exDat to combine_exons
# and view the results
combine_exons(exDat)

