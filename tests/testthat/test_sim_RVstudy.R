library(testthat)
library(Matrix)
library(SimRVPedigree)
context("sim_RVstudy")

#from simRVsequences
data("study_peds")
data("EXhaps")
data("EXmuts")

#from simRVPedigree
data("AgeSpecific_Hazards")


#SMALL TOY DATASETS FOR TESTING
toy_haps <- sparseMatrix(i = seq(1:10), j = seq(1:10), x = rep(1, 10))
toy_muts <- data.frame(colID = seq(1:10),
                       chrom = rep(1, 10),
                       position = round(seq(1001, 2000001, length.out = 10)*1000),
                       pathwaySNV = sample(x = c(TRUE, FALSE), size = 10,
                                           replace = TRUE, prob = c(0.2, 0.8)),
                       is_CRV = sample(c(rep(FALSE, 9), TRUE), size = 10))

toy_muts$marker <- paste0(toy_muts$chrom, sep = "_", toy_muts$position)

#---------------------#
# Errors and Warnings #
#---------------------#

test_that("Warning: some pedigrees do not contain affecteds and affected_only = TRUE", {

  red_peds <- study_peds
  red_peds$affected[red_peds$FamID %in% sample(1:5, size = 2)] = FALSE


  expect_warning(sim_RVstudy(ped_files = red_peds,
                             SNV_data = SNVdata(Haplotypes = EXhaps,
                                                Mutations = EXmuts)))
})


test_that("Error: pedigrees do not contain any affecteds and affected_only = TRUE", {

  red_peds <- study_peds
  red_peds$affected = FALSE

  expect_error(sim_RVstudy(ped_files = red_peds,
                           SNV_data = SNVdata(Haplotypes = EXhaps,
                                              Mutations = EXmuts)))
})

test_that("Warning if is_CRV is missing from SNV_map", {

  expect_warning(sim_RVstudy(ped_files = study_peds,
                             SNV_data = SNVdata(Haplotypes = toy_haps,
                                                Mutations = toy_muts[, -5]),
                             remove_wild = FALSE,
                             affected_only = TRUE))
})

test_that("Error: nrow(SNV_map) != ncol(haplos)", {

  expect_error(sim_RVstudy(ped_files = study_peds,
                           SNV_data = SNVdata(Haplotypes = EXhaps,
                                              Mutations = EXmuts[-1, ])))
})

test_that("Error if de novo mutations detected", {

  #choose a family from study peds
  test_fam = sample(x = c(1:5), size = 1)
  red_peds <- study_peds[study_peds$FamID == test_fam, ]
  rownames(red_peds) = NULL

  nonfounder_locs <- which(!is.na(red_peds$dadID) & red_peds$DA1 == 0)
  de_novo_subject <- sample(nonfounder_locs, size = 1)

  #set DA1 and DA2 to 1 for sampled subject
  red_peds$DA1[de_novo_subject] <- 1
  red_peds$DA2[de_novo_subject] <- 1

  expect_error(sim_RVstudy(ped_files = red_peds,
                           SNV_data = SNVdata(Haplotypes = toy_haps,
                                              Mutations = toy_muts),
                           remove_wild = FALSE,
                           affected_only = FALSE))
})


#----------------#
# Output testing #
#----------------#
toy_haps <- sparseMatrix(i = seq(1:100),
                         j = seq(1:100),
                         x = rep(1, 100))
toy_muts <- data.frame(colID = seq(1:100),
                       chrom = rep(1, 100),
                       position = round(seq(1001, 2000001, length.out = 100)*1000),
                       is_CRV = sample(c(rep(FALSE, 90), rep(TRUE, 10)), size = 100))

toy_muts$marker <- paste0(toy_muts$chrom, sep = "_", toy_muts$position)
#toy_muts$afreq <- round(runif(100, min = 0, max = 0.005), digits = 6)
toy_muts$afreq <- 1/100#round(runif(100, min = 0, max = 0.005), digits = 6)

 test_that("When mutiple RV founders in pedigree disease-locus genotypes are correct", {

  #choose a family from study peds
  test_fam = sample(x = 1:5, size = 1)
  red_peds <- study_peds[study_peds$FamID == test_fam, ]
  rownames(red_peds) = NULL

  founder_locs <- which(is.na(red_peds$dadID) & red_peds$DA1 == 0 & red_peds$DA2 == 0)

  #set DA1 to 1 for at least 2 founders
  red_peds$DA1[sample(founder_locs, size = 2)] <- 1

  study_dat = sim_RVstudy(ped_files = red_peds,
                          SNV_data = SNVdata(Haplotypes = toy_haps,
                                             Mutations = toy_muts),
                          remove_wild = TRUE,
                          affected_only = FALSE)

  for(i in 1:nrow(study_dat$ped_files)){
    #pull alleles from pedigree
    ref_alleles <- c(study_dat$ped_files$DA1[i], study_dat$ped_files$DA2[i])

    #pull appropriate rows from ped_haplos to compare
    test_alleles <- study_dat$ped_haplos[study_dat$haplo_map$ID == study_dat$ped_files$ID[i],
                                         study_dat$SNV_map$marker == study_dat$haplo_map$FamCRV[1], ]

    #test for equality
    expect_equal(ref_alleles, test_alleles)
  }
})


test_that("All pedigree members have the correct genotypes at the crv locus", {

  #Simulate a pedigree using simRVPedigree
  RVped <- sim_RVped(hazard_rates = hazard(AgeSpecific_Hazards),
                     GRR = 50, carrier_prob = 0.002,
                     RVfounder = TRUE,
                     FamID = 1,
                     num_affected = 2,
                     recall_probs = c(1),
                     founder_byears = c(1900, 1910),
                     ascertain_span = c(1970, 2015))[[2]]

  #get summary information
  RVsum <- summary(RVped)


  #simulate sequence data for pedigree
  study_dat = sim_RVstudy(ped_files = RVped,
                          SNV_data = SNVdata(Haplotypes = toy_haps,
                                             Mutations = toy_muts),
                          remove_wild = FALSE,
                          affected_only = TRUE)

  for(i in 1:nrow(study_dat$ped_files)){
    #pull alleles from pedigree
    ref_alleles <- c(study_dat$ped_files$DA1[i], study_dat$ped_files$DA2[i])

    #pull appropriate rows from ped_haplos to compare
    test_alleles <- study_dat$ped_haplos[study_dat$haplo_map$ID == study_dat$ped_files$ID[i],
                                         study_dat$SNV_map$marker == study_dat$haplo_map$FamCRV[1], ]

    #test for equality
    expect_equal(ref_alleles, test_alleles)
  }
})


test_that("rows of haplo_map are equal to rows ped_haplos", {

  study_seq <- sim_RVstudy(ped_files = study_peds,
                           SNV_data = SNVdata(Haplotypes = toy_haps,
                                              Mutations = toy_muts),
                           remove_wild = FALSE,
                           affected_only = TRUE)

  expect_equal(nrow(study_seq$ped_haplos), nrow(study_seq$haplo_map))
})



test_that("sproadic families to do carry ANY crvs in pool", {
  #create a fully sporadic sample
  ##Simulate a pedigree using simRVPedigree
  RVped <- sim_RVped(hazard_rates = hazard(AgeSpecific_Hazards),
                     GRR = 50, carrier_prob = 0.002,
                     RVfounder = TRUE,
                     FamID = 1,
                     num_affected = 2,
                     recall_probs = c(1),
                     founder_byears = c(1900, 1910),
                     ascertain_span = c(1970, 2015))[[2]]

  RVped$DA1 <- 0
  RVped$DA2 <- 0

  #sample SNVs to be in the causal pool
  toy_muts2 <- toy_muts
  toy_muts2$is_CRV = FALSE
  while (sum(toy_muts2$is_CRV) == 0) {
    toy_muts2$is_CRV <- sample(x = c(TRUE, FALSE),
                               size = 100,
                               prob = c(0.2, 0.8),
                               replace = TRUE)
  }

  #simulate sequence data for study
  study_seq <- sim_RVstudy(ped_files = RVped,
                           SNV_data = SNVdata(Haplotypes = toy_haps,
                                              Mutations = toy_muts2),
                           remove_wild = FALSE,
                           affected_only = TRUE)

  cRV_colSums <- colSums(study_seq$ped_haplos)[which(toy_muts2$is_CRV == TRUE)]

  expect_equal(cRV_colSums, rep(0, length(cRV_colSums)))

})


test_that("affecteds from genetic families carry the correct number of cRVs at famlial locus", {
  #sample SNVs to be in the causal pool
  toy_muts2 <- toy_muts
  toy_muts2$is_CRV = FALSE
  while (sum(toy_muts2$is_CRV) == 0) {
    toy_muts2$is_CRV <- sample(x = c(TRUE, FALSE),
                               size = 100,
                               prob = c(0.2, 0.8),
                               replace = TRUE)
  }

  #Simulate a pedigree using simRVPedigree
  RVped <- sim_RVped(hazard_rates = hazard(AgeSpecific_Hazards),
                     GRR = 50, carrier_prob = 0.002,
                     RVfounder = TRUE,
                     FamID = 1,
                     num_affected = 2,
                     recall_probs = c(1),
                     founder_byears = c(1900, 1910),
                     ascertain_span = c(1970, 2015))[[2]]


  #get summary information
  num_fam_RVcarriers <- sum(RVped$DA1[RVped$affected], RVped$DA2[RVped$affected])

  #simulate sequence data for study
  study_seq <- sim_RVstudy(ped_files = RVped,
                           SNV_data = SNVdata(Haplotypes = toy_haps,
                                              Mutations = toy_muts2),
                           remove_wild = FALSE,
                           affected_only = TRUE)


  fam_RVcount <- colSums(study_seq$ped_haplos[study_seq$haplo_map$affected, ])[which(toy_muts2$marker == study_seq$haplo_map$FamCRV[1])]

  expect_equal(fam_RVcount, num_fam_RVcarriers)
})
