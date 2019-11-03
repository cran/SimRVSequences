library(testthat)
library(SimRVPedigree)
context("assign_gen")

data("AgeSpecific_Hazards")

#SMALL TOY DATASETS FOR TESTING

for (i in 1:10){
  test_that("assign_gen assigns gen properly in one family", {

    #going to test this function against the true generation numbers
    #simulated by sim_ped,  NOTE: this function does not always produce a
    #full pedigree (may only consist of one person),
    #so will need to incorporate a while statement
    mem1_ped = TRUE
    while (mem1_ped) {
      ex_ped <- sim_ped(hazard_rates = hazard(hazardDF = AgeSpecific_Hazards),
                        RVfounder = TRUE,
                        GRR = 10,
                        FamID = 1,
                        founder_byears = c(1900, 1910),
                        stop_year = 2015)

      if (nrow(ex_ped) > 1){
        mem1_ped = FALSE
      }
    }

    expect_equal(assign_gen(ex_ped), ex_ped$Gen)

  })
}


for (i in 1:5){
  test_that("assign_gen assigns gen properly in multiple families", {

    #generate a list of pedigrees, each with at least 2 generations
    ped_list = list()
    for(k in 1:5){
      #going to test this function against the true generation numbers
      #simulated by sim_ped,  NOTE: this function does not always produce a
      #full pedigree (may only consist of one person),
      #so will need to incorporate a while statement
      mem1_ped = TRUE
      while (mem1_ped) {
        ex_ped <- sim_ped(hazard_rates = hazard(hazardDF = AgeSpecific_Hazards),
                          RVfounder = TRUE,
                          GRR = 10,
                          FamID = k,
                          founder_byears = c(1900, 1910),
                          stop_year = 2015)

        if (nrow(ex_ped) > 1){
          mem1_ped = FALSE
          ped_list[[k]] <- ex_ped
        }
      }
    }

    ex_peds <- do.call(rbind, ped_list)


    expect_equal(unlist(lapply(unique(ex_peds$FamID), function(x){assign_gen(ex_peds[ex_peds$FamID == x, ])})),
                 ex_peds$Gen)

  })
}


test_that("assign_gen assigns gen properly in one family for pedigrees not of class ped", {

  #going to test this function against the true generation numbers
  #simulated by sim_ped,  NOTE: this function does not always produce a
  #full pedigree (may only consist of one person),
  #so will need to incorporate a while statement
  mem1_ped = TRUE
  while (mem1_ped) {
    ex_ped <- sim_ped(hazard_rates = hazard(hazardDF = AgeSpecific_Hazards),
                      RVfounder = TRUE,
                      GRR = 10,
                      FamID = 1,
                      founder_byears = c(1900, 1910),
                      stop_year = 2015)

    if (nrow(ex_ped) > 1){
      mem1_ped = FALSE
      class(ex_ped) = "data.frame"
    }
  }

  expect_equal(assign_gen(ex_ped), ex_ped$Gen)

})

test_that("assign_gen assigns gen properly in mutliple families for pedigrees not of class ped", {

  #generate a list of pedigrees, each with at least 2 generations
  ped_list = list()
  for(k in 1:5){
    #going to test this function against the true generation numbers
    #simulated by sim_ped,  NOTE: this function does not always produce a
    #full pedigree (may only consist of one person),
    #so will need to incorporate a while statement
    mem1_ped = TRUE
    while (mem1_ped) {
      ex_ped <- sim_ped(hazard_rates = hazard(hazardDF = AgeSpecific_Hazards),
                        RVfounder = TRUE,
                        GRR = 10,
                        FamID = k,
                        founder_byears = c(1900, 1910),
                        stop_year = 2015)

      if (nrow(ex_ped) > 1){
        mem1_ped = FALSE
        ped_list[[k]] <- ex_ped
      }
    }
  }

  ex_peds <- do.call(rbind, ped_list)
  class(ex_peds) = "data.frame"


  expect_equal(unlist(lapply(unique(ex_peds$FamID), function(x){assign_gen(ex_peds[ex_peds$FamID == x, ])})),
               ex_peds$Gen)

})
