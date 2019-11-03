library(testthat)
context("get_parOffInfo")

data("study_peds")

test_that("Returns correct dads", {
  red_ped <- study_peds[study_peds$FamID == sample(1:5, size = 1), ]
  poinf <- get_parOffInfo(red_ped)

  expect_equal(red_ped$dadID[which(red_ped$ID %in% poinf$offspring_ID)],
               poinf$parent_ID[poinf$parent == "dadID"])
})


test_that("Returns correct moms", {
  red_ped <- study_peds[study_peds$FamID == sample(1:5, size = 1), ]
  poinf <- get_parOffInfo(red_ped)

  expect_equal(red_ped$momID[which(red_ped$ID %in% poinf$offspring_ID)],
               poinf$parent_ID[poinf$parent == "momID"])
})


test_that("Returns correct offspring RV status for given parent", {
  red_ped <- study_peds[study_peds$FamID == sample(1:5, size = 1), ]
  #plot(red_ped)
  poinf <- get_parOffInfo(red_ped)

  #get offspring IDs to test
  offIDs <- poinf$offspring_ID[poinf$parent == "dadID"]
  #Offspring RVstatus returned by function
  offRVstat_dad <- poinf$Off_RVstatus[poinf$parent == "dadID"]
  offRVstat_mom <- poinf$Off_RVstatus[poinf$parent == "momID"]

  #check inherited from dad
  expect_equal(red_ped[which(red_ped$ID %in% offIDs), c("DA1")],
               offRVstat_dad)

  #check inherited from dad
  expect_equal(red_ped[which(red_ped$ID %in% offIDs), c("DA2")],
               offRVstat_mom)

})


