master.code <- function() {
  rm(list = ls())
  source("ImmRes_source.R")
  ######################################################################################
  print("1. Generating de-novo cell-type signatures.")
  # The code is provided in "ImmRes1_denovoCellTypeSig.R"
  cell.type.de <- get.cell.type.sig()

  ######################################################################################
  print("2. Identifying immune resistance programs in malignant melanoma cells.")
  # The code is provided in "ImmRes2_immuneResistanceProgram.R"
  r.sc <- readRDS("../Data/scData/Mel.malignant.rds")
  r.tcga <- readRDS("../Data/PublicData/TCGA_SKCM.rds")
  load("../Results/CellTypes/cell.type.sig.RData")
  load("../Results/Resistance/PostTreatment/post.treatment.subsample.de.RData")

  print("2.1. Generating the T cell exclusion signatures.")
  exc.de <- mal.t.cell.exclusion(rB = r.tcga, r.sc = r.sc, cell.sig = cell.sig)
  print("2.2. Generating the post-treatment signatures.")
  trt.de <- get.post.trt.sig(r = r.sc, subs.trt = subs.trt)
  print("2.3. Generating the functional signatures.")
  fnc.de <- get.fnc.res.sig(r = r.sc)
  print("2.4. Combining the signatures into the immune resistance program.")
  res.sig <- get.res.program()

  ######################################################################################
  print("3. Longitudinal analysis.")
  # The code is provided in "ImmRes3_longitudinal.R"
  print("3.1 Analyzing Validation Cohort 1.")
  # valCo1<-set.ValCo1()
  valCo1.results <- test.ValCo1(r = valCo1)
  print("3.2 Analyzing a MAPKi resistance cohort (Hugo et al. 2015).")
  mapkiCo <- set.matched.MAPKi.Hugo()
  mapki.results <- test.matched.MAPKi.Hugo(r = mapkiCo)

  ######################################################################################
  print("4. Predicting OS and ICB responses in published cohorts.")
  # The code is provided in "ImmRes4_predictICBresponses.R"
  r.tcga <- set.TCGA(r.tcga = r.tcga)
  R <- set.public.ICB.cohorts()

  ######################################################################################
  print("5. Predicting ICB responses in Validation Cohort #2.")
  # The code is provided in "ImmRes5_valCohort2.R"
  print("5.1 Use different signatures to predict ICB responses.")
  # r.pd1<-set.aPD1()
  print("5.2 Testing ICB response predictions.")
  aPD1.val <- prd.aPD1()
  generate.fig5() # Provided in "ImmRes_output.R"

  ######################################################################################
  print("6. Prioritizing drugs to repress the immune resistance program.")
  find.sensitizing.drugs()

  ######################################################################################
  print("7. Exploring the impact of CDK4/6 inhibition on melanoma cells.")
  print("See the interactive data and plots in the single-cell portal")
  print("To login use the email account icr.review1@gmail.com, password icrreview2.")
  sc.url <- "https://portals.broadinstitute.org/single_cell/study/melanoma-immunotherapy-resistance#study-visualize"
  browseURL(sc.url, browser = getOption("browser"), encodeIfNeeded = FALSE)

  return()
}
