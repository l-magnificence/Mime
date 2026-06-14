# Full pipeline integration test with example data
library(testthat)
library(Mime1)

cat("Loading example data...\n")
load("D:/桌面/mime/Mime/External data/Example.cohort.Rdata")
load("D:/桌面/mime/Mime/External data/genelist.Rdata")

cat("Dataset names:", paste(names(list_train_vali_Data), collapse=", "), "\n")
cat("Training dim:", dim(list_train_vali_Data[[1]]), "\n")
cat("Genelist length:", length(genelist), "\n")

# ============================================================
# Test 1: ML.Dev.Prog.Sig mode='single' with RSF
# ============================================================
test_that("ML.Dev.Prog.Sig single mode RSF works", {
  cat("\n--- Test: ML.Dev.Prog.Sig single RSF ---\n")
  expect_no_error({
    res_single <<- ML.Dev.Prog.Sig(
      train_data = list_train_vali_Data$Dataset1,
      list_train_vali_Data = list_train_vali_Data,
      unicox.filter.for.candi = TRUE,
      unicox_p_cutoff = 0.05,
      candidate_genes = genelist,
      mode = "single",
      single_ml = "RSF",
      nodesize = 5,
      seed = 5201314
    )
  })
  expect_true(!is.null(res_single$ml.res))
  expect_true(!is.null(res_single$riskscore))
  expect_true(!is.null(res_single$Sig.genes))
  cat("PASS: ML.Dev.Prog.Sig single RSF\n")
})

# ============================================================
# Test 2: cal_AUC_ml_res with single RSF result
# ============================================================
test_that("cal_AUC_ml_res works with RSF result", {
  cat("\n--- Test: cal_AUC_ml_res single RSF ---\n")
  expect_no_error({
    auc_res <<- cal_AUC_ml_res(
      res.by.ML.Dev.Prog.Sig = res_single,
      train_data = list_train_vali_Data$Dataset1,
      inputmatrix.list = list_train_vali_Data,
      mode = "single",
      single_ml = "RSF",
      AUC_time = 1,
      auc_cal_method = "KM"
    )
  })
  expect_true(is.list(auc_res))
  # Check AUC values are in [0, 1]
  for (ds_name in names(auc_res)) {
    auc_vals <- auc_res[[ds_name]]$AUC
    valid_aucs <- auc_vals[!is.na(auc_vals)]
    if (length(valid_aucs) > 0) {
      expect_true(all(valid_aucs >= 0 & valid_aucs <= 1),
                  info = paste0("AUC out of range in ", ds_name))
    }
  }
  cat("PASS: cal_AUC_ml_res single RSF\n")
})

# ============================================================
# Test 3: ML.Dev.Prog.Sig mode='all' (full pipeline)
# ============================================================
test_that("ML.Dev.Prog.Sig all mode works", {
  cat("\n--- Test: ML.Dev.Prog.Sig mode='all' ---\n")
  expect_no_error({
    res_all <<- ML.Dev.Prog.Sig(
      train_data = list_train_vali_Data$Dataset1,
      list_train_vali_Data = list_train_vali_Data,
      unicox.filter.for.candi = TRUE,
      unicox_p_cutoff = 0.05,
      candidate_genes = genelist,
      mode = "all",
      nodesize = 5,
      seed = 5201314
    )
  })
  expect_true(length(res_all$ml.res) > 0)
  expect_true(length(res_all$riskscore) > 0)
  cat("Models trained:", length(res_all$ml.res), "\n")
  cat("Models:", paste(names(res_all$ml.res), collapse=", "), "\n")
  cat("PASS: ML.Dev.Prog.Sig all mode\n")
})

# ============================================================
# Test 4: cal_AUC_ml_res with all mode result
# ============================================================
test_that("cal_AUC_ml_res all mode works and AUC in range", {
  cat("\n--- Test: cal_AUC_ml_res all mode ---\n")
  expect_no_error({
    auc_all <<- cal_AUC_ml_res(
      res.by.ML.Dev.Prog.Sig = res_all,
      train_data = list_train_vali_Data$Dataset1,
      inputmatrix.list = list_train_vali_Data,
      mode = "all",
      AUC_time = 1,
      auc_cal_method = "KM"
    )
  })
  # Check all AUC values in [0, 1]
  for (model_name in names(auc_all)) {
    for (ds_name in names(auc_all[[model_name]])) {
      auc_vals <- auc_all[[model_name]][[ds_name]]$AUC
      valid_aucs <- auc_vals[!is.na(auc_vals)]
      if (length(valid_aucs) > 0) {
        expect_true(all(valid_aucs >= 0 & valid_aucs <= 1),
                    info = paste0("AUC > 1 in ", model_name, "/", ds_name))
      }
    }
  }
  cat("PASS: cal_AUC_ml_res all mode\n")
})

# ============================================================
# Test 5: Risk scores have no Inf/NaN
# ============================================================
test_that("risk scores contain no Inf or NaN", {
  cat("\n--- Test: risk score integrity ---\n")
  for (model_name in names(res_all$riskscore)) {
    for (ds_name in names(res_all$riskscore[[model_name]])) {
      rs <- res_all$riskscore[[model_name]][[ds_name]]$RS
      expect_true(!any(is.infinite(rs)),
                  info = paste0("Inf in RS: ", model_name, "/", ds_name))
      expect_true(!any(is.nan(rs)),
                  info = paste0("NaN in RS: ", model_name, "/", ds_name))
    }
  }
  cat("PASS: no Inf/NaN in risk scores\n")
})

# ============================================================
# Test 6: C-index values are valid
# ============================================================
test_that("C-index values are in valid range", {
  cat("\n--- Test: C-index validity ---\n")
  cc <- res_all$Cindex.res
  expect_true(all(cc$Cindex >= 0 & cc$Cindex <= 1, na.rm = TRUE),
              info = "C-index out of [0, 1]")
  cat("PASS: C-index valid\n")
})

cat("\n===================================\n")
cat("ALL PIPELINE TESTS PASSED\n")
cat("===================================\n")
