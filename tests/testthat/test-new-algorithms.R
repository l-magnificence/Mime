test_that("ML.Dev.Prog.Sig mode=all includes new algorithms", {
  # Skip if example data not available
  skip_if(!file.exists("External data/Example.cohort.Rdata"))

  load("External data/Example.cohort.Rdata")

  # Run with a small seed and all algorithms
  res <- ML.Dev.Prog.Sig(
    train_data = Example.cohort[[1]],
    list_train_vali_Data = Example.cohort,
    mode = "all",
    seed = 42,
    ensemble_top_n = 5,
    algo_timeout = 300
  )

  # Check that at least some new algorithms appear in results
  model_names <- unique(res$Cindex.res$Model)
  new_algos <- c("XGBoost-Surv", "MCP", "SCAD", "Alasso", "ORSF", "CIF", "mboost", "XGB-Cox")
  found_new <- intersect(model_names, new_algos)
  expect_true(length(found_new) > 0,
    info = paste("Expected at least one new algorithm, found:", paste(found_new, collapse=", ")))

  # Check ensemble results exist
  expect_true(!is.null(res$Ensemble_riskscore),
    info = "Ensemble_riskscore should be present when ensemble_top_n > 0")

  # Check backward compatibility: existing keys present
  expect_true("Cindex.res" %in% names(res))
  expect_true("ml.res" %in% names(res))
  expect_true("riskscore" %in% names(res))
  expect_true("Sig.genes" %in% names(res))
})
