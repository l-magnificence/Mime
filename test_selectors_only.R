#!/usr/bin/env Rscript
# 快速测试新选择器 (StabSel 和 MRMR 替代实现)

cat("====== 快速测试新选择器 ======\n\n")

# 设置用户库路径
user_lib <- file.path(Sys.getenv("USERPROFILE"), "R", "win-library",
                       paste0(R.version$major, ".", as.integer(R.version$minor)))
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
.libPaths(c(user_lib, .libPaths()))

# 加载 Mime 包
pkg_path <- "D:/test/Mime-main/Mime-main"
r_files <- list.files(file.path(pkg_path, "R"), pattern = "\\.R$", full.names = TRUE, recursive = FALSE)
for (f in r_files) {
  tryCatch(source(f, local = FALSE), error = function(e) {})
}

# 兼容层
if (!exists("var.select", envir = .GlobalEnv)) {
  var.select <- function(object, conservative = "high", ...) {
    ms <- max.subtree(object, conservative = (conservative == "high"), ...)
    depth <- ms$topvars
    if (is.null(depth)) depth <- names(sort(ms$order[, 1]))
    list(topvars = depth)
  }
  environment(var.select) <- .GlobalEnv
  assign("var.select", var.select, envir = .GlobalEnv)
}

library(survival)

# 加载数据
cat("[1] 加载数据...\n")
data_path <- "D:/桌面/大二下/combined_ComBat_corrected.csv"
raw_data <- read.csv(data_path, row.names = 1, check.names = FALSE)
raw_data$ID <- rownames(raw_data)
raw_data <- raw_data[, c("ID", "OS.time", "OS", setdiff(colnames(raw_data), c("ID", "OS.time", "OS")))]

# 特征预筛选
gene_vars <- apply(raw_data[, -(1:3)], 2, var, na.rm = TRUE)
gene_vars <- gene_vars[!is.na(gene_vars) & gene_vars > 0]
top_var_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(2000, length(gene_vars))]

cox_pvals <- numeric(length(top_var_genes))
names(cox_pvals) <- top_var_genes
for (i in seq_along(top_var_genes)) {
  cox_pvals[i] <- tryCatch({
    fit <- coxph(Surv(OS.time, OS) ~ ., data = raw_data[, c("OS.time", "OS", top_var_genes[i])])
    summary(fit)$coefficients[, "Pr(>|z|)"]
  }, error = function(e) 1)
}

sig_genes <- names(cox_pvals[cox_pvals < 0.001 & !is.na(cox_pvals)])
if (length(sig_genes) < 50) sig_genes <- names(cox_pvals[cox_pvals < 0.01 & !is.na(cox_pvals)])
if (length(sig_genes) < 50) sig_genes <- names(cox_pvals[cox_pvals < 0.05 & !is.na(cox_pvals)])
if (length(sig_genes) < 20) sig_genes <- names(sort(cox_pvals[!is.na(cox_pvals)]))[1:min(200, length(cox_pvals))]
if (length(sig_genes) > 200) sig_genes <- names(sort(cox_pvals[sig_genes]))[1:200]

raw_data <- raw_data[, c("ID", "OS.time", "OS", sig_genes)]
cat("  最终基因数:", length(sig_genes), "\n")

# 划分训练/验证集
set.seed(42)
n <- nrow(raw_data)
train_idx <- sample(1:n, size = floor(0.7 * n))
train_data <- raw_data[train_idx, ]
vali_data <- raw_data[-train_idx, ]
rownames(train_data) <- NULL
rownames(vali_data) <- NULL

list_train_vali_Data <- list(Train = train_data, Validation = vali_data)

# 测试单个新算法
cat("\n[2] 测试 mode='single' 新算法...\n")
new_algos <- c("Alasso", "CIF", "mboost", "MCP", "SCAD")
for (algo in new_algos) {
  cat(sprintf("  测试 %s ...", algo))
  tryCatch({
    res_i <- ML.Dev.Prog.Sig(
      train_data = train_data,
      list_train_vali_Data = list_train_vali_Data,
      mode = "single",
      single_ml = algo,
      seed = 42,
      algo_timeout = 120,
      cores_for_parallel = 1
    )
    if (!is.null(res_i$Cindex.res) && nrow(res_i$Cindex.res) > 0) {
      cat(sprintf(" OK C-index: %.3f / %.3f\n",
          res_i$Cindex.res$Cindex[res_i$Cindex.res$ID == "Train"],
          res_i$Cindex.res$Cindex[res_i$Cindex.res$ID == "Validation"]))
    } else {
      cat(" 无结果\n")
    }
  }, error = function(e) cat(sprintf(" 失败: %s\n", e$message)))
  gc()
}

# 测试 StabSel 和 MRMR 选择器（通过 mode="double" 模拟）
cat("\n[3] 测试 ENet Stability 选择器 (原 StabSel)...\n")
tryCatch({
  # 使用 mode="all" 但只测试 StabSel 部分 - 通过检查输出
  # 直接测试选择器本身
  est_dd <- train_data
  set.seed(42)

  cat("  运行 ENet Stability 选择器...\n")
  x_mat <- as.matrix(est_dd[, -c(1, 2)])
  y_surv <- Surv(est_dd$OS.time, est_dd$OS)
  gene_freq <- rep(0, ncol(x_mat))
  names(gene_freq) <- colnames(x_mat)
  n_fits <- 0
  t1 <- Sys.time()
  for (aa in seq(0.1, 1.0, by = 0.1)) {
    for (rr in 1:3) {
      fit_e <- tryCatch({
        idx_s <- sample(1:nrow(x_mat), floor(0.8 * nrow(x_mat)), replace = FALSE)
        cv.glmnet(x_mat[idx_s, ], y_surv[idx_s], family = "cox", alpha = aa, nfolds = 5)
      }, error = function(e) NULL)
      if (!is.null(fit_e)) {
        coefs <- coef(fit_e, s = "lambda.min")
        sel <- rownames(coefs)[which(coefs[, 1] != 0)]
        if (length(sel) > 0) gene_freq[sel] <- gene_freq[sel] + 1
        n_fits <- n_fits + 1
      }
    }
  }
  t2 <- Sys.time()
  stab_genes <- if (n_fits > 0) names(which(gene_freq / n_fits >= 0.3)) else character(0)
  cat(sprintf("  完成! %d次拟合, 选中 %d 个基因, 耗时 %.1f 秒\n", n_fits, length(stab_genes), as.numeric(difftime(t2, t1, units = "secs"))))
  cat("  选中基因:", paste(head(stab_genes, 10), collapse = ", "), "\n")
}, error = function(e) cat("  StabSel 失败:", e$message, "\n"))

cat("\n[4] 测试 UniCindex 选择器 (原 MRMR)...\n")
tryCatch({
  set.seed(42)
  t1 <- Sys.time()
  genes <- colnames(est_dd[, -c(1, 2)])
  cidx_v <- sapply(genes, function(g) {
    tryCatch({
      fit_c <- coxph(Surv(OS.time, OS) ~ ., data = est_dd[, c('OS.time', 'OS', g)])
      summary(fit_c)$concordance[1]
    }, error = function(e) NA_real_)
  })
  cidx_v <- cidx_v[!is.na(cidx_v)]
  cidx_sel <- sort(cidx_v, decreasing = TRUE)[1:min(20, length(cidx_v))]
  t2 <- Sys.time()
  cat(sprintf("  完成! 选中 %d 个基因, 耗时 %.1f 秒\n", length(cidx_sel), as.numeric(difftime(t2, t1, units = "secs"))))
  cat("  Top-10 基因 C-index:\n")
  print(head(cidx_sel, 10))
}, error = function(e) cat("  UniCindex 失败:", e$message, "\n"))

cat("\n====== 快速测试完成 ======\n")
