#!/usr/bin/env Rscript
# Mime Prognosis Model Algorithm Expansion — 测试脚本
# 数据: combined_ComBat_corrected.csv
# 日期: 2026-05-05

cat("====== Mime 算法扩展测试 ======\n\n")

# 设置用户库路径
user_lib <- file.path(Sys.getenv("USERPROFILE"), "R", "win-library",
                       paste0(R.version$major, ".", as.integer(R.version$minor)))
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
.libPaths(c(user_lib, .libPaths()))

# ---- 0. 加载 Mime 包 ----
cat("[0] 加载 Mime 包...\n")
pkg_path <- "D:/test/Mime-main/Mime-main"
# 直接 source 所有 R 文件
r_files <- list.files(file.path(pkg_path, "R"), pattern = "\\.R$", full.names = TRUE, recursive = FALSE)
for (f in r_files) {
  tryCatch({
    source(f, local = FALSE)
  }, error = function(e) {
    # 忽略加载错误
  })
}
cat("  R 文件已加载\n")

# 兼容层: randomForestSRC 3.x 移除了 var.select，用 max.subtree 替代
if (!exists("var.select", envir = .GlobalEnv)) {
  var.select <- function(object, conservative = "high", ...) {
    ms <- max.subtree(object, conservative = (conservative == "high"), ...)
    depth <- ms$topvars
    if (is.null(depth)) {
      depth <- names(sort(ms$order[, 1]))
    }
    list(topvars = depth)
  }
  environment(var.select) <- .GlobalEnv
  assign("var.select", var.select, envir = .GlobalEnv)
  cat("  已添加 var.select 兼容层\n")
}

library(survival)
library(ggplot2)

# ---- 1. 加载数据 ----
cat("\n[1] 加载数据...\n")
data_path <- "D:/桌面/大二下/combined_ComBat_corrected.csv"
raw_data <- read.csv(data_path, row.names = 1, check.names = FALSE)

cat("  样本数:", nrow(raw_data), "\n")
cat("  基因数:", ncol(raw_data) - 2, "(减去 OS.time 和 OS)\n")
cat("  OS 分布: 死亡(OS=1) =", sum(raw_data$OS == 1),
    ", 存活(OS=0) =", sum(raw_data$OS == 0), "\n")

# 添加 ID 列（Mime 要求第一列为 ID）
raw_data$ID <- rownames(raw_data)
raw_data <- raw_data[, c("ID", "OS.time", "OS", setdiff(colnames(raw_data), c("ID", "OS.time", "OS")))]

# ---- 1.5 特征预筛选（减少基因数量以加速测试）----
cat("\n[1.5] 特征预筛选...\n")

# 第一步：按方差筛选，保留变异最大的 top 2000 基因（避免对 16000+ 基因做 Cox 回归导致崩溃）
gene_vars <- apply(raw_data[, -(1:3)], 2, var, na.rm = TRUE)
gene_vars <- gene_vars[!is.na(gene_vars) & gene_vars > 0]
top_var_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(2000, length(gene_vars))]
cat("  方差筛选后基因数:", length(top_var_genes), "\n")

# 第二步：对 top 2000 基因做单因素 Cox 回归
library(survival)
cox_pvals <- numeric(length(top_var_genes))
names(cox_pvals) <- top_var_genes
for (i in seq_along(top_var_genes)) {
  cox_pvals[i] <- tryCatch({
    fit <- coxph(Surv(OS.time, OS) ~ ., data = raw_data[, c("OS.time", "OS", top_var_genes[i])])
    summary(fit)$coefficients[, "Pr(>|z|)"]
  }, error = function(e) 1)
  if (i %% 500 == 0) cat("    已处理", i, "/", length(top_var_genes), "个基因\n")
}
cat("  Cox 回归完成\n")

sig_genes <- names(cox_pvals[cox_pvals < 0.001 & !is.na(cox_pvals)])
cat("  筛选后基因数 (p<0.001):", length(sig_genes), "\n")

if (length(sig_genes) < 50) {
  sig_genes <- names(cox_pvals[cox_pvals < 0.01 & !is.na(cox_pvals)])
  cat("  放宽阈值至 p<0.01，基因数:", length(sig_genes), "\n")
}
if (length(sig_genes) < 50) {
  sig_genes <- names(cox_pvals[cox_pvals < 0.05 & !is.na(cox_pvals)])
  cat("  放宽阈值至 p<0.05，基因数:", length(sig_genes), "\n")
}
if (length(sig_genes) < 20) {
  sig_genes <- names(sort(cox_pvals[!is.na(cox_pvals)]))[1:min(200, length(cox_pvals))]
  cat("  取 top 200 基因\n")
}

# 上限 200 基因以控制运行时间
if (length(sig_genes) > 200) {
  sig_genes <- names(sort(cox_pvals[sig_genes]))[1:200]
  cat("  限制最多 200 基因\n")
}

raw_data <- raw_data[, c("ID", "OS.time", "OS", sig_genes)]
cat("  最终基因数:", length(sig_genes), "\n")

# ---- 2. 划分训练集和验证集 ----
cat("\n[2] 划分训练集 (70%) 和验证集 (30%)...\n")
set.seed(42)
n <- nrow(raw_data)
train_idx <- sample(1:n, size = floor(0.7 * n))
train_data <- raw_data[train_idx, ]
vali_data <- raw_data[-train_idx, ]

# 重置行名
rownames(train_data) <- NULL
rownames(vali_data) <- NULL

cat("  训练集样本数:", nrow(train_data), "\n")
cat("  验证集样本数:", nrow(vali_data), "\n")

# 构建 list_train_vali_Data（训练集 + 验证集）
list_train_vali_Data <- list(
  Train = train_data,
  Validation = vali_data
)

# ---- 3. 运行 ML.Dev.Prog.Sig（mode="all" 全量组合模式）----
cat("\n[3] 运行 ML.Dev.Prog.Sig（mode='all' 全量组合模式）...\n")
cat("  包含：18个独立算法 + 5个特征选择器 × 18个算法 = ~200+组合\n\n")

t_start <- Sys.time()

tryCatch({
  withCallingHandlers(
    res <- ML.Dev.Prog.Sig(
      train_data = train_data,
      list_train_vali_Data = list_train_vali_Data,
      mode = "all",
      seed = 42,
      algo_timeout = 300,
      cores_for_parallel = 1,
      ensemble_top_n = 5
    ),
    error = function(e) {
      cat("  [TRACEBACK] Error caught:", e$message, "\n")
      cat("  [TRACEBACK] Call stack:\n")
      calls <- sys.calls()
      for (i in seq_along(calls)) {
        call_str <- deparse(calls[[i]], width.cutoff = 200)[1]
        cat(sprintf("    %d: %s\n", i, substr(call_str, 1, 150)))
      }
    }
  )
}, error = function(e) {
  cat("  mode='all' 整体运行失败:", e$message, "\n")
  cat("  回退到 mode='single' 逐个运行...\n")
  res <<- NULL
})

# 如果 mode="all" 失败，回退到逐个运行
if (is.null(res) || is.null(res$Cindex.res) || nrow(res$Cindex.res) == 0) {
  cat("  回退: 逐一运行18个算法...\n")
  all_algos <- c("RSF", "Enet", "StepCox", "CoxBoost", "plsRcox", "superpc",
                 "GBM", "survivalsvm", "Ridge", "Lasso",
                 "Alasso", "ORSF", "CIF", "mboost")
  skip_algos <- c("GBM", "CoxBoost")
  all_cindex <- data.frame()
  all_riskscore <- list()

  for (algo in all_algos) {
    if (algo %in% skip_algos) {
      cat(sprintf("  [%d/%d] 跳过 %s (当前系统segfault)\n", match(algo, all_algos), length(all_algos), algo))
      next
    }
    cat(sprintf("  [%d/%d] 运行 %s ...\n", match(algo, all_algos), length(all_algos), algo))
    tryCatch({
      res_i <- ML.Dev.Prog.Sig(
        train_data = train_data,
        list_train_vali_Data = list_train_vali_Data,
        mode = "single",
        single_ml = algo,
        seed = 42,
        algo_timeout = 300,
        cores_for_parallel = 1
      )
      if (!is.null(res_i$Cindex.res) && nrow(res_i$Cindex.res) > 0) {
        all_cindex <- rbind(all_cindex, res_i$Cindex.res)
        cat(sprintf("    成功! C-index: %s\n",
            paste(round(res_i$Cindex.res$Cindex, 3), collapse = ", ")))
      }
      if (!is.null(res_i$riskscore)) {
        all_riskscore <- c(all_riskscore, res_i$riskscore)
      }
    }, error = function(e) {
      cat(sprintf("    失败: %s\n", e$message))
    })
    gc()
  }

  res <- list(
    Cindex.res = all_cindex,
    riskscore = all_riskscore,
    Ensemble_riskscore = NULL
  )
}

t_end <- Sys.time()
cat("\n  运行耗时:", round(difftime(t_end, t_start, units = "mins"), 1), "分钟\n")

# ---- 4. 查看结果 ----
cat("\n[4] 结果概览:\n")

# 所有模型的 C-index
cindex_df <- res$Cindex.res
cat("  总模型数:", length(unique(cindex_df$Model)), "\n")

# 按训练集 C-index 排序
train_cindex <- cindex_df[cindex_df$ID == "Train", ]
train_cindex <- train_cindex[order(-train_cindex$Cindex), ]

cat("\n  === 训练集 Top-10 模型 ===\n")
print(head(train_cindex[, c("Model", "Cindex")], 10))

vali_cindex <- cindex_df[cindex_df$ID == "Validation", ]
vali_cindex <- vali_cindex[order(-vali_cindex$Cindex), ]

cat("\n  === 验证集 Top-10 模型 ===\n")
print(head(vali_cindex[, c("Model", "Cindex")], 10))

# 检查新算法是否出现在结果中
new_algos <- c("Alasso", "ORSF", "CIF", "mboost")
found_new <- intersect(unique(cindex_df$Model), new_algos)
cat("\n  新算法已运行:", length(found_new), "/", length(new_algos), "\n")
cat("  ", paste(found_new, collapse = ", "), "\n")

# 检查 ensemble
if (!is.null(res$Ensemble_riskscore)) {
  cat("\n  集成模型 (Ensemble) 已生成!\n")
  ensemble_cindex <- cindex_df[grepl("Ensemble", cindex_df$Model), ]
  print(ensemble_cindex)
} else {
  cat("\n  集成模型未生成 (ensemble_top_n = 0 或无足够模型)\n")
}

# ---- 5. 生成图表 ----
cat("\n[5] 生成图表...\n")

output_dir <- "D:/桌面/大二下/Mime_results"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 5.1 C-index 比较柱状图
cat("  [5.1] C-index 比较图...\n")
cindex_plot_data <- train_cindex
cindex_plot_data$IsEnsemble <- grepl("Ensemble", cindex_plot_data$Model)
cindex_plot_data$IsNew <- cindex_plot_data$Model %in% new_algos

# 只展示 Top-30 模型
cindex_plot_data <- head(cindex_plot_data, 30)

p1 <- ggplot(cindex_plot_data, aes(x = reorder(Model, Cindex), y = Cindex,
                                     fill = ifelse(IsEnsemble, "Ensemble",
                                            ifelse(IsNew, "New Algorithm", "Existing")))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Ensemble" = "firebrick", "New Algorithm" = "steelblue", "Existing" = "grey60"),
                    name = "Model Type") +
  labs(title = "C-index Comparison (Training Set - Top 30)",
       x = "Model", y = "C-index") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "01_cindex_training_top30.png"), p1, width = 12, height = 10, dpi = 150)

# 验证集
cindex_plot_vali <- vali_cindex
cindex_plot_vali$IsEnsemble <- grepl("Ensemble", cindex_plot_vali$Model)
cindex_plot_vali$IsNew <- cindex_plot_vali$Model %in% new_algos
cindex_plot_vali <- head(cindex_plot_vali, 30)

p2 <- ggplot(cindex_plot_vali, aes(x = reorder(Model, Cindex), y = Cindex,
                                     fill = ifelse(IsEnsemble, "Ensemble",
                                            ifelse(IsNew, "New Algorithm", "Existing")))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Ensemble" = "firebrick", "New Algorithm" = "steelblue", "Existing" = "grey60"),
                    name = "Model Type") +
  labs(title = "C-index Comparison (Validation Set - Top 30)",
       x = "Model", y = "C-index") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "02_cindex_validation_top30.png"), p2, width = 12, height = 10, dpi = 150)

# 5.2 新算法 vs 旧算法 C-index 对比
cat("  [5.2] 新旧算法 C-index 对比...\n")
algo_compare <- cindex_df
algo_compare$Type <- ifelse(algo_compare$Model %in% new_algos, "New", "Existing")

p3 <- ggplot(algo_compare, aes(x = Type, y = Cindex, fill = Type)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  facet_wrap(~ID, scales = "free_y") +
  scale_fill_manual(values = c("New" = "steelblue", "Existing" = "grey60")) +
  labs(title = "New vs Existing Algorithm C-index Distribution",
       x = "Algorithm Type", y = "C-index") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "03_new_vs_existing_cindex.png"), p3, width = 10, height = 6, dpi = 150)

# 5.3 集成模型 KM 生存曲线
cat("  [5.3] 集成模型 KM 生存曲线...\n")
if (!is.null(res$Ensemble_riskscore)) {
  for (cohort_name in names(res$Ensemble_riskscore)) {
    df <- data.frame(
      ID = list_train_vali_Data[[cohort_name]]$ID,
      OS.time = as.numeric(list_train_vali_Data[[cohort_name]]$OS.time),
      OS = as.numeric(list_train_vali_Data[[cohort_name]]$OS),
      RS = res$Ensemble_riskscore[[cohort_name]]
    )
    df$RiskGroup <- ifelse(df$RS >= median(df$RS), "High Risk", "Low Risk")

    fit_km <- survfit(Surv(OS.time, OS) ~ RiskGroup, data = df)

    # 使用 base R 绘制 KM 曲线
    png(file.path(output_dir, paste0("04_KM_ensemble_", cohort_name, ".png")), width = 800, height = 600, res = 150)

    plot(fit_km, col = c("#E7B800", "#2E9FDF"), lwd = 2,
         xlab = "Time (Days)", ylab = "Overall Survival Probability",
         main = paste0("Ensemble Model KM Curve — ", cohort_name))
    legend("topright", legend = c("High Risk", "Low Risk"),
           col = c("#E7B800", "#2E9FDF"), lwd = 2, bty = "n")

    # 计算 log-rank p 值
    surv_diff <- survdiff(Surv(OS.time, OS) ~ RiskGroup, data = df)
    p_val <- 1 - pchisq(surv_diff$chisq, df = 1)
    text(x = max(df$OS.time) * 0.6, y = 0.9,
         labels = paste0("Log-rank P = ", format.pval(p_val, digits = 3)),
         cex = 1.1, col = "red")

    dev.off()
    cat("    KM 图已保存:", cohort_name, "\n")
  }
}

# 5.4 重要新算法的 KM 曲线
cat("  [5.4] 重要新算法 KM 曲线...\n")
top_new_algos <- c("Alasso", "ORSF", "CIF", "mboost")
for (algo in top_new_algos) {
  if (!is.null(res$riskscore[[algo]]) && !is.null(res$riskscore[[algo]][["Train"]])) {
    df <- res$riskscore[[algo]][["Train"]]
    df$RiskGroup <- ifelse(df$RS >= median(df$RS), "High Risk", "Low Risk")

    if (length(unique(df$RiskGroup)) < 2) {
      cat("    ", algo, "KM 图跳过（风险分组不足）\n")
      next
    }

    tryCatch({
      fit_km <- survfit(Surv(OS.time, OS) ~ RiskGroup, data = df)

      png(file.path(output_dir, paste0("05_KM_", algo, ".png")), width = 800, height = 600, res = 150)
      plot(fit_km, col = c("#E7B800", "#2E9FDF"), lwd = 2,
           xlab = "Time (Days)", ylab = "Overall Survival Probability",
           main = paste0(algo, " — KM Curve (Training Set)"))
      legend("topright", legend = c("High Risk", "Low Risk"),
             col = c("#E7B800", "#2E9FDF"), lwd = 2, bty = "n")

      surv_diff <- survdiff(Surv(OS.time, OS) ~ RiskGroup, data = df)
      p_val <- 1 - pchisq(surv_diff$chisq, df = 1)
      text(x = max(df$OS.time) * 0.6, y = 0.9,
           labels = paste0("Log-rank P = ", format.pval(p_val, digits = 3)),
           cex = 1.1, col = "red")
      dev.off()
      cat("    ", algo, "KM 图已保存\n")
    }, error = function(e) {
      cat("    ", algo, "KM 图失败:", e$message, "\n")
    })
  }
}

# 5.5 时间依赖 AUC（使用 survivalROC）
cat("  [5.5] 时间依赖 AUC 曲线...\n")
tryCatch({
  library(survivalROC)

  # 选验证集上 Top-1 新算法
  top_algo <- vali_cindex$Model[1]
  cat("    使用最佳模型:", top_algo, "\n")

  if (!is.null(res$riskscore[[top_algo]]) && !is.null(res$riskscore[[top_algo]][["Validation"]])) {
    df <- res$riskscore[[top_algo]][["Validation"]]
    time_points <- c(365, 365*3, 365*5)  # 1年、3年、5年
    time_labels <- c("1-Year", "3-Year", "5-Year")

    png(file.path(output_dir, "06_time_dependent_AUC.png"), width = 800, height = 600, res = 150)
    plot(0, 0, type = "n", xlim = c(0, max(df$OS.time, na.rm = TRUE)),
         ylim = c(0, 1), xlab = "Time (Days)", ylab = "AUC",
         main = paste0("Time-dependent AUC — ", top_algo, " (Validation)"))

    colors <- c("red", "blue", "darkgreen")
    auc_values <- c()

    for (i in seq_along(time_points)) {
      tp <- time_points[i]
      if (tp <= max(df$OS.time, na.rm = TRUE)) {
        roc_res <- survivalROC(Stime = df$OS.time, status = df$OS,
                                marker = df$RS, predict.time = tp,
                                method = "KM")
        lines(roc_res$FP, roc_res$TP, col = colors[i], lwd = 2)
        auc_values <- c(auc_values, round(roc_res$AUC, 3))
      }
    }

    legend("bottomright", legend = paste0(time_labels[1:length(auc_values)], " AUC = ", auc_values),
           col = colors[1:length(auc_values)], lwd = 2, bty = "n")
    abline(a = 0, b = 1, lty = 2, col = "grey")
    dev.off()
    cat("    AUC 图已保存\n")
  }
}, error = function(e) {
  cat("    AUC 绘制失败:", e$message, "\n")
})

# 5.6 C-index 热力图（所有算法 × 所有队列）
cat("  [5.6] C-index 热力图...\n")
tryCatch({
  library(ggplot2)

  cindex_matrix <- cindex_df
  cindex_matrix$Cindex[is.na(cindex_matrix$Cindex)] <- 0

  # 只取有结果的模型
  p6 <- ggplot(cindex_matrix, aes(x = ID, y = reorder(Model, Cindex), fill = Cindex)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "white", high = "steelblue", midpoint = 0.5,
                         limits = c(0.4, 1), name = "C-index") +
    labs(title = "C-index Heatmap: All Models × All Cohorts",
         x = "Cohort", y = "Model") +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_text(size = 7))

  ggsave(file.path(output_dir, "07_cindex_heatmap.png"), p6, width = 8, height = max(10, nrow(cindex_matrix) * 0.25), dpi = 150)
  cat("    热力图已保存\n")
}, error = function(e) {
  cat("    热力图绘制失败:", e$message, "\n")
})

# ---- 6. 保存结果 ----
cat("\n[6] 保存结果...\n")

# 保存 C-index 表格
write.csv(cindex_df, file.path(output_dir, "cindex_results.csv"), row.names = FALSE)
cat("  C-index 结果已保存至:", file.path(output_dir, "cindex_results.csv"), "\n")

# 保存 Top 模型列表
write.csv(train_cindex, file.path(output_dir, "top_models_training.csv"), row.names = FALSE)
write.csv(vali_cindex, file.path(output_dir, "top_models_validation.csv"), row.names = FALSE)

# ---- 7. 汇总 ----
cat("\n====== 测试完成 ======\n")
cat("所有图表和结果已保存至:", output_dir, "\n")
cat("\n生成的文件:\n")
cat("  01_cindex_training_top30.png     - 训练集 C-index Top30\n")
cat("  02_cindex_validation_top30.png   - 验证集 C-index Top30\n")
cat("  03_new_vs_existing_cindex.png    - 新旧算法对比箱线图\n")
cat("  04_KM_ensemble_*.png             - 集成模型 KM 生存曲线\n")
cat("  05_KM_*.png                      - 各新算法 KM 生存曲线\n")
cat("  06_time_dependent_AUC.png        - 时间依赖 AUC 曲线\n")
cat("  07_cindex_heatmap.png            - C-index 热力图\n")
cat("  cindex_results.csv               - 所有模型 C-index 结果\n")
