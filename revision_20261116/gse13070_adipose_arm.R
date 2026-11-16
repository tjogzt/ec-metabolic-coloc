# gse13070_adipose_arm.R — R1-11/R0-11：GSE13070 脂肪臂打分（同研究 2×2）
# 对比口径：基线（preTZD、preClamp）脂肪与肌肉；IRd(35/51) vs NRd(6/18)
# 程序分数 = mean(z 化后 UP 基因) − mean(z 化后 DOWN 基因)？不——与主文一致，
# 用 log2 表达量直接 mean(UP) − mean(DOWN)，MW 检验在供体层（每供体一芯片）。
suppressMessages({library(limma)})
setwd("data")

read_series <- function(f) {
  ln <- readLines(gzfile(f))
  i0 <- grep("^!series_matrix_table_begin", ln)
  i1 <- grep("^!series_matrix_table_end", ln)
  pd0 <- grep("^!Sample_title", ln, value = TRUE)
  ttl <- strsplit(pd0[1], "\t")[[1]][-1]
  ttl <- gsub('"', "", ttl)
  tb <- read.delim(text = ln[(i0 + 1):(i1 - 1)], check.names = FALSE, quote = "\"")
  ex <- as.matrix(tb[, -1]); rownames(ex) <- tb[[1]]
  colnames(ex) <- ttl
  list(ex = ex, title = ttl)
}

collapse <- function(ex, annot) {
  a <- annot[match(rownames(ex), annot$ID), c("ID", "Gene Symbol")]
  ok <- !is.na(a[["Gene Symbol"]]) & a[["Gene Symbol"]] != "" & a[["Gene Symbol"]] != "---"
  ex <- ex[ok, , drop = FALSE]; a <- a[ok, , drop = FALSE]
  parts <- strsplit(trimws(as.character(a[["Gene Symbol"]])), "\\s*///\\s*")
  nrep <- lengths(parts)
  if (any(nrep > 1)) {
    ex <- ex[rep(seq_len(nrow(ex)), nrep), , drop = FALSE]
    a  <- a[rep(seq_len(nrow(a)), nrep), , drop = FALSE]
    a[["Gene Symbol"]] <- unlist(parts, use.names = FALSE)
  }
  if (max(ex, na.rm = TRUE) > 50) ex <- log2(ex + 1)
  m <- rowMeans(ex, na.rm = TRUE)
  o <- order(a[["Gene Symbol"]], -m); a <- a[o, ]; ex <- ex[o, , drop = FALSE]
  first <- !duplicated(a[["Gene Symbol"]])
  ex <- ex[first, , drop = FALSE]; rownames(ex) <- a[["Gene Symbol"]][first]
  ex
}

S <- read_series("GSE13070_series_matrix.txt.gz")
annot <- read.delim("GPL570_annot.tsv", check.names = FALSE)
ex <- collapse(S$ex, annot)
ttl <- colnames(ex)

prog <- read.delim("../results/IR_state_program_v2.tsv", row.names = 1)
UP <- rownames(prog)[prog$direction == "UP"]
DN <- rownames(prog)[prog$direction == "DOWN"]
cat(sprintf("程序 UP %d / DN %d；矩阵可打分 UP %d / DN %d\n",
            length(UP), length(DN),
            sum(UP %in% rownames(ex)), sum(DN %in% rownames(ex))))

score <- function(exx) {
  u <- exx[UP[UP %in% rownames(exx)], , drop = FALSE]
  d <- exx[DN[DN %in% rownames(exx)], , drop = FALSE]
  colMeans(u) - colMeans(d)
}

run <- function(pat, label) {
  keep <- grepl(pat, ttl)
  e <- ex[, keep, drop = FALSE]; t <- ttl[keep]
  e <- e[complete.cases(e) & is.finite(rowSums(e)), , drop = FALSE]
  e <- normalizeBetweenArrays(e, method = "quantile")
  grp <- ifelse(grepl("_IRd", t), "IR", ifelse(grepl("_NRd", t), "NR", NA))
  sc <- score(e)
  cat(sprintf("\n=== %s === 样本 %d（IR %d / NR %d）\n", label, length(grp), sum(grp=="IR"), sum(grp=="NR")))
  w <- wilcox.test(sc[grp=="IR"], sc[grp=="NR"])
  cat(sprintf("分数中位 IR %.4f / NR %.4f；差 IR−NR = %+.4f；MW p = %s\n",
              median(sc[grp=="IR"]), median(sc[grp=="NR"]),
              median(sc[grp=="IR"]) - median(sc[grp=="NR"]), format(w$p.value, digits=6)))
  cat(sprintf("均值 IR %.4f / NR %.4f；均值差 %+.4f\n",
              mean(sc[grp=="IR"]), mean(sc[grp=="NR"]),
              mean(sc[grp=="IR"]) - mean(sc[grp=="NR"])))
  # 与主文梯度同口径的"方向一致性"：脂肪臂 IR vs NR 的 DE 与程序方向
  g <- factor(grp, levels = c("NR", "IR"))
  fit <- eBayes(lmFit(e, model.matrix(~ g)), trend = TRUE, robust = TRUE)
  tt <- topTable(fit, coef = "gIR", number = Inf, sort.by = "P")
  ov <- intersect(rownames(prog), rownames(tt))
  same <- sign(prog[ov, "logFC_adipose_obese"]) == sign(tt[ov, "logFC"])
  cat(sprintf("方向一致性：重叠 %d，同向 %d（%.1f%%），二项 p = %.3g\n",
              length(ov), sum(same), 100*mean(same),
              binom.test(sum(same), length(ov), 0.5)$p.value))
  rho <- cor.test(prog[ov, "logFC_adipose_obese"], tt[ov, "logFC"], method = "spearman")
  cat(sprintf("秩相关 rho = %+.3f（p = %.3g）\n", rho$estimate, rho$p.value))
  invisible(list(grp=grp, sc=sc, tt=tt))
}

run("^adipose_preTZD", "GSE13070 脂肪 基线 preTZD（IR 35 vs NR 6）")
run("^muscle_preTZD_preClamp", "GSE13070 肌肉 基线 preTZD_preClamp（IR 51 vs NR 18）")
cat("\n完成 gse13070_adipose_arm\n")

# ---- 精确 p 值（供正文引用） ----
