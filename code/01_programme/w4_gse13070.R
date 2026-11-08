## w4_gse13070.R — A2 补充检验：GSE13070 骨骼肌 基线 IRd+ vs IRd
## 预注册 A2 把 GSE13070 列为发现集 B；实际按 A2 第 2 条的分工，它只承担
## 「跨组织方向一致性」补充检验，不参与程序定义、不影响 T-A/T-B 判定。
## 对比口径：preTZD_preClamp 的肌肉样本（TZD 前、钳夹前基线）。

suppressPackageStartupMessages({
  library(limma); library(data.table)
})

ROOT <- Sys.getenv("ECIR_ROOT", unset = normalizePath("."))
setwd(file.path(ROOT, "data"))

## ---- 读 series matrix（手工解析，GEO 矩阵结构固定） ----
read_series <- function(f) {
  ln <- readLines(gzfile(f), warn = FALSE)
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
  # 多符号条目（"A /// B"）必须拆开，否则伪基因名留在索引、只以多符号形式
  # 注释的基因整体缺席（GPL570 实测丢失 2787 个符号）。
  parts <- strsplit(trimws(as.character(a[["Gene Symbol"]])), "\\s*///\\s*")
  nrep <- lengths(parts)
  if (any(nrep > 1)) {
    ex <- ex[rep(seq_len(nrow(ex)), nrep), , drop = FALSE]
    a  <- a[rep(seq_len(nrow(a)),  nrep), , drop = FALSE]
    a[["Gene Symbol"]] <- unlist(parts, use.names = FALSE)
    cat(sprintf("  [折叠] 多符号条目 %d 行 → 展开 %d 行\n", sum(nrep > 1), sum(nrep)))
  }
  if (max(ex, na.rm = TRUE) > 50) {
    cat(sprintf("  [预处理] 原始强度（max %.0f）→ log2(x+1)\n", max(ex, na.rm = TRUE)))
    ex <- log2(ex + 1)
  }
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

## ---- 选基线肌肉样本：muscle + preTZD + preClamp ----
keep <- grepl("^muscle_preTZD_preClamp", ttl)
ex <- ex[, keep, drop = FALSE]; ttl <- ttl[keep]

## 分组（标签已于 characteristics 字段核实）：
##   IRd / IRd+ / IRd- 三者的 "insulin resistance" 字段**都写 IRd（insulin resistant）**，
##   +/- 后缀指的是 "response to tzd treatment"（应答者/不应答者），不是抵抗/敏感！
##   真正的慢性 IR 对比是 IRd（合计 51）vs NRd（insulin sensitive，18）。
grp <- ifelse(grepl("_IRd", ttl), "IR",
       ifelse(grepl("_NRd", ttl), "NR", NA))
stopifnot(!any(is.na(grp)))
cat(sprintf("基因 %d，样本 %d（抵抗 %d / 敏感 %d）\n",
            nrow(ex), length(grp), sum(grp == "IR"), sum(grp == "NR")))

# 顺序：先去 NA 行再做分位数归一化（NA 会污染 normalizeBetweenArrays）
ex <- ex[complete.cases(ex) & is.finite(rowSums(ex)), , drop = FALSE]
ex <- normalizeBetweenArrays(ex, method = "quantile")
g <- factor(grp, levels = c("NR", "IR"))

fit <- eBayes(lmFit(ex, model.matrix(~ g)), trend = TRUE, robust = TRUE)
tt <- topTable(fit, coef = "gIR", number = Inf, sort.by = "P")
write.table(tt, "../DE_GSE13070_muscle_baseline_IR.tsv", sep = "\t", quote = FALSE)

sig <- tt[tt$adj.P.Val < 0.05 & abs(tt$logFC) >= 0.3, ]
cat(sprintf("\nFDR<0.05 且 |logFC|>=0.3：%d\n", nrow(sig)))
cat(sprintf("最小 FDR %.4f；名义 p<0.001：%d（期望 %.0f）\n",
            min(tt$adj.P.Val), sum(tt$P.Value < 0.001), nrow(tt) * 0.001))
cat(sprintf("logFC 中位绝对值 %.3f\n", median(abs(tt$logFC))))
# tt 按 P 排序：head 是最显著、tail 是最不显著，不能当"上调/下调"
up5 <- rownames(tt)[order(-tt$logFC)][1:5]
dn5 <- rownames(tt)[order(tt$logFC)][1:5]
cat(sprintf("logFC 最高 5：%s\n", paste(up5, collapse = ", ")))
cat(sprintf("logFC 最低 5：%s\n", paste(dn5, collapse = ", ")))

## ---- 与 v2 程序的方向一致性（A2 补充检验） ----
prog <- read.delim("../results/IR_state_program_v2.tsv", row.names = 1)
ov <- intersect(rownames(prog), rownames(tt))
same <- sign(prog[ov, "logFC_adipose_obese"]) == sign(tt[ov, "logFC"])
cat(sprintf("\n=== 与 v2 程序（%d 基因）的方向一致性 ===\n", nrow(prog)))
cat(sprintf("重叠 %d，同向 %d（%.1f%%），二项 p = %.3g\n",
            length(ov), sum(same), 100 * mean(same),
            binom.test(sum(same), length(ov), 0.5)$p.value))
cat("完成 GSE13070\n")
