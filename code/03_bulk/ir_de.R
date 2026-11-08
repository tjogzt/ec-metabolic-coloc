suppressMessages({library(GEOquery); library(limma)})
setwd("data")
options(timeout = 900)

## ---------- 工具 ----------
read_series <- function(f) {
  es <- getGEO(filename = f, getGPL = FALSE)
  ex <- exprs(es); pd <- pData(es)
  colnames(ex) <- pd[colnames(ex), "title"]     # 列名统一用样本标题
  list(ex = ex, pd = pd)
}

collapse_to_gene <- function(ex, annot, idcol = "ID", symcol = NULL) {
  if (is.null(symcol) || !(symcol %in% colnames(annot)))
    symcol <- c("Gene Symbol", "GENE_SYMBOL", "GENE")[c("Gene Symbol", "GENE_SYMBOL", "GENE") %in% colnames(annot)][1]
  a <- annot[match(rownames(ex), annot[[idcol]]), c(idcol, symcol)]
  keep <- !is.na(a[[symcol]]) & a[[symcol]] != "" & a[[symcol]] != "---"
  ex <- ex[keep, , drop = FALSE]; a <- a[keep, , drop = FALSE]
  # 每基因取平均表达最高的探针
  m <- rowMeans(ex, na.rm = TRUE)
  ord <- order(a[[symcol]], -m)
  a <- a[ord, ]; ex <- ex[ord, , drop = FALSE]
  first <- !duplicated(a[[symcol]])
  ex <- ex[first, , drop = FALSE]
  rownames(ex) <- a[[symcol]][first]
  ex
}

## ---------- A. GSE24215 骨骼肌：胰岛素 vs 基础（t0，配对） ----------
A <- read_series("GSE24215_series_matrix.txt.gz")
a_annot <- read.delim("GPL6480_annot.tsv", check.names = FALSE)
exA <- collapse_to_gene(A$ex, a_annot)
cat("GSE24215 基因数:", nrow(exA), " 样本:", ncol(exA), "\n")
cat("表达值范围:", round(range(exA, na.rm = TRUE), 2), "\n")

titles <- colnames(exA)
is_t0 <- grepl("_t0_", titles)
is_base <- grepl("basal", titles)
is_ins <- grepl("insulin", titles)
is_tech <- grepl("techrep", titles)
sel <- is_t0 & !is_tech
exA0 <- exA[, sel, drop = FALSE]
cond <- ifelse(grepl("basal", colnames(exA0)), "basal", "insulin")
subj <- paste0("s", sub(".*rep([0-9]+).*", "\\1", colnames(exA0)))
cat("t0 用于配对比较：", sum(cond == "basal"), "基础 vs", sum(cond == "insulin"), "胰岛素\n")
print(table(subj, cond))

design <- model.matrix(~ subj + cond)
fit <- eBayes(lmFit(exA0, design), trend = TRUE, robust = TRUE)
ttA <- topTable(fit, coef = "condinsulin", number = Inf, sort.by = "P")
write.table(ttA, "../DE_GSE24215_insulin_vs_basal.tsv", sep = "\t", quote = FALSE)
cat("GSE24215 DE：FDR<0.05 的基因", sum(ttA$adj.P.Val < 0.05),
    "；其中 |logFC|>=0.3 的", sum(ttA$adj.P.Val < 0.05 & abs(ttA$logFC) >= 0.3), "\n")

## ---------- B. GSE24422 脂肪细胞：胰岛素方向确认 ----------
B <- read_series("GSE24422_series_matrix.txt.gz")
b_annot <- read.delim("GPL570_annot.tsv", check.names = FALSE)
exB <- collapse_to_gene(B$ex, b_annot)
cat("\nGSE24422 基因数:", nrow(exB), " 样本:", ncol(exB), "\n")

tb <- colnames(exB)
grp <- ifelse(grepl("^Adip_Ins", tb), "ins",
       ifelse(grepl("^Adip_rep", tb), "ctl", NA))
selB <- !is.na(grp)
exB2 <- exB[, selB, drop = FALSE]; gB <- factor(grp[selB], levels = c("ctl", "ins"))
cat("脂肪细胞组分：对照", sum(gB == "ctl"), "胰岛素", sum(gB == "ins"), "\n")
d2 <- model.matrix(~ gB)
fitB <- eBayes(lmFit(exB2, d2), trend = TRUE, robust = TRUE)
ttB <- topTable(fitB, coef = "gBins", number = Inf, sort.by = "P")
write.table(ttB, "../DE_GSE24422_insulin_adipocyte.tsv", sep = "\t", quote = FALSE)
cat("GSE24422：FDR<0.05 的基因", sum(ttB$adj.P.Val < 0.05), "\n")

## ---------- C. 方向一致性 ----------
common <- intersect(rownames(ttA), rownames(ttB))
cat("\n两数据集共有基因:", length(common), "\n")
ag <- ttA[common, ]; bg <- ttB[common, ]
sigA <- ag$adj.P.Val < 0.05 & abs(ag$logFC) >= 0.3
cat("A 中显著且 |logFC|>=0.3：", sum(sigA), "\n")
if (sum(sigA) > 0) {
  s <- sign(ag$logFC[sigA]); t2 <- sign(bg$logFC[sigA])
  agree <- mean(s == t2, na.rm = TRUE)
  cat("B 中方向一致比例：", round(100 * agree, 1), "%（二项检验 p=",
      signif(binom.test(sum(s == t2, na.rm = TRUE), sum(!is.na(t2)), 0.5)$p.value, 3), "）\n")
  cat("B 中同向且 FDR<0.05 的：",
      sum(s == t2 & bg$adj.P.Val[sigA] < 0.05, na.rm = TRUE), "\n")
}
saveRDS(list(ttA = ttA, ttB = ttB, exA = exA, exB = exB),
        "../ir_de_cache.rds")
cat("\n完成。\n")
