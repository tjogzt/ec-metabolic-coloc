suppressMessages({library(GEOquery); library(limma)})
setwd("data")

read_series <- function(f) {
  es <- getGEO(filename = f, getGPL = FALSE)
  ex <- exprs(es); pd <- pData(es)
  colnames(ex) <- pd$title
  list(ex = ex, pd = pd)
}

collapse <- function(ex, annot, symcol = NULL) {
  if (is.null(symcol) || !(symcol %in% colnames(annot)))
    symcol <- c("Gene Symbol", "GENE_SYMBOL", "GENE")[
      c("Gene Symbol", "GENE_SYMBOL", "GENE") %in% colnames(annot)][1]
  a <- annot[match(rownames(ex), annot$ID), c("ID", symcol)]
  ok <- !is.na(a[[symcol]]) & a[[symcol]] != "" & a[[symcol]] != "---"
  ex <- ex[ok, , drop = FALSE]; a <- a[ok, , drop = FALSE]
  m <- rowMeans(ex, na.rm = TRUE)
  o <- order(a[[symcol]], -m); a <- a[o, ]; ex <- ex[o, , drop = FALSE]
  first <- !duplicated(a[[symcol]])
  ex <- ex[first, , drop = FALSE]; rownames(ex) <- a[[symcol]][first]
  ex
}

## ================= 发现队列：GSE17025（91 肿瘤 vs 11 正常内膜） =================
S <- read_series("GSE17025_series_matrix.txt.gz")
annot570 <- read.delim("GPL570_annot.tsv", check.names = FALSE)
ex <- collapse(S$ex, annot570)
cat("GSE17025 基因", nrow(ex), " 样本", ncol(ex), "\n")

hist <- sapply(strsplit(as.character(S$pd$characteristics_ch1), ": "), function(x) x[2])
# characteristics 的第 5 列是 histology
hcol <- grep("histology", as.character(S$pd[1, ]), value = FALSE)
histv <- as.character(S$pd[[grep("^histology", colnames(S$pd))[1]]])
histv <- sub("histology: ", "", histv)
cat("组织学分布：\n"); print(table(histv))

grp <- ifelse(grepl("endometrioid|papillary serous", histv), "tumor",
       ifelse(grepl("endometrium", histv), "normal", NA))
cat("分组："); print(table(grp))
keep <- !is.na(grp)
ex <- ex[, keep, drop = FALSE]; grp <- factor(grp[keep], levels = c("normal", "tumor"))
ex <- normalizeBetweenArrays(ex, method = "quantile")

design <- model.matrix(~ grp)
fit <- eBayes(lmFit(ex, design), trend = TRUE, robust = TRUE)
tt <- topTable(fit, coef = "grptumor", number = Inf, sort.by = "P")
write.table(tt, "../DE_GSE17025_tumor_vs_normal.tsv", sep = "\t", quote = FALSE)
cat("\n差异基因：FDR<0.05 =", sum(tt$adj.P.Val < 0.05),
    "；FDR<0.05 且 |logFC|>=1 =", sum(tt$adj.P.Val < 0.05 & abs(tt$logFC) >= 1), "\n")
cat("上调前 10：", paste(head(rownames(tt)[tt$logFC > 0], 10), collapse = ", "), "\n")
cat("下调前 10：", paste(head(rownames(tt)[tt$logFC < 0], 10), collapse = ", "), "\n")

## ================= 复现队列：GSE196033（57 癌 vs 30 内膜，769 基因面板） =================
T <- read_series("GSE196033_series_matrix.txt.gz")
exT <- T$ex
tv <- as.character(T$pd[[grep("^tisue", colnames(T$pd), ignore.case = TRUE)[1]]])
tv <- sub("^[^:]*: ", "", tv)
cat("\nGSE196033 分组："); print(table(tv))
gT <- factor(ifelse(grepl("cancer", tv), "tumor", "normal"), levels = c("normal", "tumor"))
exT <- log2(exT + 1); exT <- normalizeBetweenArrays(exT, method = "quantile")
fitT <- eBayes(lmFit(exT, model.matrix(~ gT)), trend = TRUE, robust = TRUE)
ttT <- topTable(fitT, coef = "gTtumor", number = Inf, sort.by = "P")
write.table(ttT, "../DE_GSE196033_tumor_vs_normal.tsv", sep = "\t", quote = FALSE)

## ================= 方向一致性 =================
sig <- rownames(tt)[tt$adj.P.Val < 0.05 & abs(tt$logFC) >= 1]
ov <- intersect(sig, rownames(ttT))
cat("\n发现签名", length(sig), "个基因；与 GSE196033 面板重叠", length(ov), "个\n")
if (length(ov) > 10) {
  ag <- sign(tt[ov, "logFC"]); bg <- sign(ttT[ov, "logFC"])
  cat("方向一致率：", round(100 * mean(ag == bg), 1), "%",
      "（二项 p =", signif(binom.test(sum(ag == bg), length(ov), 0.5)$p.value, 3), "）\n")
  cat("同向且在复现队列 FDR<0.05 的：", sum(ag == bg & ttT[ov, "adj.P.Val"] < 0.05), "\n")
}
write.table(data.frame(gene = sig,
                       direction = ifelse(tt[sig, "logFC"] > 0, "UP", "DOWN"),
                       logFC = tt[sig, "logFC"], FDR = tt[sig, "adj.P.Val"]),
            "../EC_signature_v1.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("\n已写 ../EC_signature_v1.tsv\n")
