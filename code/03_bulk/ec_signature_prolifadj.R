suppressMessages({library(GEOquery); library(limma)})
setwd("data")

read_series <- function(f) {
  es <- getGEO(filename = f, getGPL = FALSE)
  ex <- exprs(es); pd <- pData(es); colnames(ex) <- pd$title
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

S <- read_series("GSE17025_series_matrix.txt.gz")
annot570 <- read.delim("GPL570_annot.tsv", check.names = FALSE)
ex <- collapse(S$ex, annot570)

histv <- sub("histology: ", "", as.character(S$pd[[grep("^histology", colnames(S$pd))[1]]]))
grp <- ifelse(grepl("endometrioid|papillary serous", histv), "tumor",
       ifelse(grepl("endometrium", histv), "normal", NA))
keep <- !is.na(grp)
ex <- normalizeBetweenArrays(ex[, keep, drop = FALSE], method = "quantile")
grp <- factor(grp[keep], levels = c("normal", "tumor"))

## ---- 每样本增殖分数（35 个细胞周期基因的均值 z） ----
CYCLING <- c("MKI67","PCNA","CCNB1","CCNB2","CCNA2","CDK1","CDC20","BUB1","BUB1B",
             "TOP2A","AURKA","AURKB","PLK1","TPX2","UBE2C","RRM2","TYMS","MCM2",
             "MCM3","MCM4","MCM5","MCM6","MCM7","CENPF","KIF11","KIF23","ANLN",
             "PRC1","ECT2","NCAPG","ASPM","DLGAP5","CENPE","SMC4","NUSAP1")
cc <- intersect(CYCLING, rownames(ex))
z <- t(scale(t(ex[cc, , drop = FALSE])))
prolif <- colMeans(z, na.rm = TRUE)
cat("增殖分数：正常组均值", round(mean(prolif[grp == "normal"]), 2),
    "；肿瘤组均值", round(mean(prolif[grp == "tumor"]), 2), "\n")
cat("  组间 t 检验 p =", signif(t.test(prolif ~ grp)$p.value, 3), "\n")

## ---- 校正模型：~ 组 + 增殖分数 ----
design <- model.matrix(~ grp + prolif)
fit <- eBayes(lmFit(ex, design), trend = TRUE, robust = TRUE)
tt <- topTable(fit, coef = "grptumor", number = Inf, sort.by = "P")
write.table(tt, "../DE_GSE17025_tumor_vs_normal_prolifadj.tsv", sep = "\t", quote = FALSE)

hit <- intersect(cc, rownames(tt))
insig <- hit[tt[hit, "adj.P.Val"] > 0.05]
cat("\n校正后增殖基因仍在 FDR<0.05 的：", length(hit) - length(insig), "/", length(hit), "\n")
sig <- rownames(tt)[tt$adj.P.Val < 0.05 & abs(tt$logFC) >= 1]
cat("校正后差基因（FDR<0.05 且 |logFC|>=1）：", length(sig),
    "（UP", sum(tt[sig, "logFC"] > 0), "/ DOWN", sum(tt[sig, "logFC"] < 0), "）\n")
cat("其中仍是增殖基因的：", sum(sig %in% CYCLING), "\n")

## ---- 紧凑版：上下调各 200，按 |logFC| ----
s <- tt[sig, , drop = FALSE]
s$direction <- ifelse(s$logFC > 0, "UP", "DOWN")
s$gene <- rownames(s)
up <- head(s[s$direction == "UP", ][order(-abs(s$logFC[s$direction == "UP"])), ], 200)
dn <- head(s[s$direction == "DOWN", ][order(-abs(s$logFC[s$direction == "DOWN"])), ], 200)
comp <- rbind(up, dn)
write.table(data.frame(gene = comp$gene, direction = comp$direction,
                       logFC = comp$logFC, FDR = comp$adj.P.Val),
            "../EC_signature_v2_prolifadj_compact.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(data.frame(gene = rownames(tt)[rownames(tt) %in% sig],
                       direction = ifelse(tt[sig, "logFC"] > 0, "UP", "DOWN"),
                       logFC = tt[sig, "logFC"], FDR = tt[sig, "adj.P.Val"]),
            "../EC_signature_v2_prolifadj_full.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("\n已写 EC_signature_v2_prolifadj_{compact,full}.tsv\n")
cat("紧凑版中增殖基因数：", sum(comp$gene %in% CYCLING), "\n")
