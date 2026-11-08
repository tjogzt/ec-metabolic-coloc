suppressMessages({library(GEOquery); library(limma)})
# w4_de.R — 慢性 IR 状态两个数据集的差异分析
#
# 与预注册的偏离（在看结果之前记录）：
#   1) GSE15773 实际为配对设计：每例受试者取网膜 + 皮下两个脂肪库（IR 9 例、IS 10 例）。
#      预注册只写"单组比较"，未预见库间配对；为避免脂肪库差异混入组效应，
#      模型改为 ~ group + depot（库作协变量）。
#   2) GSE18732 含 NGT / IGT / DM 三组，预注册只规定 DM vs NGT，
#      故剔除 IGT，不做三组回归。
setwd("data")
options(timeout = 900)

read_series <- function(f) {
  es <- getGEO(filename = f, getGPL = FALSE)
  ex <- exprs(es); pd <- pData(es)
  colnames(ex) <- pd$title
  list(ex = ex, pd = pd)
}

collapse <- function(ex, annot, symcol = NULL) {
  if (is.null(symcol) || !(symcol %in% colnames(annot)))
    symcol <- c("Gene Symbol", "Gene symbol", "Symbol", "GENE_SYMBOL", "GENE")[
      c("Gene Symbol", "Gene symbol", "Symbol", "GENE_SYMBOL", "GENE") %in% colnames(annot)][1]
  a <- annot[match(rownames(ex), annot$ID), c("ID", symcol)]
  ok <- !is.na(a[[symcol]]) & a[[symcol]] != "" & a[[symcol]] != "---"
  ex <- ex[ok, , drop = FALSE]; a <- a[ok, , drop = FALSE]
  # 多符号条目（"A /// B"）必须拆开：否则这些符号以伪基因名留在索引里，
  # 而只以多符号形式注释的基因会整体缺席（GPL570 实测丢失 2787 个符号）。
  parts <- strsplit(trimws(as.character(a[[symcol]])), "\\s*///\\s*")
  nrep <- lengths(parts)
  if (any(nrep > 1)) {
    ex <- ex[rep(seq_len(nrow(ex)), nrep), , drop = FALSE]
    a  <- a[rep(seq_len(nrow(a)),  nrep), , drop = FALSE]
    a[[symcol]] <- unlist(parts, use.names = FALSE)
    cat(sprintf("  [折叠] 多符号条目 %d 行 → 展开为 %d 行\n",
                sum(nrep > 1), sum(nrep)))
  }
  # 若为原始荧光强度（最大值远大于 50），先做 log2(x+1)，否则 limma 会给出虚高 logFC
  if (max(ex, na.rm = TRUE) > 50) {
    cat(sprintf("  [预处理] 原始强度（max %.0f）→ log2(x+1)\n", max(ex, na.rm = TRUE)))
    ex <- log2(ex + 1)
  }
  m <- rowMeans(ex, na.rm = TRUE)
  o <- order(a[[symcol]], -m); a <- a[o, ]; ex <- ex[o, , drop = FALSE]
  first <- !duplicated(a[[symcol]])
  ex <- ex[first, , drop = FALSE]; rownames(ex) <- a[[symcol]][first]
  ex
}

get_annot <- function(gpl, dest) {
  if (!file.exists(dest)) {
    message("取平台注释 ", gpl)
    g <- getGEO(gpl, destdir = ".")
    tb <- Table(g)
    write.table(tb, dest, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  read.delim(dest, check.names = FALSE)
}

report <- function(tt, label) {
  cat(sprintf("  FDR<0.05 且 |logFC|>=0.3：%d\n",
              sum(tt$adj.P.Val < 0.05 & abs(tt$logFC) >= 0.3)))
  cat(sprintf("  上调前 5：%s\n", paste(head(rownames(tt)[tt$logFC > 0], 5), collapse = ", ")))
  cat(sprintf("  下调前 5：%s\n", paste(head(rownames(tt)[tt$logFC < 0], 5), collapse = ", ")))
  write.table(tt, paste0("../DE_", label, "_state.tsv"), sep = "\t", quote = FALSE)
}

## ================= GSE15773 脂肪：IR vs IS（含 depot 协变量） =================
S <- read_series("GSE15773_series_matrix.txt.gz")
ex <- collapse(S$ex, get_annot("GPL570", "GPL570_annot.tsv"))
ttl <- colnames(ex)
grp <- ifelse(grepl("Insulin resistant", ttl, ignore.case = TRUE), "IR",
       ifelse(grepl("Insulin sensitive", ttl, ignore.case = TRUE), "IS", NA))
depot <- ifelse(grepl("omental", ttl, ignore.case = TRUE), "omental",
         ifelse(grepl("subcutaneous", ttl, ignore.case = TRUE), "subcut", NA))
keep <- !is.na(grp) & !is.na(depot)
ex <- normalizeBetweenArrays(ex[, keep, drop = FALSE], method = "quantile")
g <- factor(grp[keep], levels = c("IS", "IR"))
dp <- factor(depot[keep])
cat(sprintf("\n=== GSE15773 脂肪 ===\n基因 %d，样本 %d（IR %d / IS %d；网膜 %d / 皮下 %d）\n",
            nrow(ex), sum(keep), sum(g == "IR"), sum(g == "IS"),
            sum(dp == "omental"), sum(dp == "subcut")))
print(table(g, dp))
# 19 张芯片来自 10 名受试者（每人皮下+网膜各一张），必须按受试者作区组。
# 只把 depot 当协变量会把同一人的两张芯片当独立样本（伪重复，残差 df 16 vs 应为 8）。
subj <- paste(g, sub(".*patient\\s*([0-9]+).*", "\\1", ttl[keep], ignore.case = TRUE), sep = "_")
cat(sprintf("  受试者 %d 名，每人 %s 张芯片\n", length(unique(subj)),
            paste(range(table(subj)), collapse = "-")))
design <- model.matrix(~ g + dp)
dup <- duplicateCorrelation(ex, design, block = subj)
cons <- dup$consensus
if (!is.finite(cons) || cons < 0) { cons <- 0; cat("  [警告] 区组内相关非正，按 0 处理\n") }
cat(sprintf("  区组内相关 consensus = %.3f\n", cons))
fit <- eBayes(lmFit(ex, design, block = subj, correlation = cons),
              trend = TRUE, robust = TRUE)
# 注意：duplicateCorrelation 不改变设计矩阵的残差自由度（本例 19-3=16），
# 区组内相关通过方差-covariance 进入，而不是靠消耗 df。下面按"样本数"印。
cat(sprintf("  模型 df.residual = %d（设计参数 %d，样本 %d）\n",
            fit$df.residual[1], ncol(design), ncol(ex)))
tt_a <- topTable(fit, coef = "gIR", number = Inf, sort.by = "P")
report(tt_a, "GSE15773_adipose_IR")

## ================= GSE18732 肌肉：DM vs NGT（剔除 IGT） =================
# 平台 GPL9486 = Affymetrix U133 Plus 2.0 的**转录本级自定义 CDF**（Hs133P_Hs_ENST v10），
# 探针 ID 是 Ensembl 转录本号（ENST…），不是 U133P2 探针号，故不能用 GPL570 注释。
# 映射改用本地 GENCODE v44 注释（ref2/gencode.v44.annotation.gtf.gz）建立的 ENST→符号表，
# 无需下载平台文件。
S2 <- read_series("GSE18732_series_matrix.txt.gz")
E2 <- read.delim("ENST2SYMBOL.tsv", check.names = FALSE)
.miss <- !rownames(S2$ex) %in% E2$ID
cat(sprintf("  ENST 未匹配 %d/%d（%.1f%%）——GENCODE v44 不含旧稳定号所致\n",
            sum(.miss), length(.miss), 100 * mean(.miss)))
if (mean(.miss) > 0.10)
  warning(sprintf("ENST 映射丢失 %.1f%%，该队列基因宇宙与其他队列不可比", 100 * mean(.miss)))
ex2 <- collapse(S2$ex, E2)
ex2 <- ex2[complete.cases(ex2) & is.finite(rowSums(ex2)), , drop = FALSE]
cat(sprintf("  ENST→符号折叠后 %d 基因（已去含 NA 的行）\n", nrow(ex2)))
ttl2 <- colnames(ex2)
g2 <- ifelse(grepl("diabetic", ttl2, ignore.case = TRUE), "DM",
      ifelse(grepl("normal", ttl2, ignore.case = TRUE), "NGT", NA))
keep2 <- !is.na(g2)
ex2 <- normalizeBetweenArrays(ex2[, keep2, drop = FALSE], method = "quantile")
g2 <- factor(g2[keep2], levels = c("NGT", "DM"))
cat(sprintf("\n=== GSE18732 肌肉 ===\n基因 %d，样本 %d（DM %d / NGT %d；已剔除 IGT）\n",
            nrow(ex2), sum(keep2), sum(g2 == "DM"), sum(g2 == "NGT")))
fit2 <- eBayes(lmFit(ex2, model.matrix(~ g2)), trend = TRUE, robust = TRUE)
tt_b <- topTable(fit2, coef = "g2DM", number = Inf, sort.by = "P")
report(tt_b, "GSE18732_muscle_DM")

## ================= 跨组织一致性（预注册第 3 条） =================
sa <- tt_a[tt_a$adj.P.Val < 0.05 & abs(tt_a$logFC) >= 0.3, ]
sb <- tt_b[tt_b$adj.P.Val < 0.05 & abs(tt_b$logFC) >= 0.3, ]
common <- intersect(rownames(sa), rownames(sb))
cat(sprintf("\n脂肪显著 %d，肌肉显著 %d，交集 %d\n", nrow(sa), nrow(sb), length(common)))
if (length(common) > 5) {
  ag <- sign(sa[common, "logFC"]) == sign(sb[common, "logFC"])
  cat(sprintf("方向一致率 %.1f%%（%d/%d），二项 p = %.3g\n",
              100 * mean(ag), sum(ag), length(common),
              binom.test(sum(ag), length(common), 0.5)$p.value))
}
cat("\n完成\n")

## ================= 修正案 A2：GSE25401 脂肪 肥胖 vs 非肥胖 =================
# 平台 GPL6244（Affy Human Gene 1.0 ST），注释取 GEO 的 .annot 文件
read_annot_geo <- function(f) {
  ln <- readLines(gzfile(f), warn = FALSE)
  i <- grep("^ID\t", ln)[1]
  if (is.na(i)) stop("annot 文件里找不到以 ID<TAB> 开头的表头")
  read.delim(text = ln[i:length(ln)], check.names = FALSE, quote = "")
}
S3 <- read_series("GSE25401_series_matrix.txt.gz")
ex3 <- collapse(S3$ex, read_annot_geo("GPL6244.annot.gz"))
ttl3 <- colnames(ex3)
g3 <- ifelse(grepl("Obese", ttl3, ignore.case = TRUE) & !grepl("Nonobese", ttl3, ignore.case = TRUE), "Obese",
      ifelse(grepl("Nonobese", ttl3, ignore.case = TRUE), "Nonobese", NA))
keep3 <- !is.na(g3)
ex3 <- normalizeBetweenArrays(ex3[, keep3, drop = FALSE], method = "quantile")
ex3 <- ex3[complete.cases(ex3) & is.finite(rowSums(ex3)), , drop = FALSE]
g3 <- factor(g3[keep3], levels = c("Nonobese", "Obese"))
cat(sprintf("\n=== GSE25401 脂肪（A2 发现 A） ===\n基因 %d，样本 %d（Obese %d / Nonobese %d）\n",
            nrow(ex3), sum(keep3), sum(g3 == "Obese"), sum(g3 == "Nonobese")))
fit3 <- eBayes(lmFit(ex3, model.matrix(~ g3)), trend = TRUE, robust = TRUE)
tt_c <- topTable(fit3, coef = "g3Obese", number = Inf, sort.by = "P")
report(tt_c, "GSE25401_adipose_obese")
