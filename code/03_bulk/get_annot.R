options(timeout = 900)
suppressMessages(library(GEOquery))
for (gpl in c("GPL570", "GPL6480")) {
  message("下载 ", gpl, " …")
  g <- tryCatch(getGEO(gpl, destdir = "."),
                error = function(e) { message("失败: ", conditionMessage(e)); NULL })
  if (is.null(g)) next
  tb <- Table(g)
  cat(gpl, "列:", paste(head(colnames(tb), 8), collapse = " | "), "\n")
  cat(gpl, "行数:", nrow(tb), "\n")
  keep <- intersect(c("ID", "GB_ACC", "GENE_SYMBOL", "Gene Symbol", "Symbol",
                      "GENE", "ORF", "SEQUENCE"), colnames(tb))
  cat(gpl, "可用标识列:", paste(keep, collapse = " | "), "\n\n")
  write.table(tb[, unique(c("ID", keep)), drop = FALSE],
              file = paste0(gpl, "_annot.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE)
}
