require(fgsea)

# Retrieve GO terms of interest
tab <- data.table(file= list.files("db/DNA_repair/", "IDs.txt$", full.names = T))
tab[, GO:= gsub("_IDs.txt$", "", basename(file))]
tab <- tab[, fread(file, col.names = "FBgn", header= F), GO]
pathways <- split(tab$FBgn, tab$GO)

# Import transcriptome ----
dat <- readRDS("Rdata/final_gene_features_table.rds")
dat <- dat[!is.na(log2FoldChange_PH29) & !is.na(log2FoldChange_PHD11)]
# clean to avoid similar FC (=0)
set.seed(1)
dat[, clean_PH29:= if(.N==1) log2FoldChange_PH29 else jitter(rep(log2FoldChange_PH29, .N), amount = .01), log2FoldChange_PH29]
dat[, clean_PHD11:= if(.N==1) log2FoldChange_PHD11 else jitter(rep(log2FoldChange_PHD11, .N), amount = .01), log2FoldChange_PHD11]

# Save table
GOs <- merge(dat, tab, by= "FBgn")
GOs <- GOs[, .(GO, FBgn, symbol, PcG_bound,
               log2FoldChange_PH18, log2FoldChange_PHD11, log2FoldChange_PH29, 
               padj_PH18, padj_PHD11, padj_PH29,
               diff_PH18, diff_PHD11, diff_PH29)]
setorderv(GOs, c("GO", "log2FoldChange_PH29"), c(1, -1))
fwrite(GOs,
       "db/DNA_repair/table_genes.txt",
       col.names = T,
       row.names = F,
       sep= "\t",
       quote= F)