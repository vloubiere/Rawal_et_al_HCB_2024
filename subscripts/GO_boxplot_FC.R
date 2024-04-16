require(vlfunctions)

full <- readRDS("Rdata/final_gene_features_table.rds")
dat <- fread("db/DNA_repair/table_genes.txt")

# Cc <- c("lightgrey", "cornflowerblue", "tomato")
Cc <- c("lightgrey", "tomato")

pdf("pdf/chetan_boxlplot_FC.pdf", 8, 3)
vl_par(mfrow= c(1,3),
       font.main= 1,
       lwd= .25)
dat[, {
  vl_boxplot(.(full$log2FoldChange_PH18,
               log2FoldChange_PH18,
               full$log2FoldChange_PHD11,
               log2FoldChange_PHD11,
               full$log2FoldChange_PH29,
               log2FoldChange_PH29),
             dat,
             main= paste0(GO, " (n=", .N, ")"),
             violin= T,
             tilt.names= T,
             col= "white",
             viocol= adjustcolor(Cc, .3),
             at= rep(c(1,2.5,4), each= 2)+c(-0.25, 0.25),
             compute.pval= list(c(1,2), c(3,4), c(5, 6), c(2,4), c(2,6)),
             boxwex= .075,
             viowex= .4,
             lwd= .25,
             xaxt= "n",
             ylab= "RNA-Seq fold change (log2)")
  vl_tilt_xaxis(c(1,2.5,4),
                labels= c("No ph-KD", "Transient ph-KD", "Constant ph-KD"))
  abline(h= 0,
         lwd= .25,
         lty= "11")
  vl_legend(fill= adjustcolor(Cc, .3),
            legend = c("All genes", "GO term"),
            border= "black")
  .SD
}, GO]
dev.off()