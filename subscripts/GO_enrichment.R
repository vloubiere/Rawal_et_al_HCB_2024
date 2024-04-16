setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")

# Import transcriptome
dat <- readRDS("Rdata/final_gene_features_table.rds")
genes <- list("Up constant ph-KD"= dat[diff_PH29=="up", FBgn],
              "Down constant ph-KD"= dat[diff_PH29=="down", FBgn],
              "Up transient ph-KD"= dat[diff_PHD11=="up", FBgn],
              "Down transient ph-KD"= dat[diff_PHD11=="down", FBgn])

enr <- vl_GO_enrich(genes, dat$FBgn, "Dm", "BP")

grep("replication|repair|break|damage|instability", unique(enr[padj<0.05 & log2OR>0, name]), value = T)
sel <- c("behavior",
         "cell−cell signaling",
         "cellular response to DNA damage stimulus",
         "compound eye photoreceptor cell differentiation",
         "DNA metabolic process",
         "DNA repair",
         "DNA replication",
         "DNA−templated DNA replication",
         "G protein−coupled receptor signaling pathway",
         "imaginal disc development",
         "ion transport",
         "muscle cell differentiation",
         "neuromuscular synaptic transmission",
         "regulation of secretion by cell",
         "regulation of trans−synaptic signaling",
         "regulation of transcription by RNA polymerase II",
         "regulation of transport",
         "secretion",
         "sex differentiation",
         "striated muscle cell development",
         "striated muscle cell differentiation",
         "tissue development",
         "transcription by RNA polymerase II",
         "vesicle−mediated transport in synapse",
         "muscle system process",
         "regulation of ion transport")

pdf("pdf/chetan_GO_enrichments.pdf", 7, 8)
vl_par(mai= c(1, 4.5, 1, 2))
plot(enr[name %in% c(enr[set_hit>10 & log2OR>1, name])],
     padj.cutoff= 0.01,
     cex.balloons= .35)
vl_par(mai= c(2.75, 4.6, 2.5, 2))
plot(enr[name %in% sel],
     padj.cutoff= 0.01,
     cex.balloons= .5)
dev.off()