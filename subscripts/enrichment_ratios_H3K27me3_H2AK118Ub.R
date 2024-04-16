setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)

# Import data ----
dat <- readRDS("Rdata/final_gene_features_table.rds")
dat <- dat[!is.na(log2FoldChange_PH18)]

# Define controls ----
PcG_bound <- dat[(PcG_bound)]
ctl <- vl_select_act_matched(group= dat[(PcG_bound)],
                             control= dat[!(PcG_bound)],
                             act.col.name = "baseMean")

# Define regions ----
TSS_bound <- vl_resizeBed(PcG_bound, "start", 0, 0)
TSS_unbound <- vl_resizeBed(ctl, "start", 0, 0)
window_bound <- vl_resizeBed(PcG_bound, "start", 5000, 5000)
window_unbound <- vl_resizeBed(ctl, "start", 5000, 5000)

# Plot function ----
pl <- function(tracks,
               names,
               main= "test",
               ylim1= NULL,
               ylim2= NULL)
{
  # Quantif signal
  upstream <- 25000
  downstream <- 75000
  PcG_bound <- vl_bw_average_track(TSS_bound,
                                   names= names,
                                   tracks = tracks,
                                   upstream = upstream,
                                   downstream = downstream,
                                   plot= F)
  unbound <- vl_bw_average_track(TSS_unbound,
                                 names= names,
                                 tracks = tracks,
                                 upstream = upstream,
                                 downstream = downstream,
                                 plot= F)
  # Compute ratio ----
  PcG_bound <- PcG_bound[, .(score= mean(score)), .(name, bin.x)]
  unbound <- unbound[, .(score= mean(score)), .(name, bin.x)]
  .c <- merge(PcG_bound, unbound, by= c("bin.x", "name"))
  
  # Plot ----
  vl_par(mai= c(.9, .9, .9, .4),
         mgp= c(1.25,  0.35, 0))
  plot(NA,
       xlim= range(.c$bin.x),
       ylim= range(.c[, score.x/score.y]),
       frame= F,
       ylab= "Enrich. ratio PcG bound/unbound",
       xaxt= "n",
       xlab= NA)
  axis(1,
       c(-upstream, 0, downstream),
       c(-upstream, "TSS", downstream))
  Cc <- c("lightgrey", "limegreen", "tomato")
  .c[, {
    lines(bin.x,
          score.x/score.y, col= Cc[.GRP])
  }, name]
  title(main= main)
  abline(v= c(-2500, 1500), lty= "11")
  legend("topright",
         fill= Cc,
         legend= unique(.c$name),
         bty= "n",
         cex= .7)
  
  # Quantification ----
  PcG_bound <- lapply(tracks, function(x) data.table(score= vl_bw_coverage(window_bound, x)))
  names(PcG_bound) <- names
  unbound <- lapply(tracks, function(x) data.table(ctl= vl_bw_coverage(window_unbound, x)))
  names(unbound) <- names
  m <- cbind(rbindlist(PcG_bound, idcol = T),
             rbindlist(unbound))
  m[, .id:= factor(.id, names)]
  vl_par(mai= c(.9, .4, .9, .2),
         mgp= c(1,  0.35, 0))
  vl_boxplot(score/ctl~.id,
             m,
             col= Cc,
             compute.pval= list(c(1,2), c(1,3), c(2,3)),
             names= levels(names),
             ylab= "Enrich. ratio PcG bound/unbound",
             tilt.names= T)
}

# Plot ----
pdf("pdf/chetan_average_tracks.pdf", 5, 3)
layout(matrix(1:2, ncol= 2),
       widths = c(1, .3))
pl(tracks= c("db/bw/cutnrun/H2AK118Ub_PH18_merge.bw",
             "db/bw/cutnrun/H2AK118Ub_PHD11_merge.bw",
             "db/bw/cutnrun/H2AK118Ub_PH29_merge.bw"),
   names= factor(c("PH18", "PHD11", "PH29"), c("PH18", "PHD11", "PH29")),
   main= "H2AK118Ub")
pl(tracks= c("db/bw/cutnrun/H3K27me3_PH18_merge.bw",
             "db/bw/cutnrun/H3K27me3_PHD11_merge.bw",
             "db/bw/cutnrun/H3K27me3_PH29_merge.bw"),
   names= factor(c("PH18", "PHD11", "PH29"), c("PH18", "PHD11", "PH29")),
   main= "H3K27me3")
# pl("db/bw/cutnrun/H2AK118Ub_PH18_merge.bw", "H2AK118Ub no ph-KD", ylim1= c(1.5, 4.5), ylim2 = c(0, 8))
# pl("db/bw/cutnrun/H2AK118Ub_PHD11_merge.bw", "H2AK118Ub transient ph-KD", ylim1= c(1.5, 4.5), ylim2 = c(0, 8))
# pl("db/bw/cutnrun/H2AK118Ub_PH29_merge.bw", "H2AK118Ub constant ph-KD", ylim1= c(1.5, 4.5), ylim2 = c(0, 8))
# pl("db/bw/cutnrun/H3K27me3_PH18_merge.bw", "H3K27me3 no ph-KD")
# pl("db/bw/cutnrun/H3K27me3_PHD11_merge.bw", "H3K27me3 transient ph-KD")
# pl("db/bw/cutnrun/H3K27me3_PH29_merge.bw", "H3K27me3 constant ph-KD")
dev.off()