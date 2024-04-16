setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
library(HilbertCurve)
library(circlize)
set.seed(12345)

vl_hilbert <- function(track,
                       title= "pixel mode",
                       clip= quantile(cov, .9995),
                       width= 2500)
{
  # Bin genome ----
  bins <- vl_binBSgenome("dm6",
                         bins.width = width,
                         restrict.seqnames = c("chrX", "chr2L", "chr2R", "chr3L", "chr3R", "chr4"))
  bins[, idx:= .I]
  # Compute coverage and clip extremes ----
  cov <- vl_bw_coverage(bins, track)
  cov[cov>clip] <- clip
  print(paste("Max value set to", clip))
  # Color function ----
  col_fun = colorRamp2(c(0, max(cov)), c("white", "red"))
  # Initiate curve ----
  hc = HilbertCurve(1,
                    max(bins$idx)+1,
                    level = 10,
                    mode = "pixel",
                    title = title)
  # Add signal layer ----
  x1 <- bins$idx
  x2 <- bins$idx+1
  hc_layer(hc,
           x1 = x1,
           x2 = x2,
           mean_mode = "absolute",
           col = col_fun(cov))
  # Add chromosome polygons ----
  x1 <- bins[, idx[1], seqnames]$V1
  x2 <- bins[, idx[.N], seqnames]$V1+1
  hc_polygon(hc, x1 = x1, x2 = x2)
  hc_text(hc,
          x1 = x1,
          x2 = x2, 
          labels = unique(bins$seqnames),
          gp = gpar(fontsize = 10),
          centered_by = "polygon")
}

pdf("pdf/chetan_hilbert_curvers.pdf", 3, 3)
vl_hilbert("db/bw/cutnrun/H2AK118Ub_PH18_merge.bw", title= "H2AK118Ub no ph-KD", clip = 8)
vl_hilbert("db/bw/cutnrun/H2AK118Ub_PHD11_merge.bw", title= "H2AK118Ub transient ph-KD", clip = 8)
vl_hilbert("db/bw/cutnrun/H2AK118Ub_PH29_merge.bw", title= "H2AK118Ub constant ph-KD", clip = 8)
vl_hilbert("db/bw/cutnrun/H3K27me3_PH18_merge.bw", title= "H3K27me3 no ph-KD", clip= 9)
vl_hilbert("db/bw/cutnrun/H3K27me3_PHD11_merge.bw", title= "H3K27me3 transient ph-KD", clip= 9)
vl_hilbert("db/bw/cutnrun/H3K27me3_PH29_merge.bw", title= "H3K27me3 constant ph-KD", clip= 9)
vl_hilbert("db/bw/cutnrun/PH_PH18_merge.bw", title= "PH no ph-KD", clip= 25)
vl_hilbert("db/bw/cutnrun/PH_PHD11_merge.bw", title= "PH transient ph-KD", clip= 25)
vl_hilbert("db/bw/cutnrun/PH_PH29_merge.bw", title= "PH constant ph-KD", clip= 25)
dev.off()