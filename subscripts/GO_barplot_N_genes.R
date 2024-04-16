require(vlfunctions)

dat <- fread("db/DNA_repair/table_genes.txt")
dat <- melt(dat, c("GO", "FBgn"), patterns("^diff"))

pdf("pdf/chetan_barplot_up_down_genes.pdf", 10, 3)
layout(matrix(1:3, ncol= 3))
vl_par(font.main= 1,
       omi= c(0, 0, 0, 2))
Cc <- c("tomato", "cornflowerblue")
dat[, {
  .c <- .SD[, .(up= sum(value=="up"), down= sum(value=="down")), variable]
  bar <- barplot(.c$up,
                 col= "tomato",
                 main= paste0(GO, " (", length(unique(FBgn)), ")"),
                 ylab= "Number of genes",
                 ylim= c(-max(.c$down), max(.c$up)))
  barplot(-.c$down,
          col= "cornflowerblue",
          main= GO,
          ylab= "Number of genes",
          add= T)
  vl_tilt_xaxis(bar,
                par("usr")[3]-strheight("M")*2,
                labels = c("No ph-KD", "Transient ph-KD", "Constant ph-KD"))
  text(bar, .c$up, .c$up, xpd= T, pos= 3, offset= 0.25, cex= 5/6)
  text(bar, -.c$down, .c$down, xpd= T, pos= 1, offset= 0.25, cex= 5/6)
  .SD
}, GO]
legend(par("usr")[2],
       par("usr")[4],
       fill= Cc,
       legend = c("Up-regulated",
                  "Down regulated"),
       xpd= NA,
       bty= "n")
dev.off()