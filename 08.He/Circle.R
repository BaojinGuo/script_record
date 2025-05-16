
library(circlize)

fai_file <- "NLL.plot.fai"
atac_dir <- "./"

fai <- read.table(fai_file, header = FALSE, stringsAsFactors = FALSE)
colnames(fai) <- c("chr", "length")
chr_lengths <- fai[, c("chr", "length")]
atac_files <- list.files(atac_dir, pattern = "\\.bed$", full.names = TRUE)

atac_list <- lapply(atac_files, function(file) {
  df <- read.table(file, header = FALSE)
  colnames(df) <- c("chr", "start", "end", "id", "score")
  return(df)
})
names(atac_list) <- basename(atac_files)

circos.clear()
circos.par(start.degree = 90, gap.after = c(rep(1, nrow(chr_lengths)-1), 10))

circos.initialize(factors = chr_lengths$chr,
                  xlim = cbind(rep(0, nrow(chr_lengths)), chr_lengths$length))

# 染色体名称 track
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  chr <- CELL_META$sector.index
  xlim <- CELL_META$xlim
  circos.text(mean(xlim), 0.5, chr, cex = 0.6, facing = "inside", niceFacing = TRUE)
}, bg.border = NA, track.height = 0.05)


colors <- rainbow(length(atac_list))

for (i in seq_along(atac_list)) {
  dat <- atac_list[[i]]
  col <- colors[i]
  circos.trackPlotRegion(
    ylim = c(0, max(dat$score, na.rm = TRUE)),
    panel.fun = function(region, value, ...) {
      chr <- CELL_META$sector.index
      dat_chr <- dat[dat$chr == chr, ]
      if (nrow(dat_chr) > 0) {
        circos.rect(dat_chr$start, 0, dat_chr$end, dat_chr$score,
                    col = col, border = NA)
      }
    },
    track.height = 0.1,
    bg.border = NA
  )
}

legend("right", legend = names(atac_list), fill = colors, cex = 0.6, bty = "n")
