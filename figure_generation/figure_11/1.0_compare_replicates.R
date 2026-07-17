# chrY should be removed

library(tidyverse)
library(patchwork)
library(UpSetR)

repsdir = "data/My/Integrated/bedRmod_reps/"
outdir = "data/My/Integrated/Illumina_polyA/res/compare_replicates/"

used_methods = c("GLORI", "CAMseq", "eTAMseq", "BSseq", "UMBS")

for (method in used_methods) {
  fbeds = dir(repsdir, pattern = method, full.names = T)
  lbeds = lapply(fbeds, function(x) {
    data.table::fread(x, data.table = FALSE) %>% dplyr::filter(`#chrom` != "chrY")
  })
  names(lbeds) = paste0("rep", 1:length(lbeds))
  sapply(lbeds, nrow)

  lsites = lapply(lbeds, function(x) {
    x %>% mutate(ID = paste(`#chrom`, chromEnd, strand, sep = "_")) %>% pull(ID)
  })

  set_max_n = max(sapply(lsites, length))
  pdf(paste0(outdir, "replicates_Upset_", method, ".pdf"), width = 4, height = 4)
  print(upset(fromList(lsites), order.by = "degree", decreasing = F, 
    set_size.show = TRUE, set_size.scale_max = set_max_n*1.5, sets.x.label = "Number of Sites"))
  dev.off()

  lp = apply(combn(names(lbeds), 2), 2, function(x) {
    print(paste0(x[1], x[2]))

    merged = inner_join(lbeds[[x[1]]], lbeds[[x[2]]], by = c("#chrom", "chromEnd", "strand"))

    ggplot(merged, aes(frequency.x, frequency.y)) + 
      ggrastr::rasterise(ggpointdensity::geom_pointdensity(size = 0.1, method = "neighbors"), dpi = 300) + 
      geom_abline(slope = 1, linetype = "dashed") + 
      scale_color_distiller(palette = "Blues", direction = 1) + 
      scale_x_continuous(limits = c(0, 100)) + 
      scale_y_continuous(limits = c(0, 100)) + 
      xlab(paste0("Modification Level (%) in ", x[1])) +
      ylab(paste0("Modification Level (%) in ", x[2])) + 
      ggpubr::stat_cor() + 
      theme_bw() + 
      theme(panel.grid = element_blank(), aspect.ratio = 1)


  })
  
  num_col = length(lp) %/% 2 + 1
  ggsave(wrap_plots(lp, ncol = num_col), filename = paste0(outdir, "replicates_Correlation_", method, ".pdf"), 
    width = 4*num_col, height = 4*num_col)
}
