rm(list=ls())
library(tidyverse)
library(patchwork)
library(ggseqlogo)

used_mods = c("I", "m5C", "m6A", "Y")
datadir = "data/My/Integrated/bedRmod_final/"
outdir = "data/My/Integrated/Illumina_polyA/res/Integrated/"

# ================================================
# 1. number of sites by mod
# ================================================
file = paste0(datadir, "Illumina_combined_polyARNA_tRNA_rRNA_rmchrY.bed")
final = data.table::fread(file, data.table = F, check.names = T) %>%
  mutate(
    ID = paste(X.chrom, chromEnd, strand, sep = "_"),
    RNAtype = case_when(
      grepl("hs_tRNA", X.chrom) ~ "tRNA", 
      grepl("hs_rRNA", X.chrom) ~ "rRNA", 
      .default = "polyA RNA"
    ), 
    name = ifelse(name == "mxU", "U*", name)
  )

table(final$name, final$RNAtype)

polyA = final %>% filter(!grepl("hs_", X.chrom))
table(polyA$name)

## read motif (got use getfasta)
# bedtools getfasta -s -fi ${fasta} -bed ${len9_bed} -bedOut > ${motif_bed}
motif = data.table::fread(paste0(outdir, "../motif/polyA_motif.bed"), data.table = F) %>%
    mutate(ID = paste(V1, V3-4, V6, sep = "_"))
table(polyA$ID %in% motif$ID, polyA$name)
table(duplicated(motif$ID))

polyA$motif = gsub("T", "U", toupper(motif$V7[match(polyA$ID, motif$ID)]))

lmotif = split(polyA$motif, polyA$name) %>%
  lapply(function(x) x[!is.na(x)])
sapply(lmotif, length)

base_color = c("A" = "#721817", "C" = "#001427", "G" = "#2D6E1E", "U" = "#F0A202")
lp = lapply(c("I", "m5C", "m6A", "Y"), function(mod) {
  ggseqlogo(lmotif[mod], ncol = , col_scheme = make_col_scheme(chars = names(base_color), 
    cols = base_color), method = "prob", seq_type = "rna") + 
    scale_x_continuous(expand = c(0, 0), breaks = seq(1, 9, 1), 
      label = c(paste0("-", 4:1), "0", paste0("+", 1:4))) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1), breaks = seq(0, 1, 0.5)) + 
    theme_classic() + 
    theme(aspect.ratio = 0.3) + 
    ggtitle(paste0("polyA RNA ", mod, " (n=", length(lmotif[[mod]]), ")"))
})

ggsave(wrap_plots(lp), filename = paste0(outdir, "/polyA_motif.pdf"), width = 8, height = 4)


lp = lapply(c("I", "m5C", "m6A", "Y"), function(mod) {
  ggseqlogo(lmotif[mod], ncol = , col_scheme = make_col_scheme(chars = names(base_color), 
    cols = base_color), method = "bits", seq_type = "rna") + 
    scale_x_continuous(expand = c(0, 0), breaks = seq(1, 9, 1), 
      label = c(paste0("-", 4:1), "0", paste0("+", 1:4))) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 2), breaks = seq(0, 2, 1)) + 
    theme_classic() + 
    theme(aspect.ratio = 0.3) + 
    ggtitle(paste0("polyA RNA ", mod, " (n=", length(lmotif[[mod]]), ")"))
})

ggsave(wrap_plots(lp), filename = paste0(outdir, "/polyA_motif_bits.pdf"), width = 8, height = 4)
