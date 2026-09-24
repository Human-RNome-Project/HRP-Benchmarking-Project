suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(patchwork)
  library(ggseqlogo)
})

option_list = list(
  make_option(c("-i", "--input"), type = "character",
    help = "Merged Illumina bedRmod, e.g., Illumina_combined_polyARNA_tRNA_rRNA_rmchrY.bed [required]"),
  make_option(c("-f", "--fasta", type = "character"), 
    help = "Genome fasta file used to grep motif, e.g. GRCh38.primary_assembly.genome.fa [required]"),
  make_option(c("-o", "--outdir"), type = "character",
    help = "Output directory [required]")
)

parser = OptionParser(
  usage = "%prog [options]",
  description = "Annotate modifications sites to 5' UTR, CDS and 3' UTR of canonical ensembl transcript of genes",
  option_list = option_list
)
args = parse_args(parser)


used_mods = c("I", "m5C", "m6A", "Y")

# ================================================
# 1. get motifs
# ================================================
final = data.table::fread(args$input, data.table = F, check.names = T) %>%
  mutate(
    ID = paste(X.chrom, chromEnd, strand, sep = "_"),
    RNAtype = case_when(
      grepl("hs_tRNA", X.chrom) ~ "tRNA", 
      grepl("hs_rRNA", X.chrom) ~ "rRNA", 
      .default = "polyA RNA"
    ), 
    name = ifelse(name == "mxU", "U*", name)
  )

polyA = final %>% filter(!grepl("hs_", X.chrom))

polyA %>% 
  mutate(
    chr = `X.chrom`, start = chromStart-4, end = chromEnd+4,
    score = score, name = name, strand = strand
  ) %>% 
  dplyr::select(chr, start, end, name, score, strand) %>%
  filter(start > 0 & (end - start) == 9) %>%
  data.table::fwrite(file = paste0(args$outdir, "/polyA_len9.bed"), sep = "\t", col.names = F, scipen = 999)
system(paste0("bedtools getfasta -s -fi ", args$fasta, " -bed ", 
  args$outdir, "/polyA_len9.bed", " -bedOut > ", args$outdir, "/polyA_len9_motif.bed"))

# ================================================
# 2. read motifs
# ================================================
motif = data.table::fread(paste0(args$outdir, "/polyA_len9_motif.bed"), data.table = F) %>%
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
    cols = base_color), method = "bits", seq_type = "rna") + 
    scale_x_continuous(expand = c(0, 0), breaks = seq(1, 9, 1), 
      label = c(paste0("-", 4:1), "0", paste0("+", 1:4))) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 2), breaks = seq(0, 2, 1)) + 
    theme_classic() + 
    theme(aspect.ratio = 0.3) + 
    ggtitle(paste0("polyA RNA ", mod, " (n=", length(lmotif[[mod]]), ")"))
})

ggsave(wrap_plots(lp), filename = paste0(args$outdir, "/polyA_motif_bits.pdf"), width = 8, height = 4)
