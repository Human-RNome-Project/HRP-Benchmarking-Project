suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(patchwork)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})


option_list = list(
  make_option(c("-i", "--input"), type = "character",
    help = "Merged Illumina bedRmod, e.g., Illumina_combined_polyARNA_tRNA_rRNA_rmchrY.bed [required]"),
  make_option(c("-a", "--anno", type = "character"), 
    help = "Illumina_polyA_mod_annotated.tsv from 0.0_anno_to_gene.R"),
  make_option(c("-o", "--outdir"), type = "character",
    help = "Output directory [required]")
)

parser = OptionParser(
  usage = "%prog [options]",
  description = "Barplot for top20 genes with highest average number of modifications sites",
  option_list = option_list
)
args = parse_args(parser)

mod_color = c(
   # A — m6A family (deep crimson → light blush)
   "m6A" = "#721817",
   "Am" = "#D44F3E",
   "m1A" = "#A52020",
   "mA?" = "#E8907A",
   "m6,6A" = "#F5C4B8",
   "i6A" = "#B12A27", 
   "t6A" = "#EDAA93",

   # C — m5C family (near-black navy → pale sky)
   "m5C" = "#001427",
   "Cm" = "#0D3B6E",
   "mC?" = "#1E6EB5",
   "ac4c" = "#6AAED6",
   "f5C" = "#2B5599", 
   "m3C" = "#BADAF0",

   # G — Inosine family (dark forest → pale mint)
   "I" = "#2D6E1E",
   "Gm" = "#74B354",
   "mG?" = "#4A8532",
   "m2,2,7G" = "#A8D48A",
   "m1G" = "#3B7B31",
   "m1I" = "#BBD89F", 
   "m2,2G" = "#8FC36B",
   "m2G" = "#5D993F",
   "m7G" = "#D5E7C3",

   # U — Psi family (deep amber → pale gold)
   "Y" = "#F0A202",
   "Um" = "#C47A02",
   "mU?" = "#F5BE45",
   "acp3U" = "#B87612", 
   "D" = "#A35B15", 
   "m5U" = "#F8CC59",
   "s2U" = "#FBDCA1",
   "U*" = "#F8D870")

used_mods = c("I", "m5C", "m6A", "Y")

## read sites and annotations
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

anno = data.table::fread(args$anno, data.table = F) %>% 
  filter(!(X.chrom == "chrY")) %>%
  mutate(ID = paste(X.chrom, chromEnd, strand, sep = "_"))


# ================================================
# 1. top genes (Fig SuppC)
# ================================================
anno$frequency = final$frequency[match(anno$ID, final$ID)]
mode = "Coding"

df_summary <- anno %>%
  as.data.frame() %>%
  filter(region == "Exonic") %>%
  group_by(gene_name, name, gene_type) %>%
  summarise(
    count = n(),
    freq_sum = sum(frequency, na.rm = TRUE),
    .groups = "drop" 
  ) %>%
  pivot_wider(
    names_from = name,
    values_from = c(count, freq_sum),
    values_fill = list(count = 0, freq_sum = 0)
  )
# data.table::fwrite(df_summary, sep = "\t", file = paste0(args$outdir, "/Integrated_top20_", mode, "_gene_by_mod_levels.tsv"))



df = df_summary %>% filter(gene_type == "protein_coding")

# order by mod frequency
used_stat = "freq_sum"
lp2 = lapply(used_mods, function(mod) {
  df_plot = df %>% 
    arrange(desc(.data[[paste0(used_stat, "_", mod)]])) %>%
    slice_head(n = 20) %>%
    mutate(label = factor(gene_name, levels = rev(gene_name)))
  ggplot(df_plot, aes(x = .data[[paste0(used_stat, "_", mod)]]/100, y = label)) + 
    geom_bar(stat = "identity", fill = mod_color[mod]) +
    theme_classic() + 
    scale_x_continuous(position = "top", expand = c(0, NULL)) + 
    xlab(paste("Average number of ", mod, " sites\nper transcripts")) + ylab(NULL) + 
    theme(aspect.ratio = 2) 
})

ggsave(wrap_plots(lp2, nrow = 1), 
  filename = paste0(args$outdir, "/Integrated_top20_", mode, "_gene_by_mod_levels.pdf"), width = 10, height = 5)