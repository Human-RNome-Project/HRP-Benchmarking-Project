rm(list=ls())
library(tidyverse)
library(patchwork)

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


anno = data.table::fread(paste0(outdir, "../Illumina_combined_exon_region.tsv"), data.table = F)
table(anno$region, anno$name)
anno = anno %>% 
  filter(!(X.chrom == "chrY")) %>%
  mutate(ID = paste(X.chrom, chromEnd, strand, sep = "_"))


# ================================================
# 1. top genes (Fig SxC)
# ================================================
anno$frequency = final$frequency[match(anno$ID, final$ID)]

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
data.table::fwrite(df_summary, file = paste0(outdir, ""))


for (mode in c("All", "Coding")) {
  if (mode == "Coding") {
    df = df_summary %>% filter(gene_type == "protein_coding")
  } else {
    df = df_summary
  }
  ## order by mod counts
  used_stat = "count"
  lp1 = lapply(used_mods, function(mod) {
    df_plot = df %>% 
      arrange(desc(.data[[paste0(used_stat, "_", mod)]])) %>%
      slice_head(n = 20) %>%
      mutate(label = factor(gene_name, levels = rev(gene_name)))
    ggplot(df_plot, aes(x = .data[[paste0(used_stat, "_", mod)]], y = label)) + 
      geom_bar(stat = "identity", fill = mod_color[mod]) +
      theme_classic() + 
      scale_x_continuous(position = "top", expand = c(0, NULL)) + 
      xlab(paste("Number of", mod, "sites")) + ylab(NULL) + 
      theme(aspect.ratio = 2) 
  })

  ggsave(wrap_plots(lp1, nrow = 1), 
    filename = paste0(outdir, "Integrated_top20_", mode, "_gene_by_mod_counts.pdf"), width = 10, height = 5)

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
    filename = paste0(outdir, "Integrated_top20_", mode, "_gene_by_mod_levels.pdf"), width = 10, height = 5)
}


used_stat = "freq_sum"
lp2_0 = lapply(used_mods, function(mod) {
  df_top = df_summary %>% 
    filter(gene_type == "protein_coding") %>%
    arrange(desc(.data[[paste0(used_stat, "_", mod)]])) %>%
    slice_head(n = 20)
  
  df_plot = anno %>% filter(region == "Exonic") %>%
    filter(gene_name %in% df_top$gene_name, name == mod) %>%
    mutate(label = factor(gene_name, levels = rev(df_top$gene_name)))
  
  ggplot(df_plot, aes(x = frequency, y = label)) + 
    geom_boxplot(outliers = F, fill = scales::alpha(mod_color[mod], 0.5)) + 
    geom_jitter(color = mod_color[mod], height = 0.1, size = 0.8) +
    theme_classic() + 
    scale_x_continuous(position = "top", limits = c(0, 105)) + 
    xlab(paste(mod, "levels (%)")) + ylab(NULL) + 
    theme(aspect.ratio = 2) 
})
ggsave(wrap_plots(lp2_0, nrow = 1), 
    filename = paste0(outdir, "Integrated_top20_Coding_gene_by_mod_levels_boxplot.pdf"), width = 10, height = 5)


# ================================================
# 2. GO
# ================================================
fun_go = function(gl, main, mod) {
  require(clusterProfiler)
  require(org.Hs.eg.db)
  if (mod %in% names(mod_color)) {
    used_color = mod_color[mod]
  } else {
    used_color = "#3B3B3B"
  }

  firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }

  en = enrichGO(gl, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", 
                ont = "BP", qvalueCutoff = 0.05)
  en = simplify(en)
  dfen = data.frame(en) %>%
    mutate(ID2 = sapply(Description, firstup))

  dfplot = dfen %>% 
    slice_head(n = 10) %>%
    mutate(ID2 = factor(ID2, levels = rev(ID2)))

  p = ggplot(dfplot, aes(x = -log10(p.adjust), y = ID2)) + 
    geom_bar(stat = "identity", fill = used_color) + 
    ylab(NULL) + xlab("-log10(FDR)") + 
    theme_classic() + 
    theme(aspect.ratio = 1.2/20*nrow(dfplot), 
          axis.text = element_text(colour = "black"), 
          axis.ticks = element_line(colour = "black"),
          axis.ticks.length = unit(6, "points")) + 
    ggtitle(paste0(main))
  if (nrow(dfplot) == 0) {
    return(ggplot())
  } else {return(p)}
  # ggsave(p3, file = paste0("05_m6A/Gene_m6A_pileup_", project, "/", main, "_enrich_GO_top20.pdf"), 
  #   width = 6, height = 6)
}


lp3 = lapply(used_mods, function(mod) {
  gl = df_summary %>% 
    arrange(desc(.data[[paste0("count_", mod)]])) %>%
    slice_head(n = 100) %>%
    pull(gene_name)
  fun_go(gl, paste0("GO for top 100 ", mod, " genes"), mod)
})
ggsave(wrap_plots(lp3, ncol = 1), filename = paste0(outdir, "Integrated_GO_top20_genes.pdf"), width = 10, height = 10)


lp4 = lapply(used_mods, function(mod) {
  gl = df_binary %>% 
    filter(.data[[paste0("count_", mod)]] == 1 & rowSums(df_binary[, 2:5]) == 1) %>%
    pull(gene_name) %>% unique()
  fun_go(gl, paste0("GO for ", length(gl), " ", mod, " only genes"), mod)
})
ggsave(wrap_plots(lp4, ncol = 1), filename = paste0(outdir, "Integrated_GO_specific_mod_only_genes.pdf"), width = 10, height = 10)


lp5 = lapply(used_mods, function(mod) {
  gl = df_binary %>% 
    filter(.data[[paste0("count_", mod)]] == 1) %>%
    pull(gene_name) %>% unique()
  fun_go(gl, paste0("GO for ", length(gl), " ", mod, " only genes"), mod)
})
ggsave(wrap_plots(lp5, ncol = 1), filename = paste0(outdir, "Integrated_GO_all_mod_genes.pdf"), width = 10, height = 10)

gene_multimod = df_binary$gene_name[rowSums(df_binary[, 2:5]) > 1]
p4 = fun_go(gene_multimod, "GO for genes with multitype mod", "multi")
ggsave(p4, filename = paste0(outdir, "Integrated_GO_genes_with_multitype_mod.pdf"), width = 10, height = 10)

