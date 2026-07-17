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
# 1. number of sites by mod (Fig6A, SxA, SxB)
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

Nsites_by_mod = table(final$name, final$RNAtype) %>% 
  as.data.frame.matrix() %>%
  mutate(Modification = rownames(.), .before = 1)
data.table::fwrite(Nsites_by_mod, file = paste0(outdir, "Number_of_sites_by_mod.tsv"), sep = "\t")


polyA = final %>% filter(!grepl("hs_", X.chrom))

df_nsites = final %>% 
  group_by(name, RNAtype) %>%
  summarise(n = n())

lp = lapply(c("polyA RNA", "rRNA", "tRNA"), function(x) {
  df = df_nsites %>% filter(RNAtype == x)
  maxn = max(df$n)
  ggplot(df, aes(name, n, fill = name)) + 
    geom_bar(stat = "identity", show.legend = F) + 
    scale_fill_manual(values = mod_color) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
    geom_text(aes(y = n + maxn*0.05, label = n)) + 
    xlab("Modification Type") + ylab("Number of sites") + 
    theme_classic() + 
    theme(aspect.ratio = 8/nrow(df), axis.text.x = element_text(angle = 45, hjust = 1)) + 
    ggtitle(x) 
})
ggsave(wrap_plots(lp, nrow = 1), filename = paste0(outdir, "Integrated_barplot_num_of_sites.pdf"), width = 18, height = 5)


# ================================================
# 2. metagene plot
# ================================================
anno = data.table::fread(paste0(outdir, "../Illumina_combined_exon_region.tsv"), data.table = F)
table(anno$region, anno$name)
anno = anno %>% 
  filter(!(X.chrom == "chrY")) %>%
  mutate(ID = paste(X.chrom, chromEnd, strand, sep = "_"))

## Fig6C
p2_0 = ggplot(anno %>% mutate(name = factor(name, levels = c("I", "m5C", "m6A", "Y"))), aes(x = name, fill = region)) + 
  geom_bar(stat = "count", position = "fill") + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0), label = scales::percent) + 
  scale_fill_manual(values = ggsci::pal_nejm()(3)) + 
  xlab(NULL) + ylab("Proportion of sites (%)") + 
  theme_bw() + 
  theme(aspect.ratio = 2, legend.position = "top")
ggsave(p2_0, filename = paste0(outdir, "Integrated_region_proportion.pdf"), width = 4, height = 5)

## Fig6D
p2 = ggplot(anno %>% filter(region == "Exonic"), aes(rel_location, fill = name, color = name)) + 
  geom_density() + 
  geom_vline(xintercept = c(1,2), color = "gray") + 
  scale_x_continuous(breaks = seq(0.5, 2.5, 1), labels = c("5'UTR", "CDS", "3'UTR")) +
  scale_color_manual(values = mod_color) + 
  labs(color = "Modification", fill = "Modification") + 
  scale_fill_manual(values = scales::alpha(mod_color, 0.2)) + 
  xlab(NULL) + ylab("Density") + 
  theme_classic() + 
  theme(aspect.ratio = 1) + 
  ggtitle("polyA RNA")
ggsave(p2, filename = paste0(outdir, "Integrated_metagene.pdf"), width = 4, height = 5)


## peak around stop codon
anno = anno %>% mutate(dis2stop = end - utr5_len - cds_len)
ggplot(anno %>% filter(region == "Exonic"), aes(dis2stop, fill = name, color = name)) + 
  geom_bar(stat = "count") + 
  scale_color_manual(values = mod_color) + 
  coord_cartesian(xlim = c(-250, 250)) + 
  labs(color = "Modification", fill = "Modification") + 
  scale_fill_manual(values = scales::alpha(mod_color, 0.2)) + 
  xlab("Distance to stop codon") + ylab("Number of sites") + 
  theme_classic() + 
  ggtitle("polyA RNA") + 
  facet_wrap(~name, ncol = 2, scales = "free_y")



## plot mean frequency (Fig6F)
anno = anno %>% mutate(frequency = final$frequency[match(ID, final$ID)])
ggplot(anno, aes(rel_location, frequency, color = name)) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 20), linewidth = 1) + 
  geom_vline(xintercept = c(1,2), color = "gray") + 
  scale_x_continuous(breaks = seq(0.5, 2.5, 1), labels = c("5'UTR", "CDS", "3'UTR")) +
  scale_y_continuous(limits = c(0, 100)) + 
  scale_color_manual(values = mod_color) + 
  labs(color = "Modification") + 
  xlab(NULL) + ylab("Mean Modification Level (%)") + 
  theme_classic() + 
  theme(aspect.ratio = 1, legend.position = "top") + 
  ggtitle("polyA RNA")
ggsave(filename = paste0(outdir, "Mean_ratio_along_transripts.pdf"), width = 4, height = 5)
