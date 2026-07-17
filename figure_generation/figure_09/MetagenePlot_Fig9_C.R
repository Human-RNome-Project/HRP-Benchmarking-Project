# ============================================================
# Note: 
# 1. Use Ensembl_canonical transcript as representative transcript for each gene
# 2. Duplicated annotations will exist if multiple genes overlapped
# =============================================================

rm(list=ls()); gc()
library(tidyverse)
library(GenomicFeatures)
library(qs)
library(dplyr)

datadir = "/home/johannes/RNome/FINAL_BEDRMOD_FINAL/"
outdir = "/home/johannes/RNome/MetagenePlot/"

# ===========================================================
# 1. read gtf 
# ===========================================================
## get ensemble canonical transcripts as representative transcript by tag in GTF
# zcat gencode.v49.primary_assembly.annotation.gtf.gz | \
#   awk -v OFS="\t" '$3=="transcript"' | \
#   grep "Ensembl_canonical" | \
#   > gencode.v49.primary_assembly.annotation.transcripts.ensembl.canonical.gtf
fcano_trans = "/home/aaron/genome/RNome_hg38/my.gencode.v49/gencode.v49.primary_assembly.annotation.transcripts.ensembl.canonical.gtf"
cano_trans = rtracklayer::import(fcano_trans) %>%
  as.data.frame() %>%
  dplyr::select(transcript_id, gene_id, gene_name, gene_type)

## read gtf (slow!)
txdb <- makeTxDbFromGFF("/home/aaron/genome/RNome_hg38/my.gencode.v49/gencode.v49.primary_assembly.annotation.gtf.gz")
# Total: 509650 transcripts; 78899 genes

## get tx length, utr5, ytr3, cds length
txlen = transcriptLengths(txdb, with.utr5_len = T, with.cds_len = T, with.utr3_len = T) %>% 
  as.data.frame() %>% 
  filter(tx_name %in% cano_trans$transcript_id)
cano_trans = inner_join(cano_trans, txlen[, -1], 
  by = join_by("transcript_id" == "tx_name", "gene_id" == "gene_id"))

## get genomic position of transcript
gr_tx = transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id"))
gr_tx = gr_tx[gr_tx$tx_name %in% cano_trans$transcript_id]

## get exons of transcripts
exon_by_tx = exonsBy(txdb, by = "tx", use.names = T)
table(names(exon_by_tx) %in% cano_trans$transcript_id)
exon_by_tx = exon_by_tx[cano_trans$transcript_id]

ANNO = list(cano_trans = cano_trans, exon_by_tx = exon_by_tx, gr_tx = gr_tx)
qs::qsave(ANNO, file = paste0(outdir, "Anno_database_canonical_transcripts.qs"))
rm(cano_trans, txlen, exon_by_tx, gr_tx, txdb)

ANNO = qs::qread("~/RNome/MetagenePlot/Anno_database_canonical_transcripts.qs")

# =======================================================================
# 2. read mod sites
# =======================================================================
valid_chr = as.character(unique(unlist(seqnames(ANNO$exon_by_tx))))

file = "~/RNome/FINAL_BEDRMOD_FINAL/ONT_polyARNA_rRNA_tRNA_combined.filtered_rmchrY.bed"
bed = data.table::fread(file, data.table = F, check.names = T) %>% 
  filter(X.chrom %in% valid_chr) %>%
  mutate(ID = paste(X.chrom, chromEnd, strand, sep = "_")) %>%
  dplyr::select(X.chrom, chromEnd, name, strand, ID)

gr_bed = makeGRangesFromDataFrame(bed, seqnames.field = "X.chrom", 
  start.field = "chromEnd", end.field = "chromEnd", strand.field = "strand")

## anno to transcripts
raw_anno = findOverlaps(gr_bed, ANNO$gr_tx) %>%
  as.data.frame() %>%
  mutate(ID = bed$ID[queryHits], transcript_id = ANNO$gr_tx$tx_name[subjectHits])
df_anno = left_join(bed, raw_anno[, c("ID", "transcript_id")], by = "ID") %>%
  mutate(region = ifelse(is.na(transcript_id), "Intergenic", "Intronic"))

## anno to exon
df_exon = mapToTranscripts(gr_bed, transcripts = ANNO$exon_by_tx) %>% 
  as.data.frame()
df_exon = cbind(df_exon[, c("seqnames", "end")], ID = bed$ID[df_exon$xHits])

## calculate relative location: utr5 (0-1), cds (1-2), utr3 (2-3)
df_exon = left_join(df_exon, ANNO$cano_trans, by = join_by("seqnames" == "transcript_id")) %>% 
  mutate(
    rel_location = case_when(
      end <= utr5_len ~ end / utr5_len, 
      (end > utr5_len & end <= utr5_len + cds_len) ~ (end - utr5_len) / cds_len + 1,
      (end > (utr5_len + cds_len) & end <= tx_len) ~ (end - utr5_len - cds_len) / utr3_len + 2, 
      .default = NaN 
    )
  )

## combine exon annotation to transcript annotation
df_anno = left_join(
    df_anno, df_exon, 
    by = join_by("transcript_id" == "seqnames", "ID" == "ID")
  ) %>% 
  mutate(
    region = ifelse(!is.na(tx_len), "Exonic", region)
  )

## fix gene id for intron
idx = match(df_anno$transcript_id, ANNO$cano_trans$transcript_id)
df_anno = df_anno %>% mutate(
  gene_id = ANNO$cano_trans$gene_id[idx], 
  gene_name = ANNO$cano_trans$gene_name[idx], 
  gene_type = ANNO$cano_trans$gene_type[idx]
)

table(df_anno$name, df_anno$region)
table(df_anno$rel_location > 3)  # non-coding genes? (tx_len > 0; utr5,cds,utr3=0)
table(is.na(df_anno$rel_location))

data.table::fwrite(df_anno, file = paste0(outdir, "Illumina_combined_exon_region.tsv"), sep = "\t")


# =======================================================================
# 3. metagene plot (Nature / Human RNome style)
# =======================================================================
library(patchwork)
# df_anno = data.table::fread(paste0(outdir, "Illumina_combined_exon_region.tsv"), data.table = F)

## Official MODIFICATION_COLORS (Human RNome / rnome_style.py). Do not improvise.
## NB: I belongs to the A family -> official colour #98343E (was wrongly green).
mod_color = c(
  "m6A" = "#721817",
  "m5C" = "#001427",
  "I"   = "#98343E",
  "Y"   = "#F0A202",
  "Am"  = "#D44F3E",
  "Cm"  = "#0D3B6E",
  "Gm"  = "#74B354",
  "Um"  = "#C47A02"
)

## Two modification groups, each ordered canonically by base (A,C,G,U)
grp_base    = c("I", "m6A", "m5C", "Y")     # base modifications
grp_methyl  = c("Am", "Cm", "Gm", "Um")     # 2'-O-methylations

## Nature / Human RNome theme: Arial, 5-7 pt text, black axes, no coloured text.
theme_rnome = function(base_size = 6) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text         = element_text(colour = "black", size = base_size),
      axis.text    = element_text(colour = "black", size = base_size),
      axis.line    = element_line(colour = "black", linewidth = 0.3),
      axis.ticks   = element_line(colour = "black", linewidth = 0.3),
      axis.title   = element_text(colour = "black", size = 7),
      plot.title   = element_text(colour = "black", size = 7, hjust = 0),
      legend.title = element_text(colour = "black", size = 7),
      legend.text  = element_text(colour = "black", size = 6),
      legend.key.size = unit(3, "mm"),
      plot.tag     = element_text(size = 8, face = "bold"),
      aspect.ratio = 1
    )
}

make_metagene = function(df, mods) {
  d = df %>% filter(region == "Exonic", name %in% mods) %>%
    mutate(name = factor(name, levels = mods))
  ggplot(d, aes(rel_location, color = name)) +
    geom_density(fill = NA, linewidth = 0.5) +
    geom_vline(xintercept = c(1, 2), color = "gray70", linewidth = 0.3) +
    scale_x_continuous(breaks = seq(0.5, 2.5, 1),
                       labels = c("5'UTR", "CDS", "3'UTR")) +
    scale_color_manual(values = mod_color, limits = mods) +
    labs(color = "Modification", x = NULL, y = "Density") +
    theme_rnome()
}

p1 = make_metagene(df_anno, grp_base)    # base modifications
p2 = make_metagene(df_anno, grp_methyl)  # 2'-O-methylations

## combined two-panel figure (side by side, 183 mm double column)
fig = p1 + p2

## Vector output for submission (TrueType-42 fonts embedded via cairo_pdf).
ggsave(paste0(outdir, "metagenplot.pdf"), fig,
       device = cairo_pdf, width = 183, height = 80, units = "mm")
## PNG preview only (Nature does not accept png for submission).
ggsave(paste0(outdir, "metagenplot.png"), fig,
       width = 183, height = 80, units = "mm", dpi = 600)

## individual panels (single-column width)
ggsave(paste0(outdir, "metagenplot_base.pdf"), p1,
       device = cairo_pdf, width = 89, height = 70, units = "mm")
ggsave(paste0(outdir, "metagenplot_methyl.pdf"), p2,
       device = cairo_pdf, width = 89, height = 70, units = "mm")
ggsave(paste0(outdir, "metagenplot_base.png"), p1,
       width = 89, height = 70, units = "mm", dpi = 600)
ggsave(paste0(outdir, "metagenplot_methyl.png"), p2,
       width = 89, height = 70, units = "mm", dpi = 600)