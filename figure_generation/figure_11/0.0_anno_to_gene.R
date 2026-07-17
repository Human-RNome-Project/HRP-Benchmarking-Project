# ============================================================
# Note: 
# 1. Use Ensembl_canonical transcript as representative transcript for each gene
# 2. Duplicated annotations will exist if multiple genes overlapped
# =============================================================

rm(list=ls()); gc()
library(tidyverse)
library(GenomicFeatures)

datadir = "data/My/Integrated/bedRmod_final/"
outdir = "data/My/Integrated/Illumina_polyA/res/"

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

# =======================================================================
# 2. read mod sites
# =======================================================================
valid_chr = as.character(unique(unlist(seqnames(ANNO$exon_by_tx))))

file = paste0(datadir, "Illumina_combined_polyARNA_tRNA_rRNA.bed")
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
