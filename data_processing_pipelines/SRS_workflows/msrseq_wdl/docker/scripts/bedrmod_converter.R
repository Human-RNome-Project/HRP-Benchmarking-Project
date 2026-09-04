#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(openxlsx)
  library(Biostrings)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  cat("Usage: bedrmod_converter.R <data_cleaned_zip> <chromosomal_fasta> <isodecoder_table> <score_script> <out_dir>\n")
  quit(status = 1)
}

input_zip          <- args[1]
ref_seq_file       <- args[2]
mod_positions_file <- args[3]
score_script       <- args[4]
out_dir            <- args[5]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(score_script)
`%nin%` = Negate(`%in%`)

# Read reference sequence for chromosomal data
ref_seq <- read.table(
  ref_seq_file,
  sep = "",
  stringsAsFactors = FALSE,
  fill = TRUE
) %>%
  as_tibble() %>%  
  filter(startsWith(V1, ">Homo_sapiens_tRNA")) %>%
  mutate(V1 = sub(".*?-", "", V1)) %>%
  select(-V2,-V3,-V4,-V5,-V6,-V8,-V9,-V13) %>%
  separate(V11, into = c("chrom", "pos"), sep = ":") %>%
  separate(pos, into = c("chromStart_gene", "chromEnd_gene"), sep = "-") %>%
  dplyr::rename(
    gene       = V1,
    bp         = V7,
    scan_score = V10,
    strand     = V12
  )

# Define modification map
mod_map <- c(
  `"` = "m1A",
  I   = "I",
  O   = "m1I",
  M   = "ac4C",
  `'` = "m3C",
  `>` = "f5C",
  K   = "m1G",
  L   = "m2G",
  `#` = "Gm",
  R   = "m2,2G",
  `7` = "m7G",
  Q   = "Q",
  W   = "o2yW",
  `2` = "s2U",
  X   = "acp3U",
  P   = "Y"
)

# Read modification positions table
mod_positions <- read.xlsx(mod_positions_file) %>%
  dplyr::rename(gene = Anticodon) %>%
  dplyr::select(gene, sequence) %>%
  mutate(
    gene = paste0(
      substr(gene, 1, 3), "-",        
      substr(gene, 4, nchar(gene)),   
      "-1"                            
    )
  ) %>%
  mutate(
    pos = str_split(sequence, "", simplify = TRUE) %>% as.data.frame()
  ) %>%
  unnest_wider(pos) %>%
  dplyr::rename_with(~ gsub("^pos\\$V", "", .x), -c(gene, sequence)) %>%
  dplyr::rename_with(~ gsub("V", "", .x), -c(gene, sequence)) %>%
  dplyr::select(-sequence) %>%
  pivot_longer(
    cols = -gene,
    names_to = "position",
    values_to = "nucleotide"
  ) %>%
  mutate(
    position   = as.integer(position),
    nucleotide = na_if(nucleotide, "")
  ) %>%
  na.omit() %>%
  mutate(nucleotide = dplyr::recode(nucleotide, !!!mod_map)) %>%
  mutate(gene = ifelse(gene == "Sec-TCA-1-1", "SeC-TCA-1-1", gene),
         gene = ifelse(gene == "iMe-tCAT-1-1", "iMet-CAT-1-1", gene),
         gene = ifelse(gene == "iMe-tCAT-2-1", "iMet-CAT-2-1", gene),
         changed = ifelse((nucleotide == "Y" & position %in% c(55,54)), TRUE, FALSE),
         position = ifelse((changed == TRUE), position - 1, position))

write_bedrmod <- function(df, file) {
  metadata <- c(
    "#fileformat=bedRModv2",
    "#organism=9606",
    "#modification_type=RNA",
    "#modification_names=m1A:m1A:A,I:I:A,m1I:m1I:A,ac4C:ac4C:C,m3C:m3C:C,f5C:f5C:C,m1G:m1G:G,m2G:m2G:G,Gm:Gm:G,m22G:m22G:G,m7G:m7G:G,Q:Q:G,o2yW:o2yW:G,s2U:s2U:U,acp3U:acp3U:U,Y:Y:U",
    "#assembly=GRCh38",
    "#annotation_source=tRNAscan-SE",
    "#annotation_version=1",
    "#sequencing_platform=Illumina NovaSeqX",
    "#basecalling=",
    "#bioinformatics_workflow=https://github.com/Luke-F1875/MSRseq_data_processing_pipeline",
    "#experiment=https://doi.org/10.1038/s41467-022-30261-3",
    "#external_source="
  )

  col_header <- paste0("#", paste(colnames(df), collapse = "\t"))
  writeLines(c(metadata, col_header), con = file)

  write.table(
    df,
    file = file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    append = TRUE
  )
}

message("Reading modification data from: ", input_zip)
# Extract inner csv name if zip or read directly if csv
if (grepl("\\.zip$", input_zip)) {
  zip_files <- unzip(input_zip, list = TRUE)$Name
  target_csv <- grep("\\.csv$", zip_files, value = TRUE)[1]
  modification_data <- read.csv(unz(input_zip, target_csv))
} else {
  modification_data <- read.csv(input_zip)
}

process_rep <- function(target_rep) {
  message("Processing Replicate: ", target_rep)
  
  modification_data_fil <- modification_data %>%
    filter(pileup >= 50) %>%
    filter(bin_start == "60") %>%
    mutate(deletion_rate = deletion/pileup) %>%
    filter(rep == target_rep) %>%
    mutate(mutation_raw = round(pileup * mutation)) %>%
    select(-X, -matches("^[TCGAN]$"), -stop, -insertion, -bin_stop, -bin_start, -rep)

  # Mutation scores
  mod_mut_score <- modification_data_fil %>%
    dplyr::rename(n_mod = mutation_raw, n_total = pileup)
  mod_data_scores_score <- calculate_bedrmod_score_stoichiometry(mod_mut_score) %>%
    dplyr::rename(mut_score = score)
  mut_scores_min <- mod_data_scores_score %>%
    select(gene, position, treatment, mut_score)

  # Deletion scores
  mod_del_score <- modification_data_fil %>%
    dplyr::rename(n_mod = deletion, n_total = pileup)
  del_data_scores_score <- calculate_bedrmod_score_stoichiometry(mod_del_score) %>%
    dplyr::rename(del_score = score)
  del_scores_min <- del_data_scores_score %>%
    select(gene, position, treatment, del_score)

  # Merge
  final_data <- modification_data_fil %>%
    left_join(mut_scores_min, by = c("gene", "position", "treatment")) %>%
    left_join(del_scores_min, by = c("gene", "position", "treatment")) %>%
    filter(mutation >= .05 | deletion_rate >= .05)

  df_merged <- left_join(final_data, ref_seq, by = "gene")

  df_all <- left_join(mod_positions, df_merged, by = c("gene", "position")) %>%
    na.omit() %>%
    filter(!nucleotide %in% c("A", "G", "U", "C")) %>%
    mutate(position = ifelse((changed == TRUE), position + 1, position)) %>%
    dplyr::rename(name = nucleotide) %>%
    mutate(chrom = gsub("chr", "", chrom)) %>%
    mutate(itemRgb = '0,0,0') %>%
    dplyr::rename(coverage = pileup) %>%
    mutate(chromStart = as.numeric(chromStart_gene) + as.numeric(position) - 1,
           chromEnd = as.numeric(chromStart_gene) + as.numeric(position)) %>%
    dplyr::select(-chromEnd_gene, -chromStart_gene, -source, -AA, -anticodon, -file_name) %>%
    mutate(thickStart = chromStart,
           thickEnd = chromEnd,
           score = coverage) %>%
    select(chrom, chromStart, chromEnd, name, mut_score, del_score, strand, thickStart, thickEnd, itemRgb, coverage, mutation, deletion_rate, treatment)

  # BS
  bseq <- df_all %>%
    filter(treatment == "HRPC_BS", deletion_rate >= .05) %>%
    dplyr::select(-mut_score, -mutation, -treatment) %>%
    dplyr::rename(score = del_score, frequency = deletion_rate) %>%
    mutate(frequency = frequency * 100) %>%
    dplyr::select(chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb, coverage, frequency)

  # CBH
  cbhseq <- df_all %>%
    filter(treatment == "HRPC_CBH", mutation >= .05) %>%
    dplyr::select(-del_score, -deletion_rate, -treatment) %>%
    dplyr::rename(score = mut_score, frequency = mutation) %>%
    mutate(frequency = frequency * 100) %>%
    dplyr::select(chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb, coverage, frequency)

  # ctrl mutation
  ctrlseq_mut <- df_all %>%
    filter(treatment == "HRPC_ctrl", mutation >= .05) %>%
    dplyr::select(-del_score, -deletion_rate, -treatment) %>%
    dplyr::rename(score = mut_score, frequency = mutation) %>%
    mutate(frequency = frequency * 100) %>%
    dplyr::select(chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb, coverage, frequency)

  # ctrl deletion
  ctrlseq_del <- df_all %>%
    filter(treatment == "HRPC_ctrl", deletion_rate >= .05) %>%
    dplyr::select(-mut_score, -mutation, -treatment) %>%
    dplyr::rename(score = del_score, frequency = deletion_rate) %>%
    mutate(frequency = frequency * 100) %>%
    dplyr::select(chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb, coverage, frequency)

  # Write files
  write_bedrmod(bseq, file.path(out_dir, paste0("HRPC_BS_", target_rep, ".bedRMod")))
  write_bedrmod(cbhseq, file.path(out_dir, paste0("HRPC_CBH_", target_rep, ".bedRMod")))
  write_bedrmod(ctrlseq_mut, file.path(out_dir, paste0("HRPC_ctrl_mut_", target_rep, ".bedRMod")))
  write_bedrmod(ctrlseq_del, file.path(out_dir, paste0("HRPC_ctrl_del_", target_rep, ".bedRMod")))
}

for (r in 1:3) {
  process_rep(r)
}

message("Done converting all 12 bedRMod files to: ", out_dir)
