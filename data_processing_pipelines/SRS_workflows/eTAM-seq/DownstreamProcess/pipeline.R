#####################################################################
# Description: Differential m6A Analysis - Rnome mrna vs ivt        #
# Date: 2026/02                                                     #
# Author: Frederic MURISIER                                         #
#                                                                   #
# Load and prepare the m6A methylation data                         #
# Identify the sequences motifs associated to m6A                   #
# Do a differential analysis between ivt and ctrl with edgeR        #
# Identify the statistically significant different methylated sites #
# Export the results in txt and tsv                                 #
#####################################################################

# ==================
# Libraries
# ==================
library(GenomicRanges)
library(dplyr)
library(stringr)
library(Biostrings)
library(tidyverse)
library(rtracklayer)
library(purrr)
library(edgeR)
library(limma)
library(ggplot2)
library(txdbmaker)
library(data.table)
library(Rsamtools)
library(GenomicFeatures)
library(org.Dm.eg.db)
library(AnnotationDbi)
library(writexl)

# ==================
# Resolution conflicts
# ==================
count <- dplyr::count
select <- dplyr::select
filter <- dplyr::filter

# ==================
# Parameters
# ==================
p_local <- '/home/h/OneDrive/Unil/eTam/Data/rnome/to_transfert/empty_folder_to_fill/FOLDER_JY/JYR/'
PATH <- "/FOLDER/JYR/"
PATH <- p_local
IVT <- 'IVT' #IVT use to generate the delta_IVT
nIVT <- 2
CTRL <- 'mRNA'
nCTRL <- 3
TYPE <- 'ETAM_RNA_mrna/metadata/'
TYPE <- paste0(p_local, 'ETAM_RNA_mrna/metadata/')
TYPE2 <- paste0(p_local, 'ETAM_RNA_IVT/metadata/')
EXT <- '.tsv' #
pattern_to_remove <- 'metadata_form_'
MIN_COVERAGE <- 10                   # Min coverage of read
MIN_SAMPLES <- min(nIVT,nCTRL)       # Nombre minimal d'échantillons
FDR_THRESHOLD <- 0.05                # FDR Threshold (False Discovery Rate)
DELTA_THRESHOLD <- 0.1               # Minimal methylation difference

# Motif patterns
# DRACH = consensus motif for m6A (D=A/G/T, R=A/G, H=A/C/T)
DRACH <- c("AAACA","AAACC","AAACT","AGACA","AGACC","AGACT",
           "GAACA","GAACC","GAACT","GGACA","GGACC","GGACT",
           "TAACA","TAACC","TAACT","TGACA","TGACC","TGACT")
# DRAC = short variant (4 nucléotides)
DRAC <- c("AAAC","AGAC","GAAC","GGAC","TAAC","TGAC")
# RAC = shorter motif  (3 nucléotides)
RAC <- c("AAC", "GAC")
# RRAT = alternatif motif
RRAT <- c("AAAT", "AGAT", "GAAT", "GGAT")

# Patterns for detection
DRACH_PATTERN <- paste0("^(", paste(DRACH, collapse = "|"), ")")
DRAC_PATTERN <- paste0("^(", paste(DRAC, collapse = "|"), ")")  # Pattern starting with DRAC
RAC_PATTERN <- paste0("^.(", paste(RAC, collapse = "|"), ")")   # Pattern with RAC in position 2
RRAT_PATTERN <- "^(A|G)(A|G)AT"                                 # Pattern for RRAT

# ==================
# PARAMETRES DE STYLE
# ==================
x_font <- 14
FONT_FAMILY <- "sans"
FONT_SIZE_TITLE <- x_font
FONT_SIZE_AXIS_TITLE <- x_font-2
FONT_SIZE_AXIS_TEXT <- x_font-3
FONT_SIZE_LEGEND <- x_font-4
FONT_SIZE_GEOM_TEXT <- x_font-10
FONT_COLOR <- "black"

# ==================
# Theme custom
# ==================
theme_custom <- function() {
  theme_bw() +                     # Theme classic black and white
    theme(
      text = element_text(family = FONT_FAMILY, color = FONT_COLOR),
      plot.title = element_text(size = FONT_SIZE_TITLE, face = "bold",
                                color = FONT_COLOR, hjust = 0.5),  # Title center and bold
      axis.title = element_text(size = FONT_SIZE_AXIS_TITLE, color = FONT_COLOR),
      axis.text = element_text(size = FONT_SIZE_AXIS_TEXT, color = FONT_COLOR),
      legend.text = element_text(size = FONT_SIZE_LEGEND, color = FONT_COLOR),
      legend.title = element_blank(),                   # No legend title
      legend.key.size = unit(0.8, 'cm'),                # Legend key size
      panel.grid.major = element_blank(),               # Remove major panel grid lines
      panel.grid.minor = element_blank(),               # Remove minor panel grid lines
      axis.line = element_line(colour = "black"),       # Axe line black
      panel.border = element_blank()                    # No panel border
    )
}

# ==================
# Setup directories
# ==================
path <- paste0(PATH)
path
setwd(path)

# Create output directories
dir.create("plots", showWarnings = FALSE)
dir.create("input_DE", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

cat("\n=== Working directory ===\n")
cat("Path:", path, "\n")

# ==================
# Helper Functions
# ==================

# Rename files with hyphens
rename_files_with_hyphens <- function(folder) {
  # List all files .annot.txt
  files <- list.files(path = folder, pattern = EXT, full.names = TRUE)
  files
  # Select the ones with hyphen
  files_with_hyphen <- files[grep("-", files)]

  # Rename every files (remplace - par _)
  for (file in files_with_hyphen) {
    new_name <- gsub("-", "_", file)
    file.rename(file, new_name)
    cat("Renamed:", basename(file), "->", basename(new_name), "\n")
  }
  # Return the updated list
  list.files(path = folder, pattern = EXT, full.names = TRUE)
}

# Import annotation files
import_annot_files <- function(files) {
  for (file in files) {
    name <- gsub(EXT, "", basename(file))  # name without extension
    data <- read.table(file, header = TRUE, quote = "", sep = "\t")  # reading file
    assign(name, data, envir = .GlobalEnv)  # Assigne in the global environnement
  }
  cat("\n=== Imported", length(files), "annotation files ===\n")
}

# Add motif columns
add_motif_columns <- function(df) {
  df %>%
    mutate(
      protected = nbA / nbA.G,  # Protection ratio (méthylation)
      motif = substr(context, 2, nchar(context) - 1),  # Extraction of the motif from the context
      # Detection of the differents motifs
      DRACH = str_detect(motif, DRACH_PATTERN),
      DRAC = str_detect(motif, DRAC_PATTERN),
      RAC = str_detect(motif, RAC_PATTERN),
      RRAT = str_detect(motif, RRAT_PATTERN),
      coord = paste0(chr, "-", pos),  # Coordonninates chromosome-position
      strand = case_when(  # Identification of the strand
        nt == "A" ~ "+",
        nt == "T" ~ "-",
        TRUE ~ ""
      )
    )
}

# Merge protected columns across samples
merge_protected_columns <- function(samples) {
  merged_df <- NULL

  for (sample in samples) {
    df <- get(sample, envir = .GlobalEnv)  # Get the dataframe
    temp_df <- df[, c("coord", "protected")]  # Select coord et protected
    colnames(temp_df)[2] <- sample  # Rename the column protected

    if (is.null(merged_df)) {
      merged_df <- temp_df  # First Sample
    } else {
      merged_df <- merge(merged_df, temp_df, by = "coord", all = FALSE)  # Fusion
    }
  }

  rownames(merged_df) <- merged_df$coord
  merged_df[, -1]  # Return without the column coord
}

# Export for edgeR
export_for_edgeR <- function(sample, output_dir) {
  df <- get(sample, envir = .GlobalEnv)

  # Split the coordinate chromosome-position
  coord_split <- strsplit(as.character(df$coord), "-")
  new_df <- data.frame(
    V1 = sapply(coord_split, `[`, 1),  # Chromosome
    V2 = as.numeric(sapply(coord_split, `[`, 2)),  # Position start
    V3 = as.numeric(sapply(coord_split, `[`, 2)) + 1,  # Position end
    V4 = df$protected,  # Ratio méthylation (frequence)
    V5 = df$nbA,  # Number of A
    V6 = df$nbG  # Number of G
  )

  output_file <- file.path(output_dir, paste0(gsub(pattern_to_remove, "", sample), ".txt"))
  write.table(new_df, file = output_file, sep = "\t",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Write file in bedRMod format
write_bedRMod <- function(df, output_file) {

  # Vérifications minimales
  required_cols <- c("coord", "chr", "pos", "pvalue", "protected", "strand","nbA.G")
  missing_cols <- setdiff(required_cols, colnames(df))

  if (length(missing_cols) > 0) {
    stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
  }

  # Construction du dataframe bedRMod
  bedmod_df <- data.frame(
    chrom = df$chr,
    chromStart = df$pos - 1,
    chromEnd = df$pos,
    name = df$coord,
    score = df$pvalue,
    strand = df$strand,
    thickStart = df$pos - 1,
    thickEnd = df$pos,
    itemRgb = "0,0,0",
    coverage = df$nbA.G,
    frequency = df$protected,
    stringsAsFactors = FALSE
  )

  # Header metadata
  header_lines <- c(
    "#fileformat=bedRModv2",
    "#organism=9606",
    "#modification_type=RNA",
    "#modification_names=21891:m6A:A",
    "#assembly=GRCh38",
    "#annotation_source=Ensembl",
    "#annotation_version=93",
    "#sequencing_platform=Illumina NovaSeq 6000",
    "#chrom chromStart chromEnd name score strand thickStart thickEnd itemRgb coverage frequency"
  )

  # writing file
  con <- file(output_file, open = "wt")

  writeLines(header_lines, con)

  write.table(
    bedmod_df,
    con,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  close(con)

  message("bedRMod file successfully written to ", output_file)
}

# ==================
# Load and process data
# ==================

# Rename and import files
files1 <- rename_files_with_hyphens(TYPE)  # Rename files
import_annot_files(files1)  # Importe all files
files1
# Rename and import files
files2 <- rename_files_with_hyphens(TYPE2)  # Rename files
import_annot_files(files2)  # Importe all files
files2

files <- append(files1,files2)

# Get sample names
samples <- gsub(EXT, "", basename(files))
samples2 <- gsub(pattern_to_remove, "", samples)

cat("\n=== Samples loaded ===\n")
print(samples2)

# ==================
# Add motif annotations
# ==================

for (sample_name in samples) {
  sample_df <- get(sample_name, envir = .GlobalEnv)
  sample_df <- add_motif_columns(sample_df)  # Add motif columns
  print(summary(sample_df$coverage))
  assign(sample_name, sample_df, envir = .GlobalEnv)  # Reassigne
}



# ==================
# Add motif coverage,conv rate, freq, pval and fdr
# ==================

p_err <- 0.001

for (sample_name in samples) {
  sample_df <- get(sample_name, envir = .GlobalEnv)
  sample_df$coverage <- sample_df$nbA + sample_df$nbG
  sample_df$conv_rate <- sample_df$nbG / sample_df$coverage
  sample_df$pvalue <- mapply(function(k, N) {
    binom.test(k, N, p = p_err, alternative = "greater")$p.value
  }, sample_df$nbG, sample_df$coverage)

  sample_df$fdr <- p.adjust(sample_df$pvalue, method = "BH")

  assign(sample_name, sample_df, envir = .GlobalEnv)  # Reassigne
}



cat("\n=== Motif annotations added ===\n")

# ==================
# Quality control plots
# ==================

# Count A sites per sample
summary_df <- data.frame(
  Sample = samples2,
  nbAsites = sapply(samples, function(x) nrow(get(x, envir = .GlobalEnv)))
)

test1 <- get(row.names(summary_df)[1])
test1 <- metadata_form_mRNA_3
min(test1$nbA.G)# 10  which means they have been filtered for 10 reads of coverage
max(test1$protected) # shouldn always be =1
summary(test1$protected)
test1_df_sub <- test1[test1$protected >= 0.5,]
nrow(test1)
nrow(test1_df_sub)
colnames(test1)

test2 <- metadata_form_IVT_3
summary(is.na(test1$coverage))
test2_df_sub <- test2[test2$protected >= 0.5,]
nrow(test2)
nrow(test2_df_sub)
colnames(test2)

# Plot number of A sites
p_sites <- ggplot(summary_df, aes(x = Sample, y = nbAsites / 1e6)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = round(nbAsites / 1e6, 2)), vjust = -0.5,
            size = FONT_SIZE_GEOM_TEXT) +
  labs(title = "A sites with coverage ≥ 10 reads",
       x = "", y = "A sites (M)") +
  theme_custom() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_sites
x_size <- 15
ggsave(paste0(path,"/plots/nb_A_sites_per_sample.pdf"), p_sites,
       width = x_size, height = x_size-2, bg = "transparent", dpi = 300)


for (sample_name in samples) {
  sample_df <- get(sample_name, envir = .GlobalEnv)
  sample_name_all <- sub('metadata_form_','',sample_name)
  write_bedRMod(df = sample_df,output_file = paste0(sample_name_all,"_all_bedRMod.bed"))
  sample_df_sub <- sample_df[sample_df$pvalue <= 0.05,]
  write_bedRMod(df = sample_df_sub,output_file = paste0(sample_name_all,"_pval0.05_bedRMod.bed"))
}


# ==================
# Merge data across samples
# ==================

merged_data <- merge_protected_columns(samples)
cat("\n=== Merged data dimensions ===\n")
cat("Sites:", nrow(merged_data), "Samples:", ncol(merged_data), "\n")
head(merged_data)

# ==================
# Create pdata
# ==================
grp <- strsplit(summary_df$Sample,'_')
grp_names <- sapply(grp, `[`, 1)
grp_names

pdata <- summary_df %>%
  mutate(
    Group = factor(grp_names, levels = sort(unique(grp_names)))
  )

write.csv(pdata, file = "pdata.csv", row.names = FALSE)

# PCA plot
color_plot <- c("black", "maroon1","red")
col.cell <- color_plot[pdata$Group]
plotMDS(merged_data, gene.selection = "common", col = col.cell, cex = 1.2)
legend("bottomleft", fill = color_plot, legend = levels(pdata$Group))
png("plots/PCA_plot.png", width = 7, height = 7, res = 350, units = "in")
dev.off()

# ==================
# Export for edgeR
# ==================

cat("\n=== Exporting for edgeR ===\n")
lapply(samples, export_for_edgeR, output_dir = "input_DE")


# ==================
# Differential analysis with edgeR
# ==================

cat("\n=== Running differential analysis ===\n")

# Prepare data
pdata_sub <- pdata %>%
  filter(Group %in% c(CTRL,IVT)) %>%
  droplevels() %>%
  mutate(
    Group = relevel(Group, CTRL),
    File = file.path(paste0("input_DE/", Sample, EXT))
  )
pdata_sub

# ==================
# Merge results with annotations
# ==================

cat("\n=== Merging results with annotations ===\n")

# Get annotation from first control sample
annot_sample <- get(
  paste0('/home/h/OneDrive/Unil/eTam/Data/rnome/to_transfert/empty_folder_to_fill/FOLDER_JY/JYR/ETAM_RNA_mrna/raw/',
  'NextSeq2000_mrna_260211_R1/',pdata_sub$Sample[1]), envir = .GlobalEnv)
annot <- annot_sample[, c("coord", "chr", "pos", "nt", "kind", "gene",
                          "motif", "DRACH", "DRAC", "RAC", "RRAT", "strand")]

# Split the samples by groups
samples_ctrl <- pdata_sub %>% filter(Group == CTRL) %>% pull(Sample)
samples_ivt <- pdata_sub %>% filter(Group == IVT) %>% pull(Sample)

# merge the data of each groups
merged_ctrl <- merge_protected_columns(samples_ctrl)
merged_ivt <- merge_protected_columns(samples_ivt)
head(merged_ivt)
head(merged_ctrl)
nrow(top_ivt)



# Merge all data
qlf_df_2 <- top_ivt %>%
  tibble::rownames_to_column("coord") %>%
  left_join(annot, by = "coord") %>%
  left_join(merged_ctrl %>% tibble::rownames_to_column("coord"), by = "coord") %>%
  left_join(merged_ivt %>% tibble::rownames_to_column("coord"), by = "coord") %>%
  mutate(
    mean_ctrl = rowMeans(select(., starts_with(CTRL)), na.rm = TRUE),  # mean ctrl
    mean_IVT = rowMeans(select(., starts_with(IVT)), na.rm = TRUE),  # mean IVT
    delta_ivt = mean_ctrl - mean_IVT,  # Difference of methylation
  ) %>%
  arrange(desc(FDR))  # order by FDR descendent

head(qlf_df_2)
nrow(qlf_df_2)
colnames(qlf_df_2)

# ==================
# Merge data for each cdt
# ==================

qlf_df_ctrl <- top_ivt %>%
  tibble::rownames_to_column("coord") %>%
  left_join(annot, by = "coord") %>%
  left_join(merged_ctrl %>% tibble::rownames_to_column("coord"), by = "coord") %>%
  mutate(
    mean_ctrl = rowMeans(select(., starts_with(CTRL)), na.rm = TRUE),  # mean ctrl
  ) %>%
  arrange(desc(FDR))  # order by FDR descendent
nrow(qlf_df_ctrl)
colnames(qlf_df_ctrl)


qlf_df_IVT <- top_ivt %>%
  tibble::rownames_to_column("coord") %>%
  left_join(annot, by = "coord") %>%
  left_join(merged_ivt %>% tibble::rownames_to_column("coord"), by = "coord") %>%
  mutate(
    mean_IVT = rowMeans(select(., starts_with(IVT)), na.rm = TRUE),  # mean ivt
  ) %>%
  arrange(desc(FDR))  # order by FDR descendent
nrow(qlf_df_IVT)
colnames(qlf_df_IVT)

# ==================
# Filter significant sites
# ==================
qlf_sub_2 <- qlf_df_2 %>%
  filter(delta_ivt >= DELTA_THRESHOLD & PValue <= FDR_THRESHOLD) # Delta ≥ 0.1 et p ≤ 0.05
head(qlf_sub_2)
nrow(qlf_sub_2)
nrow(top_ivt)

#---------
# by cdt
#---------

qlf_sub_ctrl <- qlf_df_ctrl %>%
  filter(PValue <= FDR_THRESHOLD) # p ≤ 0.05
head(qlf_sub_ctrl)
nrow(qlf_sub_ctrl)

qlf_sub_ivt <- qlf_df_IVT %>%
  filter(PValue <= FDR_THRESHOLD) # p ≤ 0.05
head(qlf_sub_ivt)
nrow(qlf_sub_ivt)
head(qlf_df_IVT)
nrow(qlf_df_IVT)
summary(qlf_df_IVT$PValue)


#-------------
# by cdt end
#-------------

cat("\n=== Significant sites IVT ===\n")
cat("Total:", nrow(qlf_sub_2), "\n")
cat("Hypermethylated (", CTRL, " > ", IVT, "):", sum(qlf_sub_2$delta_ivt > 0), "\n")
cat("Hypomethylated (", CTRL, " < ", IVT, "):", sum(qlf_sub_2$delta_ivt < 0), "\n")

# ========================
# Export results bed Rmod
# ========================

write_bedRMod(
  df = qlf_df_ctrl,
  output_file = "mRNA_all_bedRMod.bed",
  mean_column = "mean_ctrl"
)
write_bedRMod(
  df = qlf_sub_ctrl,
  output_file = "mRNA_pV0.05_bedRMod.bed",
  mean_column = "mean_ctrl"
)


write_bedRMod(
  df = qlf_df_IVT,
  output_file = "IVT_all_bedRMod.bed",
  mean_column = "mean_IVT"
)
write_bedRMod(
  df = qlf_sub_ivt,
  output_file = "IVT_pV0.05_bedRMod.bed",
  mean_column = "mean_IVT"
)


cat("\n=== Pipeline complete! ===\n")


