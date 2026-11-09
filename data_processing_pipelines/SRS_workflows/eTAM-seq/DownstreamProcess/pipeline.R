#####################################################################
# Description: Differential m6A Analysis - Rnome mrna vs ivt        #
# Date: 2026/02 (rev. 2026/11)                                      #
#                                                                   #
# Load and prepare the m6A methylation data                         #
# Identify the sequences motifs associated to m6A                   #
# Do a differential analysis between ivt and ctrl with edgeR        #
# Identify the statistically significant different methylated sites #
# Export the results in tsv and bedRMod v2                          #
#                                                                   #
# Usage:                                                            #
#   Rscript pipeline.R <data_dir> [out_dir]                         #
#     data_dir : folder containing ETAM_RNA_mrna/metadata/ and      #
#                ETAM_RNA_IVT/metadata/                             #
#     out_dir  : where results are written (default: data_dir)      #
#####################################################################

# ==================
# Libraries
# ==================
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(data.table)
  library(edgeR)
  library(limma)
})

# ==================
# Resolution conflicts
# ==================
count <- dplyr::count
select <- dplyr::select
filter <- dplyr::filter

# ==================
# Arguments
# ==================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  DATA_DIR <- args[1]
  OUT_DIR <- if (length(args) >= 2) args[2] else DATA_DIR
} else if (interactive() && exists("DATA_DIR")) {
  # Interactive use: set DATA_DIR (and optionally OUT_DIR) before sourcing
  if (!exists("OUT_DIR")) OUT_DIR <- DATA_DIR
} else {
  stop("Usage: Rscript pipeline.R <data_dir> [out_dir]")
}
if (!dir.exists(DATA_DIR)) stop("data_dir not found: ", DATA_DIR)

# ==================
# Parameters
# ==================
CTRL <- 'mRNA'                       # Biological (methylated) condition
IVT <- 'IVT'                         # Unmethylated control, used to generate delta_ivt
CTRL_DIR <- file.path(DATA_DIR, 'ETAM_RNA_mrna', 'metadata')
IVT_DIR <- file.path(DATA_DIR, 'ETAM_RNA_IVT', 'metadata')
EXT <- '.tsv'
pattern_to_remove <- 'metadata_form_'
MIN_COVERAGE <- 10                   # Min coverage of read
FDR_THRESHOLD <- 0.05                # Significance threshold (applied to PValue)
DELTA_THRESHOLD <- 0.1               # Minimal methylation difference
P_NONCONV <- NA                      # Background non-conversion rate (A read as A without m6A)
                                     # for the per-site binomial test on nbA.
                                     # NA = estimated from the IVT samples (pooled nbA / nbA.G)

# bedRMod header (annotation = RefSeq GenBank records of GRCh38.p12, see PreProcess)
BEDRMOD_ORGANISM <- 9606
BEDRMOD_ASSEMBLY <- 'GRCh38.p12'
BEDRMOD_ANNOT_SOURCE <- 'RefSeq'
BEDRMOD_ANNOT_VERSION <- 'GCF_000001405.38_GRCh38.p12'
BEDRMOD_PLATFORM <- 'Illumina NextSeq 2000'

# Motif patterns
# DRACH = consensus motif for m6A (D=A/G/T, R=A/G, H=A/C/T)
DRACH <- c("AAACA","AAACC","AAACT","AGACA","AGACC","AGACT",
           "GAACA","GAACC","GAACT","GGACA","GGACC","GGACT",
           "TAACA","TAACC","TAACT","TGACA","TGACC","TGACT")
# DRAC = short variant (4 nucleotides)
DRAC <- c("AAAC","AGAC","GAAC","GGAC","TAAC","TGAC")
# RAC = shorter motif (3 nucleotides)
RAC <- c("AAC", "GAC")
# RRAT = alternative motif
RRAT <- c("AAAT", "AGAT", "GAAT", "GGAT")

# Patterns for detection
DRACH_PATTERN <- paste0("^(", paste(DRACH, collapse = "|"), ")")
DRAC_PATTERN <- paste0("^(", paste(DRAC, collapse = "|"), ")")  # Pattern starting with DRAC
RAC_PATTERN <- paste0("^.(", paste(RAC, collapse = "|"), ")")   # Pattern with RAC in position 2
RRAT_PATTERN <- "^(A|G)(A|G)AT"                                 # Pattern for RRAT

# ==================
# Style parameters
# ==================
x_font <- 14
FONT_FAMILY <- "sans"
FONT_SIZE_TITLE <- x_font
FONT_SIZE_AXIS_TITLE <- x_font-2
FONT_SIZE_AXIS_TEXT <- x_font-3
FONT_SIZE_LEGEND <- x_font-4
FONT_SIZE_GEOM_TEXT <- x_font-10
FONT_COLOR <- "black"

theme_custom <- function() {
  theme_bw() +
    theme(
      text = element_text(family = FONT_FAMILY, color = FONT_COLOR),
      plot.title = element_text(size = FONT_SIZE_TITLE, face = "bold",
                                color = FONT_COLOR, hjust = 0.5),
      axis.title = element_text(size = FONT_SIZE_AXIS_TITLE, color = FONT_COLOR),
      axis.text = element_text(size = FONT_SIZE_AXIS_TEXT, color = FONT_COLOR),
      legend.text = element_text(size = FONT_SIZE_LEGEND, color = FONT_COLOR),
      legend.title = element_blank(),
      legend.key.size = unit(0.8, 'cm'),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.border = element_blank()
    )
}

# ==================
# Setup directories
# ==================
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
for (d in c("plots", "input_DE", "tables", "bedRMod")) {
  dir.create(file.path(OUT_DIR, d), showWarnings = FALSE)
}

cat("\n=== Directories ===\n")
cat("Data:  ", DATA_DIR, "\n")
cat("Output:", OUT_DIR, "\n")

# ==================
# Helper Functions
# ==================

# List the metadata files of one folder
list_metadata_files <- function(folder) {
  if (!dir.exists(folder)) stop("Folder not found: ", folder)
  files <- list.files(folder, pattern = paste0("\\", EXT, "$"), full.names = TRUE)
  if (length(files) == 0) stop("No *", EXT, " file in ", folder)
  files
}

# Clean sample name: drop extension and prefix, hyphens -> underscores
clean_sample_name <- function(file) {
  name <- sub(paste0("\\", EXT, "$"), "", basename(file))
  name <- sub(pattern_to_remove, "", name, fixed = TRUE)
  gsub("-", "_", name)
}

# Add motif columns
add_motif_columns <- function(df) {
  df %>%
    mutate(
      protected = nbA / nbA.G,                          # Protection ratio (methylation)
      motif = substr(context, 2, nchar(context) - 1),   # Motif extracted from the context
      DRACH = str_detect(motif, DRACH_PATTERN),
      DRAC = str_detect(motif, DRAC_PATTERN),
      RAC = str_detect(motif, RAC_PATTERN),
      RRAT = str_detect(motif, RRAT_PATTERN),
      coord = paste0(chr, "-", pos),                    # chromosome-position
      strand = case_when(
        nt == "A" ~ "+",
        nt == "T" ~ "-",
        TRUE ~ "."
      )
    )
}

# Add coverage, conversion rate, binomial p-value and FDR
# The test is on nbA (A protected from deamination = m6A candidate): is the fraction of A
# reads higher than the background non-conversion rate p0?
add_site_stats <- function(df, p0) {
  df$coverage <- df$nbA + df$nbG
  df$conv_rate <- df$nbG / df$coverage
  # Identical to binom.test(nbA, coverage, p = p0, alternative = "greater")$p.value
  df$pvalue <- pbinom(df$nbA - 1, df$coverage, p0, lower.tail = FALSE)
  df$fdr <- p.adjust(df$pvalue, method = "BH")
  df
}

# Merge one column across samples (sites present in all samples)
merge_sample_column <- function(sample_names, column = "protected") {
  merged_df <- NULL
  for (s in sample_names) {
    temp_df <- sites[[s]][, c("coord", column)]
    colnames(temp_df)[2] <- s
    merged_df <- if (is.null(merged_df)) temp_df else
      merge(merged_df, temp_df, by = "coord", all = FALSE)
  }
  rownames(merged_df) <- merged_df$coord
  merged_df[, -1, drop = FALSE]
}

# Export for edgeR (Bismark coverage-like format)
export_for_edgeR <- function(sample_name, output_file) {
  df <- sites[[sample_name]]
  new_df <- data.frame(
    V1 = df$chr,             # Chromosome
    V2 = df$pos,             # Position start
    V3 = df$pos + 1,         # Position end
    V4 = df$protected,       # Methylation ratio (not used by readBismark2DGE)
    V5 = df$nbA,             # Methylated count (A kept)
    V6 = df$nbG              # Unmethylated count (A -> G)
  )
  fwrite(new_df, output_file, sep = "\t", col.names = FALSE, quote = FALSE)
}

# p-value -> bedRMod score: integer 0..1000, phred-like (-10 * log10 p), capped at 1000
pvalue_to_score <- function(p) {
  s <- round(-10 * log10(pmax(p, 1e-300)))
  as.integer(pmin(1000, pmax(0, ifelse(is.na(s), 0, s))))
}

# Write file in bedRMod v2 format
#   pvalue_col    : column turned into the score (0..1000)
#   frequency_col : methylation ratio (0..1), written as a percentage (0..100)
#   coverage_col  : read coverage at the site
write_bedRMod <- function(df, output_file, pvalue_col, frequency_col, coverage_col) {

  required_cols <- c("chr", "pos", "strand", pvalue_col, frequency_col, coverage_col)
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
  }

  keep <- !is.na(df[[frequency_col]]) & !is.na(df[[coverage_col]])
  if (any(!keep)) {
    message(sum(!keep), " sites without frequency/coverage dropped from ", basename(output_file))
  }
  df <- df[keep, , drop = FALSE]

  bedmod_df <- data.frame(
    chrom = df$chr,
    chromStart = df$pos - 1,
    chromEnd = df$pos,
    name = "m6A",
    score = pvalue_to_score(df[[pvalue_col]]),
    strand = ifelse(df$strand %in% c("+", "-"), df$strand, "."),
    thickStart = df$pos - 1,
    thickEnd = df$pos,
    itemRgb = "0,0,0",
    coverage = as.integer(round(df[[coverage_col]])),
    frequency = as.integer(round(100 * df[[frequency_col]])),
    stringsAsFactors = FALSE
  )

  header_lines <- c(
    "#fileformat=bedRModv2",
    paste0("#organism=", BEDRMOD_ORGANISM),
    "#modification_type=RNA",
    "#modification_names=21891:m6A:A",
    paste0("#assembly=", BEDRMOD_ASSEMBLY),
    paste0("#annotation_source=", BEDRMOD_ANNOT_SOURCE),
    paste0("#annotation_version=", BEDRMOD_ANNOT_VERSION),
    paste0("#sequencing_platform=", BEDRMOD_PLATFORM),
    paste0("#", paste(colnames(bedmod_df), collapse = "\t"))
  )

  con <- file(output_file, open = "wt")
  writeLines(header_lines, con)
  write.table(bedmod_df, con, sep = "\t",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  close(con)

  message("bedRMod file written: ", output_file)
}

# ==================
# Load and process data
# ==================
files <- c(list_metadata_files(CTRL_DIR), list_metadata_files(IVT_DIR))
sample_names <- clean_sample_name(files)
if (anyDuplicated(sample_names)) stop("Duplicated sample names: ",
                                      paste(sample_names[duplicated(sample_names)], collapse = ", "))

sites <- setNames(
  lapply(files, function(f) read.table(f, header = TRUE, quote = "", sep = "\t")),
  sample_names
)

cat("\n=== Samples loaded ===\n")
print(sample_names)

# ==================
# Add motif annotations and per-site statistics
# ==================
sites <- lapply(sites, add_motif_columns)

if (is.na(P_NONCONV)) {
  ivt_samples <- sample_names[sub("_.*$", "", sample_names) == IVT]
  if (length(ivt_samples) == 0) stop("No ", IVT, " sample to estimate P_NONCONV; set it by hand.")
  P_NONCONV <- sum(vapply(sites[ivt_samples], function(d) sum(d$nbA), numeric(1))) /
    sum(vapply(sites[ivt_samples], function(d) sum(d$nbA.G), numeric(1)))
}
cat(sprintf("\nBackground non-conversion rate (P_NONCONV): %.4f\n", P_NONCONV))

sites <- lapply(sites, add_site_stats, p0 = P_NONCONV)
cat("\n=== Motif annotations and per-site statistics added ===\n")

# ==================
# Sample sheet
# ==================
pdata <- data.frame(
  Sample = sample_names,
  nbAsites = vapply(sites, nrow, integer(1)),
  Group = sub("_.*$", "", sample_names),
  stringsAsFactors = FALSE
)
unknown <- setdiff(unique(pdata$Group), c(CTRL, IVT))
if (length(unknown) > 0) stop("Samples with a group other than ", CTRL, "/", IVT, ": ",
                              paste(pdata$Sample[pdata$Group %in% unknown], collapse = ", "))
pdata$Group <- factor(pdata$Group, levels = c(CTRL, IVT))
nCTRL <- sum(pdata$Group == CTRL)
nIVT <- sum(pdata$Group == IVT)
if (nCTRL < 2 || nIVT < 2) stop("At least 2 replicates per group are needed (", CTRL, ": ",
                                nCTRL, ", ", IVT, ": ", nIVT, ")")
MIN_SAMPLES <- min(nIVT, nCTRL)      # Minimal number of samples reaching MIN_COVERAGE

fwrite(pdata, file.path(OUT_DIR, "pdata.csv"))
print(pdata)

# ==================
# Quality control
# ==================
for (s in sample_names) {
  cat(sprintf("%-15s sites=%d  min coverage=%d  median protected=%.3f\n",
              s, nrow(sites[[s]]), min(sites[[s]]$nbA.G), median(sites[[s]]$protected)))
}

p_sites <- ggplot(pdata, aes(x = Sample, y = nbAsites / 1e6)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = round(nbAsites / 1e6, 2)), vjust = -0.5,
            size = FONT_SIZE_GEOM_TEXT) +
  labs(title = paste0("A sites with coverage ≥ ", MIN_COVERAGE, " reads"),
       x = "", y = "A sites (M)") +
  theme_custom() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

x_size <- 15
ggsave(file.path(OUT_DIR, "plots", "nb_A_sites_per_sample.pdf"), p_sites,
       width = x_size, height = x_size-2, bg = "transparent", dpi = 300)

# ==================
# Per-sample bedRMod export
# ==================
for (s in sample_names) {
  df <- sites[[s]]
  write_bedRMod(df, file.path(OUT_DIR, "bedRMod", paste0(s, "_all_bedRMod.bed")),
                pvalue_col = "pvalue", frequency_col = "protected", coverage_col = "nbA.G")
  write_bedRMod(df[df$pvalue <= 0.05, ],
                file.path(OUT_DIR, "bedRMod", paste0(s, "_pval0.05_bedRMod.bed")),
                pvalue_col = "pvalue", frequency_col = "protected", coverage_col = "nbA.G")
}

# ==================
# Merge data across samples + MDS plot
# ==================
merged_data <- merge_sample_column(sample_names, "protected")
cat("\n=== Merged data dimensions ===\n")
cat("Sites:", nrow(merged_data), "Samples:", ncol(merged_data), "\n")

color_plot <- c("black", "red")
png(file.path(OUT_DIR, "plots", "PCA_plot.png"), width = 7, height = 7, res = 350, units = "in")
plotMDS(as.matrix(merged_data), gene.selection = "common",
        col = color_plot[pdata$Group], cex = 1.2)
legend("bottomleft", fill = color_plot, legend = levels(pdata$Group))
dev.off()

# ==================
# Export for edgeR
# ==================
cat("\n=== Exporting for edgeR ===\n")
pdata$File <- file.path(OUT_DIR, "input_DE", paste0(pdata$Sample, ".txt"))
for (i in seq_len(nrow(pdata))) export_for_edgeR(pdata$Sample[i], pdata$File[i])

# ==================
# Differential analysis with edgeR
# ==================
cat("\n=== Running differential analysis ===\n")

y <- readBismark2DGE(pdata$File, sample.names = pdata$Sample)

# Filter by coverage
Methylation <- gl(2, 1, ncol(y), labels = c("Me", "Un"))  # Methylated / Unmethylated
Me <- y$counts[, Methylation == "Me"]
Un <- y$counts[, Methylation == "Un"]
Coverage <- Me + Un
HasCoverage <- rowSums(Coverage >= MIN_COVERAGE) >= MIN_SAMPLES
cat("Sites with coverage >=", MIN_COVERAGE, "in >=", MIN_SAMPLES, "samples:", sum(HasCoverage), "\n")
y <- y[HasCoverage, , keep.lib.sizes = FALSE]

# Normalize library sizes (same size for the Me and Un columns of a sample)
TotalLibSize <- 0.5 * y$samples$lib.size[Methylation == "Me"] +
  0.5 * y$samples$lib.size[Methylation == "Un"]
y$samples$lib.size <- rep(TotalLibSize, each = 2)

# Design matrix
designSL <- model.matrix(~0 + Group, data = pdata)
design <- modelMatrixMeth(designSL)
colnames(design) <- c(colnames(design)[1:nrow(pdata)], levels(pdata$Group))

cat("\n=== Estimating dispersion (this may take several minutes) ===\n")
y <- estimateDisp(y, design)

cat("\n=== Fitting model (this may take several minutes) ===\n")
fit <- glmQLFit(y, design)
save(fit, file = file.path(OUT_DIR, paste0("fit_", CTRL, "_", IVT, ".RData")))

# Contrast IVT - CTRL
contr_ivt <- makeContrasts(contrasts = paste(IVT, "-", CTRL), levels = design)
qlf_ivt <- glmQLFTest(fit, contrast = contr_ivt)

cat("\n=== Summary of differential sites ===\n")
print(summary(decideTests(qlf_ivt)))

top_ivt <- as.data.frame(topTags(qlf_ivt, n = Inf))
top_ivt <- top_ivt[, setdiff(colnames(top_ivt), c("Chr", "Locus")), drop = FALSE]

# ==================
# Merge results with annotations
# ==================
cat("\n=== Merging results with annotations ===\n")

samples_ctrl <- pdata$Sample[pdata$Group == CTRL]
samples_ivt <- pdata$Sample[pdata$Group == IVT]

# Annotation: union of all samples, one row per site
annot <- bind_rows(lapply(sites, function(df)
  df[, c("coord", "chr", "pos", "nt", "kind", "gene",
         "motif", "DRACH", "DRAC", "RAC", "RRAT", "strand")])) %>%
  distinct(coord, .keep_all = TRUE)

group_means <- function(sample_set, column) {
  m <- merge_sample_column(sample_set, column)
  data.frame(coord = rownames(m), value = rowMeans(as.matrix(m), na.rm = TRUE))
}

res <- top_ivt %>%
  rownames_to_column("coord") %>%
  left_join(annot, by = "coord") %>%
  left_join(rename(group_means(samples_ctrl, "protected"), mean_ctrl = value), by = "coord") %>%
  left_join(rename(group_means(samples_ivt, "protected"), mean_IVT = value), by = "coord") %>%
  left_join(rename(group_means(samples_ctrl, "nbA.G"), cov_ctrl = value), by = "coord") %>%
  left_join(rename(group_means(samples_ivt, "nbA.G"), cov_IVT = value), by = "coord") %>%
  mutate(delta_ivt = mean_ctrl - mean_IVT) %>%
  arrange(FDR)

# ==================
# Filter significant sites
# ==================
res_sig <- res %>% filter(delta_ivt >= DELTA_THRESHOLD & PValue <= FDR_THRESHOLD)

cat("\n=== Significant sites ===\n")
cat("Total:", nrow(res_sig), "\n")
cat("Hypermethylated (", CTRL, " > ", IVT, "):", sum(res_sig$delta_ivt > 0), "\n")

fwrite(res, file.path(OUT_DIR, "tables", paste0(CTRL, "_vs_", IVT, "_all.tsv")), sep = "\t")
fwrite(res_sig, file.path(OUT_DIR, "tables",
                          paste0(CTRL, "_vs_", IVT, "_delta", DELTA_THRESHOLD,
                                 "_p", FDR_THRESHOLD, ".tsv")), sep = "\t")

# ========================
# Export results bedRMod
# ========================
res_p <- res %>% filter(PValue <= FDR_THRESHOLD)
out_bed <- function(name) file.path(OUT_DIR, "bedRMod", name)

write_bedRMod(res, out_bed(paste0(CTRL, "_all_bedRMod.bed")),
              pvalue_col = "PValue", frequency_col = "mean_ctrl", coverage_col = "cov_ctrl")
write_bedRMod(res_p, out_bed(paste0(CTRL, "_pV", FDR_THRESHOLD, "_bedRMod.bed")),
              pvalue_col = "PValue", frequency_col = "mean_ctrl", coverage_col = "cov_ctrl")
write_bedRMod(res, out_bed(paste0(IVT, "_all_bedRMod.bed")),
              pvalue_col = "PValue", frequency_col = "mean_IVT", coverage_col = "cov_IVT")
write_bedRMod(res_p, out_bed(paste0(IVT, "_pV", FDR_THRESHOLD, "_bedRMod.bed")),
              pvalue_col = "PValue", frequency_col = "mean_IVT", coverage_col = "cov_IVT")

cat("\n=== Pipeline complete! ===\n")



