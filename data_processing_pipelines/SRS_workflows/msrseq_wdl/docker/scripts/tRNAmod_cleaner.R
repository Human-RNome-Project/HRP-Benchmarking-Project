#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: tRNAmod_cleaner.R <tsv_dir> <out_csv> <out_zip>\n")
  quit(status = 1)
}

input_path <- args[1]
out_csv    <- args[2]
out_zip    <- args[3]

`%nin%` = Negate(`%in%`)

data_cleaner <- function(PATH) {
  message("Reading files from: ", PATH)
  FILES <- data.frame(file_name = sort(list.files(PATH, pattern = "\\.tsv$"))) %>%
    filter(!grepl("unassigned", file_name),
           file_name %nin% c("kallisto.out", "kallisto.err"))
  
  FILES <- FILES %>%
    separate(file_name, 
             sep="_", c("library", "junk1", "junk2", "junk3",
                        "barcode", "junk4",
                        "bin_start", "bin_stop",
                        "junkmore"),
             fill="right", remove=FALSE) %>%
    select(-matches("junk"))
  
  namelist_master <- data.frame(
    Index = c(rep("TP-AB-19s-R-Con", 3), rep("TP-AB-19s-R-BS", 3), rep("TP-AB-19s-R-CBH", 3)),
    Barcode = rep(c("bc8", "bc9", "bc10"), 3),
    Treatment = c(rep("HRPC_ctrl", 3), rep("HRPC_BS", 3), rep("HRPC_CBH", 3)),
    Rep = rep(1:3, 3)
  )
  
  FILES <- FILES %>%
    inner_join(namelist_master, 
               by = c("library" = "Index",
                      "barcode" = "Barcode")) %>%
    select(-library, -barcode) %>%
    rename(treatment = Treatment,
           rep = Rep)
  
  FILES$file_name <- as.character(FILES$file_name)
  FILES <- FILES %>%
    filter(rep != "remove")
  
  FILES <- FILES %>%
    mutate(bin_start = ifelse(is.na(bin_start), -3, bin_start),
           bin_stop = ifelse(is.na(bin_stop), -3, bin_stop))
  
  message("Found ", nrow(FILES), " matching files to process")
  
  read_in_one <- function(row) {
    fpath <- file.path(PATH, row$file_name)
    output <- read.csv(fpath, header=TRUE, sep="\t") %>%
      mutate(rep = row$rep,
             treatment = row$treatment,
             bin_start = row$bin_start,
             bin_stop  = row$bin_stop)
    
    output$gene <- as.character(output$gene)
    output$base <- as.character(output$base)
    return(output)
  }
  
  Ecoli_charging_counts_data <- FILES %>%
    group_by(file_name) %>%
    do(read_in_one(.)) %>%
    ungroup()
  
  final_data <- Ecoli_charging_counts_data %>%
    mutate(source = ifelse(grepl("mt", gene), "mitochondrial", "cytosolic")) %>%
    separate(gene, sep="-", c("junk1", "AA", "anticodon", "junk2", "junk3"), fill="right", remove=FALSE) %>%
    select(-junk1, -junk2, -junk3) %>%
    mutate(
      AA = ifelse(source == "mitochondrial", str_extract(gene, "(?<=mt)\\w{3}"), AA),
      anticodon = ifelse(source == "mitochondrial", str_sub(gene, -3), anticodon)
    )
  
  final_data <- final_data %>%
    mutate(gene = ifelse(is.na(anticodon), str_replace(gene, ".*_", ""), gene),
           AA = ifelse(is.na(anticodon), "Meso", AA),
           anticodon = ifelse(is.na(anticodon), gene, anticodon))
  
  final_data <- final_data %>%
    mutate(gene = str_replace(gene, ".*tRNA-", "")) %>%
    na.omit(select = c(AA, anticodon))
  
  return(final_data)
}

data_cleaned <- data_cleaner(input_path)

message("Writing CSV to: ", out_csv)
write.csv(data_cleaned, out_csv, row.names = TRUE)

message("Zipping to: ", out_zip)
out_csv_dir  <- dirname(normalizePath(out_csv))
out_csv_base <- basename(out_csv)
out_zip_full <- normalizePath(out_zip, mustWork = FALSE)

old_wd <- getwd()
setwd(out_csv_dir)
zip(out_zip_full, out_csv_base)
setwd(old_wd)

message("Done tRNA modification data aggregation.")
