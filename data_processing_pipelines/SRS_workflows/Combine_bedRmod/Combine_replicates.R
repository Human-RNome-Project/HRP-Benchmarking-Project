suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(data.table)
})

# function to merge sites from different replciates
#' @param fins A list of input bedRmod files to be merged as replicates
#' @param fout Output file path, e.g., Merged_replicates.bed
#' @param min_required_reps Minimum number of replicates required to output a modification site as a confident one, default is 2
fun_merge_rep <- function(fins, fout, min_required_reps = 2) {
  # read bedRmod files of all replicates as a list of data.frames
  ldf = lapply(fins, function(f) {
    data.table::fread(f, data.table = F) %>%
      dplyr::rename(chrom = `#chrom`)
  })
  
  # combine replicates
  df = do.call(rbind, ldf) %>%
    group_by(chrom, chromStart, chromEnd, name, strand, thickStart, thickEnd) %>%
    summarise(
      score = round(max(score, na.rm = T), 2), 
      itemRgb = "0,0,0",
      coverage = sum(coverage, na.rm = T),
      frequency = round(mean(frequency, na.rm = T), 2),
      nRep = sum(nRep, na.rm = T),
      repName = paste(repName, collapse = ";"),
      repScore = paste(repScore, collapse = ";"),
      repCov = paste(repCov, collapse = ";"),
      repFreq = paste(repFreq, collapse = ";"),
      method = paste(method, collapse = ";")
    )
  
  # write out as bedRmod
  all_col = c("chrom", "chromStart", "chromEnd", "name", "score", "strand", 
              "thickStart", "thickEnd", "itemRgb", "coverage", "frequency",
              "nRep", "repName", "repScore", "repCov", "repFreq", "method")
  df = df[, all_col] %>% filter(nRep >= min_required_reps)
  colnames(df)[1] = paste0("#", colnames(df)[1])
  data.table::fwrite(df, fout, sep = "\t")
}


option_list = list(
  make_option(c("-i", "--input"), action = "store_true",
              help = "Input bedRmod files for each replicate [required]"),
  make_option(c("-o", "--output"), type = "character",
              help = "Combined bedRmod with pooled coverage and mean frequency across replicates [required]"),
  make_option(c("--min_reps"), type = "integer", default = 2,
              help = "Minimum number of replicates required to keep a modification site, default: 2")
)


parser = OptionParser(
  usage = "%prog [options]",
  description = "Combine bedRmod files for replicates.",
  option_list = option_list
)

args = parse_args(parser, positional_arguments = TRUE)

fun_merge_rep(
  fins = args$args,
  fout = args$options$output,
  min_required_reps = args$options$min_reps
)