suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(data.table)
})


# function to merge sites from different methods detecting the same modification, e.g. GLORI, CAM-seq, eTAM-seq
#' @param fins A list of input bedRmod files to be merged as different methods for the same modification
#' @param fout Output file path, e.g., Merged_methods.bed
fun_merge_method <- function(fins, fout) {
  ldf = lapply(fins, function(f) {
    data.table::fread(f, data.table = F) %>%
      dplyr::rename(chrom = `#chrom`)
  })
  
  ldf = lapply(ldf, function(x) {
    x %>% 
      arrange(score, coverage) %>%
      mutate(
        newscore = (1:n())/n()
      )
  })
  
  nmeth = length(ldf)
  df = do.call(rbind, ldf) %>%
    group_by(chrom, chromStart, chromEnd, name, strand, thickStart, thickEnd) %>%
    summarise(
      score = round(sum(newscore, na.rm = T)/nmeth*10000, 2),
      itemRgb = "0,0,0",
      coverage = sum(coverage, na.rm = T),
      frequency = round(mean(frequency, na.rm = T), 2),
      nRep = sum(nRep, na.rm = T),
      repName = paste(repName, collapse = ";"),
      repScore = paste(repScore, collapse = ";"),
      repCov = paste(repCov, collapse = ";"),
      repFreq = paste(repFreq, collapse = ";"),
      repNewScore = paste(newscore, collapse = ";"),
      method = paste(unique(method), collapse = ";")
    )
  
  # write out as bedRmod
  all_col = c("chrom", "chromStart", "chromEnd", "name", "score", "strand", 
              "thickStart", "thickEnd", "itemRgb", "coverage", "frequency",
              "nRep", "repName", "repScore", "repCov", "repFreq", "method")
  df = df[, all_col]
  colnames(df)[1] = paste0("#", colnames(df)[1])
  data.table::fwrite(df, fout, sep = "\t")
}


option_list = list(
  make_option(c("-i", "--input"), action = "store_true",
              help = "Input bedRmod files for different methods that detect the same modification type [required]"),
  make_option(c("-o", "--output"), type = "character",
              help = "Combined bedRmod with pooled coverage and mean frequency across replicates [required]")
)


parser = OptionParser(
  usage = "%prog [options]",
  description = "Combine bedRmod files for different methods that detect the same modification type.",
  option_list = option_list
)

args = parse_args(parser, positional_arguments = TRUE)

fun_merge_method(
  fins = args$args,
  fout = args$options$output
)