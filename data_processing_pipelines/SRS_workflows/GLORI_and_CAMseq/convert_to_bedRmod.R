suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(data.table)
})

# function to convert output of GLORI-DUO-tools to bedRmod
#' @param fin Input file path, e.g., {sampleName}.totalm6A.FDR.csv.gz, sampleName will be inferred from the input file name
#' @param outdir Output file directory, the output bedRmod will be {outdir}/{sampleName}.bed
#' @param seqMethod Sequencing methods used which will be output in the "method" column, e.g., GLORI, CAM-seq
#' @param modification Modification type which will be output in the "name" column, default is m6A
#' @param min_AGcov minimal reads coverage for a site, default is 15
#' @param min_Acov minimal unconverted A coverage for a site, default is 5
#' @param min_Signal_Ratio minimum ratio of signal reads (eg. reads with unconverted As less than 3), default is 0.8
#' @param min_Ratio minimum m6A level, default is 0.1
#' @param max_FDR maximum FDR, default is 0.05
fun_convert_to_bedRmod = function(
  fin, outdir, seqMethod, modification = "m6A", 
  min_AGcov = 15, min_Acov = 5, min_Signal_Ratio = 0.8, 
  min_Ratio = 0.1, max_FDR = 0.05
  ) {
  # read data and generate columns of bedRmod
  sampleName = sub("\\.totalm6A\\.FDR\\.csv(\\.gz)?$", "", basename(fin))

  if (sampleName == basename(fin)) {
    stop(
      paste0(
        "Cannot infer sample name from input file: ", basename(fin),
        ". Expected a file ending with .totalm6A.FDR.csv or .totalm6A.FDR.csv.gz"
      ),
      call. = FALSE
    )
  }

  if (!dir.exists(args$outdir)) {
    dir.create(args$outdir, recursive = TRUE)
  }

  df = data.table::fread(fin, data.table = F) %>% 
    filter(
      AGcov >= min_AGcov, Acov >= min_Acov, Signal_Ratio >= min_Signal_Ratio, 
      Ratio >= min_Ratio, P_adjust < max_FDR
    ) %>%
    mutate(chrom = Chr, chromStart = Sites-1, chromEnd = Sites, name = modification, 
      score = round(-log10(P_adjust), 2), 
      strand = Strand,
      thickStart = Sites-1, thickEnd = Sites, itemRgb = "0,0,0",
      coverage = AGcov, frequency = round(Ratio * 100, 2), 
      nRep = 1, repScore = score, repCov = coverage, 
      repName = sampleName, repFreq = frequency, method = seqMethod)

  # write out as bedRmod
  all_col = c("chrom", "chromStart", "chromEnd", "name", "score", "strand", 
    "thickStart", "thickEnd", "itemRgb", "coverage", "frequency",
    "nRep", "repName", "repScore", "repCov", "repFreq", "method")

  
  df = df[, all_col]
  colnames(df)[1] = paste0("#", colnames(df)[1])
  data.table::fwrite(df, paste0(outdir, "/", sampleName, ".bed"), sep = "\t")
}

option_list = list(
  make_option(c("-i", "--input"), type = "character",
    help = "Input GLORI-DUO-tools CSV or CSV.GZ file [required]"),
  make_option(c("-o", "--outdir"), type = "character",
    help = "Output directory [required]"),
  make_option(c("--seqMethod"), type = "character", 
    help = "Sequencing method [required]"),
  make_option(c("-m", "--modification"), type = "character", default = "m6A",
    help = "Modification type [default: %default]"),
  make_option(c("--min_AGcov"),  type = "integer", default = 15, 
    help = "minimal reads coverage for a site [default: %default]"), 
  make_option(c("--min_Acov"), type = "integer", default = 5,
    help = "minimal unconverted A coverage for a site [default: %default]"), 
  make_option(c("--min_Signal_Ratio"), type = "double", default = 0.8, 
    help = "minimum ratio of signal reads (eg. reads with unconverted As less than 3) [default: %default]"), 
  make_option(c("--min_Ratio"), type = "double", default = 0.1, 
    help = "minimum m6A level [default: %default]"), 
  make_option(c("--max_FDR"), type = "double", default = 0.05, 
    help = " maximum FDR [default: %default]")
)

parser = OptionParser(
  usage = "%prog [options]",
  description = "Convert GLORI-DUO-tools output to bedRmod format.",
  option_list = option_list
)
args = parse_args(parser)

required_args = c("input", "outdir", "seqMethod")
missing_args = required_args[vapply(args[required_args], is.null, logical(1))]

if (length(missing_args) > 0) {
  print_help(parser)
  stop(
    paste0("Missing required argument(s): ", paste(missing_args, collapse = ", ")),
    call. = FALSE
  )
}

if (!file.exists(args$input)) {
  stop(paste0("Input file does not exist: ", args$input), call. = FALSE)
}

fun_convert_to_bedRmod(
  fin = args$input,
  outdir = args$outdir,
  seqMethod = args$seqMethod,
  modification = args$modification,
  min_AGcov = args$min_AGcov,
  min_Acov = args$min_Acov,
  min_Signal_Ratio = args$min_Signal_Ratio,
  min_Ratio = args$min_Ratio,
  max_FDR = args$max_FDR
)