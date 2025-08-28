#!/usr/bin/env Rscript


if (!requireNamespace("argparse", quietly = TRUE)) {
  install.packages("argparse", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("rtracklayer", quietly = TRUE)) {
  install.packages("rtracklayer", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table", repos = "https://cloud.r-project.org")
}

library(argparse)
library(rtracklayer)
library(data.table)


parser <- ArgumentParser(description = "Split UTR regions into 5' UTR and 3' UTR in a GTF file.")

parser$add_argument("-i", "--input", type = "character", required = TRUE, 
                    help = "Path to the input GTF file (required)")
parser$add_argument("-o", "--output", type = "character", required = TRUE,
                    help = "Path to the output GTF file (required)")
parser$add_argument("-v", "--version", action = "version", version = "1.0.0",
                    help = "Show script version and exit")


args <- parser$parse_args()

# split UTE in to 5UTR and 3UTR
split_utr_regions <- function(gtf_obj) {

  dt <- as.data.table(copy(gtf_obj))
  

  if (is.factor(dt$type)) {
    dt[, type := as.character(type)]
  }
  
  # get UTR idx
  utr_idx <- which(dt$type == "UTR")
  
  if (length(utr_idx) == 0) {
    warning("No UTR features found in the GTF object")
    return(dt)
  }
  
  # get CDS row
  cds_dt <- dt[type == "CDS", .(transcript_id, start, end, strand)]
  
  # calculate CDS location
  cds_ranges <- cds_dt[, .(
    cds_min = min(start),
    cds_max = max(end)
  ), by = .(transcript_id, strand)]
  
  # get UTR set
  utr_dt <- dt[utr_idx, ]
  
  # merge CDS with UTR
  utr_dt <- merge(utr_dt, cds_ranges, by = c("transcript_id", "strand"), all.x = TRUE)
  
  # process group by tx ids
  utr_dt[, utr_type := {

    if (unique(strand) == "+") {
      ifelse(end < cds_min, "5UTR", 
             ifelse(start > cds_max, "3UTR", "unknown"))
    } else {
      ifelse(start > cds_max, "5UTR", 
             ifelse(end < cds_min, "3UTR", "unknown"))
    }
  }, by = .(transcript_id, strand)]
  
  # process UTR overlap
  utr_dt[utr_type == "unknown", utr_type := ifelse(
    (strand == "+" & start < cds_min & end > cds_min) | 
      (strand == "-" & end > cds_max & start < cds_max),
    "overlap", "unknown"
  )]
  
  # Update UTR information in original dt
  dt[utr_idx, type := utr_dt$utr_type]
  
  # remove level
  if (is.factor(dt$type)) {
    dt[, type := as.character(type)]
  }
  out <- GenomicRanges::makeGRangesFromDataFrame(df = dt, keep.extra.columns = TRUE)
  return(out)
}

main <- function() {
  if (!file.exists(args$input)) {
    stop(paste("Input file does not exist:", args$input))
  }
  
  # import gtf
  cat("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] Reading GTF file...\n")
  gtf_obj <- rtracklayer::import(args$input)
  

  cat("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] Splitting UTR regions...\n")
  result_gr <- split_utr_regions(gtf_obj)

  cat("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] Writing output file...\n")
  rtracklayer::export(result_gr, args$output)
  
  cat("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] Processing complete! Output file:", args$output, "\n")
}


if (!interactive()) {
  main()
}