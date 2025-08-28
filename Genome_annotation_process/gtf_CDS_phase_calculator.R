#!/usr/bin/env Rscript

# load packages
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

# parameters
parser <- ArgumentParser(description = "Calculate phase scores for CDS features in a GTF file.")

parser$add_argument("-i", "--input", type = "character", required = TRUE, 
                    help = "Path to the input GTF file (required, can be gziped)")
parser$add_argument("-o", "--output", type = "character", required = TRUE,
                    help = "Path to the output GTF file (required)")
parser$add_argument("-v", "--version", action = "version", version = "1.0.0",
                    help = "Show script version and exit")


args <- parser$parse_args()

# CDS phase calculater
calculate_cds_phase <- function(gtf_obj) {
  # use data.table 
  dt <- data.table::as.data.table(copy(gtf_obj))
  
  # only change CDS
  cds_idx <- which(dt$type == "CDS")
  
  if (length(cds_idx) == 0) {
    warning("No CDS features found in the GTF object")
    return(dt)
  }
  
  # create CDS table
  cds_dt <- dt[cds_idx, ]
  
  # get CDS length
  cds_dt[, cds_length := end - start + 1]
  
  # group by transcripts
  cds_dt[, order_idx := .I] 
  
  # sort with transcript id for different strands
  cds_dt[strand == "+", rank := frank(start), by = transcript_id]
  cds_dt[strand == "-", rank := frank(-start), by = transcript_id]
  
  # calculate the phase score with order
  setorder(cds_dt, transcript_id, rank)
  
  # calculate the cumsum and phase score
  cds_dt[, cum_len := shift(cumsum(cds_length), fill = 0), by = transcript_id]
  
  cds_dt[, phase := (3 - (cum_len %% 3)) %% 3]
  
  # recover the original order of CDS table
  setorder(cds_dt, order_idx)
  
  # add phase score to the original table
  dt[cds_idx, phase := cds_dt$phase]
  out <- GenomicRanges::makeGRangesFromDataFrame(df = dt, keep.extra.columns = TRUE)
  
  return(out)
}

# main function
main <- function() {
  
  if (!file.exists(args$input)) {
    stop(paste("Input file does not exist:", args$input))
  }
  
  
  cat("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] Reading GTF file...\n")
  gtf_obj <- rtracklayer::import(args$input)
  
  
  cat("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] Calculating CDS phases...\n")
  result_gr <- calculate_cds_phase(gtf_obj)
  
  cat("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] Writing output file...\n")
  rtracklayer::export(result_gr, file = args$output)
  
  
  cat("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] Processing complete! Output file:", args$output, "\n")
}


if (!interactive()) {
  main()
}
