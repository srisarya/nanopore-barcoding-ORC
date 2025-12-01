#!/usr/bin/env Rscript

# Amplicon Summary Script
# Usage: Rscript amplicon_summary.R <base_dir> <amplicon_keyword> [gene_keyword]
# Example: Rscript amplicon_summary.R Lakes_day1/ COI
# Example: Rscript amplicon_summary.R Lakes_day1 rRNA 18S
# Example: Rscript amplicon_summary.R Lakes_day1 rRNA 28S

suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript amplicon_summary.R <base_dir> <amplicon_keyword> [gene_keyword]\n")
  cat("Example: Rscript amplicon_summary.R Lakes_day1/ COI\n")
  cat("Example: Rscript amplicon_summary.R Lakes_day1 rRNA 18S\n")
  cat("Example: Rscript amplicon_summary.R Lakes_day1 rRNA 28S\n")
  quit(status = 1)
}

base_dir <- args[1]
amplicon_keyword <- args[2]
gene_keyword <- if (length(args) >= 3) args[3] else NULL

# Remove trailing slash if present
base_dir <- sub("/$", "", base_dir)

# Check if base directory exists
if (!dir.exists(base_dir)) {
  cat("Error: Base directory", base_dir, "does not exist\n")
  quit(status = 1)
}

# Find subdirectory containing the amplicon keyword
all_subdirs <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)
matching_dirs <- all_subdirs[grepl(amplicon_keyword, all_subdirs, ignore.case = TRUE)]

if (length(matching_dirs) == 0) {
  cat("Error: No subdirectory found containing '", amplicon_keyword, "' in ", base_dir, "\n", sep = "")
  cat("Available subdirectories:\n")
  cat(paste("  -", all_subdirs), sep = "\n")
  quit(status = 1)
}

if (length(matching_dirs) > 1) {
  cat("Warning: Multiple subdirectories match '", amplicon_keyword, "':\n", sep = "")
  cat(paste("  -", matching_dirs), sep = "\n")
  cat("Using the first match: ", matching_dirs[1], "\n", sep = "")
}

amplicon_folder <- matching_dirs[1]

# Construct output prefix automatically
base_name <- basename(base_dir)
if (!is.null(gene_keyword)) {
  output_prefix <- paste0(base_name, "_", gene_keyword, "_summary")
} else {
  output_prefix <- paste0(base_name, "_", amplicon_keyword, "_summary")
}

# Set data directory
data_dir <- file.path(base_dir, amplicon_folder)

cat("Base directory:", base_dir, "\n")
cat("Amplicon folder:", amplicon_folder, "\n")
cat("Data directory:", data_dir, "\n")
if (!is.null(gene_keyword)) {
  cat("Gene filter:", gene_keyword, "\n")
}
cat("Output prefix:", output_prefix, "\n\n")

# Get all sample directories
sample_dirs <- list.dirs(data_dir, recursive = FALSE, full.names = TRUE)
sample_names <- basename(sample_dirs)

# Initialize results list
results_list <- list()

# Process each sample
for (i in seq_along(sample_dirs)) {
  sample_name <- sample_names[i]
  sample_dir <- sample_dirs[i]
  
  # Find all fasta/fa files in this directory
  all_fasta_files <- list.files(sample_dir, pattern = "\\.(fasta|fa)$", full.names = TRUE)
  
  # Filter by gene keyword if provided
  if (!is.null(gene_keyword)) {
    fasta_files <- all_fasta_files[grepl(gene_keyword, basename(all_fasta_files), ignore.case = TRUE)]
    
    # If no files match the gene keyword, check if there are any fasta files at all
    if (length(fasta_files) == 0 && length(all_fasta_files) > 0) {
      # Gene-specific file not found, but other fasta files exist
      results_list[[i]] <- data.frame(
        sample = sample_name,
        amplicon_found = "no",
        num_hits = 0L,
        best_hit_readcount = NA_integer_,
        best_hit_header = NA_character_,
        stringsAsFactors = FALSE
      )
      next
    }
  } else {
    fasta_files <- all_fasta_files
  }
  
  if (length(fasta_files) == 0) {
    # No amplicon found
    results_list[[i]] <- data.frame(
      sample = sample_name,
      amplicon_found = "no",
      num_hits = 0L,
      best_hit_readcount = NA_integer_,
      best_hit_header = NA_character_,
      stringsAsFactors = FALSE
    )
  } else {
    # Amplicon found - process all fasta files
    all_readcounts <- integer(0)
    all_headers <- character(0)
    total_sequences <- 0L
    
    for (fasta_file in fasta_files) {
      # Read sequences using Biostrings
      seqs <- tryCatch({
        readDNAStringSet(fasta_file)
      }, error = function(e) {
        # If DNA fails, try reading as any type
        readAAStringSet(fasta_file)
      })
      
      total_sequences <- total_sequences + length(seqs)
      
      # Extract readcounts from sequence names
      seq_names <- names(seqs)

      if (length(seq_names) > 0) {
        # Pattern to match _readcount_ followed by digits
        readcounts <- sapply(seq_names, function(name) {
          match_result <- regmatches(name, regexpr("_readcount_([0-9]+)", name, perl = TRUE))
          if (length(match_result) > 0) {
            num_str <- sub(".*_readcount_([0-9]+).*", "\\1", match_result)
            return(as.integer(num_str))
          } else {
            return(NA_integer_)
          }
        }, USE.NAMES = FALSE)
        
        # Collect valid readcounts and corresponding headers
        valid_indices <- which(!is.na(readcounts))
        if (length(valid_indices) > 0) {
          all_readcounts <- c(all_readcounts, readcounts[valid_indices])
          all_headers <- c(all_headers, seq_names[valid_indices])
        }
      }
    }
    
    # Find best hit readcount and corresponding header
    if (length(all_readcounts) > 0) {
      best_idx <- which.max(all_readcounts)
      best_readcount <- all_readcounts[best_idx]
      best_header <- all_headers[best_idx]
    } else {
      best_readcount <- NA_integer_
      best_header <- NA_character_
    }
    
    results_list[[i]] <- data.frame(
      sample = sample_name,
      amplicon_found = "yes",
      num_hits = total_sequences,
      best_hit_readcount = best_readcount,
      best_hit_header = best_header,
      stringsAsFactors = FALSE
    )
  }
}

# Combine results into data frame
results <- bind_rows(results_list)

# Sort by sample name
results <- results %>% arrange(sample)

# Write TSV output
tsv_file <- paste0(output_prefix, ".tsv")
write.table(results, file = tsv_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("TSV output written to:", tsv_file, "\n")

# Print summary statistics
cat("\nSummary Statistics:\n")
cat("Total samples:", nrow(results), "\n")
cat("Samples with amplicons:", sum(results$amplicon_found == "yes"), "\n")
cat("Samples without amplicons:", sum(results$amplicon_found == "no"), "\n")
cat("Total hits:", sum(results$num_hits), "\n")