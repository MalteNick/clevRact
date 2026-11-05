lookup <- function(chrom, start, end) {
  gene_data <- read.csv(
    .annotation_file("GRCh37.csv"),
    sep = ";",
    stringsAsFactors = FALSE
  )
  subset_rows <- gene_data[
    gene_data$Chromosome == chrom &
      gene_data$Begin <= end &
      gene_data$End >= start,
  ]
  return(unique(subset_rows$Symbol))
}

exact_deletions_to_genes <- function(df) {
  new_df <- df[0, ]
  blocks <- split(df, df$ID)
  queue <- blocks

  find_first_cnv_del_row <- function(variants) {
    pattern <- "^CNV_del_([0-9XY]+)_(\\d+)_(\\d+)$"
    matches <- regexec(pattern, variants)
    parts <- regmatches(variants, matches)

    for (i in seq_along(parts)) {
      if (length(parts[[i]]) == 4) {
        chrom <- parts[[i]][2]
        start <- as.integer(parts[[i]][3])
        end <- as.integer(parts[[i]][4])
        if (end > start) {
          return(list(index = i, chrom = chrom, start = start, end = end))
        }
      }
    }
    return(NULL)
  }

  while (length(queue) > 0) {
    block <- queue[[1]]
    queue <- queue[-1]
    result <- find_first_cnv_del_row(block$Variant)

    if (is.null(result)) {
      new_df <- rbind(new_df, block)
    } else {
      chars <- lookup(result$chrom, result$start, result$end)
      for (ch in chars) {
        modified_block <- block
        modified_block$Variant[result$index] <- ch
        queue <- c(queue, list(modified_block))
      }
    }
  }

  return(new_df)
}
