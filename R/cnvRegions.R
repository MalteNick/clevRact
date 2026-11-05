# Closed genomic intervals. Segments are alternative descriptions of one event;
# they are never inserted into a patient tree as additional events.
.cnv_exact_interval <- function(label) {
  compact <- gsub("\\s+", "", label)
  parts <- regmatches(compact, regexec(
    "^(CNV_(?:del|dup|loh|post_del|unknown))_([0-9]+|X|Y)_([0-9]+)_([0-9]+)$",
    compact, perl = TRUE))[[1]]
  if (length(parts) != 5L) return(NULL)
  start <- as.numeric(parts[4]); end <- as.numeric(parts[5])
  if (!is.finite(start) || !is.finite(end) || end < start) {
    stop("Invalid exact CNV interval: ", label)
  }
  data.frame(original = label, type = parts[2], chromosome = parts[3],
    start = start, end = end, stringsAsFactors = FALSE)
}

cnv_region_partition <- function(labels) {
  intervals <- Filter(Negate(is.null), lapply(unique(as.character(labels)), .cnv_exact_interval))
  if (!length(intervals)) stop("CNV region postprocessing requires exact CNV coordinate identifiers.")
  intervals <- do.call(rbind, intervals)
  regions <- list(); mapping <- list()
  for (chromosome in unique(intervals$chromosome)) {
    chr <- intervals[intervals$chromosome == chromosome, , drop = FALSE]
    # End+1 distinguishes overlapping endpoints and keeps adjacent CNVs separate.
    cuts <- sort(unique(c(chr$start, chr$end + 1)))
    for (i in seq_len(length(cuts) - 1L)) {
      start <- cuts[i]; end <- cuts[i + 1L] - 1
      covering <- chr[chr$start <= start & chr$end >= end, , drop = FALSE]
      if (!nrow(covering)) next
      for (type in unique(covering$type)) {
        subset <- covering[covering$type == type, , drop = FALSE]
        label <- paste(type, "region", chromosome,
          format(start, scientific = FALSE, trim = TRUE),
          format(end, scientific = FALSE, trim = TRUE), sep = "_")
        regions[[length(regions) + 1L]] <- data.frame(
          region = label, type = type, chromosome = chromosome,
          start = start, end = end, width = end - start + 1,
          source_intervals = nrow(subset), stringsAsFactors = FALSE)
        mapping[[length(mapping) + 1L]] <- data.frame(
          original = subset$original, region = label, stringsAsFactors = FALSE)
      }
    }
  }
  list(regions = do.call(rbind, regions), mapping = do.call(rbind, mapping),
    intervals = intervals, coordinate_convention = "closed inclusive")
}

# Accept the native signature result, including the full pre-support-filter hash.
# This permits individually rare exact CNVs to pool after signature computation.
postprocess_cnv_regions <- function(signature_information, minimum_support = 2L,
                                   max_expansions_per_signature = 100000L,
                                   region_partition = NULL) {
  stopifnot(length(minimum_support) == 1L, is.finite(minimum_support), minimum_support >= 1)
  signatures <- signature_information[[1]]
  side <- signature_information[[2]]
  if (length(side) >= 2L && inherits(side[[2]], "hash")) {
    signatures <- hash_table_to_dataframe(side[[2]])
  }
  data <- side[[1]]
  partition <- if (is.null(region_partition)) cnv_region_partition(data$Variant) else region_partition
  stopifnot(is.data.frame(partition$mapping), is.data.frame(partition$regions))
  alternatives <- split(partition$mapping$region, partition$mapping$original)
  cols <- grep("^Variant [0-9]+$", names(signatures), value = TRUE)
  records <- new.env(hash = TRUE, parent = emptyenv())
  provenance <- list()
  for (i in seq_len(nrow(signatures))) {
    labels <- as.character(unlist(signatures[i, cols, drop = FALSE]))
    labels <- labels[!is.na(labels) & labels != "NA" & nzchar(labels)]
    if (!any(labels %in% names(alternatives))) next
    # Native patient serialization is parsed, but only a literal c(...) call is
    # accepted. Never evaluate arbitrary imported R code.
    expression <- str2lang(as.character(signatures$Patients[i]))
    if (!is.call(expression) || !identical(expression[[1]], as.name("c"))) {
      stop("Expected a literal patient vector in computed signatures.")
    }
    values <- as.list(expression)[-1L]
    if (!all(vapply(values, is.character, logical(1)))) stop("Invalid patient vector.")
    patients <- unique(unlist(values, use.names = FALSE))
    choices <- lapply(labels, function(label) {
      if (label %in% names(alternatives)) alternatives[[label]] else label
    })
    combinations <- prod(lengths(choices))
    if (combinations > max_expansions_per_signature) {
      stop("CNV region expansion exceeds the limit for one signature; restrict its length or catalogue.")
    }
    grid <- expand.grid(choices, stringsAsFactors = FALSE)
    for (j in seq_len(nrow(grid))) {
      substituted <- as.character(unlist(grid[j, , drop = FALSE]))
      key <- paste(substituted, collapse = ";")
      previous <- records[[key]]
      if (is.null(previous)) previous <- list(labels = substituted, patients = character(), sources = integer())
      previous$patients <- union(previous$patients, patients)
      previous$sources <- union(previous$sources, i)
      records[[key]] <- previous
    }
  }
  rows <- list()
  for (key in sort(ls(records))) {
    record <- records[[key]]
    if (length(record$patients) < minimum_support) next
    row <- as.list(setNames(rep(NA_character_, length(cols)), cols))
    row[seq_along(record$labels)] <- as.list(record$labels)
    row[["Number of occurrences"]] <- length(record$patients)
    row$Patients <- paste(capture.output(dput(sort(record$patients))), collapse = "")
    row[["Source signatures"]] <- length(record$sources)
    rows[[length(rows) + 1L]] <- as.data.frame(row, check.names = FALSE, stringsAsFactors = FALSE)
    for (index in record$sources) provenance[[length(provenance) + 1L]] <- data.frame(
      region_signature = key, source_signature = paste(
        as.character(unlist(signatures[index, cols, drop = FALSE])), collapse = ";"),
      source_patients = as.character(signatures$Patients[index]), stringsAsFactors = FALSE)
  }
  result <- if (length(rows)) do.call(rbind, rows) else signatures[0, cols, drop = FALSE]
  # Source p-values and relative frequencies cannot be inherited after pooling.
  list(signatures = result, regions = partition$regions, mapping = partition$mapping,
    provenance = if (length(provenance)) do.call(rbind, provenance) else data.frame(),
    coordinate_convention = partition$coordinate_convention,
    scope = "Alternative CNV-region signatures computed after extraction; no new patient-tree events or inherited significance")
}
