compute_signatures <- function(prepared_data, input) {
  prepared_dataframe <- prepared_data[[1]]
  # Split dataframe by ID and init hashtables
  patient_split_dataframe <- split(
    prepared_dataframe,
    prepared_dataframe$Patient.ID
  )
  signature_patients_map <- hash()
  signature_amounts_per_patient <- hash()
  # Handle signatures separately for every patient
  for (id in names(patient_split_dataframe)) {
    adj_matrix <- compute_adjacency_matrix(
      patient_split_dataframe[[id]],
      input$radio_multiple_variants_per_clone
    )
    if (any(adj_matrix != 0)) {
      variant_count <- as.integer(dim(adj_matrix)[1])
      pathways_matrix <- adj_matrix
      for (signature_length in seq.int(2, input$signature_length_setter)) {
        if (signature_length != 2) {
          pathways_matrix <- compute_n_walk_matrix(pathways_matrix, adj_matrix)
        }
        valid_paths <- nzwhich(pathways_matrix, arr.ind = TRUE)
        signature_amounts_per_patient[[id]][[signature_length]] <-
          nrow(valid_paths)
        for (i in seq_len(nrow(valid_paths))) {
          variant_values <-
            patient_split_dataframe[[id]]$Variant[valid_paths[i, ]]
          key <- paste(variant_values, collapse = ";")
          if (is.null(signature_patients_map[[key]])) {
            signature_patients_map[[key]] <- c(id)
          } else {
            if (!(id %in% signature_patients_map[[key]])) {
              signature_patients_map[[key]] <- append(
                signature_patients_map[[key]],
                id
              )
            }
          }
        }
      }
    }
  }
  signature_dataframe <- hash_table_to_dataframe(signature_patients_map)
  signature_dataframe <- signature_dataframe[
    signature_dataframe$`Number of occurrences` >=
      input$minimum_support_of_signatures,
  ]
  return(list(
    signature_dataframe,
    list(
      prepared_dataframe,
      signature_patients_map,
      signature_amounts_per_patient
    )
  ))
}

compute_alt_statistics <- function(
  signature_dataframe,
  prepared_dataframe,
  signature_patients_map,
  signature_amounts_per_patient,
  input
) {
  signature_dataframe$`p value` <- NA
  variant_cols <- grep(
    "^Variant \\d+$",
    colnames(signature_dataframe),
    value = TRUE
  )
  p_values <- hash()
  possible_occurences <- hash()
  # Create hashmap linking every variant to observed patient IDs.
  variant_observation <- hash()
  patient_split_dataframe <- split(
    prepared_dataframe,
    prepared_dataframe$Patient.ID
  )
  for (id in names(patient_split_dataframe)) {
    observed_variants <- patient_split_dataframe[[id]]$Variant
    for (observed_variant in observed_variants) {
      if (is.null(variant_observation[[observed_variant]])) {
        variant_observation[[observed_variant]] <- c(id)
      } else {
        variant_observation[[observed_variant]] <- c(
          variant_observation[[observed_variant]],
          id
        )
      }
    }
  }
  for (i in seq_len(nrow(signature_dataframe))) {
    variant_key <- gsub(
      "(;NA)+$",
      "",
      paste(signature_dataframe[i, variant_cols], collapse = ";")
    )
    variants <- unlist(strsplit(variant_key, ";"))
    # Get common patients that all variants appear in
    id_lists <- lapply(variants, function(v) variant_observation[[v]])
    id_lists <- Filter(Negate(is.null), id_lists)
    common_patients <- if (length(id_lists) > 0) {
      Reduce(intersect, id_lists)
    } else {
      c()
    }
    # Prepare significance computation
    value <- signature_patients_map[[variant_key]]
    trials <- c()
    probabilites <- c()
    possible_occurences[[variant_key]] <- 0
    for (id in common_patients) {
      observed_variants <- patient_split_dataframe[[id]]$Variant
      observed_counts <- table(observed_variants)
      required_counts <- table(variants)
      observed_counts <- observed_counts[names(required_counts)]
      observed_counts[is.na(observed_counts)] <- 0
      if (all(required_counts <= observed_counts)) {
        possible_occurences[[variant_key]] <-
          possible_occurences[[variant_key]] + 1
        if (input$significance_estimation_method == "Sandmann Nick") {
          multiplier <- 1
          for (variant in variants) {
            multiplier <- multiplier * observed_counts[[variant]]
          }
          probability <- multiplier /
            (length(unique(observed_variants))^length(variants))
          probabilites <- c(probability, probabilites)
        } else {
          probabilites <- c(
            1 / (length(unique(observed_variants))^length(variants)),
            probabilites
          )
        }
        trials <- c(
          signature_amounts_per_patient[[id]][[length(variants)]],
          trials
        )
      }
    }
    patient_probabilities <- 1 - (1 - probabilites)^trials
    p_values[[variant_key]] <- .signature_enrichment_tail(
      length(value), patient_probabilities
    )
  }
  return(list(p_values, possible_occurences))
}

hash_table_to_dataframe <- function(signature_patients_map) {
  keys_vec <- keys(signature_patients_map)

  if (length(keys_vec) == 0) {
    return(data.frame(
      `Variant 1` = character(0),
      `Number of occurrences` = numeric(0),
      Patients = character(0),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }

  max_parts <- max(vapply(keys_vec, function(key) {
    length(unlist(strsplit(key, ";")))
  }, integer(1)))
  rows <- list()

  for (key in keys_vec) {
    split_key <- unlist(strsplit(key, ";"))
    if (length(split_key) < max_parts) {
      split_key <- c(split_key, rep(NA, max_parts - length(split_key)))
    }

    support <- length(signature_patients_map[[key]])
    row <- c(
      split_key,
      support,
      paste0(
        'c("',
        paste(signature_patients_map[[key]], collapse = '", "'),
        '")'
      )
    )
    rows <- append(rows, list(row))
  }

  df <- do.call(rbind, rows)
  col_names <- c(
    paste0("Variant ", seq_len(max_parts)),
    "Number of occurrences",
    "Patients"
  )
  colnames(df) <- col_names
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  df$`Number of occurrences` <- as.numeric(df$`Number of occurrences`)

  return(df)
}

add_stats_to_dataframe <- function(df, stats) {
  p_values <- stats[[1]]
  possible_occurences <- stats[[2]]
  df <- as.data.frame(t(apply(df, 1, function(row) {
    key <- gsub(
      "(;NA)+$",
      "",
      paste(row[grep("^Variant", colnames(df), value = TRUE)], collapse = ";")
    )
    row[["Relative occurrence among patients with identical variants"]] <-
      as.numeric(row[["Number of occurrences"]]) /
      as.numeric(possible_occurences[[key]])
    row[["p value"]] <- p_values[[key]]
    return(row)
  })))
  return(df)
}

.signature_enrichment_tail <- function(observed, probabilities) {
  # Exact Bernoulli convolution gives P(X >= observed), including equality.
  mass <- 1
  for (probability in probabilities) {
    mass <- c(mass * (1 - probability), 0) + c(0, mass * probability)
  }
  sum(mass[seq_along(mass) - 1 >= observed])
}
