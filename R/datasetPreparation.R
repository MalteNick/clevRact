prepareDataset <- function(df, input) {
    df_valid_trees <- df %>%
        group_by(.data$Patient.ID) %>%
        filter(all(.data$Parent %in% .data$Clone | .data$Parent == 0)) %>%
        ungroup()
    number_of_invalid_trees <-
        length(unique(df$Patient.ID)) -
        length(unique(df_valid_trees$Patient.ID))
    df <- df_valid_trees

    if (grepl("Yes", input$dropdown_variant_generation)) {
        preprocessing <- generate_variants_from_locations(
            df[, c(
                "Patient.ID", "Variant", "Clone", "Start", "End", "Parent",
                "Chromosome", "Variant.Type"
            )],
            input
        )
        if (is.null(preprocessing)) {
            return(NULL)
        }
        generated_dataset <- preprocessing[[1]]
        invalid_location_counter <- preprocessing[[2]]
        undescribed_region_counter <- preprocessing[[3]]
        number_of_unique_variants <- length(unique(generated_dataset$Variant))
    } else {
        generated_dataset <- df[, c(
            "Patient.ID", "Variant", "Clone", "Parent", "Chromosome",
            "Variant.Type"
        )]
        invalid_location_counter <- 0
        undescribed_region_counter <- 0
        number_of_unique_variants <- nrow(unique(
            df[, c("Chromosome", "Start", "End")]
        ))
    }

    single_occurence_counter <- 0
    cleaned_dataframe <- generated_dataset
    variant_id_lookup <- setNames(
        seq_len(nrow(generated_dataset)),
        generated_dataset$Variant
    )

    if (input$dropdown_cnv_clone_summary) {
        generated_dataset <- summarize_cnvs(generated_dataset)
    }

    if (
        grepl("Enumeration", input$dropdown_handler_reoccurrence, fixed = TRUE)
    ) {
        if (grepl("SNV", input$dropdown_handler_reoccurrence, fixed = TRUE)) {
            generated_dataset <- enumerate_dataframe(
                generated_dataset,
                "SNV/INDEL"
            )
        } else if (
            input$dropdown_handler_reoccurrence == "Enumeration for CNVs only"
        ) {
            generated_dataset <- enumerate_dataframe(generated_dataset, "CNV")
        } else {
            generated_dataset <- enumerate_dataframe(generated_dataset, "CNV")
            generated_dataset <- enumerate_dataframe(
                generated_dataset,
                "SNV/INDEL"
            )
        }
    }

    generated_dataset <- generated_dataset[
        order(
            generated_dataset$Patient.ID,
            generated_dataset$Clone,
            generated_dataset$Variant
        ),
    ]

    return(list(
        prepared_data = generated_dataset,
        cleaned_data = cleaned_dataframe,
        variant_id_lookup = variant_id_lookup,
        il_counter = invalid_location_counter,
        ur_counter = undescribed_region_counter,
        so_counter = single_occurence_counter,
        nou_variants = number_of_unique_variants,
        noi_trees = number_of_invalid_trees
    ))
}

generate_variants_from_locations <- function(dataset, input) {
    invalid_location_counter <- 0
    undescribed_region_counter <- 0
    generated_df <- data.frame(
        Patient.ID = character(),
        Variant = character(),
        Chromosome = character(),
        Clone = character(),
        Parent = character(),
        Variant.Type = character(),
        stringsAsFactors = FALSE
    )

    reference_genome <- if (
        input$dropdown_variant_generation == "Yes, using GRCh37"
    ) {
        "GRCh37"
    } else {
        "GRCh38"
    }

    gene_data_set <- fread(
        .annotation_file(paste0(reference_genome, ".csv"))
    )
    p_q_arm_positions <- fread(
        .annotation_file("arm_positions.csv")
    )

    for (i in seq_len(nrow(dataset))) {
        row <- dataset[i, ]
        chromosome <- trimws(row$Chromosome)
        start_pointer <- as.integer(trimws(row$Start))
        end_pointer <- as.integer(trimws(row$End))
        if (is.na(chromosome) || is.na(start_pointer) || is.na(end_pointer)) {
            next()
        }
        valid_chromosome <- grepl("^(1[0-9]|2[0-2]|[1-9]|X|Y)$", chromosome)
        valid_start <- !is.null(row$Start) && grepl("^-?\\d+$", row$Start)
        valid_end <- !is.null(row$End) && grepl("^-?\\d+$", row$End)
        if (
            !valid_chromosome ||
                !valid_start ||
                !valid_end ||
                start_pointer > end_pointer
        ) {
            invalid_location_counter <- invalid_location_counter + 1
            next
        }

        if (start_pointer == end_pointer) {
            if (input$dropdown_summary_snvs == "Use exact SNVs") {
                new_variant_identifier <- paste(
                    as.character(c(row$Chromosome, row$Start)),
                    collapse = "_"
                )
            } else {
                gene_matches <- gene_data_set[
                    tolower(trimws(gene_data_set$Chromosome)) == tolower(
                        chromosome
                    ) &
                        start_pointer >= as.integer(gene_data_set$Begin) &
                        start_pointer <= as.integer(gene_data_set$End),
                ]
                new_variant_identifier <- if (nrow(gene_matches) > 0) {
                    gene_matches$Symbol[1]
                } else {
                    undescribed_region_counter <- undescribed_region_counter + 1
                    "INVALID"
                }
            }
        } else {
            if (input$dropdown_summary_cnvs == "Use exact CNVs") {
                new_variant_identifier <- paste(
                    row$Variant.Type,
                    "_",
                    paste(
                        as.character(c(row$Chromosome, row$Start, row$End)),
                        collapse = "_"
                    )
                )
            } else if (
                input$dropdown_summary_cnvs == "Summarize on chromosome level"
            ) {
                new_variant_identifier <- gsub(
                    "\\s+",
                    "",
                    paste(row$Variant.Type, "_", row$Chromosome)
                )
            } else {
                begin_q <- as.integer(p_q_arm_positions[
                    p_q_arm_positions$Reference == reference_genome &
                        trimws(p_q_arm_positions$Chromosome) == chromosome
                ]$BeginQ)
                end_q <- as.integer(p_q_arm_positions[
                    p_q_arm_positions$Reference == reference_genome &
                        trimws(p_q_arm_positions$Chromosome) == chromosome
                ]$EndQ)
                new_variant_identifier <- if (row$Start == 0 &&
                    row$End == end_q) {
                    paste(row$Variant.Type, "_", row$Chromosome)
                } else if (row$Start == 0 && row$End == begin_q - 2) {
                    paste(row$Variant.Type, "_", row$Chromosome, "p", sep = "")
                } else if (row$Start == begin_q && row$End == end_q) {
                    paste(row$Variant.Type, "_", row$Chromosome, "q", sep = "")
                } else if (row$Start > 0 && row$End < begin_q) {
                    paste(row$Variant.Type, "_in_", row$Chromosome, "p",
                        sep = ""
                    )
                } else if (row$Start > begin_q) {
                    paste(row$Variant.Type, "_in_", row$Chromosome, "q",
                        sep = ""
                    )
                } else {
                    paste(row$Variant.Type, "_in_", row$Chromosome)
                }
            }
        }

        generated_df <- rbind(generated_df, data.frame(
            Patient.ID = row$Patient.ID,
            Variant = new_variant_identifier,
            Chromosome = chromosome,
            Clone = row$Clone,
            Parent = row$Parent,
            Variant.Type = row$Variant.Type
        ))
    }

    if (nrow(generated_df) == 0) {
        return(NULL)
    }
    rm(gene_data_set)

    generated_df_split <- split(generated_df, generated_df$Patient.ID)
    generated_df_split <- lapply(generated_df_split, function(df_part) {
        if (!any(df_part$Variant == "INVALID")) {
            return(df_part)
        }

        clones_to_remove <- df_part %>%
            group_by(.data$Clone) %>%
            summarize(all_invalid = all(.data$Variant == "INVALID")) %>%
            filter(.data$all_invalid) %>%
            pull("Clone")

        if (length(clones_to_remove) == 0) {
            return(df_part)
        }

        for (clone_to_remove in clones_to_remove) {
            parent_value <- unique(df_part$Parent[
                df_part$Clone == clone_to_remove
            ])
            if (length(parent_value) != 1) next
            df_part$Parent[df_part$Parent == clone_to_remove] <- parent_value
        }
        return(df_part)
    })

    final_df <- do.call(rbind, generated_df_split)
    final_df <- final_df %>% filter(.data$Variant != "INVALID")

    generated_df_split <- split(final_df, final_df$Patient.ID)
    generated_df_split <- lapply(generated_df_split, function(df_part) {
        unique_clones <- sort(unique(df_part$Clone))
        clone_mapping <- setNames(seq_along(unique_clones), unique_clones)
        df_part <- df_part %>%
            mutate(
                Clone = clone_mapping[as.character(.data$Clone)],
                Parent = ifelse(
                    as.character(.data$Parent) %in% names(clone_mapping),
                    clone_mapping[as.character(.data$Parent)],
                    .data$Parent
                )
            )
        return(df_part)
    })

    final_df <- do.call(rbind, generated_df_split)

    return(list(
        generated_dataframe = final_df,
        il_counter = invalid_location_counter,
        ur_counter = undescribed_region_counter
    ))
}

summarize_cnvs <- function(df) {
    required_columns <- c(
        "Patient.ID", "Chromosome", "Variant.Type", "Clone", "Variant", "Parent"
    )
    if (!all(required_columns %in% colnames(df))) {
        stop(
            "The dataframe must contain the columns: Patient.ID, Chromosome, ",
            "Variant.Type, Clone, Parent, and Variant"
        )
    }

    split_dfs <- split(df, df$Patient.ID)

    process_sub_df <- function(sub_df) {
        df_counts <- sub_df %>%
            group_by(.data$Chromosome, .data$Clone, .data$Variant.Type) %>%
            mutate(count = n())

        df_grouped <- df_counts %>%
            filter(.data$count > 1, grepl("CNV", .data$Variant.Type)) %>%
            group_by(.data$Chromosome, .data$Clone, .data$Variant.Type) %>%
            summarize(
                Patient.ID = .data$Patient.ID[1],
                Parent = .data$Parent[1],
                Variant = paste(.data$Variant[1], "_complex", sep = ""),
                .groups = "drop"
            )

        df_unique <- df_counts %>%
            filter(.data$count == 1) %>%
            select(-all_of("count"))
        df_non_unique_SNVs <- df_counts %>%
            filter(.data$count > 1 & !grepl("CNV", .data$Variant.Type)) %>%
            select(-all_of("count"))

        df_final <- bind_rows(df_grouped, df_unique, df_non_unique_SNVs)

        return(df_final)
    }

    processed_dfs <- lapply(split_dfs, process_sub_df)
    result_df <- do.call(rbind, processed_dfs)
    reordered_df <- result_df %>%
        select(all_of(c(
            "Patient.ID", "Variant", "Chromosome", "Variant.Type", "Clone",
            "Parent"
        )))

    return(reordered_df)
}

enumerate_dataframe <- function(data_frame, variant_type) {
    data_frame$Variant_alt <- data_frame$Variant
    patient_list <- split(data_frame, data_frame$Patient.ID)

    patient_list <- lapply(patient_list, function(patient_data_frame) {
        variant_counts <- table(patient_data_frame$Variant)
        variants_to_enumerate <- names(variant_counts[variant_counts > 1])

        for (variant in variants_to_enumerate) {
            for (i in seq_len(nrow(patient_data_frame))) {
                if (
                    patient_data_frame$Variant[i] == variant &&
                        grepl(variant_type, patient_data_frame$Variant.Type[i])
                ) {
                    parent <- patient_data_frame$Parent[i]
                    clone <- patient_data_frame$Clone[i]

                    while (parent != 0) {
                        filtered_df <- patient_data_frame[
                            patient_data_frame$Parent == parent &
                                patient_data_frame$Clone == clone &
                                patient_data_frame$Variant == variant,
                        ]
                        enumerated_values <- filtered_df$Variant_alt[
                            grepl(
                                paste0("^", variant, "_[a-z]$"),
                                filtered_df$Variant_alt
                            )
                        ]
                        enumerations <- sub(
                            paste0("^", variant, "_([a-z])$"),
                            "\\1",
                            enumerated_values
                        )

                        if (length(enumerations) > 0) {
                            highest_letter <- max(enumerations)
                            next_letter <- rawToChar(
                                as.raw(as.integer(charToRaw(
                                    highest_letter
                                )) + 1)
                            )
                            patient_data_frame$Variant_alt[i] <- paste0(
                                patient_data_frame$Variant[i],
                                "_",
                                next_letter
                            )
                            break
                        } else {
                            clone <- parent
                            parent <- patient_data_frame %>%
                                filter(.data$Clone == parent) %>%
                                select(all_of("Parent")) %>%
                                distinct() %>%
                                pull("Parent")
                        }
                    }

                    if (parent == 0) {
                        filtered_df <- patient_data_frame[
                            patient_data_frame$Parent == parent &
                                patient_data_frame$Clone == clone &
                                patient_data_frame$Variant == variant,
                        ]
                        enumerated_values <- filtered_df$Variant_alt[
                            grepl(
                                paste0("^", variant, "_[a-z]$"),
                                filtered_df$Variant_alt
                            )
                        ]
                        enumerations <- sub(
                            paste0("^", variant, "_([a-z])$"),
                            "\\1",
                            enumerated_values
                        )
                        if (length(enumerations) > 0) {
                            highest_letter <- max(enumerations)
                            next_letter <- rawToChar(
                                as.raw(as.integer(charToRaw(
                                    highest_letter
                                )) + 1)
                            )
                            patient_data_frame$Variant_alt[i] <- paste0(
                                patient_data_frame$Variant[i],
                                "_",
                                next_letter
                            )
                        } else {
                            patient_data_frame$Variant_alt[i] <- paste0(
                                patient_data_frame$Variant[i],
                                "_a"
                            )
                        }
                    }
                }
            }
        }
        return(patient_data_frame)
    })

    df_with_alterations <- do.call(rbind, patient_list)
    df_with_alterations$Variant <- df_with_alterations$Variant_alt
    enumerated_df <- df_with_alterations[
        setdiff(names(df_with_alterations), "Variant_alt")
    ]
    enumerated_df <- enumerated_df %>%
        mutate(
            Variant = if_else(
                str_detect(.data$Variant.Type, variant_type) &
                    !str_detect(.data$Variant, "_[a-z]$"),
                paste0(.data$Variant, "_a"),
                .data$Variant
            )
        )

    return(enumerated_df)
}
