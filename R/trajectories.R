compute_trajectories <- function(
    signature_data,
    patient_data,
    input,
    signature_type
) {
    signatures <- expand_variants(signature_data)
    signatures <- signatures[
        ,
        !(names(signatures) %in% c("p.value", "Number.of.occurrences"))
    ]
    if (!nrow(signatures)) {
        signatures <- data.frame(
            ID = integer(), Clone = integer(),
            Parent = integer(), Variant = character(), Patients = character()
        )
    }
    # Algorithm 3: recurrent independent branches joined at the empty root.
    # These bases are discovered from the cohort, even when no linear seeds
    # exist.
    # Unset inputs preserve Algorithm 3 for existing programmatic callers.
    root_bases <- if (is.null(input$compute_independent_root_branches) ||
        isTRUE(as.logical(input$compute_independent_root_branches))) {
        process_indi_patient_branches(
            patient_data[[1]], signatures, signature_type,
            minimum_support = input$minimum_support_of_trajectories,
            allow_repeated = isTRUE(as.logical(
                input$multiple_identical_variants_for_branches
            ))
        )
    } else {
        list()
    }
    trajectories <- if (length(root_bases)) {
        rbind(signatures, do.call(rbind, root_bases))
    } else {
        signatures
    }
    trajectories_bases <- trajectories
    split_signatures <- split(signatures, signatures$ID)
    global_id_counter <- if (nrow(trajectories)) max(trajectories$ID) + 1 else 1
    for (x in seq_len(input$recombination_number_setter)) {
        if (is.null(trajectories_bases)) {
            break
        }
        # Generate pairing of previous trajectories and signatures
        split_trajectories <- split(trajectories_bases, trajectories_bases$ID)
        id_pairs <- expand.grid(
            traj_id = names(split_trajectories),
            sig_id = names(split_signatures),
            stringsAsFactors = FALSE
        )
        # Generate Combinations if valid
        results <- lapply(seq_len(nrow(id_pairs)), function(i) {
            trajectory_df <- split_trajectories[[id_pairs$traj_id[i]]]
            signature_df <- split_signatures[[id_pairs$sig_id[i]]]
            combine_trajectories(
                trajectory_df, signature_df, patient_data,
                input
            )
        })
        results <- Filter(Negate(is.null), results)
        if (length(results) == 0) break
        # Renumber IDs
        for (i in seq_along(results)) {
            df <- results[[i]]
            n_ids <- length(unique(df$ID))
            old_ids <- sort(unique(df$ID))
            new_ids <- global_id_counter + seq_len(n_ids)
            id_map <- setNames(new_ids, old_ids)
            df$ID <- id_map[as.character(df$ID)]
            results[[i]] <- df
            global_id_counter <- max(df$ID)
        }
        new_trajectories <- do.call(rbind, results)

        # Check if new trajectories actually occur in candidates
        split_trajectories <- split(new_trajectories, new_trajectories$ID)
        filtered <- lapply(split_trajectories, function(df_part) {
            patient_list <- eval(parse(text = df_part$Patients[1]))
            keep_patients <- c()
            for (patient in patient_list) {
                trajectory_variant_counts <- table(df_part$Variant)
                current_patient_data <- patient_data[[1]][
                    patient_data[[1]]$Patient.ID == patient, ,
                    drop = FALSE
                ]
                patient_variant_counts <- table(current_patient_data$Variant)
                patient_variant_counts_aligned <- patient_variant_counts[
                    names(trajectory_variant_counts)
                ]
                patient_variant_counts_aligned[
                    is.na(patient_variant_counts_aligned)
                ] <- 0
                if (all(
                    trajectory_variant_counts <= patient_variant_counts_aligned
                )) {
                    if (
                        is_trajectory_in_patient(
                            current_patient_data,
                            df_part,
                            signature_type
                        )
                    ) {
                        keep_patients <- c(keep_patients, patient)
                    }
                }
            }
            if (length(keep_patients) < input$minimum_support_of_trajectories) {
                return(NULL)
            }
            df_part$Patients <- paste0(
                "c(",
                paste(shQuote(keep_patients), collapse = ", "),
                ")"
            )
            df_part
        })
        filtered_non_null <- Filter(Negate(is.null), filtered)
        trajectories_bases <- if (length(filtered_non_null) > 0) {
            do.call(rbind, filtered_non_null)
        } else {
            NULL
        }
        trajectories <- rbind(trajectories_bases, trajectories)
    }

    # Remove doubles with different clone numbering
    split_trajectories <- split(trajectories, trajectories$ID)
    ids <- names(split_trajectories)
    ids_to_remove <- c()
    if (length(ids) > 1) {
        for (i in seq_len(length(ids) - 1)) {
            for (j in seq.int(i + 1, length(ids))) {
                id1 <- ids[i]
                id2 <- ids[j]
                if (id1 %in% ids_to_remove || id2 %in% ids_to_remove) next
                df1 <- split_trajectories[[id1]]
                df2 <- split_trajectories[[id2]]
                if (nrow(df1) != nrow(df2)) next
                if (!identical(df1$Patients, df2$Patients)) next
                if (compare_trajectories(df1, df2)) {
                    ids_to_remove <- c(ids_to_remove, id2)
                }
            }
        }
    }
    trajectories <- trajectories[!(trajectories$ID %in% ids_to_remove), ]
    unique_ids <- unique(trajectories$ID)
    new_ids <- setNames(seq_along(unique_ids), unique_ids)
    trajectories$ID <- new_ids[as.character(trajectories$ID)]
    trajectories$`Number of occurrences` <- vapply(
        trajectories$Patients,
        function(p) {
            length(eval(parse(text = p)))
        },
        integer(1)
    )
    return(trajectories)
}

estimate_significance <- function(
    trajectories_df,
    patient_df,
    input,
    signature_type
) {
    shapes <- assign_shapes(trajectories_df)
    patient_topology_statistics <- generate_shape_statistics(
        patient_df,
        shapes$shape_lookup,
        max_value <- max(trajectories_df$Clone, na.rm = TRUE),
        signature_type
    )
    trajectories <- shapes$trajectories
    #   trajectories <- exact_deletions_to_genes(trajectories)
    # }
    trajectories <- add_p_values(
        trajectories,
        patient_topology_statistics,
        patient_df
    )
    if (input$p_value_adjustment_trajectories == "Benjamini-Hochberg") {
        trajectories$`p value` <- p.adjust(trajectories$`p value`,
            method =
                "BH"
        )
    }
    trajectories$Shape <- NULL
    return(trajectories)
}

filter_duplicates <- function(trajectories_df) {
    return(trajectories_df)
}

expand_variants <- function(df) {
    variant_cols <- grep("^Variant", names(df), value = TRUE)
    id_counter <- 1
    result_list <- list()
    for (i in seq_len(nrow(df))) {
        row_values <- as.vector(unlist(df[i, variant_cols], use.names = FALSE))
        first_na <- which(is.na(row_values))
        variant_values <- if (length(first_na)) {
            row_values[seq_len(first_na[1] - 1)]
        } else {
            row_values
        }
        variant_values <- variant_values[variant_values != ""]
        n <- length(variant_values)
        if (n == 0) next
        rows <- data.frame(
            ID = rep(id_counter, n),
            Clone = seq_len(n),
            Parent = seq(0, n - 1),
            Variant = variant_values,
            `Number of occurrences` = rep(df$`Number of occurrences`[i], n),
            `p value` = rep(df$`p value`[i], n),
            `Patients` = rep(df$`Patients`[i], n),
            stringsAsFactors = FALSE
        )
        result_list[[length(result_list) + 1]] <- rows
        id_counter <- id_counter + 1
    }
    if (length(result_list) == 0) {
        return(data.frame())
    }
    result_df <- do.call(rbind, result_list)
    rownames(result_df) <- NULL
    return(result_df)
}

combine_trajectories <- function(
    trajectory_df,
    signature_df,
    patient_data,
    input
) {
    if (all(apply(
        trajectory_df[, !(names(trajectory_df) %in% "ID")],
        1,
        function(row) {
            any(apply(
                signature_df[, !(names(signature_df) %in% "ID")],
                1,
                function(x) identical(row, x)
            ))
        }
    ))) {
        return(NULL)
    }
    common_patients <- intersect(
        eval(parse(text = trajectory_df$Patients[1])),
        eval(parse(text = signature_df$Patients[1]))
    )
    if (length(common_patients) < input$minimum_support_of_trajectories) {
        return(NULL)
    } else {
        overlap <- intersect(trajectory_df$Variant, signature_df$Variant)
        if (length(overlap) == 0) {
            return(NULL)
        } else {
            root_variant <- signature_df$Variant[signature_df$Parent == 0]
            if (root_variant %in% overlap) {
                ID <- 1
                clones_with_root <- trajectory_df$Clone[
                    trajectory_df$Variant == root_variant
                ]
                new_trajectories_df <- data.frame(
                    Clone = integer(0),
                    ID = integer(0),
                    Variant = character(0),
                    Parent = integer(0),
                    Patients = character(0)
                )
                sig_to_add <- signature_df[-1, ]
                clone_offset <- max(trajectory_df$Clone)
                sig_to_add$Clone <- sig_to_add$Clone + clone_offset - 1
                if (input$multiple_identical_variants_for_branches) {
                    for (clone in clones_with_root) {
                        sig_to_add$Parent <- seq(
                            from = clone,
                            length.out = nrow(sig_to_add)
                        )
                        new_trajectory <- trajectory_df
                        new_trajectory <- rbind(new_trajectory, sig_to_add)
                        new_trajectory$Patients <- paste0(
                            'c("',
                            paste(common_patients, collapse = '", "'),
                            '")'
                        )
                        new_trajectory$ID <- ID
                        ID <- ID + 1
                        new_trajectories_df <- rbind(
                            new_trajectories_df,
                            new_trajectory
                        )
                    }
                } else {
                    for (clone in clones_with_root) {
                        clone_id <- clone
                        new_trajectory <- trajectory_df
                        signature_anchor <- 1
                        while (signature_anchor < max(signature_df$Clone)) {
                            children <- new_trajectory$Parent == clone_id
                            if (any(
                                rows <- new_trajectory$Variant[children] ==
                                    signature_df$Variant[
                                        signature_df$Parent == signature_anchor
                                    ][1]
                            )) {
                                clone_id <- new_trajectory$Clone[which(rows)[1]]
                                signature_anchor <- signature_anchor + 1
                            } else {
                                if (
                                    nrow(signature_df[
                                        signature_df$Clone > signature_anchor,
                                    ]) ==
                                        0
                                ) {
                                    break()
                                } else {
                                    new_nodes <- signature_df[
                                        signature_df$Clone > signature_anchor,
                                    ]
                                    new_nodes$Clone <-
                                        new_nodes$Clone +
                                        max(new_trajectory$Clone) -
                                        min(new_nodes$Clone) +
                                        1
                                    new_nodes$Parent <-
                                        new_nodes$Parent -
                                        min(new_nodes$Parent) +
                                        clone_id
                                    new_trajectory <- rbind(
                                        new_trajectory,
                                        new_nodes
                                    )
                                    new_trajectory$Patients <- paste0(
                                        'c("',
                                        paste(common_patients,
                                            collapse =
                                                '", "'
                                        ),
                                        '")'
                                    )
                                    break()
                                }
                            }
                        }
                        if (nrow(new_trajectory) == nrow(trajectory_df)) {
                            return(NULL)
                        }
                        new_trajectory$ID <- ID
                        ID <- ID + 1
                        new_trajectories_df <- rbind(
                            new_trajectories_df,
                            new_trajectory
                        )
                    }
                }
                return(new_trajectories_df)
            } else {
                return(NULL)
            }
        }
    }
}

is_trajectory_in_patient <- function(patient, trajectory, signature_type) {
    trajectory_graph <- df_to_igraph(trajectory)
    patient_graph <- df_to_igraph(patient)
    patient_graph <- prune_labels(patient_graph, unique(trajectory$Variant))
    possible_sub_graphs <- generate_possible_pruned_graphs(
        patient_graph,
        c(trajectory[["Variant"]], "g")
    )
    for (graph in possible_sub_graphs) {
        g_clean <- remove_empty_nodes(graph)
        if (signature_type == 2) {
            extensions <- dissolve_multivariant_nodes(g_clean)
            for (g in extensions) {
                all_labels <- unique(c(V(g)$label, V(trajectory_graph)$label))
                label_to_color <- setNames(seq_along(all_labels), all_labels)
                color1 <- label_to_color[V(g)$label]
                color2 <- label_to_color[V(trajectory_graph)$label]
                if (graph.isomorphic.vf2(
                    g,
                    trajectory_graph,
                    vertex.color1 = color1,
                    vertex.color2 = color2
                )$iso) {
                    return(TRUE)
                }
            }
        } else {
            if (vcount(g_clean) == vcount(trajectory_graph)) {
                all_labels <- unique(c(V(g_clean)$label, V(
                    trajectory_graph
                )$label))
                label_to_color <- setNames(seq_along(all_labels), all_labels)

                color1 <- label_to_color[V(g_clean)$label]
                color2 <- label_to_color[V(trajectory_graph)$label]

                if (graph.isomorphic.vf2(
                    g_clean,
                    trajectory_graph,
                    vertex.color1 = color1,
                    vertex.color2 = color2
                )$iso) {
                    return(TRUE)
                }
            }
        }
    }
    return(FALSE)
}

df_to_igraph <- function(df) {
    node_labels <- df %>%
        group_by(.data$Clone) %>%
        summarise(
            label = paste(.data$Variant, collapse = ";"), .groups =
                "drop"
        )
    if (!0 %in% node_labels$Clone) {
        node_labels <- bind_rows(
            data.frame(Clone = 0, label = "g"),
            node_labels
        )
    }
    edges <- df %>%
        distinct(.data$Parent, .data$Clone)
    g <- graph_from_data_frame(edges, directed = TRUE, vertices = node_labels)
    return(g)
}

remove_empty_nodes <- function(g) {
    empty_nodes <- V(g)[V(g)$label == ""]
    if (length(empty_nodes) > 0) {
        for (node in empty_nodes) {
            parent <- neighbors(g, node, mode = "in")
            children <- neighbors(g, node, mode = "out")
            if (length(parent) > 0 && length(children) > 0) {
                p_id <- as.numeric(parent[1])
                for (child in children) {
                    c_id <- as.numeric(child)
                    if (!are_adjacent(g, p_id, c_id)) {
                        g <- add_edges(g, c(p_id, c_id))
                    }
                }
            }
        }
        empty_ids <- as.numeric(empty_nodes)
        g <- delete_vertices(g, empty_ids)
    }
    return(g)
}

prune_labels <- function(g, keep_list) {
    V(g)$label <- vapply(V(g)$label, function(label) {
        fragments <- unlist(strsplit(label, ";"))
        filtered_fragments <- fragments[
            fragments %in% keep_list | fragments == "g"
        ]
        new_label <- paste(filtered_fragments, collapse = ";")
        return(new_label)
    }, character(1))
    g <- remove_empty_nodes(g)
    return(g)
}

generate_possible_pruned_graphs <- function(g, variant_list) {
    graph_fragments <- unlist(strsplit(V(g)$label, ";"))
    graph_fragments_table <- table(graph_fragments)
    target_fragments_table <- table(variant_list)
    fragments_to_remove <- list()
    for (fragment in names(graph_fragments_table)) {
        if (graph_fragments_table[fragment] > target_fragments_table[
            fragment
        ]) {
            fragments_to_remove[[fragment]] <-
                graph_fragments_table[fragment] -
                target_fragments_table[fragment]
        }
    }
    if (length(fragments_to_remove) == 0) {
        return(list(g))
    }

    get_removal_assignments <- function(required, max_removals) {
        if (length(max_removals) == 0) {
            if (required == 0) {
                return(list(numeric(0)))
            } else {
                return(list())
            }
        }
        assignments <- list()
        for (i in 0:min(max_removals[1], required)) {
            sub_assignments <- get_removal_assignments(
                required - i,
                max_removals[-1]
            )
            if (length(sub_assignments) > 0) {
                for (sa in sub_assignments) {
                    assignments[[length(assignments) + 1]] <- c(i, sa)
                }
            }
        }
        return(assignments)
    }

    remove_fragments_from_node <- function(node_label,
                                            fragment,
                                            count_to_remove) {
        fragments <- unlist(strsplit(node_label, ";"))
        to_remove <- rep(FALSE, length(fragments))
        removed <- 0
        for (j in seq_along(fragments)) {
            if (fragments[j] == fragment && removed < count_to_remove) {
                to_remove[j] <- TRUE
                removed <- removed + 1
            }
        }
        new_fragments <- fragments[!to_remove]
        new_label <- paste(new_fragments, collapse = ";")
        return(new_label)
    }

    generate_combinations <- function(graph,
                                        fragments_to_remove,
                                        current_index = 1) {
        if (current_index > length(fragments_to_remove)) {
            return(list(graph))
        }
        fragment <- names(fragments_to_remove)[current_index]
        remove_count <- fragments_to_remove[[fragment]]
        nodes_with_fragment <- V(graph)[vapply(V(graph), function(node) {
            sum(unlist(strsplit(V(graph)$label[node], ";")) == fragment) > 0
        }, logical(1))]
        if (length(nodes_with_fragment) == 0) {
            return(generate_combinations(
                graph, fragments_to_remove, current_index + 1
            ))
        }
        available_counts <- vapply(nodes_with_fragment, function(node) {
            sum(unlist(strsplit(V(graph)$label[node], ";")) == fragment)
        }, integer(1))
        assignments <- get_removal_assignments(remove_count, available_counts)
        graphs <- list()
        for (assignment in assignments) {
            graph_copy <- graph_from_data_frame(
                d = igraph::as_data_frame(graph, what = "edges"),
                vertices = igraph::as_data_frame(graph, what = "vertices"),
                directed = TRUE
            )
            for (j in seq_along(nodes_with_fragment)) {
                node_index <- as.numeric(nodes_with_fragment[j])
                num_to_remove <- assignment[j]
                if (num_to_remove > 0) {
                    old_label <- V(graph_copy)$label[node_index]
                    new_label <- remove_fragments_from_node(
                        old_label,
                        fragment,
                        num_to_remove
                    )
                    V(graph_copy)$label[node_index] <- new_label
                }
            }
            graphs <- c(graphs, generate_combinations(
                graph_copy, fragments_to_remove, current_index + 1
            ))
        }
        graphs
    }
    generate_combinations(g, fragments_to_remove)
}

dissolve_multivariant_nodes <- function(g) {
    variant_split <- strsplit(V(g)$label, ";")
    offsets <- c(0L, head(cumsum(lengths(variant_split)), -1L))
    unique_variant_split <- lapply(seq_along(variant_split), function(i) {
        variants <- variant_split[[i]]
        setNames(
            paste0(variants, "::::::", offsets[i] + seq_along(variants)),
            variants
        )
    })
    V(g)$label <- vapply(unique_variant_split, function(v) {
        paste(v, collapse = ";")
    }, character(1))
    multi_nodes <- which(vapply(unique_variant_split, length, integer(1)) > 1)
    if (length(multi_nodes) == 0) {
        V(g)$label <- gsub("::::::\\d+", "", V(g)$label)
        return(list(g))
    }
    fragment_perms <- lapply(multi_nodes, function(i) {
        x <- unique_variant_split[[i]]
        perms <- gtools::permutations(
            n = length(x),
            r = length(x),
            v = x,
            repeats.allowed = FALSE
        )
        split(perms, seq_len(nrow(perms)))
    })
    all_perm_combos <- expand.grid(
        fragment_perms,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    all_graphs <- list()
    # Graph edits invalidate vertex sequences and renumber numeric IDs. Keep
    # stable names for nodes yet to be expanded and plain neighbour names before
    # deleting the original node.
    multi_node_names <- V(g)$name[multi_nodes]
    for (combo_idx in seq_len(nrow(all_perm_combos))) {
        g_new <- g
        combo <- all_perm_combos[combo_idx, ]
        for (j in seq_along(multi_nodes)) {
            node_name <- multi_node_names[j]
            permutation <- unlist(combo[[j]])
            parents <- igraph::as_ids(neighbors(g_new, node_name, mode = "in"))
            children <- igraph::as_ids(neighbors(g_new, node_name,
                mode =
                    "out"
            ))
            g_new <- delete_vertices(g_new, node_name)
            new_node_names <- paste0(node_name, "_", seq_along(permutation))
            g_new <- add_vertices(
                g_new,
                length(permutation),
                name = new_node_names,
                label = permutation
            )
            for (k in seq_len(length(new_node_names) - 1)) {
                g_new <- add_edges(g_new, c(
                    new_node_names[k],
                    new_node_names[k + 1]
                ))
            }
            if (length(parents) > 0) {
                for (p in parents) {
                    g_new <- add_edges(g_new, c(p, new_node_names[1]))
                }
            }
            if (length(children) > 0) {
                for (cn in children) {
                    g_new <- add_edges(
                        g_new,
                        c(new_node_names[length(new_node_names)], cn)
                    )
                }
            }
        }
        V(g_new)$label <- gsub("::::::\\d+", "", V(g_new)$label)
        all_graphs[[length(all_graphs) + 1]] <- g_new
    }
    return(lapply(all_graphs, function(gx) {
        V(gx)$label <- gsub("::::::\\d+", "", V(gx)$label)
        return(gx)
    }))
}

compare_trajectories <- function(df1, df2) {
    # Preserve mutation labels, edge direction and root branches. Comparing only
    # unordered non-root edge pairs incorrectly merges distinct root stars and
    # opposite evolutionary orders.
    g1 <- df_to_igraph(df1[order(df1$Variant), , drop = FALSE])
    g2 <- df_to_igraph(df2[order(df2$Variant), , drop = FALSE])
    if (vcount(g1) != vcount(g2) || igraph::ecount(g1) != igraph::ecount(g2)) {
        return(FALSE)
    }
    labels <- unique(c(V(g1)$label, V(g2)$label))
    graph.isomorphic.vf2(g1, g2,
        vertex.color1 = match(V(g1)$label, labels),
        vertex.color2 = match(V(g2)$label, labels)
    )$iso
}

generate_shape_statistics <- function(
    df,
    existing_shape_lookup,
    n,
    signature_type
) {
    patient_dfs <- split(df, df$Patient.ID)
    patient_shapes <- list()
    all_graph_types <- unique(unlist(lapply(
        existing_shape_lookup,
        function(x) x$id
    )))
    shape_statistics <- data.frame(
        Patient.ID = integer(0),
        matrix(ncol = length(all_graph_types), nrow = 0)
    )
    colnames(shape_statistics)[-1] <- all_graph_types
    x <- 1
    for (patient_id in names(patient_dfs)) {
        x <- x + 1
        patient_df <- patient_dfs[[patient_id]]
        all_nodes <- unique(as.character(patient_df$Clone))
        root_node <- setdiff(unique(as.character(patient_df$Parent)), all_nodes)
        if (length(root_node) != 1) {
            stop(
                "Could not uniquely identify root node."
            )
        }
        root_node <- root_node[1]

        all_graphs <- list()
        graph_weights <- list()
        patient_shape_count <- setNames(
            rep(0, length(all_graph_types)),
            all_graph_types
        )
        for (i in seq.int(2, n + 1)) {
            if (signature_type == 2) {
                extended_combinations <- combinations(
                    n = length(all_nodes),
                    r = i,
                    v = all_nodes,
                    repeats.allowed = TRUE
                )
            } else {
                if (i <= length(all_nodes)) {
                    extended_combinations <- combinations(
                        n = length(all_nodes),
                        r = i,
                        v = all_nodes,
                        repeats.allowed = FALSE
                    )
                } else {
                    break
                }
            }
            extended_combinations <- split(
                extended_combinations,
                seq(nrow(extended_combinations))
            )
            for (chosen_clones in extended_combinations) {
                clone_counts <- table(chosen_clones)
                chosen_set <- unique(chosen_clones)

                valid <- all(vapply(names(clone_counts), function(clone) {
                    sum(patient_df$Clone == clone) >= clone_counts[[clone]]
                }, logical(1)))
                if (!valid) next

                clone_occurrences <- rep(0, length(chosen_clones))
                names(clone_occurrences) <- chosen_clones
                weight <- 1
                for (j in seq_along(chosen_clones)) {
                    clone <- chosen_clones[[j]]
                    clone_occurrences[clone] <- clone_occurrences[clone] + 1
                    available <- sum(patient_df$Clone == clone)
                    weight <- weight * (available - (clone_occurrences[
                        clone
                    ] - 1))
                }
                g_sub <- make_empty_graph(directed = TRUE)
                clone_map <- list()
                for (clone in chosen_set) {
                    count <- clone_counts[[clone]]
                    if (count == 1) {
                        clone_map[[clone]] <- clone
                        g_sub <- add_vertices(g_sub, 1, name = clone)
                    } else {
                        seq_names <- paste0(clone, "_", seq_len(count))
                        clone_map[[clone]] <- seq_names
                        g_sub <- add_vertices(g_sub, count, name = seq_names)
                        g_sub <- add_edges(g_sub, as.vector(t(embed(
                            seq_names, 2
                        )[, 2:1])))
                    }
                }

                g_sub <- add_vertices(g_sub, 1, name = root_node)

                for (clone in chosen_set) {
                    count <- clone_counts[[clone]]
                    first_clone_node <- if (count == 1) {
                        clone_map[[clone]]
                    } else {
                        head(clone_map[[clone]], 1)
                    }
                    parent <- as.character(
                        patient_df$Parent[patient_df$Clone == clone][1]
                    )
                    if (is.na(parent)) next
                    current <- parent
                    while (!(current %in% names(clone_map)) &&
                        current != root_node) {
                        current <- as.character(
                            patient_df$Parent[patient_df$Clone == current][1]
                        )
                    }
                    parent_node <- if (current == root_node) {
                        root_node
                    } else {
                        mapped <- clone_map[[current]]
                        if (length(mapped) == 1) mapped else tail(mapped, 1)
                    }
                    g_sub <- add_edges(g_sub, c(parent_node, first_clone_node))
                }

                all_graphs[[length(all_graphs) + 1]] <- g_sub
                graph_weights[[length(graph_weights) + 1]] <- weight
            }
        }

        shape_assignment <- assign_unlabeled_topologies(
            all_graphs,
            existing_shape_lookup
        )
        shape_assignment$weights <- graph_weights
        for (i in seq_along(shape_assignment$weights)) {
            shape_id <- shape_assignment$shape_id[i]
            weight <- shape_assignment$weights[[i]]

            if (!shape_id %in% names(patient_shape_count)) {
                patient_shape_count[shape_id] <- 0
            }
            patient_shape_count[shape_id] <- patient_shape_count[shape_id] +
                weight
        }
        patient_shapes[[patient_id]] <- patient_shape_count
    }
    empty_df <- data.frame(
        Patient.ID = character(0),
        matrix(0, nrow = 0, ncol = length(names(existing_shape_lookup)))
    )
    for (patient_id in names(patient_shapes)) {
        new_row <- c(patient_id)
        for (shape_id in names(existing_shape_lookup)) {
            if (shape_id %in% names(patient_shapes[[patient_id]])) {
                new_row <- c(new_row, patient_shapes[[patient_id]][shape_id])
            } else {
                new_row <- c(new_row, 0)
            }
        }
        empty_df <- rbind(empty_df, new_row)
    }
    colnames(empty_df)[-1] <- names(existing_shape_lookup)
    return(empty_df)
}

assign_unlabeled_topologies <- function(graph_list, shape_lookup) {
    shape_ids <- character(length(graph_list))
    if (length(shape_lookup) == 0) {
        next_id <- 1
    } else {
        existing_ids <- as.integer(sub("shape_", "", names(shape_lookup)))
        next_id <- max(existing_ids) + 1
    }
    for (i in seq_along(graph_list)) {
        g <- graph_list[[i]]

        matched <- FALSE
        for (j in seq_along(shape_lookup)) {
            if (graph.isomorphic(g, shape_lookup[[j]])) {
                shape_ids[i] <- names(shape_lookup)[j]
                matched <- TRUE
                break
            }
        }

        if (!matched) {
            shape_name <- paste0("shape_", next_id)
            shape_lookup[[shape_name]] <- g
            shape_ids[i] <- shape_name
            next_id <- next_id + 1
        }
    }

    return(list(
        shape_ids = shape_ids,
        shape_lookup = shape_lookup
    ))
}

.trajectory_shape_key <- function(edges) {
    # Canonical rooted-tree encoding is independent of row order and node IDs.
    parents <- as.character(edges$Parent)
    clones <- as.character(edges$Clone)
    roots <- setdiff(unique(parents), clones)
    if (length(roots) != 1L) {
        stop(
            "Could not uniquely identify trajectory root."
        )
    }
    children <- split(clones, parents)
    encode <- function(node) {
        branches <- vapply(children[[node]], encode, character(1))
        paste0("(", paste(sort(branches), collapse = ""), ")")
    }
    encode(roots[[1]])
}

assign_shapes <- function(trajectories_df) {
    graphs <- split(trajectories_df, trajectories_df$ID)
    shape_lookup <- list()
    no_of_graphs <- length(graphs)
    shape_ids <- character(no_of_graphs)
    id_counter <- 1
    shape_key_map <- list()
    for (i in seq_along(graphs)) {
        g_df <- graphs[[i]]
        edges <- g_df %>% select(all_of(c("Parent", "Clone")))

        edge_key <- .trajectory_shape_key(edges)

        if (!(edge_key %in% names(shape_lookup))) {
            shape_id <- paste0("shape_", id_counter)
            g <- graph_from_data_frame(edges, directed = TRUE)
            shape_lookup[[edge_key]] <- list(
                id = shape_id,
                graph = g
            )
            id_counter <- id_counter + 1
        }

        shape_ids[i] <- shape_lookup[[edge_key]]$id
        shape_key_map[[as.character(g_df$ID[1])]] <- edge_key
    }

    trajectories_df$Shape <- unlist(lapply(
        trajectories_df$ID,
        function(id) shape_lookup[[shape_key_map[[as.character(id)]]]]$id
    ))

    shape_lookup_named <- setNames(
        lapply(shape_lookup, function(x) x$graph),
        vapply(shape_lookup, function(x) x$id, character(1))
    )

    return(list(
        trajectories = trajectories_df,
        shape_lookup = shape_lookup_named
    ))
}

add_p_values <- function(
    trajectory_dataframe,
    topology_dataframe,
    patient_data
) {
    unique_rows <- unique(trajectory_dataframe[, c("ID", "Patients", "Shape")])
    total_trajectory_count <- nrow(unique_rows)
    for (i in seq_len(total_trajectory_count)) {
        id <- unique_rows$ID[i]
        patients <- unique_rows$Patients[i]
        shape <- unique_rows$Shape[i]
        matching_rows <- trajectory_dataframe$ID == id &
            trajectory_dataframe$Patients == patients &
            trajectory_dataframe$Shape == shape
        patient_list <- unlist(strsplit(patients, ","))
        a <- length(patient_list)
        topology_counts <- topology_dataframe[[shape]]
        target_variants <- trajectory_dataframe$Variant[matching_rows]
        b <- numeric(length(topology_counts))
        for (i in seq_along(unique(patient_data$Patient.ID))) {
            patient_id <- unique(patient_data$Patient.ID)[i]
            observed_variants <- patient_data$Variant[
                patient_data$Patient.ID == patient_id
            ]
            observed_counts <- table(observed_variants)
            required_counts <- table(target_variants)
            observed_counts <- observed_counts[names(required_counts)]
            observed_counts[is.na(observed_counts)] <- 0
            if (all(required_counts <= observed_counts)) {
                b[i] <- 1 / (length(unique(observed_variants))^sum(
                    matching_rows
                ))
            } else {
                b[i] <- 0
            }
        }
        keep_indices <- b != 0 & topology_counts != 0
        b <- b[keep_indices]
        topology_counts <- topology_counts[keep_indices]
        if (sum(as.numeric(topology_counts)) != 0) {
            patient_probabilities <- 1 - (1 - b)^as.integer(topology_counts)
            x <- .signature_enrichment_tail(a, patient_probabilities)
        } else {
            x <- 2
        }
        trajectory_dataframe$`p value`[matching_rows] <- x
        occurrence_column <- "Number of occurrences"
        relative_occurrence_column <-
            "Relative occurrence among patients with identical variants"

        trajectory_dataframe[[occurrence_column]][matching_rows] <-
            as.character(unique(
                trajectory_dataframe[[occurrence_column]][matching_rows]
            ))
        trajectory_dataframe[[relative_occurrence_column]][matching_rows] <-
            as.numeric(unique(
                trajectory_dataframe[[occurrence_column]][matching_rows]
            )) /
                sum(b != 0)
    }
    return(trajectory_dataframe)
}

# Enumerate unordered root stars, pruning each level by native patient support.
# Canonical variant indices avoid factorial permutations; repeated labels are
# bounded by observed event multiplicity and permitted only when requested.
process_indi_patient_branches <- function(
    df,
    trajectories_bases,
    signature_type,
    minimum_support = 1L,
    allow_repeated = FALSE
) {
    stopifnot(
        length(minimum_support) == 1L, is.finite(minimum_support),
        minimum_support >= 1
    )
    if (!nrow(df)) {
        return(list())
    }
    patients <- split(df, df$Patient.ID)
    variants <- sort(unique(as.character(df$Variant)))
    if (!length(variants)) {
        return(list())
    }
    max_size <- max(vapply(patients, nrow, integer(1)))
    observed_counts <- lapply(patients, function(patient) {
        table(
            patient$Variant
        )
    })
    support <- lapply(variants, function(variant) {
        names(patients)[vapply(patients, function(patient) {
            variant %in% patient$Variant
        }, logical(1))]
    })
    bases <- lapply(which(lengths(support) >= minimum_support), function(i) {
        list(indices = i, patients = support[[i]])
    })
    results <- list()
    next_id <- if (nrow(trajectories_bases)) max(trajectories_bases$ID) else 0
    while (length(bases) && length(bases[[1]]$indices) < max_size) {
        new_bases <- list()
        for (base in bases) {
            first <- tail(base$indices, 1) + if (allow_repeated) 0L else 1L
            if (first > length(variants)) next
            for (i in seq.int(first, length(variants))) {
                indices <- c(base$indices, i)
                labels <- variants[indices]
                candidates <- intersect(base$patients, support[[i]])
                if (length(candidates) < minimum_support) next
                trajectory <- data.frame(
                    ID = 0L, Clone = seq_along(labels),
                    Parent = 0L, Variant = labels, Patients = "",
                    stringsAsFactors = FALSE
                )
                required <- table(labels)
                keep <- candidates[vapply(candidates, function(patient_id) {
                    counts <- observed_counts[[patient_id]][names(required)]
                    counts[is.na(counts)] <- 0L
                    all(counts >= required) && is_trajectory_in_patient(
                        patients[[patient_id]], trajectory, signature_type
                    )
                }, logical(1))]
                if (length(keep) < minimum_support) next
                next_id <- next_id + 1L
                trajectory$ID <- next_id
                # dput supplies a correctly escaped R character vector for
                # downstream
                # native support intersections, including quotes in patient
                # identifiers.
                trajectory$Patients <- paste(capture.output(dput(keep)),
                    collapse = ""
                )
                results[[length(results) + 1L]] <- trajectory
                new_bases[[length(new_bases) + 1L]] <- list(
                    indices = indices, patients = keep
                )
            }
        }
        bases <- new_bases
    }
    results
}
