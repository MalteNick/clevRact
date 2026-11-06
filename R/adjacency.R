compute_adjacency_matrix <- function(dataframe, edges_within) {
    # Innitialize adjacency matrix
    sparse_array <- SVT_SparseArray(dim = c(nrow(dataframe), nrow(dataframe)))
    # Get the original tumor-tree reachability matrix as a helper structure.
    dataframe <- dataframe[order(dataframe$Clone, dataframe$Variant), ]
    tree_structure_dataframe <- unique(dataframe[, c("Parent", "Clone")])
    graph <- graph_from_data_frame(tree_structure_dataframe, directed = TRUE)
    graph <- delete_vertices(graph, V(graph)[which(V(graph)$name == 0)])
    distance_matrix <- distances(graph, mode = "out")
    clone_order <- as.character(unique(dataframe$Clone))
    if (length(clone_order) > 1) {
        distance_matrix <- distance_matrix[clone_order, clone_order]
    }
    reachability_matrix <- ifelse(is.finite(distance_matrix), 1, 0)
    unique_values_df <- as.vector(unique(tree_structure_dataframe[, c(
        "Clone"
    )]))
    colnames(reachability_matrix) <- unique(tree_structure_dataframe$Clone)
    rownames(reachability_matrix) <- unique(tree_structure_dataframe$Clone)
    diag(reachability_matrix) <- 0
    # Setup adjacency matrix of exntended tumor trees by blocks
    clone_counts <- table(dataframe$Clone)
    first_occurence <- setNames(
        rep(NA, length(unique(dataframe$Clone))),
        unique(dataframe$Clone)
    )
    for (i in seq_len(nrow(dataframe))) {
        clone_value <- dataframe$Clone[i]
        if (is.na(first_occurence[clone_value])) {
            first_occurence[clone_value] <- i
        }
    }
    index_pairs <- which(reachability_matrix == 1, arr.ind = TRUE)
    for (pair in seq_len(nrow(index_pairs))) {
        col_index <- index_pairs[pair, "row"]
        row_index <- index_pairs[pair, "col"]
        col_name <- rownames(reachability_matrix)[row_index]
        row_name <- colnames(reachability_matrix)[col_index]
        m <- clone_counts[row_name]
        n <- clone_counts[col_name]
        start_row <- first_occurence[row_name]
        start_col <- first_occurence[col_name]
        sparse_array[
            seq.int(start_row, length.out = m),
            seq.int(start_col, length.out = n)
        ] <- 1
    }
    # Add diagonal blocks for type-2 extended tumor trees
    if (edges_within) {
        start_index <- 1
        for (clone in unique(dataframe$Clone)) {
            clone_rows <- which(dataframe$Clone == clone)
            n <- length(clone_rows)
            if (n > 0) {
                indices <- expand.grid(
                    i = seq.int(start_index, length.out = n),
                    j = seq.int(start_index, length.out = n)
                )
                sparse_array[as.matrix(indices)] <- 1
                start_index <- start_index + n
            }
        }
        sparse_array[
            cbind(seq_len(nrow(dataframe)), seq_len(nrow(dataframe)))
        ] <- 0
    }
    return(sparse_array)
}

compute_n_walk_matrix <- function(prev_array, adjacency) {
    non_zero_entries <- nzwhich(prev_array, arr.ind = TRUE)
    dim <- length(dim(prev_array)) + 1
    n <- dim(prev_array)[[1]]
    new_pathways_matrix <- SVT_SparseArray(dim = c(rep(n, dim)))
    for (entry in seq_len(nrow(non_zero_entries))) {
        prev_indices <- non_zero_entries[entry, seq_len(dim - 1)]
        last_index <- prev_indices[dim - 1]
        for (new_index in seq_len(n)) {
            if (adjacency[last_index, new_index] == 1) {
                if (length(unique(c(prev_indices, new_index))) == dim) {
                    coordinates <- rev(c(
                        non_zero_entries[entry, seq_len(dim - 1)],
                        new_index
                    ))
                    linear_coordinate <- 1 + sum(
                        (coordinates - 1) *
                            rev(cumprod(c(1, rep(n, length(coordinates) - 1))))
                    )
                    new_pathways_matrix[linear_coordinate] <- 1
                }
            }
        }
    }
    return(new_pathways_matrix)
}
