createSeaObjectFromDataset <- function(dataset) {
    if (!("T1" %in% colnames(dataset))) {
        stop("Cannot generate plots without timepoint data")
    }
    dataset <- mark_dead_branches(dataset)
    if (!"Clone" %in% colnames(dataset) || !"Parent" %in% colnames(dataset)) {
        stop("Dataset must include 'Clone' and 'Parent' columns.")
    }
    clones <- unique(dataset$Clone)
    parents <- vapply(clones, function(clone) {
        unique(dataset$Parent[dataset$Clone == clone])
    }, numeric(1))
    fracTable <- matrix(NA, nrow = length(parents), ncol = 1)
    fracTable[1, 1] <- 100
    fracTable <- compute_tree_matrix(dataset)
    timepoints <- seq(10, 10 * length(grep("^T\\d+$", names(dataset))), by = 10)
    seaObject <- createSeaObject(
        fracTable = fracTable,
        parents = parents,
        timepoints = timepoints,
        timepointInterpolation = TRUE,
        originTimepoint = 0
    )
    return(seaObject = seaObject)
}

distribute_ccf <- function(parent, parent_ccf, parents, fracTable, level = 1) {
    children <- which(parents == parent)
    if (length(children) > 0) {
        child_ccf <- parent_ccf * 0.8 / length(children)
        for (child in children) {
            fracTable[child, 1] <- child_ccf
            fracTable <- distribute_ccf(
                child,
                child_ccf,
                parents,
                fracTable,
                level + 1
            )
        }
    }
    return(fracTable)
}

compute_tree_matrix <- function(df) {
    # Select relevant columns and remove duplicates
    t_cols <- grep("^T\\d+$", names(df), value = TRUE)
    dead_cols <- grep("^Dead_T\\d+$", names(df), value = TRUE)
    df_trimmed <- df %>%
        select(all_of(c("Clone", "Parent", t_cols, dead_cols))) %>%
        distinct()

    # Step 2: Build tree structure
    tree <- split(df_trimmed$Clone, df_trimmed$Parent)
    clones <- setdiff(unique(df_trimmed$Clone), 0)
    result_matrix <- matrix(0,
        nrow = length(clones), ncol = length(t_cols),
        dimnames = list(as.character(clones), t_cols)
    )

    # Return each updated matrix explicitly from the recursive traversal.
    propagate <- function(node, value, t_index, result_matrix) {
        children <- tree[[as.character(node)]]
        if (is.null(children)) {
            return(result_matrix)
        }
        if (node == 0) {
            parent_active <- TRUE
        } else {
            parent_row <- df_trimmed %>% filter(.data$Clone == node)
            parent_active <- parent_row[[t_cols[t_index]]]
            parent_active <- ifelse(is.na(parent_active), TRUE, parent_active)
        }
        alive_children <- Filter(function(child) {
            row <- df_trimmed %>% filter(.data$Clone == child)
            is_dead <- row[[dead_cols[t_index]]]
            is.na(is_dead) || !is_dead
        }, children)
        if (length(alive_children) == 0) {
            return(result_matrix)
        }
        total_value <- if (parent_active) value * 0.8 else value
        each_value <- total_value / length(alive_children)
        for (child in alive_children) {
            result_matrix[as.character(child), t_index] <- each_value
            result_matrix <- propagate(
                child, each_value, t_index, result_matrix
            )
        }
        result_matrix
    }
    for (i in seq_along(t_cols)) {
        result_matrix <- propagate(0, 100, i, result_matrix)
    }
    return(result_matrix)
}

mark_dead_branches <- function(df) {
    unique_nodes <- df %>% distinct(.data$Clone, .data$Parent, .keep_all = TRUE)
    t_cols <- grep("^T\\d+$", names(df), value = TRUE)
    for (t_col in t_cols) {
        dead_flags <- setNames(rep(NA, nrow(unique_nodes)), unique_nodes$Clone)
        is_dead_branch <- function(clone) {
            rows <- unique_nodes[unique_nodes$Clone == clone, ]
            if (nrow(rows) == 0) {
                return(TRUE)
            }
            if (any(as.logical(rows[[t_col]], na.rm = TRUE))) {
                return(FALSE)
            }
            children <- unique_nodes$Clone[unique_nodes$Parent == clone]
            all_children_dead <- all(vapply(
                children, is_dead_branch,
                logical(1)
            ))
            return(all_children_dead)
        }
        for (clone in unique_nodes$Clone) {
            dead_flags[as.character(clone)] <- is_dead_branch(clone)
        }
        unique_nodes[[paste0("Dead_", t_col)]] <- dead_flags
    }
    df_with_flags <- df %>%
        left_join(
            unique_nodes %>%
                select(all_of(c("Clone", "Parent")), starts_with("Dead_")),
            by = c("Clone", "Parent")
        )
    return(df_with_flags)
}
