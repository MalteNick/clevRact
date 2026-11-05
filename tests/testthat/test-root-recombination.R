root_cohort <- function(labels, parents, ids = c("p1", "p2")) {
  do.call(rbind, lapply(ids, function(id) data.frame(
    Patient.ID = id, Variant = labels, Clone = seq_along(labels),
    Parent = parents, stringsAsFactors = FALSE)))
}
root_input <- function(repeats = FALSE, rounds = 1L) list(
  minimum_support_of_trajectories = 2L,
  multiple_identical_variants_for_branches = repeats,
  recombination_number_setter = rounds)
root_seed <- function(a, b, ids = c("p1", "p2")) {
  data.frame(`Variant 1` = a, `Variant 2` = b,
    `Number of occurrences` = length(ids), `p value` = 0.5,
    Patients = paste(capture.output(dput(ids)), collapse = ""),
    check.names = FALSE)
}
root_groups <- function(result) split(result, result$ID)

test_that("independent mutations are found without any linear signatures", {
  cohort <- root_cohort(c("IDH2", "NRAS"), c(0L, 0L),
    c("patient'one", 'patient"two'))
  result <- compute_trajectories(data.frame(), list(cohort), root_input(), 1L)
  expect_equal(length(root_groups(result)), 1L)
  expect_setequal(result$Variant, c("IDH2", "NRAS"))
  expect_true(all(result$Parent == 0))
  expect_true(all(result$`Number of occurrences` == 2L))
  expect_setequal(eval(parse(text = result$Patients[1])), unique(cohort$Patient.ID))
})

test_that("root initialization enumerates combinations once and respects support", {
  cohort <- root_cohort(c("A", "B", "C"), c(0L, 0L, 0L))
  # A rare label must not produce a recurrent root base.
  cohort <- rbind(cohort, data.frame(Patient.ID = "p1", Variant = "rare",
    Clone = 4L, Parent = 0L))
  result <- compute_trajectories(data.frame(), list(cohort), root_input(), 1L)
  groups <- root_groups(result)
  expect_equal(unname(sort(lengths(lapply(groups, function(x) x$Variant)))), c(2L, 2L, 2L, 3L))
  expect_false("rare" %in% result$Variant)
  keys <- vapply(groups, function(x) paste(sort(x$Variant), collapse = ";"), "")
  expect_equal(anyDuplicated(keys), 0L)
})

test_that("same-clone and linear events are not mistaken for independent branches", {
  linear <- root_cohort(c("A", "B"), c(0L, 1L))
  expect_equal(nrow(compute_trajectories(data.frame(), list(linear), root_input(), 1L)), 0L)
  same_clone <- root_cohort(c("A", "B"), c(0L, 0L))
  same_clone$Clone <- 1L
  expect_equal(nrow(compute_trajectories(data.frame(), list(same_clone), root_input(), 2L)), 0L)
})

test_that("independent bases are extended by the ordinary linear signatures", {
  cohort <- root_cohort(c("A", "B", "C"), c(0L, 0L, 1L))
  result <- compute_trajectories(root_seed("A", "C"), list(cohort), root_input(), 1L)
  groups <- root_groups(result)
  target <- Filter(function(x) setequal(x$Variant, c("A", "B", "C")), groups)
  expect_true(length(target) >= 1L)
  expect_true(any(vapply(target, function(x) {
    a <- x$Clone[x$Variant == "A"]
    all(x$Parent[x$Variant %in% c("A", "B")] == 0) &&
      x$Parent[x$Variant == "C"] == a
  }, logical(1))))
  expect_true(all(result$`Number of occurrences` >= 2))
})

test_that("repeated root labels require opt-in and observed multiplicity", {
  cohort <- root_cohort(c("A", "A"), c(0L, 0L))
  expect_equal(nrow(compute_trajectories(data.frame(), list(cohort), root_input(), 1L)), 0L)
  result <- compute_trajectories(data.frame(), list(cohort), root_input(TRUE), 1L)
  expect_equal(nrow(result), 2L)
  expect_true(all(result$Variant == "A"))
  expect_true(all(result$`Number of occurrences` == 2L))
})


test_that("deduplication preserves root labels and directed evolutionary order", {
  star <- data.frame(Clone=c(1L,2L),Parent=c(0L,0L),Variant=c("A","B"))
  other <- star; other$Variant <- c("A","C")
  expect_false(compare_trajectories(star,other))
  forward <- star; forward$Parent <- c(0L,1L)
  reverse <- forward; reverse$Variant <- rev(forward$Variant)
  expect_false(compare_trajectories(forward,reverse))
  renumbered <- forward; renumbered$Clone <- c(7L,9L);renumbered$Parent <- c(0L,7L)
  expect_true(compare_trajectories(forward,renumbered[2:1,]))
})


test_that("independent-root computation can be switched off without disabling ordinary branches", {
  input <- root_input()
  input$compute_independent_root_branches <- FALSE
  independent <- root_cohort(c("A", "B"), c(0L, 0L))
  expect_equal(nrow(compute_trajectories(data.frame(), list(independent), input, 1L)), 0L)
  input$compute_independent_root_branches <- TRUE
  expect_equal(nrow(compute_trajectories(data.frame(), list(independent), input, 1L)), 2L)
  input$compute_independent_root_branches <- FALSE
  ordinary <- root_cohort(c("A", "B", "C"), c(0L, 1L, 1L))
  seeds <- rbind(root_seed("A", "B"), root_seed("A", "C"))
  result <- compute_trajectories(seeds, list(ordinary), input, 1L)
  expect_true(any(vapply(root_groups(result), function(x) {
    if (!setequal(x$Variant, c("A", "B", "C"))) return(FALSE)
    a <- x$Clone[x$Variant == "A"]
    all(x$Parent[x$Variant %in% c("B", "C")] == a)
  }, logical(1))))
})


test_that("independent-root selection accepts Shiny radio-button string values", {
  cohort <- root_cohort(c("A", "B"), c(0L, 0L))
  input <- root_input()
  input$compute_independent_root_branches <- "FALSE"
  expect_equal(nrow(compute_trajectories(data.frame(), list(cohort), input, 1L)), 0L)
  input$compute_independent_root_branches <- "TRUE"
  expect_equal(nrow(compute_trajectories(data.frame(), list(cohort), input, 1L)), 2L)
})

test_that("repeated-root selection accepts Shiny radio-button string values", {
  cohort <- root_cohort(c("A", "A"), c(0L, 0L))
  input <- root_input()
  input$compute_independent_root_branches <- "TRUE"
  input$multiple_identical_variants_for_branches <- "TRUE"
  expect_equal(nrow(compute_trajectories(data.frame(), list(cohort), input, 1L)), 2L)
  input$multiple_identical_variants_for_branches <- "FALSE"
  expect_equal(nrow(compute_trajectories(data.frame(), list(cohort), input, 1L)), 0L)
})

test_that("the GUI computes independent branches with an empty linear seed table", {
  shiny::testServer(.build_server(), {
    cohort <- root_cohort(c("A", "A"), c(0L, 0L))
    state$signature_data_set <- data.frame(`Variant 1` = character(),
      `Number of occurrences` = numeric(), Patients = character(),
      `p value` = numeric(),
      `Relative occurrence among patients with identical variants` = numeric(),
      check.names = FALSE)
    state$preprocessed_data <- list(cohort)
    state$signature_type <- 1L
    session$setInputs(absolute_frequency_signatures = 2,
      signature_search = "", min_p_value = 1, relative_frequency = 0,
      minimum_support_of_trajectories = 2L, recombination_number_setter = 1L,
      compute_independent_root_branches = "TRUE",
      multiple_identical_variants_for_branches = "TRUE",
      compute_trajectories = 1L)
    expect_equal(nrow(state$trajectory_data_set), 2L)
    expect_true(all(state$trajectory_data_set$Variant == "A"))
  })
})
