test_that("enrichment uses the inclusive upper tail", {
  expect_equal(.signature_enrichment_tail(2, c(0.1, 0.2)), 0.02)
  expect_equal(.signature_enrichment_tail(1, c(0.1, 0.2)), 0.28)
  expect_equal(.signature_enrichment_tail(0, c(0.1, 0.2)), 1)
  expect_equal(.signature_enrichment_tail(3, c(0.1, 0.2)), 0)
  expect_equal(.signature_enrichment_tail(0, numeric()), 1)
  expect_equal(.signature_enrichment_tail(1, numeric()), 0)
  expect_equal(.signature_enrichment_tail(1, c(0, 1)), 1)
  expect_equal(.signature_enrichment_tail(2, c(0, 1)), 0)
})

test_that("enrichment agrees with binomial and enumerated probabilities", {
  for (observed in 0:21) {
    expect_equal(.signature_enrichment_tail(observed, rep(0.2, 20)),
      pbinom(observed - 1, 20, 0.2, lower.tail = FALSE))
  }
  probabilities <- c(0.02, 0.3, 0.6, 0.9)
  outcomes <- as.matrix(expand.grid(rep(list(0:1), length(probabilities))))
  masses <- apply(outcomes, 1, function(x) {
    prod(ifelse(x == 1, probabilities, 1 - probabilities))
  })
  for (observed in 0:5) {
    expect_equal(.signature_enrichment_tail(observed, probabilities),
      sum(masses[rowSums(outcomes) >= observed]))
  }
})

test_that("signature and trajectory scoring return the same enrichment tail", {
  patients <- data.frame(Patient.ID = rep(c("p1", "p2"), each = 3),
    Variant = rep(c("A", "B", "C"), 2), stringsAsFactors = FALSE)
  signatures <- data.frame(`Variant 1` = "A", `Variant 2` = "B",
    check.names = FALSE)
  observed <- hash::hash()
  observed[["A;B"]] <- c("p1", "p2")
  counts <- hash::hash()
  counts[["p1"]] <- c(0, 1)
  counts[["p2"]] <- c(0, 1)
  statistics <- compute_alt_statistics(signatures, patients, observed, counts,
    list(significance_estimation_method = "Sandmann"))
  signature_p <- statistics[[1]][["A;B"]]
  # Two independent patients, each with q = 1/3^2 and one candidate position.
  expect_equal(signature_p, 1 / 81)
  trajectory <- data.frame(ID = c(1, 1), Patients = "p1,p2", Shape = "1",
    Variant = c("A", "B"), `Number of occurrences` = "2",
    check.names = FALSE, stringsAsFactors = FALSE)
  scored <- add_p_values(trajectory, list(`1` = c(1, 1)), patients)
  expect_equal(scored$`p value`, rep(signature_p, 2))
})

test_that("equivalent rooted shapes share topology counts", {
  first <- data.frame(ID = 1, Clone = 1:4, Parent = c(0, 1, 1, 2),
    Variant = LETTERS[1:4], Patients = "p1,p2",
    `Number of occurrences` = "2", check.names = FALSE)
  second <- first[c(4, 3, 2, 1), ]
  second$ID <- 2
  mapping <- c(`0` = 0, `1` = 8, `2` = 5, `3` = 9, `4` = 2)
  second$Clone <- unname(mapping[as.character(second$Clone)])
  second$Parent <- unname(mapping[as.character(second$Parent)])
  chain <- first
  chain$ID <- 3
  chain$Parent <- c(0, 1, 2, 3)
  shapes <- assign_shapes(rbind(first, second, chain))
  expect_equal(length(shapes$shape_lookup), 2L)
  ids <- vapply(split(shapes$trajectories$Shape,
    shapes$trajectories$ID), function(x) unique(x), character(1))
  expect_equal(unname(ids[1]), unname(ids[2]))
  expect_equal(unname(ids[1] == ids[3]), FALSE)
  patients <- do.call(rbind, lapply(c("p1", "p2"), function(id) {
    data.frame(Patient.ID = id, Clone = 1:4, Parent = c(0, 1, 1, 2),
      Variant = LETTERS[1:4])
  }))
  shapes <- assign_shapes(rbind(first, second))
  counts <- generate_shape_statistics(patients, shapes$shape_lookup, 4, 1L)
  scored <- add_p_values(shapes$trajectories, counts, patients)
  expect_equal(scored$`p value`, rep((1 / 4^4)^2, 8))
})
