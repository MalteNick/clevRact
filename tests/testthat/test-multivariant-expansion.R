test_that('expansion preserves ancestry across multiple multivariant clones', {
  df <- data.frame(Patient.ID='p',Variant=c('A','B','C','D'),
    Clone=c(1L,1L,2L,2L),Parent=c(0L,0L,1L,1L))
  graphs <- dissolve_multivariant_nodes(df_to_igraph(df))
  expect_length(graphs,4L)
  for(g in graphs) {
    expect_setequal(igraph::V(g)$label,c('g','A','B','C','D'))
    expect_true(igraph::is_dag(g))
    labels <- igraph::V(g)$label
    distance <- igraph::distances(g,mode='out')
    expect_true(all(is.finite(distance[match(c('A','B'),labels),match(c('C','D'),labels)])))
  }
})

test_that('same-clone expansion handles repeated labels and stable neighbour IDs', {
  patient <- data.frame(Patient.ID='p',Variant=c('A','B','A','C'),
    Clone=c(1L,1L,2L,2L),Parent=c(0L,0L,1L,1L))
  target <- data.frame(Variant=c('B','A','C'),Clone=1:3,Parent=c(0L,1L,2L))
  expect_true(is_trajectory_in_patient(patient,target,2L))
  expect_false(is_trajectory_in_patient(patient,target,1L))
})
