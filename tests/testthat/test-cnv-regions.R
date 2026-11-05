test_that("overlaps form inclusive atomic regions without joining adjacent intervals", {
  result <- cnv_region_partition(c("CNV_del_9_100_200", "CNV_del_9_150_250", "CNV_del_9_251_300", "CNV_loh_9_150_250"))
  del <- result$regions[result$regions$type == "CNV_del", ]
  expect_equal(del$start, c(100,150,201,251))
  expect_equal(del$end, c(149,200,250,300))
  expect_equal(del$source_intervals, c(1L,2L,1L,1L))
  expect_true("CNV_loh_region_9_150_200" %in% result$regions$region)
})

test_that("postprocessing pools patients and treats segments as alternatives", {
  signatures <- data.frame(`Variant 1`=c("CNV_del_9_100_200","CNV_del_9_150_250"),
    `Variant 2`="NOTCH1",`Number of occurrences`=c(2,2),
    Patients=c('c("p1", "p2")','c("p2", "p3")'),`p value`=.001,check.names=FALSE)
  data <- data.frame(Variant=signatures$`Variant 1`)
  original <- data
  result <- postprocess_cnv_regions(list(signatures,list(data)),minimum_support=1)
  target <- result$signatures[result$signatures$`Variant 1`=="CNV_del_region_9_150_200",]
  expect_equal(target$`Number of occurrences`,3)
  expect_equal(target$`Source signatures`,2L)
  expect_equal(nrow(result$signatures),3L)
  expect_true(all(result$signatures$`Variant 2`=="NOTCH1"))
  expect_false("Variant 3" %in% names(result$signatures))
  expect_false("p value" %in% names(result$signatures))
  expect_identical(data,original)
})

test_that("full computed catalogue pools singleton CNVs after extraction", {
  map <- hash::hash();map[["CNV_del_9_100_200;NOTCH1"]]<-"p1"
  map[["CNV_del_9_150_250;NOTCH1"]]<-"p2"
  data<-data.frame(Variant=c("CNV_del_9_100_200","CNV_del_9_150_250"))
  result<-postprocess_cnv_regions(list(data.frame(),list(data,map)),minimum_support=2)
  expect_equal(nrow(result$signatures),1L)
  expect_equal(result$signatures$`Number of occurrences`,2)
  expect_equal(result$signatures$`Variant 1`,"CNV_del_region_9_150_200")
})

test_that("malformed coordinates and patient expressions fail explicitly", {
  expect_error(cnv_region_partition("CNV_del_9_200_100"),"Invalid")
  expect_error(cnv_region_partition("CNV_del_in_9p"),"exact")
  sig<-data.frame(`Variant 1`="CNV_del_9_1_2",`Variant 2`="A",Patients='system("false")',check.names=FALSE)
  expect_error(postprocess_cnv_regions(list(sig,list(data.frame(Variant=sig$`Variant 1`)))),"literal")
})

test_that("subgroups reuse the complete cohort's globally defined regions", {
  partition<-cnv_region_partition(c("CNV_del_9_100_200","CNV_del_9_150_250"))
  sig<-data.frame(`Variant 1`="CNV_del_9_100_200",`Variant 2`="NOTCH1",
    Patients='c("p1", "p2")',check.names=FALSE)
  result<-postprocess_cnv_regions(list(sig,list(data.frame(Variant=sig$`Variant 1`))),
    minimum_support=2,region_partition=partition)
  expect_equal(result$signatures$`Variant 1`,c("CNV_del_region_9_100_149","CNV_del_region_9_150_200"))
  expect_identical(result$regions,partition$regions)
  expect_true(all(result$signatures$`Number of occurrences`==2))
})

test_that("the optional Shiny action preserves original signatures and tree data", {
  shiny::testServer(.build_server(), {
    sig<-data.frame(`Variant 1`="CNV_del_9_100_200",`Variant 2`="NOTCH1",
      `Number of occurrences`=2,Patients='c("p1", "p2")',check.names=FALSE)
    original<-data.frame(Patient.ID=c("p1","p1","p2","p2"),
      Variant=rep(c("CNV_del_9_100_200","NOTCH1"),2),Clone=rep(1:2,2),Parent=rep(0:1,2))
    state$signature_data_set<-sig
    state$signature_side_information<-list(original)
    state$preprocessed_data<-list(original)
    session$setInputs(minimum_support_of_signatures=2,postprocess_cnv_regions=1)
    expect_equal(nrow(state$cnv_region_result$signatures),1L)
    expect_identical(state$signature_data_set,sig)
    expect_identical(state$preprocessed_data[[1]],original)
  })
})
