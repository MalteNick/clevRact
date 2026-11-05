test_that('reference annotation retains X chromosome small variants', {
  for(reference in c('GRCh37','GRCh38')) {
    genes <- read.csv(.annotation_file(paste0(reference,'.csv')),sep=';')
    rows <- genes[genes$Symbol %in% c('PHF6','KDM6A') & genes$Chromosome=='X',]
    rows <- rows[!duplicated(rows$Symbol),]
    expect_equal(nrow(rows),2L)
    data <- data.frame(Patient.ID=c('p1','p2'),Variant=rows$Symbol,
      Clone=1L,Parent=0L,Chromosome='X',Start=rows$Begin,End=rows$Begin,
      Variant.Type='SNV/INDEL')
    input <- list(dropdown_variant_generation=paste('Yes, using',reference),
      dropdown_summary_snvs='Summarize on gene level',dropdown_summary_cnvs='Use exact CNVs')
    result <- generate_variants_from_locations(data,input)
    expect_equal(nrow(result[[1]]),2L)
    expect_setequal(result[[1]]$Variant,c('PHF6','KDM6A'))
    expect_equal(result[[3]],0L)
  }
})
