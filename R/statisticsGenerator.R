generateStatistics <- function(dataset) {
  patient_ids <- unique(dataset$Patient.ID)
  rowCounts <- numeric(length(patient_ids))
  uniqueCloneCounts <- numeric(length(patient_ids))
  i <- 1
  for (id in patient_ids) {
    patientData <- dataset[dataset$Patient.ID == id, , drop = FALSE]
    rowCounts[i] <- nrow(patientData)
    uniqueCloneCounts[i] <- length(unique(patientData$Clone))
    i <- i + 1
  }

  rowCountFrequency <- table(rowCounts)
  uniqueCloneCountFrequency <- table(uniqueCloneCounts)

  variantPatientMap <- tapply(
    dataset$Patient.ID,
    dataset$Variant,
    function(x) length(unique(x))
  )
  variantPatientFrequencyDistribution <- table(variantPatientMap)

  variantFrequency <- table(dataset$Variant)
  sortedVariants <- sort(variantFrequency, decreasing = TRUE)

  top10Variants <- head(sortedVariants, 10)
  top10VariantsPercentage <- (top10Variants / sum(variantFrequency)) * 100

  chromosomeFrequency <- table(dataset$Chromosome)
  if (length(chromosomeFrequency) == 0) {
    chromosomePercentage <- NULL
  } else {
    chromosomePercentage <- (chromosomeFrequency / nrow(dataset)) * 100
  }

  percentage_snv_indel <- mean(
    dataset$Variant.Type == "SNV/INDEL",
    na.rm = TRUE
  ) * 100
  percentage_NA <- mean(dataset$Variant.Type == "NA", na.rm = TRUE) * 100
  percentage_cnv <- 100 - percentage_snv_indel - percentage_NA

  return(list(
    RowCountFrequency = rowCountFrequency, # Number of Variants per ID
    UniqueCloneCountFrequency = uniqueCloneCountFrequency,
    UniqueVariants = length(variantPatientMap),
    Top10Variants = top10Variants, # 10 most frequent variant labels
    Top10VariantsPercentage = top10VariantsPercentage,
    ChromosomePercentage = chromosomePercentage,
    SNVINDELPercentage = percentage_snv_indel, # Percentage of local variants
    CNVPercentage = percentage_cnv, # Percentage of copy number alterations
    NAPercentage = percentage_NA # Percentage of undescribed alterations
  ))
}
