##' Prepares a data frame ready to be transformed to a matrix for fishplot.
##' 
##' Prepares a data frame to be transformed to a matrix for fishplot. As imput 'data' it needs a data frame
##' with COVERAGE, VAF and SAMPLE_NAME columns
##' @param data The data frame with starting data
##' @param diagnosis The name of the diagnosis in the SAMPLE_NAME column
##' @param relapse The name of the relapse in the SAMPLE_NAME column
##' @param relapse2 The name of the second relapse in the SAMPLE_NAME column (default = FALSE, second relaps). 
##' @param coverage The minimum coverage of the called variants (default = 30)
##' @param min.vaf The minimum VAF value in % (dafault = 5)
##' @param relapse2 If the data of a secondo relapse is present (default = FALSE)
##' @return A data frame with these columns: VAF_DIAGNOSIS, GENE, VAF_RELAPSE
##' 
##' 
##' @export
data.for.fishplot <- function(data, diagnosis = "C01_ES_C03_REM", relapse = "C02_RIC_C03_REM", relapse2 = FALSE, 
                              coverage = 30, min.vaf = 5){
  # Filter data by coverage >= 30
  data.filtered.by.coverage <- data %>% dplyr::filter(COVERAGE >= coverage)
  # Obtain VAF as 100%
  data.filtered.by.coverage$VAF <- data.filtered.by.coverage$VAF * 100
  # VAF and GENE data from diagnosis
  diagnosis.data <- data.filtered.by.coverage %>% dplyr::filter(SAMPLE_NAME == diagnosis) %>% 
    dplyr::select(VAF, GENE)
  diagnosis.data <- diagnosis.data %>% dplyr::group_by(GENE) %>% dplyr::summarise(VAF_DIAGNOSIS = mean(VAF))
  # VAF and GENE data from relapse
  relapse.data <- data.filtered.by.coverage %>% dplyr::filter(SAMPLE_NAME == relapse) %>% 
    dplyr::select(VAF, GENE)
  relapse.data <- relapse.data %>% dplyr::group_by(GENE) %>% dplyr::summarise(VAF_RELAPSE = mean(VAF))
  # VAF and GENE data from relapse II
  if (relapse2 != FALSE){
    relapse2.data <- data.filtered.by.coverage %>% dplyr::filter(SAMPLE_NAME == relapse2) %>% 
      dplyr::select(VAF, GENE)
    relapse2.data <- relapse2.data %>% dplyr::group_by(GENE) %>% dplyr::summarise(VAF_RELAPSE2 = mean(VAF))
  }
  # Join the samples by GENE column
  diagnosis.and.relapse.data <- dplyr::full_join(diagnosis.data, relapse.data, by = "GENE")
  if (relapse2 != FALSE){
    diagnosis.and.relapse.data <- dplyr::full_join(diagnosis.and.relapse.data, relapse2.data, by = "GENE")
  }
  
  # Change names of VAF.x and VAF.y
  if (relapse2 != FALSE) {
    colnames(diagnosis.and.relapse.data) <- c("GENE", "VAF_DIAGNOSIS", "VAF_RELAPSE", "VAF_RELAPSE2")
  } 
  else {
    colnames(diagnosis.and.relapse.data) <- c("GENE", "VAF_DIAGNOSIS", "VAF_RELAPSE")
  }
  
  # Change NA to 0
  diagnosis.and.relapse.data[is.na(diagnosis.and.relapse.data)] <- 0
  
  # Filter diagnosis.and.relapse.data with VAF > 5 in at least one sample
  filtered.diagnosis.and.relapse.data <-  diagnosis.and.relapse.data %>% dplyr::filter(VAF_DIAGNOSIS > min.vaf | VAF_RELAPSE > min.vaf)

  filtered.diagnosis.and.relapse.data
}


##' Get a matrix with VAF_DIAGNOSIS, VAF_RELAPSE and eventualy VAF_RELAPSE2 columns and GENES as row names to ben passed todbscan::dbscan
##' 
##' Starting from the data.frame from data.for.fishplot function
##' obtains a matrix with VAF_DIAGNOSIS, VAF_RELAPSE and eventually VAF_RELAPSE and GENES as row names
##' to be passed to dbscan::dbscan clustering function
##' 
##' @param data The data frame obtained from data.for.fishplot function
##' @param relapse2 If the data from a second relapse is present (default = FALSE)
##' @return A mtrix with VAF_DIAGNOSIS, VAF_RELAPSE and eventually VAR_RELAPSE2 columns and GENES as row names
##' 
##' @export
matrix.for.fishplot <- function(data, relapse2 = FALSE){
  if (relapse2 != FALSE) {
    matrix.for.fishplot <- as.matrix(data[,c(1,3,4)])
  }
  else {
    matrix.for.fishplot <- as.matrix(data[,c(1,3)])
  }
  row.names <- data$GENE
  rownames(matrix.for.fishplot) <- row.names
  matrix.for.fishplot
}


##' Build a matrix ready to be passed to createFishObject 
##' 
##' Starting with data from data.for.fishplot and clustering results from dbscan::dbscan
##' build a matrix to be used in createFishObject
##' 
##' @param fishplot.data.frame data from data.for.fishplot funtion
##' @param clustering.results Clustering results from dbscan::dbscan function
##' @param relapse2 If the data from a second relapse is present (default = FALSE)
##' @return A matrix to be used in createFishObject function
##' 
##' @export
build.fishplot.matrix <- function(fishplot.data.frame, clustering.results, relapse2 = FALSE){
  # Add cluster comlumn to fishplot.data.frame
  fishplot.data.frame$CLUSTER <- clustering.results$cluster
  fishplot.data.frame
  # Obtain median of each CLUSTER from VAF_DIAGNOSIS and VAF_RELAPSE
  # First step: group by CLUSTER
  grouped.fishplot.data.frame <- fishplot.data.frame %>% dplyr::group_by(CLUSTER)
  diagnosis.clusters.medians  <-  grouped.fishplot.data.frame %>% dplyr::summarise(VAF = median(VAF_DIAGNOSIS))
  relapse.clusters.medians    <-  grouped.fishplot.data.frame %>% dplyr::summarise(VAF = median(VAF_RELAPSE))
  if (relapse2 != FALSE) {
    relapse2.clusters.medians   <-  grouped.fishplot.data.frame %>% dplyr::summarise(VAF = median(VAF_RELAPSE2))
  }
  
  # Set fish maximum value
  max.value <- 95
  
  # Adjust the VAF column based on the max.values
  diagnosis.clusters.medians.adjusted <- diagnosis.clusters.medians$VAF / sum(diagnosis.clusters.medians$VAF) * max.value
  relapse.clusters.medians.adjusted   <- relapse.clusters.medians$VAF / sum(relapse.clusters.medians$VAF) * max.value
  if (relapse2 != FALSE){
    relapse2.clusters.medians.adjusted  <- relapse2.clusters.medians$VAF / sum(relapse2.clusters.medians$VAF) * max.value
  }
  
  # Create the central columns
  # Add the 2 central columns
  central.column1 <- rep(0, length(diagnosis.clusters.medians.adjusted))
  central.column2 <- rep(0, length(relapse.clusters.medians.adjusted))
  
  if (relapse2 != FALSE){
    central.column3 <- rep(0, length(relapse2.clusters.medians.adjusted))
    central.column4 <- rep(0, length(relapse2.clusters.medians.adjusted))
    fishplot.matrix.data <- c(diagnosis.clusters.medians.adjusted, central.column1, central.column2, 
                              relapse.clusters.medians.adjusted, relapse.clusters.medians.adjusted,
                              central.column3, central.column4,
                              relapse2.clusters.medians.adjusted, relapse2.clusters.medians.adjusted)
  }
  else {
    fishplot.matrix.data <- c(diagnosis.clusters.medians.adjusted, central.column1, central.column2, 
                              relapse.clusters.medians.adjusted, relapse.clusters.medians.adjusted)
  }
  
  fishplot.matrix <- matrix(fishplot.matrix.data, nrow = length(central.column1), byrow = F)
  
  if (relapse2 != FALSE){
    fishplot.matrix <- rbind(c(99, 1, 1, 99, 99, 1, 1, 99, 99), fishplot.matrix)
  }
  else {
    fishplot.matrix <- rbind(c(99, 1, 1, 99, 99), fishplot.matrix)
  }
  
  fishplot.matrix
}



fish.colors = c("#C4C3C3", "#F9F8F7", "#F26BFB", "#86EFFF", "#FFC386", 
                "#C3FF86", "#FFA086", "#FCFFC1", "#CDCDFF", "#0000C2",
                "#B60000", "#9B00B6", "#8C2424", "#24818C", "#8C2481",
                "#F77CBD", "#F77C7C", "#F7BA7C", "#7CF7BA", "#7C98F7",
                "#000000", "#8AA876", "#BD93B1", "#FFCDF3", "#09178E")

