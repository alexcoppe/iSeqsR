##' Reads a tab separated file with variants information
##' 
##' Reads a tab separated file with variants information produced by vcf_2_table.py script from iSeqs suite https://github.com/alexcoppe/iSeqs
##' Can load data from any tsv file, it only requires a column called "Patient"
##' 
##' 
##' @param path The path to the directory containing your coverage statistics file
##' @param filename The name of the file produced by vcf_2_table.py
##' @return A dplyr tbl_df with variants information
##'
##' @export
read.variants.tsv <- function (filename="variants.txt", path="") {
  variants.data <- dplyr::tbl_df(read.table(file.path(path, filename), header = T, sep = "\t", stringsAsFactors = F, row.names = NULL)) %>% unique()
  variants.data <- dplyr::mutate(variants.data, key=paste(chr,pos,ref,alt,sep=","))
  variants.data$Patient <- as.character(variants.data$Patient)
  variants.data
}


##' Summarises the number of variants per patient
##' 
##' Given variants data obtained with read.variants.tsv, returns a data frame with the number of variats found in each patient
##' The returned data.frame has 3 columns: total number of variants, number of variants with MODERATE predicted impact and 
##' number of variants with HIGH predicted impact for each patient.
##' 
##' 
##' @param variants.table A tbl_df obtained with read.variants.tsv function or a data frame with impact, Patient and key columns. Key column should be a column with the following string: chr,pos,alt,ref
##' @return A data frame with one line for each patient. It has 3 columns: total number of variants, number of variants with MODERATE impact and number of variants with HIGH impact
##'
##' @export
summarise.variants.by.patient <- function (variants.table) {
  data.with.selected.columns <- variants.table %>% dplyr::select(key, Patient, impact)
  
  table.with.all.variants <- data.with.selected.columns %>% dplyr::group_by(Patient) %>% dplyr::summarise(number.of.variants = dplyr::n_distinct(key))
  table.with.moderate.variants <- data.with.selected.columns %>% dplyr::filter(impact == "MODERATE")  %>% dplyr::group_by(Patient) %>% dplyr::summarise(moderate.variants=dplyr::n_distinct(key))
  table.with.high.variants <- data.with.selected.columns %>% dplyr::filter(impact == "HIGH")  %>% dplyr::group_by(Patient) %>% dplyr::summarise(high.variants=dplyr::n_distinct(key))
  list.of.dfs <- list(table.with.all.variants, table.with.moderate.variants, table.with.high.variants)
  
  merged.df <- Reduce(function(...) merge(..., by="Patient", all=T), list.of.dfs)
  merged.df[is.na(merged.df)] <- 0
  merged.df
}

##' Summarises the number of variants by impact
##' 
##' Given variants data obtained with read.variants.tsv, returns a data frame with the number of variats per impact
##' The returned data.frame has column for each impact specified by the impacts parameter
##' 
##' @param variants.table A tbl_df obtained with read.variants.tsv function or a data frame with impact, Patient and key columns. Key column should be a column with the following string: chr,pos,alt,ref
##' @param impacts a character vector with impacts that should be summarised
##' @return A named vector with the number of variants for each impact
##'
##' @export
summarise.variants.by.impact <- function (variants.table, impacts=c("MODERATE", "HIGH")) {
  data.with.selected.columns <- variants.table %>% dplyr::select(key,impact) 
  data.with.selected.columns %>% unique() %>% nrow()
  
  variants.by.impact <- sapply(impacts, function(x) {
    dplyr::filter(data.with.selected.columns, impact==x) %>% unique() %>% nrow()
  }
  )
  variants.by.impact
}

