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


##' Plots a barchart of number of variants per patient
##' 
##' Given variants data obtained with read.variants.tsv, plots a barchart with the number of variants found in each patient
##' 
##' 
##' @param variants.table A tbl_df obtained with read.variants.tsv function or a data frame with impact, Patient and key columns. Key column should be a column with the following string: chr,pos,alt,ref
##' @return A barchart with the number of variants in each patient.
##'
##' @export
patients.by.number.of.variants.plot <- function (variants.table) {
  variants.by.patient <- summarise.variants.by.patient(variants.table)
  patients.by.number.of.variants.plot <- sorted.by.mutations <- variants.by.patient %>% dplyr::arrange(desc(number.of.variants))
  sorted.by.mutations$Patient <- factor(sorted.by.mutations$Patient, levels = sorted.by.mutations$Patient)
  condition <- c("a")
  ggplot2::ggplot(data=sorted.by.mutations, ggplot2::aes(x=Patient, y=number.of.variants, fill=condition)) + ggplot2::geom_bar(stat="identity") +
    ggplot2::theme(axis.title.x = ggplot2::element_text(face="bold", colour="#333333", size=18)) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(angle = 0, hjust = 1, size=13,color="#333333")) +
    ggplot2::ylab("Number of variants") + ggplot2::xlab("Patient")
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


##' Summarises genes by the number of patients in which are found mutated
##' 
##' Given variants data obtained with read.variants.tsv, returns a data frame with the number of patients in which a gene is found mutated.
##' 
##' @param variants.table A tbl_df obtained with read.variants.tsv function or a data frame with impact, Patient and gene_name columns.
##' @return A tbl_df with 2 columns: gene name, number of patients in which the 
##'
##' @export
variants.recurrence.by.gene.name <- function(variants.table) {
  recurrent.variants <- variants.table %>% dplyr::group_by(gene_name) %>% summarise(in.n.samples=dplyr::n_distinct(Patient))
  ordered.recurrent.variants <- arrange(recurrent.variants, desc(in.n.samples))
  ordered.recurrent.variants
}


