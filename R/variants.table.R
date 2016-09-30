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


##' Plots a matirx plot of mutated genes vs patient
##' 
##' Given variants data obtained with read.variants.tsv, plots a matrix plot of mutated genes vs patient.
##' 
##' @param variants.table A tbl_df obtained with read.variants.tsv function or a data frame with Patient and gene_name columns.
##' @param cluster_cols boolean value determining if columns should be clustered 
##' @param cluster_rows boolean value determining if rows should be clustered
##' @param min.recurrence minimum number of patients in which should be mutated to be used 
##' @return A matrix plot with mutated genes vs patient
##'
##' @export
gene.vs.patient.matrix.plot <- function (variants.table,  cluster_cols = F,  cluster_rows = F, min.recurrence = 2, ...) {
  #Create a data table with 3 columns: gene symbol, sample in which gene is mutated, a column of 1s
  data.for.heatmap <- dplyr::select(variants.table, gene_name, Patient) %>% dplyr::mutate(Patient=as.character(Patient)) %>% dplyr::mutate(in.sample=1)
  #data.for.heatmap <- dplyr::select(filtered.data.with.symbols.no.na, symbol, sample, ratio) %>% mutate(sample=as.character(sample))
  data.for.heatmap <- unique(data.for.heatmap)
  data.for.heatmap
  
  #Create the matrix to be used for plotting the heatmap. Rows are gene symbols, columns are samples
  heatmap.data <- reshape2::dcast(data.for.heatmap, gene_name ~ Patient)
  heatmap.data[is.na(heatmap.data)] <- 0
  heatmap.data <- tbl_df(heatmap.data)
  
  heatmap.data$recurrency <- rowSums(heatmap.data[,2:ncol(heatmap.data)])
  
  #Select only genes mutated in at least 2 samples
  recurrent.genes <- dplyr::filter(heatmap.data, recurrency >= min.recurrence)
  
  #Oder rows of data frame by recurrencys
  recurrent.genes <- arrange(recurrent.genes, desc(recurrency))
  
  #Build the matrix for heatmap plotting
  heatmap.matrix <- as.matrix(recurrent.genes[,2:(ncol(recurrent.genes)-1)])
  rownames(heatmap.matrix) <- recurrent.genes$gene_name
  
  #Build the heatmap plot
  pheatmap::pheatmap(heatmap.matrix, legend=F, cluster_cols= cluster_cols, cluster_rows =  cluster_rows, color=colorRampPalette(c("white",  "#6699FF"))(50), 
                     border_color="#777777", display_numbers = FALSE, ...)
}



