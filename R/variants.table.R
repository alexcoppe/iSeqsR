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


##' Merges two variants.tables.
##' 
##' Merges two variants tables obtained from read.variants.tsv function. It uses chromosome, genomic position, reference, alternative allele as a key.
##' In the case of a variant present in both tables it uses columns for the first table.
##' 
##' @param t1 a variants table obtained by read.variants.tsv
##' @param t2 a variants table obtained by read.variants.tsv
##' @return a data frame with variants from both t1 and t2. If a variant is present in both tables it keeps column for the first table.
merge2.variant.tables <- function (t1,t2) {
  t1 <- mutect.table %>% unique()
  t2 <- mutect2.table %>% unique()
  i <- intersect(t1$key, t2$key)
  t1.diff <- setdiff(t1$key,i)
  t2.diff <- setdiff(t2$key,i)
  
  common.variants <-  t1 %>% dplyr::filter(key %in% i)
  t2.variants <-  t2 %>% dplyr::filter(key %in% t2.diff)
  t1.variants <- t1 %>% dplyr::filter(key %in% t1.diff)
  merged <- rbind(common.variants, t2.variants, t1.variants)
  merged
}

##' Merges variants.tables.
##' 
##' Merges variants tables obtained from read.variants.tsv function. It uses chromosome, genomic position, reference, alternative allele as a key.
##' In the case of a variant present in both tables it uses columns for the first table.
##' 
##' @param variants.list a  list of variants table obtained by read.variants.tsv
##' @return a data frame with variants obtained by merging all variant tables in the variants.list parameter; uses chrosome, position, alternative allele and reference allele as key.
##'
##' @export
variant.tables.merge <- function(variants.list) {
  all.variants <- Reduce(merge2.variant.tables, variants.list)
  all.variants
}




build.co.mut.data <- function (all.variants, min.occurence=2, feature="impact") {
  columns.to.select <- c("gene_name", "Patient", feature)
  i <- match(columns.to.select, names(all.variants))
  co.mut.data <- all.variants %>% dplyr::select(i)
  
  if (feature == "impact" | (feature != "impact" & feature != "annotation")) {
    fill <- co.mut.data %>% dplyr::group_by(gene_name,Patient) %>% dplyr::summarise(fill = paste(unique(impact), collapse = " & ") )
  }
  if (feature == "annotation") {
    fill <- co.mut.data %>% dplyr::group_by(gene_name,Patient) %>% dplyr::summarise(fill = paste(unique(annotation), collapse = " & ") )
  }
  co.mut.data.with.fill.col <- merge(co.mut.data, fill, by = c("gene_name", "Patient")) %>% dplyr::select(gene_name, Patient, fill) %>% unique()
  
  coocurrent.genes <- co.mut.data %>% dplyr::group_by(gene_name) %>% dplyr::summarise(n = dplyr::n_distinct(Patient)) %>% dplyr::filter(n >= min.occurence)
  filtered.co.mut.data <- merge(co.mut.data.with.fill.col, coocurrent.genes, by="gene_name") %>% dplyr::arrange(n) %>% unique()
  
  sorted.genes <- dplyr::arrange(filtered.co.mut.data, n)$gene_name %>% unique
  filtered.co.mut.data$gene_name <- factor(filtered.co.mut.data$gene_name, levels = sorted.genes)

  #The following block of code is needed by cases in which there are patients with no mutation
  all.patients <- all.variants$Patient %>% unique()
  patients.with.mutations <- filtered.co.mut.data$Patient %>% unique()
  patients.without.mutations <- setdiff(all.patients, patients.with.mutations)
  if (length(patients.without.mutations) != 0)  {
    mutated.genes <- filtered.co.mut.data$gene_name %>% unique()
    not.present.gene.patient.pairs <- expand.grid(patients.without.mutations, mutated.genes)
    colnames(not.present.gene.patient.pairs) <- c("Patient", "gene_name")
    #not.present.gene.patient.pairs$impact <- "NA"
    not.present.gene.patient.pairs$fill <- "NA"
    
    i <- match(not.present.gene.patient.pairs$gene_name, filtered.co.mut.data$gene_name)
    n <- filtered.co.mut.data$n[i]
    not.present.gene.patient.pairs$n <- n
    
    all.data <- rbind(filtered.co.mut.data, not.present.gene.patient.pairs)
    all.data[ all.data == "NA" ] <- NA
    all.data
  } else {
    filtered.co.mut.data
  }
}


chroma <- c(60)
luminance <- c(70)
hue <- c(360,0)


##' Produces a very basic coMut plot.
##' 
##' Given a variants table data frame, builds a very basic coMut plot (https://www.broadinstitute.org/blog/visualizing-cancer-genome).
##' The plot shows mutation type at gene level. It also shows if a gene in a sample is hit by mutation with different types.
##' 
##' @param variants.table A tbl_df obtained with read.variants.tsv function or a data frame with impact, Patient and gene_name columns
##' @param min.occurence Minimun number of patient in which the gene should be mutated to be included in the plot
##' @param legend.position The position where to draw the legend, "none" for no legend
##' @return A ggplot
##'
##' @export
coMut.plot.main <- function(all.variants, min.occurence=2, legend.position="right", feature="impact") {

  data <- build.co.mut.data(all.variants, min.occurence = min.occurence, feature=feature)
    
  legend <- ifelse(legend.position=="none", "none", legend.position)
  
  theme.for.plot <- ggplot2::theme(
    axis.ticks=ggplot2::element_blank(),
    panel.grid.major=ggplot2::element_blank(),
    panel.background=ggplot2::element_rect(fill='#EEEEEE',colour='#EEEEEE'),
    axis.text.x=ggplot2::element_text(angle=45, hjust=1),
    axis.title.x=ggplot2::element_blank(),
    axis.title.y=ggplot2::element_blank(),
    axis.text.y=ggplot2::element_text(face = "italic"),
    legend.position = legend
  )

  ggplot2::ggplot(data, ggplot2::aes(x = Patient, y = gene_name,  fill=factor(fill))) +
    ggplot2::geom_tile(colour = "#EEEEEE") + theme.for.plot +
    ggplot2::scale_fill_discrete(name="Mutation Type", na.value="transparent",  h = hue, l=luminance, c=chroma )
}


##' ##' Produces the right plot in a coMut plot: a barchart with the number of patients in which a gene.
##' 
coMut.plot.hit.genes <- function(all.variants, min.occurence=2, legend.position="right", feature="impact") {
  
  data <- build.co.mut.data(all.variants, min.occurence = min.occurence, feature = feature)
  
  t <- ggplot2::theme(axis.line=ggplot2::element_blank(),
                      axis.text.y=ggplot2::element_blank(),
                      axis.ticks.y=ggplot2::element_blank(),
                      axis.title.x=ggplot2::element_blank(),
                      axis.title.y=ggplot2::element_blank(),
                      legend.position="right",
                      panel.background=ggplot2::element_blank(),
                      panel.border=ggplot2::element_blank(),
                      panel.grid.major=ggplot2::element_blank(),
                      panel.grid.minor=ggplot2::element_blank(),
                      plot.background=ggplot2::element_blank()
  )
  
  p2 <- ggplot2::ggplot(data, ggplot2::aes(gene_name, fill=fill)) + ggplot2::geom_bar() + ggplot2::coord_flip() +  t + 
      ggplot2::scale_fill_discrete(name="Mutation Type", na.value="transparent", h = hue, l=luminance, c=chroma)
  p2
}

##' ##' Produces the right plot in a coMut plot: a barchart with the number of patients in which a gene.
##' 
coMut.plot.mutations.in.sample <- function(all.variants, min.occurence=2, feature="impact") {
  
  data <- build.co.mut.data(all.variants, min.occurence = min.occurence, feature = feature)
  
  t <- ggplot2::theme(axis.line=ggplot2::element_blank(),
                      axis.text.x=ggplot2::element_blank(),
                      axis.ticks.x=ggplot2::element_blank(),
                      axis.title.y=ggplot2::element_blank(),
                      axis.title.x=ggplot2::element_blank(),
                      legend.position="none",
                      panel.background=ggplot2::element_blank(),
                      panel.border=ggplot2::element_blank(),
                      panel.grid.major=ggplot2::element_blank(),
                      panel.grid.minor=ggplot2::element_blank(),
                      plot.background=ggplot2::element_blank()
  )
  
  #p <- ggplot2::ggplot(filtered.co.mut.data, ggplot2::aes(Patient, fill=fill)) + ggplot2::geom_bar() + t + 
   #     ggplot2::scale_fill_discrete(name="Mutation Type")
  
  p <- ggplot2::ggplot(data, ggplot2::aes(Patient, fill=fill)) + ggplot2::geom_bar() +  t + 
    ggplot2::scale_fill_discrete(name="Mutation Type", na.value="transparent",  h = hue, l= luminance, c=chroma)
  p
}


##' Produces a coMut plot.
##' 
##' Given a variants table data frame, builds a very basic coMut plot (https://www.broadinstitute.org/blog/visualizing-cancer-genome).
##' The plot shows mutation type at gene level and on the right a bar chart with the number of patients in which a gene is.
##' It also shows if a gene in a sample is hit by mutation with different types.
##' 
##' @param variants.table A tbl_df obtained with read.variants.tsv function or a data frame with impact, Patient and gene_name columns
##' @param min.occurence Minimun number of patient in which the gene should be mutated to be included in the plot
##' @return A ggplot
##'
##' @export
coMut.plot <- function(all.variants, min.occurence=2, feature="impact") {
  upper.plot <- coMut.plot.mutations.in.sample(all.variants, min.occurence = min.occurence, feature = feature)
  main.plot <- coMut.plot.main(all.variants, legend.position = "none", min.occurence = min.occurence, feature = feature)
  right.plot <- coMut.plot.hit.genes(all.variants, min.occurence = min.occurence,  feature = feature)

  g.main <- ggplot2::ggplotGrob(main.plot)
  g.right <- ggplot2::ggplotGrob(right.plot)
  g.upper <- ggplot2::ggplotGrob(upper.plot)
  
  maxwidth = grid::unit.pmax(g.upper$widths[1:3],
                             g.main$widths[1:3])
  
  g.upper$widths[1:3] <- maxwidth
  g.main$widths[1:3] <- maxwidth
  
  maxheight <- grid::unit.pmax(g.main$heights[1:5],
                               g.right$heights[1:5])
  
  g.main$heights[1:5] <- maxheight
  g.right$heights[1:5] <- maxheight
  
  blankPanel <- grid::grid.rect(gp=grid::gpar(col="white"))
  
  gridExtra::grid.arrange(g.upper, blankPanel, g.main, g.right, ncol=2, heights=c(1,5), widths=c(2,1))
}






