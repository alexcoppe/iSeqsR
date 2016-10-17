mutation.group <- new.env(hash = T)
mutation.group$"C>T" <- "C>T/G>A" 
mutation.group$"G>A" <- "C>T/G>A"
mutation.group$"C>G" <- "C>G/G>C"
mutation.group$"G>C" <- "C>G/G>C"
mutation.group$"C>A" <- "C>A/G>T"
mutation.group$"G>T" <- "C>A/G>T"
mutation.group$"T>A" <- "T>A/A>T"
mutation.group$"A>T" <- "T>A/A>T"
mutation.group$"T>C" <- "T>C/A>G"
mutation.group$"A>G" <- "T>C/A>G"
mutation.group$"T>G" <- "T>G/A>C"
mutation.group$"A>C" <- "T>G/A>C"



build.mutation.categories.table <- function (variants.table) {
  mutation.pattern.data <- dplyr::select(variants.table, chr, pos, ref, alt, Patient) %>% unique()
  only.point.mutations <- dplyr::filter(mutation.pattern.data, nchar(ref) == 1 & nchar(alt) == 1) %>% mutate(mut=paste(ref, alt, sep=">"))
  only.point.mutations <- only.point.mutations %>% mutate(mut.type=unlist(mget(mut, mutation.group) ))
  
  mutations.in.patient <- only.point.mutations %>% group_by(Patient) %>% summarise(mutations.in.patient=length(Patient))
  mutation.pattern.data <- only.point.mutations %>% group_by(Patient, mut.type) %>% dplyr::summarise(number.of.mutations=length(mut.type), fraction.of.mutation=length(mut.type)/  unlist(mutations.in.patient[mutations.in.patient$Patient == Patient[1],2]) )
  mutation.pattern.data
}

choice.axis.text.size <- function (samples.number) {
  if(samples.number <= 10){
    axis.text.size <- 28
  } else if(samples.number > 10 & samples.number <= 20) {
    axis.text.size <- 20
  } else {
    axis.text.size <- 15
  }
  axis.text.size
}


##' Build a barchart with mutational spectrum per patient
##' 
##' Build a barchart with mutational spectrum of six transition and transversion categories for each patient.
##' 
##' @param variants.table A tbl_df obtained with read.variants.tsv function or a data frame with chr, pos, ref, alt and Patient columns.
##' @return A ggplot object
##'
##' @export
mutations.spectrum.barchart <- function (variants.table) {
  mutation.categories.table <- build.mutation.categories.table(variants.table)
  
  #This functions needs refactoring, should create a singe ggplot theme
  p <- ggplot2::ggplot(mutation.categories.table, ggplot2::aes(Patient, fill=mut.type)) + ggplot2::geom_bar(ggplot2::aes(weight=fraction.of.mutation), width=0.95 ) +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = 'white', colour = 'white')) +
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::ggtitle("Mutation spectrum of six transition (Ti) and transversion (Tv) categories") +
    ggplot2::theme(plot.title=ggplot2::element_text(colour = "#444444", face = "bold")) +
    ggplot2::xlab("Patient") +
    ggplot2::theme(axis.title.x = ggplot2::element_text(colour = "#666666")) +
    ggplot2::ylab("Transition/Transversion Frequency") + 
    ggplot2::theme(axis.title.y = ggplot2::element_text(colour = "#666666")) +
    ggplot2::theme(axis.text=ggplot2::element_text(size=12, angle=45), axis.title=ggplot2::element_text(size=14, color="#444444"))
  p
}

##' Build a barchart with mutational spectrum per patient
##' 
##' Build a barchart with mutational spectrum of six transition and transversion categories for each patient.
##' 
##' @param variants.table A tbl_df obtained with read.variants.tsv function or a data frame with chr, pos, ref, alt and Patient columns.
##' @return A ggplot object
##'
##' @export
mutations.spectrum.barchart2 <- function (variants.table) {
  mutation.categories.table <- build.mutation.categories.table(variants.table)
  samples.number <- dplyr::select(variants.table,Patient) %>% unique() %>% nrow()
  axis.text.size <- choice.axis.text.size(samples.number)
  axis.title.size <- 30
  legend.text.size <- 22
  p <- ggplot2::ggplot(mutation.categories.table, ggplot2::aes(Patient, fill=mut.type), ylim=c(0,1)) + 
    ggplot2::geom_bar(ggplot2::aes(weight=fraction.of.mutation), width=0.8, position = "dodge" ) + 
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = 'white', colour = '#BBBBBB')) +
    ggplot2::scale_y_continuous(expand = c(0,0), minor_breaks = seq(0, 1, 0.05) ) +
    ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = "#DDDDDD")) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_line(colour = "#AAAAAA", linetype = "dotted")) +
    ggplot2::coord_cartesian(ylim=c(0,1)) +
    ggplot2::xlab("Patient") + 
    ggplot2::theme(axis.title.x = ggplot2::element_text(colour = "#666666")) +
    ggplot2::ylab("Transition/Transversion Frequency") + 
    ggplot2::theme(axis.title.y = ggplot2::element_text(colour = "#666666")) +
    ggplot2::theme(axis.text=ggplot2::element_text(size=axis.text.size), axis.text.x=ggplot2::element_text(angle= 45, hjust= 1), axis.title=ggplot2::element_text(size=axis.title.size, color="#999999")) +
    ggplot2::theme(plot.title=ggplot2::element_text(colour = "#444444", face = "bold")) +
    ggplot2::scale_fill_discrete(name="") +
    ggplot2::theme(legend.text = ggplot2::element_text(colour="#444444", size = legend.text.size, face = "bold"))
  p
}







