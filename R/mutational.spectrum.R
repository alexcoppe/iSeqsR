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


##' Build a barchart with mutational spectrum per patient
##' 
##' Build a barchart with mutational spectrum of six transition and transversion categories for each patient.
##' 
##' @param variants.table A tbl_df obtained with read.variants.tsv function or a data frame with chr, pos, ref, alt and Patient columns.
##' @return A ggplot object
##'
##' @export
mutations.spectrum.barchart <- function (variants.table) {
  mutation.pattern.data <- dplyr::select(variants.table, chr, pos, ref, alt, Patient) %>% unique()
  only.point.mutations <- dplyr::filter(mutation.pattern.data, nchar(ref) == 1 & nchar(alt) == 1) %>% mutate(mut=paste(ref, alt, sep=">"))
  only.point.mutations <- only.point.mutations %>% mutate(mut.type=unlist(mget(mut, mutation.group) ))
  
  mutations.in.patient <- only.point.mutations %>% group_by(Patient) %>% summarise(mutations.in.patient=length(Patient))
  mutation.pattern.data <- only.point.mutations %>% group_by(Patient, mut.type) %>% dplyr::summarise(a=length(mut.type), b=length(mut.type)/  unlist(mutations.in.patient[mutations.in.patient$Patient == Patient[1],2]) )
  mutation.pattern.data
  
  library(ggplot2)
  p <- ggplot(mutation.pattern.data, aes(Patient, fill=mut.type)) + geom_bar(aes(weight=b), width=0.95 ) +
    theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
    scale_y_continuous(expand = c(0,0)) +
    ggtitle("Mutation spectrum of six transition (Ti) and transversion (Tv) categories") +
    theme(plot.title=element_text(colour = "#444444", face = "bold")) +
    xlab("Patient") +
    theme(axis.title.x = element_text(colour = "#666666")) +
    ylab("Transition/Transversion Frequency") + 
    theme(axis.title.y = element_text(colour = "#666666")) +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14, color="#444444"))
  p
}
