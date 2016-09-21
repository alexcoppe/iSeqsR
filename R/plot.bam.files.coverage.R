
##' Plots a barchart of bam files coverage.
##' 
##' Plot a barchart of bam files coverage. It needs coverage statistic files generated with the command:
##' bedtools coverage -hist -abam sample.bam -b target_regions.bed | grep ^all > sample-coverage-hist.txt
##' 
##' 
##' @param path The path to the directory containing your coverage statistics file
##' @param intervals A numeric vector with coverage intervals
##' @param pattern A regular expression pattern. Only file names which match the regular expression will be returned
##' @return A ggplot object to be printed
##' @examples
##' fpath <- system.file("extdata", "sample2-coverage-hist.txt", package="iSeqsR")
##' directory <- dirname(fpath)
##' bam.files.coverage.plot(path = directory, pattern="-coverage-hist.txt$", intervals = c(0,10,20,30,50))
##'
##' @export
bam.files.coverage.plot <- function (path=".", intervals=c(0,5,10,20,30, 50), pattern="-coverage-hist.txt$") {
  files <- list.files(path = path, pattern=pattern)
  intervals <- c(intervals, 1000)
  names <- c("coverage", "bases", "exome", "perc")
  patient.t <- c()
  bases.t <- c()
  categories.t <- c()
  exome.length <- 0
  i <- 1
  for (i in 1:length(files)){
    f <- files[i]
    data <- read.table(file.path(path, f), sep="\t", header=F)
    cols <- ncol(data)
    data <- data[,2:cols]
    colnames(data) <- names
    bases.in.interval.of.coverage <- zoo::rollapply(intervals, 2, function(x) {
      sum(subset(data, coverage >= x[1] & coverage < x[2])$bases)
    } )
    bases.t <- c(bases.t, bases.in.interval.of.coverage)
    patient <- as.character(f)
    patients <- rep(patient, length(intervals) - 1)
    exome.length <- data[["exome"]][1]
    patient.t <- c(patient.t, patients)
    categories <- rollapply(intervals, 2, function(x) {
      paste(x[1], x[2], sep="-")
    } )
    categories[length(categories)] <- paste(">", intervals[length(intervals) - 1], sep="")
    categories.t <- c(categories.t, categories)
    
  }
  
  categories.t <- factor(categories.t, levels=categories, ordered=T)
  
  d <- data.frame(patients=patient.t, bases=bases.t, coverage=categories.t)
  
  labs <- gsub("-coverage-hist.txt", "", files)
  axis.text.size <- 26
  axis.title.size <- 30
  legend.text.size <- 22
  p <- ggplot2::ggplot(d, ggplot2::aes(patients, fill=coverage)) + ggplot2::geom_bar(ggplot2::aes(weight=bases / exome.length, order=rev(categories.t))) +
    ggplot2::xlab("Patient") + 
    ggplot2::ylab("% of targeted bases covered") +
    ggplot2::theme(axis.title.y = ggplot2::element_text(colour = "#666666") ) +
    ggplot2::theme(axis.title.x = ggplot2::element_text(colour = "#666666")) +
    ggplot2::theme(axis.text=ggplot2::element_text(size=axis.text.size), axis.title=ggplot2::element_text(size=axis.title.size, color="#343123")) +
    ggplot2::theme(legend.text = ggplot2::element_text(colour="#444444", size = legend.text.size, face = "bold")) +
    ggplot2::scale_x_discrete(labels=labs)
  
  p
}
