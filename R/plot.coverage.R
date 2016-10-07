
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
##' coverage.barchart(path = directory, pattern="-coverage-hist.txt$", intervals = c(0,10,20,30,50))
##'
##' @export
coverage.barchart <- function (path=".", intervals=c(0,5,10,20,30, 50), pattern="-coverage-hist.txt$") {
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
    categories <- zoo::rollapply(intervals, 2, function(x) {
      paste(ifelse(x[1] == 0, x[1], x[1] + 1), x[2], sep="-")
    } )
    categories[length(categories)] <- paste(">", intervals[length(intervals) - 1], sep="")
    categories.t <- c(categories.t, categories)
    
  }
  
  categories.t <- factor(categories.t, levels=categories, ordered=T)
  
  d <- data.frame(patients=patient.t, bases=bases.t, coverage=categories.t)
  
  labs <- gsub(pattern, "", files)
  axis.text.size <- 14
  axis.title.size <- 18
  legend.text.size <- 18
  
  theme.for.plot <- ggplot2::theme(
    axis.title.y = ggplot2::element_text(colour = "#666666"),
    axis.title.x = ggplot2::element_text(colour = "#666666"),
    axis.text=ggplot2::element_text(size=axis.text.size),
    axis.text.x=ggplot2::element_text(angle=45, hjust=1),
    axis.title=ggplot2::element_text(size=axis.title.size, color="#343123"),
    legend.text = ggplot2::element_text(colour="#444444", size = legend.text.size, face = "bold")
  )
  
  p <- ggplot2::ggplot(d, ggplot2::aes(patients, fill=coverage)) + 
    ggplot2::geom_bar(ggplot2::aes(weight=bases / exome.length, order=rev(categories.t))) +
    ggplot2::xlab("Patient") + 
    ggplot2::ylab("% of targeted bases covered") +
    ggplot2::scale_x_discrete(labels=labs) + theme.for.plot
  
  p
}



##' Plot the cumulative distribution of the fraction of targeted bases covered by n reads
##' 
##' Plot the cumulative distribution of the fraction of targeted bases covered by n reads. It needs coverage statistic files generated with the command:
##' bedtools coverage -hist -abam sample.bam -b target_regions.bed | grep ^all > sample-coverage-hist.txt
##' It is based on the code from this the post in Getting Things Done in Genetics: http://www.gettinggeneticsdone.com/2014/03/visualize-coverage-exome-targeted-ngs-bedtools.html
##' 
##' 
##' @param path The path to the directory containing your coverage statistics file
##' @param pattern A regular expression pattern. Only file names which match the regular expression will be returned
##' @examples
##' fpath <- system.file("extdata", "sample2-coverage-hist.txt", package="iSeqsR")
##' directory <- dirname(fpath)
##' coverage.cumulative.distribution.plot(path = directory, pattern="-coverage-hist.txt$")
##'
##' @export
coverage.cumulative.distribution.plot <- function(path=".", pattern="-coverage-hist.txt$"){
  files <- list.files(path = path, pattern=pattern)
  labs <- gsub(pattern, "", files)
  cov <- list()
  cov_cumul <- list()
  for (i in 1:length(files)) {
    cov[[i]] <- read.table(file.path(path, files[i]))
    cov_cumul[[i]] <- 1-cumsum(cov[[i]][,5])
  }
  

  palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(11,"Paired"))(length(cov))
  
  plot(cov[[1]][2:401, 2], cov_cumul[[1]][1:400], type='n', xlab="Depth", ylab="Fraction of capture targeted bases with coverage \u2265 depth", ylim=c(0,1.0), main="Targeted Regions Coverage")
  abline(v = 10, col = "gray60")
  abline(v = 20, col = "gray60")
  abline(v = 30, col = "gray60")
  abline(v = 40, col = "gray60")
  abline(v = 50, col = "gray60")
  abline(v = 100, col = "gray60")
  abline(h = 0.50, col = "gray60")
  axis(1, at=c(10,20,30,40,50,100), labels=c(10,20,30,40,50,100))
  axis(2, at=c(0.90), labels=c(0.90))
  axis(2, at=c(0.50), labels=c(0.50))
  
  for (i in 1:length(cov)) points(cov[[i]][2:401, 2], cov_cumul[[i]][1:400], type='l', lwd=3, col=palette[i])
  
  legend("topright", legend=labs, col=palette, lty=1, lwd=4)
  
}



##' Calculates mean coverage for all samples in a directory
##' 
##' Calculates mean coverage for all samples in a given directory. It needs coverage statistic files generated with the command:
##' bedtools coverage -hist -abam sample.bam -b target_regions.bed | grep ^all > sample-coverage-hist.txt
##'
##' 
##' @param path The path to the directory containing your coverage statistics file
##' @param pattern A regular expression pattern. Only file names which match the regular expression will be returned
##' @examples
##' fpath <- system.file("extdata", "sample2-coverage-hist.txt", package="iSeqsR")
##' directory <- dirname(fpath)
##' samples.mean.coverage(path = directory, pattern="-coverage-hist.txt$")
##'
##' @export
##' @import dplyr
samples.mean.coverage <- function (path=".", pattern="-coverage-hist.txt$") {
  files <- list.files(path = path, pattern=pattern)
  names <- c("coverage", "bases", "exome", "perc")
  
  all.data <- data.frame(coverage=numeric(0), bases=numeric(0), exome=numeric(0), perc=numeric(0), patient=character(0))
  for (i in 1:length(files)){
    f <- files[i]
    data <- read.table(file.path(path, files[i]), sep="\t", header=F)
    cols <- ncol(data)
    data <- data[,2:cols]
    colnames(data) <- names
    patient <- gsub(pattern, "", as.character(f))
    data$patient <- patient
    exome.length <- data[["exome"]][1]
    all.data <- rbind(data, all.data)
  }
  all.data$mult <- all.data$coverage * all.data$bases
  all.data
  patients.mean.coverage <- dplyr::group_by(all.data, patient) %>% dplyr::summarise(  mean.coverage=sum(as.numeric(mult)) / exome[1]   )
  patients.mean.coverage
}


