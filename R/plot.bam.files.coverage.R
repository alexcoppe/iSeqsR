
##' Plot histogram of bam files coverage
##' 
##' @param intervals A numeric vector with coverage intervals
##' @param pattern The pattern to be matched to find coverage statistics file in the current directory
##' @export
plot.bam.files.coverage <-
function (intervals=c(0,5,10,20,30, 50), pattern="-coverage-hist.txt$") {
  files <- list.files(pattern=pattern)
  intervals <- c(intervals, 1000)
  names <- c("coverage", "bases", "exome", "perc")
  patient.t <- c()
  bases.t <- c()
  categories.t <- c()
  exome.length <- 0
  i <- 1
  for (i in 1:length(files)){
    f <- files[i]
    #print(f)
    data <- read.table(f, sep="\t", header=F)
    cols <- ncol(data)
    data <- data[,2:cols]
    colnames(data) <- names
    bases.in.interval.of.coverage <- rollapply(intervals, 2, function(x) {
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
  print(d)
  p <- ggplot(d, aes(patients, fill=coverage)) + geom_bar(aes(weight=bases / exome.length, order=rev(categories.t))) +
    xlab("Patient") + 
    ylab("% of targeted bases covered") +
    theme(axis.title.y = element_text(colour = "#666666") ) +
    theme(axis.title.x = element_text(colour = "#666666")) +
    theme(axis.text=element_text(size=axis.text.size), axis.title=element_text(size=axis.title.size, color="#343123")) +
    theme(legend.text = element_text(colour="#444444", size = legend.text.size, face = "bold")) +
    scale_x_discrete(labels=labs)
  
  p
}
