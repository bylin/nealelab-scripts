library("ggplot2")
library("grid")

addBarGraph <- function(plot) {
  #default colors, see http://colorbrewer2.org/
  data <- read.table(file="/Users/Student/introns.txt", head=TRUE, sep="\t")
  plot <- plot + geom_bar(data=data, aes(x=Dataset, y=Percentage, fill=Classification, width=0.62 ),stat="identity")
  classifColors <- c("#E41A1C","#377EB8","#984EA3","#4DAF4A", "#FF7F00")
  classifs <- c("Simple repeats", "Uncategorized", "DNA transposon", "Non-LTR retrotransposon", "LTR retrotransposon")
  plot <- plot + scale_fill_manual(name = "", values=classifColors, breaks=classifs)
  return(plot)
}


addLineGraph <- function(plot) {
  data <- read.table(file="/Users/Student/gc.txt", head=TRUE, sep="\t")
  plot <- plot + geom_line(data=data,aes(x=Dataset,y=Percentage,group=1,linetype="GC Content"),colour="black", size=1)
  plot <- plot + geom_point(data=data, aes(x=Dataset,y=Percentage))
  plot <- plot + scale_linetype_manual(values=1, name="")
  return(plot)
}

formatPlot <- function(plot) {
  plot <- plot + theme(axis.title.x = element_text(vjust=0,face="bold"), axis.title.y = element_text(vjust=0,face="bold"),
                       axis.text.y = element_text(colour="black"), axis.text.x = element_text(colour="black", angle=45, hjust=1),
                       plot.margin = unit(c(1,0,1,1), "cm"))
  plot <- plot + xlab("Intron set") + ylab("Percentage of dataset")
}

tiff("plot.tiff", height=600, width=960)
theme_set(theme_gray(base_size=20))
plot <- ggplot()
plot <- addBarGraph(plot)
plot <- addLineGraph(plot)
plot <- formatPlot(plot)
plot
dev.off()