library("ggplot2")
library("grid")
library("scales")

addBarGraph <- function(plot) {
  #default colors, see http://colorbrewer2.org/
  data <- read.table(file="/Users/Student/introns.txt", head=TRUE, sep="\t")
  datasets <- c("Genome", "All", "First", "First 20-49 Kb", "First 50-99 Kb", "First 100+ Kb", "Distal", "Distal 20-49 Kb", "Distal 50-99 Kb", "Distal 100+ Kb")
  data$Dataset <- factor(data$Dataset, levels=datasets, ordered=TRUE)
  plot <- plot + geom_bar(data=data, aes(x=Dataset, y=Percentage, fill=Classification, width=0.62),stat="identity")
  classifColors <- c( "#2B8CBE","#FEC44F",  "#045A8D",  "#41AB5D","#74A9CF" ,"#78C679"  , "#D7301F","#FF7F00" )
  classifs <- c("Gypsy LTR", "Copia LTR", "Other LTR",  "LINE", "Other retrotransposon","DNA transposon","Uncategorized", "Simple repeats")
  plot <- plot + scale_fill_manual(name = "", values=classifColors, breaks=classifs)
  return(plot)
}


addGCGraph <- function(plot) {
  data <- read.table(file="/Users/Student/gc.txt", head=TRUE, sep="\t")
  plot <- plot + geom_line(data=data,aes(x=Dataset,y=GC.Content,group=1, linetype="GC Content"),colour="black", size=1)
  plot <- plot + geom_point(data=data, aes(x=Dataset,y=GC.Content), colour="black")
  plot <- plot + scale_linetype_manual(values=1, name="")
  return(plot)
}


formatPlot <- function(plot) {
  plot <- plot + theme(axis.title.x = element_text(vjust=0,face="bold"), axis.title.y = element_text(vjust=0,face="bold"),
                       axis.text.y = element_text(colour="black"), axis.text.x = element_text(colour="black", angle=45, hjust=1),
                       plot.margin = unit(c(1,0,1,1), "cm"))
  plot <- plot + scale_y_continuous(labels=percent, breaks=c(0, 0.15, 0.30, 0.45, 0.60, 0.75))
  #plot <- plot + guides(fill = guide_legend(nrow = 2)) + theme(legend.position="bottom")
  plot <- plot + xlab("Intron set") + ylab("Percentage of dataset")
}

main <- function() {
  #theme_set(theme_gray(base_size=20))
  plot <- ggplot()
  plot <- addBarGraph(plot)
  plot <- addGCGraph(plot)
  plot <- formatPlot(plot)
  plot
  ggsave(plot, file="intronRepeats.tiff", dpi=400)
}
main()