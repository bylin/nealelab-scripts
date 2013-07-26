library("ggplot2")
library("grid")

addBarGraph <- function(plot) {
  #default colors, see http://colorbrewer2.org/
  data <- read.table(file="/Users/Student/introns.txt", head=TRUE, sep="\t")
  datasets <- c("Genome", "All", "First", "First 20-49 Kb", "First 50-99 Kb", "First 100+ Kb", "Mod 20-49 Kb", "Mod 50-99 Kb", "Mod 100+ Kb", "High 20-49 Kb", "High 50-99 Kb", "High 100+ Kb")
  data$Dataset <- factor(data$Dataset, levels=datasets, ordered=TRUE)
  plot <- plot + geom_bar(data=data, aes(x=Dataset, y=Percentage, fill=Classification, width=0.62),stat="identity")
  #classifColors <- c("#FF7F00","#984EA3", "#377EB8", "#4DAF4A",  "#E41A1C", "#FFFF33")
  classifColors <- c("#E41A1C" ,"#4DAF4A" ,  "#377EB8", "#984EA3" ,"#FFFF33", "#FF7F00")
  classifs <- c("LTR Gypsy", "LTR Copia","LINE", "Other retrotransposon","DNA transposon","Simple repeats","Uncategorized")
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
  plot <- plot + xlab("Intron set") + ylab("Percentage of dataset")
}

main <- function() {
  theme_set(theme_gray(base_size=20))
  plot <- ggplot()
  plot <- addBarGraph(plot)
  plot <- addGCGraph(plot)
  plot <- formatPlot(plot)
  plot
}
main()