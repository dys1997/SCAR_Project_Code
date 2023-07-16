library(ggplot2)
library(shiny)
library(ggsci)
library(RColorBrewer)
library(Seurat)
library(optparse)
library(openxlsx)

op_list <- list(make_option(c("-i", "--inputfile1"), type = "character", default = F,action = "store", help = "Singlecell rds file",metavar="character"),
                make_option(c("-I", "--inputfile2"), type = "character", default = F,action = "store", help = "Pathway matrix",metavar="character"),
                make_option(c("-p", "--plot_dir"), type = "character", default = F,action = "store", help = "Output dir of plots",metavar="character"),
                make_option(c("-t", "--table"), type = "character", default = F, action = "store", help = "An csv file with pathway information",metavar="character"))

parser <- OptionParser(option_list = op_list)
opt = parse_args(parser)
rds = readRDS(opt$inputfile1)
pathway = read.csv(opt$inputfile2,row.names = 1)
pathway = pathway[which(rownames(pathway)!=""),]
getpalette =colorRampPalette(brewer.pal(12,"Set3"))
for (i in rownames(pathway)){
  p= Seurat::VlnPlot(rds,features = i)
  p = p+scale_fill_manual(values = getpalette(length(unique(rds$celltype))))+
    xlab("")+ylab("")+
    theme(plot.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1,size = 8),
          legend.text = element_text(size = 10))
  name = i
  name<-gsub(" ","_",name)
  name<-gsub("/","_or_",name)
  ggsave(paste(opt$plot_dir,paste(name,".png",sep = ""),sep = "/"),p,height = 7,width = 10,dpi=100)
}
pathway.list = rownames(pathway)
for (name in length(pathway.list)){
    pathway.list[name]<-gsub(" ","_",pathway.list[name])
    pathway.list[name]<-gsub("/","_or_",pathway.list[name])
}
Pathway = as.data.frame(pathway.list)
colnames(Pathway) = "Pathway"
write.csv(Pathway,paste(opt$table,"Pathway.csv",sep = "/"),row.names = FALSE)
