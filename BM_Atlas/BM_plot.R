library(ggplot2)
library(shiny)
library(ggsci)
library("RColorBrewer")
library(Seurat)
library(openxlsx)
library(optparse)
op_list <- list(make_option(c("-i", "--inputfile"), type = "character", default = 5,action = "store", help = "Singlecell rds file",metavar="character"),
                make_option(c("-o", "--outputfile"), type = "character", default = F, action = "store", help = "outputfile ",metavar="character"))

parser <- OptionParser(option_list = op_list)
opt = parse_args(parser)
print(opt$inputfile)
print(opt$outputfile)
par(pin = c(50,30))
diz = readRDS(opt$inputfile)####
ff = length(unique(diz$celltype))
getpalette =colorRampPalette(brewer.pal(12,"Set3"))
#读取rds文件

#读取marker列表
mk_data = read.csv("/mnt/alamo01/users/dengys/Project/SCAR/scar_bm/scripts/BM_Atlas.csv")
#读取数据集中表达的基因与所有biomarker之间的交集
gene_list =intersect(mk_data$Biomarker,rownames(diz@assays$RNA@data))
for (i in gene_list){
  p= Seurat::VlnPlot(diz,features = i)
  p=p+scale_fill_manual(values = getpalette(ff))+
    xlab("")+ylab("")+
    theme(plot.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1,size = 15),
          legend.text = element_text(size = 10))
  ggsave(plot = p,filename = paste(opt$outputfile,paste(as.character(i),"_violin.png",sep = ""),sep="/"),height = 7,width = 10,dpi=100)
  # data = p$data
  # colnames(data)<- c("count","cluster")
  # p2 = ggplot(data,aes(x=cluster,y=count,fill=cluster))+
  #   geom_boxplot()+
  #   scale_fill_manual(values = getpalette(ff))+
  #   xlab("")+ylab("")+
  #   theme_classic()+
  #   theme(plot.title = element_text(size = 30),
  #         axis.text.x = element_text(angle = 45, hjust = 1,size = 15),
  #         legend.title=element_blank(),
  #         legend.text = element_text(size = 10)
  #   )
  # ggsave(plot = p2,filename = paste(opt$outputfile,paste(as.character(i),"_box.png",sep = ""),sep="/"),height = 10,width = 15)
}
data = as.data.frame(gene_list)
colnames(data)="Gene"
write.xlsx(x = data, file = paste(opt$outputfile,"Gene.xlsx",sep="/"),
           sheetName = "Sheet1", row.names = FALSE)
