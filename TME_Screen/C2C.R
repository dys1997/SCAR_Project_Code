#DizEED chat第二版
library(Seurat)
library(SeuratData)
library(CellChat)
library(patchwork)
library(ggplot2)
library(Matrix)
library(cowplot)
library(optparse)
op_list <- list(make_option(c("-i", "--inputfile"), type = "character", default = 5,action = "store", help = "Singlecell rds file",metavar="character"),
                make_option(c("-o", "--outputfile"), type = "character", default = F, action = "store", help = "outputfile ",metavar="character")
                )
parser <- OptionParser(option_list = op_list)
opt = parse_args(parser)



setwd(dir = opt$outputfile)
SCAR_Atlas_0489<- readRDS(file = opt$inputfile)

####
CellChatDB <- CellChatDB.human
setwd(dir = opt$outputfile)
dir.create('CT')
dir.create('HP')
dir.create('TSV')
####若cluster idents包含0，则运行
plusC_cluster_ids = paste0("C", levels(SCAR_Atlas_0489))
names(plusC_cluster_ids) <- levels(SCAR_Atlas_0489)
SCAR_Atlas_0489 <- RenameIdents(SCAR_Atlas_0489, plusC_cluster_ids)
#
data.input <- GetAssayData(SCAR_Atlas_0489,assay = 'RNA',slot = 'data')
labels <- Idents(SCAR_Atlas_0489)
identity <- data.frame(group = labels, row.names = names(labels))
SCAR_Atlas_0489chat <- createCellChat(object = data.input, group.by = "labels", assay = "RNA")
SCAR_Atlas_0489chat <- addMeta(SCAR_Atlas_0489chat, meta = identity, meta.name = "labels")
SCAR_Atlas_0489chat <- setIdent(SCAR_Atlas_0489chat, ident.use = "labels") # set "labels" as default cell identity
showDatabaseCategory(CellChatDB)
SCAR_Atlas_0489chat@DB <- CellChatDB
info <- dplyr::glimpse(CellChatDB$interaction)

SCAR_Atlas_0489chat <- subsetData(SCAR_Atlas_0489chat)
SCAR_Atlas_0489chat@data.signaling
future::plan("multisession", workers = 4)
SCAR_Atlas_0489chat <- identifyOverExpressedGenes(SCAR_Atlas_0489chat)
SCAR_Atlas_0489chat <- identifyOverExpressedInteractions(SCAR_Atlas_0489chat)
SCAR_Atlas_0489chat <- projectData(SCAR_Atlas_0489chat, PPI.human)
SCAR_Atlas_0489chat <- computeCommunProb(SCAR_Atlas_0489chat,raw.use = T)
SCAR_Atlas_0489chat <- filterCommunication(SCAR_Atlas_0489chat, min.cells = 10)
SCAR_Atlas_0489chat <- computeCommunProbPathway(SCAR_Atlas_0489chat)
SCAR_Atlas_0489chat <- aggregateNet(SCAR_Atlas_0489chat)



#vcell2cell 可视化
setwd(dir = paste(opt$outputfile,'/CT/',sep=""))
groupSize <- as.numeric(table(SCAR_Atlas_0489chat@idents))
#p4
png(file = paste('CTall','.png',sep = ''),
    width = 2200,height = 2200,units = 'px',res = 150)
netVisual_circle(SCAR_Atlas_0489chat@net$count,vertex.weight = groupSize,weight.scale = T,label.edge = F,
                 title.name = 'Interaction weights')
dev.off()
#each celltype
mat <- SCAR_Atlas_0489chat@net$weight
head(Idents(SCAR_Atlas_0489))
test <- as.data.frame(table(SCAR_Atlas_0489@active.ident))
test.list <- as.character(test$Var1)
#each cluster的png格式-cell2cell
setwd(dir = paste(opt$outputfile,'/CT/',sep=""))
for (
  i in 1:nrow(mat)) {
  png(file = paste('CT',i,'.png',sep = ''),width = 800,height = 800,units = 'px',res = 80)
  #png(file = paste(rownames(mat)[i],'.png',sep = ''),width = 800,height = 800,units = 'px',res = 80)
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  dev.off()
}

##heatmap of pathway
setwd( dir = paste(opt$outputfile,'/HP/',sep=""))
SCAR_Atlas_0489chat@netP[["pathways"]]
pathway.list <- SCAR_Atlas_0489chat@netP[["pathways"]]
pathway.list

for (i in 1:length(pathway.list)) {
  png(file = paste('HP',i,'.png',sep = ''),width = 1500,height = 1000,units = 'px',res = 120)
  #png(file = paste(pathway.list[i],' signaling network.png',sep = ''),width = 1500,height = 1000,units = 'px',res = 120)
  print(netVisual_heatmap(SCAR_Atlas_0489chat, signaling = pathway.list[i], color.heatmap = "Reds"))
  dev.off()
}
###csv
setwd(dir = paste(opt$outputfile,'/TSV/',sep=""))
slot.name = 'netP'
df.net <- subsetCommunication(SCAR_Atlas_0489chat)
write.csv(df.net, file = paste(opt$outputfile,'/TSV/CellChat.tsv',sep=""),quote = F,sep = ',')
#saveRDS(object = SCAR_Atlas_0489chat,file = 'D:/YAN/XIEHE/SCAR/C2C/SCAR_Atlas_0489/SCAR_Atlas_0489.rds')




printer = file(paste(opt$outputfile,"/good_job.txt",sep=""),"w")
writeLines("goodjob!!",con=printer,sep=" ")
#writeLines("The same line.",con=printer)
close(printer)