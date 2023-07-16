library(Seurat)
data <- Read10X('path_to_cellranger_result')
SCAR_Atlas_0316 <- CreateSeuratObject(counts = data,project = 'SCAR_Atlas_0316',min.cells = 3,min.features = 200)
####QC

grep(pattern = "^MT-",rownames(SCAR_Atlas_0316),value = T)
SCAR_Atlas_0316[["percent.MT"]] <- PercentageFeatureSet(SCAR_Atlas_0316,pattern = "^MT-")
SCAR_Atlas_0316 <- subset(SCAR_Atlas_0316,subset = nFeature_RNA>200 & 
                            nFeature_RNA <6000 & 
                            percent.MT < 10)
SCAR_Atlas_0316
#
###STEP5
SCAR_Atlas_0316 <- NormalizeData(SCAR_Atlas_0316, normalization.method = "LogNormalize",scale.factor = 10000)
SCAR_Atlas_0316 <- FindVariableFeatures(SCAR_Atlas_0316,selection.method = "vst",nfeatures = 2000)

all.SCAR_Atlas_0316.gene <- rownames(SCAR_Atlas_0316)
SCAR_Atlas_0316 <- ScaleData(SCAR_Atlas_0316,features = all.SCAR_Atlas_0316.gene)
SCAR_Atlas_0316 <- RunPCA(SCAR_Atlas_0316 ,features = VariableFeatures(object = SCAR_Atlas_0316))
SCAR_Atlas_0316 <- JackStraw(SCAR_Atlas_0316,num.replicate = 100,dims = 30)
SCAR_Atlas_0316 <- ScoreJackStraw(SCAR_Atlas_0316,dims = 1:30)
JackStrawPlot(SCAR_Atlas_0316,dims = 1:30)

nPCs <- 30
SCAR_Atlas_0316 <- FindNeighbors(SCAR_Atlas_0316,dims = 1:nPCs)
SCAR_Atlas_0316 <- FindClusters(SCAR_Atlas_0316, resolution = 1)
head(Idents(SCAR_Atlas_0316))
SCAR_Atlas_0316 <- RunUMAP(SCAR_Atlas_0316, dims = 1:nPCs)
SCAR_Atlas_0316 <- RunTSNE(SCAR_Atlas_0316, dims = 1:nPCs)
DimPlot(SCAR_Atlas_0316, reduction = "tsne", label = T, label.size = 4,pt.size = 1.5,repel = T)
###STEP7
saveRDS(SCAR_Atlas_0316, file = "path_to_save_rds")

#/mnt/alamo01/users/dengys/lmr/SCAR_Atlas_0316/SCAR_Atlas_0316lmr.rds
SCAR_Atlas_0316 <- readRDS("path_to_save_rds")

#'MZB1','LYZ','FCGR3B','CD79A','FOXP3','KLRF1','SLC4A10','CD8A','CD4'
SCAR_Atlas_0316features <- c('MZB1','LYZ','FCGR3B','CD79A','FOXP3','KLRF1','SLC4A10','CD8A','CD4')

DotPlot(object = SCAR_Atlas_0316, 
        cols = c("white", "red"),
        features= SCAR_Atlas_0316features, 
        cluster.idents=T) + 
  theme(axis.text.x = element_text("control",size = 10,face="bold",angle = 90)
  )

new.cluster.ids <- c("CD8_T_cell-0",   
                     "CD4_T_cell-1",    
                     "CD4_T_cell-2",  
                     "Natural_killer_cell-3",            
                     " Regulatory T cell-4", 
                     "CD8_T_cell-5",
                     "Natural_killer_cell-6",
                     "Plasma-7",
                     "Myeloid-8",
                     "Mucosal_associated_invariant_T_cell-9",
                     'CD4_T_cell-10',
                     'CD8_T_cell-11',
                     'CD8_T_cell-12',
                     'B_cell-13',
                     'B_cell-14',
                     "CD8_T_cell-15",    
                     "Mucosal_associated_invariant_T_cell-16",
                     'CD4_T_cell-17',
                     'B_cell-18',
                     'Myeloid-19',
                     'CD4_T_cell-20',
                     'Myeloid-21',
                     "CD4_T_cell-22",    
                     "Neutrophil-23",  
                     'Mucosal_associated_invariant_T_cell-24',
                     'CD8_T_cell-25',
                     'Neutrophil-26',
                     ' Regulatory T cell-27',
                     'Plasma-28',
                     'CD4_T_cell-29',
                     "Plasma-30",    
                     "Myeloid-31",
                     'CD8_T_cell-32',
                     'Mucosal_associated_invariant_T_cell-33',
                     'Natural_killer_cell-34',
                     'B_cell-35')


names(new.cluster.ids) <- levels(SCAR_Atlas_0316)
SCAR_Atlas_0316<- RenameIdents(SCAR_Atlas_0316, new.cluster.ids)
Idents(SCAR_Atlas_0316)
SCAR_Atlas_0316[['celltype']] <- Idents(SCAR_Atlas_0316)
unique(SCAR_Atlas_0316$celltype)
DimPlot(SCAR_Atlas_0316, reduction = "tsne", label = T, label.size = 4,pt.size = 1.5,repel = T)


saveRDS(SCAR_Atlas_0316, file = "path_to_save_rds")

SCAR_Atlas_0316$Annotation     


