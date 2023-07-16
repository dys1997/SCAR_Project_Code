# .libPaths() 
library(optparse)
op_list <- list(make_option(c("-i", "--inputfile"), type = "character", default = 5,action = "store", help = "Singlecell rds file",metavar="character"),
                make_option(c("-o", "--outputfile"), type = "character", default = F, action = "store", help = "outputfile ",metavar="character"))

parser <- OptionParser(option_list = op_list)
opt = parse_args(parser)
print(opt$inputfile)
print(opt$outputfile)

library(CytoTRACE)

##marrow_10x_expr_results <- CytoTRACE(mat = marrow_10x_expr)
#phenotype为每个barcodes的
#plotCytoGenes(marrow_10x_expr_results, numOfGenes = 10)
#plotCytoTRACE(marrow_10x_expr_results, phenotype = marrow_10x_pheno, gene = "Top2a")
#ggplot2::ggsave('marrow_10x_expr_results_plotCytoTRACE.png',width = 12)
library(Seurat)
rds = readRDS(opt$inputfile)

if(length(rownames(rds@meta.data)) > 3000){
    set.seed(1)
    select = sample(c(1:length(rownames(rds@meta.data))),3000)  
    phe <- Idents(rds)[select]
    phe = as.character(phe)
    rownames(rds@meta.data)[select]
    names(phe) <- rownames(rds@meta.data)[select]
    mat_3k <- as.matrix(rds@assays$RNA@counts)[,select]
    mat_3k[1:4,1:4]
    results <- CytoTRACE(mat = mat_3k)
}else{
    phe <- Idents(rds)
    phe = as.character(phe)
    names(phe) <- rownames(rds@meta.data)
    mat_3k <- as.matrix(rds@assays$RNA@counts)
    mat_3k[1:4,1:4]
    results <- CytoTRACE(mat = mat_3k)
}


plotCytoGenes(results, numOfGenes = 10,outputDir = opt$outputfile)
plotCytoTRACE(results, phenotype = phe,outputDir = opt$outputfile)
#ggplot2::ggsave(p,'marrow_10x_expr_results_plotCytoTRACE.png',width = 12)


# Use the CytoTRACE function to read the expression matrix of single-cell data, predict the stemness of the cell through the expression matrix information, and save the prediction result as an object that CytoTRACE can recognize. The prediction process of this function can be divided into two steps. In the first step, Cytotrace will count the genes in the expression matrix and convert them into TPM or CPM, and then calculate the Pearson correlation between the normalized expression of each gene and the gene count, and get the Gene counts signature (GCS) by correlation sorting . In the second step, CytoTRACE converts the normalized expression matrix into a Markov matrix, based on which the non-negative least squares regression (NNLS) is applied to the GCS, and finally the GCS is adjusted iteratively using a diffusion process based on the probability structure of the Markov process . The resulting values are ranked and scaled between 0 and 1, representing the predicted order of cells by their relative differentiation status (0, more differentiated; 1, less differentiated). Finally, the object and the cell type information of the single cell are simultaneously passed as parameters to the plotCytoTRACE function for visual display of the cell stemness prediction results.