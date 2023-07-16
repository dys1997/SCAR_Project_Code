#.libPaths("/home/dengys/R/x86_64-pc-linux-gnu-library/4.2")
library(Seurat)
#library(SeuratData)
library(optparse)
op_list <- list(make_option(c("-i", "--inputfile"), type = "character", default = 5,action = "store", help = "Singlecell rds file",metavar="character"),
                make_option(c("-o", "--outputfile"), type = "character", default = F, action = "store", help = "outputfile ",metavar="character")
                )

parser <- OptionParser(option_list = op_list)
opt = parse_args(parser)
print(opt$inputfile)
print(opt$outputfile)


object<-readRDS(opt$inputfile)

subset<-SplitObject(object = object, split.by = "celltype")

dir.create(opt$outputfile)
setwd(opt$outputfile)
celltype<-unique(as.character(object$celltype))
j<-0
for(i in 1:length(celltype)){
  tmp<-subset[[i]]
  name<-gsub(" ","_",celltype[i])
  name<-gsub("/","_or_",name)
  #print(paste0(SP,"----",name,"------cell number:",ncol(tmp)))
  if(ncol(tmp)>100){
    tmp.data<-tmp@assays$RNA@data
    tmp.data.t<-t(as.matrix(tmp.data))
    j<-j+1
    dir.create(paste0("CT",j))
    setwd(paste0("CT",j))
    write.csv(tmp.data.t,"data_t_exp.csv",quote=F)
    write.csv(tmp.data.t,paste0(name,"_data_t_exp.csv"),quote=F)
    n<- ncol(tmp) %/% 10
    set.seed(777)
    sam<-sample(1:ncol(tmp),n)
    tmp.data.t.tmp<-tmp.data.t[sam,]
    write.csv(tmp.data.t.tmp,"data_t_exp_10.csv",quote=F)
    setwd("./../")
  }
}
write.csv(tmp.data.t.tmp,"good_job.csv",quote=F)