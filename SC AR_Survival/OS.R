rm(list=ls())
library(cgdsr) 
library(DT)
library(survival)
library(dplyr)
library (survminer)
mycgds <- CGDS("http://www.cbioportal.org/")
test(mycgds)
all_TCGA_studies <-getCancerStudies(mycgds)
DT::datatable(all_TCGA_studies)
blca_tcga <- "blca_tcga"
getCaseLists(mycgds,blca_tcga)[,c(1,2)]
getGeneticProfiles(mycgds,	blca_tcga)[,c(1,2)]
choose_genes_human <- read.table("C:/Users/86131/Desktop/human.txt")
gene1 <- choose_genes_human[c(1:1000),]
gene2 <- choose_genes_human[c(1001:2000),]
gene3 <- choose_genes_human[c(2001:3000),]
gene4 <- choose_genes_human[c(3001:4000),]
gene5 <- choose_genes_human[c(4001:5000),]
gene6 <- choose_genes_human[c(5001:6000),]
gene7 <- choose_genes_human[c(6001:7000),]
gene8<- choose_genes_human[c(7001:8000),]
gene9 <- choose_genes_human[c(8001:9000),]
gene10 <- choose_genes_human[c(9001:10000),]
gene11 <- choose_genes_human[c(10001:11000),]
gene12 <- choose_genes_human[c(11001:12000),]
gene13 <- choose_genes_human[c(12001:13000),]
gene14 <- choose_genes_human[c(13001:14000),]
gene15 <- choose_genes_human[c(14001:15000),]
gene16 <- choose_genes_human[c(15001:16000),]
gene17 <- choose_genes_human[c(16001:17000),]
gene18 <- choose_genes_human[c(17001:18000),]
gene19 <- choose_genes_human[c(18001:19433),]
GENE <- list(gene1,gene2,gene3,gene4,gene5,gene6,gene7,gene8,gene9,gene10,gene11,gene12,gene13,gene14,gene15,gene16,gene17,gene18,gene19)
#
getCaseLists(mycgds,blca_tcga)[,c(1,2)]
getGeneticProfiles(mycgds,blca_tcga)[,c(1,2)]
mycaselist <-getCaseLists(mycgds,blca_tcga)[6,1]
mygeneticprofile <- getGeneticProfiles(mycgds,blca_tcga)[7,1]

for(i in 6:19) {
expr <- getProfileData(mycgds,GENE[[i]],
                       mygeneticprofile,mycaselist)
#View(expr)
mycaselist <- getCaseLists(mycgds,blca_tcga)[1,1] 
myclinicaldata <- getClinicalData(mycgds,mycaselist)
#View(myclinicaldata)
choose_columns=c('OS_STATUS','OS_MONTHS')
choose_clinicaldata <- myclinicaldata[,choose_columns]
dat=cbind(choose_clinicaldata[,c('OS_STATUS','OS_MONTHS')],expr[rownames(choose_clinicaldata),])
dat=dat[dat$OS_MONTHS > 0,]
write.csv(dat,"C:/Users/86131/Desktop/COAD2.csv")
library(survival)
library(dplyr)
library (survminer)
dat <- read.csv("C:/Users/86131/Desktop/COAD2.csv")
removeRowsAllNa  <- function(x){x[apply(x, 1, function(y) any(!is.na(y))),]}
removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}
dat1 <- removeColsAllNa(dat)
dat2 <- removeRowsAllNa(dat1)
dat3 <- na.omit(dat2)
library(tidyr)
dat4 <- separate(dat3,col = OS_STATUS,into = c("OS_STATUS","B"),sep = ':')
dat <- dat4[,-3]
dat$OS_STATUS <- as.numeric(dat$OS_STATUS)

for(gene in colnames(dat)[4:ncol(dat)]){
  group <- ifelse(dat[[gene]] > median(dat[[gene]]),
                  "high","low")
  if(length(table(group))==1) next
  fit <- survfit(Surv(OS_MONTHS,OS_STATUS) ~ group,data = dat)
  p=ggsurvplot(fit,palette = c("#FF3030", "#0000FF"),
               #title = "Lung cancer (n=40)",
               plot.btitle = element_text(hjust = 0.5),
               font.main = 18,
               xlab ="OS(months)", 
               usurv.median.line = "hv", 
               #risk.table = TRUE,
               #risk.table.col = "strata",
               pval = TRUE,
               legend.title = "",
               legend = c(0.15,0.1), 
               legend.labs = c("A High expression", "A Low  expression") ,
               break.x.by = 25,
               break.y.by = 0.25,
               xlim = c(0,200),
               ggtheme =theme_survminer())  
  png(file=paste0("F:/SCAR/SCAR_survial/BL/",gene,"_BLCA.png"),width = 600,height = 600,units = "px")
  print(p)
  dev.off()
}
 }