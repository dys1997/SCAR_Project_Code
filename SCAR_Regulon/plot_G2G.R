library(optparse)
library(igraph)
#library(Cairo)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(getopt)
op_list <- list(make_option(c("-i", "--inputfile"), type = "character", default = 5,action = "store", help = "Singlecell rds file",metavar="character"))

parser <- OptionParser(option_list = op_list)
opt = parse_args(parser)
print(opt$inputfile)
setwd(dir = opt$inputfile)

dir.create("iGraph")
dir.create("GO")

exp<-read.table(paste(opt$inputfile,"/expr_mat.adjacencies.tsv",sep=""),header = T,sep = "\t")
colnames(exp)<-c("regulatoryGene","targetGene","weight")
#j<-0
set.seed(777)
reguG<-unique(as.character(exp$regulatoryGene))
if(length(reguG) > 32){
  reguG<-reguG[1:32]
  exp<-subset(exp, subset = regulatoryGene %in% reguG)
}
for (j in 1:length(reguG)) {
  print(j)
  i<-reguG[j]
  tmp<-subset(exp,subset = regulatoryGene %in% i)
  if(nrow(tmp)> 32){
    object<-tmp[1:32,]
  }else{
    object<-tmp
  }
  
  relations <- data.frame(from=object$regulatoryGene,
                          to=object$targetGene,
                          #paris=f$pair,
                          weight=object$weight*10
                          #mode=f$mode
  )

  g <- graph_from_data_frame(relations, directed=TRUE,vertices=unique(union(object$targetGene,object$regulatoryGene)))
  #x<-length(unique(V(g)))

  #print(x)
  V(g)[as.character(unique(object$targetGene))]$color<-"#6CB92D"
    V(g)[as.character(unique(object$regulatoryGene))]$color<-"#3DBBC3"
      #V(g)$count<-count_source[V(g),2] 
    V(g)[as.character(unique(object$targetGene))]$size<-8
    V(g)[as.character(unique(object$regulatoryGene))]$size <- 15
    V(g)[as.character(unique(object$targetGene))]$label.cex<-0.8
    V(g)[as.character(unique(object$targetGene))]$label.color<-"black"
      V(g)[as.character(unique(object$targetGene))]$label.font<-3
      V(g)[as.character(unique(object$regulatoryGene))]$label.cex <- 1.5
      V(g)[as.character(unique(object$regulatoryGene))]$label.color <- "red"
        E(g)$arrow.size <- 0
        #E(g)$edge.color <- "gray80"
        E(g)$width <- E(g)$weight*0.01
        

cat("GOOD HERE!")
        
        l<-layout_nicely(g)

        
        png(paste0("./iGraph/iGraph",j,".png"),width = 700,height = 700,units = "px",res = 100)
        plot(g,layout=l,
             #vertex.label.color="red",
             #vertex.label.cex=0.6,
             vertex.frame.color="white",
             vertex.color=V(g)$color,
             vertex.shape="circle",
             vertex.size=V(g)$size,
             edge.width=E(g)$width,
             edge.arrow.width=0,
             edge.arrow.size=0,
             arrow.mode="0",
             edge.color="gray50",
             #mark.groups=unique(V(g)),
             #mark.border=0
        )
        #plot(p)
        dev.off()
        pdf(paste0("./iGraph/iGraph",j,".pdf"),width = 6,height = 6)
        plot(g,layout=l,
             #vertex.label.color="red",
             #vertex.label.cex=0.6,
             vertex.frame.color="white",
             vertex.color=V(g)$color,
             vertex.shape="circle",
             vertex.size=V(g)$size,
             edge.width=E(g)$width,
             edge.arrow.width=0,
             edge.arrow.size=0,
             arrow.mode="0",
             edge.color="gray50",
             #mark.groups=unique(V(g)),
             #mark.border=0
        )
        #plot(p)
        dev.off()
        
        
        if(nrow(tmp) > 100){
          gene<- as.character(tmp$targetGene)
          mygene<-select(org.Hs.eg.db,columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL",keys=as.character(gene))
          mygene<-mygene$ENTREZID
          mygene<-na.omit(mygene)
          ego<-enrichGO(OrgDb = "org.Hs.eg.db",gene=mygene,ont="BP",pvalueCutoff=0.05,readable=TRUE)
          #ego$TF<-i
          #write.table(ego,paste0("./GO/GO",j,".txt"),sep="\t",row.names = F,quote = F)
          ego<-data.frame(ego,stringsAsFactors = F)
          if(nrow(ego) > 0){
            ego$TF<-i
            write.table(ego,paste0("./GO/GO",j,".tsv"),sep="\t",row.names = F,quote = F)
            
            if(nrow(ego) > 32){
              selected_go<-ego[1:32,]
            }else {
              selected_go<-ego
            }
            
            tt <- factor(selected_go$Description, levels=rev(unique(selected_go$Description)))
            pp = ggplot(selected_go,aes(-1*log10(p.adjust),tt))
            pbubble = pp +geom_point(aes(size=Count,color=-1*log10(p.adjust)))+
              scale_colour_gradient(low="blue",high="red") +
              theme_bw() +scale_size(range=c(2, 10)) +
              labs(title = paste(unique(tmp$regulatoryGene),"_GO",sep = ""))+
              ylab("")+xlim(min(-1*log10(selected_go$p.adjust))*0.9, max(-1*log10(selected_go$p.adjust))*1.1)+
              theme(axis.text.y = element_text(size = 20, family = "Helvetica", color = "black",  angle = 0),
                    axis.text.x = element_text(size = 15, family = "Helvetica", color = "black", angle = 0))+
              theme(title = element_text(color="black",size = 30,family = "Helvetica"))+
              theme(legend.title = element_text(color="black", size=20, family = "Helvetica"))+
              theme(legend.text = element_text(color="azure4", size = 15,  family = "Helvetica"))+
              theme(axis.title = element_text(color="black", size=16, family = "Helvetica"))+
              theme(legend.key.size=unit(1.1,'cm'))
            #ggsave(file="HVG_shared_GO.pdf", plot=pbubble, width=20, height=10)
            pdf(paste0("./GO/GO",j,".pdf"),width=20, height=20)
            plot(pbubble)
            dev.off()
            
            png(paste0("./GO/GO",j,".png"),width = 700,height = 700,units = "px",res = 40)
            plot(pbubble)
            dev.off()
          }
        }
        
}


printer = file(paste(opt$inputfile,"/goodjob.txt",sep=""),"w")
writeLines("goodjob!!",con=printer,sep=" ")
#writeLines("The same line.",con=printer)
close(printer)
