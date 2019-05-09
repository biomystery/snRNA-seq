# Mouse data compared to druncseq (Habib2017, mouse dataset)

## pariwise correlation of the average expression for the genes
## in each cell-type signature defined by our data and cell types defined by Drunc-seq

source('./libs.R')
# load signature genes for cell type defined by us  -----------------------


anno.sigGenes <- read.table('../data/new_analysis/mouse/intron/mouse.marker.genes.txt',
                            stringsAsFactors = F,header = F)$V1
length(anno.sigGenes)#4090

# load druncseq avg expression data  --------------------------------------

# expression

dat.exp.log <- fread("../data/mouse/Mouse_Processed_GTEx_Data.DGE.log-UMI-Counts.txt")
dat.exp.log <- dat.exp.log[GENE%in% anno.sigGenes]
dat.exp.log.t <- dcast(melt(dat.exp.log,id.vars = "GENE"),variable ~ GENE)

# cell id
dat <- read.csv(file = "../data/mouse/Table3_Data_Info_mouse.csv",skip = 27,stringsAsFactors = F)[,-8]
colnames(dat) <-sub("X.","",colnames(dat))

# merge data
dat.exp.log.m <- merge(dat.exp.log.t,dat[,c("Cell.ID","Cluster.name")],by.y = "Cell.ID", by.x = "variable")
colnames(dat.exp.log.m)[1]<- "Cell.ID"

# change to long
dat.exp.log.l <- melt(dat.exp.log.m,id.vars = c("Cell.ID","Cluster.name"))
head(dat.exp.log.l)

# avg
dat.exp.log.avg <- dat.exp.log.l %>% group_by(variable,Cluster.name) %>%
  summarise(mean.exp=mean(value,na.rm=T))

dat.exp.log.avg <- dat.exp.log.l %>% group_by(variable,Cluster.name) %>%
  summarise(mean.exp=transfunc(value))

transfunc <- function(x) log2(sum(2^x-1)/length(x)+1)

dat.exp.log.avg.2<- dat.exp.log.avg;
dat.exp.log.avg.2$Cluster.name <- sub("[0-9]","",dat.exp.log.avg$Cluster.name)
dat.exp.log.avg.2<- dat.exp.log.avg.2 %>% group_by(variable,Cluster.name) %>%
  summarise(mean.exp=mean(mean.exp,na.rm=T))


# spread
dat.exp.log.avg.s <- dat.exp.log.avg%>% spread(key = Cluster.name,value = mean.exp)
druncseq.mouse <- dat.exp.log.avg.s[,-ncol(dat.exp.log.avg.s)]

druncseq.mouse <- as.data.frame(druncseq.mouse)
rownames(druncseq.mouse) <- druncseq.mouse$variable; druncseq.mouse$variable <- NULL

druncseq.mouse.2 <- dat.exp.log.avg.2%>% spread(key = Cluster.name,value = mean.exp)
druncseq.mouse.2<- druncseq.mouse.2[,-ncol(druncseq.mouse.2)]

druncseq.mouse.2 <- as.data.frame(druncseq.mouse.2)
rownames(druncseq.mouse.2) <- druncseq.mouse.2$variable; druncseq.mouse.2$variable <- NULL


# load our expression data  -----------------------------------------------

# expression
dat.exp <- lapply(1:14,function(x){
  fn <- paste0("../data/new_analysis/mouse/intron/average_expression/AFB_merged_intronic_clusteringcluster",x,"averageumi.txt")
  read.table(fn,header = T,stringsAsFactors = F,col.names = paste0("C",x))
})

dat.exp <- do.call(cbind,dat.exp)

our.mouse <- dat.exp[ rownames(druncseq.mouse),]
setDF(our.mouse)

all.equal(rownames(our.mouse),rownames(druncseq.mouse))


# Our cluster anno --------------------------------------------------------
our.mouse.clust <- read.csv('../data/new_analysis/mouse/intron/mouse.cluster.anno.txt',
                              header = F,stringsAsFactors = F,col.names = c("clust","cell.type"))
rownames(our.mouse.clust)<- paste0("C",our.mouse.clust$clust)

# correlation matrix  -----------------------------------------------------
our.mouse <- log2(our.mouse+1)
calcCormat <- function(a,b) apply(a,2,function(x) apply(b,2,function(y) cor(x,y,method = "spearman")))
#calcCormat <- function(a,b) apply(a,2,function(x) apply(b,2,function(y) cor(x,y)))
cor.mat <- calcCormat(our.mouse,druncseq.mouse)
cor.mat.2 <- calcCormat(our.mouse,druncseq.mouse.2)

pd <- (cor.mat%>% tbl_df()%>% rownames_to_column("Druncseq")%>% gather(ours,value,- "Druncseq"))
ggplot(pd,aes(Druncseq,ours))+
  geom_tile(aes(fill=value))+
  geom_text(aes(label=round(value,1)))+
  scale_fill_gradient(low = "white",high = "red")+theme_bw()



if(T){
    pdf(file =   "../figs/fig1d.pdf",width = 7,height = 7)
  cols=colorRampPalette(c(                        "white",
                          rgb(203/255,72/255,85/255)))(11)
  pheatmap(cor.mat,scale = "none",cluster_cols = T,cluster_rows = T,colorRampPalette(brewer.pal(9,"Reds"))(11),
           cellwidth = 12,cellheight = 12,breaks = seq(0.3,0.85,by = 0.05),
           main = "Mouse pair-wised Spearman's correlation \nof avg expression between our data (row) and \n drunc-seq (column) from our signature genes")
  pheatmap(cor.mat.2,scale = "none",cluster_cols = T,cluster_rows = T,colorRampPalette(brewer.pal(9,"Reds"))(11),
           cellwidth = 12,cellheight = 12,breaks = seq(0.3,0.85,by = 0.05),
           main = "Mouse pair-wised Spearman's correlation \nof avg expression between our data (row) and \n drunc-seq (column) from our signature genes (reduced)")
  colnames(cor.mat) <- our.mouse.clust[colnames(cor.mat),'cell.type']
  colnames(cor.mat.2) <- our.mouse.clust[colnames(cor.mat.2),'cell.type']
  pheatmap(cor.mat,scale = "none",cluster_cols = T,cluster_rows = T,colorRampPalette(brewer.pal(9,"Reds"))(11),
           cellwidth = 12,cellheight = 12,breaks = seq(0.3,0.85,by = 0.05),
           main = "Mouse pair-wised Spearman's correlation \nof avg expression between our data (row) and \n drunc-seq (column) from our signature genes")
  pheatmap(cor.mat.2,scale = "none",cluster_cols = T,cluster_rows = T,colorRampPalette(brewer.pal(9,"Reds"))(11),
           cellwidth = 12,cellheight = 12,breaks = seq(0.3,0.85,by = 0.05),
           main = "Mouse pair-wised Spearman's correlation \nof avg expression between our data (row) and \n drunc-seq (column) from our signature genes (reduced)")


  dev.off()

}
range(cor.mat.2)
range(cor.mat)
