# Mouse data compared to druncseq 

## pariwise correlation of the average expression for the genes 
## in each cell-type signature defined by our data and cell types defined by Drunc-seq


# load signature genes for cell type defined by us  -----------------------

anno.cell <- read.table("./data/new_analysis/mouse/AFB_new_analysis/afb_merged_svdmeta_data")
anno.sigGenes <- lapply(1:11,function(x){
  fn <- paste0("./data/new_analysis/mouse/AFB_new_analysis/afb_merged_final.markergenes",x,".csv") 
  data.frame(genes=read.csv(fn,header = T,stringsAsFactors = F)$X,clust=x)
})
anno.sigGenes <- do.call(rbind,anno.sigGenes)
dim(anno.sigGenes)#5707 

# load druncseq avg expression data  --------------------------------------

# expression 
require(data.table)
dat.exp.log <- fread("./data/mouse/Mouse_Processed_GTEx_Data.DGE.log-UMI-Counts.txt")
dat.exp.log <- dat.exp.log[GENE%in% anno.sigGenes$genes]
dat.exp.log.t <- dcast(melt(dat.exp.log,id.vars = "GENE"),variable ~ GENE)

# cell id 
dat <- read.csv(file = "./data/mouse/Table3_Data_Info_mouse.csv",skip = 27,stringsAsFactors = F)[,-8]
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
require(data.table)
dat.exp.log <- fread("./data/new_analysis/mouse/AFB_new_analysis/afb_merged_final_umi_allgenes.txt")
our.mouse <- dat.exp.log[V1 %in% rownames(druncseq.mouse)]
setDF(our.mouse)
rownames(our.mouse) <- our.mouse$V1;our.mouse$V1 <- NULL 

our.mouse <- our.mouse[rownames(druncseq.mouse),]

all.equal(rownames(our.mouse),rownames(druncseq.mouse))


# correlation matrix  -----------------------------------------------------
our.mouse <- log2(our.mouse+1)
calcCormat <- function(a,b) apply(a,2,function(x) apply(b,2,function(y) cor(x,y,method = "spearman")))
calcCormat <- function(a,b) apply(a,2,function(x) apply(b,2,function(y) cor(x,y)))
cor.mat <- calcCormat(our.mouse,druncseq.mouse)
cor.mat.2 <- calcCormat(our.mouse,druncseq.mouse.2)

pd <- (cor.mat%>% tbl_df()%>% rownames_to_column("Druncseq")%>% gather(ours,value,- "Druncseq"))
ggplot(pd,aes(Druncseq,ours))+
  geom_tile(aes(fill=value))+
  geom_text(aes(label=round(value,1)))+
  scale_fill_gradient(low = "white",high = "red")+theme_bw()

if(T){
  pdf(file =   "fig1d.pdf",width = 7,height = 7)
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
