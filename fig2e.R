# Human data compared to druncseq 

## pariwise correlation of the average expression for the genes 
## in each cell-type signature defined by our data and cell types defined by Drunc-seq

# load our expression data  -----------------------------------------------

# expression 
our.human <- read.table("./data/new_analysis/human/fc_hc_latest_svd_umi_upregmarkers_nodup.txt")
dim(our.human)#3508 

all.equal(rownames(our.mouse),rownames(druncseq.mouse))

# load druncseq expression data  --------------------------------------

# expression 
require(data.table)
require(Seurat)
require(tidyverse)


dat.exp <- fread("./data/human/GTEx_droncseq_hip_pcf.umi_counts.txt")
setDF(dat.exp)
dat.exp <- dat.exp%>% column_to_rownames("V1")

dim(dat.exp) # 32111 by 14963 

# seurat object 
droncseq.human <- CreateSeuratObject(raw.data = dat.exp, min.cells = 3, min.genes = 200, 
                           project = "snRNAseq")

#  mito
mito.genes <- grep(pattern = "^MT-", x = rownames(x = droncseq.human@data), value = TRUE)
percent.mito <- Matrix::colSums(droncseq.human@raw.data[mito.genes, ])/Matrix::colSums(droncseq.human@raw.data)
droncseq.human<- AddMetaData(object = droncseq.human, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = droncseq.human, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

# normalization 
droncseq.human <- NormalizeData(object = droncseq.human, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

# variable genes
droncseq.human <- FindVariableGenes(object = droncseq.human, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = droncseq.human@var.genes)

# scale data 
droncseq.human <- ScaleData(object = droncseq.human, vars.to.regress = c("nUMI", "percent.mito"))
dim(droncseq.human@scale.data) #26805 14963
droncseq.human <- droncseq.human@scale.data
genes <- intersect(rownames(our.human),rownames(droncseq.human))
our.human <- our.human[genes,]
droncseq.human <- droncseq.human[genes,]
saveRDS(droncseq.human,file="droncseq.scaled.log.umi.Rds")

# avg  druncseq expression data  --------------------------------------
# cell id 
dat <- read.csv(file = "./data/human/Table7_data_info_human.csv",skip = 20,stringsAsFactors = F)[,-6]
colnames(dat) <-sub("X.","",colnames(dat))
all(dat$Cell.ID %in% colnames(droncseq.human)) # True 
dat <- dat %>% filter(! Cluster.ID %in% c(16,17,18))
dat <- dat %>% mutate(Cluster.Name.simple = sub("[1-3]","",Cluster.Name))
dat <- dat %>% column_to_rownames("Cell.ID")
dat <- dat[genes,]

# merge data
tmp <- as.data.frame(t(droncseq.human) )
tmp <- cbind(tmp,dat[rownames(tmp),c("Cluster.Name","Cluster.Name.simple")])
droncseq.human<- tmp

droncseq.human.l <- droncseq.human %>% gather(key="gene",value = "val",1:3287)

# avg
transfunc <- function(x) log2(sum(2^x-1)/length(x)+1)
droncseq.human.avg  <- droncseq.human.l %>% group_by(gene,Cluster.Name) %>% 
  summarise(mean.exp=mean(val,na.rm=T))


# spread 
droncseq.human.avg.s <- droncseq.human.avg%>% spread(key = Cluster.Name,value = mean.exp)
druncseq.mouse <- dat.exp.log.avg.s[,-ncol(dat.exp.log.avg.s)]

druncseq.mouse <- as.data.frame(druncseq.mouse)
rownames(druncseq.mouse) <- druncseq.mouse$variable; druncseq.mouse$variable <- NULL 

druncseq.mouse.2 <- dat.exp.log.avg.2%>% spread(key = Cluster.name,value = mean.exp)
druncseq.mouse.2<- druncseq.mouse.2[,-ncol(druncseq.mouse.2)]

druncseq.mouse.2 <- as.data.frame(druncseq.mouse.2)
rownames(druncseq.mouse.2) <- druncseq.mouse.2$variable; druncseq.mouse.2$variable <- NULL 





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
