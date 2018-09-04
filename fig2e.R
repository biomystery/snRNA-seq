# Human data compared to druncseq 

## pariwise correlation of the average expression for the genes 
## in each cell-type signature defined by our data and cell types defined by Drunc-seq
source('./libs.R')
# load our expression data  -----------------------------------------------

# expression 
our.human <- fread("./data/intro_clustered_1809/fc_hc_intronclustered_umi_upregmarkers.txt")
setDF(our.human)

dim(our.human)# 4186 genes 
sum(duplicated(our.human$V1)) #1557 duplicates

our.human <- unique.array(our.human)

# load druncseq expression data  --------------------------------------

## process droncseq data
if(F){
  # expression 
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
  }
droncseq.human <- readRDS(file = "./data/droncseq.scaled.log.umi.Rds")

# avg  druncseq expression data  --------------------------------------
# cell id 
droncseq.human.id <- read.csv(file = "./data/human/Table7_data_info_human.csv",skip = 20,stringsAsFactors = F)[,-6]
colnames(droncseq.human.id) <-sub("X.","",colnames(droncseq.human.id))
all(droncseq.human.id$Cell.ID %in% colnames(droncseq.human)) # True 
droncseq.human.id <- droncseq.human.id %>% filter(! Cluster.ID %in% c(16,17,18))%>%
  mutate(Cluster.Name.simple = sub("[1-3]","",Cluster.Name))
head(droncseq.human.id)

# add cell id to gene expression 
droncseq.human.long <- as.data.frame(droncseq.human) %>% 
  rownames_to_column("gene") %>% 
  gather(key = Cell.ID,value = scaled.log.umi,-"gene") 
  
df1 <- data.table(droncseq.human.long,key = "Cell.ID")
df2 <- data.table(droncseq.human.id,key = "Cell.ID")

droncseq.human.long <- merge(df1,df2)  
setDF(droncseq.human.long)
range(droncseq.human.long$scaled.log.umi) #-5,10 
rm(df1);rm(df2)

# avg
droncseq.human.avg <- droncseq.human.long %>% 
  dplyr::select(-one_of("Genes","Transcripts","Cell.ID","Cluster.ID"))%>%
  group_by(gene,Cluster.Name) %>% 
  summarise(avg.exp = mean(scaled.log.umi))

droncseq.human.avg.simple <- droncseq.human.long %>% 
    dplyr::select(-one_of("Genes","Transcripts","Cell.ID","Cluster.ID"))%>%
  group_by(gene,Cluster.Name.simple) %>% 
    summarise(avg.exp = mean(scaled.log.umi))

transfunc <- function(x) log2(sum(2^x-1)/length(x)+1)


# spread 
droncseq.human.avg.s <- droncseq.human.avg%>%
  spread(key = Cluster.Name,value = avg.exp) %>% 
  as.data.frame %>% 
  column_to_rownames("gene")

droncseq.human.avg.ss <- droncseq.human.avg.simple%>%
  spread(key = Cluster.Name.simple,value = avg.exp)%>%
  as.data.frame %>% 
  column_to_rownames("gene")

rm(droncseq.human.avg.simple);rm(droncseq.human.avg)
rm(droncseq.human.id)

# correlation matrix  -----------------------------------------------------
our.human <- our.human[rownames(droncseq.human.avg.s),]
calcCormat <- function(a,b) apply(a,2,function(x) apply(b,2,function(y) cor(x,y,method = "spearman")))

cor.mat <- calcCormat(our.human,droncseq.human.avg.s)
cor.mat.s <- calcCormat(our.human,droncseq.human.avg.ss)

pd <- (cor.mat%>% tbl_df()%>% rownames_to_column("Druncseq")%>% gather(ours,value,- "Druncseq"))
ggplot(pd,aes(Druncseq,ours))+
  geom_tile(aes(fill=value))+
  geom_text(aes(label=round(value,1)))+
  scale_fill_gradient2(low = rgb(14/255,135/255,182/255),mid ="white",high = rgb(203/255,72/255,85/255))+theme_bw()

if(T){
  pdf(file =   "fig2e.pdf",width = 7,height = 7)
  cols=colorRampPalette(c(rgb(14/255,135/255,182/255),
                          "white",
                          rgb(203/255,72/255,85/255)))(20)
  pheatmap(cor.mat,scale = "none",cluster_cols = T,cluster_rows = T,cols,
           cellwidth = 12,cellheight = 12,breaks = seq(-0.5,0.5,by = 0.05),
           main = "Pair-wised Spearman's correlation \nof avg expression between our data (row) and \n drunc-seq (column) from our signature genes (Human)")
  pheatmap(cor.mat.s,scale = "none",cluster_cols = T,cluster_rows = T,cols,
           cellwidth = 12,cellheight = 12,breaks = seq(-0.5,0.5,by = 0.05),
           main = "Pair-wised Spearman's correlation \nof avg expression between our data (row) and \n drunc-seq (column) from our signature genes (reduced,human)")
  dev.off()
}

saveRDS(list(druncseq.human.avg= droncseq.human.avg.s,
             druncseq.human.avg.s=droncseq.human.avg.ss,
             druncseq.human = droncseq.human.long,
             our.human=our.human),file = "fig2e.Rds")
