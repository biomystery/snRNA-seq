# mouse data  -------------------------------------------------------------
dat.tsne <- read.delim("./data/mouse/Mouse_Coordinates2.txt",stringsAsFactors = F)[-1,]
dim(dat.tsne); head(dat.tsne)
dat.cluster <- read.delim("./data/mouse/Mouse_Meta_Data_with_cluster.txt",stringsAsFactors = F)[-1,]
head(dat.cluster)

dat <- dat.tsne %>% right_join(dat.cluster,by="NAME")
dat <- type_convert(as.tibble(dat)) %>% unite(Cluster.ID,c("ClusterID","Cluster"),sep = ".")
dat <- dat %>% separate(Cluster.ID,into=c("ClusterID","Cluster"),sep="[.]")
dat.sum <- dat%>% group_by(Cluster.ID) %>% summarise(num=n())

# replot 
require(ggplot2)
require(tidyverse)

filter.clustersID <- 22:27
p.fig1a <- ggplot(dat %>% 
                    filter(!ClusterID %in%filter.clustersID)%>%
                    unite(Cluster.ID,c("ClusterID","Cluster"),sep = "."),aes(X,Y)) + 
  geom_point(aes(colour=Cluster.ID),size=1) 


p.figs3c_1<- ggplot(dat,aes(X,Y)) + geom_point(aes(colour= Genes<400),size=1) +
  scale_color_manual(name="<400 genes",values = setNames(c("cyan","grey"),c(T,F)))

p.figs3c_2<- ggplot(dat,aes(X,Y)) + geom_point(aes(colour= Genes<300),size=1) +
  scale_color_manual(name="<400 genes",values = setNames(c("cyan","grey"),c(T,F)))


saveRDS(list(dat,dat.sum),file = "mouse_tsne.Rds")

# replot-2 gene expression 
require(data.table)
dat.exp <- fread("./data/mouse/Mouse_Processed_GTEx_Data.DGE.log-UMI-Counts.txt")
dim(dat.exp)
dat.exp[1:2,1:2]
colnames(dat.exp)

dat.exp.cell <-unique(substr(colnames(dat.exp)[-1],1,4))
dat.cell <- unique(substr(dat$NAME,1,4))
all.equal(dat.exp.cell,dat.cell)

dat.exp.t <- dcast(melt(dat.exp,id.vars = "GENE"),variable ~ GENE)
head(dat.exp.t[1:10,1:10])
head(dat)
?left_join
dat.all <- (dat%>% left_join(dat.exp.t,by=c("NAME"="variable")))

data.all.long <- dat.all %>% gather(value="log2_UMI_counts",key="Gene",8:ncol(dat.all))

p.figs3d <- ggplot(data.all.long,aes(ClusterID,log2_UMI_counts)) + geom_violin()
pd.figs3d <- dat.all %>% group_by(ClusterID)%>% 
