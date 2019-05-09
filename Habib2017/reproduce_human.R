

# tsne-plot ---------------------------------------------------------------
dat.tsne <- read.delim("./data/human/GTEx_droncseq_hip_pcf.tsne.txt")
dim(dat.tsne)

dat.cluster <- read.delim("./data/human/GTEx_droncseq_hip_pcf.clusters.txt",header = F)
dim(dat.cluster)

# replot 
require(ggplot2)
require(tidyverse)


dat <- right_join(dat.tsne,dat.cluster,by=c("X"="V1"))
dat.region <- read.delim("./data/human/Human_DroNc-seq_Regions.txt",stringsAsFactors = F)[-1,]
dat <- (right_join(dat,dat.region,by=c("X"="NAME")))

dat[with(dat,which(V2!= as.integer(Cluster))),]
head(dat)
dat$V2 <- NULL 
colnames(dat) <- sub("X","Neucli",colnames(dat))

ggplot(dat,aes(tSNE_1,tSNE_2)) + geom_point(aes(colour=Cluster),size=1) 
ggplot(dat%>% filter(!Cluster%in%c(15,17,18,19)),aes(tSNE_1,tSNE_2)) + geom_point(aes(colour=Cluster),size=1) 

p.gaba <- ggplot(dat%>% filter(Cluster%in%c(5,6)),aes(tSNE_1,tSNE_2)) + geom_point(aes(colour=Cluster),size=1) 
dat.sum <- dat %>% group_by(Cluster)%>% summarise(num=n())

saveRDS(list(dat,dat.sum),file = "human_tsne.Rds")



# mouse data  -------------------------------------------------------------
dat.tsne <- read.delim("./data/mouse/Mouse_Coordinates2.txt",stringsAsFactors = F)[-1,]
dim(dat.tsne); head(dat.tsne)
dat.cluster <- read.delim("./data/mouse/Mouse_Meta_Data_with_cluster.txt",stringsAsFactors = F)[-1,]
head(dat.cluster)

dat <- dat.tsne %>% right_join(dat.cluster,by="NAME")
dat <- type_convert(as.tibble(dat)) %>% unite(Cluster.ID,c("ClusterID","Cluster"),sep = ".")
dat.sum <- dat%>% group_by(Cluster.ID) %>% summarise(num=n())

# replot 
require(ggplot2)
require(tidyverse)

p.fig1a <- ggplot(dat,aes(X,Y)) + geom_point(aes(colour=Cluster.ID),size=1) 
p.figs3c_1<- ggplot(dat,aes(X,Y)) + geom_point(aes(colour= Genes<400),size=1) +
  scale_color_manual(name="<400 genes",values = setNames(c("cyan","grey"),c(T,F)))
p.figs3c_2<- ggplot(dat,aes(X,Y)) + geom_point(aes(colour= Genes<300),size=1) +
  scale_color_manual(name="<400 genes",values = setNames(c("cyan","grey"),c(T,F)))


saveRDS(list(dat,dat.sum),file = "mouse_tsne.Rds")

# replot-2 gene expression 
dat.exp <- read.delim("./data/mouse/Mouse_Processed_GTEx_Data.DGE.log-UMI-Counts.txt",stringsAsFactors = F)
