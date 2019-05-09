source('../libs.R')
# TSNE  -------------------------------------------------------------
dat.tsne <- read.delim("../data/mouse/Mouse_Coordinates2.txt",stringsAsFactors = F)[-1,]
dim(dat.tsne); head(dat.tsne)
dat.cluster <- read.delim("../data/mouse/Mouse_Meta_Data_with_cluster.txt",stringsAsFactors = F)[-1,]
head(dat.cluster)

dat <- dat.tsne %>% right_join(dat.cluster,by="NAME")
dat <- type_convert(as.tibble(dat)) %>% unite(Cluster.ID,c("ClusterID","Cluster"),sep = ".")
rm(dat.tsne);rm(dat.cluster)

# replot
filter.clustersID <- 22:27
p.fig1a <- ggplot(dat %>%
                    filter(!ClusterID %in%filter.clustersID)%>%
                    unite(Cluster.ID,c("ClusterID","Cluster"),sep = "."),aes(X,Y)) +
  geom_point(aes(colour=Cluster.ID),size=1)

p.figs3c_1<- ggplot(dat,aes(X,Y)) + geom_point(aes(colour= Genes<400),size=1) +
  scale_color_manual(name="<400 genes",values = setNames(c("cyan","grey"),c(T,F)))

p.figs3c_2<- ggplot(dat,aes(X,Y)) + geom_point(aes(colour= Genes<300),size=1) +
  scale_color_manual(name="<400 genes",values = setNames(c("cyan","grey"),c(T,F)))


saveRDS(list(dat,dat.sum),file = "../dat/mouse_tsne.Rds")


# replot-2 gene expression  ---------------------------------------------
dat.exp.log <- fread("../data/mouse/Mouse_Processed_GTEx_Data.DGE.log-UMI-Counts.txt")
dim(dat.exp.log) #] 17308 13314
dat.exp.log[1:4,1:10]
colnames(dat.exp.log)[1:10]

dat.exp.log.cell <-unique(substr(colnames(dat.exp.log)[-1],1,4))
dat.cell <- unique(substr(dat$NAME,1,4))
all.equal(dat.exp.log.cell,dat.cell)

#rownames(dat.exp.log) <- dat.exp.log$GENE; dat.exp.log$GENE <- NULL
dat.exp.log.t <- dcast(melt(dat.exp.log,id.vars = "GENE"),variable ~ GENE)

head(dat.exp.log.t[1:10,1:10])
head(dat)


dat.all <- (dat%>% left_join(dat.exp.log.t,by=c("NAME"="variable")))
data.all.long <- dat.all %>% gather(value="log2_UMI_counts",key="Gene",8:ncol(dat.all))

saveRDS(data.all.long,file = "../data/mouse_all.Rds")


# # transcript per cluster ------------------------------------------------
dat <- read.csv(file = "../data/mouse/Table3_Data_Info_mouse.csv",skip = 27,stringsAsFactors = F)[,-8]
dat.2 <- read.csv(file = "../data/mouse/Table3_Data_Info_mouse.csv",skip =3,nrows = 22,stringsAsFactors = F)[,1:4]


colnames(dat) <-sub("X.","",colnames(dat))
dat <- dat%>% unite(col = "Clusters",c("Cluster","Cluster.name"),remove = F)

p.figs3d <- ggplot(dat,aes(Cluster,Transcripts)) + geom_violin(aes(fill=Clusters))
dat.filter <- (dat %>% filter(Cell.Type %in% 1:12))
f.idx <- dat$Cell.ID[dat$Cell.Type %in% 1:12]

p.figs3d %+% dat.filter

all.equal(colnames(dat.exp.log)[-1],dat$Cell.ID) # F
(all(colnames(dat.exp.log)[-1] %in% dat$Cell.ID)) # TRUE


# average gene expression  ---------------------------------------------------
dat.exp.log.filter <- dat.exp.log[,f.idx,with=F];
dat.exp.log.t.filter <- dat.exp.log.t[variable%in% f.idx]
dat.exp.log.t.filter <- merge(dat.exp.log.t.filter,dat.filter[,c("Cell.ID","Cell.Type")],by.y = "Cell.ID", by.x = "variable")
dat.exp.log.t.filter[1:6,17305:17310]
dat.exp.log.t.filter[1:6,1:6]
colnames(dat.exp.log.t.filter)[1]<- "Cell.ID"

dat.exp.log.filter.long <- melt(dat.exp.log.t.filter,id.vars = c("Cell.ID","Cell.Type"))
head(dat.exp.log.filter.long)
dat.exp.log.avg <- dat.exp.log.filter.long %>% group_by(variable,Cell.Type) %>%
  summarise(mean.exp=mean(value,na.rm=T))

dat.exp.log.avg.spread <- dat.exp.log.avg%>% spread(key = Cell.Type,value = mean.exp)
pd <- as.data.frame(dat.exp.log.avg.spread[,-1])
rownames(pd) <- as.character(dat.exp.log.avg.spread$variable)
sum(rowMeans(pd)==0)
pd <- pd[rowMeans(pd)!=0,]

cord <- c(1,3,4,5,2,6,7,10,8,9,11,12)
pd <- pd[,cord]
gene.id.dic <- apply(pd,1,which.max)
rord <- order(gene.id.dic)

cell.type.dic <- unique(dat.2$Name.1)[-13];
cell.type.dic<- colnames(pd)<- cell.type.dic[cord]

png(filename = "./fig1b.png",width = 700,height = 900)
pheatmap(pd[rord,],scale = "row",cluster_cols = F,show_rownames = F,cluster_rows = F,color = colorRampPalette(c("deepskyblue4", "white", "red"))(50))
dev.off()



# avg expression 2 (use table 4) ------------------------------------------
druncseq.cell.marker <- read.csv(file = "../data/mouse/Table4_Mouse_Clusters.csv",skip = 29,stringsAsFactors = F)
tail(druncseq.cell.marker);dim(druncseq.cell.marker); #3473
dat.exp.log.filter.long <- dat.exp.log.filter.long%>% mutate(variable = as.character(variable))
all(druncseq.cell.marker$Gene.ID %in% dat.exp.log.filter.long$variable)
sum(druncseq.cell.marker$Gene.ID%in% unique(dat.exp.log.filter.long$variable)) # 2042
grep("Snhg",unique(dat.exp.log.filter.long$variable),value = T)
druncseq.cell.marker$Gene.ID[!(druncseq.cell.marker$Gene.ID%in% unique(dat.exp.log.filter.long$variable))] # 2042
length(unique(dat.exp.log.filter.long$Cell.ID))
str(dat.exp.log.filter.long)
druncseq.cell.marker.avgLog

dat.exp.log.filter.long %>% group_by(variable,Cell.Type) %>%
  summarise(mean.exp=mean(value,na.rm=T))


