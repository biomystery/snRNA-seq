require(tidyverse)
require(ggplot2)
require(data.table)
# compare Habib2016 -------------------------------------------------------
snucseq.dat.exp.log <- read.csv(file="./data/Table_S1_marker_genes_hip_habib2016.csv",skip = 2,stringsAsFactors = F)
head(snucseq.dat.exp.log)
tail(snucseq.dat.exp.log)
dim(snucseq.dat.exp.log)

snucseq.dat.exp.log<- rbind(snucseq.dat.exp.log[,1:8] %>% gather(key = Cell.type,value = AvglogTPM,2:8),
                            snucseq.dat.exp.log[,10:ncol(snucseq.dat.exp.log)] %>% gather(key = Cell.type,value = AvglogTPM,2:7)%>%
                              rename(Gene.names=Gene.names.1),stringsAsFactors=F) %>% distinct(.)



druncseq.gaba.info <- read.csv(file = "./data/mouse/Table5_Mouse_CGABAergic.csv",skip = 2,stringsAsFactors = F)
head(druncseq.gaba.info)
tail(druncseq.gaba.info); dim(druncseq.gaba.info) # 544 genes

druncseq.gaba.cell <- dat %>% filter(grepl("GABA",Cluster))
dim(druncseq.gaba.cell) #816 nuclei == Fig. S6


# Reproduce-subclustering GABA -------------------------------------------------
dat.exp <- read.delim("./data/mouse/Mouse_Processed_GTEx_Data.DGE.log-UMI-Counts.txt",stringsAsFactors = F)
druncseq.gaba.DMG<- dat.exp[,c("GENE",druncseq.gaba.cell$NAME),with=F]
dim(druncseq.gaba.DMG) #17308 (g) * 816 (n)
druncseq.gaba.DMG[1:3,1:3]

druncseq.gaba.DMG.long <- (druncseq.gaba.DMG %>% gather(key = Cell.ID,value="logUMI",-1))

# 1. find variable genes
#druncseq.gaba.var <- transpose(druncseq.gaba.DMG)
druncseq.gaba.var <- druncseq.gaba.DMG.long %>%
  group_by(GENE) %>%
  summarise(avg= mean(logUMI),
            avg_2 = mean(exp(logUMI)-1),
            var_2 = var(exp(logUMI)-1),
            var = var(logUMI))%>%
  mutate(cv2 = var/(avg*avg),
         cv2_2 = var_2/(avg_2*avg_2))
head(druncseq.gaba.var)

write.csv(druncseq.gaba.var, "02.Expression-variation_stat.csv",row.names = F)

# check avg expression
ggplot(druncseq.gaba.var %>% gather(cv_type,val,c(cv2_2,cv2)),aes(x=cv_type,val))+ geom_violin()
ggplot(druncseq.gaba.var %>% gather(cv_type,val,c(avg,avg_2)),aes(x=cv_type,val))+ geom_violin()
ggplot(druncseq.gaba.var,aes(avg_2,cv2_2))+geom_point()

quantile(druncseq.gaba.var$avg_2,probs = c(0,.25,.5,0.75,.95,1))
#0%         25%         50%         75%         95%        100%
# 0.000000000 0.003405705 0.022464497 0.067741485 0.215173325 4.765387307

druncseq.gaba.var.threshold <- 0.005
druncseq.gaba.var.filter <- druncseq.gaba.var %>% filter(complete.cases(.))%>%
  filter(avg>=druncseq.gaba.var.threshold)

inv_avg = 1 / druncseq.gaba.var.filter$avg
fit <- lm(druncseq.gaba.var.filter$cv2 ~ inv_avg) #fitting CV^2 and mean expression to an inverse distribution
druncseq.gaba.var.filter$l <- fit$coefficients[2] / druncseq.gaba.var.filter$avg + fit$coefficients[1]
druncseq.gaba.var.filter$se <- summary(fit)$sigma
druncseq.gaba.var.filter$se_1 <- druncseq.gaba.var.filter$l + druncseq.gaba.var.filter$se

ggplot(druncseq.gaba.var.filter,aes(avg,cv2,colour=(cv2>=se_1))) + geom_point()+
  geom_line(aes(avg,l),linetype=2,linewidth=1,colour="red")+
  geom_point(aes(avg,se_1), colour="red",shape=1)

quantile(druncseq.gaba.var.variable$avg,probs = c(0,.25,.5,.75,.95,1))

#         0%         25%         50%         75%         95%        100%
#0.005004546 0.006507657 0.008329588 0.011710077 0.019001383 0.050264507

druncseq.gaba.var.variable <- druncseq.gaba.var.filter %>% filter(cv2>se_1)
quantile(druncseq.gaba.var.variable$avg_2,probs = c(0,.25,.5,.75,.95,1))
write.table(rownames(m), "04.Over-dispersed_gene-list.txt", quote = F, sep = "\t", row.names = F, col.names = F)

# 2. PCA

require(rsvd)
head(druncseq.gaba.var.variable)
druncseq.gaba.DMG.var <- druncseq.gaba.DMG.long%>% filter(GENE %in% druncseq.gaba.var.variable$GENE)

druncseq.gaba.DMG.var.mat <- druncseq.gaba.DMG.var %>% spread(Cell.ID,logUMI)
dim(druncseq.gaba.DMG.var.mat)
druncseq.gaba.DMG.var.mat[1:3,1:3]
rownames(druncseq.gaba.DMG.var.mat) <- druncseq.gaba.DMG.var.mat$GENE;druncseq.gaba.DMG.var.mat$GENE <- NULL
druncseq.gaba.DMG.var.mat.scaled <- t(scale(t(druncseq.gaba.DMG.var.mat)))

pd <- druncseq.gaba.DMG.var.mat.scaled; pd[pd>3] <- 3; pd[pd < -3 ] <-3
pheatmap(pd,scale = "none",show_rownames = F,show_colnames = F)

druncseq.gaba.var.rpca<- rpca(t(druncseq.gaba.DMG.var.mat))
druncseq.gaba.var.rpca.summary <- summary(druncseq.gaba.var.rpca)
druncseq.gaba.var.rpca.summary[,1:5]

ggscreeplot(druncseq.gaba.var.rpca,type = "cum")
dp <- druncseq.gaba.var.rpca.summary[4,-1]-druncseq.gaba.var.rpca.summary[4,-ncol(druncseq.gaba.var.rpca.summary)]
plot(1:length(dp),dp,type="l")


ggcorplot(druncseq.gaba.var.rpca,alpha = .3,top.n = 10)
ggindplot(druncseq.gaba.var.rpca,ellipse = T)


# 3.

# seurat pacakge ----------------------------------------------------------
require(Seurat)
rownames(druncseq.gaba.DMG) <- druncseq.gaba.DMG$GENE; druncseq.gaba.DMG$GENE <- NULL
setDF(druncseq.gaba.DMG,rownames = druncseq.gaba.DMG$GENE)
class(druncseq.gaba.DMG);rownames(druncseq.gaba.DMG)

druncseq.gaba.DMG.filter <- druncseq.gaba.DMG[apply(druncseq.gaba.DMG,1,max)>0,]
range(apply(druncseq.gaba.DMG.filter,1,max)) # 0.11 to 6.15
range(apply(druncseq.gaba.DMG.filter,2,max)) # 3.73 to 6.15

tissue <- CreateSeuratObject(raw.data = druncseq.gaba.DMG.filter,do.scale = T)
tissue <- FindVariableGenes(object = tissue, mean.function = ExpMean, dispersion.function = LogVMR,
                            x.low.cutoff = 0.005)
tissue <- RunPCA(tissue,do.print = FALSE, pc.genes = tissue@var.genes,pcs.compute = 20)

VizPCA(object = tissue, pcs.use = 1:2)
PCAPlot(object = tissue, dim.1 = 1, dim.2 = 2)


# singificatn pc
tissue <- JackStraw(object = tissue, num.replicate = 100, do.print = FALSE)
(p.jsp<-JackStrawPlot(object = tissue, PCs = 1:20))
p.score <- unlist(strsplit(levels(p$data$PC.Score),split = " "))
p.score.df <- data.frame(pc=p.score[seq(1,40,by=2)],
                         pval= as.numeric(p.score[seq(2,40,by=2)]))
pcs <- which(p.score.df$pval<0.01)
length(pcs)
PCElbowPlot(object = tissue)

tissue=FindClusters(tissue,reduction.type = "pca",dims.use = pcs,
                    resolution = 1.2,print.output = 0,k.param = 20,
                    save.SNN = T,force.recalc = T)



# tSNE
tissue <- RunTSNE(tissue, dims.use = pcs, do.fast = T)

save(tissue, file = "../data/mouse_gaba.Robj")

# biomarker
if(T){

  #
  fcon= file("../data/gaba.txt")
  sink(fcon)
  PrintFindClustersParams(tissue)
  cat("# variable genes:",length(tissue@var.genes))
  cat("range(rowMax(tissue@scale.data)):",range(rowMax(tissue@scale.data)))
  sink()
  close(fcon)


  #
  gaba.glist <- c("Pvalb", "Cck","Cnr1","Sst","Cxcl14",  "Ndnf", "Npy","Vcan",
                  "Vip","Htr3a","Rgs12")

  plist <- list()

  (plist[[1]]<-TSNEPlot(tissue,do.label = TRUE))
  plist[[2]]<-VlnPlot(object = tissue, features.plot =gaba.glist)

  pdf(file = "gaba.pdf",onefile = T,width = 10,height = 7)
  for(i in 1:2)
    print(plist[[i]])
  FeaturePlot(object = tissue, features.plot = gaba.glist, cols.use = c("grey", "blue"),
              reduction.use = "tsne")
  print(p.jsp)
  PCHeatmap(object = tissue, pc.use = 1:12, cells.use = 500, do.balanced = TRUE,
            label.columns = FALSE, use.full = FALSE)

  dev.off()
}
