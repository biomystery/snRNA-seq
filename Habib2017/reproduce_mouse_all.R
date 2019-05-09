require(tidyverse)
require(ggplot2)
require(data.table)

# Reproduce-all clustering -------------------------------------------------

require(Seurat)
dat.exp[1:3,1:3]
all(apply(dat.exp[,-1],MARGIN = 1,max)>0)
setDF(dat.exp)
class(dat.exp)
rownames(dat.exp) <- dat.exp$GENE;dat.exp$GENE <- NULL 

range(apply(dat.exp,1,max)) # 0.67 to 6.15 
range(apply(dat.exp,2,max)) # 2.47 to 6.5


tissue <- CreateSeuratObject(raw.data = dat.exp,do.scale = T)
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
pcs <- which(as.numeric(p.score[seq(2,40,by=2)])<0.01)
length(pcs)
PCElbowPlot(object = tissue)
PCHeatmap(object = tissue, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = tissue, pc.use = 13:20, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = tissue, pc.use = pcs, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)


tissue=FindClusters(tissue,reduction.type = "pca",dims.use = pcs,
                    resolution = 2,print.output = 0,k.param = 50,
                    save.SNN = T,force.recalc = T)

if(T){
  # 
  #fcon= file("gaba_3.txt")
  #sink(fcon)
  PrintFindClustersParams(tissue)
  cat("# variable genes:",length(tissue@var.genes),"\n")
  #cat("range(rowMax(tissue@scale.data)):",range(rowMax(tissue@scale.data)),"\n")
  cat ("#clusters: " ,max(as.numeric(levels(tissue@ident))),"\n")
  cat("#nuclei per cluster:", table(tissue@ident),"\n")
  cat("total nuclei",length(tissue@ident))
  #sink()
  #close(fcon)
  
}
# tSNE
tissue <- RunTSNE(tissue, dims.use = pcs, do.fast = T,dim.embed = 3)

save(tissue, file = "./mouse_gaba_3.Robj")

# biomarker 
if(T){
  
  
  
  #
  gaba.glist <- c("Pvalb", "Cck","Cnr1","Sst","Cxcl14",  "Ndnf", "Npy","Vcan",
                  "Vip","Htr3a","Rgs12")
  library(gridExtra)
  plist <- list()
  
  (plist[[1]]<-TSNEPlot(tissue,do.label = TRUE))
  plist[[2]]<-VlnPlot(object = tissue, features.plot =gaba.glist)
  
  pdf(file = "gaba_3.pdf",onefile = T,width = 10,height = 7)
  for(i in 1:2)
    print(plist[[i]])
  FeaturePlot(object = tissue, features.plot = gaba.glist, cols.use = c("grey", "blue"), 
              reduction.use = "tsne")
  print(p.jsp)
  
  dev.off()
}
