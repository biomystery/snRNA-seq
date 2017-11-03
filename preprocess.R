# read data ---------------------------------------------------------------

#wc -l GTEx_droncseq_hip_pcf.umi_counts.txt
# 32112 genes 
# systems ("head -n1 GTEx_droncseq_hip_pcf.umi_counts.txt | awk  '{print NF; exit}' ") 
# 14693 nucs 

# wc -l GTEx_droncseq_hip_pcf.clusters.txt 
# 14963 

# https://pubpeer.com/publications/E38493AECAF3189FC0708887F9EC39


require(Seurat)

# create obj --------------------------------------------------------------
file ="./data/GTEx_droncseq_hip_pcf.umi_counts.txt"
file<- "./data/human/Human_Processed_GTEx_Data.DGE.UMI-Counts.txt"
name ="hip_pcf"
require(data.table)
tissue.data[1:10,1:10]

process_umi<-function(file,name,min.cell=10,min.genes=200,...){
  tissue.data <- fread(file);  setDF(tissue.data,rownames = tissue.data$V1);  tissue.data$V1 <- NULL
  tissue <- CreateSeuratObject(raw.data = tissue.data,project=name, normalization.method = "LogNormalize",do.scale=TRUE,...)
  mito.genes <- grep("^mt-", rownames(tissue@data), value = T)
  percent.mito <- colSums(as.array(expm1(tissue@data[mito.genes, ]))) / colSums(as.array(expm1(tissue@data)))
  tissue <- AddMetaData(tissue, percent.mito, "percent.mito")
  VlnPlot(tissue, c("nGene", "nUMI", "percent.mito"), nCol = 3)
  GenePlot(tissue, "nUMI", "nGene")
  return (tissue)
}

tissue <- process_umi(file,name)


# find variable genes  ----------------------------------------------------
tissue <- FindVariableGenes(tissue,do.plot = T)

length(tissue@var.genes) # 1249 genes 


