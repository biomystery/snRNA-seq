rm(list=ls())
require(caret)
require(ggplot2)
require(tidyverse)
require(data.table)


# prepare the testing data ------------------------------------------------
snucseq.gaba.genemarker <- read.csv(file = "./data/habib2016/Table_S2_marker_genes_GABAergic_habib2016.csv",
                                    stringsAsFactors = F,skip = 2)
snucseq.gaba.genemarker.2 <- read.delim(file = "./data/habib2016/GABAergic_edit_subcluster_marker_gene.txt",
                                    stringsAsFactors = F,skip = 0)

head(snucseq.gaba.genemarker.2)
snucseq.gaba.genemarker<-snucseq.gaba.genemarker.2 %>% gather(key="gaba_region",value = "logTPM",-1)

glists <- c("Pvalb", "Tac1", "Sst", "Npy", "Reln", "Vcan", "Nos1", "Cck", 
             "Cnr1", "Htr3a", "Calb1", "Vip", "Calb2", "Cxcl14")
ngene <- length(glists)

p<- ggplot(snucseq.gaba.genemarker%>% filter(GENE.NAMES %in% glists),aes(gaba_region,logTPM)) + 
  geom_bar(stat = "identity",position = "dodge",aes(fill=GENE.NAMES))
ggplotly(p) # not exactly 


# all cluster data  -----------------------------------------------------
snucseq.gaba.cellID  <- read.delim(file = "./data/habib2016/Coordinates_GABAergic.txt",
                                 stringsAsFactors = F,skip = 1)
head(snucseq.gaba.cellID)
ggplot(snucseq.gaba.cellID,aes(numeric,numeric.1))+geom_point()
snucseq.cellID <- read.delim(file="./data/habib2016/Coordinates_Major_cell_types.txt",
                             stringsAsFactors = F,skip = 1)
dim(snucseq.cellID);rm(snucseq.cellID)


# cell- iD --------------------------------------------------------------
snucseq.cellID <- read.delim(file = "./data/habib2016/CLUSTER_AND_SUBCLUSTER_INDEX.txt",
                             stringsAsFactors = T)[-1,]
table(snucseq.cellID$CLUSTER)
table(snucseq.cellID$SUB.CLUSTER)


# expression matrix  ------------------------------------------------------

snucseq.exp <- fread(file = "./data/habib2016/DATA_MATRIX_LOG_TPM.txt")
snucseq.exp[1:3,1:3]
all(snucseq.gaba$cellID$CellID %in% colnames(snucseq.exp))
setDF(snucseq.exp)
# Prepare data for RF model -----------------------------------------------

snucseq.gaba <- list()
snucseq.gaba$cellID <- snucseq.cellID%>% 
  filter(CLUSTER%in% "GABAergic") %>% 
  select("NAME","SUB.CLUSTER") %>% 
  rename("CellID"="NAME","Cluster"="SUB.CLUSTER")%>%
  mutate(CellID = as.character(CellID))

snucseq.gaba$expr = snucseq.exp[,c("GENE",snucseq.gaba$cellID$CellID)]
dim(snucseq.gaba$expr)
save(file = "snucseq.gaba.Rdata",snucseq.gaba)
