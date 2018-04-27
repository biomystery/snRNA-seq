rm(list=ls())
require(ggplot2)
require(tidyverse)
require(data.table)
druncseq.gaba <- list()
# CellID annotation -----------------------------------------------------
druncseq.gaba.cellID <- read.delim(file = "./data/mouse/Mouse.GABAergic.tSNE.And.Sub-CLusters.txt",skip = 1,
  stringsAsFactors = F)

druncseq.gaba$cellID <- druncseq.gaba.cellID %>% select("TYPE","group")%>% 
  rename("CellID"="TYPE","Cluster"="group")


# expression data  --------------------------------------------------------
druncseq.gaba.logUMI <- fread("./data/mouse/Mouse_Processed_GTEx_Data.DGE.log-UMI-Counts.txt")

all(druncseq.gaba$cellID$CellID %in% colnames(druncseq.gaba.logUMI))
druncseq.gaba.logUMI[1:3,1:3]
setDF(druncseq.gaba.logUMI)
druncseq.gaba$logUMI <- druncseq.gaba.logUMI[,c("GENE",druncseq.gaba$cellID$CellID)]
dim(druncseq.gaba$logUMI)



# variable genes from Seraut pacakge --------------------------------------
load(file = "./mouse_gaba_3.Robj")
write.csv(tissue@var.genes,file = "variable_gene_list.txt",quote = F,row.names = F,col.names = F)


# find virable genes  -----------------------------------------------------
#curve(dgamma(x, scale=1.5, shape=2),from=0, to=15, main="Gamma distribution")
druncseq.gaba.sum <- data.frame( gene = rownames(druncseq.gaba$logUMI),
                                 m = apply(2^druncseq.gaba$logUMI-1,1,mean),
                                cv = apply(2^druncseq.gaba$logUMI-1,1,sd)) %>% 
  mutate(cv = cv/m)



model.gamma <- with(druncseq.gaba.sum,glm(cv ~ m, family = "Gamma"))
model.gamma.sum <- summary(model.gamma)
y.fit <- predict(model.gamma)

plot(druncseq.gaba.sum$m,(druncseq.gaba.sum$cv))
points(1/druncseq.gaba.sum[complete.cases(druncseq.gaba.sum),]$m,y.fit,col=2,pch=1)

model.gamma.final <- data.frame(druncseq.gaba.sum[complete.cases(druncseq.gaba.sum),],
                                cv.predict =y.fit,stringsAsFactors = F)%>%
  mutate(cv.excess = cv -cv.predict)%>%
  arrange(desc(cv.excess)) %>%
  filter(m >= 0.005 & abs(cv.excess)>=0.2)

head(model.gamma.final,n = 100)

# Save data  --------------------------------------------------------------
rownames(druncseq.gaba$logUMI)<- druncseq.gaba$logUMI$GENE;druncseq.gaba$logUMI$GENE <- NULL  
druncseq.gaba$variableGenes <- tissue@var.genes
all(tissue@var.genes %in% rownames(druncseq.gaba$logUMI))
saveRDS(file = "druncseq.gaba.Rds",object = druncseq.gaba)



# all cells  --------------------------------------------------------------
druncseq.sum <- data.frame( gene = druncseq.gaba.logUMI$GENE,
                            m = apply(2^druncseq.gaba.logUMI[,-1]-1,1,mean),
                            cv = apply(2^druncseq.gaba.logUMI[,-1]-1,1,sd)) %>% 
  mutate(cv = cv/m)
druncseq.gamma <- with(druncseq.sum,glm(cv ~ I(1/m), family = Gamma(link = "log")))

druncseq.gamma.final <- data.frame(druncseq.sum[complete.cases(druncseq.sum),],
                                cv.predict =predict(druncseq.gamma),stringsAsFactors = F)%>%
  mutate(cv.excess = cv -cv.predict)%>%
  arrange(desc(cv.excess)) 

with(druncseq.gamma.final,expr = {
  par(mfrow=c(1,2))
  plot(m,cv)
  lines(m,cv.predict,col=2)
  plot(m,cv.excess)
})

druncseq.gamma.final.2 <- druncseq.gamma.final %>% 
  filter(m >= 0.005 & cv.excess >0.2)

# test --------------------------------------------------------------------

glmGamma <- glm(response ~ I(1/x1, family = Gamma(link = "identity"))
library(MASS)
myshape <- gamma.shape(glmGamma)
gampred <- predict(glmGamma , type = "response", se = T, dispersion = 1/myshape$alpha) 
summary(glmGamma, dispersion = 1/myshape$alpha)

# others ------------------------------------------------------------------
ggplot(druncseq.gaba.sum[s,],aes(1/m,cv))+geom_point() + geom_smooth(method = rlm) +
  geom_abline(intercept = fit$coefficients[1], slope = fit$coefficients[2])+
  geom_abline(intercept = fit.2$coefficients[1], slope = fit.2$coefficients[2],colour='red')+
  geom_abline(intercept = fit.3$coefficients[1], slope = fit.3$coefficients[2],colour='green')
fit <- with(druncseq.gaba.sum[],rlm(1/cv~m))
fit.2 <- with(druncseq.gaba.sum[],lm(1/cv~m))

s <- quantile(druncseq.gaba.sum$m,probs = c(0.01,.99))
s <- (druncseq.gaba.sum$m > s[1] & druncseq.gaba.sum$m <s[2])
fit.3 <- with(druncseq.gaba.sum[s,],lm(1/cv~m))

require(MASS)
require(car)
with(druncseq.gaba.sum,hist(m,breaks = seq(0,5,by = .2)))
# other -------------------------------------------------------------------
# reproduce the heatmap 
ggplot(druncseq.gaba.cellID,aes(numeric,numeric.1)) + geom_point(aes(colour=as.factor(group)))


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


