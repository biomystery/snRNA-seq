rm(list=ls())
require(caret)
set.seed(100)
# read data  --------------------------------------------------------------
druncseq.gaba <- readRDS(file = "./druncseq.gaba.Rds")

ls(druncseq.gaba)
attach(druncseq.gaba)

# partition data  ---------------------------------------------------------
y <- as.factor(cellID$Cluster)
round(table(y)/length(y),2)
train_idx <- createDataPartition(y,p = .6)$Resample1
round(table(y[train_idx])/length(train_idx),2)

# check name
all.equal(cellID$CellID,colnames(logUMI))

y_train <- y[train_idx]
y_test <- y[-train_idx]

x_train<- t(logUMI[variableGenes,train_idx]) # samples - row; features - column
x_test<- t(logUMI[variableGenes,-train_idx])



# construct the model  ----------------------------------------------------
trainctrl <- trainControl(verboseIter = TRUE)
sink("rf_training.log.txt",append = T)
rf_model <- train(x_train,y_train,method = "rf",trControl = trainctrl)
rf_model
sink()
saveRDS(rf_model,file = "../data/mouse_gaba_rf_model.Rds")

# ROC curve
plot(rf_model)

# predict the test  -------------------------------------------------------
rf_model.predict <- predict(rf_model,newdata = x_test)
rf_model.predict.prob <- predict(rf_model,newdata = x_test,type = "prob")

rf_model.predict.final <- apply(rf_model.predict.prob,1,which.max)
rf_model.predict.final[apply(rf_model.predict.prob,1,max)<0.25]<- 9
confusionMatrix(data=as.factor(rf_model.predict.final),factor(y_test,levels=1:9))

cat("====")
confusionMatrix(data = rf_model.predict,y_test)


# predict the sNuc-seq ---------------------------------------------------
load(file = "snucseq.gaba.Rdata")
ls(snucseq.gaba)
dim(snucseq.gaba$cellID) # 133,2
dim(snucseq.gaba$expr) # 25392 134
sum(variableGenes %in% snucseq.gaba$expr$GENE)
length(variableGenes)
1062-991
variableGenes[!variableGenes %in% snucseq.gaba$expr$GENE]
"Sp3os" %in% snucseq.gaba$expr$GENE
"Sp3os" %in% rownames(logUMI)
rownames(snucseq.gaba$expr) <- snucseq.gaba$expr$GENE; snucseq.gaba$expr$GENE <- NULL
x_snucseq <- snucseq.gaba$expr[variableGenes,]
sum(is.na(x_snucseq))/133
x_snucseq[is.na(x_snucseq)] <- 0
x_snucseq <- t(x_snucseq)
colnames(x_snucseq) <- colnames(x_train)
y_predict <- predict(object = rf_model,newdata = x_snucseq)

y_snucseq <- droplevels(snucseq.gaba$cellID$Cluster)
y_snucseq.2 <- as.factor(as.numeric(y_snucseq))

res.c <- confusionMatrix(data = y_predict,y_snucseq.2)
pd <- as.data.frame(res.c$table) %>% group_by(Reference) %>%
  mutate(Percentage= Freq/sum(Freq)*100)%>%
  mutate(Percentage.2 = cut(Percentage,breaks = c(-0.01,25,50,75,100)))%>%
  mutate(Reference.n = levels(y_snucseq)[Reference])

ggplot(pd,aes(Prediction,Reference.n)) + geom_point(aes(size=Percentage.2,colour=Percentage.2))

# try tuneRF from Randomforest package ------------------------------------
# Algorithm Tune (tuneRF)
sink("rf_training.log.txt",append = T)
cat("tuneRF(x_train, y_train, stepFactor=1.2, improve=1e-5, ntree=1000,trace = T,plot = T) \n")
set.seed(100)
bestmtry <- tuneRF(x_train, y_train, stepFactor=1.2, improve=1e-5,ntree=1000,trace = T,plot = T)
print(bestmtry)
sink()





# mannual search best parameters  -----------------------------------------
# ref: https://machinelearningmastery.com/tune-machine-learning-algorithms-in-r/

manual_rf <- function(){
  customRF <- list()
  customRF$parameters <- data.frame(parameter = c("mtry", "ntree"), class = rep("numeric", 2), label = c("mtry", "ntree"))
  customRF$grid <- function(x, y, len = NULL, search = "grid") {}
  customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
  }
  customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
    predict(modelFit, newdata)
  customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
    predict(modelFit, newdata, type = "prob")
  customRF$sort <- function(x) x[order(x[,1]),]
  customRF$levels <- function(x) x$classes
  return(customRF)
}

# train model
control <- trainControl(method="repeatedcv", number=10, repeats=3)
tunegrid <- expand.grid(.mtry=c(1:15), .ntree=c(1000, 1500, 2000, 2500))
set.seed(seed)
custom <- train(Class~., data=dataset, method=customRF, metric=metric, tuneGrid=tunegrid, trControl=control)
summary(custom)
plot(custom)


# Z-score transform x & y  ------------------------------------------------
sink("rf_training.log.txt",append = T)
cat("===\n NOW use z-value of expression (instead of logUMI) \n")

x_z <- t(scale(t(logUMI[variableGenes,])))
x_train_z<- t(x_z[,train_idx]) # samples - row; features - column
x_test_z<- t(x_z[,-train_idx])

sink("rf_training.log.txt",append = T)
cat("tuneRF(x_train_z, y_train_z, stepFactor=1.5, improve=1e-5, ntree=1000,trace = T,plot = T) \n")
set.seed(100)
bestmtry <- tuneRF(x_train_z, y_train, stepFactor=1.5, improve=1e-5,ntree=1000,trace = T,plot = T)
print(bestmtry)
sink()

# Z-score transform x per cell  ------------------------------------------------
sink("rf_training.log.txt",append = T)
cat("===\n NOW use z-value of expression (normalized in each cell) \n")
x_z <- scale(logUMI[variableGenes,])
x_train_z<- t(x_z[,train_idx]) # samples - row; features - column
x_test_z<- t(x_z[,-train_idx])

set.seed(100)
?tuneRF
bestmtry <- tuneRF(x_train_z, y_train, stepFactor=1.1, improve=1e-5,ntreeTry=1000,trace = T,plot = T,mtryStart = 72)
print(bestmtry)
sink()



detach(druncseq.gaba)



