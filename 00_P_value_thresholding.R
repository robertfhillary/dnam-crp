####### Load requisite libraries ###########
library(glmnet)
library(tidyverse)
library(tidyr)

# Mean Impute 
meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)

######## Create functions for tidying data #######
# M2Beta
m_to_beta <- function (val) 
{
  beta <- 2^val/(2^val + 1)
  return(beta)
}

# Mean Impute 
meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)

# Read in Training Set Methylation 
Train=readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/normalised-betas.rds")
# Extract rownames from test set
data=load("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/Beta_3525_norm_bgcorrect_0.001BetaThreshold_probefilter.RObject")
Test=dat

# We want to loop through different thresholds 
cpgs.loop=list.files("/Cluster_Filespace/Marioni_Group/Rob/CRP/Thresholds/Predictors/",".")
cpgs.loop=rev(cpgs.loop) # reversing order so smallest ones can go first as test runs 

for(probes in cpgs.loop){ 
# Read in CpG subsets 
cpgs=read.csv(paste0("/Cluster_Filespace/Marioni_Group/Rob/CRP/Thresholds/Predictors/",probes))

# Tidy up methylation dataframes 
datMethTrain=Train[which(row.names(Train) %in% row.names(Test)),]
datMethTest=Test[which(row.names(Test) %in% row.names(Train)),]
datMethTrain=as.data.frame(t(datMethTrain))
datMethTest=as.data.frame(t(datMethTest))

# Extract threshold of interest 
nm=gsub("Predictors_","", probes)

# Read in phenotypes 
datPhenoTrain=read.csv(paste0("/Cluster_Filespace/Marioni_Group/Rob/CRP/GS_EWAS/crp_resids.csv"))
PhenoTest=read.csv("/Cluster_Filespace/Marioni_Group/Rob/CRP/Projected_Scores/scores_LBC36W2.csv")

datPhenoTrain=datPhenoTrain[,c("Sample_Sentrix_ID","crp")]
datPhenoTest=PhenoTest[,c("Basename","CRP")]

datPhenoTrain=datPhenoTrain[which(!is.na(datPhenoTrain$crp)),]
datPhenoTest=datPhenoTest[which(!is.na(datPhenoTest$CRP)),]

datPhenoTrain=datPhenoTrain[which(datPhenoTrain$Sample_Sentrix_ID %in% row.names(datMethTrain)),]
datPhenoTest=datPhenoTest[which(datPhenoTest$Basename %in% row.names(datMethTest)),]

## Match order 
ids1=datPhenoTrain$Sample_Sentrix_ID
ids2=datPhenoTest$Basename

datMethTrain=datMethTrain[match(ids1,row.names(datMethTrain)),]
datMethTest=datMethTest[match(ids2,row.names(datMethTest)),]

## Tidy up data 
row.names(datPhenoTrain)=datPhenoTrain$Sample_Sentrix_ID
row.names(datPhenoTest)=datPhenoTest$Basename

datPhenoTrain$Sample_Sentrix_ID=NULL
datPhenoTest$Basename=NULL

names(datPhenoTrain)[1] <- "CRP"
names(datPhenoTest)[1] <- "CRP"

#Subset methylation data to your set of CpGs
CpGs <- cpgs$CpG
datMethTrain1 <- datMethTrain[,which(colnames(datMethTrain) %in% CpGs)]
datMethTest1 <- datMethTest[,which(colnames(datMethTest) %in% CpGs)]
  
## Match order of CpGs
CpGs=CpGs[which(CpGs %in% colnames(datMethTrain1))]
ids=colnames(datMethTrain1)
CpGs=CpGs[match(ids,CpGs)]
datMethTest1=datMethTest1[,match(ids,colnames(datMethTest1))]
  
#Impute missing values if needed. You can use a different imputation method of your choice but we have not found
#     this makes a significant difference.
meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)
datMethTrain1 <- apply(datMethTrain1,2,meanimpute)
datMethTest1 <- apply(datMethTest1,2,meanimpute)
  
datMethTrain1<-if((range(datMethTrain1,na.rm=T)> 1)[[2]] == "TRUE") { message("Suspect that M Values are present. Converting to Beta Values");m_to_beta(datMethTrain1) } else { message("Suspect that Beta Values are present");datMethTrain1}
datMethTest1<-if((range(datMethTest1,na.rm=T)> 1)[[2]] == "TRUE") { message("Suspect that M Values are present. Converting to Beta Values");m_to_beta(datMethTest1) } else { message("Suspect that Beta Values are present");datMethTest1}
  
#check order of data is correct
if(all(colnames(datMethTrain1) == CpGs)){
message("CpGs are all in order")
  } else(message(paste("Only",sum(colnames(datMethTrain1) == CpGs),"CpGs are in order")))
  
  if(all(rownames(datMethTrain1) == rownames(datPhenoTrain))){ #may need to change rownames(datPhenoTrain) to a column name
    message("Samples are all in order")
  } else(message(paste("Only",sum(rownames(datMethTrain1) == rownames(datPhenoTrain)),"Samples are in order")))
  
# Set seed for projections   
set.seed(1234)
#Perform PCA and projections. 
PCA = prcomp(datMethTrain1,scale.=F)
TrainPCData = PCA$x[,1:(dim(PCA$x)[2])]
TestPCData = predict(PCA,as.matrix(datMethTest1))[,1:(dim(PCA$x)[2])]
  
#Select phenotype to be predicted
TrainCRP = datPhenoTrain$CRP
TestCRP = datPhenoTest$CRP
  
#Train PC predictor. Can test different models using different alpha and lambda parameters (see glmnet documentation)
cv = cv.glmnet(TrainPCData, TrainCRP, nfolds=8,alpha=0.5, family="gaussian")
fit = glmnet(TrainPCData, TrainCRP, family="gaussian", alpha=0.5, nlambda=100)

#Most likely your final model will only use a small subset of PCs. Thus you can compress your model:
CalcPCCRP <- vector(mode = "list",length = 0)
temp = as.matrix(coef(cv,s = cv$lambda.min))
CalcPCCRP$model = temp[temp!=0,][-1]
CalcPCCRP$intercept = temp[1,1]
CalcPCCRP$center = PCA$center
CalcPCCRP$rotation = PCA$rotation[,names(CalcPCCRP$model)]
PC.CRP <- sweep(as.matrix(datMethTest1),2,CalcPCCRP$center) %*% CalcPCCRP$rotation %*% CalcPCCRP$model + CalcPCCRP$intercept

# Run elnet on original training data - not PCA transformed 
lasso.cv <- cv.glmnet(datMethTrain1, TrainCRP, family="gaussian", alpha = 0.5, ncores = 8, nfolds = 20) # cross validation to get best lambda
fit2 <- glmnet(datMethTrain1, TrainCRP, family = "gaussian", alpha = 0.5, ncores = 8, lambda = lasso.cv$lambda.min) # model fit 

# Project into test data 
# gather data
ElnetCRP <- vector(mode = "list",length = 0)
temp = as.matrix(coef(lasso.cv,s = cv$lambda.min))
ElnetCRP$model = temp[temp!=0,][-1]
weights=as.data.frame(ElnetCRP$model)
names(weights)[1]="Beta"
ElnetCRP$intercept = temp[1,1]
# project 
prediction=datMethTest1[,which(colnames(datMethTest1)%in% names(ElnetCRP$model))]
ids=names(ElnetCRP$model)
prediction=t(prediction[,match(ids,colnames(prediction))])
tmp=prediction*(weights$Beta)
Elnet.CRP=as.data.frame(colSums(tmp))

# Project using original Wielscher weights 
CpGs1=cpgs[which(cpgs$CpG %in% colnames(datMethTest1)),]
ids=CpGs1$CpG
prediction2=t(datMethTest1[,match(ids,colnames(datMethTest1))])
tmp2=prediction2*(CpGs1$Beta)
Wielscher.CRP=as.data.frame(colSums(tmp2))

# Tidy up CRP objects and combine 
names(PC.CRP)[1]="PC_CRP"
names(Elnet.CRP)[1]="Elnet_CRP"
names(Wielscher.CRP)[1]="Wielscher_CRP"
scores=cbind(PC.CRP,Elnet.CRP,Wielscher.CRP)
scores$Basename=row.names(scores)

# Merge within original phenotype dataframe and write out 
out=merge(PhenoTest[,c("Basename","CRP")],scores,by="Basename")
write.csv(out, paste0("/Cluster_Filespace/Marioni_Group/Rob/CRP/Thresholds_v2/Scores/",probes), row.names=F)

## Remove temporary objects and clean up
rm(prediction)
rm(prediction2)
rm(tmp)
rm(tmp2)
rm(cv)
rm(fit)
rm(lasso.cv)
rm(fit2)
gc()

# Print to denote completion 
print(probes)
}




