####### Load requisite libraries ###########
library(glmnet)
library(tidyverse)
library(data.table)

######## Create functions for tidying data #######
# M2Beta
m_to_beta <- function (val) 
{
  beta <- 2^val/(2^val + 1)
  return(beta)
}

# Mean Impute 
meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)

# Read in CpG subsets 
cpgs=read.csv("/Cluster_Filespace/Marioni_Group/Rob/CRP/Predictors_by_groups_CRP_calc.csv")

# Read in Training Set Methylation 
datMethTrain=readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/normalised-betas.rds")

# Read in CpGs within target methylation dataset - this is so PCA is performed on same CpG set 
cpgnames=readRDS("/Cluster_Filespace/Marioni_Group/Rob/CRP/Elnet/CpGs_combined_commontoall.rds")

# Find intersection and save out for prediction analyses - next script 
inter=intersect(cpgnames,row.names(datMethTrain))
# Find intersection with wielscher cpgs for this experiment 
cpgs1=cpgs[cpgs$Predictor%in%"wielscher",]
inter1=intersect(inter,cpgs1$CpG)

# Tidy up methylation dataframes 
datMethTrain=datMethTrain[which(row.names(datMethTrain)%in%inter1),]
datMethTrain=as.data.frame(t(datMethTrain))
datMethTrain=datMethTrain[,match(inter1,names(datMethTrain))]

# Read in phenotypes 
datPhenoTrain=read.csv(paste0("/Cluster_Filespace/Marioni_Group/Rob/CRP/GS_EWAS/crp_resids_rel.csv"))
# Format phenotypes 
datPhenoTrain=datPhenoTrain[,c("Sample_Sentrix_ID","crp")]
names(datPhenoTrain)[1]="Basename"
datPhenoTrain=datPhenoTrain[which(!is.na(datPhenoTrain$crp)),]
datPhenoTrain=datPhenoTrain[which(datPhenoTrain$Basename %in% row.names(datMethTrain)),]

## Tidy up data 
row.names(datPhenoTrain)=datPhenoTrain$Basename
datPhenoTrain$Basename=NULL
names(datPhenoTrain)[1] <- "CRP"
datMethTrain1=datMethTrain[which(row.names(datMethTrain)%in%row.names(datPhenoTrain)),]
ids=row.names(datPhenoTrain)
datMethTrain1=datMethTrain1[match(ids,row.names(datMethTrain1)),]

# Clean methylation dataframe 
meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)
datMethTrain1 <- apply(datMethTrain1,2,meanimpute)
datMethTrain1<-if((range(datMethTrain1,na.rm=T)> 1)[[2]] == "TRUE") { message("Suspect that M Values are present. Converting to Beta Values");m_to_beta(datMethTrain1) } else { message("Suspect that Beta Values are present");datMethTrain1}

# Ensure samples are in order 
if(all(rownames(datMethTrain1) == rownames(datPhenoTrain))){ #may need to change rownames(datPhenoTrain) to a column name
  message("Samples are all in order")
} else(message(paste("Only",sum(rownames(datMethTrain1) == rownames(datPhenoTrain)),"Samples are in order")))
  
# Perform PCA and projections
set.seed(1234)
PCA = prcomp(datMethTrain1,scale.=F)
TrainPCData = PCA$x[,1:(dim(PCA$x)[2])]

# Select phenotype to be predicted
TrainCRP = datPhenoTrain$CRP
# Train PC predictor. Can test different models using different alpha and lambda parameters (see glmnet documentation)
cv = cv.glmnet(TrainPCData, TrainCRP, nfolds=20,alpha=0.5, family="gaussian")
fit = glmnet(TrainPCData, TrainCRP, family="gaussian", alpha=0.5, lambda=cv$lambda.min)
# Examine full model
cor(TrainCRP,predict(fit,TrainPCData,s = cv$lambda.min),use="complete.obs")

# Most likely your final model will only use a small subset of PCs. Thus you can compress your model:
CalcPCCRP <- vector(mode = "list",length = 0)
temp = as.matrix(coef(cv,s = cv$lambda.min))
# Save version for plotting later that keeps all PCs 
CalcPCCRP.plot$model=temp[-1,]
CalcPCCRP$model = temp[temp!=0,][-1]
CalcPCCRP$intercept = temp[1,1]
CalcPCCRP$center = PCA$center
CalcPCCRP$rotation = PCA$rotation[,names(CalcPCCRP$model)]

  
# Save out objects necessary for calculation in test set 
saveRDS(inter1, "/Cluster_Filespace/Marioni_Group/Rob/CRP/GS_LBC_Prediction/cpgs_GS_LBC_overlap.rds")
saveRDS(PCA, "/Cluster_Filespace/Marioni_Group/Rob/CRP/GS_LBC_Prediction/PCA_object.rds")
saveRDS(CalcPCCRP, "/Cluster_Filespace/Marioni_Group/Rob/CRP/GS_LBC_Prediction/prediction_object.rds")


datMethTrain2=as.data.frame(datMethTrain1)
fit2=plsr(TrainCRP~.,data=datMethTrain2,validation="CV",ncomp=15,scale=F)
plsCV <- RMSEP(fit2,estimate='CV')
#plot(plsCV,main='')
saveRDS(fit2,"/Cluster_Filespace/Marioni_Group/Rob/CRP/GS_LBC_Prediction/pls_object.rds")

## What do PCs correlate with in training set? 
TrainPCs=TrainPCData[,colnames(TrainPCData) %in% names(CalcPCCRP.plot$model)]
# saveRDS(TrainPCs, "/Cluster_Filespace/Marioni_Group/Rob/CRP/GS_Training_PCs.rds")
# Read in covariate data
cov=read.csv("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/covariates.csv")
# Merge in CRP data
datPhenoTrain$Sample_Sentrix_ID=row.names(datPhenoTrain)
cov=merge(cov,datPhenoTrain,by="Sample_Sentrix_ID")
# Aligns IDs
ids=row.names(TrainPCs)
cov=cov[match(ids,cov$Sample_Sentrix_ID),]
# Extract covariates of interest 
cov1=cov[,c("Sample_Sentrix_ID","age","sex","rank","bmi","units","smokingScore","years","Bcell","Gran","Mono","NK","CD4T","CD8T","Batch", "CRP")]

# Run individual linear regressions with each PC and predictor to obtain strength of association 
# Loop through each predictor (continuous first)
cont=c("CRP","age","rank","bmi","units","smokingScore","years","Bcell","Gran","Mono","NK","CD4T","CD8T")
# Set up list to store results for each predictor
output.lm=list()
# Loop stage 
for(j in cont){
# Set up list to store results for each PC (within each predictor)
output.lm.tmp=list()
# Loop through each PC 
for(i in 1:ncol(TrainPCs)){
# Run lm and store output 
tmp=signif(summary(lm(scale(TrainPCs[,i]) ~ scale(cov1[,j])))$coefficients[2,c(1,2,4)],2)
names(tmp)=c("Beta","SE","P")
output.lm.tmp[[i]]=tmp
print(paste("Completed",i,"of",length(1:ncol(TrainPCs)),"PCs"))
}
# Set list for this predictor as a dataframe and store with other predictors
output.lm[[j]]=as.data.frame(do.call("rbind",output.lm.tmp))
# Show progress
print(paste("Completed",which(cont %in% j),"of",length(cont),"Predictors"))
} 
# Format results 
names(output.lm)=c("CRP","Age","SIMD","BMI","Alcohol Consumption","Smoking Score", "Education","B cells","Granzymes","Monocytes","Natural Killer Cells","CD4T+ Cells","CD8T+ Cells")
# Combine continuous results 
results1 = as.data.frame(do.call("cbind", output.lm))
results1$PC=colnames(TrainPCs)
results1=results1[,c(ncol(results1),(1:(ncol(results1)-1)))]

# Loop through each predictor (categorical second)
cat=c("sex")
# Set up list to store results for each predictor
output.lm=list()
# Loop stage 
for(j in cat){
  # Set up list to store results for each PC (within each predictor)
  output.lm.tmp=list()
  # Loop through each PC 
  for(i in 1:ncol(TrainPCs)){
    # Run lm and store output 
    tmp=signif(summary(lm(scale(TrainPCs[,i]) ~ cov1[,j]))$coefficients[2,c(1,2,4)],2)
    names(tmp)=c("Beta","SE","P")
    output.lm.tmp[[i]]=tmp
    print(paste("Completed",i,"of",length(1:ncol(TrainPCs)),"PCs"))
  }
  # Set list for this predictor as a dataframe and store with other predictors
  output.lm[[j]]=as.data.frame(do.call("rbind",output.lm.tmp))
  # Show progress
  print(paste("Completed",which(cat %in% j),"of",length(cat),"Predictors"))
} 
# Format results 
names(output.lm)=c("Sex")
# Combine continuous results 
results2 = as.data.frame(do.call("cbind", output.lm))
# Combine with previous results 
results=cbind(results1,results2)

# Loop through each predictor (anova third)
anov=c("Batch")
# Set up list to store results for each predictor
output.lm=list()
# Loop stage 
for(j in anov){
  # Set up list to store results for each PC (within each predictor)
  output.lm.tmp=list()
  # Loop through each PC 
  for(i in 1:ncol(TrainPCs)){
    # Run lm and store output 
    tmp = summary(aov(TrainPCs[,i] ~ factor(cov1[,j])))[[1]][1,c(4,5)]
    names(tmp)=c("F","P")
    output.lm.tmp[[i]]=tmp
    print(paste("Completed",i,"of",length(1:ncol(TrainPCs)),"PCs"))
  }
  # Set list for this predictor as a dataframe and store with other predictors
  output.lm[[j]]=as.data.frame(do.call("rbind",output.lm.tmp))
  # Show progress
  print(paste("Completed",which(anov %in% j),"of",length(anov),"Predictors"))
} 
# Format results 
names(output.lm)=c("Batch")
# Combine continuous results 
results3 = as.data.frame(do.call("cbind", output.lm))
# Combine with previous results 
results=cbind(results,results3)

#################
### PLOT 1 ######
#################

# Get betas for continuous variables for plotting 
betas=results[,c(1,grep("Beta",names(results)))]
betas$Sex.Beta=NULL
# Subset to first 20 PCs
betas1=betas[1:20,]
# Remove .Beta portion from colnames
names(betas1)=gsub(".Beta","",names(betas1))
# Replace rownames with PC names
row.names(betas1)=betas1$PC
betas1$PC=NULL
# Plot first 20 PCs 
pdf("/Cluster_Filespace/Marioni_Group/Rob/CRP/Plots/Corrplot_PCs_first20PCs.pdf",width=15,height=15)
corrplot(t(betas1),tl.col = 'black',tl.cex=1.7)
dev.off()


#################
### PLOT 2 ######
#################

# Get betas for continuous variables for plotting 
betas=results[,c(1,grep("Beta",names(results)))]
betas$Sex.Beta=NULL
# Find Top 20 PCs - largest coefficients 
coefs=temp[-1,] #remove intercept 
# Get absolute coefficients and rank 
coefs=abs(coefs)
coefs1=coefs[rev(order(coefs))]
# Get names of top 20 PCs 
top=names(coefs1)[1:20]
betas1=betas[which(betas$PC %in% top),]
# Match order according to top PC 
betas1=betas1[match(top,betas1$PC),]
# Remove .Beta portion from colnames
names(betas1)=gsub(".Beta","",names(betas1))
# Replace rownames with PC names
row.names(betas1)=betas1$PC
betas1$PC=NULL
# Plot first 20 PCs 
pdf("/Cluster_Filespace/Marioni_Group/Rob/CRP/Plots/Corrplot_PCs_top20PCs.pdf",width=15,height=15)
corrplot(t(betas1),tl.col = 'black',tl.cex=1.7)
dev.off()