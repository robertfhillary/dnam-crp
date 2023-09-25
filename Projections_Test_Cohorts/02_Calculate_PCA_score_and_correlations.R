####################################################################################
###### SCRIPT 2 : CALCULATE SCORE #5 - PCA+ELNET SCORE based on Wielscher CpGs  ####
####################################################################################

# Set working directory
#setwd("filepath/")

###############################################
####### 1. Load requisite libraries ###########
###############################################

# Check if missing and install if so 
list.of.packages <- c("data.table", "tidyverse","corrplot")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load packages
library(data.table)
library(tidyverse)
library(corrplot)

#####################################################
######## 2. Create functions for tidying data #######
#####################################################

# M2Beta
m_to_beta <- function (val) 
{
  beta <- 2^val/(2^val + 1)
  return(beta)
}

# Mean Impute 
meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)


#################################################################
#### 3. Read in scores from Script 1 - to merge in at end  ######
#################################################################

# Read in previous scores
scores=read.csv("crp_scores_lesspca.csv")

##################################################
#### 4. Read in PCA objects for prediction  ######
##################################################

# Common CpG sites - will need later 
inter1=readRDS("cpgs_overlap.rds")

# Prediction object
CalcPCCRP=readRDS("prediction_object.rds")

########################################################
#### 5. Read in methylation and phenotype objects ######
########################################################

# Read in methylation dataframe 
datMethTest=readRDS("methylation_example.rds")

## Of note, you will need IDs as rows and CpGs as columns before progressing
    # this is the opposite form of what is expected in script_01 # 
    # however, script should handle this for you # 
# Check if Data needs to be Transposed
if(nrow(datMethTest) > ncol(datMethTest)){
  message("It seems that CpGs are rows - data will be transposed")
  datMethTest<-as.data.frame(t(datMethTest))
}

# Subset to CpGs used for PCA-based prediction 
datMethTest=datMethTest[,which(colnames(datMethTest) %in% inter1)]

# Sanity check - are all CpGs present? 
print(paste0("There are ", length(which(inter1 %in% colnames(datMethTest))), " of ", length(inter1), " sites present"))

# Align CpGs 
inter1=inter1[which(inter1 %in% colnames(datMethTest))]
datMethTest=datMethTest[,match(inter1,colnames(datMethTest))]

# Read in phenotype data 
datPhenTest=read.csv("pheno_example.csv")

# Format phenotype file 
datPhenTest=datPhenTest[,c("Basename","CRP")]
datPhenTest=datPhenTest[which(!is.na(datPhenTest$CRP)),]
row.names(datPhenTest)=datPhenTest$Basename
TestCRP=datPhenTest$CRP

# Align files - if order of IDs differ between files 
datPhenTest=datPhenTest[which(row.names(datPhenTest)%in%row.names(datMethTest)),]
datMethTest=datMethTest[which(row.names(datMethTest)%in%row.names(datPhenTest)),]
ids=row.names(datPhenTest)
datMethTest=datMethTest[match(ids,row.names(datMethTest)),]

#########################################
#### 6. Prepare files for analysis ######
#########################################

# M values to beta values if present 
datMethTest<-if((range(datMethTest,na.rm=T)> 1)[[2]] == "TRUE") { message("Suspect that M Values are present. Converting to Beta Values");m_to_beta(datMethTest) } else { message("Suspect that Beta Values are present");datMethTest}

# Mean impute if present 
datMethTest <- apply(datMethTest,2,meanimpute)

# Ensure samples are in order 
if(all(rownames(datMethTest) == rownames(datPhenTest))){ #may need to change rownames(datPhenoTest) to a column name
  message("Samples are all in order")
} else(message(paste("Only",sum(rownames(datMethTest) == rownames(datPhenTest)),"Samples are in order")))


#########################################
######### 7. Prediction step ############
#########################################

# Set seed 
set.seed(1234)

# Get projected score 
PC.CRP <- sweep(as.matrix(datMethTest),2,CalcPCCRP$center) %*% CalcPCCRP$rotation %*% CalcPCCRP$model + CalcPCCRP$intercept

# Store scores as dataframe 
PC.CRP=as.data.frame(PC.CRP)
names(PC.CRP)[1]="Wielscher_PCA" 

# Tidy up data frame and merge with other four scores (along with CRP measure)
PC.CRP$Basename=row.names(PC.CRP)
scores1=merge(scores,PC.CRP,by="Basename",all.x=T)

# Format file before proceeding 
scores1=scores1[,c("Basename","CRP","Linear_Hillary","Wielscher","Wielscher_PCA","BayesPR","Elnet")]

write.csv(scores1, "crpscores.csv")

###############################################
######### 8. Run correlation tests ############
###############################################

# Create output dataframe 
mat=as.data.frame(matrix(nrow=5,ncol=5))

# Name columns 
names(mat)=c("Score","r","LCI","HCI","p")

# Store score names 
mat[,1]=names(scores1)[3:ncol(scores1)]

# Store their correlation with CRP 
mat[,2]=apply(scores1[,3:ncol(scores1)], 2, function(x) cor.test(x,scores1[,2])$estimate)

# Get 2.5% CI and 97.5% CI 
mat[,3]=apply(scores1[,3:ncol(scores1)], 2, function(x) cor.test(x,scores1[,2])$conf.int[1])
mat[,4]=apply(scores1[,3:ncol(scores1)], 2, function(x) cor.test(x,scores1[,2])$conf.int[2])

# Get p-value 
mat[,5]=apply(scores1[,3:ncol(scores1)], 2, function(x) cor.test(x,scores1[,2])$p.value)

