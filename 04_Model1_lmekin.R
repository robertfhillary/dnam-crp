#######################################################################
##### 4.0 MODEL TO RUN LMEKIN MODELS ACROSS 22 CHROMOSOMES   ##########
#######################################################################

##### NOTE, THIS IS THE RSCRIPT THAT IS USED BY A BASH (.SH) SCRIPT THAT RUNS THE CHROMOSOMES IN PARALLEL #####
###### COPY THIS SCRIPT INTO THE CLUSTER AND USE 'bash model1_lmekin.sh' TO RUN ALL MODELS IN PARALLEL BEHIND SCREENS #####
###### REMEMBER TO TAKE NOTE OF THE APPROPRIATE DIRECTORIES AND MODIFY IF NECESSARY ################

## THE CURRENT DIRECTORY I HAVE PLACED THESE FILES IS /Cluster_Filespace/Marioni_Group/Rob/CRP/GS_EWAS/

args = commandArgs(trailingOnly=TRUE)
variable=args[1]
variable=as.numeric(variable)
# Round one

# Round one

# Load requisite libraries
library(coxme)
library(kinship2)
library(data.table)


###################################################################
##### STEP 1. PREPARE KINSHIP MATRIX FOR RELATEDNESS ANALYSES #####
###################################################################

ped1 <- read.csv("/Cluster_Filespace/Marioni_Group/Hannah/Biomarkers_brain_health_LBC_GS/GenScot_input_data/pedigree_data/2023-03-20_pedigree.csv")
ped2 = data.frame(famid=c(4091,4384), volid=c(103027, 144865), father=c(0,0), mother=c(0,0), sex= c("F", "F"))
ped = rbind(ped1, ped2)
kin <- with(ped, pedigree(volid, father, mother, sex, famid=famid))
kin_model <- kinship(kin) 

############################################################
##### STEP 2. CREATE FUNCTION TO EXTRACT MODEL RESULTS #####
############################################################

# Create function
extract_coxme_table <- function (mod){
  beta <- fixef(mod)
  nvar <- length(beta)
  nfrail <- nrow(mod$var) - nvar
  se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
  z<- round(beta/se, 2)
  p<- signif(1 - pchisq((beta/se)^2, 1), 2)
  table=data.frame(cbind(beta,se,z,p))
  return(table)
}


# Load requisite libraries
library(coxme)
library(kinship2)
library(data.table)

# Set working directory
setwd("/Cluster_Filespace/Marioni_Group/Rob/CRP/GS_EWAS/Inputs/")

# # Read in methylation file
# tmp=readRDS(paste0("chr",variable,".rds"))
# # Prepare for analysis
# tmp=t(tmp)
# tmp=as.data.frame(tmp)
# #tmp$Sample_Sentrix_ID=row.names(tmp)
# #tmp=tmp[,c(ncol(tmp),1:(ncol(tmp)-1))]
# names(tmp)=gsub("value.", "", names(tmp))
# cpgs.sig=read.csv("/Cluster_Filespace/Marioni_Group/Rob/CRP/Tables/basic_ewas_signif.csv")
# tmp=tmp[,which(names(tmp) %in% cpgs.sig$Probe)]
# 
# # Read in phenotype and covariate file
# phen=read.csv("../crp_resids.csv")
# tmp=tmp[which(row.names(tmp)%in%phen$Sample_Sentrix_ID),]
# ids=phen$Sample_Sentrix_ID
# tmp=tmp[match(ids,row.names(tmp)),]
# 
# # Prepare matrix to store results
# mat.res=as.data.frame(matrix(nrow=length(2:ncol(tmp)),ncol=4))
# names(mat.res)=c("CpG","coef","se","pvalue")
# 
# # Loop through CpGs
# for(i in 1:ncol(tmp)){
#   # Prepare model
#   mod = lmekin(crp ~ tmp[,i] + age + factor(sex) + NK + Gran + Mono + CD4T + CD8T + Batch + (1|phen$Sample_Name), varlist = kin_model*2, na.action="na.exclude", data=phen)
#   # Extract model results
#   mod1 <- as.data.frame(extract_coxme_table(mod))
#   # Store CpG
#   mat.res[i,1] <- as.character(names(tmp)[i])
#   # Store Beta and SE
#   mat.res[i,c(2,3)] <- mod1[2,1:2]
#   # Store P
#   mat.res[i,4] <- mod1[2,4]
#   # Print to denote completion
#   print(i)
# }
# # Save out final file
# fwrite(mat.res,paste0("../Model1_lmekin/chr",variable,".txt"),row.names=F)
# 

# Round two

# Load requisite libraries
library(coxme)
library(kinship2)
library(data.table)

# Set working directory
setwd("/Cluster_Filespace/Marioni_Group/Rob/CRP/GS_EWAS/Inputs/")
var2=(variable+1)
# Read in methylation file
tmp=readRDS(paste0("chr",as.character(var2),".rds"))
# Prepare for analysis
tmp=t(tmp)
tmp=as.data.frame(tmp)
#tmp$Sample_Sentrix_ID=row.names(tmp)
#tmp=tmp[,c(ncol(tmp),1:(ncol(tmp)-1))]
names(tmp)=gsub("value.", "", names(tmp))
cpgs.sig=read.csv("/Cluster_Filespace/Marioni_Group/Rob/CRP/Tables/basic_ewas_signif.csv")
tmp=tmp[,which(names(tmp) %in% cpgs.sig$Probe)]

# Read in phenotype and covariate file
phen=read.csv("../crp.csv")
tmp=tmp[which(row.names(tmp)%in%phen$Sample_Sentrix_ID),]
ids=phen$Sample_Sentrix_ID
tmp=tmp[match(ids,row.names(tmp)),]

# Prepare matrix to store results
mat.res=as.data.frame(matrix(nrow=length(2:ncol(tmp)),ncol=4))
names(mat.res)=c("CpG","coef","se","pvalue")

# Loop through CpGs
for(i in 1:ncol(tmp)){
  # Prepare model
  mod = lmekin(crp ~ tmp[,i] + age + factor(sex) + NK + Gran + Mono + CD4T + CD8T + Batch + (1|phen$id), varlist = kin_model*2, na.action="na.exclude", data=phen)
  # Extract model results
  mod1 <- as.data.frame(extract_coxme_table(mod))
  # Store CpG
  mat.res[i,1] <- as.character(names(tmp)[i])
  # Store Beta and SE
  mat.res[i,c(2,3)] <- mod1[2,1:2]
  # Store P
  mat.res[i,4] <- mod1[2,4]
  # Print to denote completion
  print(i)
}
# Save out final file
fwrite(mat.res,paste0("../Model1_lmekin/chr",as.character(var2),".txt"),row.names=F)
