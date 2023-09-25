###################################################
#### 10. BAYESR+ ANALYSIS OF CRP ##################
###################################################

## Set working directory 
setwd("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/Chromosomes/") 

## Load requisite libraries
library(data.table)
library(limma)
library(lumi)

## Create function for mean imputation
meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)


#########################################################
######## STEP 1 - PREPARATION OF PHENOTYPE FILE ########
#########################################################

# Read in phenotype from OSCA (already prepared)
phen=read.csv("/Cluster_Filespace/Marioni_Group/Rob/CRP/GS_EWAS/crp.csv")

# Read in with methylation basenames  
samps=readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")

## Subset to those with complete genetic data 
gen=as.data.frame(fread("/Cluster_Filespace/Marioni_Group/GS/GS_GWAS/GS_GWAS.fam"))
gen=gen[which(gen$V2 %in% samps$Sample_Name),]
samps=samps[which(samps$Sample_Name %in% gen$V2),]

# Merge phenotype data in with remaining basenames 
phenos=merge(phen,samps[,c("Sample_Name","Sample_Sentrix_ID")],by.x="id",by.y="Sample_Name")

## Make additional IDs file where we have FID and IID 
fam=read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/clinical/pedigree.csv")
fam1=merge(fam,phenos,by.x="volid", by.y="id",all.y=T)
gwas.ids <- data.frame(FID =fam1$famid,
                       IID = fam1$volid)
write.table(gwas.ids,"/Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/IDs/crp_ids.list",row.names=F,quote=F,sep="\t")

# Residualise phenotype for EWAS 
phenos$crp_resid=scale(resid(lm(phenos$crp ~ phenos$age + phenos$sex)))
ids=phenos$Sample_Sentrix_ID.x

# Extract premade one for prediction model 
crp2=read.csv("/Cluster_Filespace/Marioni_Group/Rob/CRP/GS_EWAS/crp_resids_rel.csv")
crp2=crp2[match(ids,crp2$Sample_Sentrix_ID),]

# Write out phenotype 
# Unresidualised version 
write.table(x = t(as.matrix(as.numeric(scale(phenos[,2])))),file="/Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/Phenotypes/crp.csvphen",quote = F, sep = ",", row.names = F, col.names = F)
# Residualised version 
write.table(x = t(as.matrix(as.numeric(phenos[,13]))),file="/Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/Phenotypes/crp_resid.csvphen",quote = F, sep = ",", row.names = F, col.names = F)
# For prediction (common sites across cohorts)
write.table(x = t(as.matrix(as.numeric(scale(crp2[,41])))),file="/Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/Phenotypes/crp_prediction.csvphen",quote = F, sep = ",", row.names = F, col.names = F)


##########################################################
######## STEP 2 - PREPARATION OF METHYLATION FILE ########
##########################################################

# Read in methylation file # by chromosome 
for(i in 1:22){ 
  meth=readRDS(paste0("GS20k_chr", i, "_mvals.rds"))
  
  # Subset to those with phenotype in question 
  meth=meth[,which(colnames(meth) %in% phenos$Sample_Sentrix_ID.x)]
  
  # Subset to probes passing QC 
  probes=read.table("/Cluster_Filespace/Marioni_Group/Elena/gs_osca/data/cpgs_tokeep.txt", header=F)
  meth=meth[which(row.names(meth) %in% probes$V1),]
  
  # Subset to 450k only 
  meth=meth[which(row.names(meth)%in%anno$Name),]
  
  # Match order of IDs in phenotype and methylation file 
  ids=phenos$Sample_Sentrix_ID.x
  meth=meth[,match(ids,colnames(meth))]
  
  # Check order of IDs match between phenotype and methylation files 
  table(colnames(meth)==phenos$Sample_Sentrix_ID.x)
  
  # Convert to beta values 
  meth=m2beta(meth)
  
  # Mean impute - cannot have NAs in final file 
  meth <- apply(meth,1,meanimpute)
  
  # Transpose back original format - as apply changes the format
  meth=t(meth)
  
  # Prepare covariate matrix for regressions 
  
  # Match order of IDs with other files 
  samps=samps[match(ids,samps$Sample_Sentrix_ID),]
  
  # Check order of IDs match with other files 
  table(samps$Sample_Sentrix_ID==phenos$Sample_Sentrix_ID.x)
  table(samps$Sample_Sentrix_ID==colnames(meth))
  
  ## Regression step - residualise for age, sex and batch 
  design.resid <- model.matrix(~sex + age + Batch, data=samps)
  fit.resid <- limma::lmFit(meth, design.resid)
  gc()
  meth <- limma::residuals.MArrayLM(fit.resid, meth)
  meth <- meth[!is.infinite(rowSums(meth)), ]
  rm(fit.resid)
  gc()
  
  ## Scale 
  meth=t(apply(meth,1,scale))
  
  ## Write out CpGs 
  cpgs=as.data.frame(row.names(meth))
  names(cpgs)[1]="CpG"
  fwrite(cpgs, paste0("/Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/CpGs/GS20k_chr", i, "_subset_cpgs.txt"),row.names=F)
  
  # Save out residualised file 
  fwrite(meth, paste0("/Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/Methylation_unresid/GS20k_chr", i, "_resid.txt"),row.names=F)
  
  ## Remove methylation object and clean up environment 
  rm(meth)
  gc()
  print(i)
} 

## Write out ID order for later
saveRDS(ids, "/Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/IDs/crp_meth_ids.rds")

###############################################################
######## STEP 3 - PREPARE GENETIC FILES  ######################
###############################################################

####### Recode step - not in R ######

plink19 \
--bfile /Cluster_Filespace/Marioni_Group/GS/GS_GWAS/GS_GWAS \
--keep /Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/IDs/crp_ids.list \
--recodeA \
--out /Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/Genetics_Raw/crp_recode \


### Back in R #####

## Read in resulting file in R for preparation step 
snp=as.data.frame(fread("/Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/Genetics_Raw/crp_recode.raw"))

## Match order with other files 
ids1=phenos$id
snp=snp[match(ids1,snp$IID),]

# Remove redundant columns
snp[,1:6]=NULL

# Extract SNP names 
snps=as.data.frame(colnames(snp))
names(snps)[1]="Marker"
fwrite(x =as.matrix(snps), "/Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/SNPs/crp_snp_list.csv", sep = ",", row.names = F, col.names = T, quote = F) 

# Mean impute any missing values 
snp <- apply(snp,2,meanimpute)

# Transpose back original format - as apply changes the format
# snp=t(snp)

# Scale each value 
snp <- apply(snp,2,scale)

snp=t(snp)

# Save out file 
fwrite(x =as.matrix(snp), "/Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/Genetics/crp_17904.csv", sep = ",", row.names = F, col.names = F, quote = F) 



###############################################################
######## STEP 4 - PREPARE FIXED EFFECT COVARIATE FILES ########
###############################################################

#### HERE YOU MAY WANT TO TRUNCATE TO JUST WBCS - I JUST ADAPTED THIS TO MOST COVARIATES AS I HAD VARIOUS MODELS TO CHECK 

# Read in WBCs
wbc=read.table("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/wbc_quant.qcov", header=T)

# Combine with sample ID information - to help merge in with other covariates 
wbc1=merge(wbc,samps[,c("Sample_Sentrix_ID", "Sample_Name", "age","sex","Batch")],by.x="FID",by.y="Sample_Sentrix_ID")

# Match order of IDs with other files 
wbc1=wbc1[match(ids,wbc1$FID),]

## Remove Sample ID information 
wbc1$FID=NULL
wbc1$IID=NULL 
wbc1$Sample_Name=NULL 
wbc1$sex=ifelse(wbc1$sex%in%"M",0,1)
wbc1$Batch=as.numeric(as.factor(wbc1$Batch))
wbc1[1:ncol(wbc1)]=apply(wbc1[,1:ncol(wbc1)],2,scale)

# Write out covariate file
write.table(x = as.matrix(wbc1),file = "/Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/crp_full_covariates.csv" ,quote = F, sep = ",", row.names = F, col.names = F)



##############################################################
######## STEP 5 - COMBINE INDIVIDUAL CHROMOSOME FILES ########
##############################################################

# Change working directory 
setwd("/Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/Methylation/") 

# Extract files 
files=list.files(".", ".txt")
files=files[order(files)]

# Read in and rbind all methylation files 
data <- rbindlist(lapply(files,fread))

## Write out final file
#fwrite(x =as.matrix(data), "/Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/Methylation/crp_unresid_17904.csv", sep = ",", row.names = F, col.names = F, quote = F) 
fwrite(x =as.matrix(data), "/Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/Methylation/crp_pred_17904.csv", sep = ",", row.names = F, col.names = F, quote = F) 

# Extract CpGs - will need this for processing final results files as the row.names get lost in methylation file in BayesR
cpgs=list.files("../CpGs/", ".txt")
cpgs=cpgs[order(cpgs)]

# Ensure that order is same as methylation file just created 
# methylation
ids.cpg=gsub("GS20k_chr", "", files)
ids.cpg=gsub("_.*", "", ids.cpg)
# cpgs
ids1.cpg=gsub("GS20k_chr", "", cpgs)
ids1.cpg=gsub("_.*", "", ids1.cpg)
# is order same? 
table(ids.cpg==ids1.cpg)

# Read in cpg lists 
# Change working directory 
setwd("/Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/CpGs/") 
cg <- rbindlist(lapply(cpgs,fread))
names(cg)[1]="Marker"
# Save out file
fwrite(x =as.matrix(cg), "/Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/CpGs/crp_cpg_subset_list.csv", sep = ",", row.names = F, col.names = T, quote = F) 



##############################################################
######## STEP 6 - COMBINE GENETICS + METHYLATION FILE ########
##############################################################

## Combine with genetics file to make multi-omics file 
new.data=rbind(data,snp)

# Remove other files to clear up environment
rm(data)
rm(snp)
gc()

## Write out final file
fwrite(x =as.matrix(new.data), "/Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/Combined/crp_17904.csv", sep = ",", row.names = F, col.names = F, quote = F) 

# Remove object and clear up space
rm(new.data)
gc()

# Combine markers to make groups file 
cg$Group=0
snps$Group=1
new.markers=rbind(cg,snps)
new.markers[,1]=as.character(new.markers[,1])
# Save out file 
write.table(x=new.markers, "/Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/Groups/crp_2.txt",  sep=' ', row.names = F, col.names = F) 


#####################################
#### STEP 7 - BAYESR  ###############
#####################################

## Epigenetics - no residualisation - mirrors linear EWAS  
cd /Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/

../../../BayesRRcmd/src/brr --data-file Methylation/crp_unresid_17904.csv --pheno Phenotypes/crp.csvphen --analysis-type preprocess --fixed_effects crp_full_covariates.csv --fixedEffectNumber 8 --thread 12 --thread-spawned 12 --marker-cache --seed 1 

../../../BayesRRcmd/src/brr --data-file Methylation/crp_unresid_17904.csv --pheno Phenotypes/crp.csvphen --fixed_effects crp_full_covariates.csv --fixedEffectNumber 8 --analysis-type ppbayes --chain-length 10000 --burn-in 5000 --thin 5 --S "0.001,0.01,0.1" --mcmc-samples Outputs/crp_unresid.csv --thread 12 --thread-spawned 12 --marker-cache --seed 1



## Epigenetics - standard way  
cd /Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/
  
../../../BayesRRcmd/src/brr --data-file Methylation/crp_17904.csv --pheno Phenotypes/crp_resid.csvphen --analysis-type preprocess --fixed_effects crp_basic_covariates.csv --fixedEffectNumber 5 --thread 12 --thread-spawned 12 --marker-cache --seed 1 

../../../BayesRRcmd/src/brr --data-file Methylation/crp_17904.csv --pheno Phenotypes/crp_resid.csvphen --fixed_effects crp_basic_covariates.csv --fixedEffectNumber 5 --analysis-type ppbayes --chain-length 10000 --burn-in 5000 --thin 5 --S "0.001,0.01,0.1" --mcmc-samples Outputs/crp_resid.csv --thread 12 --thread-spawned 12 --marker-cache --seed 1



## Epigenetics - prediction  
cd /Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/
  
../../../BayesRRcmd/src/brr --data-file Methylation/crp_pred_17904.csv --pheno Phenotypes/crp_prediction.csvphen --analysis-type preprocess --thread 12 --thread-spawned 12 --marker-cache --seed 1 

../../../BayesRRcmd/src/brr --data-file Methylation/crp_pred_17904.csv --pheno Phenotypes/crp_prediction.csvphen --analysis-type ppbayes --chain-length 10000 --burn-in 5000 --thin 5 --S "0.001,0.01,0.1" --mcmc-samples Outputs/crp_prediction.csv --thread 12 --thread-spawned 12 --marker-cache --seed 1


## Genetics 
cd /Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/
  
../../../BayesRRcmd/src/brr --data-file Genetics/crp_17904.csv --pheno Phenotypes/crp_resid.csvphen --analysis-type preprocess --thread 24 --thread-spawned 24 --marker-cache --seed 1 

../../../BayesRRcmd/src/brr --data-file Genetics/crp_17904.csv --pheno Phenotypes/crp_resid.csvphen  --analysis-type ppbayes --chain-length 10000 --burn-in 5000 --thin 5 --S "0.01,0.1" --mcmc-samples Outputs/crp_genetics.csv --thread 12 --thread-spawned 12 --marker-cache --seed 1



## Combined 
cd /Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/
  
../../../BayesRRcmd/src/brr --data-file Combined/crp_17904.csv --pheno Phenotypes/crp_resid.csvphen --analysis-type preprocess --fixed_effects crp_basic_covariates.csv --fixedEffectNumber 5 --thread 24 --thread-spawned 24 --marker-cache --seed 1 

../../../BayesRRcmd/src/brr --data-file Combined/crp_17904.csv --pheno Phenotypes/crp_resid.csvphen --fixed_effects crp_basic_covariates.csv --fixedEffectNumber 5 --analysis-type ppbayes --chain-length 5000 --burn-in 2500 --thin 5 --group Groups/crp_2.txt --S "0.001,0.01,0.1;0.001,0.01,0.1" --mcmc-samples Outputs/crp_combined_lowergeneticpriors.csv --thread 12 --thread-spawned 12 --marker-cache --seed 1

