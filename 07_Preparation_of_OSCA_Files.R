################################################################
####### 7.0 OSCA EWAS - PREPARATION OF FILES ###################
################################################################

# Phenotype - first take unresidualised CRP as we will use fixed-effect covariates in OSCA 
# Read in phenotype
crp=read.csv("/Cluster_Filespace/Marioni_Group/Rob/CRP/GS_EWAS/crp_resids.csv")
# Read in pedigree data to extract IDs 
ped=read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/clinical/pedigree.csv")
# Identify common IDs to both files and merge 
ped1=ped[which(ped$volid %in% crp$Sample_Name),]
names(ped1)[1:2]=c("FID","IID")
crp1=merge(crp,ped1[,c("FID","IID")],by.x="Sample_Name",by.y="IID")

###########################
## Write out phenotypes ###
###########################

#### HERE, We will write out CRP data using GS IDs as this will help us link to files used in the EWAS stage 
pheno <- data.frame(FID = crp1$FID, IID = crp1$Sample_Name, phen=crp1$crp)
write.table(pheno, "/Cluster_Filespace/Marioni_Group/Rob/CRP/GS_EWAS/crp.phen", row.names=F, sep=' ')
# Here, we will write out using different IDs to support the REML analyses - the files used for REML currently rely on meth IDs 
pheno.meth <- data.frame(FID = crp1$Sample_Sentrix_ID, IID = crp1$Sample_Sentrix_ID, phen=crp1$crp)
write.table(pheno.meth, "/Cluster_Filespace/Marioni_Group/Rob/CRP/GS_EWAS/crp_reml.phen", row.names=F, sep=' ')


#############################
### Write out covariates ####
#############################

###############################
#### BASIC COVARIATES #########
###############################

# Quant covariates - age, WBCs + PCs 
# This version excludes PCs but we would ideally like to include them 
#qcov1 <- data.frame(FID = crp1$FID, IID = crp1$Sample_Name, age=crp1$age, c1=crp1$NK, c2=crp1$Bcell, c3=crp1$Mono,c4=crp1$CD4T,c5=crp1$CD8T)
#write.table(qcov1, file="/Cluster_Filespace/Marioni_Group/Rob/CRP/GS_EWAS/crp_quant.qcov", row.names=F, sep=' ', quote = T)
# This version includes PCs - lets proceed with these 
qcov1 <- data.frame(FID = crp1$FID, IID = crp1$Sample_Name, age=crp1$age, c1=crp1$NK, c2=crp1$Bcell, c3=crp1$Mono,c4=crp1$CD4T,c5=crp1$CD8T,v1=crp1$V3, v2=crp1$V4, v3=crp1$V5, v4=crp1$V6, v5=crp1$V7, v6=crp1$V8, v7=crp1$V9, v8=crp1$V10, v9=crp1$V11, v10=crp1$V12, v11=crp1$V13, v12=crp1$V14, v13=crp1$V15, v14=crp1$V16, v15=crp1$V17, v16=crp1$V18, v17=crp1$V19, v18=crp1$V20, v19=crp1$V21, v20=crp1$V22)
write.table(qcov1, file="/Cluster_Filespace/Marioni_Group/Rob/CRP/GS_EWAS/crp_quant_PCs.qcov", row.names=F, sep=' ', quote = T)
# This version only includes age - useful for adjusting BOD file if needed (WBCs not wanted there)
qcov1 <- data.frame(FID = crp1$FID, IID = crp1$Sample_Name, age=crp1$age)
write.table(qcov1, file="/Cluster_Filespace/Marioni_Group/Rob/CRP/GS_EWAS/crp_quant_age.qcov", row.names=F, sep=' ', quote = T)

# Fact covariates - sex + batch
crp1$sex=ifelse(crp1$sex%in%"M",0,1)
crp1$Batch=as.numeric(as.factor(crp1$Batch))
cov1 <- data.frame(FID = crp1$FID, IID = crp1$Sample_Name, sex=crp1$sex, batch=crp1$Batch)
write.table(cov1, file="/Cluster_Filespace/Marioni_Group/Rob/CRP/GS_EWAS/crp_fact.cov", row.names=F, sep=' ', quote = T)


#################################
### FULLY-ADJUSTED COVARIATES ###
#################################

# quantitative file 
qcov1 <- data.frame(FID = crp1$FID, IID = crp1$Sample_Name, smk=crp1$smokingScore, alc=crp1$units, bmi=crp1$bmi, ea=crp1$years, rank=crp1$rank, age=crp1$age, c1=crp1$NK, c2=crp1$Gran, c3=crp1$Mono,c4=crp1$CD4T,c5=crp1$CD8T,v1=crp1$V3, v2=crp1$V4, v3=crp1$V5, v4=crp1$V6, v5=crp1$V7, v6=crp1$V8, v7=crp1$V9, v8=crp1$V10, v9=crp1$V11, v10=crp1$V12, v11=crp1$V13, v12=crp1$V14, v13=crp1$V15, v14=crp1$V16, v15=crp1$V17, v16=crp1$V18, v17=crp1$V19, v18=crp1$V20, v19=crp1$V21, v20=crp1$V22)
write.table(qcov1, file="/Cluster_Filespace/Marioni_Group/Rob/CRP/GS_EWAS/crp_quant_full.qcov", row.names=F, sep=' ', quote = T)

# factor file  
cov1 <- data.frame(FID = crp1$FID, IID = crp1$Sample_Name, sex=crp1$sex, batch=crp1$Batch, usual=crp1$usual)
write.table(cov1, file="/Cluster_Filespace/Marioni_Group/Rob/CRP/GS_EWAS/crp_fact_full.cov", row.names=F, sep=' ', quote = T)


##############################################################################
##### IN TERMINAL - ADJUST EXISTING GS BOD FILE AND MAKE ORM FOR REML ########
##############################################################################

### ADJUST ####
osca_Linux --befile bvals-norm20k-18413-831733 --covar crp_fact.cov --qcovar crp_quant_age.qcov --adj-probe --make-bod --out bvals-norm20k-18413-831733-adj
#### MAKE ORM #### 
osca_Linux --befile bvals-norm20k-18413-831733-adj --make-orm --out bvals-norm20k-adj


