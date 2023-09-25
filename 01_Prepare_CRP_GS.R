#############################################################
##################### 1.0 CRP PREPARATION ###################
#############################################################

# Load requisite libraries 
library(data.table)
library(coxme)
library(kinship)

# Set working directory 
setwd("/Cluster_Filespace/Marioni_Group/Rob/CRP/")

# Read in pedigree data 
ped = read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/clinical/pedigree.csv")

# Tidy up dataframe 
ped$father <- as.numeric(ped$father)
ped$mother <- as.numeric(ped$mother)
ped$father[ped$father==0] <- NA
ped$mother[ped$mother==0] <- NA
ped$sex <- as.numeric(as.factor(ped$sex))
ped$sex[ped$sex==2] <- 0
ped$sex <- ped$sex+1

# Prepare matrix 
kin <- with(ped, pedigree(volid, father, mother, sex, famid=famid))
kin_model <- kinship(kin) 


# Read in CRP data 
crp=read.csv("2022-11-01_gs19335_crp.csv")
# Extract outliers 
med=median(crp$crp,na.rm=T)
sde=sd(crp$crp,na.rm=T)
upper=med+(4*sde)
lower=med-(4*sde)
crp=crp[which(crp$crp<upper),]
crp=crp[which(crp$crp>lower),] 
crp=crp[order(crp$crp),]
# Prepare data 
crp$crp=log(crp$crp+0.01)
#crp$crp=log(crp$crp+1)
crp=crp[which(!is.na(crp$crp)),]
crp=crp[which(!is.infinite(crp$crp)),]

# Prepare covariate files  - age, sex, batch, WBCs 
cov=read.csv("../EWAS_Disease_GS/covariates.csv")
# Extract basic covariates 
cov1=cov
# Subset to those with CRP data 
cov2=cov1[which(cov1$Sample_Name%in%crp$id),]
# Subset CRP data to those with covariate information 
crp1=crp[which(crp$id%in%cov2$Sample_Name),]
# Merge in methylation IDs 
crp1=merge(crp1[,c("id","crp")],cov2,by.x="id",by.y="Sample_Name")
names(crp1)[1]=c("Sample_Name")
# Save out file 
write.csv(crp1,"GS_EWAS/crp_zeroes.csv",row.names=F)

# residualise 
crp1=crp1[which(!is.na(crp1$V4)),]
crp1$resids=resid(lm(crp ~ age + sex + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + V20 + V21 + V22,data=crp1,na.action="na.exclude"))
write.csv(crp1, "GS_EWAS/crp_resids.csv",row.names=F)

# relatedness
crp1$resids_rel=resid(lmekin(crp ~ age + sex + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + V20 + V21 + V22 + (1|Sample_Name), varlist = kin_model*2, data=crp1))
write.csv(crp1, "GS_EWAS/crp_resids_rel.csv",row.names=F)
