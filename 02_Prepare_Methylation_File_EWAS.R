########################################################################
###### 2.0 PREPARATION OF METHYLATION FILE - EWAS  #####################
########################################################################

# Load requisite libraries 
library(data.table)
library(lumi)

# Set working  directory 
setwd("/Cluster_Filespace/Marioni_Group/Rob/CRP/GS_EWAS/Inputs/") 

# Loop through each chromosome 
for(j in 1:22){ 
# Read in methylation file - by chromosome 
tmp=readRDS(paste0("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/Chromosomes/GS20k_chr",j,"_mvals.rds"))
# Convert to betas 
tmp1=m2beta(tmp)
# Function to define and trim outliers 
saveRDS(tmp1,paste0("chr",j,".rds"),compress=F)
# Remove temporary dataframes 
rm(tmp)
rm(tmp1)
gc()
# Print to denote completion 
print(j)
} 