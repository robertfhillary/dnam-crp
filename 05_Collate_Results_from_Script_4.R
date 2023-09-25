#############################################################
###### 5.0 COLLATE EWAS RESULTS FROM MANUAL MODELS IN R #####
#############################################################

# Load requisite libraries 
library(data.table)

# Set working directory 
setwd("/Cluster_Filespace/Marioni_Group/Rob/CRP/GS_EWAS/")

# Set dictionary of models to read in 
models=c("Model1_lmekin") # can extend to more if you need to run sensitivity analyses 

# Loop through models 
for(j in models){ 
  # Extract all results files (i.e. chromosomes) for that model
  files=list.files(j,".txt")
  # Read in results files 
  res=lapply(files,function(x) as.data.frame(fread(paste0(j,"/",x))))
  # Combine results files 
  res=as.data.frame(do.call("rbind",res))
  # Subset to 450k only 
  # Save out file 
  fwrite(res, paste0("Outputs/CRP_",j,".txt"),row.names=F)
  # Print lambda 
  chisq <- qchisq(1-res$pvalue,1)
  lambda <- median(chisq)/qchisq(0.5,1)
  print(paste0(j, " ", lambda))
  # Print to denote completion
  print(j)
}
