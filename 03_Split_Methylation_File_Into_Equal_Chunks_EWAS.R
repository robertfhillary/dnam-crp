#########################################################################################################
######## 3.0 TAKE METHYLATION FILES FROM 2.0 AND SPLIT INTO EQUAL CHUNKS  ###############################
#########################################################################################################

# Set working directory 
setwd("/Cluster_Filespace/Marioni_Group/Rob/Diabetes_MetaEWAS/Inputs/")
# Prepare list to store results 
list=list()
# Loop through each chromosome 
for(j in 1:22){ 
  # Read in methylation file - by chromosome 
  tmp=readRDS(paste0("chr",j,".rds"))
  # Store in list 
  list[[j]]=tmp 
  # Print file 
  print(j)
} 

# Split into 22 equal bins 
list1=as.data.frame(do.call("rbind", list))

start = 1  ## start pos 
stop = 17563 ## stop pos 
for(i in 0:20){ 
  ## sample information file 
  tmp = list1[(start+(stop*i)):(stop*(i+1)),]
  ## write out samps file 
  saveRDS(tmp,file=paste0("chr",(i+1),"_alt.rds"),compress=F)
  # Print to denote completion
  print(i)
} 

# Do chromosome 22 separately 
a=which(row.names(list1) %in% row.names(tmp))
tmp1=list1[(a[length(a)]+1):nrow(list1),]
saveRDS(tmp1,file=paste0("chr22","_alt.rds"),compress=F)

  