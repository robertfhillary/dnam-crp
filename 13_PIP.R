## Calculate Posterior Inclusion Probability of CpGs - over 95% inclusion means significant i.e.  PIP > 0.95 

## EPIC ARRAY ###

setwd("/Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/PIP/") 
loop = list.files("../Comp/", pattern = ".csv") 
loop=loop[-grep("prediction",loop)]
library(data.table) 
for(i in loop){
  # Extract name 
  nm = gsub(".csv", "", i)
  # Read in CpG list in order 
  names=read.csv("../CpGs/crp_cpg_list.csv",header=T)
  # Read in components file 
  comp <- fread(paste("../Comp/", i, sep=""))  
  comp<-as.data.frame(comp) 
  # Calculate PIPs 
  pip <- lapply(comp,function(x){length(x[which(x>0)])/length(x)})
  pip <- as.data.frame(reshape2::melt(unlist(pip)))
  pip <- setDT(pip, keep.rownames = TRUE) 
  # Name CpGs 
  names(pip) <- c("Marker", "PIP") 
  pip$Marker <- names$Marker
  # Give trait name in column 
  pip$Trait <- nm 
  # Save out file 
  write.csv(pip, file = paste0(nm, "_pip.csv", sep = ""), row.names = F) 
} 


## Find All CpGs that had PIP > 0.95 and PIP > 0.2 

L <- list.files(".", ".csv")

pip_files <- lapply(L, read.csv, header = T) 
names <- as.character(L) 
names <- gsub(".csv*", "", names) 
names(pip_files) <- names 
pip_files = Map(cbind, pip_files) 

pip_top <- do.call(rbind, lapply(pip_files, function(x)x[x$PIP > 0.95,]))
pip_top1 <- do.call(rbind, lapply(pip_files, function(x)x[x$PIP > 0.2,]))


## Next, query all CpGs that are 2.5 kb < CpGs that have PIP > 0.2 
# Set up list to store dataframe 
list1=list()
# Read in annotation dataframe 
anno=readRDS("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds")
# Establish CpGs to be queried 
cpgs1=unique(pip_top1$Marker)
# Loop of interest 
for(i in 1:length(cpgs1)){ 
# Extract the CpG of interest 
tmp=anno[anno$Name %in% cpgs1[i],]
# Subset annotation to chr and basepairs +/- within 2.5kb 
anno.tmp=anno[anno$chr %in% tmp$chr,]
anno.tmp=anno.tmp[which(anno.tmp$pos <= (tmp$pos+1250)),]
anno.tmp=anno.tmp[which(anno.tmp$pos >= (tmp$pos-1250)),]
anno.tmp$group=as.character(cpgs1[i])
# Store output 
list1[[i]]=anno.tmp
# Print to denote completion 
print(i)
} 
# Combine 
list2=as.data.frame(do.call("rbind",list1))
# Establish if they are the lead CpG being queried 
list2$lead=ifelse(list2$Name == list2$group, "lead", "other")

# Next, ask if those CpGs have correlation coefficient >0.5 with the lead CpG in the group 
meth=readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/mvals.rds")
# Subset to samples passing QC 
samps=readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")
meth=meth[,which(colnames(meth) %in% samps$Sample_Sentrix_ID)]
# Subset to probes of interest 
meth1=meth[which(row.names(meth) %in% list2$Name),]
# Set up list to store new outputs 
keeps=list()
# Loop through each group to get correlation matrix 
for(i in unique(list2$group)){ 
# Subset to that group 
tmp=list2[which(list2$group %in% i),]  
# Subset methylation data to probes 
meth.tmp=meth1[which(row.names(meth1) %in% tmp$Name),]
# Correlation matrix
cor.tmp=cor(as.matrix(t(meth.tmp)))
# Subset correlation matrix to lead probe 
cor.tmp1=cor.tmp[which(row.names(cor.tmp)%in%tmp[tmp$lead%in%"lead","Name"]),]
# Remove the probe being queried as it will have a correlation of 1 with itself 
cor.tmp2=cor.tmp1[which(!names(cor.tmp1)%in%tmp[tmp$lead%in%"lead","Name"])]
# Identify those (if any) with absolute correlation coefficient > 0.5 with probe being queried 
highcor=cor.tmp2[abs(cor.tmp2)>=0.5]
lowcor=cor.tmp2[abs(cor.tmp2)<0.5]
# Separate probes into three groups - the lead CpG, those to be kept based on correlation and those to be removed
tmp$keep=NA
tmp[which(tmp$Name %in% names(highcor)),"keep"]="keep"
tmp[which(tmp$Name %in% names(lowcor)),"keep"]="remove"
tmp[tmp$lead%in%"lead","keep"]="lead"
# Store the outputs in a new file so that we can use it again for the PIP stage 
keeps[[i]]<-tmp
# Print to denote completion
print(which(unique(list2$group) %in% i))
} 

# Combine outputs and remove those that failed the correlation stage 
keeps1=as.data.frame(do.call("rbind",keeps))
keeps1=keeps1[which(keeps1$keep %in% c("lead","keep")),]

names(pip_files)=gsub("_pip", "", names(pip_files))
## Now, we need to get the group PIP for all of the assigned groups 
# Loop through groups 
for(i in unique(keeps1$group)){ 
# Find traits which this was significant for 
traits=pip_top1[pip_top1$Marker %in% i,"Trait"]
traits=gsub("_pip", "", traits)
# Define which files we will loop through for this CpG 
ind=which(names(pip_files) %in% traits)
for(j in ind){ 
# Get file 
file=pip_files[[j]]
# Subset to that group 
tmp=keeps1[which(keeps1$group %in% i),] 
# Calculate group PIP 
tmp$group_PIP=sum(file[which(file$Marker %in% tmp$Name),"PIP"])
tmp$trait=names(pip_files)[j]
# Tidy up file 
tmp1=tmp[,c("Name", "chr", "pos", "UCSC_RefGene_Name", "group", "lead", "keep", "group_PIP", "trait")]
write.csv(tmp1, paste0("/Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/Group_PIP/", i, "_", names(pip_files)[[j]],".csv"),row.names=F)
}
# Print to denote completion 
print(which(unique(keeps1$group) %in% i))
} 


## Combine group PIP results 
setwd("/Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/Group_PIP/")
# Combine files 
pips=list.files(".",".")
pips1=as.data.frame(do.call("rbind", lapply(pips, read.csv, header = T)))
# Only subset to those where group PIP > 0.8
pips2=pips1[pips1$group_PIP >= 0.8,]
pips3=pips2[pips2$lead%in%"lead",]
# Tabulate results 
table(pips3$trait)
