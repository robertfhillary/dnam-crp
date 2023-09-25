#########################################################################
################### SCRIPT 1 : CALCULATE SCORES #1-4   ##################
#########################################################################

# Set working directory
#setwd("filepath/")

#####################################################
######## 1. Create functions for tidying data #######
#####################################################

# M2beta 
m_to_beta <- function (val) 
{
  beta <- 2^val/(2^val + 1)
  return(beta)
}

# NA to mean (mean imputation)
na_to_mean <-function(methyl){
  methyl[is.na(methyl)]<-mean(methyl,na.rm=T)
  return(methyl)
}


#####################################
### 2. Read in methylation object ###
#####################################

# Methylation dataframe 
data=readRDS("methylation_example.rds") # Of note, this example is a very truncated dataframe and most CpGs will be missing, please input your typical methylation dataframe here (all available CpGs)


#####################################
### 3. Read in phenotype object #####
#####################################

# Read in CRP (outliers removed based on 4SD in original paper, and data were log-transformed)
pheno=read.csv("pheno_example.csv")


#################################################
### 4. Read in predictors and their weights #####
#################################################

# Read in weights for scores 
cpgs <- read.csv("Predictors_by_groups.csv", header = T) 


#################################################
### 5. Quality control prior to projections #####
#################################################

## Of note, you will need CpGs as rows and IDs as columns before progressing ## 
    # this is the opposite form of what is expected in script_02 # 
    # however, script should handle this for you # 
# Check if Data needs to be Transposed 
if(ncol(data) > nrow(data)){
  print(message("It seems that individuals are rows - data will be transposed!"))
  data<-as.data.frame(t(data))
}

# Subset methylation dataframe to sites for prediction 
coef=data[intersect(rownames(data), cpgs$CpG),]

# Convert M values to Beta values 
coef<-if((range(coef,na.rm=T)> 1)[[2]] == "TRUE") { message("Suspect that M Values are present. Converting to Beta Values");m_to_beta(coef) } else { message("Suspect that Beta Values are present");coef}

# Convert NAs to mean values 
coef <- t(apply(coef,1,function(x) na_to_mean(x)))


#################################################
### 6. Loop stage for projecting scores  ########
#################################################

# Loop through each predictor 
loop = unique(cpgs$Predictor)
# Set up output dataframe 
out <- data.frame()
# Initiate loop 
for(i in loop){ 
  # Find CpGs for predictor i 
  tmp=coef[intersect(row.names(coef),cpgs[cpgs$Predictor %in% i,"CpG"]),]
  tmp_coef = cpgs[cpgs$Predictor %in% i, ]
  # If there is more than one CpG present for the predictor (hopefully will be), then align CpGs across files
  if(nrow(tmp_coef) > 1) { 
    # Align CpGs in your methylation dataframe with those in the predictor file 
    tmp_coef = tmp_coef[match(row.names(tmp),tmp_coef$CpG),]
    # Multiply methylation beta values by weight for each CpG - and take sum for each individual 
    out[colnames(coef),i]=colSums(tmp_coef$Beta*tmp)
  } else {
    # If there is just one CpG, then it will already be aligned with dataframe, so can just calculate and take score 
    tmp2 = as.matrix(tmp)*tmp_coef$Beta 
    out[colnames(coef),i] = colSums(tmp2)
  }
} 


#############################################################
### 7. Prepare final file and merge with pheno file  ########
#############################################################

# Prepare final file 
out$Basename=row.names(out)
# Format 
out=out[,c(ncol(out),(1:(ncol(out)-1)))]
# Merge with phenotype file
out=merge(out, pheno,by="Basename",all.y=T)
# Write out to pass to next script 
write.csv(out, "crp_scores_lesspca.csv", row.names=F)
