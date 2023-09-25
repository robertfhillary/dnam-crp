# Load packages and functions for use in models 
library(bigmemory)
library(optparse)
library(biglasso)
library(data.table)
library(tidyverse)
library(foreign)
library(lumi)  # might have to install this or other packages 


methx=fread("/Cluster_Filespace/Marioni_Group/Rob/CRP/Elnet/gs20k-betas-scaled-cleaned.csv")
methx=as.data.frame(methx)
methx=t(methx)
methx=as.data.frame(methx)
cpgs=readRDS("/Cluster_Filespace/Marioni_Group/Rob/CRP/Elnet/CpGs_combined.rds")
names(methx)=cpgs

# Read in CpGs to keep 
cpgs1=readRDS("/Cluster_Filespace/Marioni_Group/Rob/CRP/Elnet/CpGs_combined_keep.rds")
methx=methx[,which(names(methx) %in% cpgs1)]

# subset to wielscher probes at p<3.6e-8
# wielscher=as.data.frame(fread("/Cluster_Filespace/Marioni_Group/Rob/CRP/SummStats/wieschler.txt"))
# wielscher=wielscher[wielscher$adj.P<0.05,]

crp=read.csv("/Cluster_Filespace/Marioni_Group/Rob/CRP/GS_EWAS/crp_resids_rel.csv")

# Align IDs 
#methx=methx[,which(names(methx)%in%wielscher$MarkerName)]
#methx=methx[which(row.names(methx) %in% crp$Sample_Sentrix_ID),]
ids=crp$Sample_Sentrix_ID
methx=methx[match(ids,row.names(methx)),]


################################################################################

# Split big lasso matricies using cbindBM - load relevant functions for this

################################################################################

# x needs to be a list of big matrices, and they have to have the same rownames
cbindBM_list <- function(x, binding="right", 
                         z=NULL, type=NULL, separated=NULL,
                         backingfile=NULL, backingpath=NULL,
                         descriptorfile=NULL, binarydescriptor=FALSE,
                         shared=TRUE, erase = TRUE)
{
  
  if (is.null(type)) type <- typeof(x[[1]])
  if (is.big.matrix(x[[1]])) {
    if (is.null(separated)) separated <- is.separated(x[[1]])
  } else {
    separated <- FALSE
  }
  
  cols_list <- list()
  total_cols <- 0
  for (i in 1:length(x)) {
    cols <- cleanupcols(NULL, ncol(x[[i]]), colnames(x[[i]]))
    cols_list <- append(cols_list, list(cols))
    total_cols <- total_cols + ncol(x[[i]])
  }    
  
  if (is.null(z)) {
    z <- big.matrix(nrow=nrow(x[[1]]), ncol=total_cols, type=type, init=NULL,
                    dimnames=dimnames(x[[1]]), separated=separated,
                    backingfile=backingfile, backingpath=backingpath,
                    descriptorfile=descriptorfile,
                    binarydescriptor=binarydescriptor, shared=shared)
  }
  
  counter <- 0
  for (i in 1:length(cols_list)) {
    print(i)
    if (i == 1) {
      z[, 1:length(cols_list[[i]])] <- x[[i]][,cols_list[[i]]]
    } else {
      z[, (counter + 1):(counter + length(cols_list[[i]]))] <- x[[i]][,cols_list[[i]]]
    }
    counter <- counter + length(cols_list[[i]])
    print(counter)
    
    if (erase == TRUE) {
      cat("\nErasing chunk and liberating memory...\n\n")
      x[[i]] <- "Replacement"
      gc()
    }
  }
  return(z)
}


cleanupcols <- function(cols=NULL, nc=NULL, colnames=NULL) {
  if (is.null(cols)) cols <- 1:nc
  else {
    if (!is.numeric(cols) & !is.character(cols) & !is.logical(cols))
      stop("column indices must be numeric, logical, or character vectors.")
    if (is.character(cols))
      if (is.null(colnames)) stop("column names do not exist.")
    else cols <- mmap(cols, colnames)
    if (is.logical(cols)) {
      if (length(cols) != nc)
        stop(paste("column vector length must match the number of",
                   "columns of the matrix."))
      cols <- which(cols)
    }
    tempj <- .Call("CCleanIndices", as.double(cols), as.double(nc), PACKAGE="bigmemory")
    if (is.null(tempj[[1]])) stop("Illegal column index usage in extraction.\n")
    if (tempj[[1]]) cols <- tempj[[2]]
  }
  return(cols)
}

meth=methx
div <- 30 # Number of chunks to divide OG methylation dataframe

por <- ceiling(length(colnames(meth))/div)
chunk_list <- list()

for (i in 1:div) {
  cat(paste0("\nWorking on chunk: ", i, " of ", div))
  if (i == 1) {
    chunk <- as.big.matrix(meth[,1:(por-1)])
  } else if (i == div) {
    chunk <- as.big.matrix(meth[,(por*(i-1)):length(colnames(meth))])
  } else {
    chunk <- as.big.matrix(meth[,(por*(i-1)):((por*i)-1)])
  }
  cat("\nMade chunk. Appending to chunk list...\n")
  chunk_list <- append(chunk_list, list(chunk))
  gc()
}

# Saving names prior to chunk fusing
names <- colnames(meth)
rm(meth)

cat("\nRAM clean up...\n\n")
gc()

cat("\nFusing chunks!\n\n")
methx <- cbindBM_list(x = chunk_list)
rm(chunk, chunk_list)

# Set CpG names
options(bigmemory.allow.dimnames=TRUE)
colnames(methx)<- names



set.seed(1783) # set seed to ensure fold variation minimised 

identical(rownames(methx), crp$Sample_Sentrix_ID) # TRUE

# Assign ytrain
y <- as.numeric(crp$resids_rel) # Create a numeric list for this variable to feed in as y to the model

# Run training 
lasso.cv <- cv.biglasso(methx, y, family="gaussian", alpha = 0.5, ncores = 8, nfolds = 20) # cross validation to get best lambda
# Calculate lambda.1se
cve <- lasso.cv$cve
lambda_values <- lasso.cv$lambda
# Find the index of the lambda value that corresponds to the minimum cross-validated error
min_cve_index <- which.min(cve)
# Calculate the standard error of the cross-validated errors
nfolds=20 #note that i had 20, you may have different 
se_cve <- sd(cve)/sqrt(20)
# Find lambda.1se
lambda_1se <- lambda_values[max(which(cve <= cve[min_cve_index] + se_cve))]
# Select coefficients 
fit <- biglasso(methx, y, family = "gaussian", alpha = 0.5, ncores = 8, lambda = lambda_1se) # model fit 
coefs <- coef(fit) # Get betas
coefs1=as.data.frame(as.matrix(coefs))
coefs1$CpG=row.names(coefs1)
names(coefs1)[1]="Coefficient"
coefs2=coefs1[which(coefs1$Coefficient!=0),]
coefs2=coefs2[-1,]
coefs2$Rank=abs(coefs2$Coefficient)
coefs2$Rank=rank(rev(coefs2$Rank))
coefs2=coefs2[order(coefs2$Rank),]
# Save coefficients
saveRDS(coefs2, "//Cluster_Filespace/Marioni_Group/Rob/CRP/Elnet/CRP_predictor_lambda1se.rds") 



### Troubleshooting - useful points to remember when finished 

# cd /dev/shm
# rm -rf /dev/shm/*
# pkill -u danni

# sessionInfo()


