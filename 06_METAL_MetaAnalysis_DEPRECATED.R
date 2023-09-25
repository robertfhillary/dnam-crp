######################################################################
########  6.0 META-ANALYSIS WITH PREVIOUS EWAS USING METAL ###########
######################################################################


# Load requisite libraries 
library(data.table)

# Set working directory 
setwd("/Cluster_Filespace/Marioni_Group/Rob/CRP/GS_EWAS/")

# Read in Wielscher results and meta-analyze 
tmp1=as.data.frame(fread("../SummStats/wieschler.txt"))
tmp1=tmp1[,1:4]
names(tmp1)=c("probe","coef.pat","se.pat","p.pat") 

# Read in GS meta-analyze 
tmp2=as.data.frame(fread("Outputs/CRP_Model1.txt"))
cpgs=as.data.frame(fread("/Cluster_Filespace/Marioni_Group/Elena/gs_osca/data/cpgs_tokeep.txt",header=F))
tmp2=tmp2[which(tmp2$CpG %in% cpgs$V1),]
names(tmp2)=c("probe","coef.pat","se.pat","p.pat") 
write.csv(tmp2,"../METAL/gs.csv",row.names=F,quote=F)

common=intersect(tmp1$probe,tmp2$probe)
tmp1=tmp1[which(tmp1$probe%in%common),]
tmp2=tmp2[which(tmp2$probe%in%common),]

tmp1$N=22774 
tmp2$N=17904
  
write.csv(tmp1,"../METAL/wielscher.csv",row.names=F,quote=F)
write.csv(tmp2,"../METAL/gs.csv",row.names=F,quote=F)


##################################
##### PROCESS RESULTS ############
##################################

a=as.data.frame(fread("/Cluster_Filespace/Marioni_Group/Rob/CRP/METAL/full.pat_sens1.txt"))
a2=as.data.frame(fread("/Cluster_Filespace/Marioni_Group/Rob/CRP/METAL/full.pat1.txt"))

names(a)[6]="p"
names(a2)[6]="p"

a$p=as.numeric(a$p)
a2$p=as.numeric(a2$p)

a=a[order(a$p),]
a2=a2[order(a2$p),]


###############################
#### METAL STEP ###############
###############################

## STDERR ##

cd /Cluster_Filespace/Marioni_Group/Rob/CRP/METAL/

metal
REMOVEFILTERS
SCHEME STDERR
USESTRAND OFF
GENOMICCONTROL OFF
AVERAGEFREQ OFF
MINMAXFREQ OFF
COLUMNCOUNTING LENIENT
WEIGHTLABEL     DONTUSECOLUMN
DEFAULTWEIGHT   20000
SEPARATOR COMMA
PVALUELABEL p.pat
MARKER probe
EFFECTLABEL coef.pat
STDERRLABEL se.pat
OUTFILE full.pat .txt

PROCESSFILE      gs.csv
PROCESSFILE      wielscher.csv

ANALYZE HETEROGENEITY
CLEAR



## SAMPLESIZE ##         ###### THIS IS THE PREFERRED METHOD BUT THE OTHER IS INCLUDED IN CASE IT IS USEFUL 

cd /Cluster_Filespace/Marioni_Group/Rob/CRP/METAL/
  
  metal
REMOVEFILTERS
SCHEME SAMPLESIZE
USESTRAND OFF
GENOMICCONTROL OFF
AVERAGEFREQ OFF
MINMAXFREQ OFF
COLUMNCOUNTING LENIENT
WEIGHTLABEL N
SEPARATOR COMMA
PVALUELABEL p.pat
MARKER probe
EFFECTLABEL coef.pat
STDERRLABEL se.pat
OUTFILE full.pat_sens .txt

PROCESSFILE      gs.csv
PROCESSFILE      wielscher.csv

ANALYZE HETEROGENEITY
CLEAR