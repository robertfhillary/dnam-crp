#####################################
#### INCREMENTAL R2 MODELLING #######
#####################################

# Set working directory 
setwd("/Cluster_Filespace/Marioni_Group/Rob/CRP/Projected_Scores/")

# Load requisite libraries
library(data.table)
library(corrplot)
library(RColorBrewer)
library(artyfarty)

# Read in dataset
w2=read.csv("scores_LBC36W2.csv")
# Merge in genetic data 
prs=read.table("../PGRS/lbc1936_genomewide_crp.all.score",header=T)
tar=read.csv("../../../LBC/LBC_methylation/target_QC_age_sex_date.csv")
tar1=tar[tar$cohort=="LBC36" & tar$WAVE ==2,]
prs1=merge(prs,tar1[,c("ID","Basename","age","sex")],by.x="IID",by.y="ID")
names(prs1)[3]="genetic_score"
names(prs1)[5]="age_w2"
w2=merge(w2,prs1,by="Basename")

# get variance by age and sex and also by age, sex and generations
hs.rg=summary(lm(w2$CRP~w2$age_w2+w2$sex+w2$genetic_score))$r.squared
hs.r=summary(lm(w2$CRP~w2$age_w2+w2$sex))$r.squared

# Loop through each predictor 
loop=c("Wielscher","Linear_Hillary","Wielscher_PCA","Elnet","BayesPR")
# Loop stage to get variance explained by epigenetic predictor and then also epigenetic+genetic predictor 
output=as.data.frame(matrix(nrow=length(loop),ncol=3))
names(output)=c("Predictor","R2_over_basic","R2_with_genetics")
for(i in loop){
  ind=which(loop%in%i)
  tmp.r=summary(lm(w2$CRP~w2$age_w2+w2$sex+w2[,i]))$r.squared
  tmp.r1=summary(lm(w2$CRP~w2$age_w2+w2$sex+w2$genetic_score+w2[,i]))$r.squared
  output[ind,1]=as.character(i)
  output[ind,2]=tmp.alone=tmp.r-hs.r  
  output[ind,3]=tmp.over.g=tmp.r1-hs.r
}
# Format output 
output[,1]=c("Wielscher-EWAS","Hillary-EWAS","Wielscher-PCA","Elastic Net", "Bayesian PR")


# Now we will tidy up the output for plotting 
varplot <- gather(output, type, measurement, R2_over_basic:R2_with_genetics, factor_key=TRUE)
varplot$measurement=varplot$measurement*100
varplot[11,1]="Polygenic"
varplot[11,2]="R2_over_basic"
varplot[11,3]=(hs.rg*100)
varplot[1:10,4]="DNA Methylation and Polygenic Scores"
varplot[11,4]="Polygenic Score Alone"
names(varplot)[4]="facet"

# Set up output for plot some more
varplot$Predictor=factor(varplot$Predictor,levels=c("Polygenic","Hillary-EWAS","Wielscher-EWAS","Wielscher-PCA","Bayesian PR","Elastic Net"))
varplot$type=ifelse(varplot$type %in% "R2_over_basic","DNAm Alone","DNAm and Polygenic")
varplot$type=factor(varplot$type,levels=rev(c("DNAm Alone","DNAm and Polygenic")))
t=brewer.pal(n=5,"Set3")
t=t[c(2,1,3,4,5)]
t=c("#D3D3D3",t)

# Plotting stage
pdf("/Cluster_Filespace/Marioni_Group/Rob/CRP/Plots/incrementalr2.pdf",width=10,height=8)
p1=ggplot(data=varplot,aes(x=type,y=measurement,fill=varplot$Predictor))+geom_bar(stat="identity", position = "dodge",aes(alpha=factor(varplot$type))) +   scale_alpha_manual(values = c(1,0.6)) + scale_y_continuous(expand=c(0,0), limits=c(0,30),breaks=c(0,5,10,15,20,25,30))+ ylab(expression("Incremental"~italic(R)^2~"over null model"))+coord_flip()+theme_scientific()
p2=p1+ xlab("") +geom_text(aes(label=paste(sprintf('%.1f',measurement),"%")), position=position_dodge(width=0.9), hjust=-0.25) +facet_grid(Predictor~.,space="free",scales="free") + theme(legend.position="none") + scale_fill_manual(values=t)
print(p2+ theme(strip.placement="outside",panel.margin=unit(0, "lines"),                     # Place facet labels outside x axis labels.
                strip.text=element_text(size=11),strip.background = element_rect(fill = "white"),  # Make facet label background white.
                axis.title = element_blank())+theme(axis.text=element_text(size=13.5),axis.title=element_text(size=14)))
dev.off()                                                                                                                                                                   
