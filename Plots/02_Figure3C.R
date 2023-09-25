# Load requisite libraries 
library(artyfarty)
library(ggplot2)
library(tidyverse)
library(shades)
library(RColorBrewer)

# Collate results 
setwd("/Cluster_Filespace/Marioni_Group/Rob/CRP/Thresholds_v2/Scores/")
# Extract file 
loop=list.files(".", ".")
# Set up dataframe to store correlation coefficients 
mat=as.data.frame(matrix(nrow=length(loop),ncol=14))
names(mat)=c("Threshold","n","cor_pred_CRP","cor_PC_CRP","cor_Elnet_CRP","LCI_pred_CRP","LCI_PC_CRP","LCI_Elnet_CRP","HCI_pred_CRP","HCI_PC_CRP","HCI_Elnet_CRP","p_pred_CRP","p_PC_CRP","p_Elnet_CRP")
# Loop through each file to store correlation coefficients 
for(i in loop){ 
  # Read in each and calculate correlations files and extract key information 
  res=read.csv(i)
  # Store predictor as the name of the dataframe 
  nm=gsub("Predictors_", "", i)
  nm=as.numeric(gsub(".csv","",nm))
  # Get index of file in loop 
  ind=which(loop %in% i)
  # Store p-value threshold 
  mat[ind,1]=nm
  # Store the number of CpGs 
  pred=read.csv(paste0("../../Thresholds/Predictors/",i))
  mat[ind,2]=length(unique(pred$CpG))
  # Store correlation coefficient between EWAS-predictor and CRP 
  mat[ind,3]=cor.test(res$CRP,res$Wielscher_CRP)$estimate[[1]]
  # Store correlation coefficient between PCA-predictor and CRP 
  mat[ind,4]=cor.test(res$CRP,res$PC.CRP)$estimate[[1]]
  # Store correlation coefficient between Elnet-predictor and CRP 
  mat[ind,5]=cor.test(res$CRP,res$Elnet_CRP)$estimate[[1]]
  # Store LCI for EWAS-predictor and CRP 
  mat[ind,6]=cor.test(res$CRP,res$Wielscher_CRP)$conf.int[[1]]
  # Store LCI for PCA-predictor and CRP 
  mat[ind,7]=cor.test(res$CRP,res$PC.CRP)$conf.int[[1]]
  # Store LCI for Elnet-predictor and CRP 
  mat[ind,8]=cor.test(res$CRP,res$Elnet_CRP)$conf.int[[1]]
  # Store HCI for EWAS-predictor and CRP 
  mat[ind,9]=cor.test(res$CRP,res$Wielscher_CRP)$conf.int[[2]]
  # Store HCI for PCA-predictor and CRP 
  mat[ind,10]=cor.test(res$CRP,res$PC.CRP)$conf.int[[2]]
  # Store HCI for Elnet-predictor and CRP
  mat[ind,11]=cor.test(res$CRP,res$Elnet_CRP)$conf.int[[2]]
  # Store p value for EWAS-predictor and CRP 
  mat[ind,12]=cor.test(res$CRP,res$Wielscher_CRP)$p.value[1]
  # Store p value for PCA-predictor and CRP 
  mat[ind,13]=cor.test(res$CRP,res$PC.CRP)$p.value[1]
  # Store p value for Elnet-predictor and CRP 
  mat[ind,14]=cor.test(res$CRP,res$Elnet_CRP)$p.value[1]
  # Print to denote completion
  print(paste0("Completed ", ind," of ", length(loop)))
}   


# Order thresholds by p-value
mat=mat[order(mat$Threshold),]
# Convert data wide to long
# correlation coefficient 
correlations <- gather(mat[,c("Threshold","n","cor_pred_CRP","cor_PC_CRP","cor_Elnet_CRP")], Predictor, r, cor_Elnet_CRP:cor_pred_CRP, factor_key=TRUE)
# LCI 
lci <- gather(mat[,c("Threshold","LCI_pred_CRP","LCI_PC_CRP","LCI_Elnet_CRP")], Predictor, LCI, LCI_Elnet_CRP:LCI_pred_CRP, factor_key=TRUE)
# HCI 
hci <- gather(mat[,c("Threshold","HCI_pred_CRP","HCI_PC_CRP","HCI_Elnet_CRP")], Predictor, HCI, HCI_Elnet_CRP:HCI_pred_CRP, factor_key=TRUE)
# p value
p <-  gather(mat[,c("Threshold","p_pred_CRP","p_PC_CRP","p_Elnet_CRP")], Predictor, p, p_Elnet_CRP:p_pred_CRP, factor_key=TRUE)
# Combine metrics 
gs.plot1=cbind(correlations,lci[,c("LCI")],hci[,"HCI"],p[,"p"])
# Correct names 
names(gs.plot1)=c("Threshold","n","Predictor","r","LCI","HCI", "p")
# Format table 
gs.plot1$p=signif(gs.plot1$p,1)
# Write out results
write.csv(gs.plot1, "/Cluster_Filespace/Marioni_Group/Rob/CRP/Tables/p_value_thresholding.csv",row.names=F)

# Extract colors from Figure B 
t=brewer.pal(n=5,"Set3")
t=t[c(2,1,3,4,5)]
t=c("#D3D3D3",t)
# Subset to predictors being examined here 
t=t[c(3,4,6)]
# Set order of predictors 
gs.plot1$Predictor=factor(gs.plot1$Predictor,levels=c("cor_pred_CRP","cor_PC_CRP","cor_Elnet_CRP"))

# Plot data 
# without confidence intervals 
pdf("/Cluster_Filespace/Marioni_Group/Rob/CRP/Plots/pvalue_thresholds_withintervals.pdf",width=11.8,height=8)
ggplot(data=gs.plot1, aes(x= factor(Threshold), y= r)) +
  geom_point(aes(color=factor(Predictor)), size=4,position=position_dodge(width = 0.90)) +
  theme(legend.position="right") + 
  labs(y=expression(paste("Correlation with log-CRP levels (", italic("r"), ")")), x=expression(paste(italic(p),"-value threshold (increasing significance)"))) + scale_color_manual("Predictor Type",labels=rev(c("Elastic Net","Wielscher-PCA","Wielscher-EWAS")), values=brightness(t,0.8)) + 
  scale_x_discrete(limits = rev(levels(factor(gs.plot1$Threshold))), labels = paste0(rev(levels(factor(gs.plot1$Threshold))), "\n", "(n =", rev(gs.plot1$n), ")")) +
  scale_y_continuous(limits = c(-0.1,0.55),breaks=c(-0.10,-0.05,0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55))+  theme_scientific()+theme(axis.text.x = element_text(size=10.4),axis.title.x=element_text(size=14), axis.text.y=element_text(size=12),axis.title.y=element_text(size=14),legend.text=element_text(size=12),legend.title=element_text(size=13.5))
dev.off()



# with confidence intervals 
pdf("/Cluster_Filespace/Marioni_Group/Rob/CRP/Plots/pvalue_thresholds_withintervals.pdf",width=11.8,height=8)
ggplot(data=gs.plot1, aes(x= factor(Threshold), y= r)) +
  geom_point(aes(color=factor(Predictor)), size=4,position=position_dodge(width = 0.90)) +
  theme(legend.position="right") + 
  geom_errorbar(aes(ymin=LCI,ymax=HCI,color=factor(Predictor)),width=0,position=position_dodge(width = 0.90)) + 
  labs(y=expression(paste("Correlation with log-CRP levels (", italic("r"), ")")), x=expression(paste(italic(p),"-value threshold (increasing significance)"))) + scale_color_manual("Predictor Type",labels=rev(c("Elastic Net","Wielscher-PCA","Wielscher-EWAS")), values=brightness(t,0.8)) + 
  scale_x_discrete(limits = rev(levels(factor(gs.plot1$Threshold))), labels = paste0(rev(levels(factor(gs.plot1$Threshold))), "\n", "(n =", rev(gs.plot1$n), ")")) +
  scale_y_continuous(limits = c(-0.1,0.55),breaks=c(-0.10,-0.05,0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55))+  theme_scientific()+theme(axis.text.x = element_text(size=10.4),axis.title.x=element_text(size=14), axis.text.y=element_text(size=12),axis.title.y=element_text(size=14),legend.text=element_text(size=12),legend.title=element_text(size=13.5))
dev.off()




