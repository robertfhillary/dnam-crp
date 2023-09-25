# Load requisite libraries 
library("foreign")
library("data.table")
library("survival")
library("cowplot")
library("RColorBrewer")
library("ggplot2")
library("artyfarty")

# Read in LBC data 
ph36=read.spss("/Cluster_Filespace/Marioni_Group/Rob/CRP/LBC1936_EWAS_BldBasedChonicLowGradeInflammation_RH_27FEB2023.sav",to.data.frame=T)
# Prepare data 
ph36$sbp=(ph36$sbp1sit_w2+ph36$sbp2sit_w2)/2
ph36$dbp=(ph36$dbp1sit_w2+ph36$dbp2sit_w2)/2
ph36$age_w2=ph36$agedays_w2/365.25

# Select key variables for testing 
vars=c("lbc36no","age_w2","sex","smokcurr_w2","height_w2","fev_w2","fer_w2","pef_w2","fvc_w2","sbp","dbp","bld_hdlchol_w2","bld_choles_w2","bld_hba1c_DCCT_w2","bld_iron_w2","bld_creat_w2", "bmi_w2","alcunitwk_w2","bld_triglyc_w2","bld_apttr_w2","bld_ptrat_w2","bld_hb_w2","bld_rcc_w2","bld_hct_w2","bld_mcv_w2","bld_pltlt_w2","bld_prothr_w2","bld_aptt_w2","bld_fibrin_w2","bld_ferritin_w2","adl_w2","yrsedu_w1","depind_w1","sixmwk_w2","griplh_w2","griprh_w2","diab_w2","cvdhist_w2","stroke_w2","hichol_w2","hibp_w2")
ph36=ph36[,which(names(ph36)%in%vars),]
ph36$alcunitwk_w2=as.numeric(ph36$alcunitwk_w2)
ph36$smokcurr_w2=as.numeric(as.factor(ph36$smokcurr_w2))
names(ph36)[1]="ID"

# Merge with w2 scores 
scores=read.csv("/Cluster_Filespace/Marioni_Group/Rob/CRP/Projected_Scores/scores_LBC36W2.csv")
tar=readRDS("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/targets_3489_bloodonly.rds")
tar1=tar[tar$cohort=="LBC36" & tar$WAVE==2,]
tar1=tar1[,c("ID","Basename")]
scores=merge(scores,tar1,by="Basename",all.x=T)
df=merge(ph36,scores,by="ID")

# Merge in genetic scores
gen=read.table("/Cluster_Filespace/Marioni_Group/Rob/CRP/PGRS/lbc1936_genomewide_crp.all.score",header=T)
gen$FID=NULL
names(gen)=c("ID","genetic_score")
df=merge(df,gen,by="ID",all.x=T)

################################################################
### RUN LINEAR REGRESSIONS AGAINST CONTINUOUS OUTCOMES #########
################################################################

# Define the continuous variables
cont=c("sbp","dbp","bld_hdlchol_w2","bld_choles_w2","bld_hba1c_DCCT_w2","bld_creat_w2","bmi_w2","alcunitwk_w2","bld_triglyc_w2","bld_pltlt_w2","bld_prothr_w2","bld_aptt_w2","adl_w2","yrsedu_w1","depind_w1","sixmwk_w2","griplh_w2","griprh_w2")

# define lists to save output 
list1=list()
# Loop through each outcome and predictor 
for(j in c("genetic_score","Linear_Hillary","Wielscher","Wielscher_PCA","Elnet","BayesPR","CRP")){ 
  for(var in cont){ 
    list1[[j]][[var]]<- summary(lm(scale(df[,var]) ~ scale(df[,j]) + df$age_w2 + factor(df$sex)))$coefficients[2,]
  }
} 
# combine results 
l1 <- do.call('c', list1)
# tidy up results for convenience later in plotting 
l1=lapply(split(l1,sub('.*\\.', '', names(l1))),function(x) do.call(rbind, x))
l1=as.data.frame(do.call("rbind",l1))
l1$pred=gsub("\\..*", "", row.names(l1))
l1$trait=gsub(".*\\.", "", row.names(l1))



## FEV, FVC, PEF, FER will be done separately as height is needed as an additional covariate 
lung=c("fev_w2","fer_w2","pef_w2","fvc_w2")
# Set up list to store results 
list2=list()
# Loop through each outcome and predictor predictor 
for(j in c("genetic_score","Linear_Hillary","Wielscher","Wielscher_PCA","Elnet","BayesPR","CRP")){ 
  for(var in lung){ 
    list2[[j]][[var]]<- summary(lm(scale(df[,var]) ~ scale(df[,j]) + df$age_w2 + df$sex + df$height_w2))$coefficients[2,]
  }
} 
# combine results 
l2 <- do.call('c', list2)
# format results for later 
l2=lapply(split(l2,sub('.*\\.', '', names(l2))),function(x) do.call(rbind, x))
l2=as.data.frame(do.call("rbind",l2))
l2$pred=gsub("\\..*", "", row.names(l2))
l2$trait=gsub(".*\\.", "", row.names(l2))


###################################
#### PLOT CONTINUOUS OUTCOMES #####
###################################

# Combine results for continuous outcomes 
l3=as.data.frame(rbind(l1,l2))
names(l3)=c("Beta","SE","t","P","Predictor","Trait")
l3[l3$Predictor%in%"genetic_score","Predictor"]="Genetic Score"
l3[l3$Predictor%in%"Linear_Hillary","Predictor"]="Hillary-EWAS"
l3[l3$Predictor%in%"CRP","Predictor"]="Measured CRP"
l3[l3$Predictor%in%"Wielscher","Predictor"]="Wielscher-EWAS"
l3[l3$Predictor%in%"Wielscher_PCA","Predictor"]="Wielscher-PCA"
l3[l3$Predictor%in%"Elnet","Predictor"]="Elastic Net"
l3[l3$Predictor%in%"BayesPR","Predictor"]="Bayesian PR"
l3$FDR=p.adjust(l3$P,method="BH")
# Write out 
write.csv(l3, "/Cluster_Filespace/Marioni_Group/Rob/CRP/Results/continuous_phewas_09052023.csv",row.names=F)

###########################################################
### We will just plot Bayesian PR + Measured CRP for now ########
########################################################### 

# Plot Betas against each other for Bayesian PR vs CRP 
l3.p=l3 #for supp figure 
l3=l3[which(l3$Predictor %in% c("Measured CRP","Bayesian PR")),]
t=brewer.pal(n=7,"Dark2")
# Set variable names 
nms=c("Townsend Disability Scale","Alcohol Consumption", "Activated Partial Thrombopl. Time","Total Cholesterol", "Creatinine","HbA1c (DCCT)","HDL Cholesterol","Platelet Count","Prothrombin Time","Triglycerides","Body Mass Index","Diastolic Blood Pressure","Deprivation Index","Grip Strength - Left", "Grip Strength - Right","Systolic Blood Pressure","Six Metre Walk","Education-Years","Forced Expiratory Rate","Forced Expiratory Volume-1sec","Forced Vital Capacity","Peak Expiratory Flow")
l3$Trait=rep(nms,each=2)
l3$Predictor=factor(l3$Predictor,levels=c("Measured CRP","Bayesian PR"))
l3$Trait = factor(l3$Trait, levels=unique(l3$Trait[order(l3$Beta)]))
# Calculate Confidence Intervals 
l3$HCI=l3[,1]+(1.96*l3[,2])
l3$LCI=l3[,1]-(1.96*l3[,2])

###########################################################################
### FIGURE 3A - Comparison with Measured CRP for continuous outcomes ######
###########################################################################

pdf("/Cluster_Filespace/Marioni_Group/Rob/CRP/Plots/cont_betas.pdf",width=9,height=9) 
l3$Graph="Continuous Outcomes"
g=ggplot(l3,aes(y=Beta, x=Trait, group = factor(Predictor),color=factor(Predictor))) + 
  geom_point(position=position_dodge(width=0.7),shape=17,size = 3,alpha=0.9)+
  geom_errorbar(position=position_dodge(width=0.7),aes(ymin = LCI, ymax = HCI), width = 0.05,
                colour = "gray50") + 
  ylab("Standardised Beta [95% Confidence Interval]")+ xlab ("")+ theme_scientific() + scale_color_manual(values=t[c(1,6)],name="Predictor") + 
   scale_y_continuous(limits = c(-0.6,0.6), breaks = c(-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6))+ theme(axis.text=element_text(size=11),axis.title=element_text(size=12))+coord_flip()+
  geom_hline(yintercept = 0, linetype = "dotted")+facet_wrap(~Graph)+theme(strip.background = element_rect(fill = "gray95"),strip.text=element_text(size=11))+  guides(color = "none")
print(g)
dev.off()


############################
### SUPPLEMENTARY FIGURE ###
############################
# Prep data and colours 

t1=brewer.pal(n=7,"Dark2")
l3.p$log=-log10(l3.p$FDR)
l3.p$Trait = factor(l3.p$Trait, levels=unique(l3.p$Trait[order(l3.p$log)]))
# Set names of variables
# Set order of predictors for plot 
l3.p$Predictor=factor(l3.p$Predictor,levels=c("Measured CRP","Genetic Score","Hillary-EWAS","Wielscher-EWAS","Wielscher-PCA","Bayesian PR","Elastic Net"))
pdf("/Cluster_Filespace/Marioni_Group/Rob/CRP/Plots/p_values.pdf",width=7.5,height=13)
g.supp=ggplot(l3.p,aes(x=l3.p$Trait,y=l3.p$log,color=factor(Predictor)))+scale_y_continuous(expand=c(0.01,0.01),limits=c(0,20),breaks=c(0,5,10,15,20)) + ylab(expression(-log[10](italic(p)[FDR])))+ geom_hline(yintercept=1.30103,linetype="dashed",size=0.5,color="gray55")+
scale_x_discrete(expand = c(.05, .05),name="",labels=rev(c("Body Mass Index","HDL Cholesterol","Six Metre Walk","HbA1c (DCCT)","Triglycerides","Grip Strength - Right","Peak Expiratory Flow","Townsend Disability Scale","Grip Strength - Left","Creatinine","Activated Partial Thrombopl. Time","Education-Years","Forced Expiratory Volume-1sec","Forced Vital Capacity", "Platelet Count","Deprivation Index","Alcohol Consumption","Prothrombin Time","Forced Expiratory Rate","Total Cholesterol","Diastolic Blood Pressure","Systolic Blood Pressure"))) +
geom_point(size=3,shape=18,alpha=0.73) + coord_flip()+scale_color_manual(values=t1)+theme_scientific() +theme(legend.position="none") + facet_grid(Predictor~.)+ theme(strip.placement="outside",panel.spacing=unit(0.2, "lines"),                     # Place facet labels outside x axis labels.
                                                                          strip.text=element_text(size=9.7),strip.background = element_rect(fill = NA),  # Make facet label background white.
                                                               axis.title = element_blank())+theme(axis.text.y=element_text(size=6),axis.text.x=element_text(size=12),axis.title=element_text(size=12))
print(g.supp)
dev.off()

###################################################################
### RUN LOGISTIC REGRESSIONS AGAINST CATEGORICAL OUTCOMES #########
###################################################################
# Select categorical variables for testing 
cat=c("diab_w2","cvdhist_w2","stroke_w2","hichol_w2","hibp_w2") 

# define lists to save output 
df$stroke_w2=ifelse(df$stroke_w2%in%"Yes",1,0)
df$cvdhist_w2=ifelse(df$cvdhist_w2%in%"Yes",1,0)
df$hibp_w2=ifelse(df$hibp_w2%in%"Yes",1,0)
df$hichol_w2=ifelse(df$hichol_w2%in%"Yes",1,0)
df$diab_w2=ifelse(df$diab_w2%in%"Yes",1,0)

# Set up list to store output 
list3=list()
# Loop through each predictor and outcome 
for(j in c("genetic_score","Linear_Hillary","Wielscher","Wielscher_PCA","Elnet","BayesPR","CRP")){ 
  for(var in cat){ 
    list3[[j]][[var]]<- summary(glm(df[,var] ~ scale(df[,j]) + df$age_w2 + df$sex,family="binomial"))$coefficients[2,]
  }
} 
# combine 
l4 <- do.call('c', list3)
l4=lapply(split(l4,sub('.*\\.', '', names(l4))),function(x) do.call(rbind, x))
l4=as.data.frame(do.call("rbind",l4))
l4$pred=gsub("\\..*", "", row.names(l4))
l4$trait=gsub(".*\\.", "", row.names(l4))

# Extract stats 
l4$HCI=exp(l4[,1]+(1.96*l4[,2]))
l4$LCI=exp(l4[,1]-(1.96*l4[,2]))
l4$Odds=exp(l4[,1])

names(l4)=c("Beta","SE","Z","P","Predictor","Trait","HCI","LCI","Odds")
l4[l4$Predictor%in%"genetic_score","Predictor"]="Genetic Score"
l4[l4$Predictor%in%"Linear_Hillary","Predictor"]="Hillary-EWAS"
l4[l4$Predictor%in%"CRP","Predictor"]="Measured CRP"
l4[l4$Predictor%in%"Wielscher","Predictor"]="Wielscher-EWAS"
l4[l4$Predictor%in%"Wielscher_PCA","Predictor"]="Wielscher-PCA"
l4[l4$Predictor%in%"Elnet","Predictor"]="Elastic Net"
l4[l4$Predictor%in%"BayesPR","Predictor"]="Bayesian PR"
l4$FDR=p.adjust(l4$P,method="BH")
# Write out 
write.csv(l4, "/Cluster_Filespace/Marioni_Group/Rob/CRP/Results/categorical_phewas_09052023.csv",row.names=F)

# Subset to Bayesian PR and Measured CRP
l4=l4[l4$Predictor %in% c("Bayesian PR","Measured CRP"),]
l4$Trait=factor(l4$Trait,levels=c("hichol_w2","hibp_w2","cvdhist_w2","stroke_w2","diab_w2"))
l4$Predictor=factor(l4$Predictor,levels=c("Measured CRP","Bayesian PR"))
# Graph stage 
pdf("/Cluster_Filespace/Marioni_Group/Rob/CRP/Plots/logistic.pdf",width=8,height=5)
l4$Graph="Cross-Sectional Disease Outcomes"
g2=ggplot(l4,aes(y=Odds, x=Trait, group = factor(Predictor),color=factor(Predictor))) + 
  geom_point(position=position_dodge(width=0.5),shape = 17, size = 3,alpha=0.9)+
  geom_errorbar(position=position_dodge(width=0.5),aes(ymin = LCI, ymax = HCI), width = 0.05,
                colour = "gray50")+ scale_x_discrete(labels=c("High Cholesterol", "Hypertension", "Cardiovascular\n Disease", "Stroke","Type 2 Diabetes")) + 
  ylab("Odds Ratio [95% Confidence Interval]")+ xlab ("")+ theme_scientific() + scale_color_manual(values=t[c(1,6)],name="Predictor") + 
  geom_hline(yintercept = 1, linetype = "dotted")+ scale_y_continuous(limits = c(0,3), breaks = c(0,0.5,1,1.5,2,2.5,3))+ theme(axis.text=element_text(size=11),axis.title=element_text(size=12))+coord_flip()+
   facet_wrap(~Graph)+theme(strip.background = element_rect(fill = "gray95"),strip.text=element_text(size=11)) +  guides(color = "none")
dev.off() 



###############################################
########## ALL-CAUSE MORTALITY MODELS #########
###############################################
# Hazard Ratios 
## Read in death information
dead=read.spss("/Cluster_Filespace/Marioni_Group/Rob/GrimAge/LBC1936_DNA_MethylationBasedEstimatorOfTeloLength_RM_28FEB2019.sav",to.data.frame=T)
# Code death events 
dead$event = 0
dead$event[dead$dead=="DEAD"] <- 1
# Set age information 
dead$age_event = ifelse(dead$event==0, dead$AgedaysApx_LastCensor/365.25, dead$agedays_death/365.25)
dead$age2 = dead$agedays_w2/365.25 # would need to change this to agedays_w1 if proceeding with wave1 scores for longer followup period
# Merge in information with scores 
df=merge(df,dead[,c("event","age_event","lbc36no","age2")],by.x="ID",by.y="lbc36no",all.y=T)
# Read in Wave 1 scores for convenience (if needed later)
#df1=read.csv("../Projected_Scores/scores_LBC36W1.csv")
#names(df1)[2:ncol(df1)]=paste0(names(df1)[2:ncol(df1)],"_w1")
# Merge in W1 IDs so we can link them to their wave2 samples 
#tar2=tar[tar$cohort=="LBC36" & tar$WAVE ==1,]
#tar2=tar2[,c("Basename","ID")]
#df1=merge(tar2,df1,by="Basename")
# Combine with Wave2 scores 
#df=merge(df,df1,by="ID",all.y=T) 
# Calculate information 
df$tte <- df$age_event - df$age2

# ALTHOUGH WE MERGED IN WAVE1 DATA - WE WILL PROCEED WITH WAVE2 - THE PATTERNS OF RESULTS ARE SAME IN BOTH WAVES AND WAVE2 HAS HSCRP
# Loop for hazard ratios 
list4=list()
# Loop through each predictor and outcome 
for(j in c("genetic_score","Linear_Hillary","Wielscher","Wielscher_PCA","Elnet","BayesPR","CRP")){ 
    list4[[j]]<- summary(coxph(Surv(tte, event) ~ age2 + factor(sex) + scale(df[,j]), data=df))$coefficients[3,]
  }
# combine 
l5 <- as.data.frame(do.call('rbind', list4))
l5$Predictor=gsub("\\..*", "", row.names(l5))
l5$Trait="All-Cause Mortality"

# Extract stats 
l5$HCI=exp(l5[,1]+(1.96*l5[,3]))
l5$LCI=exp(l5[,1]-(1.96*l5[,3]))
l5$Odds=exp(l5[,1])


l5[,2]=NULL
names(l5)=c("Beta","SE","Z","P","Predictor","Trait","HCI","LCI","HR")
l5[l5$Predictor%in%"genetic_score","Predictor"]="Genetic Score"
l5[l5$Predictor%in%"Linear_Hillary","Predictor"]="Hillary-EWAS"
l5[l5$Predictor%in%"CRP","Predictor"]="Measured CRP"
l5[l5$Predictor%in%"Wielscher","Predictor"]="Wielscher-EWAS"
l5[l5$Predictor%in%"Wielscher-PCA","Predictor"]="Wielscher-PCA"
l5[l5$Predictor%in%"Elnet","Predictor"]="Elastic Net"
l5[l5$Predictor%in%"BayesPR","Predictor"]="Bayesian PR"
l5$FDR=p.adjust(l5$P,method="BH")
# Write out 
write.csv(l5, "/Cluster_Filespace/Marioni_Group/Rob/CRP/Results/survival_phewas_09052023.csv",row.names=F)

# Subset to Bayesian PR and Measured CRP
l5=l5[l5$Predictor %in% c("Bayesian PR","Measured CRP"),]
l5$Predictor=factor(l5$Predictor,levels=c("Measured CRP","Bayesian PR"))
# Graph stage 
pdf("/Cluster_Filespace/Marioni_Group/Rob/CRP/Plots/hazards.pdf",width=8,height=5)
l5$Graph="All-Cause Mortality (Longitudinal)"
g3=ggplot(l5,aes(y=HR, x=Trait, group = factor(Predictor),color=factor(Predictor))) + 
  geom_point(position=position_dodge(width=0.5),shape = 17, size = 3,alpha=0.9)+
  geom_errorbar(position=position_dodge(width=0.5),aes(ymin = LCI, ymax = HCI), width = 0.05,
                colour = "gray50")+
  ylab("Hazard Ratio [95% Confidence Interval]")+ xlab ("")+ theme_scientific() + scale_color_manual(values=t[c(1,6)],name="Predictor") + 
  geom_hline(yintercept = 1, linetype = "dotted")+ scale_y_continuous(limits = c(0,3), breaks = c(0,0.5,1,1.5,2,2.5,3))+ theme(axis.text=element_text(size=11),axis.title=element_text(size=12))+coord_flip()+
  facet_wrap(~Graph)+theme(strip.background = element_rect(fill = "gray95"),strip.text=element_text(size=11)) +  guides(color = "none")
print(g3)
dev.off() 


# Get legend 
graph_legend=ggplot(l5,aes(y=HR, x=Trait, group = factor(Predictor),color=factor(Predictor))) + 
  geom_point(position=position_dodge(width=0.5),shape = 17, size = 3,alpha=0.9)+
  geom_errorbar(position=position_dodge(width=0.5),aes(ymin = LCI, ymax = HCI), width = 0.05,
                colour = "gray50")+
  ylab("Hazard Ratio [95% Confidence Interval]")+ xlab ("")+ theme_scientific() + scale_color_manual(values=t,name="Predictor") + 
  geom_hline(yintercept = 1, linetype = "dotted")+ scale_y_continuous(limits = c(0,3), breaks = c(0,0.5,1,1.5,2,2.5,3))+ theme(axis.text=element_text(size=11),axis.title=element_text(size=12))+coord_flip()+
  facet_wrap(~Graph)+theme(strip.background = element_rect(fill = "gray95"),strip.text=element_text(size=11)) +  guides(color = guide_legend(reverse=TRUE)) + theme(legend.title = element_text(hjust = 0.5))
legend=get_legend(graph_legend)


# Combined plot
pdf("/Cluster_Filespace/Marioni_Group/Rob/CRP/Plots/Figure4_combined_bayespr.pdf",width=14,height=9)
plot = ggdraw(plot_grid(g, plot_grid(g2, g3, nrow = 2, rel_widths = c(1/4,1/4)),legend, ncol = 3,rel_widths = c(0.75,1/2,0.3)))
print(plot)
dev.off()
