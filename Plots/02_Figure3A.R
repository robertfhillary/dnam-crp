# Set working directory 
setwd("/Cluster_Filespace/Marioni_Group/Rob/CRP/Correlations/")

# Load requisite libraries
library(data.table)
library(corrplot)
library(ggplot2)

# Read in results file 
cor=as.data.frame(matrix(nrow=70,ncol=5)) 

# Store Group Names 
cor[,1]=rep(c("Age 0","Age 7","Age 15", "Age 24", "Mothers","Age 40-70 (EUR)","Age 40-70 (SA)","Age 70","Age 73 - ls","Age 73 - hs","Age 76","Age 79","Age 87","Age 90"),each=5)

# Store Predictor Names 
cor[,2]=rep(c("Hillary-EWAS","Wielscher-EWAS","Wielscher-PCA","BayesR","Elnet"),14)

# ALSPAC
# age 0 (cord blood)
alspac.cord=read.csv("CORRELATIONS_ALSPAC_cord_2023-05-24.csv")
cor[1:5,3]=alspac.cord$r
cor[1:5,4]=alspac.cord$p
cor[1:5,5]="ALSPAC"

# age 7-9 
alspac.7=read.csv("CORRELATIONS_ALSPAC_age9_2023-05-24.csv")
cor[6:10,3]=alspac.7$r
cor[6:10,4]=alspac.7$p
cor[6:10,5]="ALSPAC"

# age 15 
alspac.15=read.csv("CORRELATIONS_ALSPAC_age15.450_2023-05-24.csv")
cor[11:15,3]=alspac.15$r
cor[11:15,4]=alspac.15$p
cor[11:15,5]="ALSPAC"

# age 24 
alspac.24=read.csv("CORRELATIONS_ALSPAC_age24_2023-05-24.csv")
cor[16:20,3]=alspac.24$r
cor[16:20,4]=alspac.24$p
cor[16:20,5]="ALSPAC"

# Mothers 
alspac.mothers=read.csv("CORRELATIONS_ALSPAC_mom_2023-05-24.csv")
cor[21:25,3]=alspac.mothers$r
cor[21:25,4]=alspac.mothers$p
cor[21:25,5]="ALSPAC"

# SABRE 
# Europeans
sabre.euro=read.csv("CORRELATIONS_SABRE_EUR_22062023.csv")
cor[26:30,3]=sabre.euro$r
cor[26:30,4]=sabre.euro$p
cor[26:30,5]="SABRE"

# South Asians
sabre.sa=read.csv("CORRELATIONS_SABRE_SA_22062023.csv")
cor[31:35,3]=sabre.sa$r
cor[31:35,4]=sabre.sa$p
cor[31:35,5]="SABRE"

## LBC1936 
# Wave 1
lbc36.w1=read.csv("CORRELATIONS_LBC36_W1_08052023.csv")
cor[36:40,3]=lbc36.w1$r
cor[36:40,4]=lbc36.w1$p
cor[36:40,5]="LBC1936"

# Wave 2 - lsCRP 
lbc36.w2.ls=read.csv("CORRELATIONS_LBC36_W2_lsCRP_08052023.csv")
cor[41:45,3]=lbc36.w2.ls$r
cor[41:45,4]=lbc36.w2.ls$p
cor[41:45,5]="LBC1936"

# Wave 2 - hsCRP (need to replace)
lbc36.w2.hs=read.csv("CORRELATIONS_LBC36_W2_08052023.csv")
cor[46:50,3]=lbc36.w2.hs$r
cor[46:50,4]=lbc36.w2.hs$p
cor[46:50,5]="LBC1936"

# Wave 3 
lbc36.w3=read.csv("CORRELATIONS_LBC36_W3_08052023.csv")
cor[51:55,3]=lbc36.w3$r
cor[51:55,4]=lbc36.w3$p
cor[51:55,5]="LBC1936"

# Wave 4
lbc36.w4=read.csv("CORRELATIONS_LBC36_W4_08052023.csv")
cor[56:60,3]=lbc36.w4$r
cor[56:60,4]=lbc36.w4$p
cor[56:60,5]="LBC1936"

## LBC1921
# Wave 3
lbc21.w3=read.csv("CORRELATIONS_LBC21_W3_08052023.csv")
cor[61:65,3]=lbc21.w3$r
cor[61:65,4]=lbc21.w3$p
cor[61:65,5]="LBC1921"

# Wave 4
lbc21.w4=read.csv("CORRELATIONS_LBC21_W4_08052023.csv")
cor[66:70,3]=lbc21.w4$r
cor[66:70,4]=lbc21.w4$p
cor[66:70,5]="LBC1921"


# Plot Stage
cor$Text=as.character(sprintf('%.2f',cor[1:70,3]))
#cor$Text=as.character(round(cor$Corr,2))
names(cor)=c("Test","Predictor","Corr","P","Cohort","Text")
cor$Predictor=factor(cor$Predictor,levels=c("Hillary-EWAS","Wielscher-EWAS","Wielscher-PCA","BayesR","Elnet"))
cor$Test=factor(cor$Test,levels=rev(c("Age 0","Age 7","Age 15","Age 24","Mothers","Age 40-70 (EUR)", "Age 40-70 (SA)", "Age 70", "Age 73 - ls", "Age 73 - hs", "Age 76", "Age 79", "Age 87", "Age 90")))
cor$Cohort=factor(cor$Cohort,levels=c("ALSPAC","SABRE","LBC1936","LBC1921"))

pdf("/Cluster_Filespace/Marioni_Group/Rob/CRP/Plots/heatmap_20062023.pdf",width=9,height=8)
p <- ggplot(cor, aes(x=Test, y=Predictor,fill=Corr)) +
  theme_bw() +
  theme(axis.line = element_line(colour = "gray40"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +  geom_tile(color="gray50") +
  scale_x_discrete(expand=c(0,0), name="Dataset") +
  scale_y_discrete(expand=c(0,0),name="Predictor",labels=c("Hillary-EWAS","Wielscher-EWAS","Wielscher-PCA","Bayesian PR","Elastic Net")) + geom_text(aes(label=Text,fontface=ifelse(cor$P<0.05,2,NA)))  
p1= p + scale_fill_gradient2(na.value="white",low="#D7191C", mid="white", high="#2C7BB6",limits=c(-1,1),name="Correlation Coefficient")  ## change the color theme ##
p1=p1 + theme(axis.text.y=element_text(size=11),axis.title=element_text(size=12), axis.text.x=element_text(size=11,angle=45,hjust=1))
p1=p1+coord_flip()+facet_grid(Cohort~.)+theme(legend.position="top")

print(p1 + facet_grid(Cohort~.,scales = "free", space = "free_y",switch="y")+ 
        theme(strip.placement="outside",panel.margin=unit(0, "lines"),                     # Place facet labels outside x axis labels.
              strip.text=element_text(size=11),strip.background = element_rect(fill = "white"),  # Make facet label background white.
              axis.title = element_blank()))
dev.off()

