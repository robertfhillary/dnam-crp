# Load requisite libraries 
library(ggplot2)
library(wesanderson)
library(data.table)
library(qqman)
library(artyfarty) 
library(RColorBrewer)

# Set working directory 
setwd("/Cluster_Filespace/Marioni_Group/Rob/CRP/")

#### Run on p12 for package compatability #####

######################
###### VC Plot #######
######################

# Read in variance components # I manually made this file from REML and BayesR results 
data1=read.csv("variance.csv")

# Variance component plot - Figure 2B (but easier to plot first)
data1$Trait=factor(data1$Trait,levels=c("Genetics","Methylation","Combined"))
pdf("Plots/variance_components.pdf",width=8,height=6)
ggplot(data1, aes(fill=Trait, y=(Variance*100), x=Trait)) + 
  geom_bar(position="stack", stat="identity") +
  geom_errorbar(ymin = (data1$LCI*100), ymax=(data1$HCI*100), width = 0.25) + 
  scale_y_continuous(expand=c(0,0),breaks = c(0,10,20,30,40,50,60,70,80,90,100), limits = c(0,100))  +
  scale_fill_manual(values = wes_palette("Cavalcanti1", n = 3))+ xlab("Omics Data Type") + ylab("Proportion of Variance Captured % \n [95% Confidence or Credible Interval]") + theme(legend.title.align = 0.5) + labs(fill = "Data Type") +
  theme_scientific() +  theme(axis.title=element_text(size=13),axis.text=element_text(size=12),legend.position="none") + facet_wrap(~Group) +  theme(strip.background = element_rect(fill="gray97"),strip.text.x = element_text(size = 10.5))
dev.off()



######################
## Manhattan Plot ####
######################
# Read in data 
tmp=as.data.frame(fread("/Cluster_Filespace/Marioni_Group/Rob/CRP/GS_EWAS/OSCA_Outputs/crp_standard_full.linear"))
# Fix gene information 
tmp$Gene=gsub(";.*","",tmp$Gene)
# Extract required columns
tmp=tmp[,c("Probe","Chr","bp","p")]
names(tmp)=c("SNP","CHR","BP","P")

pdf("/Cluster_Filespace/Marioni_Group/Rob/CRP/Plots/manhattan.pdf",width=8,height=6)
print(manhattan(tmp, main = "Fully-Adjusted Model", 
          col = c("steelblue1", "mediumvioletred"),suggestiveline=-log10(3.6e-8), genomewideline = F,annotatePval = 1e-60,annotateTop = T))
dev.off()


# Read in data 
tmp=as.data.frame(fread("/Cluster_Filespace/Marioni_Group/Rob/CRP/GS_EWAS/OSCA_Outputs/crp_standard_full.linear"))
# Fix gene information 
tmp$Gene=gsub(";.*","",tmp$Gene)
# Label column names as appropriate 
names(tmp)[1:3]=c("CHR","SNP","BP")
names(tmp)[8]="P"
# Subset to probes passing QC 
keep=as.data.frame(fread("/Cluster_Filespace/Marioni_Group/Elena/gs_osca/data/cpgs_tokeep.txt",header=F))
tmp=tmp[which(tmp$SNP %in% keep$V1),]

# Set up for plot - distance between CpGs 
tmp1=tmp[order(tmp$CHR,tmp$BP),]
don <- tmp1 %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(tmp1, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  
  # Add highlight and annotation information  %>%
  mutate( is_annotate=ifelse(-log10(P)>60, "yes", "no")) 
axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# Prepare label information
don[which(is.na(don$Gene)),"Gene"]="Intergenic"
don$Label=paste(don$SNP,don$Gene,sep=",")

tiff("/Cluster_Filespace/Marioni_Group/Rob/CRP/Plots/manhattan2.tiff",res=300,width=9.5,height=6.5,unit="in")
print(ggplot(don, aes(x=BPcum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("#4682B4", "#B47846"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0),limits=c(0,round((-log10(min(tmp$P))+3)))) +     # remove space between plot area and x axis
    xlab("CpG Position") +  ylab(expression(-log[10](italic(p)))) + ggtitle("Fully-Adjusted Model") + 
  # Add highlighted points

  # Add label using ggrepel to avoid overlapping
  geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=Label), size=4) +
  geom_hline(yintercept=7.44,linetype="dashed",size=1,color="seagreen3")+
  # Custom the theme:
  theme_scientific() + theme(axis.text=element_text(size=13),axis.title=element_text(size=14),plot.title = element_text(hjust = 0.5))+
  theme(legend.position="none"))
dev.off()