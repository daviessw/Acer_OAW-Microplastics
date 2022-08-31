# source("http://bioconductor.org/biocLite.R")
# biocLite("DESeq2")

#set your working directory
setwd("~/Dropbox/Acer_Microplastics") 

###conduct array quality metrics to detect and remove outliers
library(DESeq2)
library(vegan)
library(ggplot2)

#read in counts 
countData <- read.table("acer_genome_2022counts.txt")
head(countData)
length(countData[,1])
#21569

conditions <- read.csv("samples_renamed.csv")
row.names(conditions) <- conditions$sample

#### making conditions data frame ####
names(countData)=sub("_L001_R1_001.fastq.trim.sam.counts","",names(countData))
names(countData)=sub("X","",names(countData))#get rid of the X
names(countData)=substr(names(countData),1,nchar(names(countData))-4)
names(countData)

#### count data ####
totalCounts=colSums(countData)
totalCounts
barplot(totalCounts, ylab="raw counts", xlab="sample number")
#     10     11     12     14     15     16     17     19     20     22     23     24     25 
# 834588 846094 609277 828973 648174 880982 414581 800716 863849 929270 780221 761130 656062 
# 26     28     29      2     32     34     35      3      6      8      9 
# 664303 604654 166828 688099 844935 758436 521313 818571 685167 966574 674017 
min(totalCounts) #166828
max(totalCounts)  #966574

#### remove sample 29 - bad # of reads ####
#skip these next 2 lines if you want to look at the data with sample 29
countData <- countData[,c(1:15,17:24)]#sample 29 is in column 16
conditions <- subset(conditions,sample!=29)

##remove the one replicate from each treatment of the clonal pair
countData <- countData[,c(1:7,9:13,15,17:20,22:23)]
remove <- c("2","19","26","6")
conditions <- conditions[!conditions$sample %in% remove,]

#### re-making g data frame ####
g=data.frame(conditions$genet, conditions$treat)
g
#colData=g

#### DESeq time! ####
dds<-DESeqDataSetFromMatrix(countData=countData, colData=g, design=~conditions.genet+conditions.treat) #can only test for the main effects of site, pco2, temp

#one step DESeq
dds<-DESeq(dds)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing

head(dds)
res<- results(dds)

########### PCA ###############
rlog=rlogTransformation(dds, blind=TRUE) 
rld=assay(rlog)
rld_t=t(rld)
colnames(rld_t)
head(rld_t)

pca <- prcomp(rld_t,center = TRUE)
head(pca)
li <- pca$sdev^2 / sum(pca$sdev^2)
pc1v <- round(li[1] * 100, 1)
pc2v <- round(li[2] * 100, 1)
pca_s <- as.data.frame(pca$x)
head(pca_s)
pca_s <- pca_s[,c(1,2)]
pca_s$Samples = row.names(pca_s)
pca_s$treatment=g$conditions.treat
pca_s$genet=g$conditions.genet
head(pca_s)

cbPalette <- c("dodgerblue", "darkorange","firebrick2", "firebrick4")
pdf("PCA_allgenes_rlog.pdf",height=4,width=5, useDingbats = FALSE)
ggplot(pca_s, aes(PC1, PC2, fill=treatment, shape=genet, group=treatment, label = rownames(pca_s))) +
  guides(color=guide_legend(title="Treatment"))+
  guides(pch=guide_legend(title="Genotype"))+
  geom_point(color = "black", size=4) +
  scale_fill_manual(values=cbPalette)+
  scale_colour_manual(values=cbPalette)+
  scale_shape_manual(values=c(21,22,23,24,25))+
  theme_bw() +
  stat_ellipse(aes(x=PC1, y=PC2, group=treatment, colour = treatment))+
  xlab(paste0("PC1: ",pc1v,"% variance")) +
  ylab(paste0("PC2: ",pc2v,"% variance")) + 
  guides(fill = "none")
dev.off()

adonis(pca_s[,1:2] ~ treatment+genet, data = pca_s, method='eu', na.rm = TRUE)
##only genet significant when using all genes
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# treatment  3     836.1  278.71  2.0342 0.15687  0.115    
# genet      4    2986.9  746.73  5.4500 0.56038  0.001 ***
#   Residuals 11    1507.1  137.01         0.28276           
# Total     18    5330.2                 1.00000

# #######Establishing top DEGs#######
top<-head(res[order(res$pvalue),], 1000)
top=data.frame(top)
nrow(top)
top.counts <- countData[c(row.names(countData)) %in% c(row.names(top)),]
nrow(top)
top.m <- as.matrix(top.counts)
rlog.top=rlogTransformation(top.m, blind=TRUE) 
rld_t=t(rlog.top)
colnames(rld_t)
head(rld_t)
length(rld_t)

pca <- prcomp(rld_t,center = TRUE)
head(pca)

li <- pca$sdev^2 / sum(pca$sdev^2)
pc1v <- round(li[1] * 100, 1)
pc2v <- round(li[2] * 100, 1)
pca_s <- as.data.frame(pca$x)
head(pca_s)
pca_s <- pca_s[,c(1,2)]
pca_s$Samples = row.names(pca_s)
pca_s$treatment=g$conditions.treat
pca_s$genet=g$conditions.genet
head(pca_s)

cbPalette <- c("dodgerblue", "darkorange","firebrick2", "firebrick4")
pdf("PCA_top1000_rlog.pdf",height=4,width=5, useDingbats = FALSE)
ggplot(pca_s, aes(PC1, PC2, fill=treatment, shape=genet, group=treatment, label = rownames(pca_s))) +
  guides(color=guide_legend(title="Treatment"))+
  guides(pch=guide_legend(title="Genotype"))+
  geom_point(color = "black", size=4) +
  scale_fill_manual(values=cbPalette)+
  scale_colour_manual(values=cbPalette)+
  scale_shape_manual(values=c(21,22,23,24,25))+
  theme_bw() +
  stat_ellipse(aes(x=PC1, y=PC2, group=treatment, colour = treatment))+
  xlab(paste0("PC1: ",pc1v,"% variance")) +
  ylab(paste0("PC2: ",pc2v,"% variance")) + 
  guides(fill = "none")
dev.off()

adonis(pca_s[,1:2] ~ treatment+genet, data = pca_s, method='eu', na.rm = TRUE)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# treatment  3    554.63 184.876 17.1988 0.57879  0.001 ***
#   genet      4    285.39  71.347  6.6374 0.29782  0.001 ***
#   Residuals 11    118.24  10.749         0.12339           
# Total     18    958.26                 1.00000  

####GE plasticity analysis
library(dplyr) # necessary for PCA function
library(ggfortify) # plotting the PCA
library(vegan) # running the PERMANOVA (adonis2())
library(ggpubr) # for arranging multiple plots into single figure
source("PCAplast_function.R") # source the plasticity function
conds=pca_s[,3:5]
str(conds)
## To run the plasticity function, enter the following objects:
plast_out <- PCAplast(pca = pca, # the PCA dataframe containing the PCA eigenvalues
                      data = conds, # the condition/treatment data corresponding to samples
                      sample_ID = "Samples", # the name of column that provide unique ID per sample (if blank, will pull rownames for this)
                      num_pca =  "2", # the number of PCAs to include in analysis (default is 'all', but you can specify another number with a minimum of 2 PCAs)
                      control_lvl = "AMB_Control", # control level of the treatment. If blank, a control mean per control level is assumed
                      group = "genet", # the grouping column (i.e., colony). If blank, will assume control level grouping only!
                      control_col = "treatment") # what the 'treatment' column is called
head(plast_out)
## Plot the plasticity (PC distances): overlay mean and 1 standard deviation
cbPalette2 <- c("darkorange","firebrick2", "firebrick4")
pdf("plasticity_top1000.pdf",height=4,width=3, useDingbats = FALSE)
plast_plot <- ggplot(data = plast_out, aes(x = treatment, y = dist, fill=treatment, shape=genet, group=treatment)) + 
  geom_point(color = "black", size=4, alpha = 0.3, position = position_jitter(width = 0.15)) +
  scale_fill_manual(values=cbPalette2)+
  scale_colour_manual(values=cbPalette2)+
  scale_shape_manual(values=c(21,22,23,24,25))+
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0, colour = c("darkorange","firebrick2", "firebrick4")) +
  stat_summary(fun = "mean", size = 0.5, colour = c("darkorange","firebrick2", "firebrick4")) +
  theme_bw() +
  ylab("Plasticity +/- SD")+
  xlab("Treatment")+
  theme(legend.position = "none")
plast_plot
dev.off()

aov1=aov(dist~treatment, data=plast_out)
summary(aov1)
# Df Sum Sq Mean Sq F value Pr(>F)  
# treatment    2  98.85   49.43   4.349 0.0406 *
#   Residuals   11 125.03   11.37  
TukeyHSD(aov1)
# $treatment
# diff        lwr       upr     p adj
# OAW_Control-AMB_MP -0.5304615 -6.6386956  5.577773 0.9702090
# OAW_MP-AMB_MP       5.2923110 -0.4665873 11.051209 0.0723525
# OAW_MP-OAW_Control  5.8227725 -0.2854615 11.931007 0.0619245


#################### OAW pairwise comparisons
g
g$conditions.treat<-factor(g$conditions.treat, levels=c("OAW_Control","AMB_Control"))
##second term is the "control"
resOAW <- results(dds, contrast=c("conditions.treat","OAW_Control","AMB_Control"))
#how many FDR < 10%?
table(resOAW$padj<0.05)
# 0.1=36
# 0.05=20
# 0.01=6
summary(resOAW)
# out of 21327 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 18, 0.084%
# LFC < 0 (down)     : 18, 0.084%
# outliers [1]       : 0, 0%
# low counts [2]     : 764, 3.6%
# (mean count < 0)

nrow(resOAW[resOAW$padj<0.05 & !is.na(resOAW$padj),])  # Num significantly differentially expressed genes excluding the no/low count genes   #228
#20
plotMA(resOAW, main="OAW_Control vs AMB_Control")

results <- as.data.frame(resOAW)
head(results)

nrow(resOAW[resOAW$padj<0.1 & resOAW$log2FoldChange > 0 & !is.na(resOAW$padj),])
nrow(resOAW[resOAW$padj<0.1 & resOAW$log2FoldChange < 0 & !is.na(resOAW$padj),])
#UP in OAW 18
#DOWN in OAW 18

write.table(resOAW, file="OAW.txt", quote=F, sep="\t")

resOAW=data.frame(resOAW)
head(resOAW)
go_OAW = resOAW %>%
  tibble::rownames_to_column(var = "iso") %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(iso, mutated_p_updown)
head(go_OAW)
colnames(go_OAW) <- c("gene", "pval")
head(go_OAW)
write.csv(go_OAW, file="go_OAW.csv", quote=F, row.names=FALSE)

##########################MP Pairwise Comparison
conditions <- read.csv("samples_renamed.csv")
row.names(conditions) <- conditions$sample
conditions <- subset(conditions,sample!=29)
remove <- c("2","19","26","6")
conditions <- conditions[!conditions$sample %in% remove,]
g=data.frame(conditions$genet, conditions$treat)
g
g=data.frame(conditions$genet, conditions$treat)
g$conditions.treat<-factor(g$conditions.treat, levels=c("AMB_MP","AMB_Control"))
##second term is the "control"
resMP <- results(dds, contrast=c("conditions.treat","AMB_MP","AMB_Control"))
#how many FDR < 10%?
table(resMP$padj<0.05)
# 0.1=57
# 0.05=17
# 0.01=5
summary(resMP)
# out of 21327 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 48, 0.23%
# LFC < 0 (down)     : 9, 0.042%
# outliers [1]       : 0, 0%
# low counts [2]     : 19835, 93%
# (mean count < 102)

nrow(resMP[resMP$padj<0.05 & !is.na(resMP$padj),])  # Num significantly differentially expressed genes excluding the no/low count genes   #228
#17
plotMA(resMP, main="AMB_Microplastics vs AMB_Control")

results <- as.data.frame(resMP)
head(results)

nrow(resMP[resMP$padj<0.1 & resMP$log2FoldChange > 0 & !is.na(resMP$padj),])
nrow(resMP[resMP$padj<0.1 & resMP$log2FoldChange < 0 & !is.na(resMP$padj),])
#UP in MP 48
#DOWN in MP 9

write.table(resMP, file="MP.txt", quote=F, sep="\t")

resMP=data.frame(resMP)
head(resMP)
go_MP = resMP %>%
  tibble::rownames_to_column(var = "iso") %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(iso, mutated_p_updown)
head(go_MP)
colnames(go_MP) <- c("gene", "pval")
head(go_MP)
write.csv(go_MP, file="go_MP.csv", quote=F, row.names=FALSE)

#################Pairwise Comparison OAW+MP
conditions <- read.csv("samples_renamed.csv")
row.names(conditions) <- conditions$sample
conditions <- subset(conditions,sample!=29)
remove <- c("2","19","26","6")
conditions <- conditions[!conditions$sample %in% remove,]
g=data.frame(conditions$genet, conditions$treat)
g
g=data.frame(conditions$genet, conditions$treat)
g$conditions.treat<-factor(g$conditions.treat, levels=c("OAW_MP","AMB_Control"))
##second term is the "control"
resALL <- results(dds, contrast=c("conditions.treat","OAW_MP","AMB_Control"))
#how many FDR < 10%?
table(resALL$padj<0.05)
# 0.1=138
# 0.05=100
# 0.01=55
summary(resALL)
# out of 21327 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 80, 0.38% 
# LFC < 0 (down)   : 58, 0.27% 
# outliers [1]     : 0, 0% 
# low counts [2]   : 12397, 58% 
# (mean count < 9)

nrow(resALL[resALL$padj<0.05 & !is.na(resALL$padj),])  # Num significantly differentially expressed genes excluding the no/low count genes   #228
#100
plotMA(resALL, main="OAW_Microplastics vs AMB_Control")

results <- as.data.frame(resALL)
head(results)

nrow(resALL[resALL$padj<0.1 & resALL$log2FoldChange > 0 & !is.na(resALL$padj),])
nrow(resALL[resALL$padj<0.1 & resALL$log2FoldChange < 0 & !is.na(resALL$padj),])
#UP in ALL 80
#DOWN in ALL 58

write.table(resALL, file="OAW_MP.txt", quote=F, sep="\t")

resALL=data.frame(resALL)
head(resALL)
go_OAW_MP = resALL %>%
  tibble::rownames_to_column(var = "iso") %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(iso, mutated_p_updown)
head(go_OAW_MP)
colnames(go_OAW_MP) <- c("gene", "pval")
head(go_OAW_MP)
write.csv(go_OAW_MP, file="go_OAW_MP.csv", quote=F, row.names=FALSE)

###############################################################################################
##############################################################################
#--------------get pvals
valOAW=cbind(resOAW$pvalue, resOAW$padj)
head(valOAW)
colnames(valOAW)=c("pval.OAW", "padj.OAW")
length(valOAW[,1])
table(complete.cases(valOAW))

valMP=cbind(resMP$pvalue, resMP$padj)
head(valMP)
colnames(valMP)=c("pval.MP", "padj.MP")
length(valMP[,1])
table(complete.cases(valMP))

valOAW_MP=cbind(resALL$pvalue, resALL$padj)
head(valOAW_MP)
colnames(valOAW_MP)=c("pval.OAW_MP", "padj.OAW_MP")
length(valOAW_MP[,1])
table(complete.cases(valOAW_MP))

######-------------make rlogdata and pvals table
rlog=rlogTransformation(dds, blind=TRUE) 
rld=assay(rlog)
head(rld)
colnames(rld)=paste(g$conditions.treat)
head(rld)
length(rld[,1])

rldpvals=cbind(rld,valOAW, valMP, valOAW_MP)
head(rldpvals)
dim(rldpvals)
# [1] 21569    25
table(complete.cases(rldpvals))
# FALSE  TRUE 
# 20077  1492 

write.csv(rldpvals, "AcerMicro_RLDandPVALS.csv", quote=F)

#################################################################################
# VENN Diagram to include both up and down regulated genes in common for PC02
library(VennDiagram)

pOAW_up=row.names(resOAW[resOAW$padj<0.1 & !is.na(resOAW$padj) & resOAW$log2FoldChange>0,])
length(pOAW_up) #18
pOAW_down=row.names(resOAW[resOAW$padj<0.1 & !is.na(resOAW$padj) & resOAW$log2FoldChange<0,])
length(pOAW_down) #18
pMP_up=row.names(resMP[resMP$padj<0.1 & !is.na(resMP$padj) & resMP$log2FoldChange>0,])
length(pMP_up) #48
pMP_down=row.names(resMP[resMP$padj<0.1 & !is.na(resMP$padj) & resMP$log2FoldChange<0,])
length(pMP_down) #9
pALL_up=row.names(resALL[resALL$padj<0.1 & !is.na(resALL$padj) & resALL$log2FoldChange>0,])
length(pALL_up) #80
pALL_down=row.names(resALL[resALL$padj<0.1 & !is.na(resALL$padj) & resALL$log2FoldChange<0,])
length(pALL_down) #58

pOAW=row.names(resOAW[resOAW$padj<0.1 & !is.na(resOAW$padj),])
length(pOAW)
pMP=row.names(resMP[resMP$padj<0.1 & !is.na(resMP$padj),])
pALL=row.names(resALL[resALL$padj<0.1 & !is.na(resALL$padj),])

#UP
pdegs05p_up=union(pOAW_up,pMP_up)
pdegs05_up=union(pdegs05p_up,pALL_up)
length(pdegs05_up)
#129

#DOWN
pdegs05p_down=union(pOAW_down,pMP_down)
pdegs05_down=union(pdegs05p_down,pALL_down)
length(pdegs05_down)
#80

#ALL
pdegs05p=union(pOAW,pMP)
pdegs05=union(pdegs05p,pALL)
length(pdegs05)
#209

###do separately for UP, DOWN, ALL and save as PDF
library(VennDiagram)
candidates=list("OAW"=pOAW_down, "MP"=pMP_down, "OAW+MP"=pALL_down)
prettyvenn=venn.diagram(
  x = candidates,
  filename=NULL,
  col = "transparent",
  fill = c("firebrick2", "darkorange","firebrick4"),
  alpha = 0.5,
  # label.col = c("darkred", "white", "darkgreen", "white", "white", "white", "blue4"),
  cex = 2.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col = c("firebrick3", "darkorange2","black"),
  cat.cex = 2.5,
  cat.fontfamily = "sans",
  cat.dist = c(0.08, 0.08,0.08),
  cat.pos = 1
);
grid.draw(prettyvenn)

###################################heatmaps for genes AMB vs MP
rldpvals <- read.csv(file="AcerMicro_RLDandPVALS.csv", row.names=1)
head(rldpvals)
rld_MP= rldpvals[,c(3:8,10,12,17,18,22,23)]
head(rld_MP)
gg=read.table("amil_gene2geneName.tab",sep="\t")
head(gg)

p.val=0.10 # FDR cutoff
conds=rld_MP[rld_MP$padj.MP<=p.val & !is.na(rld_MP$padj.MP),]
length(conds[,1])
#57
head(conds)

exp=conds[,1:10]
head(exp)
means=apply(exp,1,mean) # means of rows
explc=exp-means # subtracting them
head(explc)
head(gg)

library(tidyverse)
df_all_iso <- explc %>%
  rownames_to_column("V1") %>%
  left_join(gg) %>%
  mutate(V2 = gsub(" OS=.*", "", V2))
head(df_all_iso)
unanno=df_all_iso[,2:11]
head(unanno)

df_only_anno <- df_all_iso %>%
  filter(!is.na(V2))
rownames(df_only_anno) <- make.unique(df_only_anno$V2)
head(df_only_anno)
anno=df_only_anno[,2:11]
head(anno)

df_unanno <- df_all_iso %>%
  filter(is.na(V2))

##color schemes
ccol=colorRampPalette(rev(c("red","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
col0=colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
#dataframe of the samples
colnames(anno)

my_sample_col <- data.frame(treatment = c("MP","AMB","MP","AMB","MP","MP","AMB","MP","AMB", "AMB"))
row.names(my_sample_col)= colnames(anno)
my_colour = list(
  Treatment = c(MP = "darkorange", AMB = "dodgerblue"))

# big heat map of all annotated genes
pdf("MP_AMB_DEGS_0.10_annotated.pdf",height=5,width=25, onefile=F)
pheatmap(anno,cluster_cols=T,scale="row", color=col0, annotation_col= my_sample_col, annotation_colors =my_colour, show_rownames = T, show_colnames =F, border_color = "NA")
dev.off()
# big heat map of all DEGs, no annotation too
pdf("MP_AMB_DEGS_0.10_unannotated.pdf",height=5,width=4, onefile=F)
pheatmap(unanno,cluster_cols=T,scale="row",color=col0, annotation_col= my_sample_col, annotation_colors =my_colour, show_rownames = F,show_colnames = F, border_color = "NA")
dev.off()

###################################heatmaps for genes AMB vs OAW
head(rldpvals)
rld_OAW= rldpvals[,c(2,4,6,10,14,16,17,18,19,20,21)]
head(rld_OAW)
gg=read.table("amil_gene2geneName.tab",sep="\t")
head(gg)

p.val=0.10 # FDR cutoff
conds=rld_OAW[rld_OAW$padj.OAW<=p.val & !is.na(rld_OAW$padj.OAW),]
length(conds[,1])
#36
head(conds)

exp=conds[,1:9]
head(exp)
means=apply(exp,1,mean) # means of rows
explc=exp-means # subtracting them
head(explc)
head(gg)

library(tidyverse)
df_all_iso <- explc %>%
  rownames_to_column("V1") %>%
  left_join(gg) %>%
  mutate(V2 = gsub(" OS=.*", "", V2))
head(df_all_iso)
unanno=df_all_iso[,2:10]
head(unanno)

df_only_anno <- df_all_iso %>%
  filter(!is.na(V2))
rownames(df_only_anno) <- make.unique(df_only_anno$V2)
head(df_only_anno)
anno=df_only_anno[,2:10]
head(anno)

df_unanno <- df_all_iso %>%
  filter(is.na(V2))

##color schemes
ccol=colorRampPalette(rev(c("red","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
col0=colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
#dataframe of the samples
colnames(anno)

#cbPalette2 <- c("darkorange","firebrick2", "firebrick4")
my_sample_col <- data.frame(treatment = c("OAW","AMB","AMB","AMB","OAW","OAW","AMB","AMB","OAW"))
row.names(my_sample_col)= colnames(anno)
my_colour = list(
  Treatment = c(OAW = "firebrick2", AMB = "dodgerblue"))

# big heat map of all annotated genes
pdf("OAW_AMB_DEGS_0.10_annotated.pdf",height=5,width=10, onefile=F)
pheatmap(anno,cluster_cols=T,scale="row", color=col0, annotation_col= my_sample_col, annotation_colors =my_colour, show_rownames = T, show_colnames =F, border_color = "NA")
dev.off()
# big heat map of all DEGs, no annotation too
pdf("OAW_AMB_DEGS_0.10_unannotated.pdf",height=4,width=4, onefile=F)
pheatmap(unanno,cluster_cols=T,scale="row",color=col0, annotation_col= my_sample_col, annotation_colors =my_colour, show_rownames = F,show_colnames = F, border_color = "NA")
dev.off()

###################################heatmaps for genes AMB vs OAW+MP
head(rldpvals)
rld_OAW_MP= rldpvals[,c(1,4,6,9,10,11,13,15,17,18,24,25)]
head(rld_OAW_MP)
gg=read.table("amil_gene2geneName.tab",sep="\t")
head(gg)

p.val=0.10 # FDR cutoff
conds=rld_OAW_MP[rld_OAW_MP$padj.OAW_MP<=p.val & !is.na(rld_OAW_MP$padj.OAW_MP),]
length(conds[,1])
#138
head(conds)

exp=conds[,1:10]
head(exp)
means=apply(exp,1,mean) # means of rows
explc=exp-means # subtracting them
head(explc)
head(gg)

library(tidyverse)
df_all_iso <- explc %>%
  rownames_to_column("V1") %>%
  left_join(gg) %>%
  mutate(V2 = gsub(" OS=.*", "", V2))
head(df_all_iso)
unanno=df_all_iso[,2:11]
head(unanno)

df_only_anno <- df_all_iso %>%
  filter(!is.na(V2))
rownames(df_only_anno) <- make.unique(df_only_anno$V2)
head(df_only_anno)
anno=df_only_anno[,2:11]
head(anno)

df_unanno <- df_all_iso %>%
  filter(is.na(V2))

##color schemes
ccol=colorRampPalette(rev(c("red","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
col0=colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
#dataframe of the samples
colnames(anno)

#cbPalette2 <- c("darkorange","firebrick2", "firebrick4")
my_sample_col <- data.frame(treatment = c("OAW_MP","AMB","AMB","OAW_MP","AMB","OAW_MP","OAW_MP","OAW_MP","AMB","AMB"))
row.names(my_sample_col)= colnames(anno)
my_colour = list(
  Treatment = c(OAW_MP = "firebrick4", AMB = "dodgerblue"))

# big heat map of all annotated genes
pdf("OAW_MP_AMB_DEGS_0.10_annotated.pdf",height=10,width=20, onefile=F)
pheatmap(anno,cluster_cols=T,scale="row", color=col0, annotation_col= my_sample_col, annotation_colors =my_colour, show_rownames = T, show_colnames =F, border_color = "NA")
dev.off()
# big heat map of all DEGs, no annotation too
pdf("OAW_MP_AMB_DEGS_0.10_unannotated.pdf",height=4,width=4, onefile=F)
pheatmap(unanno,cluster_cols=T,scale="row",color=col0, annotation_col= my_sample_col, annotation_colors =my_colour, show_rownames = F,show_colnames = F, border_color = "NA")
dev.off()

####now go do GO enrichment analysis
#this analysis was conducted separatelyusing GO_MWU.R script

###now plotting the interesting genes belonging to GO categories
library(dplyr)
library(stringr)
iso2go=read.table("amil_gene2go.tab", row.names=1) %>%
  dplyr::rename("GO_ID" = "V2")
head(iso2go)

OAW_MP_res= read.table("OAW_MP.txt")
head(OAW_MP_res)

read.csv("AcerMicro_RLDandPVALS.csv")
rldpval=readr::read_csv("AcerMicro_RLDandPVALS.csv") %>%
  select(gene=1, everything())
head(rldpval)

col0=colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)

##pulling all of the immunity terms and making a data frame
immunity=iso2go %>%
  mutate(gene = rownames(iso2go)) %>%
  filter(grepl("GO:0051092|GO:0051091|GO:0070534|GO:0000209|GO:0002697|GO:0002637|GO:0043900|GO:0002831|GO:0031347|GO:0002753|GO:0039528|GO:0098586|GO:0045089|GO:0002833|GO:0031349|GO:0045088|GO:0002230|GO:0050691|GO:0002218|GO:0002221|GO:0002758|GO:0002253|GO:0002757|GO:0002764|GO:0071360|GO:0032608|GO:0001819|GO:0001817|GO:0032606|GO:0002756|GO:0034138|GO:0035666|GO:0002755|GO:0034130|GO:0034134|GO:0002224|GO:0008063|GO:0034142|GO:0032103|GO:0050778|GO:0038061|GO:0007250|GO:1901222|GO:0042088|GO:0006955|GO:0034341|GO:0071346|GO:0006952|GO:0098542|GO:0045087|GO:0034097|GO:0071345|GO:0030522", GO_ID)) %>%
  left_join(rldpval)
head(immunity)

immuno=immunity[,c(2,3,6,8,11,12,13,15,17,19,20,26,27)]
head(immuno)
row.names(immuno)=immuno$gene
immuno1=immuno[,2:13]
head(immuno1)

##gene table
gg=read.table("amil_gene2geneName.tab",sep="\t")
head(gg)

p.val=0.10 # raw pvalue for GO enriched
conds=immuno1[immuno1$pval.OAW_MP<=p.val & !is.na(immuno1$pval.OAW_MP),]
length(conds[,1])
#113
head(conds)

exp=conds[,1:10]
head(exp)
means=apply(exp,1,mean) # means of rows
explc=exp-means # subtracting them
head(explc)
head(gg)

library(tidyverse)
df_only_anno <- df_all_iso %>%
  filter(!is.na(V2))
rownames(df_only_anno) <- make.unique(df_only_anno$V2)
head(df_only_anno)
anno=df_only_anno[,2:11]
head(anno)

#dataframe of the samples
colnames(anno)

##write out and edit manually a few gene names that have weird descriptions
write.csv(anno, "immunity_info.csv", quote=F)
anno2=read.csv("immunity_info.csv", row.names=1)

#cbPalette2 <- c("darkorange","firebrick2", "firebrick4")
my_sample_col <- data.frame(treatment = c("OAW_MP","AMB","AMB","OAW_MP","AMB","OAW_MP","OAW_MP", "OAW_MP","AMB","AMB"))
row.names(my_sample_col)= colnames(anno)
my_colour = list(
  treatment = c(OAW_MP = "firebrick4", AMB = "dodgerblue"))

# big heat map of all annotated genes
library(pheatmap)
pdf("OAW_AMB_immunity.pdf",height=8.5,width=9, onefile=F)
pheatmap(anno2,cluster_cols=T,scale="row", color=col0, annotation_col= my_sample_col, annotation_colors =my_colour, show_rownames = T, show_colnames =F, border_color = "NA")
dev.off()
