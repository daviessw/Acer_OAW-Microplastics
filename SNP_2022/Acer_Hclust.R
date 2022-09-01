setwd("~/Dropbox/BU/BU_Research/Katie_Microplastics/hclust")
bams=read.table("bams")[,1] # list of bam files, in the same order as you did your analysis
goods=c(1:length(bams))
#-------------
# clustering / PCoA based on identity by state (IBS) based on single read resampling
# (for low and/or uneven coverage)
#for 1
ma = as.matrix(read.table("acer1.ibsMat"))
ma=ma[goods,goods]
dimnames(ma)=list(bams[goods],bams[goods])
hc=hclust(as.dist(ma),"ave")
pdf("hclust1.pdf",height=5,width=6)
plot(hc,cex=0.5,xlab="")
dev.off()

#for 2
ma = as.matrix(read.table("acer2.ibsMat"))
ma=ma[goods,goods]
dimnames(ma)=list(bams[goods],bams[goods])
hc=hclust(as.dist(ma),"ave")
pdf("hclust2.pdf",height=5,width=6)
plot(hc,cex=0.5,xlab="")
dev.off()

#for 3
ma = as.matrix(read.table("acer3.ibsMat"))
ma=ma[goods,goods]
dimnames(ma)=list(bams[goods],bams[goods])
hc=hclust(as.dist(ma),"ave")
pdf("hclust3.pdf",height=5,width=6)
plot(hc,cex=0.5,xlab="")
dev.off()

#for 4
ma = as.matrix(read.table("acer4.ibsMat"))
ma=ma[goods,goods]
dimnames(ma)=list(bams[goods],bams[goods])
hc=hclust(as.dist(ma),"ave")
pdf("hclust4.pdf",height=5,width=6)
plot(hc,cex=0.5,xlab="")
dev.off()

#for 5
ma = as.matrix(read.table("acer5.ibsMat"))
ma=ma[goods,goods]
dimnames(ma)=list(bams[goods],bams[goods])
hc=hclust(as.dist(ma),"ave")
pdf("hclust5.pdf",height=5,width=6)
plot(hc,cex=0.5,xlab="")
dev.off()