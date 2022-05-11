if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

library(clusterProfiler)
library(biomaRt)
library(EnsDb.Hsapiens.v79)
library(openxlsx)
library(gplots)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(tidyverse)
library(factoextra)
library(EnsDb.Hsapiens.v79)
library(DOSE)
library(EnhancedVolcano)
library(GSVA)
library(pheatmap)
library(RColorBrewer)
library(DESeq2)
library(ComplexHeatmap)


#### Importing and merging the raw data ####
pten_ctl_2w<-read.csv("~/RNAseq/DEG/PTEN_0hr-vs-PTEN_2w/counts/raw_counts.csv")
head(pten_ctl_2w)

pten_ctl_24<-read.csv("~/RNAseq/DEG/PTEN_0hr-vs-PTEN_24hr/counts/raw_counts.csv")
head(pten_ctl_24)

## Format the raw data
colnames(pten_ctl_24)[1]<-'Gene'
colnames(pten_ctl_2w)[1]<-'Gene'

pten_ctl<-pten_ctl_24[,c(1,6:9)]
head(pten_ctl)
pten_24<-pten_ctl_24[,1:5]
head(pten_24)
pten_2w<-pten_ctl_2w[,1:5]
head(pten_2w)

## Merge the raw data
pten_ctl_24_2w<-merge(pten_2w,pten_ctl_24,by="Gene") 

## Reorder the raw data
pten_ctl_24_2w<-pten_ctl_24_2w[c(1,10:13,6:9,2:5)]

## Put Genenames as row names
RC<-data.frame(pten_ctl_24_2w[,-1], row.names = pten_ctl_24_2w[,1])

## Rename the data
colnames(RC)<-c(paste0(rep('Ctl_'), 1:4), paste0(rep('24H_'), 1:4), paste0(rep('2W_'), 1:4))

## Export the file
write.csv(RC, 'RNAseq_RawCounts.csv')


##### Standardising the data with DESeq2 #####
cts <- as.matrix(RC)
coldata<-matrix(c(rep('DOXneg',4),rep('DOX24h',4),rep('DOX2w',4)), ncol = 2, nrow = 12)
rownames(coldata)=colnames(RC)
colnames(coldata)<-c('condition', 'replicate')
coldata<-data.frame(coldata)
coldata[,'replicate']=c(1:4,1:4,1:4)
coldata$condition <- factor(coldata$condition)
coldata$replicate <- factor(coldata$replicate)


## Load data into correct format
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~0 + condition)

## Remove meaningless rows	***57500 > 26505***
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

## Conduct DEseq2
dds<-DESeq(dds)

## Tidy up the normalised gene count data
NC<-counts(dds, normalized=T)
write.csv(NC, 'NC.csv')

## Plot data of each gene individually
par(mfrow=c(2,3))
plotCounts(dds, gene="ENSG00000152583", intgroup="condition")
plotCounts(dds, gene="ENSG00000000003", intgroup="condition")
plotCounts(dds, gene="ENSG00000189221", intgroup="condition")
plotCounts(dds, gene="ENSG00000120129", intgroup="condition")
plotCounts(dds, gene="ENSG00000171862", intgroup="condition")# <-*PTEN*
plotCounts(dds, gene="ENSG00000148175", intgroup="condition")

## Significant expression analysis for each gene
PTEN<-data.frame(NC['ENSG00000171862',])
PTENm<-data.frame(mean(PTEN[1:4,]),
                  mean(PTEN[5:8,]), 
                  mean(PTEN[9:12,]))
PTENsd<-data.frame(sd(PTEN[1:4,]),
                   sd(PTEN[5:8,]), 
                   sd(PTEN[9:12,]))
PTEN<-data.frame(t(PTENsd),t(PTENm))
rownames(PTEN)=c('Ctl','24H','2W')
colnames(PTEN)=c('sd','mean')

write.csv(PTEN, 'PTEN_exp_DEseq2(1).csv')

pdf('PTEN_Expression.pdf', height=6, width=6)
ggplot(PTEN, aes(x=rownames(PTEN), y=mean)) +
  geom_point() +
  geom_errorbar(aes(ymin=PTEN[,2]-PTEN[,1], 
                    ymax=PTEN[,2]+PTEN[,1]), 
                width=.2)
dev.off()


## Further data standardisation
# Variance Stabilising Transformation (vst)

vsd <-vst(dds, blind=FALSE)
head(assay(vsd), 3)

write.csv(assay(vsd), 'NC(VST).csv')


####DEseq2 Statistical Analysis####
#CONTROL VS 24H
res_DOXneg_DOX24 <- results(dds,
                            contrast=c("condition", "DOXneg", "DOX24h"),
                            independentFiltering = F)
ctl_24_res <- res_DOXneg_DOX24[order(res_DOXneg_DOX24$padj),]

vsd_24ctl<-data.frame(assay(vsd[,c(1:8)]))

vsd_24<-merge(vsd_24ctl,ctl_24_res)
write.csv(vsd_24, 'Ctl_24_VST.csv')

vsd_24_sig<-vsd_24[vsd_24$padj<0.05,]
write.csv(vsd_24_sig, 'vsd_24_sig.csv')


#CONTROL VS 2W
res_DOXneg_DOX2w <- results(dds,
                            contrast=c("condition","DOXneg", "DOX2w"),
                            independentFiltering = F)
ctl_2w_res <- res_DOXneg_DOX2w[order(res_DOXneg_DOX2w$padj),]

vsd_2wctl<-data.frame(assay(vsd[,c(1:4, 9:12)]))

vsd_2w<-merge(vsd_2wctl,ctl_2w_res)
write.csv(vsd_2w, 'Ctl_2W_VST.csv')

vsd_2w_sig<-vsd_2w[vsd_2w$padj<0.05,]
write.csv(vsd_2w_sig, 'vsd_2w_sig.csv')


### Combining Significant Genes
sig2w<-data.frame(vsd_2w_sig[,'X'])
sig24<-data.frame(vsd_24_sig[,'X'])
colnames(sig24)='ENSG'
colnames(sig2w)='ENSG'
sigall<-rbind(sig24,sig2w)
sigall<-data.frame(sigall[!duplicated(sigall), ])
colnames(sigall)<-'ENSG'
VSD <- read.csv('NC(VST).csv')
colnames(VSD)[1]<-'ENSG'
sigvsd<-merge(VSD,sigall,by='ENSG')
sigvsd<-data.frame(sigvsd[,-1], row.names = sigvsd[,1])
colnames(sigvsd)<-c(paste0(rep('Ctl_'), 1:4), paste0(rep('24H_'), 1:4), paste0(rep('2W_'), 1:4))
write.csv(sigvsd, 'vsd_sig.csv')

#### Heatmap to display significant genes ####
## Run for 'sigvsd', vsd_2w[,9:16]', 'vsd_24[,9:16]'
a3<-as.matrix(data.frame(vsd2w[,10:17]))
colorbar<-bluered(100)
distCor <- function(a3) as.dist(1-cor(t(a3)))
hclustAvg <- function(a3) hclust(a3, method="average") #hclustfun=hclustAvg???distfun=distCor
condition_colors<-unlist(lapply(colnames(a3), function(x){
  if (grepl('Ctl',x)) {
    'blue'}
  else if (grepl('24H',x)){
    "green"}
  else{
    'red'
  }
}))
#label<-row.names(a3)
pdf('sigvsd2w_heatmap.pdf', height=8, width=12)
heatmap.2(a3, trace="none", density='none',scale="row",cexRow = 1,
          cexCol =1.5, zlim=c(-10,10),srtCol=0,adjCol=c(0.5,1),srtRow = 0,adjRow = c(0.5,1),
          Colv = F,Rowv = T,labRow = NA,
          distfun=distCor, hclustfun=hclustAvg, col=colorbar, key = T,keysize = 1,
          ColSideColors = condition_colors)
dev.off()


#### Cluster Heatmap ####
data_formatted<-as.matrix(dist(t(sigvsd), method = 'maximum'))
colors <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)

pdf("vsd_sig_clusterheatmap 2.pdf", width=10, height=10)

heatmap.2(data_formatted, trace=NULL, tracecol = NULL, col=colors)

dev.off()

####### PHeatmap with Clusters ########
select <- order(rowMeans(sigvsd),
                decreasing=TRUE)
data<-as.matrix(sigvsd[select,])

df <- data.frame(condition=colData(dds)[,c("condition")], row.names = rownames(coldata))

annot_colors=list(condition=c(DOXneg="blue",
                              DOX24h="green",
                              DOX2w='red'))

pdf("vsd_sig_clusterheatmap.pdf", width=12, height=8)
pheatmap(data, 
         cluster_rows=FALSE, 
         show_rownames=FALSE, 
         annotation_col=df,
         annotation_colors = annot_colors)
dev.off()


##### PCA Analysis of All Samples #####
pca<-prcomp(t(sigvsd), scale=FALSE)
head(pca)
plot(pca$x[,1],pca$x[,2])
pca.var<-pca$sdev^2
pca.var.per<-round(pca.var/sum(pca.var)*100,1)

# Scree Plot
pdf("Scree Plot.pdf", width=5, height=4)
barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation", ylim=c(0,100))
dev.off()

# PCA Plot
pdf("PCA.pdf",width=5, height=5)

pca.data<-data.frame(Sample=rownames(pca$x),
                     X=pca$x[,1],
                     Y=pca$x[,2])

ggplot(data=pca.data, aes(x=X, y=Y, label=Sample))+
  geom_point(size=1)+
  geom_text(color=c('blue','blue','blue','blue','green','green','green','green','red','red','red','red'))+
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep=""))+
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep=""))+
  theme_test()+
  theme_bw()+
  theme_linedraw()+
  ggtitle("PCA Graph for sig vst")

dev.off()
