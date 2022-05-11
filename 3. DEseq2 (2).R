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
##Import the raw data
RC<-read.csv("RNAseq_RawCounts.csv")
RC<-data.frame(RC[,-1], row.names = RC[,1])

## Rename the data
colnames(RC)<-c(paste0(rep('Ctl_'), 1:4), paste0(rep('24H_'), 1:4), paste0(rep('2W_'), 1:4))

## Remove samples: Ctl_4, 24H_4, 2W_4
RC<-RC[,c(1:3,5:7,9:11)]
write.csv(RC, '(2) RNAseq_RawCounts.csv')

##### Standardising the data with DESeq2 #####
cts <- as.matrix(RC)
coldata<-matrix(c(rep('DOXneg',3),rep('DOX24h',3),rep('DOX2w',3)), ncol = 2, nrow = 9)
rownames(coldata)=colnames(RC)
colnames(coldata)<-c('condition', 'replicate')
coldata<-data.frame(coldata)
coldata[,'replicate']=c(1:3,1:3,1:3)
coldata$condition <- factor(coldata$condition)
coldata$replicate <- factor(coldata$replicate)


## Load data into correct format
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~0 + condition)

## Remove meaningless rows	***57500 > 25070***
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

## Conduct DEseq2
dds<-DESeq(dds)

## Tidy up the normalised gene count data
NC<-counts(dds, normalized=T)
write.csv(NC, '(2) NC.csv')

## Plot data of each gene individually
par(mfrow=c(2,3))
plotCounts(dds, gene="ENSG00000152583", intgroup="condition")
plotCounts(dds, gene="ENSG00000000003", intgroup="condition")
plotCounts(dds, gene="ENSG00000189221", intgroup="condition")
plotCounts(dds, gene="ENSG00000120129", intgroup="condition")
plotCounts(dds, gene="ENSG00000171862", intgroup="condition")# <-*PTEN*
plotCounts(dds, gene="ENSG00000148175", intgroup="condition")
dev.off()

## Significant expression analysis for each gene
PTEN<-data.frame(NC['ENSG00000171862',])
PTENm<-data.frame(mean(PTEN[1:3,]),
                  mean(PTEN[4:6,]), 
                  mean(PTEN[7:9,]))
PTENsd<-data.frame(sd(PTEN[1:3,]),
                   sd(PTEN[4:6,]), 
                   sd(PTEN[7:9,]))
PTEN<-data.frame(t(PTENsd),t(PTENm))
rownames(PTEN)=c('Ctl','24H','2W')
colnames(PTEN)=c('sd','mean')

pdf('(2) PTEN_Expression.pdf', height=6, width=6)
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

write.csv(assay(vsd), '(2) NC(VST).csv')


####DEseq2 Statistical Analysis####
#CONTROL VS 24H
res_DOXneg_DOX24 <- results(dds,
                            contrast=c("condition", "DOX24h", "DOXneg"),
                            independentFiltering = F)
ctl_24_res<-data.frame(res_DOXneg_DOX24[order(res_DOXneg_DOX24$padj),])

vsd_24ctl<-data.frame(assay(vsd[,c(1:6)]))

vsd_24<-merge(ctl_24_res,vsd_24ctl,by=0)
write.csv(vsd_24, '(2)Ctl_24_VST.csv')

vsd_24_sig<-vsd_24[vsd_24$padj<0.05,]
write.csv(vsd_24_sig, '(2)vsd_24_sig.csv')


#CONTROL VS 2W
res_DOXneg_DOX2w <- results(dds,
                            contrast=c("condition","DOX2w", "DOXneg"),
                            independentFiltering = F)
ctl_2w_res<-data.frame(res_DOXneg_DOX2w[order(res_DOXneg_DOX2w$padj),])

vsd_2wctl<-data.frame(assay(vsd[,c(1:3, 7:9)]))

vsd_2w<-merge(ctl_2w_res,vsd_2wctl,by=0)
write.csv(vsd_2w, '(2)Ctl_2W_VST.csv')

vsd_2w_sig<-vsd_2w[vsd_2w$padj<0.05,]
write.csv(vsd_2w_sig, '(2)vsd_2w_sig.csv')


### Combining Significant Genes
sig2w<-data.frame(vsd_2w_sig[,'Row.names'])
sig24<-data.frame(vsd_24_sig[,'Row.names'])
colnames(sig24)='ENSG'
colnames(sig2w)='ENSG'
sigall<-rbind(sig24,sig2w)
sigall<-data.frame(sigall[!duplicated(sigall), ])
colnames(sigall)<-'ENSG'
VSD <- read.csv('(2) NC(VST).csv')
colnames(VSD)[1]<-'ENSG'
sigvsd<-merge(VSD,sigall,by='ENSG')
sigvsd<-data.frame(sigvsd[,-1], row.names = sigvsd[,1])
colnames(sigvsd)<-c(paste0(rep('Ctl_'), 1:3), paste0(rep('24H_'), 1:3), paste0(rep('2W_'), 1:3))
write.csv(sigvsd, '(2) vsd_sig.csv')

#### Heatmap to display significant genes ####
## Run for 'sigvsd', vsd_2w_sig[,9:14]', 'vsd_24_sig[,9:14]'#
a2<-data.frame(sigvsd)
colnames(a2)<-c(paste0(rep('Ctl_'),1:3), paste0(rep('24H_'),1:3), paste0(rep('2W_'),1:3))
a3<-as.matrix(a2)
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
pdf('(2)sigvsd24_heatmap.pdf', height=8, width=12)
heatmap.2(a3, trace="none", density='none',scale="row",cexRow = 1,
          cexCol =1.5, zlim=c(-10,10),srtCol=0,adjCol=c(0.5,1),srtRow = 0,adjRow = c(0.5,1),
          Colv = F,Rowv = T,labRow = NA,
          distfun=distCor, hclustfun=hclustAvg, col=colorbar, key = T,keysize = 1,
          ColSideColors = condition_colors)
dev.off()


#### Cluster Heatmap ####
data_formatted<-as.matrix(dist(t(sigvsd), method = 'maximum'))
colors <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)

pdf("(2) vsd_sig_clusterheatmap2.pdf", width=10, height=10)

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

pdf("(2) vsd_sig_clusterheatmap.pdf", width=12, height=8)
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
pdf("(2)Scree Plot.pdf", width=5, height=4)
barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation", ylim=c(0,100))
dev.off()

# PCA Plot
pdf("(2)PCA.pdf",width=5, height=5)

pca.data<-data.frame(Sample=rownames(pca$x),
                     X=pca$x[,1],
                     Y=pca$x[,2])

ggplot(data=pca.data, aes(x=X, y=Y, label=Sample))+
  geom_point(size=1)+
  geom_text(color=c('blue','blue','blue','green','green','green','red','red','red'))+
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep=""))+
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep=""))+
  theme_test()+
  theme_bw()+
  theme_linedraw()+
  ggtitle("PCA Graph for sig vst (2)")

dev.off()


##### MAKING A  VOLCANO PLOT #####
## I need to merge up and down genes for 2 weeks and 24 hours
## Then I need to merge these with the l2fc and padj for the genes

# 24 Hours
vc24<-data.frame(vsd_24[,c(1,3,7)])
vc24<-data.frame(vc24[,-1], row.names = vc24[,1])
colnames(vc24)<-c('log2FC','pvalue')

#**Add a new column called 'DEG' with values 'NO'
vc24$DEG<-"NO"			
# if log2Foldchange > 1 set as "UP"
vc24$DEG[vc24$log2FC > 1 & vc24$pvalue<0.05] <- "UP"
# if log2Foldchange < -1 set as "DOWN"
vc24$DEG[vc24$log2FC < -1 & vc24$pvalue<0.05] <- "DOWN"
head(vc24)

pdf('VolcanoPlot_24H.pdf', height=4, width=6)
p<- ggplot(data=vc24, aes(x=log2FC, y=-log10(pvalue))) +geom_point()
p2 <- p + geom_vline(xintercept=c(-1, 1), col="blue") +
  geom_hline(yintercept=-log10(0.05), col="blue")
p2
dev.off()

png('VolcanoPlot_24H(2).png', height=300, width=400)
p <- ggplot(data=vc24, aes(x=log2FC, y=-log10(pvalue), col=DEG)) + 
  geom_point()+theme_light()+
  geom_vline(xintercept=c(-1, 1), col="blue", linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), col="blue", linetype="dashed")+
  scale_color_manual(values=c("red","dark gray", "green"))
p+ font("xlab", size=14)+
  font("ylab", size=14)+
  font("xy.text", size=12)
dev.off()


# 2weeks
vc2w<-data.frame(vsd_2w[,c(1,3,7)])
vc2w<-data.frame(vc2w[,-1], row.names = vc2w[,1])
colnames(vc2w)<-c('log2FC','pvalue')

#**Add a new column called 'DEG' with values 'NO'
vc2w$DEG<-"NO"			
# if log2Foldchange > 1 set as "UP"
vc2w$DEG[vc2w$log2FC > 1 & vc2w$pvalue<0.05] <- "UP"
# if log2Foldchange < -1 set as "DOWN"
vc2w$DEG[vc2w$log2FC < -1 & vc2w$pvalue<0.05] <- "DOWN"
head(vc2w)

pdf('VolcanoPlot_2W.pdf', height=4, width=6)
p<- ggplot(data=vc2w, aes(x=log2FC, y=-log10(pvalue))) +geom_point()
p2 <- p + geom_vline(xintercept=c(-1, 1), col="blue") +
  geom_hline(yintercept=-log10(0.05), col="blue")
p2
dev.off()

png('VolcanoPlot_2W(2).png', height=300, width=400)
p <- ggplot(data=vc2w, aes(x=log2FC, y=-log10(pvalue), col=DEG)) + 
  geom_point()+theme_light()+
  geom_vline(xintercept=c(-1, 1), col="blue", linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), col="blue", linetype="dashed")+
  scale_color_manual(values=c("red", "dark gray", "green"))
p+ font("xlab", size=14)+
  font("ylab", size=14)+
  font("xy.text", size=12)
dev.off()



###### DEGs - printing the genes #####
## Other Gene IDs
EID<-as.character(rownames(vc24))

# Create a dataframe which identifies the GENEIDs and their associated Gene names from a gene database
Name<- ensembldb::select(EnsDb.Hsapiens.v79, keys= EID, keytype="GENEID", columns=c('SYMBOL','ENTREZID','GENEID'))
head(Name)


# 24 hour upreg DEGs 
upreg24<-data.frame(GENEID=rownames(vc24[vc24$log2FC > 1 & vc24$pvalue<0.05,]), vc24[vc24$log2FC > 1 & vc24$pvalue<0.05,])
upreg24<-merge(Name, upreg24,  by='GENEID')
upreg24<-upreg24[!duplicated(upreg24[,1]),]
write.csv(upreg24, file='24_UP_DEG.csv')

# 24 hour dnreg DEGs
downreg24<-data.frame(GENEID=rownames(vc24[vc24$log2FC < -1 & vc24$pvalue<0.05,]), vc24[vc24$log2FC < -1 & vc24$pvalue<0.05,])
downreg24<-merge(Name, downreg24, by='GENEID')
downreg24<-downreg24[!duplicated(downreg24[,1]),]
write.csv(downreg24, file='24_DN_DEG.csv')


# 2 weeks upreg DEGs 
upreg2w<-data.frame(GENEID=rownames(vc2w[vc2w$log2FC > 1 & vc2w$pvalue<0.05,]), vc2w[vc2w$log2FC > 1 & vc2w$pvalue<0.05,])
upreg2w<-merge(Name, upreg2w, by='GENEID')
upreg2w<-upreg2w[!duplicated(upreg2w[,1]),]
write.csv(upreg2w, file='2W_UP_DEG.csv')

# 2 weeks dnreg DEGs
downreg2w<-data.frame(GENEID=rownames(vc2w[vc2w$log2FC < -1 & vc2w$pvalue<0.05,]), vc2w[vc2w$log2FC < -1 & vc2w$pvalue<0.05,])
downreg2w<-merge(Name, downreg2w, by='GENEID')
downreg2w<-downreg2w[!duplicated(downreg2w[,1]),]
write.csv(downreg2w, file='2W_DN_DEG.csv')


### Identifying variable DEGs
DEG_UD<-merge(upreg24,downreg2w, by='SYMBOL')
DEG_DU<-merge(downreg24,upreg2w, by='SYMBOL')
DEG_DD<-merge(downreg24,downreg2w, by='SYMBOL')
DEG_UU<-merge(upreg24,upreg2w, by='SYMBOL')

Var_DEG<-list('Up/Dn'=DEG_UD$SYMBOL,
              'Dn/Up'=DEG_DU$SYMBOL,
              'Up/Up'=DEG_UU$SYMBOL,
              'Dn/Dn'=DEG_DD$SYMBOL)

Var_DEG<-lapply(Var_DEG, `length<-`, max(lengths(Var_DEG)))
Var_DEG<-as.data.frame(Var_DEG)

write.csv(Var_DEG, 'Alternating_GOI.csv')
