library(edgeR)
library(limma)
library(openxlsx)
library(clusterProfiler)
library(biomaRt)
library(EnsDb.Hsapiens.v79)
library(openxlsx)
library(gplots)
library(ggplot2)
library(ggpubr)
library(factoextra)
library(EnsDb.Hsapiens.v79)
library(DOSE)
library(EnhancedVolcano)
library(GSVA)
library(pheatmap)
library(RColorBrewer)


##Downloading and Merging 2 data
clinicaldata<-read.csv("(TCGA) Clinical Data.csv", header=TRUE)   

##Put the Sample.ID to row names
clinicaldata <- data.frame(clinicaldata[,-c(1:3)], row.names=clinicaldata[,3])  

##Isolate all TNBC cancer types
clinicaldata2<-split(clinicaldata,clinicaldata$Subtype)
clinicaldataB<-clinicaldata2$BRCA_Basal

write.csv(clinicaldataB, '(TCGA) Basal_ClinicalData.csv')

#### Merge cBioPortal with Xenabrowser ####
## Data collected from Xenabrowser 'TCGA Breast Cancer (BRCA)'
Gene<-read.csv("XenaData.csv", header=TRUE)
Gene1<-data.frame(Gene[,-1], row.names=Gene[,1])
Gene1<-t(Gene1)
Gene1<-data.frame('Sample.ID'=row.names(Gene1), Gene1)

##Gene data for Basal(TNBC) samples
GeneList<-row.names(clinicaldataB)
GeneList.<-data.frame('Sample.ID'=gsub('-','.',GeneList))
Gene_TNBC<-merge(Gene1,GeneList., by='Sample.ID')
Gene_TNBC<-data.frame(Gene_TNBC[,-1], row.names=Gene_TNBC[,1])
Gene_TNBC<-t(Gene_TNBC)

##Correlation Analysis for Gene_TNBC
pv<-array(1,c(20530,1))
cor<-array(1,c(20530,1))
for (i in 1:20530)
{
  VL<-cor.test(as.numeric(Gene_TNBC[i,]),as.numeric(Gene_TNBC['PTEN',]),method = "spearman")
  pv[i,1]<-VL$p.value              
  cor[i,1]<-VL$estimate}   
padj<-p.adjust(pv, method = 'fdr', n = length(pv)) ## fdr

Gene_TNBC_Stat<-data.frame('Row'=row.names(Gene_TNBC),pv,padj,cor,Gene_TNBC)
Gene_TNBC_Stat<-Gene_TNBC_Stat[order(Gene_TNBC_Stat$Row),]


## Identify all the GENEIDs found in the RNAseq
colnames(Gene_TNBC_Stat)[1]<-'SYMBOL'
SYMBOL<-Gene_TNBC_Stat[,1]

## Create a dataframe which identifies the GENEIDs and their associated Gene names from a gene database
entrezALL<- ensembldb::select(EnsDb.Hsapiens.v79, keys= SYMBOL, keytype="SYMBOL", columns=c("SYMBOL","ENTREZID","GENEID"))
head(entrezALL)

## Merge Gene names with RNAseq data
TNBC<-merge(entrezALL,Gene_TNBC_Stat, by='SYMBOL')
TNBC<-data.frame(TNBC[!duplicated(TNBC[,1]), ])
write.csv(TNBC, '(TCGA) All_Genes_TNBC.csv')


#### Isolate significant genes ####
TNBC_ord<-TNBC[order(TNBC$padj, decreasing=F),]
sig_TNBC<-data.frame(TNBC_ord[TNBC_ord$padj < 0.05,])
sig_TNBC<-sig_TNBC[complete.cases(sig_TNBC[,1]),]

sig_TNBCpos<-sig_TNBC[sig_TNBC$cor >= 0.3,]
sig_TNBCpos<-sig_TNBCpos[order(sig_TNBCpos$cor, decreasing=T),]
write.csv(sig_TNBCpos, '(TCGA) Sig_Genes_TNBC_+ve (87).csv')

sig_TNBCneg<-sig_TNBC[sig_TNBC$cor <= -0.3,]
sig_TNBCneg<-sig_TNBCneg[order(sig_TNBCneg$cor, decreasing=F),]
write.csv(sig_TNBCneg, '(TCGA) Sig_Genes_TNBC_-ve (84).csv')

ssig_TNBC<- rbind(sig_TNBCneg,sig_TNBCpos)
ssig_TNBC<-ssig_TNBC[order(ssig_TNBC$cor,decreasing = T),]
write.csv(ssig_TNBC, '(TCGA) SSig_Genes_TNBC (171).csv')


##### Heatmaps #####
### Heatmap for all basal genes
sigTNBC<-sig_TNBC[,c(1,7:177)]
sigTNBC<-data.frame(sigTNBC[,-1], row.names = sigTNBC[,1])
a1<-as.matrix(sigTNBC)

### REMOVE THE OUTLYING SAMPLES ###
Norm<-as.matrix(a1[!a1 %in% boxplot.stats(a1)$out])
a1[which(a1[]>max(Norm))]<-NA
a1[which(a1[]<min(Norm))]<-NA
a1<-a1[complete.cases(a1),]

tiff("TCGAbreastcorrelation_Norm(All).tiff", units="in", width=12, height=8, res=300)
distCor <- function(a1) as.dist(1-cor(t(a1)))
hclustAvg <- function(a1) hclust(a1, method="average")
heatmap.2(a1, trace='none', scale='row', margins=c(10,12), 
          cexRow=1, cexCol=.8, srtCol=45, adjCol=c(1,0), 
          Colv = T, Rowv = T, 
          dendrogram='none', distfun = distCor, hclustfun = hclustAvg, 
          col=colorRampPalette(rev(brewer.pal(10,"RdBu")))(250), 
          symbreaks = F, key = T, keysize = 1, main = 'TCGA', xlab = 'Samples', ylab = 'Genes of Interest')
dev.off()


### Heatmap for significant Basal genes
sigTNBC<-ssig_TNBC[,c(1,7:177)]
sigTNBC<-as.matrix(data.frame(sigTNBC[,-1], row.names = sigTNBC[,1]))
sigTNBC<-sigTNBC[,order(sigTNBC['PTEN',], decreasing=TRUE)]
a1<-as.matrix(sigTNBC)

### REMOVE THE OUTLYING SAMPLES ###
Norm<-as.matrix(a1[!a1 %in% boxplot.stats(a1)$out])
a1[which(a1[]>max(Norm))]<-NA
a1[which(a1[]<min(Norm))]<-NA
a1<-a1[complete.cases(a1),]

tiff("TCGAbreastCorrelation_Norm(Sig).tiff", units="in", width=12, height=8, res=300)
distCor <- function(a1) as.dist(1-cor(t(a1)))
hclustAvg <- function(a1) hclust(a1, method="average")
heatmap.2(a1, trace='none', scale='row', margins=c(10,12), 
          cexRow=1, cexCol=.8, srtCol=45, adjCol=c(1,0), 
          Colv = NA, Rowv = NA, 
          dendrogram='none', distfun = distCor, hclustfun = hclustAvg, 
          col=colorRampPalette(rev(brewer.pal(10,"RdBu")))(250), 
          symbreaks = F, key = T, keysize = 1, main = 'TCGA', xlab = 'Samples', ylab = 'Genes of Interest')
dev.off()


#### Pathway analysis ####
library(EnsDb.Hsapiens.v79)
library(clusterProfiler)
library(biomaRt)
library(openxlsx)
library(gplots)
library(ggplot2)
library(ggpubr)
library(factoextra)
library(DOSE)
library(EnhancedVolcano)
library(GSVA)
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)
library(cowplot)
library(DEGreport)
library(org.Hs.eg.db)
library(pathview)
library(tximport)
library(AnnotationHub)
library(ensembldb)
library(GOSemSim)
library(dplyr)
library(enrichplot)
library(GSEABase)
library(GSVAdata)

##Gathering the datasets
TNBC_A<-as.matrix(data.frame(row.names = TNBC$SYMBOL, TNBC[,7:177]))

TNBCN<-sig_TNBCneg$ENTREZID
##
TNBC_N<-as.matrix(data.frame(row.names=sig_TNBCneg$SYMBOL, sig_TNBCneg[7:177]))
TNBC_N<-TNBC_N[,order(TNBC_N[1,], decreasing=TRUE)]

TNBCP<-sig_TNBCpos$ENTREZID
##
TNBC_P<-as.matrix(data.frame(row.names=sig_TNBCpos$SYMBOL, sig_TNBCpos[7:177]))
TNBC_P<-TNBC_P[,order(TNBC_P[1,], decreasing=TRUE)]


#### Over-representation analysis ####
# GO sig pos genes
edo <- enrichGO(TNBCP, OrgDb = org.Hs.eg.db, ont = 'ALL')
edor <- edo@result
edor$GeneRatio<-sapply(edor$GeneRatio, function(x) eval(parse(text=x)))

png('GORA_UP.png', height=350, width=800)
P<-ggplot(edor,aes(GeneRatio,Description)) +
  geom_point(aes(size=Count,colour=`p.adjust`)) + 
  scale_colour_gradient(high='blue',low='red',n.breaks=10) +
  theme_bw() + 
  facet_grid(rows = vars(ONTOLOGY), scales = 'free', space = 'free_y') +
  theme(axis.title=element_text(face="bold", size=15,colour = 'black'), 
        axis.text=element_text(face="bold", size=15,colour = 'black'),
        strip.text = element_text(face="bold",size=15))
P
dev.off()

# GO sig neg genes
edo <- enrichGO(TNBCN, OrgDb = org.Hs.eg.db, ont = 'ALL')
edor <- edo@result
edor$GeneRatio<-sapply(edor$GeneRatio, function(x) eval(parse(text=x)))

png('GORA_DN.png', height=700, width=800)
P<-ggplot(edor,aes(GeneRatio,Description)) +
  geom_point(aes(size=Count,colour=`p.adjust`)) + 
  scale_colour_gradient(high='blue',low='red',n.breaks=10) +
  theme_bw() + 
  facet_grid(rows = vars(ONTOLOGY), scales = 'free', space = 'free_y') +
  theme(axis.title=element_text(face="bold", size=15,colour = 'black'), 
        axis.text=element_text(face="bold", size=15,colour = 'black'),
        strip.text = element_text(face="bold",size=15))
P
dev.off()



####GSEA####
gseTNBC <- as.numeric(TNBC$cor)
names(gseTNBC) <- TNBC$ENTREZID
gseTNBC <- sort(gseTNBC, decreasing = TRUE)
head(gseTNBC)


##KEGG
gseKEGG<-gseKEGG(gseTNBC, 
                 organism = "hsa", 
                 nPermSimple = 10000, 
                 eps = 0)

KEGG<-gseKEGG@result

png('gseaKEGG_TCGA.png', height=700, width=600)
gseaplot2(gseKEGG, geneSetID = c(9,20,35,89))
dev.off()


#### GSVA
## Hallmark Pathways
GSVA_hall<-read.gmt('~/Youme/Year 4/Projects/BIOL 6013 - Advanced Research/Results/Bioinformatics/3 Correlation Analysis/Gene Lists/h.all.v7.5.1.symbols.gmt')
GSVA_hall_list<-split(GSVA_hall[,2], GSVA_hall[,1])
GSVA_hall_list

# Conducting the GSVA & ordering results based on highest variance between groups
HALL<-data.frame(gsva(TNBC_A,GSVA_hall_list))
pv<-array(1,c(nrow(HALL),1))
cor<-array(1,c(nrow(HALL),1))
for (i in 1:nrow(HALL))
{
  VL<-cor.test(as.numeric(HALL[i,]),as.numeric(TNBC_A['PTEN',]),method = "spearman")
  pv[i,1]<-VL$p.value              
  cor[i,1]<-VL$estimate}   
padj<-p.adjust(pv, method = 'fdr', n = length(pv)) ## fdr
HALL_<-data.frame('Row'=row.names(HALL),pv,padj,cor,HALL)
HALL_<-HALL_[abs(HALL_$cor) > 0.3,]
HALL_<-HALL_[HALL_$padj < 0.05,]
HALL_<-HALL_[order(HALL_$cor, decreasing=FALSE),]
HALL_<-as.matrix(rbind(HALL_[,5:175],data.frame(row.names='PTEN', sig_TNBC[1,c(7:177)])))
HALL_<-HALL_[,order(HALL_[9,], decreasing=TRUE)]
a<-HALL_[c(2,5:9),]
png('HALL.png', height=400, width=2000)
heatmap(a, Rowv=NA, Colv = NA,
        col=redblue(100),
        margins=c(15,2))
dev.off()


##C6- Oncogenic Genesets
GSVA_c6<-read.gmt('~/Youme/Year 4/Projects/BIOL 6013 - Advanced Research/Results/Bioinformatics/3 Correlation Analysis/Gene Lists/c6.all.v7.5.1.symbols.gmt')
GSVA_c6_list<-split(GSVA_c6[,2], GSVA_c6[,1])
GSVA_c6_list

# GSVA & Visualisation
# Positively PTEN-correlated
C6<-data.frame(gsva(TNBC_A,GSVA_c6_list))
pv<-array(1,c(nrow(C6),1))
cor<-array(1,c(nrow(C6),1))
for (i in 1:nrow(C6))
{
  VL<-cor.test(as.numeric(C6[i,]),as.numeric(TNBC_A['PTEN',]),method = "spearman")
  pv[i,1]<-VL$p.value              
  cor[i,1]<-VL$estimate}   
padj<-p.adjust(pv, method = 'fdr', n = length(pv))
C6_<-data.frame('Row'=row.names(C6),pv,padj,cor,C6)
C6_<-C6_[abs(C6_$cor) > 0.3,]
C6_<-C6_[abs(C6_$padj) < 0.05,]
C6_<-C6_[order(C6_$cor, decreasing=FALSE),]
C6_<-as.matrix(rbind(C6_[,5:175],data.frame(row.names='PTEN', sig_TNBC[1,c(7:177)])))
C6_<-C6_[,order(C6_[7,], decreasing=TRUE)]
a<-C6_[]

C6_<-C6_[c('TGFB_UP.V1_DN','TGFB_UP.V1_UP','KRAS.KIDNEY_UP.V1_UP','STK33_SKM_DN','PTEN_DN.V1_DN','MTOR_UP.N4.V1_DN'),]

png('C6.png', height=300, width=800)
heatmap(a, Rowv=NA, Colv = NA,
        col=redblue(100),
        margins = c(10,5))
dev.off()
