library(clusterProfiler)
library(biomaRt)
library(EnsDb.Hsapiens.v79)
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
library(org.Hs.eg.db)
library(pathview)
library(AnnotationHub)
library(ensembldb)
library(GOSemSim)
library(dplyr)
library(enrichplot)
library(GSEABase)
library(GSVAdata)
library(Biobase)


##### Formatting the data #####
### TCGA Data
tcga_pos<-read.csv("(TCGA) Sig_Genes_TNBC_+ve (87).csv")
tcga_pos<-data.frame(tcga_pos[,-1])

tcga_neg<-read.csv("(TCGA) Sig_Genes_TNBC_-ve (84).csv")
tcga_neg<-data.frame(tcga_neg[,-1])


### RNAseq Sig Genes
VST_24<-read.csv("(2)Ctl_24_VST.csv")
colnames(VST_24)[2]<-'GENEID'
VST24<-VST_24[,c(2,4,8)]

VST_2W<-read.csv("(2)Ctl_2W_VST.csv")
colnames(VST_2W)[2]<-'GENEID'
VST2W<-VST_2W[,c(2,4,8)]


### RNAseq DEGs
W2_UP_DEG<-read.csv("2W_UP_DEG.csv")
deg2wup<-W2_UP_DEG[,-1]
deg2wup<-deg2wup[,1:3]

W2_DN_DEG<-read.csv("2W_DN_DEG.csv")
deg2wdn<-W2_DN_DEG[,-1]
deg2wdn<-deg2wdn[,1:3]

H24_UP_DEG<-read.csv("24_UP_DEG.csv")
deg24up<-H24_UP_DEG[,-1]
deg24up<-deg24up[,1:3]

H24_DN_DEG<-read.csv("24_DN_DEG.csv")
deg24dn<-H24_DN_DEG[,-1]
deg24dn<-deg24dn[,1:3]


#### Combining the data
##TCGA Genelists
#Negative PTEN correlation
tcga_ng<-data.frame(tcga_neg[,c(1:3)])
#Postitive PTEN correlation
tcga_ps<-data.frame(tcga_pos[,c(1:3)])


###Merging DEGs
#Upregulated - 24H
up24<-data.frame('GENEID'=deg24up[,1])
THGup24<-merge(up24, tcga_ps, by='GENEID')
THGu24<-merge(THGup24, VST24, by= 'GENEID')
write.csv(THGu24, 'TopHitGenes_up24.csv')
# 2 THG

#Downregulated - 24H
dn24<-data.frame('GENEID'=deg24dn[,1])
THGdn24<-merge(dn24, tcga_ng, by='GENEID')
THGd24<-merge(THGdn24, VST24, by= 'GENEID')
write.csv(THGd24, 'TopHitGenes_dn24.csv')
# 1 THG

#Upregulated - 2W
up2w<-data.frame('GENEID'=deg2wup[,1])
THGup2w<-merge(up2w, tcga_ps, by='GENEID')
THGu2w<-merge(THGup2w, VST2W, by= 'GENEID')
write.csv(THGu2w, 'TopHitGenes_up2w.csv')
# 4 THG


#Downregulated - 2W
dn2w<-data.frame('GENEID'=deg2wdn[,1])
THGdn2w<-merge(dn2w, tcga_ng, by='GENEID')
THGd2w<-merge(THGdn2w, VST2W, by= 'GENEID')
write.csv(THGd2w, 'TopHitGenes_dn2w.csv')
# 6 THG



#### Preparing the gene sets ####
### Alternating DEGs
Alternating_GOI <- read_csv("~/Youme/Year 4/Projects/BIOL 6013 - Advanced Research/Results/Bioinformatics/2 RNA-Seq Analysis/R Codes/(2) RNA-seq/Alternating_GOI.csv")
ADEG<-Alternating_GOI[,2:5]
LADEG <- list()                  
for(i in 1:ncol(ADEG)) {            
  LADEG[[i]] <- ADEG[ , i]
}
LADEG<-lapply(LADEG, function(x) x[!is.na(x)])
names(LADEG) <- colnames(ADEG)

UU<-merge(data.frame('Gene_ID'=LADEG$Up.Up),stat2w,by='Gene_ID')
UU<-UU[!duplicated(UU$Gene_ID),]
UU_ <- as.numeric(UU$log2FoldChange)
names(UU_) <- UU$Entrez_ID
UU_ <- sort(UU_, decreasing = TRUE)
head(UU_)
U_U<-UU$Entrez_ID

DD<-merge(data.frame('Gene_ID'=LADEG$Dn.Dn),stat2w,by='Gene_ID')
DD<-DD[!duplicated(DD$Gene_ID),]
DD_ <- as.numeric(DD$log2FoldChange)
names(DD_) <- DD$Entrez_ID
DD_ <- sort(DD_, decreasing = TRUE)
head(DD_)
D_D<-DD$Entrez_ID

DU<-merge(data.frame('Gene_ID'=LADEG$Dn.Up),stat2w,by='Gene_ID')
DU<-DU[!duplicated(DU$Gene_ID),]
DU_ <- as.numeric(DU$log2FoldChange)
names(DU_) <- DU$Entrez_ID
DU_ <- sort(DU_, decreasing = TRUE)
head(DU_)
D_U<-DU$Entrez_ID

UD<-merge(data.frame('Gene_ID'=LADEG$Up.Dn),stat2w,by='Gene_ID')
UD<-UD[!duplicated(UD$Gene_ID),]
UD_ <- as.numeric(UD$log2FoldChange)
names(UD_) <- UD$Entrez_ID
UD_ <- sort(UD_, decreasing = TRUE)
head(UD_)
U_D<-UD$Entrez_ID

### DEG Entrez_ID gene lists
DEG_24_U<-read.csv("24_UP_DEG.csv")
DEG24U<-DEG_24_U$ENTREZID

DEG_24_D<-read.csv("24_DN_DEG.csv")
DEG24D<-DEG_24_D$ENTREZID

DEG_2W_U<-read.csv("2W_UP_DEG.csv")
DEG2WU<-DEG_2W_U$ENTREZID

DEG_2W_D<-read.csv("2W_DN_DEG.csv")
DEG2WD<-DEG_2W_D$ENTREZID


#### 24 Hour ranked gene list
X_2_Ctl_24_VST<-read.csv("(2)Ctl_24_VST.csv")
stat24<-X_2_Ctl_24_VST[,c(2,4)]
colnames(stat24)[1]<-'ENSG'

## Identify all the GENEIDs found in the RNAseq
d2w<-data.frame(stat24[,1])
colnames(d2w)[1]<-'ENSG'
ENSG<-d2w[,1]

## Create a dataframe which identifies the GENEIDs and their associated Gene names from a gene database
entrezALL<- ensembldb::select(EnsDb.Hsapiens.v79, keys= ENSG, keytype="GENEID", columns=c("SYMBOL","ENTREZID","GENEID"))
head(entrezALL)
colnames(entrezALL)<-c("Gene_ID","Entrez_ID", "ENSG")

## Merge Gene names with RNAseq data
entrezALL2<-merge(entrezALL,stat24, by='ENSG')
head(entrezALL2)
stat24<-entrezALL2[!duplicated(entrezALL2$Gene_ID),]

FC24 <- as.numeric(stat24$log2FoldChange)
names(FC24) <- stat24$Entrez_ID
FC24 <- sort(FC24, decreasing = TRUE)
head(FC24)


#### 2 Weeks ranked gene list
X_2_Ctl_2W_VST<-read.csv("(2)Ctl_2W_VST.csv")
stat2w<-X_2_Ctl_2W_VST[,c(2,4)]
colnames(stat2w)[1]<-'ENSG'

## Identify all the GENEIDs found in the RNAseq
d2w<-data.frame(stat2w[,1])
colnames(d2w)[1]<-'ENSG'
ENSG<-d2w[,1]

## Create a dataframe which identifies the GENEIDs and their associated Gene names from a gene database
entrezALL<- ensembldb::select(EnsDb.Hsapiens.v79, keys= ENSG, keytype="GENEID", columns=c("SYMBOL","ENTREZID","GENEID"))
head(entrezALL)
colnames(entrezALL)<-c("Gene_ID","Entrez_ID", "ENSG")

## Merge Gene names with RNAseq data
entrezALL2<-merge(entrezALL,stat2w, by='ENSG')
head(entrezALL2)
stat2w<-entrezALL2[!duplicated(entrezALL2$Gene_ID),]

FC2W <- as.numeric(stat2w$log2FoldChange)
names(FC2W) <- stat2w$Entrez_ID
FC2W <- sort(FC2W, decreasing = TRUE)
head(FC2W)


#### Part 1 - KEGG pathway enrichment ####
#### KEGG GSEA - FC24/FC2W
gseaKEGG <- gseKEGG(geneList = FC2W, 
                    organism = "hsa", 
                    nPermSimple = 10000, 
                    eps = 0)

png('gseKEGG_2W.png', height=700, width=600)
gseaplot2(gseaKEGG, geneSetID = c(34,37,83,86,96,120))
dev.off()

gseaKEGG <- gseKEGG(geneList = FC24, 
                    organism = "hsa", 
                    nPermSimple = 10000, 
                    eps = 0)

png('gseKEGG_24.png', height=700, width=600)
gseaplot2(gseaKEGG, geneSetID = 18)
dev.off()


#### Part 2 - GO Analysis ####
## GO Over-representation analysis
# DEGs 24H + 2W
e2WU_ <- enrichGO(DEG2WU, OrgDb = org.Hs.eg.db, ont = 'ALL')
e2WU <- e2WU_@result
e2WU<-e2WU[e2WU$p.adjust < 0.05,]
e2WU<-e2WU[order(e2WU$GeneRatio, decreasing=TRUE),]
e2WU$GeneRatio<-sapply(e2WU$GeneRatio, function(x) eval(parse(text=x)))
e2WU<-data.frame(e2WU[], 'Reg'='2W')

e24U_ <- enrichGO(DEG24U, OrgDb = org.Hs.eg.db, ont = 'ALL')
e24U <- e24U_@result
e24U<-e24U[e24U$p.adjust < 0.05,]
e24U<-e24U[order(e24U$GeneRatio, decreasing=TRUE),]
e24U$GeneRatio<-sapply(e24U$GeneRatio, function(x) eval(parse(text=x)))
e24U<-data.frame(e24U[], 'Reg'='24H')

eU<-merge(e2WU[,1:3],e24U[,1:3],by='Description')
e24U<-merge(eU,e24U,by='Description')
e24U<-e24U[,-c(2:5)]
e24U<-e24U[order(e24U$GeneRatio, decreasing=TRUE),]
e2WU<-merge(eU,e2WU,by='Description')
e2WU<-e2WU[,-c(2:5)]
e2WU<-e2WU[order(e2WU$GeneRatio, decreasing=TRUE),]
eU<-rbind(e24U[1:30,],e2WU[1:30,])

UEU<-eU[c(21,23,27,29,30,37,39,45,52,55,58,59),]


e2WD_ <- enrichGO(DEG2WD, OrgDb = org.Hs.eg.db, ont = 'ALL')
e2WD <- e2WD_@result
e2WD<-e2WD[e2WD$p.adjust < 0.05,]
e2WD<-e2WD[order(e2WD$GeneRatio, decreasing=TRUE),]
e2WD$GeneRatio<-sapply(e2WD$GeneRatio, function(x) eval(parse(text=x)))
e2WD<-data.frame(e2WD[],'Reg'='2W')

e24D_ <- enrichGO(DEG24D, OrgDb = org.Hs.eg.db, ont = 'ALL')
e24D <- e24D_@result
e24D<-e24D[e24D$p.adjust < 0.05,]
e24D<-e24D[order(e24D$GeneRatio, decreasing=TRUE),]
e24D$GeneRatio<-sapply(e24D$GeneRatio, function(x) eval(parse(text=x)))
e24D<-data.frame(e24D[],'Reg'='24H')

eD<-merge(e2WD[,1:3],e24D[,1:3],by='Description')
e24D1<-merge(eD,e24D,by='Description')
e24D1<-e24D1[,-c(2:5)]
e24D1<-e24D1[order(e24D1$GeneRatio, decreasing=TRUE),]
e2WD1<-merge(eD,e2WD,by='Description')
e2WD1<-e2WD1[,-c(2:5)]
e2WD1<-e2WD1[order(e2WD1$GeneRatio, decreasing=TRUE),]
eD<-rbind(e2WD1,e24D1,e2WD[1:18,],e24D[1:18,])
DUD<-eD[c(1,3:9,13,14,29,30,35,40,42),]

png('GORA_UP.png', height=600, width=1100)
P<-ggplot(UEU,aes(GeneRatio,Description)) +
  geom_point(aes(size=Count,colour=`p.adjust`)) + 
  scale_colour_gradient(high='blue',low='red',n.breaks=10) +
  theme_bw() + 
  facet_grid(rows = vars(ONTOLOGY), scales = 'free', space = 'free_y',
             cols = vars(Reg)) +
  theme(axis.title=element_text(face="bold", size=17,colour = 'black'), 
        axis.text=element_text(face="bold", size=17,colour = 'black'),
        strip.text = element_text(face="bold",size=17))
P
dev.off()

#DEGs UU_ + DD_ + UD_ + DU_
UUU <- enrichGO(U_U, OrgDb = org.Hs.eg.db, ont = 'ALL')
uuu<-UUU@result
uuu<-uuu[uuu$p.adjust < 0.05,]
uuu<-uuu[order(uuu$GeneRatio, decreasing=TRUE),]
uuu$GeneRatio<-sapply(uuu$GeneRatio, function(x) eval(parse(text=x)))
uuu<-data.frame(uuu,'Reg'='UU')

DDD <- enrichGO(D_D, OrgDb = org.Hs.eg.db, ont = 'ALL')
ddd<-DDD@result
ddd<-ddd[ddd$p.adjust < 0.05,]
ddd<-ddd[order(ddd$GeneRatio, decreasing=TRUE),]
ddd$GeneRatio<-sapply(ddd$GeneRatio, function(x) eval(parse(text=x)))
ddd<-data.frame(ddd,'Reg'='DD')

uudd<-rbind(uuu[1:20,],ddd[1:20,])
uudd<-uudd[c(6,9,14,17,18,19,20,22,23,24,37,38),]

png('GORA_UU.png', height=600, width=1000)
P<-ggplot(uudd,aes(GeneRatio,Description)) +
  geom_point(aes(size=Count,colour=`p.adjust`)) + 
  scale_colour_gradient(high='blue',low='red',n.breaks=10) +
  theme_bw() + 
  facet_grid(rows = vars(ONTOLOGY),
             cols = vars(Reg), scales='free',space='free_y') +
  theme(axis.title=element_text(face="bold", size=17,colour = 'black'), 
        axis.text=element_text(face="bold", size=17,colour = 'black'),
        strip.text = element_text(face="bold",size=17))
P
dev.off()


#### GSVA ####
### Gathering matrix data for all samples in the RNAseq
## All VST
All<-read.csv("(2) NC(VST).csv")
colnames(All)[1]<-'ENSG'
All_named<-merge(statall, All, by = 'ENSG')
All_named<-All_named[!duplicated(All_named$Gene_ID),]
All<-as.matrix(data.frame(All_named[,c(5:13)], row.names=All_named$Gene_ID))
colnames(All)<-c(paste0(rep('Ctl_'),1:3), paste0(rep('24H_'),1:3), paste0(rep('2W_'),1:3))
All<-as.matrix(All)

All_Avg<-data.frame(row.names = rownames(All),
                    'Ctl'=rowMeans(All[,1:3]),
                    '24H'=rowMeans(All[,4:6]),
                    '2W'=rowMeans(All[,7:9]))

### Conducting the GSVA
## Hallmark Pathways
GSVA_hall<-read.gmt('~/Youme/Year 4/Projects/BIOL 6013 - Advanced Research/Results/Bioinformatics/3 Correlation Analysis/Gene Lists/h.all.v7.5.1.symbols.gmt')
GSVA_hall_list<-split(GSVA_hall[,2], GSVA_hall[,1])
GSVA_hall_list

# Conducting the GSVA & ordering results based on highest variance between groups
Hall<-data.frame(gsva(All,GSVA_hall_list))
#Hall<-data.frame(row.names=rownames(Hall),
#                 'Ctl'=rowMeans(Hall[,1:3]),
#                 '24H'=rowMeans(Hall[,4:6]),
#                 '2W'=rowMeans(Hall[,7:9]))

pv<-array(1,c(nrow(Hall),1))
cor<-array(1,c(nrow(Hall),1))
for (i in 1:nrow(Hall))
{
  VL<-cor.test(as.numeric(Hall[i,]),as.numeric(All['PTEN',]),method = "spearman")
  pv[i,1]<-VL$p.value              
  cor[i,1]<-VL$estimate}   
padj<-p.adjust(pv, method = 'fdr', n = length(pv)) ## fdr
Hall_<-data.frame('Row'=row.names(Hall),pv,padj,cor,Hall)
Hall_<-Hall_[abs(Hall_$cor) > 0.3,]
Hall_<-Hall_[Hall_$pv < 0.05,]
Hall_<-Hall_[order(Hall_$cor, decreasing=FALSE),]
All<-data.frame(All)
Hall_<-rbind(Hall_[,5:13],data.frame(row.names='PTEN', All['PTEN',]))

a<-as.matrix(Hall_)
png('Hall.png', height=400, width=1000)
heatmap(a, Rowv = NA, Colv = NA,
        col=bluered(500),scale = 'row')
dev.off()


##C6- Oncogenic Genesets
GSVA_c6<-read.gmt('~/Youme/Year 4/Projects/BIOL 6013 - Advanced Research/Results/Bioinformatics/3 Correlation Analysis/Gene Lists/c6.all.v7.5.1.symbols.gmt')
GSVA_c6_list<-split(GSVA_c6[,2], GSVA_c6[,1])
GSVA_c6_list

# Conducting the GSVA
C6<-data.frame(gsva(All,GSVA_c6_list))
#C6<-data.frame('Ctl'=rowMeans(C6[,1:3]),
#               '24H'=rowMeans(C6[,4:6]),
#               '2W'=rowMeans(C6[,7:9]))

pv<-array(1,c(nrow(C6),1))
cor<-array(1,c(nrow(C6),1))
for (i in 1:nrow(C6))
{
  VL<-cor.test(as.numeric(C6[i,]),as.numeric(All['PTEN',]),method = "spearman")
  pv[i,1]<-VL$p.value              
  cor[i,1]<-VL$estimate}   
padj<-p.adjust(pv, method = 'fdr', n = length(pv)) ## fdr
C6_<-data.frame('Row'=row.names(C6),pv,padj,cor,C6)
C6_<-C6_[abs(C6_$cor) > 0.3,]
C6_<-C6_[C6_$padj < 0.05,]
C6_<-C6_[order(C6_$cor, decreasing=FALSE),]
All<-data.frame(All)
C6_<-as.matrix(rbind(C6_[,5:13],data.frame(All['PTEN',])))

a<-as.matrix(C6_)
png('C6_.png', height=500, width=1500)
heatmap(a, Rowv = NA, Colv = NA,
        col=bluered(500),
        scale='row',
        margins = c(20,2))
dev.off()

