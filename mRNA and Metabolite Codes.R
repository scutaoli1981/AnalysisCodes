################################################
##Codes for differential expression analysis of mRNA data
setwd("...")#set the mRNA data file location
library(openxlsx)
library(limma)
library(qvalue)
library(DESeq2)
library(ggplot2)
library(patchwork)
htseq_counts<-read.csv("mRNA data",sep = "\t",row.names = 1)
zeroindex<-apply(htseq_counts,1,function(x){sum(x==0)})
htseq_counts<-htseq_counts[zeroindex<30,]
head(htseq_counts)
dim(htseq_counts)
col_data<-data.frame(Samples=colnames(htseq_counts),condition=rep(c("HA","PS","S"),c(20,20,20)))
col_data$condition <- factor(col_data$condition)
col_data
##Parameters
volcanofcbig<-log2(1.5)
volcanofcsmall<-log2(1/1.5)
volcanopval<-0.1
##PS vs. HA
htseq_counts1<-htseq_counts[,c(21:40,1:20)]
col_data1<-col_data[c(21:40,1:20),]
col_data1
#write.csv(col_data1,file = "col_data1.csv",row.names = F)
dds <- DESeqDataSetFromMatrix(countData = htseq_counts1,
                              colData = col_data1,
                              design = ~ condition)
# Normalize the count data
dds <- estimateSizeFactors(dds)
# To get the normalized counts
normalized_counts <- counts(dds, normalized=TRUE)
# Differential expression analysis
dds <- DESeq(dds)
# Plot a volcano plot
res1 <- results(dds, contrast=c("condition","PS","HA"))
#res1$padj<-p.adjust(res1$pvalue,method = "BH")
res1$padj<-qvalue(res1$pvalue)$qvalues
volcanorawdata1<-res1
volcanorawdata1$p.adjust<-volcanorawdata1$padj
volcanorawdata1$Threshold<-c("NoChange")
volcanorawdata1$Threshold[volcanorawdata1$log2FoldChange>=volcanofcbig & volcanorawdata1$p.adjust <= volcanopval]<-c("UP")
volcanorawdata1$Threshold[volcanorawdata1$log2FoldChange<=volcanofcsmall & volcanorawdata1$p.adjust <= volcanopval]<-c("DOWN")
head(volcanorawdata1)
table(volcanorawdata1$Threshold)
pp1<-ggplot(data=volcanorawdata1,aes(x=log2FoldChange,y=-log10(p.adjust)))+theme_bw()+
  theme(panel.grid.minor = element_blank() ,panel.grid.major = element_blank())+
  geom_hline(yintercept =-log10(volcanopval),col="grey",lty=2,lwd=1)+
  geom_vline(xintercept =volcanofcsmall,col="grey",lty=2,lwd=1)+
  geom_vline(xintercept =volcanofcbig,col="grey",lty=2,lwd=1)+
  geom_point(aes(colour=Threshold),alpha=0.9,size=2)+theme(legend.position ="none")+
  scale_color_manual(values = c("#0072B9","grey","#D90605"))+
  xlab("Fold Change (log2)") + ylab("Adjusted P Value (-log10)") +ggtitle("PS vs. HA")
pp1
##S vs. HA
htseq_counts2<-htseq_counts[,c(41:60,1:20)]
col_data2<-col_data[c(41:60,1:20),]
col_data2
#write.csv(col_data1,file = "col_data1.csv",row.names = F)
dds <- DESeqDataSetFromMatrix(countData = htseq_counts2,
                              colData = col_data2,
                              design = ~ condition)
# Normalize the count data
dds <- estimateSizeFactors(dds)
# To get the normalized counts
normalized_counts <- counts(dds, normalized=TRUE)
# Differential expression analysis
dds <- DESeq(dds)
# Plot a volcano plot
res2 <- results(dds, contrast=c("condition","S","HA"))
#res2$padj<-p.adjust(res2$pvalue,method = "BH")
res2$padj<-qvalue(res2$pvalue)$qvalues
volcanorawdata2<-res2
volcanorawdata2$p.adjust<-volcanorawdata2$padj
volcanorawdata2$Threshold<-c("NoChange")
volcanorawdata2$Threshold[volcanorawdata2$log2FoldChange>=volcanofcbig & volcanorawdata2$p.adjust <= volcanopval]<-c("UP")
volcanorawdata2$Threshold[volcanorawdata2$log2FoldChange<=volcanofcsmall & volcanorawdata2$p.adjust <= volcanopval]<-c("DOWN")
head(volcanorawdata2)
table(volcanorawdata2$Threshold)
pp2<-ggplot(data=volcanorawdata2,aes(x=log2FoldChange,y=-log10(p.adjust)))+theme_bw()+
  theme(panel.grid.minor = element_blank() ,panel.grid.major = element_blank())+
  geom_hline(yintercept =-log10(volcanopval),col="grey",lty=2,lwd=1)+
  geom_vline(xintercept =volcanofcsmall,col="grey",lty=2,lwd=1)+
  geom_vline(xintercept =volcanofcbig,col="grey",lty=2,lwd=1)+
  geom_point(aes(colour=Threshold),alpha=0.9,size=2)+theme(legend.position ="none")+
  scale_color_manual(values = c("#0072B9","grey","#D90605"))+
  xlab("Fold Change (log2)") + ylab("Adjusted P Value (-log10)") +ggtitle("S vs. HA")
pp2
##
##S vs. PS
htseq_counts3<-htseq_counts[,c(41:60,21:40)]
col_data3<-col_data[c(41:60,21:40),]
col_data3
#write.csv(col_data1,file = "col_data1.csv",row.names = F)
dds <- DESeqDataSetFromMatrix(countData = htseq_counts3,
                              colData = col_data3,
                              design = ~ condition)
# Normalize the count data
dds <- estimateSizeFactors(dds)
# To get the normalized counts
normalized_counts <- counts(dds, normalized=TRUE)
# Differential expression analysis
dds <- DESeq(dds)
# Plot a volcano plot
res3 <- results(dds, contrast=c("condition","S","PS"))
#res3$padj<-p.adjust(res3$pvalue,method = "BH")
res3$padj<-qvalue(res3$pvalue)$qvalues
volcanorawdata3<-res3
volcanorawdata3$p.adjust<-volcanorawdata3$padj
volcanorawdata3$Threshold<-c("NoChange")
volcanorawdata3$Threshold[volcanorawdata3$log2FoldChange>=volcanofcbig & volcanorawdata3$p.adjust <= volcanopval]<-c("UP")
volcanorawdata3$Threshold[volcanorawdata3$log2FoldChange<=volcanofcsmall & volcanorawdata3$p.adjust <= volcanopval]<-c("DOWN")
head(volcanorawdata3)
table(volcanorawdata3$Threshold)
pp3<-ggplot(data=volcanorawdata3,aes(x=log2FoldChange,y=-log10(p.adjust)))+theme_bw()+
  theme(panel.grid.minor = element_blank() ,panel.grid.major = element_blank())+
  geom_hline(yintercept =-log10(volcanopval),col="grey",lty=2,lwd=1)+
  geom_vline(xintercept =volcanofcsmall,col="grey",lty=2,lwd=1)+
  geom_vline(xintercept =volcanofcbig,col="grey",lty=2,lwd=1)+
  geom_point(aes(colour=Threshold),alpha=0.9,size=2)+theme(legend.position ="none")+
  scale_color_manual(values = c("#0072B9","grey","#D90605"))+
  xlab("Fold Change (log2)") + ylab("Adjusted P Value (-log10)") +ggtitle("S vs. PS")
pp3
##
pp1+pp2+pp3
ggsave("RNA_Volcanoplot.qvalue.pdf",width =12,height = 8)
##GO enrichment
chayigeneall.go<-enrichGO(gene = unique(chayigeneall1$Entrez_geneID),minGSSize=1,
                          OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", ont = "BP",
                          pvalueCutoff = 1,qvalueCutoff = 1,readable = FALSE,pool = TRUE)
gores_ALL_df<-chayigeneall.go@result
head(gores_ALL_df)
dim(gores_ALL_df)
pp5<-enrichplot::dotplot(chayigeneall.go,showCategory=20)+
  theme(panel.grid.minor = element_blank() ,panel.grid.major = element_blank())
pp5
(((pp1 | pp2) / (pp3|pp4)) + plot_layout(heights =c(1,2))) | pp5+
  plot_layout(widths=c(2,1))
ggsave("RNA_Allplots.go.qvalue.pdf",width =20,height = 14)
##Heatmap for differentially expressed genes
library(pheatmap)
htseq_countsx1<-htseq_counts
htseq_countsx<-htseq_countsx1[rownames(htseq_countsx1)%in%chayigeneall,]
dim(htseq_countsx)
annotation_col_hca = data.frame(
  SampleType = factor(rep(c("HA","PS","S"),c(20,20,20)))
)
rownames(annotation_col_hca) = colnames(htseq_countsx)
ann_colors = list(
  SampleType=c(HA="#46BCE0",PS="#F49A7D",S="#DB0200")
)
htseq_countsxx<-t(scale(t(htseq_countsx)))
htseq_countsxx[abs(htseq_countsxx)>2]<-2
pheatmap(htseq_countsxx, scale = "none", cluster_rows = TRUE,cluster_cols = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100),show_rownames = F,
         annotation_col = annotation_col_hca,border_color=NA,cutree_rows=2,
         annotation_colors = ann_colors)
pdf("chayigenesheatmap.qvalue.pdf",width = 13,height = 15)
pheatmap(htseq_countsxx, scale = "none", cluster_rows = TRUE,cluster_cols = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100),show_rownames = F,
         annotation_col = annotation_col_hca,border_color=NA,cutree_rows=2,
         annotation_colors = ann_colors)
dev.off()
################################################################
###Codes for differential expression analysis of metabolomics data
setwd("...")
library(openxlsx)
library(limma)
library(qvalue)
library(DESeq2)
library(ggplot2)
library(patchwork)
library(samr)
metabodata<-read.xlsx("Metabolite data",sheet = 1)
head(metabodata)
dim(metabodata)
metabodata1<-metabodata[,-c(1:5)]
rownames(metabodata1)<-metabodata[,3]
head(metabodata1)
metabodata2<-log2(sweep(metabodata1,2,apply(metabodata1, 2, median,na.rm=T),FUN = "/"))
head(metabodata2)
col_data<-data.frame(Samples=colnames(metabodata2),condition=rep(c("HA","PS","S"),c(20,20,20)))
col_data
dfANOVA<-metabodata2
alldatalist<-list(x=as.matrix(dfANOVA),y=rep(c(1:3),times=c(20,20,20)), 
                  geneid=as.character(1:nrow(dfANOVA)),
                  genenames=rownames(dfANOVA), logged2=TRUE)
samr.obj<-samr(alldatalist, resp.type="Multiclass", nperms=1000,random.seed=1234567)
pv<-samr.pvalues.from.perms(samr.obj$tt, samr.obj$ttstar)
pv[pv>1]<-1
dfANOVA$adjPval<-pv
head(dfANOVA)
dim(dfANOVA)
dfANOVA1<-data.frame(mz=metabodata$`m/z`,adjPval=dfANOVA$adjPval)
dfANOVA1$adjPval.gr<-ifelse(dfANOVA$adjPval<0.1,"Sign","NotSign")
head(dfANOVA1)
dfANOVApp<-ggplot(data=dfANOVA1,aes(x=mz,y=-log10(adjPval)))+theme_bw()+
  theme(panel.grid.minor = element_blank() ,panel.grid.major = element_blank())+
  geom_point(aes(colour=adjPval.gr),alpha=0.9,size=2)+theme(legend.position ="none")+
  scale_color_manual(values = c("green","#F86E6E"))+
  xlab("Peaks (m/z)") + ylab("P Value (-log10)")
dfANOVApp
##Parameters
volcanofcbig<-log2(1.5)
volcanofcsmall<-log2(1/1.5)
volcanopval<-0.1
##PS vs. HA
metabores1<-metabodata2[,c(21:40,1:20)]
col_data1<-col_data[c(21:40,1:20),]
col_data1
alldatalist1<-list(x=as.matrix(metabores1),y=rep(c(1,2),times=c(20,20)), 
                  geneid=as.character(1:nrow(metabores1)),
                  genenames=rownames(metabores1), logged2=TRUE)
samr.obj1<-samr(alldatalist1, resp.type="Two class unpaired", nperms=1000,random.seed=1234567)
pv1<-samr.pvalues.from.perms(samr.obj1$tt, samr.obj1$ttstar)
pv1[pv1>1]<-1
metabores1$padj<-pv1
metabores1$log2FoldChange<-rowMeans(metabores1[,21:40])-rowMeans(metabores1[,1:20])
head(metabores1)
volcanorawdata1<-metabores1
volcanorawdata1$p.adjust<-volcanorawdata1$padj
volcanorawdata1$Threshold<-c("NoChange")
volcanorawdata1$Threshold[volcanorawdata1$log2FoldChange>=volcanofcbig & volcanorawdata1$p.adjust <= volcanopval]<-c("UP")
volcanorawdata1$Threshold[volcanorawdata1$log2FoldChange<=volcanofcsmall & volcanorawdata1$p.adjust <= volcanopval]<-c("DOWN")
head(volcanorawdata1)
table(volcanorawdata1$Threshold)
pp1<-ggplot(data=volcanorawdata1,aes(x=log2FoldChange,y=-log10(p.adjust)))+theme_bw()+
  theme(panel.grid.minor = element_blank() ,panel.grid.major = element_blank())+
  geom_hline(yintercept =-log10(volcanopval),col="grey",lty=2,lwd=1)+
  geom_vline(xintercept =volcanofcsmall,col="grey",lty=2,lwd=1)+
  geom_vline(xintercept =volcanofcbig,col="grey",lty=2,lwd=1)+
  geom_point(aes(colour=Threshold),alpha=0.9,size=2)+theme(legend.position ="right")+
  scale_color_manual(values = c("DOWN"="#0072B9","NoChange"="grey","UP"="#D90605"))+
  xlab("Fold Change (log2)") + ylab("Adjusted P Value (-log10)") +ggtitle("PS vs. HA")
pp1
##S vs. HA
metabores2<-metabodata2[,c(41:60,1:20)]
col_data2<-col_data[c(41:60,1:20),]
col_data2
alldatalist2<-list(x=as.matrix(metabores2),y=rep(c(1,2),times=c(20,20)), 
                   geneid=as.character(1:nrow(metabores2)),
                   genenames=rownames(metabores2), logged2=TRUE)
samr.obj2<-samr(alldatalist2, resp.type="Two class unpaired", nperms=1000,random.seed=1234567)
pv2<-samr.pvalues.from.perms(samr.obj2$tt, samr.obj2$ttstar)
pv2[pv2>1]<-1
metabores2$padj<-pv2
metabores2$log2FoldChange<-rowMeans(metabores2[,21:40])-rowMeans(metabores2[,1:20])
head(metabores2)
volcanorawdata2<-metabores2
volcanorawdata2$p.adjust<-volcanorawdata2$padj
volcanorawdata2$Threshold<-c("NoChange")
volcanorawdata2$Threshold[volcanorawdata2$log2FoldChange>=volcanofcbig & volcanorawdata2$p.adjust <= volcanopval]<-c("UP")
volcanorawdata2$Threshold[volcanorawdata2$log2FoldChange<=volcanofcsmall & volcanorawdata2$p.adjust <= volcanopval]<-c("DOWN")
head(volcanorawdata2)
table(volcanorawdata2$Threshold)
pp2<-ggplot(data=volcanorawdata2,aes(x=log2FoldChange,y=-log10(p.adjust)))+theme_bw()+
  theme(panel.grid.minor = element_blank() ,panel.grid.major = element_blank())+
  geom_hline(yintercept =-log10(volcanopval),col="grey",lty=2,lwd=1)+
  geom_vline(xintercept =volcanofcsmall,col="grey",lty=2,lwd=1)+
  geom_vline(xintercept =volcanofcbig,col="grey",lty=2,lwd=1)+
  geom_point(aes(colour=Threshold),alpha=0.9,size=2)+theme(legend.position ="right")+
  scale_color_manual(values = c("DOWN"="#0072B9","NoChange"="grey","UP"="#D90605"))+
  xlab("Fold Change (log2)") + ylab("Adjusted P Value (-log10)") +ggtitle("S vs. HA")
pp2
##S vs. PS
metabores3<-metabodata2[,c(41:60,21:40)]
col_data3<-col_data[c(41:60,21:40),]
col_data3
alldatalist3<-list(x=as.matrix(metabores3),y=rep(c(1,2),times=c(20,20)), 
                   geneid=as.character(1:nrow(metabores3)),
                   genenames=rownames(metabores3), logged2=TRUE)
samr.obj3<-samr(alldatalist3, resp.type="Two class unpaired", nperms=1000,random.seed=1234567)
pv3<-samr.pvalues.from.perms(samr.obj3$tt, samr.obj3$ttstar)
pv3[pv3>1]<-1
metabores3$padj<-pv3
metabores3$log2FoldChange<-rowMeans(metabores3[,21:40])-rowMeans(metabores3[,1:20])
head(metabores3)
volcanorawdata3<-metabores3
volcanorawdata3$p.adjust<-volcanorawdata3$padj#Pval
volcanorawdata3$Threshold<-c("NoChange")
volcanorawdata3$Threshold[volcanorawdata3$log2FoldChange>=volcanofcbig & volcanorawdata3$p.adjust <= volcanopval]<-c("UP")
volcanorawdata3$Threshold[volcanorawdata3$log2FoldChange<=volcanofcsmall & volcanorawdata3$p.adjust <= volcanopval]<-c("DOWN")
head(volcanorawdata3)
table(volcanorawdata3$Threshold)
pp3<-ggplot(data=volcanorawdata3,aes(x=log2FoldChange,y=-log10(p.adjust)))+theme_bw()+
  theme(panel.grid.minor = element_blank() ,panel.grid.major = element_blank())+
  geom_hline(yintercept =-log10(volcanopval),col="grey",lty=2,lwd=1)+
  geom_vline(xintercept =volcanofcsmall,col="grey",lty=2,lwd=1)+
  geom_vline(xintercept =volcanofcbig,col="grey",lty=2,lwd=1)+
  geom_point(aes(colour=Threshold),alpha=0.9,size=2)+theme(legend.position ="right")+
  scale_color_manual(values = c("DOWN"="#0072B9","NoChange"="grey","UP"="#D90605"))+
  xlab("Fold Change (log2)") + ylab("Adjusted P Value (-log10)") +ggtitle("S vs. PS")
pp3
##
dfANOVApp/(pp1+pp2+pp3)+
  plot_layout(heights =c(1.5,1))
ggsave("metabo_allplot.adjPval.sam.pdf",width =17 ,height = 10)


