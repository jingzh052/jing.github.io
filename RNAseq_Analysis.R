BiocManager::install(version = "3.12")

BiocManager::install('EnhancedVolcano')
BiocManager::install('apeglm')
BiocManager::install("ashr")
BiocManager::install("ComplexHeatmap")
BiocManager::install("DESeq2")
BiocManager::install("edgeR")
BiocManager::install("goseq")
BiocManager::install("pathview")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("msigdbr")
BiocManager::install("VennDetail")
library(EnhancedVolcano)
library(apeglm)
library(ashr)
library(ComplexHeatmap)
library(DESeq2)
library(edgeR)
library(pathview)
library(clusterProfiler)
library(msigdbr)
library(magrittr)
library(pheatmap)
library(gplot)
library(ggfortify)

library(org.Mm.eg.db)
library(AnnotationDbi)

library(tidyverse)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(goseq)
library(VennDetail)

## Count Data Download and Summary

counts=read.delim2("readcounts_exon_unstrand_Rmatrix.txt",header=TRUE,sep='\t')

colnames(counts) <- c('Geneid','SRR445','SRR446','SRR447','SRR448','SRR449','SRR450',
                      'SRR451','SRR452','SRR453','SRR454','SRR455','SRR456')

counts <- counts[order(counts$Geneid),]

data <- counts[,-1]

rownames(data) <- counts[,1]

counts <- data

summary(counts)

###########################

##Pre-Processing data

###########################

sampleName_full <- c('SRR445','SRR446','SRR447','SRR448','SRR449','SRR450',
                'SRR451','SRR452','SRR453','SRR454','SRR455','SRR456')

treatmentName_full <- c("AA", "AA", "AA", "AA", "Ctr","Ctr","Ctr","Ctr","HFD","HFD","HFD","HFD")

colnames(counts) <- treatmentName_full

medianCountByGroup <- t(apply(counts, 1, tapply,treatmentName_full, median))

maxMedian<-apply(medianCountByGroup, 1, max)

counts_filtered<-counts[maxMedian>=10,]

counts_full <- counts_filtered

counts <- counts_full[,5:12]

sampleName <- sampleName_full[5:12]

treatmentName <- treatmentName_full[5:12]

###########################

##Exploratory Data Analysis

###########################

## Boxplot

colnames(counts) <- sampleName

counts_long <- gather(counts, Sample, Read_Count, SRR449:SRR456, factor_key=TRUE)

boxplot(log2(counts_long$Read_Count+1) ~ counts_long$Sample, las=2, xlab="",ylab="", col=c(2,2,2,2,3,3,3,3))


## Correlogram

color.scheme <- rev(brewer.pal(8,"RdBu"))

sample.cor <- cor(counts)

melt.sample.cor <- melt(sample.cor)

ggplot(data = melt.sample.cor, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+scale_fill_gradientn(colors=color.scheme, limits=c(0.98,1))+
  theme(axis.text.x=element_text(angle=90, hjust=1), axis.title.y=element_blank(), axis.title.x=element_blank())

## Principal Component Analysis

logCount <- log2(as.matrix(counts)+1)

pcaCount <- prcomp(t(logCount))

autoplot(pcaCount,label=TRUE, size=4, shape=3)

##################################

##Differential Expression Analysis 

##################################

############# DESeq2 method ##########################

treatment <- as.factor(c("Ctr","Ctr","Ctr","Ctr","HFD","HFD","HFD","HFD"))


meta <- data.frame(row.names=sampleName,treatment)

cds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ treatment)

head(cds)

cds <- estimateSizeFactors(cds)

cds <- estimateDispersions(cds)


plotDispEsts(cds, main="DESeq: Per-gene dispersion estimates")


## Significant of treatment effects

cds1 <- DESeq(cds)

cds1 <- nbinomWaldTest(cds1)

res<- results(cds1, pAdjustMethod = "BH", alpha =0.05)

res

summary(res)

sum(res$padj<0.05,na.rm=TRUE)

##volcano


plot(res$log2FoldChange, -log10(res$padj),  panel.first=grid(),
     main="Volcano plot", xlab=" log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="green")
abline(h=-log10(0.05), col="red")



## Comparison of the treatments: volcano plot

cds2<- DESeq(cds)

res1 <- results(cds2, alpha=0.05)

res1

summary(res1)

sum(res1$padj<0.05, na.rm=TRUE)

hist(res1$padj, main = "Distribition of Adjusted P-values", xlab = "Adjusted p-value")



## HFD versus Ctr

res1 <- results(cds2, contrast=c("treatment", "HFD", "Ctr"), alpha=0.05)
res1



resultsNames(cds2)
res1.shrink <- lfcShrink(cds2, coef="treatment_HFD_vs_Ctr", type='apeglm')


EnhancedVolcano(res1.shrink,
                lab = rownames(res1.shrink),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "Contrast: HFD versus Ctr",
                pCutoff = 4e-6,
                FCcutoff = 0.5,
                pointSize = 2.0,
                labSize = 3.0)

sum(res1.shrink$padj<0.05, na.rm=TRUE)



##Heatmap differentially expressed genes

color.scheme <- rev(brewer.pal(8,"RdYlBu"))

#DESeq2-normalized counts

cds_norm <- counts(cds2,normalized=T)

summary(cds_norm)

boxplot(cds_norm,las=2, ylim=c(0,300), col=c(2,2,2,2,3,3,3,3))



## Heatmap for log2(count+1) transformation

cds.norm <- normTransform(cds2)

cds.norm_df <- as.data.frame(assay(cds.norm))

summary(cds.norm_df)

boxplot(cds.norm_df)

##select <- order(rowMeans(cds_norm), decreasing=TRUE)[1:500]

df <- as.data.frame(colData(cds2)[,c("treatment")])

rownames(df) <- c('SRR449','SRR450',
                  'SRR451','SRR452','SRR453','SRR454','SRR455','SRR456')

heatmatrix <- assay(cds.norm)##[select,]

colnames(df) <- c("treatment")

pheatmap(heatmatrix, cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

## High gene count heatmap

select <- order(rowMeans(cds_norm), decreasing=TRUE)[1:30]

heatmatrix1 <- assay(cds.norm)[select,]

pheatmap(heatmatrix1, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)


##EdgeR method

### group comparison and p-values

dge = DGEList(counts=counts, group=treatment)

dge_norm = cpm(dge, normalized.lib.sizes=FALSE)

dge = calcNormFactors(dge)

dge$samples

plotMDS(dge)

dge = estimateCommonDisp(dge)
dge = estimateTagwiseDisp(dge)


fit.edgeR = exactTest(dge)

pval.edgeR = fit.edgeR$table$PValue

sum(pval.edgeR<0.05)

hist(fit.edgeR$table$PValue, main="Distribution of P-value", xlab ="P-Values")

##MA plot with 5% differentially expressed genes

res.fit.edgeR <- topTags(fit.edgeR, n=nrow(fit.edgeR$table))

plotSmear(fit.edgeR,
          de.tags=rownames(res.fit.edgeR)[which(res.fit.edgeR$table$FDR<0.05)])

## Volcano data 

volcanoData <- cbind(res.fit.edgeR$table$logFC, -log10(res.fit.edgeR$table$FDR))
colnames(volcanoData) <- c("logFC", "negLogPval")
head(volcanoData)
plot(volcanoData, panel.first=grid(), pch=20,cex=0.6, xlab=" log2(fold-change)", ylab="-log10(p-value)", main= "Edge R Volcano")
abline(v=0)
abline(v=c(-1,1), col="green")
abline(h=-log10(0.05), col="red")


####### log transformation of counts and visualization #####

dge_log <- log2(dge$counts+1)

head(dge_log)

boxplot(dge_log, col=c(2,2,2,2,3,3,3,3))

## heatmap for log transformation

pheatmap(dge_log, cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

## heatmap for gene count with large counts

select1 <- order(rowMeans(dge_log), decreasing=TRUE)[1:30]

heatmatrix2 <- dge_log[select1,]

pheatmap(heatmatrix2, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)


## heatmap for highly expressed Gene

upregulated <- res.fit.edgeR[which(res.fit.edgeR$table$FDR<0.05 & res.fit.edgeR$table$logFC>1), ] ## there are 31 upregulated

downregulated <- res.fit.edgeR[which(res.fit.edgeR$table$FDR<0.05 & res.fit.edgeR$table$logFC < -1), ] ## there are 30 downregulated

activeGene <- rbind(upregulated$table, downregulated$table)

actGeneName <- rownames(activeGene)

actNormCounts <- subset(dge_norm,rownames(dge_norm) %in% actGeneName)

actGene_cor <- cor(t(actNormCounts))

actGene_dist <- as.dist(1-actGene_cor)


## cluster

actGene_hclust <- hclust(actGene_dist,method ="complete")

plot(actGene_hclust,cex=0.6,main="HFD and Ctr")

geneGroup <- cutree(actGene_hclust,k=2)

actGeneGroup <- cbind(actNormCounts, geneGroup)

diffExp_down <- subset(actGeneGroup, geneGroup==1)[,1:8]

diffExp_up <- subset(actGeneGroup, geneGroup==2)[,1:8]

## heatmap for active normalized counts

bk1 <- c(seq(0,36.00,by=6))
bk2 <- c(36.01, seq(100, 1000, by=100), 1500, 10000)
bk <- c(bk1,bk2)  #combine the break limits for purpose of graphing

my_palette <- c(colorRampPalette(colors = c("lightblue", "darkblue"))(n = length(bk1)-1),
                c(colorRampPalette(colors = c("darkblue", "red"))(n = length(bk2)-1)))

pheatmap(actNormCounts, col= my_palette, breaks = bk, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)

## Compare DEseq versus EdgeR

venn1 <- venndetail(list(DESeq2=rownames(res1[which(res1$padj<0.05),]),edgeR=rownames(res.fit.edgeR$table[which(res.fit.edgeR$table$FDR<0.05),])))

plot(venn1)
#############################################

###Enrichment Functional Analysis

#############################################

org.Mm.egACCNUM2EG
Lkeys(org.Mm.egENSEMBL2EG)
keytypes(org.Mm.eg.db)
engkeys = keys(org.Mm.eg.db, keytype="ENSEMBL")
ensemble2entrez = AnnotationDbi::select(org.Mm.eg.db, 
                         keys=engkeys, 
                         columns=c("ENSEMBL","ENTREZID"), 
                         keytype = "ENSEMBL")

sigGene <- res.fit.edgeR[which(res.fit.edgeR$table$FDR<0.05),]$table

geneID <- sigGeneName <- rownames(sigGene)

sigGene1 <- cbind(sigGene, geneID)

idx_all <- ensemble2entrez$ENSEMBL %in% sigGeneName

idn_all <- ensemble2entrez$ENTREZID[idx_all]

gene.df_all <- bitr(idn_all, fromType="ENTREZID", toType=c("ENSEMBL","SYMBOL"), OrgDb=org.Mm.eg.db)

non_duplicates <- which(duplicated(gene.df_all$SYMBOL) == FALSE)

gene.df_all <- gene.df_all[non_duplicates, ] 

Gene_idn_all <- inner_join(sigGene1, gene.df_all, by= c("geneID"= "ENSEMBL"))

## Enricher for Gene Ontology

ggo <- groupGO(gene  = gene.df_all$ENTREZID,
                  OrgDb = org.Mm.eg.db,
                  ont   = "CC",
                  level = 3,
                  readable = TRUE)


ego_all_CC <- enrichGO(gene         = gene.df_all$ENSEMBL,
                    OrgDb         = org.Mm.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05)

ego_all_CC <- setReadable(ego_all_CC, OrgDb = org.Mm.eg.db)

ego_all_BP <- enrichGO(gene         = gene.df_all$ENSEMBL,
                    OrgDb         = org.Mm.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05)

ego_all_BP <- setReadable(ego_all_BP, OrgDb = org.Mm.eg.db)

ego_all_MF <- enrichGO(gene         = gene.df_all$ENSEMBL,
                       OrgDb         = org.Mm.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)

ego_all_MF <- setReadable(ego_all_MF, OrgDb = org.Mm.eg.db)

cluster_summary_CC <- data.frame(ego_all_CC)

cluster_summary_BP <- data.frame(ego_all_BP)

cluster_summary_MF <- data.frame(ego_all_MF)

dotplot(ego_all_BP)

dotplot(ego_all_CC)

dotplot(ego_all_MF)

###GSEA Method Using KEGG and GO

foldchanges <- Gene_idn_all$logFC

names(foldchanges) <- Gene_idn_all$ENTREZID

foldchanges <- sort(foldchanges, decreasing = TRUE)

head(foldchanges)

gseaKEGG <- gseKEGG(geneList = foldchanges, # ordered named vector of fold changes (Entrez IDs are the associated names)
                    organism = "mmu", # supported organisms listed below
                    #nPerm = 1000, # default number permutations
                    minGSSize = 10, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                    pvalueCutoff = 1, # padj cutoff value
                    verbose = FALSE)

gseaKEGG_results <- gseaKEGG@result

head(gseaKEGG_results)

gseaGO <- gseGO(geneList = foldchanges, 
                OrgDb = org.Mm.eg.db, 
                ont = 'BP', 
                #nPerm = 1000, 
                minGSSize = 10, 
                pvalueCutoff = 1,
                verbose = FALSE) 

gseaGO_results <- gseaGO@result

head(gseaGo_results)

gseaplot(gseaGO, geneSetID="GO:0002684")


###Enricher and GSEA using MSigDb
msigdbr_species()

m_df <- msigdbr(species = "Mus musculus")
head(m_df, 2) %>% as.data.frame

m_t2g <- msigdbr(species = "Mus musculus", category = "C3") %>% 
  dplyr::select(gs_name, entrez_gene)


em_all <- enricher(gene.df_all$ENTREZID, TERM2GENE=m_t2g)

head(em_all)

em2 <- GSEA(foldchanges, TERM2GENE = m_t2g)

head(gesa)

get_kegg_plots <- function(x) {
  pathview(gene.data = foldchanges, pathway.id = gseaKEGG_results$ID[x], species = "mmu", 
           limit = list(gene = 2, cpd = 1))
}

purrr::map(1:length(gseaKEGG_results$ID), get_kegg_plots)

##upregulated


idx_up <- ensemble2entrez$ENSEMBL %in% rownames(diffExp_up)

idn_up <- ensemble2entrez$ENTREZID[idx_up]

gene.df_up <- bitr(idn_up, fromType="ENTREZID", toType=c("ENSEMBL","SYMBOL"), OrgDb=org.Mm.eg.db)

non_duplicates <- which(duplicated(gene.df_up$SYMBOL) == FALSE)

gene.df_up <- gene.df_up[non_duplicates, ]

ego_up <- enrichGO(gene         = gene.df_up$ENSEMBL,
                    OrgDb         = org.Mm.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.1,
                    qvalueCutoff  = 0.1)

ego2_up <- enrichGO(gene         = gene.df_up$ENSEMBL,
                    OrgDb         = org.Mm.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.1,
                    qvalueCutoff  = 0.1)

ego3_up <- enrichGO(gene         = gene.df_up$ENSEMBL,
                    OrgDb         = org.Mm.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.1,
                    qvalueCutoff  = 0.1)

ego_up <- setReadable(ego_up, OrgDb = org.Mm.eg.db)

ego2_up <- setReadable(ego2_up, OrgDb = org.Mm.eg.db)

ego3_up <- setReadable(ego3_up, OrgDb = org.Mm.eg.db)

cluster_summary_1 <- data.frame(ego_up)

cluster_summary_2 <- data.frame(ego2_up)

cluster_summary_MF <- data.frame(ego3_up)

cluster_summary_1

cluster_summary_2

cluster_summary_MF

## downregulated

idx_down <- ensemble2entrez$ENSEMBL %in% rownames(diffExp_down)

idn_down <- ensemble2entrez$ENTREZID[idx_down]

gene.df_down <- bitr(idn_down, fromType="ENTREZID", toType=c("ENSEMBL","SYMBOL"), OrgDb=org.Mm.eg.db)

non_duplicates <- which(duplicated(gene.df_down$SYMBOL) == FALSE)

gene.df_down <- gene.df_down[non_duplicates, ]

ego_down <- enrichGO(gene         = gene.df_down$ENSEMBL,
                   OrgDb         = org.Mm.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.1,
                   qvalueCutoff  = 0.1)

ego2_down <- enrichGO(gene         = gene.df_down$ENSEMBL,
                    OrgDb         = org.Mm.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.1,
                    qvalueCutoff  = 0.1)

ego3_down <- enrichGO(gene         = gene.df_down$ENSEMBL,
                      OrgDb         = org.Mm.eg.db,
                      keyType       = 'ENSEMBL',
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.1,
                      qvalueCutoff  = 0.1)

ego_down <- setReadable(ego_down, OrgDb = org.Mm.eg.db)

ego2_down <- setReadable(ego2_down, OrgDb = org.Mm.eg.db)

ego3_down <- setReadable(ego3_down, OrgDb = org.Mm.eg.db)

cluster_summary_3 <- data.frame(ego_down)

cluster_summary_4 <- data.frame(ego2_down)

cluster_summary_MF_2 <- data.frame(ego3_down)

cluster_summary_3

cluster_summary_4

cluster_summary_MF_2

## Venn comparison

venn2 <- venndetail(list(upregulated=cluster_summary_1$ID,
                         downregulated=cluster_summary_3$ID))

plot(venn2)


