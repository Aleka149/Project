BiocManager::install("ggplot2")
BiocManager::install("ggrepel")
BiocManager::install("gplots")
BiocManager::install("dplyr")
BiocManager::install("tidyverse")
BiocManager::install("edgeR")
BiocManager::install("gplots")
BiocManager::install("ggfortify")
BiocManager::install("clusterprofiler")
BiocManager::install("org.HS.eg.db")
BiocManager::install("pheatmap")
library(ggplot2)
library(gplots)
library(ggrepel)
library(tidyverse)
library(ggfortify)
library(pheatmap)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(gplots)
set.seed(777)
# Read the data
IDs <- read.table("IDs.tabular" ,header=TRUE)
View(IDs)
normcounts <- read.table("edgeR_normcounts.tabular" ,header=TRUE)
View(normcounts)
# Format the data
class(IDs$GeneID) #integer
class(normcounts$ENTREZID) #integer also
Format <- merge(IDs,normcounts,by.x="GeneID",by.y="ENTREZID")
View(Format)

sum(is.na(Format)) # after running we have 80 na
Format <- na.omit(Format)
sum(is.na(Format))
rownames(Format) <- Format$SYMBOL #fixes rownames now that there's no NA's
View(Format) #loads table to updated form
Format <- Format[,-c(1,10)] #removes columns named SYMVBOL and GeneID
view(Format)

# metadata
md <- data.frame(colnames(Format))
colnames(md)[1] <- "SampleID"
md$group <- rep(c("control", "cancer"),times=c(4,4))
md$study <- rep(c("Roels", "Verboom"), times=c(4,4))

  
# Create pca as asked
count_pca <- prcomp(t(Format))

autoplot(count_pca)
autoplot(count_pca, data=md, colour="group", size=4)

autoplot(count_pca, data=md, colour="group", size=4, frame=TRUE, frame.type = 'norm')+
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
# the pca plot tells if the edgeR is done and if it was upregulated or downregulated
#
#

#t-test + adjusting pvalues

pvalue <-data.frame(p_value=rep(0,nrow(Format))) # doing the t-test and for
rownames(pvalue) <-rownames(Format)
for (i in 1:nrow(Format)){
  pvalue$p_value[i] <-t.test(Format [i,1:4], Format [i,5:8], alternative = "two.sided")$p.value
  
}
View(pvalue) # 13323
pvalue$padjust <- p.adjust(pvalue$p_value, method = "fdr") # Benjamini and Hochberg

head(pvalue)

Format$pval<-pvalue$p_value
count_subset_sig <- Format[Format$pval < 0.05,]

count_subset_sig$ControlMean<-rowMeans(count_subset_sig[,1:4])
count_subset_sig$PatientMean<-rowMeans(count_subset_sig[,5:8])

head(count_subset_sig)
#
#

#log2fc 

FC_ctrl <- rowMeans(Format[,1:4])
FC_cancer <- rowMeans(Format[,5:8])
FC <- as.data.frame(FC_cancer-FC_ctrl)
colnames(FC) <- "FoldChange" #log2FC
FDR_pvalue <- pvalue
FDR_pvalue$log2FC <- FC$FoldChange # we add the logFC column to our FDR pvalue matrix
padjusted <- pvalue

# visualization 1

universe <- data.frame(gene=rownames(Format)) # create a data frame with the names of all measured genes
universe <- bitr(universe$gene, fromType = "SYMBOL",toType = "ENTREZID", OrgDb = "org.Hs.eg.db" ) # convert symbol to ENTREZID

 table(duplicated(universe$ENTREZID)) #check if duplicated entrez ID
universe <- data.frame(gene=unique(universe$ENTREZID)) # remove duplicated ids if present


# Define input parameters
sig_genes_pathway <- rownames(Format) # create a data frame with the names of all signficiant genes (FDR<0.05)

genes_entrez <- bitr(sig_genes_pathway, fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)  # convert symbol to entrez

# table(duplicated(genes_entrez$ENTREZID)) # remove duplicated genes


# KEGG pathway analysis
kegg1 <- enrichKEGG(gene= genes_entrez$ENTREZID,
                    organism = "hsa",
                    keyType = "kegg",
                    universe = universe$ENTREZID, 
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH")


kegg1 <- setReadable(kegg1, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

dotplot(kegg1, showCategory=5, font.size=15) 

summary_enrich <- kegg1@result

# GO pathway analysis
go1 <- enrichGO(gene= genes_entrez$ENTREZID, 
                OrgDb = 'org.Hs.eg.db', 
                ont = "BP", 
                pvalueCutoff = 0.05, 
                pAdjustMethod = "BH", 
                universe = universe$ENTREZID, 
                qvalueCutoff = 0.2, 
                minGSSize = 10, 
                maxGSSize = 500, readable = TRUE)

dotplot(go1, showCategory=10, font.size=15) 

summary_go <- go1@result

# visualization 2

# making the mean values into a long format
stackbar<-gather(long_subset, key="group", value="mean", -c(ind))

#plotting this with geom_col
ggplot(stackbar, aes(x=ind, y=mean, fill=group))+
  geom_col(position="dodge")+
  facet_wrap(~ind)
# I also wish to plot the mean values as stack bars
# visualization 3


p <- FDR_pvalue %>%
  ggplot(mapping = aes(x=log2FC, y=padjust))+
  geom_point()+
  theme_minimal()
p2 <- p +
  geom_hline(yintercept = c(-2,2), col="red") +
  geom_vline(xintercept= c(-2,2), col="red")
# Color the dots based on their differential expression
FDR_pvalue$diffexpressed <- "NO" # create a new column in the data frame and fill it with NO
FDR_pvalue$diffexpressed[FDR_pvalue$log2FC > 2 & FDR_pvalue$padjusted < 0.05] <- "UP"
FDR_pvalue$diffexpressed[FDR_pvalue$log2FC < -2 & FDR_pvalue$padjusted < 0.05] <- "DOWN"
p <- FDR_pvalue %>%
  ggplot(mapping=aes(x=log2FC, y=padjust), col=diffexpressed)+
  geom_point()+
 # theme_minimal()
p
p2 <- p +
  geom_hline(yintercept = log10(0.05), col="red")+
  geom_vline(xintercept= c(-2,2), col="red")
mycolors <- c("blue", "red", "black") # specifiy the colors you want to use
names(mycolors) <- c("DOWN", "UP", "NO") # add column names
p2
p3 <- p2+
scale_color_manual(values=mycolors)
