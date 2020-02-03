#The purpose of this code is to perform differential expression analysis on
#RNAseq data on RNA extracted from the sea squirt Ciona robusta using the 
#DEseq2 library. The experimental condition involved knocking out the 
#transcription factor Pax3/7 using CRISPR-Cas9 compared to a wild type condition.

#Written by Jameson Orvis, adapted from code by Elijah K. Lowe.


# 1. Load all libraries ---------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("tximportData")
BiocManager::install("tximport")

BiocManager::install("DESeq2")

if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')

BiocManager::install("DRIMSeq")


library("tidyverse")
library("tximport")
library("readr")
library("tximportData")
library("DESeq2")
library(EnhancedVolcano)
library("DRIMSeq")
library(readxl)

# 2. Read in files --------------------------------------------------------

setwd("D:/quants/") #directory with quantified read data from salmon
dir <- "data/"

#reads file containing formatted list of every sample and associated metadata
samples <- read.table(file.path("samples.txt"), header=TRUE, fill = TRUE) 
samples$condition <- factor(rep(c("WT","Pax_cpr"),each=12))
rownames(samples) <- samples$assay

#files points to the quants file for each sample
files <- file.path(paste0(samples$assay), "quant.sf")
names(files) <- samples$assay

#Verifies that files exist for every sample
file.exists(files)

#A quant file saved as .csv. Used here as basically a list of all possible transcripts
quants_gene_list <- read_csv(paste0(dir,"quants_gene_list.csv"))

#cuts off the list of gene names in quants_gene_list after the second dot. This is to allow
#mapping of all transcripts to a single gene, since there can be multiple transcripts associated
#with a single gene. Tximport imports transcript data from raw data.
tx2geneKH <- data.frame(TXNAME = quants_gene_list$Name, GENEID = gsub("\\.v.+", "", quants_gene_list$Name))
txi.g <- tximport(files, type="salmon", tx2gene=tx2geneKH)
write_csv(tx2geneKH, paste0(dir,"ciona_rob_tx2gene.csv"))

# 3. Setup differential expression dataframes --------------------------------------------------

#Imports tximport data into DESeqDataSet
dds.g <- DESeqDataSetFromTximport(txi.g, samples, ~ rep + condition)
dds <- collapseReplicates(dds.g, dds.g$sample)

#boxplot of read counts for each sample
boxplot(log10(counts(dds)+1))
dds <- estimateSizeFactors(dds)
boxplot(log10(counts(dds,normalized=TRUE)+1)) #normalized read counts

#VarianceStabilizingTransformation essentially normalizes the data and 
#it to a homoskedastic distribution. Plots PCA.
vsd.g <- varianceStabilizingTransformation(dds.g)
pcaData <- plotPCA(vsd.g, intgroup=c("condition", "rep"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=rep)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
svg("all_gene-level_pca.svg")
dev.off()

#reorders conditions for wild type to be first.
dds.g$condition <- relevel(dds.g$condition, ref = "WT")

#actually performs differential gene analysis.
dds.g <- DESeq(dds.g)
res.g <- DESeq2::results(dds.g)
summary(res.g)

#writes results to a tibble
tib <- 
  res.g %>% 
  as_tibble() %>% 
  mutate(Kh.id = rownames(res.g))

#Creates a volcano plot of results
svg("Pax37_CRISPR_volcano_plot.svg")
EnhancedVolcano(res.g,
                lab = rownames(res.g),
                FCcutoff = 0.3, 
                pCutoff = 0.05,
                x = 'log2FoldChange',
                y = 'pvalue')
dev.off()

#Appends links containing in-situ hybridization data and other relevant information
#to the results spreadsheet
insitu_data <- read_excel("D:/quants/data/KHID-UniqueName-URLs-InSitu-COMPLETE.xlsx")
tib_appended <- tib %>%
  left_join(insitu_data, by = c("Kh.id" = "KHID"))
write_csv(as_data_frame(arrange(tib_appended, log2FoldChange)), paste0(dir,"Pax37_CRISPR_all_genes.csv"))
