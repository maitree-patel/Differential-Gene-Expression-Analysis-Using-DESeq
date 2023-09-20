#install.packages("pheatmap")
library(ggplot2)
library(pheatmap)

#loading in normalized count martix
norm.counts <- read.csv("/Users/maitreepatel/Desktop/R/RNASeq/normal_countmatrix.csv",
                        header = TRUE,
                        row.names = 1)

#loading deseq results (with the stats)
deseq.results <- read.csv("/Users/maitreepatel/Desktop/R/RNASeq/deseq_results.csv",
                          header = TRUE,
                          row.names = 1)

deseq.results$sig <- ifelse(deseq.results$padj <= 0.05, 
                            "yes",
                            "no")
deseq.results <- na.omit(deseq.results)

View(deseq.results)

library(tidyverse)

#plotting base mean on x and log2fold on y with a significant p-value or not
deseq.results %>%
  ggplot(aes(x = log10(baseMean),
             y = log2FoldChange,
             color = sig)) +
  geom_point() #changing base mean to log10
#tells us there are differentially expressed genes among read counts 
#in both the untreated (points above) and treated (below zero)

#deseq.results %>%
  #ggplot(aes(x = rownames(deseq.results),
             #y = baseMean)) +
  #geom_bar(stat = "identity")

#pheatmap
significant <- subset(deseq.results,
                      padj <= 0.05)
all.significant <- merge(norm.counts,
                         significant,
                         by = 0) #0 for merging by the row name
#a merged df with normalized read counts with all rows having a significant pvalue
#for pheat map we need only read counts
View(all.significant)
signi.counts <- all.significant[,2:8] #without rownames

#getting row names
row.names(signi.counts) <- all.significant$Row.names
View(signi.counts)

pheatmap(signi.counts)
#log10 values and so cant visualize it that well

pheatmap(log2(signi.counts + 1))
#high = red
#low = blue

#looking at where the data point falls in the distribution
#median read count of the row
pheatmap(log2(signi.counts + 1), 
         scale = "row")
#in short, within the row which genes are higher and lower

#removing rownames and clustering by the dendrograms
pheatmap(log2(signi.counts + 1), 
         scale = "row",
         show_rownames = FALSE,
         treeheight_row = FALSE,
         treeheight_col = FALSE)
#all the untreated samples show different expression than treated samples