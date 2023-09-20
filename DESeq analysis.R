library(BiocManager)
#BiocManager::install("DESeq2")
#BiocManager::install("apeglm")
library(DESeq2)
library(apeglm)

#loading read count data 
count_matrix <- read.csv("/Users/maitreepatel/Desktop/R/RNASeq/count_matrix.csv",
                         header = TRUE,
                         row.names = 1)
View(count_matrix)
rownames(count_matrix)
#treatment are the columns
#genes are the rows

#loading in metadata
#if it is tab delimited, use read.table
metadata <- read.csv("/Users/maitreepatel/Desktop/R/RNASeq/design.csv",
                     header = TRUE)
#this gives information about the column
#while making this metadata file, it is important to note that the order os the sample is the same as that in data
#DESeq requires column names and the treatment information

#linking the two files using DESeq2
dds <- DESeqDataSetFromMatrix(count_matrix,
                              metadata,
                              ~Treatment) #function depends on data
#will be a large file

#filtering out genes with expression less than 10
filter.genes <- rowSums(counts(dds)) >= 10
#subsetting all rows with expressin >=10
dds <- dds[filter.genes,]

#running the main DESeq function for analysis for calculating differentially expressed genes
#performs normalization - accounting for different library size
ddsDE <- DESeq(dds)

#exporting normalized read counts
norm.counts <- counts(ddsDE,
                      normalized = TRUE)
#View(norm.counts)
write.csv(norm.counts,
          "normal_countmatrix.csv")

#getting DESeq results
result <- results(ddsDE, 
                  alpha = 0.05)
#checking overall results 
#summary(result)
#It measure a total of 11577 genes
#449 genes are upregulated
#415 genes are down regulated

#getting more detailed results
#reordering results using adjusted p-value
result.ordered <- result[order(result$padj),] 
write.csv(result.ordered,
          "deseq_results.csv")
#explanation:
#basemean is mean expression among the 7 treatments

#log2base tells us direction and intensity of change in expression

#seeing which which treatment is first:
resultsNames(ddsDE)
##gene expression is same between two means = get 1
##gene expression in untreated is greater that treated = between 0 and 1
##gene expression in treated is greater than untreated = greater than 1
##then log base two of the ratios above
##therefore lfc > 0 449 genes are expressed more in untreated
##lfc < 0 415 genes are expressed more in treated
#pvalue - each gene tested for significant gene expression (difference in mean between single end and paired end and treated and untreated?)
#corrected pvalue - adjusts for the fact that 5% of the pvalues show signifact difference by chance (so 5% of the 11000 genes will be shown to have significant expression)
#adjusted p value lowers the significance threshold, say, alpha = 0.001

#visualizing gene expression
plotMA(ddsDE,
       ylim = c(-5,5))
#+ve is more expressed for untreated
#-ve is more exressed for treated


