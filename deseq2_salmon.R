#---- set up libraries----------------------------------------------------------
library(readr)
library(tximport)
library(DESeq2)
library(edgeR)
library(dplyr)
library(io)

#---- set up parameters --------------------------------------------------------
project <- 'test' # name of project
condition_1 <- 'Smoke_KO'
condition_2 <- 'Smoke_Control'
cutoff_fdr <- 0.1
cutoff_logFC <- 1

#---- set up path --------------------------------------------------------------
path_tx2gene <- '/Users/yyang18/Google Drive/Data/tx2gene/' # path to tx2gene
path_input_data <- './input_data/' # path to input data
path_salmon <- paste0('./input_data/','salmon') # path to salmon
path_annotation <- '/Users/yyang18/Google Drive/Data/Gene_ID_Convert/' # path to gene annotation file
path_output_data <- './output_data/' # path to output data

#---- read files ---------------------------------------------------------------
tx2gene <- read.csv(paste0(path_tx2gene,'gencode.v37.annotation.csv'), row.names = 1)
group <- read.csv(paste0(path_input_data,'group.csv'), header=TRUE) # sample & condition
annotation = read.csv(paste0(path_annotation,'gencode_v37_gene_annotation.csv'))

#---- build ddsTxi -------------------------------------------------------------
rownames(group) <- group$sample
files <- file.path(path_salmon, group$sample, "quant.sf")
names(files) <- group$sample
print(all(file.exists(files)))
txi <- tximport(files, type="salmon", tx2gene = tx2gene)
print(dim(txi$counts)) # dimension of read counts
group$condition <- factor(group$condition)
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = group,
                                   design = ~ condition)
keep <- rowSums(counts(ddsTxi)) >= 10 # pre-filtering
ddsTxi <- ddsTxi[keep,]

#---- logcpm -------------------------------------------------------------------
y <- DGEList(counts = counts(ddsTxi), group = factor(group$condition)) %>%
     calcNormFactors()
logcpm <- cpm(y, log=TRUE)
row_name <- data.frame(gene_id = row.names(logcpm)) %>%
            merge(annotation)
row.names(logcpm) <- row_name$gene_name
write.csv(logcpm, paste0(path_output_data,project,'_logcpm.csv')) # logcpm matrix
logcpm_z <- t(scale(t(logcpm)))
write.csv(logcpm_z, paste0(path_output_data,project,'_logcpm_z.csv')) # logcpm z normalized matrix


#---- differential analysis ----------------------------------------------------
ddsTxi <- DESeq(ddsTxi)
res <- results(ddsTxi, contrast = c('condition',condition_1,condition_2))
res <- as.data.frame(res)
res$gene_id <- row.names(res)
res <- merge(res,annotation) %>%
       subset(select = c(1,8,9,2,3,4,5,6,7))
res$change = ifelse(res$padj < cutoff_fdr & abs(res$log2FoldChange) >= cutoff_logFC, 
                    ifelse(res$log2FoldChange > cutoff_logFC ,'Up','Down'),'Not Sig')
write.csv(res, paste0(path_output_data,project,'_deseq2.csv'), row.names = FALSE) # result of deseq2
gsea <- select(res,'gene_name','stat')
write.csv(gsea,paste0(path_output_data,project,'_gsea.csv'), row.names = FALSE) # input of gsea
res_diff = subset(res, change == 'Up' | change == 'Down')
write.csv(res_diff,paste0(path_output_data,project,'_diff.csv'), row.names=FALSE) # differentially expressed genes
heatmap = logcpm_z[res_diff$gene_name,]
write.csv(heatmap,paste0(path_output_data,project,'_heatmap.csv')) # input of heatmap 
res$gene_name[res$change == 'Not Sig'] <- ''
write.csv(res,paste0(path_output_data,project,'_volcano.csv'), row.names = FALSE) # input of volcano plot



