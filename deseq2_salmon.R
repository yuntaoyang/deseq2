#---- set up libraries----------------------------------------------------------
library(readr)
library(tximport)
library(DESeq2)
library(edgeR)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)

#---- set up parameters --------------------------------------------------------
project <- 'test' # name of project
condition_1 <- 'Smoke_KO'
condition_2 <- 'Smoke_Control'
cutoff_fdr <- 0.05
cutoff_logFC <- 1
cutoff_gsea <- 0.25

#---- set up path --------------------------------------------------------------
path_tx2gene <- '/Users/yyang18/Google Drive/Data/tx2gene/' # path to tx2gene
path_input_data <- './input_data/' # path to input data
path_salmon <- paste0('./input_data/','salmon') # path to salmon
path_annotation <- '/Users/yyang18/Google Drive/Data/Gene_ID_Convert/' # path to gene annotation file
path_output_raw <- './output_raw/' # path to output raw data
path_output_report <- './output_report/' # path to output report data

#---- create output directory --------------------------------------------------
dir.create(path_output_raw)
dir.create(path_output_report)

#---- read files ---------------------------------------------------------------
tx2gene <- read.csv(paste0(path_tx2gene,'gencode.v37.annotation.csv'), row.names = 1)
group <- read.csv(paste0(path_input_data,'group.csv'), header=TRUE) # sample & condition
annotation = read.csv(paste0(path_annotation,'gencode_v37_gene_annotation.csv')) # gencode gene annotation

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
write.csv(logcpm, paste0(path_output_raw,project,'_logcpm.csv')) # logcpm matrix for output raw
logcpm_z <- t(scale(t(logcpm)))
write.csv(logcpm_z, paste0(path_output_raw,project,'_logcpm_z.csv')) # logcpm z normalized matrix for output raw

#---- differential analysis ----------------------------------------------------
ddsTxi <- DESeq(ddsTxi)
res <- results(ddsTxi, contrast = c('condition',condition_1,condition_2))
res <- as.data.frame(res)
res$gene_id <- row.names(res)
res <- merge(res,annotation) %>%
       subset(select = c(1,8,9,2,3,4,5,6,7))
res$change = ifelse(res$padj < cutoff_fdr & res$log2FoldChange > cutoff_logFC,'Up',
                    ifelse(res$padj < cutoff_fdr & res$log2FoldChange < -cutoff_logFC,'Down','NS'))
res$gene_id <- substr(res$gene_id,1,15)
write.csv(res, paste0(path_output_raw,project,'_deseq2.csv'), row.names = FALSE) # result of deseq2 for output raw
gsea <- dplyr::select(res,'gene_name','stat')
write.table(gsea,paste0(path_output_raw,project,'_gsea.rnk'), row.names = FALSE, sep='\t') # input of gsea for output raw
res_diff = subset(res, change == 'Up' | change == 'Down')
write.csv(res_diff,paste0(path_output_raw,project,'_diff.csv'), row.names=FALSE) # differentially expressed genes for output raw
heatmap = logcpm_z[res_diff$gene_name,]
write.csv(heatmap,paste0(path_output_raw,project,'_heatmap.csv')) # input of heatmap for output raw
res_vol <- res
res_vol$gene_name[res_vol$change == 'NS'] <- ''
write.csv(res_vol,paste0(path_output_raw,project,'_volcano.csv'), row.names = FALSE) # input of volcano plot for output raw

#---- output_report_DEG --------------------------------------------------------
res_deg <- res %>%
  filter(change=='Up'|change=='Down') %>%
  dplyr::select(-c('baseMean','lfcSE','stat'))
full_name <- AnnotationDbi::select(org.Hs.eg.db, keys=res_deg$gene_id, columns=c("GENENAME"), keytype="ENSEMBL")
colnames(full_name) <- c("gene_id","description")
res_deg <- merge(res_deg,full_name)
write.csv(res_deg,paste0(path_output_report,project,'_gene.csv'), row.names = FALSE) # differentially expressed genes for output report

#---- output_report_GO ---------------------------------------------------------
gene_go <- res$stat
names(gene_go) <- as.character(res$gene_name) #named vector
gene_go <- sort(gene_go,decreasing=T) #decreasing order
gseaGO <- gseGO(gene_go, OrgDb=org.Hs.eg.db,
                ont='BP',keyType="SYMBOL",pAdjustMethod = "BH",
                minGSSize=10, maxGSSize=500,eps=0,
                pvalueCutoff=1, verbose=FALSE, by="fgsea")
write.csv(gseaGO,paste0(path_output_raw,project,'_go.csv'), row.names = FALSE) # GO GSEA output raw
egseGO_filter <- gseaGO %>%
  filter(p.adjust < cutoff_gsea) %>%
  select(-c('rank','leading_edge'))
write.csv(egseGO_filter,paste0(path_output_report,project,'_go_report.csv'), row.names = FALSE) # GO GSEA output report

#---- output_report_KEGG -------------------------------------------------------
gene_tx <- bitr(names(gene_go),fromType="SYMBOL",toType=c("ENTREZID"),
                 OrgDb = org.Hs.eg.db)
colnames(gene_tx)[1] <- "gene_name"
gene_tx <- merge(gene_tx,res,by="gene_name")
gene_kegg <- gene_tx$stat
names(gene_kegg) <- as.character(gene_tx$ENTREZID) #named vector
gene_kegg <- sort(gene_kegg,decreasing=T) #decreasing order
egseKEGG <- gseKEGG(gene_kegg,organism='hsa',keyType="kegg",
                    minGSSize=10, maxGSSize=500,eps=0,
                    pvalueCutoff=1, pAdjustMethod = "BH",
                    verbose=FALSE, by="fgsea")
write.csv(egseKEGG,paste0(path_output_raw,project,'_kegg.csv'), row.names = FALSE) # KEGG GSEA output raw
egseKEGG_filter <- egseKEGG %>%
  filter(p.adjust < cutoff_gsea) %>%
  select(-c('rank','leading_edge'))
write.csv(egseKEGG_filter,paste0(path_output_report,project,'_kegg_report.csv'), row.names = FALSE) # KEGG GSEA output report

