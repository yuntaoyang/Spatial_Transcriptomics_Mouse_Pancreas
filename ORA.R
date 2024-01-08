#---- Load Packages ------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(stringr)
#---- Set Directory ------------------------------------------------------------
in_dir <- '/Users/yyang18/My Drive (yyang18@uth.edu)/Project/McGovern_Yanan_Cao/2022_10_10_Integration_Two_Batches/out/Aging_vs_Young/'
out_dir <- '/Users/yyang18/My Drive (yyang18@uth.edu)/Project/McGovern_Yanan_Cao/2022_10_10_Integration_Two_Batches/out/pathways/'
#---- Define Functions ---------------------------------------------------------
# GO-BP
GO <- function(file)
{
  df <- read.csv(paste0(in_dir, file, '.csv'), row.names = 'X') %>%
    filter(pct.1 > 0.1, pct.2 > 0.1)
  df$gene <- row.names(df)
  # up-regulated
  df_up <- df %>%
    filter(p_val_adj < 0.05, avg_log2FC > 0.25)
  ego <- enrichGO(gene = df_up$gene,
                  keyType = "SYMBOL",
                  OrgDb = org.Mm.eg.db,
                  ont = 'BP',
                  minGSSize = 10,
                  maxGSSize = 500,
                  pAdjustMethod = "BH",
                  pvalueCutoff = 1)
  write.csv(ego, paste0(out_dir, 'file/GO/', file, '_GO_Up.csv'), row.names = FALSE)
  ego_filter <- ego %>%
    filter(qvalue < 0.05)
  write.csv(ego_filter, paste0(out_dir, 'file/GO/', file, '_GO_Up_filter.csv'), row.names = FALSE)
  # down-regulated
  df_down <- df %>%
    filter(p_val_adj < 0.05, avg_log2FC < -0.25)
  ego <- enrichGO(gene = df_down$gene,
                  keyType = "SYMBOL",
                  OrgDb = org.Mm.eg.db,
                  ont = 'BP',
                  minGSSize = 10,
                  maxGSSize = 500,
                  pAdjustMethod = "BH",
                  pvalueCutoff = 1)
  write.csv(ego, paste0(out_dir, 'file/GO/', file, '_GO_Down.csv'), row.names = FALSE)
  ego_filter <- ego %>%
    filter(qvalue < 0.05)
  write.csv(ego_filter, paste0(out_dir, 'file/GO/', file, '_GO_Down_filter.csv'), row.names = FALSE)
}
# KEGG
KEGG <- function(file)
{
  df <- read.csv(paste0(in_dir, file, '.csv'), row.names = 'X') %>%
    filter(pct.1 > 0.1, pct.2 > 0.1)
  df$gene <- row.names(df)
  ids <- bitr(df$gene, 'SYMBOL', 'ENTREZID', 'org.Mm.eg.db')
  data <- merge(df,ids, by.x='gene', by.y='SYMBOL')
  # up-regulated
  data_up <- data %>%
    filter(p_val_adj < 0.05, avg_log2FC > 0.25)
  ekegg <- enrichKEGG(gene = data_up$ENTREZID,
                      organism = "mmu",
                      keyType = "ncbi-geneid",
                      minGSSize = 10,
                      maxGSSize = 500,
                      pAdjustMethod = "BH",
                      pvalueCutoff = 1)
  ekegg <- as.data.frame(ekegg)
  ekegg$geneID <- as.character(ekegg$geneID)
  gene_name = c()
  for(i in ekegg$geneID)
  {
    gene_list = c()
    for (j in as.list(strsplit(i,'/'))[[1]])
    {
      gene_list <- c(gene_list,paste(data[data$ENTREZID==j,]$gene,collapse='/'))
    }
    gene_name <- c(gene_name,paste(gene_list,collapse='/'))
  }
  ekegg <- subset(ekegg, select = -c(geneID))
  ekegg['gene_name'] <- gene_name
  write.csv(ekegg, paste0(out_dir, 'file/KEGG/', file, '_KEGG_Up.csv'), row.names = FALSE)
  ekegg_filter <- ekegg %>%
    filter(qvalue < 0.05)
  write.csv(ekegg_filter, paste0(out_dir, 'file/KEGG/', file, '_KEGG_Up_filter.csv'), row.names = FALSE)
  # down-regulated
  data_down <- data %>%
    filter(p_val_adj < 0.05, avg_log2FC < -0.25)
  ekegg <- enrichKEGG(gene = data_down$ENTREZID,
                      organism = "mmu",
                      keyType = "ncbi-geneid",
                      minGSSize = 10,
                      maxGSSize = 500,
                      pAdjustMethod = "BH",
                      pvalueCutoff = 1)
  ekegg <- as.data.frame(ekegg)
  ekegg$geneID <- as.character(ekegg$geneID)
  gene_name = c()
  for(i in ekegg$geneID)
  {
    gene_list = c()
    for (j in as.list(strsplit(i,'/'))[[1]])
    {
      gene_list <- c(gene_list,paste(data[data$ENTREZID==j,]$gene,collapse='/'))
    }
    gene_name <- c(gene_name,paste(gene_list,collapse='/'))
  }
  ekegg <- subset(ekegg, select = -c(geneID))
  ekegg['gene_name'] <- gene_name
  write.csv(ekegg, paste0(out_dir, 'file/KEGG/', file, '_KEGG_Down.csv'), row.names = FALSE)
  ekegg_filter <- ekegg %>%
    filter(qvalue < 0.05)
  write.csv(ekegg_filter, paste0(out_dir, 'file/KEGG/', file, '_KEGG_Down_filter.csv'), row.names = FALSE)
}
#---- Run ORA ------------------------------------------------------------------
data <- c('aggr_AAP_vs_YAP', 'aggr_Acinar_AAP_vs_YAP', 'aggr_Islet_AAP_vs_YAP', 'aggr_Stroma_AAP_vs_YAP')
# GO
lapply(data, GO)
# KEGG
lapply(data, KEGG)
#---- GO - Dot Plot: top 20 up-regulated pathways ------------------------------
# Generate data
pathways <- read.csv(paste0(out_dir, 'Top_Pathways.csv'))
colnames(pathways)[1] <- "Description"
pathways$GO <- paste(seq_along(pathways$Description), pathways$Description, sep = ". ")
# Load the GO data
load_data <- function(file_name, cluster)
{
  df <- read.csv(paste0(out_dir, '/file/GO/', file_name)) %>%
    dplyr::select(Description, qvalue, Count) %>%
    filter(Description %in% pathways$Description) %>%
    left_join(pathways, by = 'Description') %>%
    rename('qvalue' = 'FDR') %>%
    arrange(factor(GO, levels = pathways$GO))
  df$type <- rep(cluster, dim(df)[1])
  return(df)
}
aggr_all <- load_data('aggr_AAP_vs_YAP_GO_Up.csv', 'All Clusters')
aggr_Acinar <- load_data('aggr_Acinar_AAP_vs_YAP_GO_Up.csv', 'Acinar')
aggr_Stroma <- load_data('aggr_Stroma_AAP_vs_YAP_GO_Up.csv', 'Stroma')
aggr_Islet <- load_data('aggr_Islet_AAP_vs_YAP_GO_Up.csv', 'Islet')
data <- rbind(aggr_all, aggr_Acinar, aggr_Stroma, aggr_Islet) %>%
  filter(FDR < 0.05)
data$type <- factor(data$type, levels = c('All Clusters', 'Acinar', 'Stroma', 'Islet'))
data$GO <- factor(data$GO, levels = rev(pathways$GO))
colnames(data)[3] <- "size"
write.csv(data, paste0(out_dir, 'figure/', 'dotplot_top20_pathways.csv'), row.names = FALSE)
# Generate the figure
p <- ggplot(data, aes(x = type, y = GO, color = -log10(FDR), size = size)) + 
  geom_point() +
  theme_bw() +
  labs(x = expression("Aging_AAP vs. Young_AAP"), y = 'GO') +
  scale_size_continuous(breaks=c(10,30,50)) +
  scale_color_gradient(low = "blue", high = "red", breaks = c(5,10,15,20,25)) +
  theme(legend.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        legend.key.size = unit(1.5, 'cm'),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x = element_text(angle = 60, vjust = 0.5, size = 20),
        axis.text.y = element_text(size = 20))
ggsave(paste0(out_dir, 'figure/', 'dotplot_top20_pathways.tiff'), width=16, height=12, compression = "lzw")
#---- GO - Gene-Concept Network ------------------------------------------------
# GO: All Clusters
df <- read.csv(paste0(in_dir, 'aggr_AAP_vs_YAP', '.csv'), row.names = 'X') %>%
  filter(pct.1 > 0.1, pct.2 > 0.1)
df$gene <- row.names(df)
# up-regulated
df_up <- df %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.25)
ego <- enrichGO(gene = df_up$gene,
                keyType = "SYMBOL",
                OrgDb = org.Mm.eg.db,
                ont = 'BP',
                minGSSize = 10,
                maxGSSize = 500,
                pAdjustMethod = "BH",
                pvalueCutoff = 1)
write.csv(head(ego, 12), paste0(out_dir, 'figure/', 'gene_concept_network_top12.csv'), row.names=FALSE)
df_up_select <- df_up %>%
  filter(avg_log2FC > 2)
geneList <- df_up_select$avg_log2FC
names(geneList) <- df_up_select$gene
# only keep genes with avg_log2FC > 2
n = 1
for(i in ego@result$geneID)
{
  genes <- ''
  for (j in as.list(strsplit(i,'/'))[[1]])
  {
    if (j %in% names(geneList))
    {
      genes <- paste0(genes,'/',j)
    }
  }
  ego@result$geneID[n] <- substring(genes,2)
  ego@result$Description[n] <- n
  n = n + 1
}
# Gene-Concept Network - top 12 pathways
p <- cnetplot(ego, color.params = list(foldChange = geneList), showCategory = 12,
              cex.params = list(category_label = 2, gene_label = 1.5)) +
  scale_colour_gradient(low = "yellow", high = "red", name = "log2FC", breaks <- c(3,4)) +
  theme_bw() +
  theme(legend.text = element_text(size = 35),
        legend.title = element_text(size = 35),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank())
ggsave(paste0(out_dir, 'figure/', 'gene_concept_network_top12.tiff'), width=18, height=14, compression = "lzw")
