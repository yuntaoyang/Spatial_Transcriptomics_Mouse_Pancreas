#---- Load Packages ------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(glmGamPoi)
#---- Set Directory ------------------------------------------------------------
root_dir <- "/Users/yyang18/My Drive (yyang18@uth.edu)/Project/McGovern_Yanan_Cao/2022_10_10_Integration_Two_Batches/data/"
out_dir <- '/Users/yyang18/My Drive (yyang18@uth.edu)/Project/McGovern_Yanan_Cao/2022_10_10_Integration_Two_Batches/out/'
#---- Load data into Seurat ----------------------------------------------------
sample_merge <- readRDS(file = paste0(root_dir, 'Seurat/', 'seurat_object.rds'))
#---- After rename clusters based on condition ---------------------------------
sample_merge <- RenameIdents(sample_merge, 
                             "C0" = 'Acinar', 
                             "C1" = 'Acinar', 
                             "C2" = 'Acinar', 
                             "C3" = 'Stroma', 
                             "C4" = 'Acinar', 
                             "C5" = 'Stroma', 
                             "C6" = 'Islet', 
                             "C7" = 'Acinar')
sample_merge$sample_type <- Idents(sample_merge)
sample_merge$sample_type.ident <- paste(Idents(sample_merge), sample_merge$orig.ident, sep = "_")
Idents(sample_merge) <- "sample_type.ident"
#---- Count spots for each cluster ---------------------------------------------
cluster_count <- as.data.frame(table(sample_merge$sample_type.ident))
colnames(cluster_count) <- c('cluster','count')
write.csv(cluster_count, paste0(out_dir, 'Aging_vs_Young/', 'cluster_count.csv'), row.names = FALSE)
#---- Identification of Spatially Variable Features ----------------------------
# differential expression between any clusters
de_markers <- function(object, cluster_1, cluster_2, file_name)
{
  markers <- FindMarkers(object, ident.1 = cluster_1, ident.2 = cluster_2,
                         logfc.threshold = 0, min.pct = 0)
  write.csv(markers, paste0(out_dir, 'Aging_vs_Young/', 'aggr_', file_name, '.csv'))
  markers_filter <- markers  %>%
    filter(p_val_adj < 0.05, avg_log2FC > 0.25 | avg_log2FC < -0.25, pct.1 > 0.1, pct.2 > 0.1)
  write.csv(markers_filter, paste0(out_dir, 'Aging_vs_Young/', 'aggr_', file_name, '_filter.csv'))
}
# All clusters in Aging vs. All clusters in Young
YAP_clustes <- c('Acinar_Young_AAP1', 'Acinar_Young_AAP2', 'Islet_Young_AAP1', 'Islet_Young_AAP2', 'Stroma_Young_AAP1', 'Stroma_Young_AAP2')
AAP_clustes <- c('Acinar_Aging_AAP1', 'Acinar_Aging_AAP2', 'Islet_Aging_AAP1', 'Islet_Aging_AAP2', 'Stroma_Aging_AAP1', 'Stroma_Aging_AAP2')
de_markers(sample_merge, AAP_clustes, YAP_clustes, 'AAP_vs_YAP')
# Cluster n in Aging (AAP_1 & AAP_2) vs. Cluster n in Young (YAP_1 & YAP_2)
de_markers(sample_merge, c('Acinar_Aging_AAP1', 'Acinar_Aging_AAP2'), c('Acinar_Young_AAP1', 'Acinar_Young_AAP2'), 'Acinar_AAP_vs_YAP')
de_markers(sample_merge, c('Islet_Aging_AAP1', 'Islet_Aging_AAP2'), c('Islet_Young_AAP1', 'Islet_Young_AAP2'), 'Islet_AAP_vs_YAP')
de_markers(sample_merge, c('Stroma_Aging_AAP1', 'Stroma_Aging_AAP2'), c('Stroma_Young_AAP1', 'Stroma_Young_AAP2'), 'Stroma_AAP_vs_YAP')
#---- Dotplot ------------------------------------------------------------------
# Generate data
degs <- read.csv(paste0(out_dir, 'Aging_vs_Young/', 'Selected_DEG.csv'))
load_data <- function(file_name,cluster)
{
  df <- read.csv(paste0(out_dir, 'Aging_vs_Young/', file_name)) %>%
    select(X, avg_log2FC, p_val_adj) %>%
    filter(X %in% degs$genes) %>%
    rename('genes' = 'X', 'FDR' = 'p_val_adj', 'Log2FC' = 'avg_log2FC') %>%
    arrange(factor(genes, levels = degs$genes))
  df$type <- rep(cluster,31)
  return(df)
}
aggr_all <- load_data('aggr_AAP_vs_YAP.csv', 'All Clusters')
aggr_Acinar <- load_data('aggr_Acinar_AAP_vs_YAP.csv', 'Acinar')
aggr_Islet <- load_data('aggr_Islet_AAP_vs_YAP.csv', 'Islet')
aggr_Stroma <- load_data('aggr_Stroma_AAP_vs_YAP.csv' ,'Stroma')
data <- rbind(aggr_all, aggr_Acinar, aggr_Islet, aggr_Stroma) %>%
  filter(FDR < 0.05, Log2FC > 0.25 | Log2FC < -0.25)
data$type <- factor(data$type, levels = c('All Clusters', 'Acinar', 'Stroma', 'Islet'))
data$genes <- factor(data$genes, levels = rev(unique(data$genes)))
write.csv(data, paste0(out_dir, 'Aging_vs_Young/', 'dotplot.csv'), row.names = FALSE)
# Generate the figure
p <- ggplot(data, aes(x = type, y = genes, color = Log2FC, size = -log10(FDR))) + 
  geom_point() +
  theme_bw() +
  labs(x = 'Aging_AAP vs. Young_AAP', y = 'Differentially Expressed Genes') +
  scale_size_continuous(breaks=c(5, 10, 50, 100)) +
  scale_color_gradient2(low = "blue4", mid = 'white', high = "red2", midpoint = 0) +
  theme(legend.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        legend.key.size = unit(1.5, 'cm'),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x = element_text(angle = 60, vjust = 0.5, size = 20),
        axis.text.y = element_text(size = 20))
ggsave(paste0(out_dir, 'Aging_vs_Young/', 'dotplot.tiff'), width=10, height=12, compression = "lzw")
