#---- Load Packages ------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(reshape2)
#---- Set Directory ------------------------------------------------------------
GEO <- 'GSE108097'
root_dir <- "/Users/yyang18/My Drive (yyang18@uth.edu)/Project/McGovern_Yanan_Cao/2022_10_10_Integration_Two_Batches/data/"
out_dir <- '/Users/yyang18/My Drive (yyang18@uth.edu)/Project/McGovern_Yanan_Cao/2022_10_10_Integration_Two_Batches/out/'
#---- Load and preprocess scRNA-seq data into Seurat ---------------------------
load_sc <- function(file)
{
  df <- read.csv(file = paste0(root_dir, GEO, '/' ,file))
  df <- df %>% 
    mutate(across('assigned_cluster', str_replace, 'Dendrtic cell', 'Dendritic cell')) # typo
  # count cell types
  cell_types <- as.data.frame(table(df$assigned_cluster))
  colnames(cell_types) <- c('cell type','number')
  cell_types <- cell_types[order(cell_types$number, decreasing = TRUE),]
  write.csv(cell_types,paste0(out_dir, 'integration/', 'GSM2906458_cell_types.csv'), row.names = FALSE)
  # merge sub-cell types
  df <- df %>% 
    mutate(across('assigned_cluster', str_replace, 'Endothelial cell_Fabp4 high', 'Endothelial cell')) %>% 
    mutate(across('assigned_cluster', str_replace, 'Endothelial cell_Lrg1 high', 'Endothelial cell')) %>% 
    mutate(across('assigned_cluster', str_replace, 'Endothelial cell_Tm4sf1 high', 'Endothelial cell')) %>% 
    mutate(across('assigned_cluster', str_replace, 'Erythroblast_Hbb-bt high', 'Erythroblast')) %>% 
    mutate(across('assigned_cluster', str_replace, 'Erythroblast_Igkc high', 'Erythroblast')) %>% 
    mutate(across('assigned_cluster', str_replace, 'Smooth muscle cell_Acta2 high', 'Smooth muscle cell')) %>% 
    mutate(across('assigned_cluster', str_replace, 'Smooth muscle cell_Rgs5 high', 'Smooth muscle cell')) %>% 
    mutate(across('assigned_cluster', str_replace, 'Stromal cell_Fn1 high', 'Stromal cell')) %>% 
    mutate(across('assigned_cluster', str_replace, 'Stromal cell_Mfap4 high', 'Stromal cell')) %>% 
    mutate(across('assigned_cluster', str_replace, 'Stromal cell_Smoc2 high', 'Stromal cell'))
  # create Seurat object
  count <- df %>% select(-c(barcode, assigned_cluster))
  row.names(count) <- df$barcode
  object <- CreateSeuratObject(counts = t(count))
  object@meta.data$clusters <- df$assigned_cluster
  return(object)
}
GSM2906458 <- load_sc('GSM2906458.csv')
sc_ref <- SCTransform(GSM2906458, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunTSNE(dims = 1:30) %>%
  RunUMAP(dims = 1:30)
#---- UMAP: GSE108097 ----------------------------------------------------------
umap <- DimPlot(sc_ref, reduction = "umap", group.by = "clusters", pt.size = 1.2,
                label = TRUE, label.size = 8, repel = TRUE) +
  ggtitle('') +
  theme(axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30), 
        legend.text=element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size=15)))
ggsave(paste0(out_dir, 'integration/figures/', 'GSM2906458_umap', '.tiff'), width=15, height=10, compression = "lzw")
#---- t-SNE: GSE108097 ---------------------------------------------------------
t_sne <- DimPlot(sc_ref, reduction = "tsne", group.by = "clusters", pt.size = 1.2,
                 label = TRUE, label.size = 8, repel = TRUE) +
  ggtitle('') +
  theme(axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30), 
        legend.text=element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size=15)))
ggsave(paste0(out_dir, 'integration/figures/', 'GSM2906458_tsne', '.tiff'), width=15, height=10, compression = "lzw")
#---- Load and preprocess spatial data into Seurat -----------------------------
visium <- readRDS(file = paste0(root_dir, 'Seurat/', 'seurat_object.rds'))
# Rename clusters
visium <- RenameIdents(visium,
                       "C0" = 'Acinar_C0',
                       "C1" = 'Acinar_C1',
                       "C2" = 'Acinar_C2',
                       "C3" = 'Stroma_C3',
                       "C4" = 'Acinar_C4',
                       "C5" = 'Stroma_C5',
                       "C6" = 'Islet_C6',
                       "C7" = 'Acinar_C7')
#---- Integrate spatial data with scRNA-seq ------------------------------------
anchors <- FindTransferAnchors(reference = sc_ref, query = visium, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = sc_ref$clusters, prediction.assay = TRUE,
                                  weight.reduction = visium[["pca"]], dims = 1:30)
visium[["predictions"]] <- predictions.assay
DefaultAssay(visium) <- "predictions"
#---- Save the predicted data --------------------------------------------------
# cellular barcode & cell type
df_cell_type <- as.data.frame(t(visium@assays$predictions@data))
df_cell_type$cell_type <- names(df_cell_type)[max.col(df_cell_type, ties.method = "first")]
df_cell_type <- df_cell_type[, c("max", "cell_type")]
# cellular barcode & cluster & sample
df_cluster <- data.frame(cluster = visium@active.ident)
df_cluster$sample <- c(rep('Young_AAP2', 768), rep('Aging_AAP2', 2096), rep('Young_AAP1', 1627), rep('Aging_AAP1', 1470))
# merge df_cell_type with df_cluster
df_predicted <- merge(df_cluster, df_cell_type, by = "row.names", all = TRUE)
colnames(df_predicted)[1] <- 'spatial_barcode'
colnames(df_predicted)[4] <- 'prediction_score'
df_predicted <- df_predicted[order(match(df_predicted$spatial_barcode, row.names(visium@meta.data))), ]
# add a column for merging clusters with same tissue types
df_predicted <- df_predicted %>%
  mutate(cluster_merge = sapply(strsplit(as.character(cluster), "_"), `[`, 1))
write.csv(df_predicted, paste0(out_dir, 'integration/', 'integration.csv'), row.names = FALSE)
#---- Bar Plot (Before Merging): percentage of cell types for all samples ------
df_predicted <- read.csv(paste0(out_dir, 'integration/', 'integration.csv'))
df_all <- melt(t(apply(as.data.frame.matrix(table(df_predicted$cluster, df_predicted$cell_type)),
                       1, function(x) round(100 * x / sum(x), 2))),
               varnames = c("cluster", "cell_type"),
               value.name = "percentage")
df_all$cluster <- factor(df_all$cluster, levels = c('Acinar_C0', 'Acinar_C1', 'Acinar_C2', 'Stroma_C3',
                                                    'Acinar_C4', 'Stroma_C5', 'Islet_C6', 'Acinar_C7'))
p <- ggplot(df_all, aes(x = cluster, y = percentage, fill = cell_type)) + 
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_manual(values=c("#EF7F49","#B79F00","#00B4F0","#C77CFF")) +
  theme(legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.key.size = unit(1, 'cm'),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(angle = 60, vjust = 0.5, size = 15),
        axis.text.y = element_text(size = 15)) +
  labs(fill = "Cell Type") +
  xlab('Cluster') +
  ylab('Percentage') 
ggsave(paste0(out_dir, 'integration/figures/', paste0('barplot_before_merging.tiff')), width=8, height=8, compression = "lzw")
#---- Bar Plot (After Merging): percentage of cell types for all samples -------
df_predicted <- read.csv(paste0(out_dir, 'integration/', 'integration.csv'))
df_all <- melt(t(apply(as.data.frame.matrix(table(df_predicted$cluster_merge, df_predicted$cell_type)),
                       1, function(x) round(100 * x / sum(x), 2))),
               varnames = c("cluster", "cell_type"),
               value.name = "percentage")
df_all$cluster <- factor(df_all$cluster, levels = c('Acinar', 'Stroma', 'Islet'))
p <- ggplot(df_all, aes(x = cluster, y = percentage, fill = cell_type)) + 
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_manual(values=c("#EF7F49","#B79F00","#00B4F0","#C77CFF")) +
  theme(legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.key.size = unit(1, 'cm'),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(angle = 60, vjust = 0.5, size = 15),
        axis.text.y = element_text(size = 15)) +
  labs(fill = "Cell Type") +
  xlab('Cluster') +
  ylab('Percentage') 
ggsave(paste0(out_dir, 'integration/figures/', paste0('barplot_after_merging.tiff')), width=8, height=8, compression = "lzw")
#---- Spatial Dim Plot ---------------------------------------------------------
visium@meta.data$cell_type <- df_predicted$cell_type
spatial_dim <- function(object)
{
  p <- SpatialDimPlot(object, group.by = 'cell_type')
  return(p)
}
p_1 <- spatial_dim(visium)[[1]] +
  scale_fill_manual(values = c("#EF7F49","#00B4F0","#C77CFF")) +
  labs(title = '', fill = "Cell Type") +
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15)) +
  guides(fill = guide_legend(override.aes = list(size=10)))
ggsave(paste0(out_dir, 'integration/', 'figures/', 'Young_AAP2_Cell_Types.tiff'), width=8, height=8, compression = "lzw")

p_2 <- spatial_dim(visium)[[2]] +
  scale_fill_manual(values=c("#EF7F49","#B79F00","#00B4F0","#C77CFF")) +
  labs(title = '', fill = "Cell Type") +
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15)) +
  guides(fill = guide_legend(override.aes = list(size=10)))
ggsave(paste0(out_dir, 'integration/', 'figures/', 'Aging_AAP2_Cell_Types.tiff'), width=8, height=8, compression = "lzw")

p_3 <- spatial_dim(visium)[[3]] +
  scale_fill_manual(values = c("#EF7F49","#B79F00","#00B4F0","#C77CFF")) +
  labs(title = '', fill = "Cell Type") +
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15)) +
  guides(fill = guide_legend(override.aes = list(size=10)))
ggsave(paste0(out_dir, 'integration/', 'figures/', 'Young_AAP1_Cell_Types.tiff'), width=8, height=8, compression = "lzw")

p_4 <- spatial_dim(visium)[[4]] +
  scale_fill_manual(values = c("#EF7F49","#B79F00","#00B4F0","#C77CFF")) +
  labs(title = '', fill = "Cell Type") +
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15)) +
  guides(fill = guide_legend(override.aes = list(size=10)))
ggsave(paste0(out_dir, 'integration/', 'figures/', 'Aging_AAP1_Cell_Types.tiff'), width=8, height=8, compression = "lzw")
