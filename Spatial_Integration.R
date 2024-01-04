#---- Load Packages ------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(glmGamPoi)
#---- Set Directory ------------------------------------------------------------
root_dir <- "/Users/yyang18/My Drive (yyang18@uth.edu)/Project/McGovern_Yanan_Cao/2022_10_10_Integration_Two_Batches/data/"
out_dir <- '/Users/yyang18/My Drive (yyang18@uth.edu)/Project/McGovern_Yanan_Cao/2022_10_10_Integration_Two_Batches/out/'
h5_mat_name <- './filtered_feature_bc_matrix.h5'
#---- Load data into Seurat and normalize --------------------------------------
load_data <- function(sample_id, slice, sample_name)
{
  setwd(paste0(root_dir,sample_id))
  object <- Seurat::Load10X_Spatial(
    data.dir = paste0(root_dir,sample_id), 
    filename = h5_mat_name,
    assay = 'Spatial',
    slice = slice, 
    filter.matrix = TRUE, 
    to.upper = FALSE)
  object <- SCTransform(object, assay='Spatial', verbose=FALSE, method="glmGamPoi", min_cells=1)
  object$orig.ident <- sample_name
  return(object)
}
YAP_3859 <- load_data('3859-YAP', 'YAP3859', 'Young_AAP2')
AAP_3956 <- load_data('3956-AAP', 'AAP3956', 'Aging_AAP2')
YAP_3858 <- load_data('3858-YAP', 'YAP3858', 'Young_AAP1')
AAP_3955 <- load_data('3955-AAP', 'AAP3955', 'Aging_AAP1')
#---- Merge samples ------------------------------------------------------------
st.list = list(YAP_3859 = YAP_3859, AAP_3956 = AAP_3956, 
               YAP_3858 = YAP_3858, AAP_3955 = AAP_3955)
options(future.globals.maxSize = 2000 * 1024^2)
st.features = SelectIntegrationFeatures(st.list, nfeatures = 3000, verbose = FALSE)
st.list <- PrepSCTIntegration(object.list = st.list, anchor.features = st.features,
                              verbose = TRUE)
int.anchors <- FindIntegrationAnchors(object.list = st.list, normalization.method = "SCT",
                                      verbose = FALSE, anchor.features = st.features)
sample_merge <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT",
                              verbose = FALSE)
#---- Dimensional reduction and clustering -------------------------------------
cluster <- function(object, resolution)
{
  object <- RunPCA(object, verbose = FALSE)
  object <- FindNeighbors(object, dims = 1:30)
  object <- FindClusters(object, verbose = FALSE, resolution = resolution)
  object <- RunUMAP(object, reduction = "pca", dims = 1:30)
  object <- RunTSNE(object, reduction = "pca", dims = 1:30)
  return(object)
}
sample_merge <- cluster(sample_merge, 0.3)
#---- Rename clusters: Reproduce -----------------------------------------------
levels(x = sample_merge) <- c("0", "3", "1", "4", "5", "2", "6", "7")
sample_merge <- RenameIdents(sample_merge, 
                             "0" = 'C0', 
                             "3" = 'C1', 
                             "1" = 'C2', 
                             "4" = 'C3', 
                             "5" = 'C4', 
                             "2" = 'C5', 
                             "6" = 'C6', 
                             "7" = 'C7')
#---- Save the Seurat object ---------------------------------------------------
saveRDS(sample_merge, file = paste0(root_dir, 'Seurat/', 'seurat_object.rds'))
#---- UMAP plots ---------------------------------------------------------------
umap <- function(object)
{
  p_1 <- DimPlot(object, reduction = "umap", group.by = c("ident"), pt.size = 1.5,
                 label = TRUE, label.size = 7) +
    ggtitle('') +
    theme(axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 30), 
          legend.text=element_text(size=20)) +
    guides(color = guide_legend(override.aes = list(size=20)))
  ggsave(paste0(out_dir, 'aggr/', 'aggr', '_umap_clusters.png'), width=13, height=10)
  
  p_2 <- DimPlot(object, reduction = "umap", group.by = c("orig.ident"), pt.size = 1.5) +
    ggtitle('') +
    theme(axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 30), 
          legend.text=element_text(size=20)) +
    guides(color = guide_legend(override.aes = list(size=20)))
  ggsave(paste0(out_dir, 'aggr/', 'aggr', '_umap_samples.png'), width=13, height=10)
}
umap(sample_merge)
#---- t-SNE plots --------------------------------------------------------------
t_sne <- function(object)
{
  p_1 <- DimPlot(object, reduction = "tsne", group.by = c("ident"), pt.size = 1.5,
                 label = TRUE, label.size = 7) +
    ggtitle('') +
    theme(axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 30), 
          legend.text=element_text(size=20)) +
    guides(color = guide_legend(override.aes = list(size=20)))
  ggsave(paste0(out_dir, 'aggr/', 'aggr', '_tsne_clusters.png'), width=13, height=10)
  
  p_2 <- DimPlot(object, reduction = "tsne", group.by = c("orig.ident"), pt.size = 1.5) +
    ggtitle('') +
    theme(axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 30), 
          legend.text=element_text(size=20)) +
    guides(color = guide_legend(override.aes = list(size=20)))
  ggsave(paste0(out_dir, 'aggr/', 'aggr', '_tsne_samples.png'), width=13, height=10)
}
t_sne(sample_merge)
#---- Spatial Dim Plot ---------------------------------------------------------
spatial_dim <- function(object)
{
  p <- SpatialDimPlot(object)
  return(p)
}
spatial_dim(sample_merge)[[1]] +
  scale_color_manual(values = c("#F8766D","#B79F00","#7CAE00","#00BA38","#00BFC4","#00B4F0","#C77CFF")) +
  labs(title = '', fill = "Clusters") +
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15)) +
  guides(fill = guide_legend(override.aes = list(size=10)))
ggsave(paste0(out_dir, 'aggr/', 'Young_AAP2_spatial_dim.png'), width=8, height=8)

spatial_dim(sample_merge)[[2]] +
  scale_color_manual(values = c("#F8766D","#B79F00","#7CAE00","#00BA38","#00BFC4","#00B4F0","#C77CFF")) +
  labs(title = '', fill = "Clusters") +
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15)) +
  guides(fill = guide_legend(override.aes = list(size=10)))
ggsave(paste0(out_dir, 'aggr/', 'Aging_AAP2_spatial_dim.png'), width=8, height=8)
spatial_dim(sample_merge)[[3]] +
  scale_color_manual(values = c("#F8766D","#B79F00","#7CAE00","#00BA38","#00BFC4","#00B4F0","#C77CFF")) +
  labs(title = '', fill = "Clusters") +
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15)) +
  guides(fill = guide_legend(override.aes = list(size=10)))
ggsave(paste0(out_dir, 'aggr/', 'Young_AAP1_spatial_dim.png'), width=8, height=8)

spatial_dim(sample_merge)[[4]] +
  scale_color_manual(values = c("#F8766D","#B79F00","#7CAE00","#00BA38","#00BFC4","#00B4F0","#C77CFF","#F564E3")) +
  labs(title = '', fill = "Clusters") +
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15)) +
  guides(fill = guide_legend(override.aes = list(size=10)))
ggsave(paste0(out_dir, 'aggr/', 'Aging_AAP1_spatial_dim.png'), width=8, height=8)
#---- Markers: one cluster vs. other clusters ----------------------------------
markers_all <- FindAllMarkers(sample_merge, logfc.threshold = 0, min.pct = 0)
for (i in split(markers_all, markers_all$cluster))
{
  write.csv(i, paste0(out_dir, 'aggr/', gsub(' ', '', paste0('markers_all_cluster_', unique(i$cluster), '.csv'))),
            row.names = FALSE)
  i_filter <- filter(i, avg_log2FC > 0.25, pct.1 > 0.1, pct.2 > 0.1, p_val_adj < 0.05)
  write.csv(i_filter, paste0(out_dir, 'aggr/', gsub(' ', '', paste0('markers_filter_cluster_', unique(i$cluster), '.csv'))),
            row.names = FALSE)
}

