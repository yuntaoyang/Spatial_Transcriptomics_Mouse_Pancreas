#---- Load Packages ------------------------------------------------------------
library(Seurat)
library(ggplot2)
#---- Set Directory ------------------------------------------------------------
root_dir <- "/Users/yyang18/My Drive (yyang18@uth.edu)/Project/McGovern_Yanan_Cao/2022_10_10_Integration_Two_Batches/data/"
out_dir <- '/Users/yyang18/My Drive (yyang18@uth.edu)/Project/McGovern_Yanan_Cao/2022_10_10_Integration_Two_Batches/out/'
#---- Load the Seurat object ---------------------------------------------------
sample_merge <- readRDS(file = paste0(root_dir, 'Seurat/', 'seurat_object.rds'))
#---- Rename clusters ----------------------------------------------------------
sample_merge <- RenameIdents(sample_merge, 
                             "C0" = 'Acinar_C0', 
                             "C1" = 'Acinar_C1', 
                             "C2" = 'Acinar_C2', 
                             "C3" = 'Stroma_C3', 
                             "C4" = 'Acinar_C4', 
                             "C5" = 'Stroma_C5', 
                             "C6" = 'Islet_C6', 
                             "C7" = 'Acinar_C7')
#---- Spatial Feature Plot -----------------------------------------------------
spatial_feature <- function(object, feature, min, max)
{
  p <- SpatialFeaturePlot(object, features = c(feature), image.alpha = 0)
  p_1 <- p[[1]] + 
    scale_fill_gradientn(colours=c('blue', 'white', 'red'),limits=c(min, max),oob=scales::squish) + 
    ggtitle('') +
    theme(legend.text=element_text(size=30), legend.title=element_text(size=30), legend.key.size = unit(2, 'cm'))
  ggsave(paste0(out_dir, 'markers/', feature, '/Young_AAP2_spatial_', feature, '.tiff'), width=8, height=8, compression = "lzw")
  
  p_2 <- p[[2]] + 
    scale_fill_gradientn(colours=c('blue', 'white', 'red'),limits=c(min, max),oob=scales::squish) + 
    ggtitle('') +
    theme(legend.text=element_text(size=30), legend.title=element_text(size=30), legend.key.size = unit(2, 'cm'))
  ggsave(paste0(out_dir, 'markers/', feature, '/Aging_AAP2_spatial_',feature,'.tiff'), width=8, height=8, compression = "lzw")
  
  p_3 <- p[[3]] + 
    scale_fill_gradientn(colours=c('blue', 'white', 'red'),limits=c(min, max),oob=scales::squish) + 
    ggtitle('') +
    theme(legend.text=element_text(size=30), legend.title=element_text(size=30), legend.key.size = unit(2, 'cm'))
  ggsave(paste0(out_dir, 'markers/', feature, '/Young_AAP1_spatial_',feature,'.tiff'), width=8, height=8, compression = "lzw")
  
  p_4 <- p[[4]] + 
    scale_fill_gradientn(colours=c('blue', 'white', 'red'),limits=c(min, max),oob=scales::squish) + 
    ggtitle('') +
    theme(legend.text=element_text(size=30), legend.title=element_text(size=30), legend.key.size = unit(2, 'cm'))
  ggsave(paste0(out_dir, 'markers/', feature, '/Aging_AAP1_spatial_',feature,'.tiff'), width=8, height=8, compression = "lzw")
}
#---- Violin plot --------------------------------------------------------------
vlnplot <- function(object, feature)
{
  VlnPlot(object, features = c(feature)) + 
    xlab('Clusters') + 
    theme(legend.text = element_text(size = 25),
          legend.key.size = unit(1.5, 'cm'),
          axis.title.x = element_text(size = 25),
          axis.title.y = element_text(size = 25),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          plot.title = element_text(size = 25))
  ggsave(paste0(out_dir, 'markers/', feature, '/aggr_violin_', feature, '.tiff'), width=10, height=6, compression = "lzw")
}
#---- Feature Plot -------------------------------------------------------------
feature_plot <- function(object,feature)
{
  FeaturePlot(object, features = c(feature), 
              label = TRUE, label.size = 6, cols = c("lightgrey" ,"#DE1F1F")) +
    theme(legend.text = element_text(size = 25),
          legend.key.size = unit(1.5, 'cm'),
          axis.title.x = element_text(size = 25),
          axis.title.y = element_text(size = 25),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          plot.title = element_text(size = 25))
  ggsave(paste0(out_dir, 'markers/', feature, '/aggr_feature_', feature, '.tiff'), width=8, height=7, compression = "lzw")
}
#---- Gene Plot ----------------------------------------------------------------
gene_plot <- function(object, feature, min, max)
{
  spatial_feature(object, feature, min, max)
  vlnplot(object, feature)
  feature_plot(object, feature)
}
# Figure 3
gene_plot(sample_merge, "Insrr",-3,7)
gene_plot(sample_merge, "Prss3",-3,7)
gene_plot(sample_merge, "Cpa1",-3,3)
gene_plot(sample_merge, "Ins1",-3,7)
gene_plot(sample_merge, "Col1a2",-3,7)
gene_plot(sample_merge, "Acta2",-2,9)
gene_plot(sample_merge, "S100a9",-3,9)
gene_plot(sample_merge, "Cd68",-3,9)
#---- Rename clusters ----------------------------------------------------------
sample_merge <- RenameIdents(sample_merge, 
                             "Acinar_C0" = 'Acinar', 
                             "Acinar_C1" = 'Acinar', 
                             "Acinar_C2" = 'Acinar', 
                             "Stroma_C3" = 'Stroma', 
                             "Acinar_C4" = 'Acinar', 
                             "Stroma_C5" = 'Stroma', 
                             "Islet_C6" = 'Islet', 
                             "Acinar_C7" = 'Acinar')
# Figure 4
gene_plot(sample_merge, "Il6",-2,8)
gene_plot(sample_merge, "Il10",-2,2)
gene_plot(sample_merge, "Ptgs2",-2,8)
gene_plot(sample_merge, "Hamp",-2,4)

# Others
gene_plot(sample_merge, "Rgs5",-3,7)
gene_plot(sample_merge, "Gcg",-3,8)
gene_plot(sample_merge, "Prss3",-3,7)
gene_plot(sample_merge, "Fgb",-2,8)
gene_plot(sample_merge, "Fgr",-2,8)
gene_plot(sample_merge, "Chil1",-2,8)
gene_plot(sample_merge, "Grem2",-2,8)

