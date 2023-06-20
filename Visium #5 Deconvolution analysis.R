# Deconvolution analysis

# setup------------------------------------------------------------------------
rm(list = ls()) ; gc() ; gc()
Sys.setenv("OPENBLAS_NUM_THREADS"=8)
easypackages::libraries(c("tidyverse","Seurat","patchwork","spacexr", "ggplot2","STdeconvolve",'ggpubr'))
wdir <- "set to your workdir"
setwd(wdir)
se <- readRDS('RDS/subclus_SPN.rds')

# RCTD------------------------------------------------------------------------
# STdata
vis_coords <- GetTissueCoordinates(se)
VisiumData <- SpatialRNA(coords = vis_coords, counts = GetAssayData(se,slot = "counts", assay = "Spatial"))
barcodes <- colnames(VisiumData@counts)
plot_puck_continuous(puck = VisiumData, barcodes = barcodes, plot_val = VisiumData@nUMI, 
                     size=1, ylimit=c(0,round(quantile(VisiumData@nUMI,0.9))), 
                     title='plot of nUMI') 

# reference data
SCdata <- readRDS("RDS/SCdata.rds")

## UMAP plot(reference data)
DimPlot(SCdata,group.by = 'cell_type')

SCreference <- Reference(counts = GetAssayData(SCdata,assay = "RNA", slot = "counts"),
                         cell_types = factor(SCdata$cell_type,levels = c('AdrenalCortex','Macrophage','Tcell')))

# run RCTD
myRCTD <- create.RCTD(VisiumData, SCreference, max_cores = 8)
myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")

# Create variables from the myRCTD object to plot results
barcodes <- colnames(myRCTD@spatialRNA@counts) # list of spatial barcodes
weights <- myRCTD@results$weights # Weights for each cell type per barcode

# Normalize per spot weights so cell type probabilities sum to 1 for each spot
norm_weights <- normalize_weights(weights) 
cell_type_names <- colnames(norm_weights) # List of cell types

# add deconvolution result to Visium meta data
se@meta.data <- merge(se@meta.data,norm_weights,by=0,all=T) %>% column_to_rownames(var="Row.names")
# saveRDS(se,'RDS/Deconvolution.rds')

# plot
se <- readRDS('RDS/Deconvolution.rds')
cell_type_names <- c('AdrenalCortex','Macrophage','Tcell')
SpatialFeaturePlot(se,features = cell_type_names, pt.size.factor = 3)
ggsave(filename = "Figure/Deconvolution_spatial.svg", height = 4, width = 12, dpi = 300)
SpatialFeaturePlot(se,features = 'Macrophage', pt.size.factor = 3)
ggsave(filename = "Figure/Deconvolution_Macrophage.svg", height = 4, width = 4, dpi = 300)
SpatialFeaturePlot(se,features = c('CD14','CD68','CD163'), ncol = 3, pt.size.factor = 3)
ggsave(filename = "Figure/Deconvolution_MacrophageMarkers.svg", height = 4, width = 12, dpi = 300)


col_pal <- c("#F8766D", "#619CFF", "#F564E3", "#00BFC4") ; names(col_pal) <- c('AdrenalCortex','SPN-F','SPN-R','CPA')
p <- lapply(cell_type_names,function(x){
  plot_df <- se@meta.data %>% dplyr::select(Proportion = all_of(x), group = annotation)

  lab <- wilcox.test(Proportion ~ group, data = plot_df) %>% .$p.value
  lab <- ifelse(lab < 0.05, '*','NS')

  plot_df %>% 
    ggplot(aes(x = group, y = Proportion)) +
    geom_violin(aes(fill = group), color = NA, alpha = 0.6) + 
    scale_fill_manual(values = col_pal) +
    geom_boxplot(width = 0.07, outlier.colour = NA, position = position_dodge(width = 0.9)) +
    ggsignif::geom_signif(comparisons = list(c("SPN-F", "SPN-R")), map_signif_level = F, annotations = lab) +
    theme_classic() + ggtitle(x) + ylab('Estimated Proportion') +
    scale_y_continuous(labels = scales::percent) +
    theme(plot.title = element_text(size = 20, face = 'bold'), 
          axis.title.x = element_blank(), axis.title.y = element_text(size = 15),
          axis.text.x = element_text(size = 20))
})
wrap_plots(p, nrow = 1) + guide_area() + plot_layout(guides = 'collect')
ggsave(filename = "Figure/Deconvolution.svg", height = 4, width = 12, dpi = 300)

