# Trajectory analysis 2
# # # # # # # # # # 
# setup--------------------------------------------------------
rm(list = ls()) ; gc() ; gc()

reticulate::use_condaenv(condaenv = "/Users/norifusa/miniforge3-intel/envs/R_reticulate_intel")
easypackages::libraries(c("tidyverse","Seurat","sctransform","glmGamPoi","patchwork",
                          "SeuratWrappers","ElPiGraph.R","igraph","destiny","phateR",
                          'slingshot','plotly','presto','ComplexHeatmap','circlize','RColorBrewer','scales','zoo',
                          "monocle3","CytoTRACE","reticulate"
))
wdir <- "set to your workdir"
setwd(wdir)
load("RDS/Trajectory.rda")


# PseudotimeDE----------------------------------------------------------------------------
library(PseudotimeDE)
for(group in c('SPN_lineage','CPA_lineage')){
  subse <- subset(se,Trajectory %in% c('Root',group) & annotation == 'SPN')
  pt_df <- subse@meta.data %>% rownames_to_column(var = 'cell') %>% dplyr::select(cell, pseudotime = Pseudotime) %>% as_tibble()
  norm_df <- GetAssayData(subse, assay = 'SCT', slot = 'data') %>% as.matrix()
  select_gene <- apply(norm_df,MARGIN = 1, FUN = var) %>% sort(decreasing = T) %>% 
    head(3000) %>% # narrow down to the top 3000 genes with the highest variance
    names()
  res <- runPseudotimeDE(gene.vec = select_gene, ori.tbl = pt_df,sub.tbl = NULL, 
                         mat = norm_df, model = "gaussian",mc.cores = 6)
  ## add correlation coefficient
  cor <- apply(t(norm_df[select_gene,]),MARGIN = 2,FUN = function(x){cor(x,subse$Pseudotime,method = 'pearson')})
  res$cor <- cor[res$gene]
  
  df <- res %>% filter(fix.pv < 0.05) %>% arrange(-cor)
  saveRDS(df,file = paste0('RDS/PT_DE_',group,'.rds'))
  p <- lapply(1:nrow(df), function(x){
    tmp <- df[x,]
    plotCurve(gene.vec = tmp$gene, ori.tbl = pt_df,mat = norm_df, model.fit = tmp$gam.fit) + ylab('Expression Level')
  }) %>% wrap_plots(ncol = 5)
  ggsave(p, filename = paste0('Figure/PT_DE_',group,'.pdf'), height = ceiling(nrow(df)/5)*3, width = 15, dpi = 300, limitsize = F)
}

DE_SPNlin <- readRDS('RDS/PT_DE_SPN_lineage.rds')
DE_CPAlin <- readRDS('RDS/PT_DE_CPA_lineage.rds')

# Heatmap(along pseudotime)---------------------------------------------------
se <- readRDS('RDS/subclus_SPN.rds')

for(lineage in c('SPN_lineage','CPA_lineage')){
tmp <- se %>% 
  subset(Trajectory %in% c("Root",lineage))
exp <- GetAssayData(tmp,slot="data",assay="SCT")

## order spots along pseudotime
cell_order <- tmp@meta.data %>% arrange(Pseudotime) %>% rownames()
heat.mat <- as.matrix(GetAssayData(tmp, slot = "data", assay = "SCT"))
## select genes
# keep.gene <- names(head(sort(apply(heat.mat,1,var),decreasing=T),n=100)) # 分散上位の遺伝子
keep.gene <- readRDS(paste0('RDS/PT_DE_',lineage,'.rds')) %>% .$gene %>% unique()
heat.mat <- heat.mat[keep.gene,cell_order]
## scale expression
heat.mat <- t(apply(heat.mat, MARGIN = 1, FUN = scale))
colnames(heat.mat) <- cell_order
## order genes(expression peaks in the order of Pseudotime)
gene.sort <- order(apply(t(rollapply(t(heat.mat), width=10, by=1, FUN=mean)), 1, which.max))
heat.mat <- heat.mat[gene.sort,]

## Gene clustering by kmeans
library(factoextra)
library(cluster)
set.seed(1234)
fviz_nbclust(heat.mat, kmeans, method = "wss")
fviz_nbclust(heat.mat, kmeans, method = "silhouette")

centers <- 2
clus_gene <- kmeans(heat.mat, centers = centers)$cluster
## save kmeans result 
# data.frame(gene = names(clus_gene), cluster = clus_gene) %>% arrange(cluster) %>% 
#   xlsx::write.xlsx(file = paste0('Result/',lineage,'_kmeans.xlsx'))
## load
clus_gene_df <- xlsx::read.xlsx(paste0('Result/',lineage,'_kmeans.xlsx'), sheetIndex = 1)
clus_gene <- clus_gene_df$cluster ; names(clus_gene) <- clus_gene_df$gene

clus_sort <- sapply(1:centers, function(cluster){
  gene <- names(clus_gene[clus_gene==cluster])
  mean_order <- apply(heat.mat[gene,],MARGIN = 1, FUN = which.max) %>% mean()
  names(mean_order) <- cluster
  return(mean_order)
}) %>% sort() %>% names()
row_split <- data.frame(gene = names(clus_gene), cluster = clus_gene, row.names = 'gene') %>% 
  mutate(cluster = factor(cluster, levels = clus_sort))
row_split <- row_split[rownames(heat.mat),]

## Heatmap
anno.df <- tmp@meta.data[cell_order,c("annotation","Pseudotime")]
col_fun <- colorRamp2(breaks = seq(from = -2, to = 2, length = 100), colors = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100))
# col_fun <- colorRamp2(breaks = seq(from = -2, to = 2, length = 100), colors = viridis::magma(n=100)) # viridis
col_pal <- c("#F8766D", "#619CFF", "#F564E3", "#00BFC4") ; names(col_pal) <- c('AdrenalCortex','SPN-F','SPN-R','CPA')
branch_col <- c("#7FC97F", "#BEAED4", "#FDC086") ; names(branch_col) <- c("Root","SPN_lineage",'CPA_lineage') # RColorBrewer::brewer.pal(name = 'Accent',n=3)
top_annotation <- HeatmapAnnotation(df = anno.df, show_legend = T,
                                    annotation_height = unit(.5, "inches"),
                                    annotation_legend_param = list(Pseudotime = list(title_position = "topcenter",
                                                                                     direction = "horizontal", 
                                                                                     legend_width = unit(2, "inches"))),
                                    col = list(annotation = col_pal,
                                               Pseudotime = colorRamp2(breaks = seq(0,max(tmp$Pseudotime),len = 100), 
                                                                       colors = viridis::viridis(n=100))),
                                    show_annotation_name = F)
ht <- Heatmap(matrix = heat.mat, col = col_fun, 
              row_split = row_split, cluster_row_slices = F,
              # cluster_rows = F,
              heatmap_width = unit(4, "inches"),
              heatmap_height = unit(7, "inches"),
              top_annotation = top_annotation,
              show_heatmap_legend = T, 
              use_raster = T,raster_by_magick = T,
              show_row_dend = F,
              cluster_columns = F, show_column_dend = F,
              show_column_names = F, show_row_names = F, 
              row_names_gp = grid::gpar(fontsize = 6), 
              heatmap_legend_param = list(title = "Scaled Expression", title_position = "topcenter",
                                          direction = "horizontal", legend_width = unit(2, "inches")))

draw(ht,merge_legends = TRUE, heatmap_legend_side = "bottom")
fig <- grid.grabExpr(draw(ht,merge_legends = TRUE, heatmap_legend_side = "bottom")) 
ggsave(fig, filename = paste0('Figure/Heatmap_',lineage,'.svg'), dpi = 300, width = 6, height = 8)
}





