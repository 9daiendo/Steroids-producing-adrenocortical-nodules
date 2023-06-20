# Subclustering of SPN
# # # # # # # # # # 
# setup--------------------------------------------------------
rm(list = ls()) ; gc() ; gc()

easypackages::libraries(c("tidyverse","Seurat","sctransform","glmGamPoi","patchwork",
                          "SeuratWrappers","ElPiGraph.R","igraph",'ggrepel','presto',
                          'hdWGCNA','WGCNA','ComplexHeatmap','circlize','RColorBrewer','scales','zoo'
))
wdir <- "set to your workdir"
setwd(wdir)
load("RDS/Trajectory.rda")

# subclustering--------------------------------------------------------------------
se <- subset(se, annotation == 'SPN')
se@meta.data <- se@meta.data[,grep('SCT_snn_res.',colnames(se@meta.data),invert = T)]
resolution <- seq(from = 0.4, to = 2, by = 0.1)
set.seed(1234)
se <- se %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = resolution)

# Decision using SC3.stability method
prefix <- "SCT_snn_res."
clt <- clustree::clustree(se, prefix = prefix, node_colour = "sc3_stability")
clt.plot <- clt$data %>% dplyr::select(!!sym(prefix), sc3_stability) %>% 
  dplyr::group_by(!!sym(prefix)) %>% dplyr::summarise(mean(sc3_stability)) %>% as.data.frame()
names(clt.plot)[2] <- "sc3_stability"
clt.plot[,prefix]<-as.numeric(as.character(clt.plot[,prefix]))
optimal.resolution <- clt.plot %>% arrange(desc(sc3_stability)) %>% head(1) %>% pull(!!sym(prefix))
val_plot <- ggplot(clt.plot, mapping=aes(x=!!sym(prefix), y=sc3_stability)) +
  geom_point() + theme_classic() + geom_line() + xlab("resolution") + 
  geom_vline(xintercept = optimal.resolution,color="red", size=2) + 
  scale_x_continuous(breaks = resolution)

ggsave(val_plot,filename = "Figure/sc3_stability_SPN.tiff",dpi = 300, height = 6, width = 6)
print(paste0("Optimal resolution: ",optimal.resolution))

# set optimal resolution
se@misc$optimal.resolution <- optimal.resolution
se$seurat_clusters <- se@meta.data[,paste0("SCT_snn_res.",optimal.resolution)]
SpatialDimPlot(se,group.by = 'seurat_clusters')

group = 'seurat_clusters' ; Idents(se) <- se@meta.data[,group]
col_pal <- c("#F8766D", "#619CFF", "#F564E3", "#00BFC4") ; names(col_pal) <- c('AdrenalCortex','SPN-F','SPN-R','CPA')
spatial <- wrap_plots(nrow = 2,
                     lapply(levels(se),function(x){
                       tmp = se ; tmp$group = factor(ifelse(Idents(se)==x,x,'other'),levels = c(x,'other'))
                       Idents(tmp) = tmp$group 
                       cols = c(col_pal[x], '#F5F5F5') ; names(cols) = c(x,'other')
                       SpatialDimPlot(tmp,image.alpha = 0, cols = cols, pt.size.factor = 2.4) + ggtitle(x) + 
                         theme(legend.position = 'none',plot.title = element_text(size = 20))
                     }))
ggsave(plot = spatial, height = 8, width = 12, dpi = 300, filename = paste0("Figure/Spatial_",group,'_SPN.svg'))

# referenceのplot(Histological annotation)
se$HistologicalAnnotation2 <- gsub(se$HistologicalAnnotation2, pattern = '_N5',replacement = '')
se$HistologicalAnnotation2 <- gsub(se$HistologicalAnnotation2, pattern = '_N6',replacement = '')
se$HistologicalAnnotation2 <- gsub(se$HistologicalAnnotation2, pattern = 'CPA-F',replacement = 'CPA')
SpatialDimPlot(se,group.by = 'HistologicalAnnotation2')
group = 'HistologicalAnnotation2' ; Idents(se) <- se@meta.data[,group]
col_pal <- c("#F8766D", "#619CFF", "#F564E3", "#00BFC4") ; names(col_pal) <- c('AdrenalCortex','SPN-F','SPN-R','CPA')
spatial <- lapply(levels(se),function(x){
                  tmp = se ; tmp$group = factor(ifelse(Idents(se)==x,x,'other'),levels = c(x,'other'))
                  Idents(tmp) = tmp$group 
                  cols = c(col_pal[x], '#F5F5F5') ; names(cols) = c(x,'other')
                  SpatialDimPlot(tmp,image.alpha = 0, cols = cols, pt.size.factor = 2.4) + ggtitle(x) + 
                  theme(legend.position = 'none',plot.title = element_text(size = 20))
                  })[c(2,3)] %>% 
  wrap_plots(nrow = 1)

ggsave(plot = spatial, height = 6, width = 12, dpi = 300, filename = paste0("Figure/Spatial_",group,'_SPN.svg'))


# Comparison of histological findings and clustering results
df1 <- se@meta.data %>% 
  mutate(HistologicalAnnotation = as.character(HistologicalAnnotation2)) %>% 
  mutate(HistologicalAnnotation = gsub(x = HistologicalAnnotation,pattern = '_N5',replacement = '')) %>% 
  mutate(HistologicalAnnotation = gsub(x = HistologicalAnnotation,pattern = '_N6',replacement = '')) %>% 
  mutate(HistologicalAnnotation = gsub(x = HistologicalAnnotation,pattern = 'CPA-F',replacement = 'CPA')) %>% 
  mutate(HistologicalAnnotation = factor(HistologicalAnnotation,levels = c('AdrenalCortex','SPN-F','SPN-R','CPA')))
df2 <- df1 %>% 
  dplyr::select(seurat_clusters) %>% group_by(seurat_clusters) %>% dplyr::summarise(count = n()) %>% 
  rename(seurat_cluster_count = count)
col_pal <- c("#F8766D", "#619CFF", "#F564E3", "#00BFC4") ; names(col_pal) <- c('AdrenalCortex','SPN-F','SPN-R','CPA')
df1 %>% 
  dplyr::select(seurat_clusters, HistologicalAnnotation) %>% 
  group_by(seurat_clusters,HistologicalAnnotation) %>% dplyr::summarise(count = n()) %>%
  left_join(df2, by = 'seurat_clusters') %>% mutate(Fraction = count/seurat_cluster_count) %>% 
  ggplot(aes(x = HistologicalAnnotation, y = seurat_clusters)) +
  geom_point(aes(size = Fraction, color = HistologicalAnnotation)) +
  scale_size(name = 'fraction of spots', range = c(0,8)) +
  ylab('Unsupervised clustering') +
  scale_color_manual(values = col_pal) +
  theme_bw()
ggsave(height = 6, width = 8, dpi = 300, filename = "Figure/Compare_TranscriptvsHistological_SPN.svg")

# annotation--------------------------------------------------------------------
new.ident <- factor(c(
  '0' = 'SPN-R',
  '1' = 'SPN-F',
  '2' = 'SPN-F',
  '3' = 'SPN-R',
  '4' = 'SPN-F',
  '5' = 'SPN-R'
),
levels = c("SPN-F","SPN-R"))
se$annotation <- new.ident[se$seurat_clusters]
group = 'annotation' ; Idents(se) <- se@meta.data[,group]
col_pal <- c("#F8766D", "#619CFF", "#F564E3", "#00BFC4") ; names(col_pal) <- c('AdrenalCortex','SPN-F','SPN-R','CPA')
spatial = wrap_plots(nrow = 1,
                     lapply(levels(se),function(x){
                       tmp = se ; tmp$group = factor(ifelse(Idents(se)==x,x,'other'),levels = c(x,'other'))
                       Idents(tmp) = tmp$group 
                       cols = c(col_pal[x], '#F5F5F5') ; names(cols) = c(x,'other')
                       SpatialDimPlot(tmp,image.alpha = 0, cols = cols, pt.size.factor = 2.4) + ggtitle(x) + 
                         theme(legend.position = 'none',plot.title = element_text(size = 20))
                     }))
ggsave(spatial,filename = 'Figure/SpatialDimClustering_SPN.svg',dpi = 300, height = 6, width = 12)

# saveRDS(se,'RDS/subclus_SPN.rds')
se <- readRDS('RDS/subclus_SPN.rds')

# DE analysis--------------------------------------------------------------------
DEG_df <- presto::wilcoxauc(se, group_by = 'annotation', seurat_assay = "SCT") %>%
  filter(group == 'SPN-F') %>% arrange(-logFC)
# xlsx::write.xlsx(DEG_df, file = 'Result/DEG_SPN.xlsx', row.names = F) # labelする遺伝子に◯をつける

## Volcano plot
toptable <- xlsx::read.xlsx(file = 'Result/DEG_SPN.xlsx',sheetIndex = 1)
up <- toptable %>% filter(logFC > 0.25 & padj < 0.05) %>% nrow() ; down <- toptable %>% filter(logFC < -0.25 & padj < 0.05) %>% nrow()
highlight_df <- toptable %>% filter(label == '◯')

toptable %>% 
  mutate(col = ifelse(logFC > 0.25 & padj < 0.05, "SPN-F", ifelse(logFC < -0.25 & padj < 0.05, "SPN-R", "n.p."))) %>%
  ggplot(aes(x = logFC, y = -log10(padj))) + geom_point(aes(color = col), alpha = 0.6) +
  scale_color_manual(values = c("#619CFF","#F564E3","lightgrey")) + 
  ggrepel::geom_label_repel(data = highlight_df,aes(label = feature), box.padding = unit(0.5, "lines"), max.overlaps = 10) +
  geom_vline(xintercept = 0.25, size = 0.2, linetype = "dashed") +
  geom_vline(xintercept = -0.25, size = 0.2, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), size = 0.2, linetype = "dashed") +
  geom_text(aes(x = max(abs(logFC)*1.1), y = 2), label = paste0("SPN-F (", up,")"), hjust = 1, vjust = 0.5) +
  geom_text(aes(x = -max(abs(logFC)*1.1), y = 2), label = paste0("SPN-R (",down,")"), hjust = 0, vjust = 0.5) +
  theme_classic() + NoLegend() + theme(axis.title = element_text(size = 15), plot.title = element_text(face = "bold", size = 20)) +
  xlab("Log fold change") + ylab(bquote(~-Log[10]~"P")) +
  xlim(c(-max(abs(toptable$logFC)*1.1),max(abs(toptable$logFC))*1.1)) +
  ylim(c(0,max(-log10(toptable$padj))*1.1))

ggsave(filename = paste0("Figure/Volcano_SPN.svg"), height = 6, width = 6, dpi=300)

## Violin plot of steroidogenic enzyme genes
steroidogenic_gene <- c("CYP11A1","CYP11B1","CYP11B2","CYP21A2","CYP17A1","CYB5A","SULT2A1")

plot_df <- GetAssayData(se,assay = 'SCT', slot = 'data') %>% as.data.frame() %>% .[steroidogenic_gene,] %>% t() %>% as.data.frame()
all.equal(rownames(plot_df),Cells(se))
plot_df <- cbind(plot_df,se@meta.data)
p <- lapply(steroidogenic_gene,function(x){
  lab <- toptable %>% filter(feature == x) %>% 
    mutate(lab = ifelse(abs(logFC) > 0.25 & padj < 0.05, '*','NS')) %>% .$lab
  plot_df %>% 
    dplyr::select(Expression = all_of(x), everything()) %>% 
    ggplot(aes(x = annotation, y = Expression)) +
    geom_violin(aes(fill = annotation), color = NA, alpha = 0.6) + 
    scale_fill_manual(values = col_pal, name = 'SPN layer') +
    geom_boxplot(width = 0.07, outlier.colour = NA, position = position_dodge(width = 0.9)) +
    ggsignif::geom_signif(comparisons = list(c("SPN-F", "SPN-R")), map_signif_level = F, annotations = lab) +
    theme_classic() + ggtitle(x) + ylab('Expression Level') +
    theme(plot.title = element_text(size = 20, face = 'bold'), 
          legend.title = element_text(size = 15),legend.text = element_text(size = 15),
          axis.title.x = element_blank(), axis.title.y = element_text(size = 15),
          axis.text.x = element_text(size = 20))
  })
wrap_plots(p, nrow = 2) + guide_area() + plot_layout(guides = 'collect')
ggsave(filename = "Figure/SPN_steroidogenic_gene.svg", height = 8, width = 16, dpi = 300)

p <- lapply(steroidogenic_gene,function(x){
  SpatialFeaturePlot(se,features = x, pt.size.factor = 3) + 
    ggtitle(x) +
    scale_fill_gradientn(colours = rev(pals::brewer.spectral(n = 100)),
                         labels = c('low','high'), name = 'Scaled\nexpression',
                         breaks = unlist(lapply(list(min,max), function(f) f(GetAssayData(se[x,],slot='data',assay='SCT'))))) +
    theme(legend.direction = 'vertical', 
          legend.title = element_text(size = 15),legend.text = element_text(size = 15),
          plot.title = element_text(size = 20, face = 'bold')) 
})
wrap_plots(p, nrow = 2) + guide_area() + plot_layout(guides = 'collect')
ggsave('Figure/SPN_steroidogenic_gene_spatial.svg',height = 8, width = 16, dpi = 300)

## Enrichment analysis
df <- lapply(c('SPN-F','SPN-R'),function(x){
  xlsx::read.xlsx(paste0('Result/Metascape/',x,'/metascape_result.xlsx'),sheetName = 'Enrichment') %>%
    mutate(Summary = sapply(strsplit(GroupID,'_'),'[[',2)) %>% 
    mutate(GroupID = as.integer(sapply(strsplit(GroupID,'_'),'[[',1))) %>% 
    mutate(Category = gsub(x = Category, pattern = 'Canonical Pathways', replacement = 'CP'),
           Category = gsub(x = Category, pattern = 'GO Biological Processes', replacement = 'GOBP'),
           Category = gsub(x = Category, pattern = 'KEGG Pathway', replacement = 'KEGG'),
           Category = gsub(x = Category, pattern = 'Reactome Gene Sets', replacement = 'REACTOME'),
           Category = gsub(x = Category, pattern = 'WikiPathways', replacement = 'WIKI')) %>% 
    mutate(type = x)
}) %>% do.call(rbind,.)
select_pathway <- df %>% 
  filter(Category %in% c('GOBP','KEGG')) %>% # GOBPに絞る場合
  filter(Summary != 'Summary') %>%
  group_by(type,GroupID) %>% arrange(LogP) %>% slice(1:1) %>%
  group_by(type) %>% top_n(n = 10, wt = -GroupID) %>% # 上位10groupのtermにしぼる
  .$Term
select_df <- df %>% filter(Term %in% select_pathway) %>% group_by(type) %>% distinct(Term, .keep_all = T) %>% ungroup() %>% 
  separate(InTerm_InList,c("count","total"),sep="/") %>% 
  mutate(GeneRatio = as.double(count)/as.double(total)) %>%
  mutate(Description = paste0(Category,': ',Description)) %>% # pathwayの頭にCategoryを追加する
  mutate(Description = stringr::str_wrap(Description,width = 50)) %>% 
  mutate(type = factor(type,levels=c('SPN-F','SPN-R'))) %>% 
  arrange(type)
orderVec <- NULL
for(i in 1:length(select_df$Description)){
  if(select_df$Description[i] %in% orderVec){next} else {
    orderVec <- c(orderVec,select_df$Description[i])}}
select_df <- select_df %>% transform(Description = factor(Description, levels = orderVec))

select_df %>% 
  ggplot(aes(x = type, y = Description)) + 
  geom_point(aes(color = -LogP, size = GeneRatio)) +
  scale_size_continuous(range=c(5, 10)) +
  # scale_x_discrete(limits = c('SPN-F','SPN-R'), label = c('AdrenalCortex_Up'='Adrenal\nCortex','CPA_Up'='CPA','CPCC_Up'='SPN')) +
  scale_color_viridis_c(option = 'D') + theme_minimal() + 
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 14, face = 'bold'),
        axis.text.y = element_text(size = 10))

ggsave('Figure/SPN_GOplot.svg',dpi=300, height = 6, width = 7)


# GSEA analysis--------------------------------------------------
library(fgsea)
library(msigdbr)
human.genes <- msigdbr(species = "Homo sapiens")
table(human.genes$gs_subcat)
genesets.interest <- human.genes %>% filter(gs_subcat %in% c("GO:BP",'CP:KEGG')) 
pathways.interest <- genesets.interest %>% split(x = .$gene_symbol, f = .$gs_name)
names(pathways.interest) %>% length()

se <- readRDS('RDS/subclus_SPN.rds')
toptable <- wilcoxauc(se,group_by = "annotation", seurat_assay = "SCT") %>%
  filter(group=="SPN-F") %>% arrange(-logFC)
stats <- toptable %>% 
  # filter(abs(logFC) > 0.25 & padj < 0.05) %>%
  arrange(desc(logFC)) %>% dplyr::select(feature, logFC) %>% deframe()

gsea_res <- fgsea(pathways = pathways.interest, stats = stats) 

gsea_res %>% 
  filter(pval < 0.05 & padj < 0.05) %>%
  arrange(-NES) %>% view()



# pathway activity inference(PROGENy)---------------------------------------
se <- readRDS('RDS/subclus_SPN.rds')
library(decoupleR)
net <- get_progeny(organism = 'human', top = 100)
data <- se 
mat <- GetAssayData(data, assay = 'SCT', slot = 'data') %>% as.matrix()
acts <- run_wmean(mat=mat, net=net, .source='source', .target='target', .mor='weight', times = 100, minsize = 5)
acts
# Extract norm_wmean and store it in pathwayswmean in data
data[['pathwayswmean']] <- acts %>%
  filter(statistic == 'norm_wmean') %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)
# Change assay
DefaultAssay(object = data) <- "pathwayswmean"
# Scale the data
data <- ScaleData(data)
data@assays$pathwayswmean@data <- data@assays$pathwayswmean@scale.data
p1 <- DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend() + ggtitle('Cell types')
p2 <- FeaturePlot(data, features = c("WNT"), min.cutoff = 0) & 
  scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')
p1 | p2
wilcoxauc(data,seurat_assay = 'pathwayswmean') %>% filter(group == 'SPN-F' & padj < 0.05)
SpatialFeaturePlot(data, rownames(data), pt.size.factor = 3)

# Extract activities from object as a long dataframe
df <- t(as.matrix(data@assays$pathwayswmean@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(data)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))
# Transform to wide matrix
top_acts_mat <- df %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  t() %>% 
  as.matrix()

library(ComplexHeatmap)
library(circlize)
col_fun <- colorRamp2(breaks = seq(from = -2, to = 2, length = 100), 
                      colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))
col_fun <- colorRamp2(breaks = seq(from = -1, to = 1, length = 100), 
                      colors = colorRampPalette(c("Darkblue", "white","red"))(100))
ht <- Heatmap(top_acts_mat, col = col_fun,column_names_rot = 0, column_names_centered = T,
              heatmap_legend_param = list(title = "Pathway\nactivity"),cluster_columns = F)
ht
p <- grid.grabExpr(draw(ht)) 
ggsave(p, filename = "Figure/PROGENy.svg", height = 6, width = 7, dpi = 300)

## DoRotheaによる転写因子推定
net <- get_dorothea(organism='human', levels=c('A', 'B', 'C'))
net$source %>% table()
acts <- run_wmean(mat=mat, net=net, .source='source', .target='target',.mor='mor', times = 100, minsize = 5)
acts
# Extract norm_wmean and store it in tfswmean in pbmc
data[['tfswmean']] <- acts %>%
  filter(statistic == 'norm_wmean') %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)
# Change assay
DefaultAssay(object = data) <- "tfswmean"
# Scale the data
data <- ScaleData(data)
data@assays$tfswmean@data <- data@assays$tfswmean@scale.data
SpatialFeaturePlot(data,'NR5A1')

# Extract activities from object as a long dataframe
df <- t(as.matrix(data@assays$tfswmean@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(data)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Get top tfs with more variable means across clusters
tfs <- wilcoxauc(data,seurat_assay = 'tfswmean') %>% 
  filter(group == 'SPN-F' & padj < 0.05) %>% 
  pull(feature)

# Subset long data frame to top tfs and transform to wide matrix
top_acts_mat <- df %>%
  filter(source %in% tfs) %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>% t() %>% 
  as.matrix()

col_fun <- colorRamp2(breaks = seq(from = -1, to = 1, length = 100), 
                      colors = colorRampPalette(c("Darkblue", "white","red"))(100))
ht <- Heatmap(top_acts_mat, col = col_fun,column_names_rot = 0, column_names_centered = T,
              heatmap_legend_param = list(title = "TF activity score"),cluster_columns = F)
ht
p <- grid.grabExpr(draw(ht)) 
ggsave(p, filename = "Figure/DoRothea.svg", height = 8, width = 7, dpi = 300)

toptable <- wilcoxauc(data,seurat_assay = 'tfswmean') %>% filter(group == 'SPN-F') 
up <- toptable %>% filter(logFC > 0 & padj < 0.05) %>% nrow() ; down <- toptable %>% filter(logFC < 0 & padj < 0.05) %>% nrow()
toptable %>% mutate(col = ifelse(logFC > 0 & padj < 0.05, "SPN-F", ifelse(logFC < 0 & padj < 0.05, "SPN-R", "n.p."))) %>%
  ggplot(aes(x = logFC, y = -log10(padj))) + geom_point(aes(color = col), alpha = 0.6) +
  scale_color_manual(values = c("#619CFF","#F564E3","lightgrey")) + 
  geom_vline(xintercept = 0, size = 0.2, linetype = "dashed") +
  # geom_vline(xintercept = -0.25, size = 0.2, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), size = 0.2, linetype = "dashed") +
  geom_text(aes(x = max(abs(logFC)), y = 1), label = paste0("SPN-F (", up,")"), hjust = 1, vjust = 0.5) +
  geom_text(aes(x = -max(abs(logFC)), y = 1), label = paste0("SPN-R (",down,")"), hjust = 0, vjust = 0.5) +
  theme_classic() + NoLegend() + theme(axis.title = element_text(size = 15), plot.title = element_text(face = "bold", size = 20)) +
  xlab("TF activity score") + ylab(bquote(~-Log[10]~"P")) +
  xlim(c(-max(abs(toptable$logFC)),max(abs(toptable$logFC)))) +
  geom_label_repel(data = toptable %>% filter(padj < 0.05) %>% group_by(logFC>0) %>% top_n(n = 10, wt = abs(logFC)),
                   aes(label = feature), box.padding = unit(0.5, "lines"), max.overlaps = 10) 
ggsave(filename = "Figure/DoRotheaVolcano.svg", height = 6, width = 6, dpi = 300)

