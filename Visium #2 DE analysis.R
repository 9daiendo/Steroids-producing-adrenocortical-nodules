#2 DE analysis, Enrichment analysis

# setup--------------------------------------------------------------------
rm(list = ls()) ; gc() ; gc()

easypackages::libraries(c("Seurat","sctransform","glmGamPoi","patchwork","reshape2",
                          "org.Hs.eg.db","clusterProfiler","presto",'Hmisc',"tidyverse"))
wdir <- "set to your workdir"
setwd(wdir)
load(file = "RDS/PostDataPreparation.rda")
se <- postQC

## modify meta data label
df <- se@meta.data %>% 
  # mutate(HistologicalAnnotation2 = gsub(x = HistologicalAnnotation,pattern = 'CPCC1_N',replacement = 'SPN_N5')) %>% 
  # mutate(HistologicalAnnotation2 = gsub(x = HistologicalAnnotation2,pattern = 'CPCC2_N',replacement = 'SPN_N6')) %>% 
  mutate(HistologicalAnnotation = sapply(str_split(HistologicalAnnotation2,'[-]'),'[[',1)) %>% 
  mutate(HistologicalAnnotation = ifelse(HistologicalAnnotation%in%c('Capsule_A','Capsule_T'),'Capsule',HistologicalAnnotation))

se$HistologicalAnnotation <- factor(df$HistologicalAnnotation,levels=c('AdrenalCortex','CPA','SPN_N5','SPN_N6','Capsule'))
se$HistologicalAnnotation2 <- df$HistologicalAnnotation2
SpatialDimPlot(se,group.by = 'HistologicalAnnotation')
ggsave('Figure/HistologicalAnnotation.svg',dpi=300,height=6,width = 8)

# Clustering--------------------------------------------------------------------
resolution <- seq(from = 0.4, to = 1.2, by = 0.1)
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
val_plot
ggsave(val_plot,filename = "Figure/sc3_stability.tiff",dpi = 300, height = 6, width = 6)
print(paste0("Optimal resolution: ",optimal.resolution))

# set optimal resolution
se@misc$optimal.resolution <- optimal.resolution
se$seurat_clusters <- se@meta.data[,paste0("SCT_snn_res.",optimal.resolution)]
SpatialDimPlot(se,group.by = 'seurat_clusters')

# Comparison of histological findings and clustering results--------------------------------------------------------------------
df1 <- se@meta.data %>% 
  mutate(HistologicalAnnotation = as.character(HistologicalAnnotation)) %>% 
  mutate(HistologicalAnnotation = factor(ifelse(HistologicalAnnotation%in%c('SPN_N5','SPN_N6'),'SPN',HistologicalAnnotation),
                                         levels = c('AdrenalCortex','SPN','CPA','Capsule')))
df2 <- df1 %>% 
  dplyr::select(seurat_clusters) %>% group_by(seurat_clusters) %>% dplyr::summarise(count = n()) %>% 
  rename(seurat_cluster_count = count)
df1 %>% 
  dplyr::select(seurat_clusters, HistologicalAnnotation) %>% 
  group_by(seurat_clusters,HistologicalAnnotation) %>% dplyr::summarise(count = n()) %>%
  left_join(df2, by = 'seurat_clusters') %>% mutate(Fraction = count/seurat_cluster_count) %>% 
  ggplot(aes(x = HistologicalAnnotation, y = seurat_clusters)) +
  geom_point(aes(size = Fraction, color = HistologicalAnnotation)) +
  scale_size(name = 'fraction of spots', range = c(0,8)) +
  ylab('Unsupervised clustering') +
  theme_bw()
ggsave(height = 6, width = 8, dpi = 300, filename = "Figure/Compare_TranscriptvsHistological.svg")


group = 'HistologicalAnnotation'
Idents(se) <- se@meta.data[,group]
col_pal = scales::hue_pal()(n = length(levels(se))) ; names(col_pal) <- levels(se)
umap = DimPlot(se,label = T) + NoLegend()
spatial = wrap_plots(nrow = 3,
                     lapply(levels(se),function(x){
                       SpatialDimPlot(se, cells.highlight = CellsByIdentities(object = se, idents = x),
                                      facet.highlight = T, cols.highlight = c(col_pal[x],'lightgrey'))
                     }))
p = umap + spatial + plot_layout(design = c('ABBB'))
ggsave(plot = p, height = 6, width = 12, dpi = 300, filename = paste0("Figure/Umap_",group,'.tiff'))

group = 'seurat_clusters' ; Idents(se) <- se@meta.data[,group]
col_pal = scales::hue_pal()(n = length(levels(se))) ; names(col_pal) <- levels(se)
spatial = wrap_plots(nrow = 3,
                     lapply(levels(se),function(x){
                       tmp = se ; tmp$group = factor(ifelse(Idents(se)==x,x,'other'),levels = c(x,'other'))
                       Idents(tmp) = tmp$group 
                       cols = c(col_pal[x], '#F5F5F5') ; names(cols) = c(x,'other')
                       SpatialDimPlot(tmp,image.alpha = 0, cols = cols) + ggtitle(x) + 
                         theme(legend.position = 'none',plot.title = element_text(size = 20))
                     }))
ggsave(plot = umap, height = 6, width = 6, dpi = 300, filename = paste0("Figure/Umap_",group,'.tiff'))
ggsave(plot = spatial, height = 12, width = 20, dpi = 300, filename = paste0("Figure/Spatial_",group,'.tiff'))


# Annotation for clusters--------------------------------------------------------------------
new.ident <- factor(c(
  '0' = 'CPA',
  '1' = 'CPA',
  '2' = 'CPA',
  '3' = 'CPA',
  '4' = 'CPA',
  '5' = 'CPA',
  '6' = 'Capsule',
  '7' = 'CPA',
  '8' = 'AdrenalCortex',
  '9' = 'SPN',
  '10' = 'Capsule',
  '11' = 'CPA',
  '12' = 'SPN',
  '13' = 'CPA',
  '14' = 'CPA',
  '15' = 'CPA'
),
levels = c("AdrenalCortex","SPN","CPA","Capsule"))
se$annotation <- new.ident[se$seurat_clusters]
group = 'annotation' ; Idents(se) <- se@meta.data[,group]
col_pal = scales::hue_pal()(n = length(levels(se))) ; names(col_pal) <- levels(se)
spatial = wrap_plots(nrow = 2,
                     lapply(levels(se),function(x){
                       tmp = se ; tmp$group = factor(ifelse(Idents(se)==x,x,'other'),levels = c(x,'other'))
                       Idents(tmp) = tmp$group 
                       cols = c(col_pal[x], '#F5F5F5') ; names(cols) = c(x,'other')
                       SpatialDimPlot(tmp,image.alpha = 0, cols = cols) + ggtitle(x) + 
                         theme(legend.position = 'none',plot.title = element_text(size = 20))
                     }))
ggsave(spatial,filename = 'Figure/SpatialDimClustering.svg',dpi = 300, height = 12, width = 12)

umap = DimPlot(se)
ggsave(umap,filename = 'Figure/UmapClustering.svg',dpi = 300, height = 6, width = 8)

# save(se,file = "RDS/PostAnnotation.rda")
load("RDS/PostAnnotation.rda")

# DE analysis--------------------------------------------------------------------
for(group in c("seurat_clusters","HistologicalAnnotation",'annotation')){
  Idents(se) <- se@meta.data[,group]
  col_pal = scales::hue_pal()(n = length(levels(se))) ; names(col_pal) <- levels(se)
  # DEG_df <- FindAllMarkers(se, logfc.threshold = 0.5) %>% filter(p_val_adj < 0.01)
  DEG_df <- presto::wilcoxauc(se, group_by = group, seurat_assay = "SCT") %>% 
    filter(abs(logFC) > 0.25 & padj < 0.05) %>% arrange(group,-logFC)
  xlsx::write.xlsx(DEG_df, file = paste0("Result/DEG_", group, ".xlsx"), row.names = F)
  features <- DEG_df %>% group_by(group) %>% top_n(wt = auc, n = 10) %>% .$feature %>% unique()
  heatmap <- AverageExpression(object = se, group.by = group, features = features, return.seurat = F, assays = "SCT", slot = "data")[[1]] %>% 
    as.data.frame() %>% pheatmap::pheatmap(scale = "row", angle_col = 45, fontsize_col = 16,cluster_rows = T)
  ggsave(plot = heatmap, height = length(features)*0.2, width = length(unique(se@meta.data[,group])), dpi = 300,
         filename = paste0("Figure/Heatmap_",group,".tiff"))
}

## DEGのheatmap
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(pals)

Idents(se) <- se$annotation
col_fun <- colorRamp2(breaks = seq(from = -2, to = 2, length = 100), 
                      colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))

DEG_df <- presto::wilcoxauc(se, group_by = 'annotation', seurat_assay = "SCT") %>%
  filter(group != 'Capsule') %>% mutate(group = factor(group,levels=c('AdrenalCortex','SPN','CPA'))) %>% 
  filter(logFC > 0.25 & padj < 0.05) %>% 
  arrange(-logFC) %>% group_by(group) %>% slice(1:100) %>% ungroup()
DEG_df %>% 
  dplyr::select(group,feature,logFC,padj,auc) %>% arrange(group,-logFC) %>% 
  xlsx::write.xlsx(., file = 'Result/top100DEG.xlsx')
features <- DEG_df$feature

marker.gene <- xlsx::read.xlsx('Result/top100DEG.xlsx',sheetIndex = 1) %>% 
  drop_na(label) %>% .$feature

subse <- subset(se,annotation!='Capsule')

DotPlot(subse,features = marker.gene, assay = 'SCT', group.by = 'annotation',scale = T)+ 
  # scale_colour_gradientn(colours = rev(pals::brewer.spectral(n = 100)))+
  scale_colour_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)) +
  theme(axis.text.x = element_text(angle=45,hjust=1))

mat_df <- GetAssayData(subse,assay='SCT',slot='data') %>% as.data.frame() %>% .[unique(features),] %>% apply(., MARGIN = 1, FUN = scale) %>% t()
colnames(mat_df) <- Cells(subse)

ident <- ifelse(as.character(Idents(subse))=='AdrenalCortex','Adrenal\nCortex',as.character(Idents(subse))) %>% 
  factor(.,levels = c('Adrenal\nCortex','SPN','CPA'))
ident.colors <- (scales::hue_pal())(n = 4)[1:3] ; names(ident.colors) <- c('Adrenal\nCortex','SPN','CPA')

top_annotation <- HeatmapAnnotation(df = data.frame(row.names = "label",
                                                    label = colnames(subse),
                                                    ident = ident), 
                                    show_legend = F,
                                    col = list(ident = ident.colors),
                                    show_annotation_name = F)
right_annotation <- rowAnnotation(foo = anno_mark(at = sapply(marker.gene,FUN = function(x){charmatch(x,rownames(mat_df))}), labels = marker.gene))

ht <- Heatmap(matrix = mat_df, col = col_fun, 
              heatmap_width = unit(10, "in"), heatmap_height = unit(8, "in"),
              top_annotation = top_annotation,
              right_annotation = right_annotation,
              show_heatmap_legend = T, use_raster = T,
              column_split = ident, column_title_rot = 0,
              cluster_rows = F, cluster_columns = F, 
              show_column_names = F, show_row_names = F, 
              show_row_dend = F, show_column_dend = F,
              row_names_gp = grid::gpar(fontsize = 14), 
              heatmap_legend_param = list(title = "Scaled Expression", 
                                          # legend_height = unit(4, "cm"), title_position = "leftcenter-rot",
                                          legend_width = unit(2, "in"), title_position = "topcenter",
                                          direction = "horizontal"))
draw(ht,heatmap_legend_side = 'bottom')
fig <- grid.grabExpr(draw(ht,heatmap_legend_side = 'bottom')) 
ggsave(fig, filename = "Figure/Heatmap.svg", dpi = 300, width = 11, height = 9)

# Correlation between clusters--------------------------------------------------------------------
library(ComplexHeatmap)
library(circlize)
col_fun <- colorRamp2(breaks = seq(from = -1, to = 1, length = 100), 
                      colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))

norm_data = AverageExpression(subset(se,annotation!='Capsule'),assay = 'SCT', slot = 'scale.data', group.by = 'annotation')[[1]] %>% as.matrix()
colnames(norm_data) <- ifelse(colnames(norm_data)=='AdrenalCortex','Adrenal\nCortex',colnames(norm_data))
cor_data = rcorr(x = norm_data, type = 'pearson')

ht = Heatmap(cor_data$r, col = col_fun,heatmap_legend_param = list(title = "Pearson\ncorrelation"))
p <- grid.grabExpr(draw(ht)) 
ggsave(p, filename = "Figure/Corheatmap.tiff", height = 6, width = 7, dpi = 300)

# Enrichment analysis(metascape)--------------------------------------------------------------------
## Gene list
cluster.genes <- presto::wilcoxauc(se,group_by="annotation",seurat_assay="SCT")
df <- cluster.genes %>% 
  filter(abs(logFC) > 0.25 & padj < 0.05) %>% 
  mutate(UpDown = ifelse(logFC>0,'Up','Down')) %>% 
  mutate(type = paste0(group,"_",UpDown)) %>% 
  group_by(type) %>% dplyr::select(feature) %>% group_split(.keep = T)
df %>% lapply(.,function(x){unlist(length(x$feature))})

label <- df %>% lapply(.,function(x){unlist(unique(as.vector(x$type)))}) %>% unlist
gene_df <- df %>% lapply(.,function(x){unlist(as.vector(x$feature))}) %>% plyr::ldply(., rbind) %>% t()
colnames(gene_df) <- label
gene_df %>% write.csv('gene_df.csv',row.names = F) # metascape用
gene_list <- df %>% lapply(.,function(x){unlist(paste0(x$feature,collapse = ', '))})
names(gene_list) <- label

# Metascape結果のplot
df <- lapply(c('AdrenalCortex','SPN','CPA'),function(x){
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
  # filter(Summary == 'Summary') %>%
  filter(GroupID <= 10) %>% 
  group_by(GroupID,type) %>% arrange(LogP) %>% 
  slice(1:1) %>%
  .$Term
select_df <- df %>% filter(Term %in% select_pathway) %>% group_by(type) %>% distinct(Term, .keep_all = T) %>% ungroup() %>% 
  separate(InTerm_InList,c("count","total"),sep="/") %>% 
  mutate(GeneRatio = as.double(count)/as.double(total)) %>%
  mutate(Description = paste0(Category,': ',Description)) %>% # pathwayの頭にCategoryを追加する
  mutate(Description = stringr::str_wrap(Description,width = 50)) %>% 
  mutate(type = factor(type,levels=c('AdrenalCortex','SPN','CPA'))) %>% 
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
  scale_x_discrete(limits = c('AdrenalCortex','SPN','CPA'), label = c('AdrenalCortex'='Adrenal\nCortex','SPN'='SPN','CPA'='CPA')) +
  scale_color_viridis_c(option = 'D') + theme_minimal() + 
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 12, face = 'bold'),
        axis.text.y = element_text(size = 10))

ggsave('Figure/GOplot.svg',dpi=300, height = 9, width = 8)





# Enrichment analysis(clusterprofiler)------------------------------------------------------------ 
## GO Over-Representation Analysis
for(group in c("seurat_clusters","HistologicalAnnotation",'annotation')){
  gene_df <- read.csv(paste0("Result/DEG_",group,".csv"),row.names = 1)
  entrez_df <- AnnotationDbi::select(x = org.Hs.eg.db, keys = unique(gene_df$feature), columns = c("SYMBOL","ENTREZID"), keytype = "SYMBOL")
  df <- gene_df %>% left_join(entrez_df, by = c("feature"="SYMBOL")) %>% 
    mutate(FC = factor(ifelse(logFC>0,"upregulated","downregulated"),levels=c("upregulated","downregulated")))
  formula_res <- clusterProfiler::compareCluster(geneClusters = ENTREZID ~ group + FC, 
                                                 data = df, fun = "enrichGO", OrgDb = org.Hs.eg.db,
                                                 ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.05, readable = TRUE)
  saveRDS(formula_res,file = paste0("RDS/formula_res_",group,".rds"))
}

for(group in c("seurat_clusters","HistologicalAnnotation",'annotation')){
  formula_res <- readRDS(file = paste0("RDS/formula_res_",group,".rds"))
  p <- clusterProfiler::dotplot(formula_res, x = "FC", showCategory = 5) + facet_grid(~group) +
    theme_bw(base_size = 12) + theme(axis.title.x = element_blank())
  ggsave(p, height = 16,width = 20, dpi = 300, limitsize = FALSE,
         filename = paste0("Figure/GOanalysis_",group,".tiff"))
}

## GSEA
human.genes <- msigdbr::msigdbr(species = "Homo sapiens")
pathways.GOBP <- human.genes %>% 
  filter(gs_subcat == "CP:KEGG") %>% 
  split(x = .$gene_symbol, f = .$gs_name)
cluster.genes <- presto::wilcoxauc(se,group_by="annotation",seurat_assay="SCT")
gsea_res <- list(0,0,0,0) ; names(gsea_res) <- levels(Idents(se))

for(i in levels(Idents(se))){
  ranks <- cluster.genes %>% filter(group == i) %>%
    arrange(desc(auc),desc(logFC)) %>% 
    dplyr::select(feature, auc) %>% 
    deframe()
  fgseaRes <- fgsea::fgsea(pathways.GOBP, stats = ranks, nPermSimple = 1000)
  fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES)) %>% 
    mutate(genelist = sapply(leadingEdge, function(x){head(x, n = 3) %>% paste0(collapse = ", ")}))
  
  p <- fgseaResTidy %>% filter(padj < 0.05) %>% group_by(NES>0) %>% 
    top_n(wt = abs(NES), n = 10) %>% 
    ggplot(aes(x = NES, y = reorder(pathway, NES))) + geom_col(aes(fill=-padj)) +
    geom_text(aes(label = genelist, hjust = ifelse(NES>0,0,1), x = ifelse(NES>0,0.25,-0.25)), color = "white", fontface = "bold") + 
    labs(y = "GO BP Pathways", x = "Normalized Enrichment Score") + theme_bw() +
    scale_fill_gradientn(colors = c("blue","red"), name = "adjusted p value")
  ggsave(p, filename = paste0("Figure/GSEA_",i,".tiff"),height = 12, width = 16, dpi=300)
  gsea_res[[i]] <- fgseaResTidy %>% mutate(group = i)
}
merged.res <- do.call(rbind,gsea_res)
select.pathway <- merged.res %>% filter(padj < 0.05) %>% 
  group_by(group,NES>0) %>% top_n(n = 5, wt = abs(NES)) %>% .$pathway
merged.res %>% 
  filter(pathway %in% select.pathway) %>% 
  dplyr::select(pathway,NES,group) %>% 
  pivot_wider(values_from = 'NES',names_from = 'group') %>% 
  # filter(pathway %in% c(
  #   'GOBP_CELL_CELL_SIGNALING_BY_WNT',
  #   'GOBP_RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA',
  #   'GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY',
  #   'GOBP_REGULATION_OF_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY',
  #   'GOBP_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY',
  #   'GOBP_POSITIVE_REGULATION_OF_CELL_DEATH',
  #   'GOBP_RESPONSE_TO_INTERFERON_ALPHA',
  #   'GOBP_REGULATION_OF_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_VIA_DEATH_DOMAIN_RECEPTORS',
  #   'GOBP_REGULATION_OF_CANONICAL_WNT_SIGNALING_PATHWAY',
  #   'GOBP_REGULATION_OF_CELLULAR_RESPONSE_TO_STRESS',
  #   'GOBP_REGULATION_OF_AUTOPHAGY'
  # )) %>%
  column_to_rownames(var = 'pathway') %>% 
  drop_na() %>%
  pheatmap::pheatmap(scale = 'none')

interest.pathway <- "KEGG_P53_SIGNALING_PATHWAY"
merged.res %>% 
  filter(pathway==interest.pathway) %>%
  ggplot(aes(x = NES, y = pathway, fill = group)) +geom_bar(stat = "identity", position = "dodge") 

