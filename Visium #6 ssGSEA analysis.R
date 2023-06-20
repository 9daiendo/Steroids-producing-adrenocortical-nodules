# ssGSEA

# setup--------------------------------------------------------
rm(list = ls()) ; gc() ; gc()
easypackages::libraries(c("tidyverse","Seurat","sctransform","glmGamPoi","patchwork",
                          "SeuratWrappers", "msigdbr", "rdist","GSVA",'escape'))
wdir <- "set to your workdir"
setwd(wdir)
se <- readRDS('RDS/Deconvolution.rds')

# ssGSEA--------------------------------------------------------
set.seed(1234)   

# prepare gene sets
## genesets
human.genes <- msigdbr(species = "Homo sapiens")
genesets.interest <- human.genes %>% filter(gs_cat %in% c('H')) 
pathways.interest <- genesets.interest %>% split(x = .$gene_symbol, f = .$gs_name)
gs <- pathways.interest
# gs <- pathways.interest[grep("_SENESCENCE", names(pathways.interest), value = T)]

## curated gene sets (DOI: 10.1126/sciadv.add0422)
gs <- c(gs,readRDS("RDS/sci_adv_gs.rds"))

# run ssGSEA 
mat <- GetAssayData(se,slot = "counts", assay = "Spatial") # raw countを入力とする
ES <- enrichIt(obj = mat, gene.sets = gs, ssGSEA.norm = T, cores = 6,  min.size = 1)
se <- AddMetaData(se, ES)
# saveRDS(se,"RDS/ssGSEA.rds")

res_df <- data.frame(Cor = apply(ES,MARGIN = 2, FUN = function(x){cor.test(x,se$Macrophage)$estimate}),
           P = apply(ES,MARGIN = 2, FUN = function(x){cor.test(x,se$Macrophage)$p.value})) %>% 
  rownames_to_column(var = 'Pathway') %>% 
  left_join(by = 'Pathway',data.frame(Gene = sapply(gs,FUN = function(x){paste0(x,collapse = ',')}), Pathway = names(gs)))

select_pathway <- res_df %>% 
  filter(abs(Cor) > 0.3 & P < 0.05) %>% 
  filter(Pathway %in% grep('HALLMARK',.$Pathway,value=T)) %>% ## HALLMARK geneset
  filter(Pathway %in% c('SENESCENCE_EGGERT','KUILMAN_OIS','SENESCENCE_OZCAN','ACOSTA_ET_AL_SENESCENCE','FRIDMAN_TAINSKY_SENESCENCE','SASP_COPPE','SenMayo','SENESCENCE_BUHL')) %>% # Senescence
  arrange(-Cor) %>% .$Pathway

select_pathway <- c('HALLMARK_IL6_JAK_STAT3_SIGNALING','HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_APOPTOSIS','SenMayo')

col_pal <- c("#F8766D", "#619CFF", "#F564E3", "#00BFC4") ; names(col_pal) <- c('AdrenalCortex','SPN-F','SPN-R','CPA')
lapply(select_pathway,function(pathway){
  cor <- paste('r =',res_df %>% filter(Pathway == pathway) %>% .$Cor %>% round(2))
  p <- res_df %>% filter(Pathway==pathway) %>% .$P ; p <- ifelse(p < 0.01,'P < 0.01', paste('P =',round(p,2)))
  se@meta.data %>% 
    ggplot(aes(x = Macrophage, y = !!sym(pathway))) +
    geom_point(aes(color = annotation), alpha = 0.7, size = 2) + 
    scale_color_manual(name = 'SPN layer',values = col_pal) +
    geom_smooth(method = 'lm', se = F, linetype = "dashed", color = '#39B600') +
    theme_classic() + 
    labs(title = gsub(x = pathway,pattern = 'HALLMARK_',replacement = '') %>% 
           gsub(x = ., pattern = '_', replacement = ' ') %>% 
           stringr::str_wrap(.,width = 22) %>% toupper(),
         subtitle = paste0(cor,', ',p),
         x = 'Macrophage Proportion', y = 'NES') +
    theme(aspect.ratio = 1, plot.title = element_text(size = 14, face = 'bold'), 
          plot.subtitle = element_text(size = 12, hjust = 0.1, vjust = 0.5),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12), legend.title = element_blank()) 
}) %>% wrap_plots(ncol = 4) + 
  # guide_area() + 
  plot_layout(guides = 'collect')
ggsave('Figure/ssGSEA_select.svg',height = ceiling((length(select_pathway)+0)/4)*4, width = 4*4, dpi = 300)
  


## plot
SpatialFeaturePlot(se,features = names(gs), crop = F)
SpatialFeaturePlot(subset(se, Branch%in%c('CPCC','CPCC_CPA') & annotation == 'CPCC'),
                   features = names(gs), crop = F)
ridgeEnrichment(se@meta.data, gene.set = "ZR_up", group = "Branch", add.rug = TRUE)

for(pathway in colnames(ES)){
  SpatialFeaturePlot(subset(se, Branch%in%c('CPCC','CPCC_CPA') & annotation == 'CPCC'),
                     features = pathway, crop = F, pt.size.factor = 1.4)
  ggsave(paste0('Figure/ssGSEA/',pathway,'.tiff'),dpi=300,height=6,width=6)
}

ggsave(wrap_plots(p, ncol = 4), filename = "Figure/ssGSEAplots.pdf", height = 60, width = 28, dpi = 300, limitsize = F)

## get significance
output <- getSignificance(se@meta.data, gene.sets = names(gs), 
                          group = "annotation", fit = "ANOVA") 
output %>% dim()
output %>% filter(p.value < 0.01 & FDR < 0.05) %>% arrange(p.value)


subset(se, Branch%in%c('CPCC','CPCC_CPA') & annotation == 'CPCC') %>% 
  SpatialDimPlot(.,cells.highlight = .@meta.data %>% top_n(n=5,wt = Macrophage) %>% rownames(),
                 crop = F, pt.size.factor = 1.4)

# GSVA------------------------------------------------------------
mat <- GetAssayData(se,slot = "data", assay = "SCT")

# run GSVA
gsva.es <- gsva(mat, gs, parallel.sz = 8)
saveRDS(gsva.es,"RDS/gsva_es.rds")
all.equal(colnames(gsva.es),colnames(se))

for(i in 1:nrow(gsva.es)){
  se@meta.data[rownames(gsva.es)[i]] <- gsva.es[i,]
}
p <- list()
for(i in 1:nrow(gsva.es)){
  p[[i]] <- SpatialFeaturePlot(se,rownames(gsva.es)[i])
}
ggsave(wrap_plots(p, ncol = 4), filename = "Figure/GSVAplots.pdf", height = 60, width = 28, dpi = 300, limitsize = F)