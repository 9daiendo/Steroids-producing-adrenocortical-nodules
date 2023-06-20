# Adjacent_CPA vs SPN vs CPA

# setup--------------------------------------------------------
rm(list = ls()) ; gc() ; gc()

easypackages::libraries(c("tidyverse",'edgeR','limma','readxl','biomaRt','sva','eulerr',
                          'ComplexHeatmap','colorRamps', 'circlize','pals','RColorBrewer',
                          'msigdbr','fgsea',
                          'AnnotationDbi','org.Hs.eg.db','patchwork'))
wdir <- "set to your workdir"
setwd(wdir)

col_pal <- scales::hue_pal()(n=4) ; names(col_pal) <- c('Adjacent_CPA','SPN','CPA','Adjacent_nonCPA')
steroidogenic <- c("CYP11A1",'STAR',"CYP11B2",'HSD3B2','CYP11B1',"CYP17A1","CYP21A2","CYB5A","SULT2A1")

load(file = 'RDS/bulkRNAseq.rda')

# significance DEG
lfc_thr <- 0.5 ; p_thr <- 0.1

# PCA--------------------------------------------------------
select <- metadata %>% 
  filter(TumorType != 'Adjacent_nonCPA') %>% .$Label
group <- factor(metadata[select,'TumorType'])
design <- model.matrix(~group)
select_d <- d[,select] ; select_d <- estimateDisp(select_d,design)
select_lcpm <- lcpm[,select]

# make normalized expression matrix
exp_df <- select_lcpm
variablefeature <- 1000 # top 1000 variable genes
select_gene <- apply(exp_df, MARGIN = 1, FUN = var) %>% sort(decreasing = T) %>% head(variablefeature) %>% names()

pca <- exp_df[select_gene,] %>% t() %>% prcomp()
summary(pca) # proportion of varianceが寄与率 cumulative proportionが累積寄与率
factoextra::fviz_screeplot(pca, ncp = 30) + theme_classic()
percentage <- round((pca$sdev)^2 / sum((pca$sdev)^2) * 100, 2)
percentage <- paste0(colnames(pca$x), "(", paste(as.character(percentage), "%", ")", sep=""))

df <- select_d$samples %>%
  left_join((pca$x %>% as.data.frame() %>% rownames_to_column(var = "Label")), by = "Label")
df <- exp_df %>% t() %>% as.data.frame() %>% 
  rownames_to_column(var = "Label") %>% 
  inner_join(df, by = "Label")

# PCA plot
df %>% 
  mutate(TumorType = factor(TumorType, levels = c('Adjacent_CPA','SPN','CPA'))) %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  xlim(c(min(df$PC1)*1.2,max(df$PC1)*1.2)) + ylim(c(min(df$PC2)*1.2,max(df$PC2)*1.2)) +
  # xlab(percentage[1]) + ylab(percentage[2]) +
  # ggplot(aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(fill = TumorType), pch = 21, size = 5) + 
  scale_fill_manual(name = 'TissueType', values = col_pal, label = c('Adjacent_CPA' = 'Adj.AC', 'SPN', 'CPA')) +
  # geom_point(aes(fill = PathologyID), pch = 21, size = 5) +
  # geom_point(aes(fill = Batch), pch = 21, size = 5) +
  # ggrepel::geom_text_repel(aes(label = Label),box.padding = 0.6, size = 3) +
  ggforce::geom_mark_ellipse(aes(fill = TumorType),color = NA, alpha = .2) +
  theme_classic() +
  theme(aspect.ratio = 1,
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_text(size = 15, face = "plain"),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12, face = "plain"))

ggsave(filename = "Figure/bulkRNAseq/PCA.svg", width = 7, height = 6, dpi = 300)
ggsave(filename = "Figure/bulkRNAseq/PCA.pdf", width = 7, height = 6, dpi = 300)

# DEG--------------------------------------------------------------------------------
group1 <- 'CPA' ; group2 <- 'SPN'
group1 <- 'CPA' ; group2 <- 'Adjacent_CPA'
group1 <- 'SPN' ; group2 <- 'Adjacent_CPA'

contrast <- paste0(group1,'-',group2)
select <- metadata %>% 
  filter(TumorType %in% c(group1,group2)) %>% .$Label
group <- factor(metadata[select,'TumorType'],levels = c(group2,group1))
design <- model.matrix(~group)

select_d <- d[,select] ; select_d <- estimateDisp(select_d,design)
select_lcpm <- lcpm[,select]

fit <- glmFit(select_d, design) ; lrt <- glmLRT(fit,coef = 2)

result <- topTags(lrt,  n = Inf, adjust.method="BH") %>% as.data.frame() %>% 
  arrange(desc(logFC)) %>% dplyr::select(logFC, padj = FDR,GENENAME,ENSEMBL,ENTREZID) %>% rownames_to_column(var = 'SYMBOL')
result %>% arrange(-logFC) %>% 
  filter(abs(logFC) > lfc_thr & padj < p_thr) %>%
  filter(SYMBOL %in% grep('HLA',SYMBOL,value=T))
result %>% arrange(-logFC) %>% 
  filter(abs(logFC) > lfc_thr & padj < p_thr) %>%
  filter(SYMBOL %in% steroidogenic)

# save result
result %>% filter(abs(logFC) > lfc_thr & padj < p_thr) %>% xlsx::write.xlsx(file = paste0('Result/bulkRNAseq/',contrast,'.xlsx'))

# Volcano plot--------------------------------------------------------------------------------
toptable <- result
up <- toptable %>% filter(logFC > lfc_thr & padj < p_thr) %>% nrow()
down <- toptable %>% filter(logFC < -lfc_thr & padj < p_thr) %>% nrow()
highlight_df <- toptable %>% filter(abs(logFC) > lfc_thr & padj < p_thr) %>% arrange(-abs(logFC)) %>% group_by(logFC > 0) %>% dplyr::slice(1:10)

vol_col_pal <- c(col_pal[c(group1,group2)], "n.p." = "lightgrey")

p <- toptable %>% 
  mutate(col = ifelse(logFC > lfc_thr & padj < p_thr, group1, ifelse(logFC < -lfc_thr & padj < p_thr, group2, "n.p."))) %>%
  ggplot(aes(x = logFC, y = -log10(padj))) + geom_point(aes(color = col), alpha = 0.6) +
  scale_color_manual(values = vol_col_pal) + 
  ggrepel::geom_text_repel(data = highlight_df,aes(label = SYMBOL), box.padding = unit(0.5, "lines"), max.overlaps = 20) +
  geom_vline(xintercept = lfc_thr, size = 0.2, linetype = "dashed") +
  geom_vline(xintercept = -lfc_thr, size = 0.2, linetype = "dashed") +
  geom_hline(yintercept = -log10(p_thr), size = 0.2, linetype = "dashed") +
  geom_text(aes(x = max(abs(logFC)*1.2), y = max(abs(-log10(padj)))), label = paste0(group1," (", up,")"), hjust = 1, vjust = 0.5) +
  geom_text(aes(x = -max(abs(logFC)*1.2), y = max(abs(-log10(padj)))), label = paste0(group2," (",down,")"), hjust = 1, vjust = 0.5) +
  theme_classic() + 
  theme(legend.position = "none", aspect.ratio = 1,
        axis.title = element_text(size = 15), 
        plot.title = element_text(face = "bold", size = 20)) +
  xlab(bquote(~Log[2]~'FC')) + ylab(bquote(~-Log[10]~"P")) +
  xlim(c(-max(abs(toptable$logFC)*1.2),max(abs(toptable$logFC))*1.2)) +
  ylim(c(0,max(-log10(toptable$padj))*1.1)) +
  coord_flip()

# save plot
ggsave(p, filename = paste0('Figure/bulkRNAseq/Volcano_',contrast,'.svg'), height = 6, width = 6, dpi=300)
ggsave(p, filename = paste0('Figure/bulkRNAseq/Volcano_',contrast,'.pdf'), height = 6, width = 6, dpi=300)

# Heatmap--------------------------------------------------------------------------------
col_fun <- colorRamp2(breaks = seq(from = -2, to = 2, length = 100), 
                      colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100))

mat <- select_lcpm[steroidogenic[steroidogenic%in%rownames(select_lcpm)],] %>%
  apply(.,MARGIN = 1,FUN = scale) %>% t()
colnames(mat) <- colnames(select_lcpm)

annotation_col <- data.frame(TissueType = select_d$samples$TumorType,
                             Sample = select_d$samples$Label, row.names = 'Sample')
top_annotation <- HeatmapAnnotation(df = annotation_col, col = list(TissueType = col_pal), 
                                    show_legend = T, show_annotation_name = F)

ht <- Heatmap(matrix = mat, col = col_fun, 
              column_title_rot = 0, column_title_gp = gpar(fontsize = 0), 
              show_heatmap_legend = T, 
              cluster_rows = F, cluster_columns = F, show_row_dend = F, show_column_dend = F,
              column_split = annotation_col$TissueType,
              # column_gap = unit(5,'mm'),
              width = ncol(mat)*unit(8, "mm"), 
              height = nrow(mat)*unit(8, "mm"),
              show_column_names = T, column_names_rot = 90,
              show_row_names = T, row_names_side = 'left',
              top_annotation = top_annotation, 
              
              heatmap_legend_param = list(direction =  'vertical', legend_height = unit(40, "mm"),
                                          title_position = 'leftcenter-rot', title = "Scaled Expression"))
draw(ht, merge_legend = T, heatmap_legend_side = "right")

# save plot
Fig <- grid.grabExpr(draw(ht, merge_legend = T, heatmap_legend_side = "right"))
ggsave(Fig, filename = paste0('Figure/bulkRNAseq/Heatmap',group1,'vs',group2,'.svg'), width = 6, height = 4, dpi = 300)
ggsave(Fig, filename = paste0('Figure/bulkRNAseq/Heatmap',group1,'vs',group2,'.pdf'), width = 6, height = 4, dpi = 300)


# venn diagram--------------------------------------------------------------------------------
# make gene list
group1 <- 'CPA' ; group2 <- 'Adjacent_CPA'
contrast <- paste0(group1,'-',group2)

result <- xlsx::read.xlsx(paste0('Result/bulkRNAseq/',contrast,'.xlsx'),sheetIndex = 1)
DEG_CPAvsAdjacent_CPA <- result %>% dplyr::filter(logFC > lfc_thr & padj < p_thr) %>% .$SYMBOL

group1 <- 'SPN' ; group2 <- 'Adjacent_CPA'
contrast <- paste0(group1,'-',group2)
result <- xlsx::read.xlsx(paste0('Result/bulkRNAseq/',contrast,'.xlsx'),sheetIndex = 1)
DEG_SPNvsAdjacent_CPA <- result %>% dplyr::filter(logFC > lfc_thr & padj < p_thr) %>% .$SYMBOL

set1 <- DEG_CPAvsAdjacent_CPA
set2 <- DEG_SPNvsAdjacent_CPA

gene_df <- data.frame(gene = unique(c(set1,set2)), row.names = 'gene',
                      set1 = unique(c(set1,set2)) %in% set1,
                      set2 = unique(c(set1,set2)) %in% set2)
colnames(gene_df) <- c('CPA', 'SPN')
unname(col_pal[colnames(gene_df)])

# venn で共通の遺伝子の可視化
ggvenn::ggvenn(gene_df, text_size = 5, set_name_size = F,
               columns = colnames(gene_df),
               show_percentage = F, fill_alpha = .4,
               # auto_scale = T,
               stroke_color = 'grey',
               fill_color = unname(col_pal[colnames(gene_df)])) +
  theme(aspect.ratio = 1.5) + coord_flip() 
ggsave('Figure/bulkRNAseq/venn.svg',width = 4, height = 6,dpi = 300)
ggsave('Figure/bulkRNAseq/venn.pdf',width = 4, height = 6,dpi = 300)


# Enrichment analysis(metascape)--------------------------------------------------------------------
gene_df %>% filter(CPA == F & SPN == T) %>% rownames() %>% cat(sep = '\n', file = 'Result/bulkRNAseq/SPNonly.txt')
gene_df %>% filter(CPA == T & SPN == T) %>% rownames() %>% cat(sep = '\n', file = 'Result/bulkRNAseq/SPN_CPA_common.txt')
gene_df %>% filter(CPA == T & SPN == F) %>% rownames() %>% cat(sep = '\n', file = 'Result/bulkRNAseq/CPAonly.txt')

group <- 'CPAonly' # SPNonly CPA_SPN_common
df <- xlsx::read.xlsx(paste0('Result/Metascape/bulkRNAseq/',group,'/metascape_result.xlsx'),sheetName = 'Enrichment') %>%
  mutate(Summary = sapply(strsplit(GroupID,'_'),'[[',2)) %>% 
  mutate(GroupID = as.integer(sapply(strsplit(GroupID,'_'),'[[',1))) %>% 
  mutate(Category = gsub(x = Category, pattern = 'Canonical Pathways', replacement = 'CP'),
         Category = gsub(x = Category, pattern = 'GO Biological Processes', replacement = 'GOBP'),
         Category = gsub(x = Category, pattern = 'KEGG Pathway', replacement = 'KEGG'),
         Category = gsub(x = Category, pattern = 'Reactome Gene Sets', replacement = 'REACTOME'),
         Category = gsub(x = Category, pattern = 'WikiPathways', replacement = 'WIKI')) %>% 
  mutate(Symbols = sapply(Symbols,function(x){
    tmp <- str_split(x,',')
    if(length(unlist(tmp))>=5){
      head(unlist(tmp),n=5) %>% paste0(collapse = ",") %>% paste0(.,'...') %>% return()
    }else{
      unlist(tmp) %>% paste0(collapse = ",") %>% return()
    }
  })) %>% 
  separate(InTerm_InList,c("count","total"),sep="/") %>% 
  mutate(GeneRatio = as.double(count)/as.double(total)) %>%
  mutate(Description = paste0(Category,': ',Description)) %>% # pathwayの頭にCategoryを追加する
  mutate(Description = stringr::str_wrap(Description,width = 35)) %>% 
  mutate(padj = Log.q.value.) %>% 
  mutate(count = as.numeric(count))

select_pathway <- df %>% 
  filter(Category %in% c('GOBP','KEGG')) %>% # GOBPに絞る場合
  filter(Summary != 'Summary') %>%
  group_by(GroupID) %>% arrange(Log.q.value.) %>% dplyr::slice(1:1) %>% ungroup() %>%
  top_n(wt = -GroupID, n = 10) %>% .$Description

# select_pathway <- c(select_pathway,'GOBP: MHC class II protein complex assembly')

df %>% filter(Summary != 'Summary') %>%
  filter(Description %in% select_pathway) %>% 
  ggplot(aes(x = -LogP, y = fct_reorder(Description,-LogP))) + 
  geom_col(fill = ifelse(group == 'SPNonly','#7CAE00',ifelse(group == 'CPAonly','#00BFC4','#3EB662FF')),alpha = .7) + 
  geom_text(aes(label = Symbols),x = 0.1, hjust = 0, size = 4) +
  DOSE::theme_dose(10) + xlab(bquote(~-Log[10]~"P")) +
  theme(axis.title.y = element_blank(),aspect.ratio = 1.5,
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 12,hjust = 0))

ggsave(paste0('Figure/bulkRNAseq/GOplot_',group,'.svg'),dpi=300, height = 6, width = 7)
ggsave(paste0('Figure/bulkRNAseq/GOplot_',group,'.pdf'),dpi=300, height = 6, width = 7)


# GSEA analysis--------------------------------------------------------------------------------
# gene sets
ZF_geneset <- c("LEF1","ATP2B3","CCN3","PCP4","PTGDS","FGF12",
                "MCOLN3","VAT1L","MINAR1","HSD3B2","ASTN1","DNAI3",
                "KCNK2","PEG10","NCAM1","VSNL1","ZNF711","NRK",
                "MCF2","CALN1","ST6GALNAC5","KCNJ5-AS1","LUZP2")
ZR_geneset <- c("FGG","SLC13A5","SLC27A2","ORM1",
                "CYB5A","PAH","TSPAN12","SGPP2",
                "CES1","CYP26B1","SERPINA3","PLEKHA5",
                "RND2","GSTA5","ENPP2","ORM2","PROK1",
                "CYSLTR2","PTGER3","MB","REEP1","CAV2")

pathways.interest <- list('ZF_geneset' = ZF_geneset, 'ZR_geneset' = ZR_geneset)


# fgsea
stats <- result %>% arrange(desc(logFC)) %>% 
  dplyr::select(SYMBOL, logFC) %>% deframe()
gsea_res <- fgsea(pathways = pathways.interest, stats = stats) 
gsea_res %>% 
  filter(pval < 0.05 & padj < 0.1) %>%
  arrange(-NES) 

# plot
lapply(names(pathways.interest),function(select_pathway){

  ES <- gsea_res %>% filter(pathway == select_pathway) %>% .$ES %>% round(.,2)
  padj <- gsea_res %>% filter(pathway == select_pathway) %>% .$padj
  padj <- ifelse(padj < 0.01, '< 0.01', paste0('= ',round(padj,2)))
  range <- gsea_res %>% filter(pathway == select_pathway) %>% .$ES*1.2 %>% c(.,-.)
  
  rnk <- rank(-stats) ; ord <- order(rnk)
  statsAdj <- stats[ord] ; statsAdj <- sign(statsAdj) * (abs(statsAdj)^1)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathways.interest[[select_pathway]], names(statsAdj)))))
  pathway <- sort(pathway)
  tmp <- calcGseaStat(statsAdj, selectedStats = pathway, returnAllExtremes = TRUE)
  top <- max(tmp$tops) ; bottom <- min(tmp$bottoms) ; range <- abs(top - bottom)
  
  label_pos <- data.frame(x = c(0,length(stats)), 
                          y = rep(bottom - range/6, 2), 
                          label = c(group1, group2)) 
  arrow_pos <- data.frame(x = c(1600, length(stats) - 1600),
                          y = rep(bottom - range/20, 2),
                          xend = c(0, length(stats)),
                          yend = rep(bottom - range/20, 2))
  plotEnrichment(pathways.interest[[select_pathway]], stats, gseaParam = 1, ticksSize = 0.2) + 
  geom_text(data = label_pos, mapping = aes(x = x, y = y, label = label), size = 5, fontface = 'bold',
            vjust = "inward", hjust = "inward") +
  geom_segment(data = arrow_pos, aes(x = x, y = y, xend = xend, yend = yend),
               lineend = "square", linejoin = "mitre",
               size = 1, arrow = arrow(length = unit(2, "mm"))) +  
  labs(title = gsub(select_pathway, pattern = '_', replacement = ' '),
         subtitle = paste0('ES = ',ES,', P ', padj)) + 
    xlab('Rank') + ylab('Enrichment Score') + 
    theme(aspect.ratio = 1, 
          plot.title = element_text(size = 20, face = 'bold'),
          plot.subtitle = element_text(size = 15),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 15)) 
}) %>% wrap_plots(ncol = 2)

ggsave(filename = paste0('Figure/bulkRNAseq/GSEA_',contrast,'.svg'), height = 4, width = 7, dpi = 300)
ggsave(filename = paste0('Figure/bulkRNAseq/GSEA_',contrast,'.pdf'), height = 4, width = 7, dpi = 300)

# Violin plot of NR5A1--------------------------------------------------------------------------------
select <- metadata %>% 
  filter(TumorType != 'Adjacent_nonCPA') %>% .$Label
select_d <- d[,select]
select_lcpm <- lcpm[,select]

plot_df <- select_lcpm[c(steroidogenic,'NR5A1'),] %>% t() %>% as.data.frame()
all.equal(rownames(plot_df),colnames(select_d))

plot_df <- cbind(plot_df,select_d$samples)

CPAvsAdj <- xlsx::read.xlsx('Result/bulkRNAseq/CPA-Adjacent_CPA.xlsx',sheetIndex = 1)
CPAvsSPN <- xlsx::read.xlsx('Result/bulkRNAseq/CPA-SPN.xlsx',sheetIndex = 1)
SPNvsAdj <- xlsx::read.xlsx('Result/bulkRNAseq/SPN-Adjacent_CPA.xlsx',sheetIndex = 1)


for(gene in c(steroidogenic,'NR5A1')){
    
  lab_CPAvsSPN <- ifelse(gene%in%CPAvsSPN$SYMBOL,'*','NS')
  lab_CPAvsAdj <- ifelse(gene%in%CPAvsAdj$SYMBOL,'*','NS')
  lab_SPNvsAdj <- ifelse(gene%in%SPNvsAdj$SYMBOL,'*','NS')
  
  plot_df %>% 
    dplyr::select(Expression = all_of(gene), everything()) %>% 
    mutate(TissueType = factor(TumorType, levels = c('Adjacent_CPA','SPN','CPA'))) %>% 
    ggplot(aes(x = TissueType, y = Expression)) +
    geom_violin(aes(fill = TissueType), color = NA, alpha = 0.6) + 
    scale_fill_manual(values = col_pal, name = 'TissueType', label = c('Adjacent_CPA' = 'Adj.AC','SPN','CPA')) +
    geom_boxplot(width = 0.07, outlier.colour = NA, position = position_dodge(width = 0.9)) +
    ggsignif::geom_signif(comparisons = list(c('SPN','Adjacent_CPA'),c('CPA','SPN'),c('CPA','Adjacent_CPA')), 
                          annotations = c(lab_SPNvsAdj,lab_CPAvsSPN,lab_CPAvsAdj),
                          step_increase = 0.1, map_signif_level = F, tip_length = 0) +
    theme_classic() + ggtitle(gene) +
    ylab(bquote(~Log[2]~"CPM")) +
    scale_x_discrete(label = c('Adjacent_CPA' = 'Adj.AC','SPN','CPA')) +
    theme(plot.title = element_text(size = 20, face = 'bold'), aspect.ratio = 1,
          legend.title = element_text(size = 15),legend.text = element_text(size = 12),
          axis.title.x = element_blank(), axis.title.y = element_text(size = 15),
          axis.text.x = element_text(size = 15))
  
  ggsave(paste0('Figure/bulkRNAseq/Violin_',gene,'.pdf'), height = 4, width = 6, dpi = 300)
  ggsave(paste0('Figure/bulkRNAseq/Violin_',gene,'.svg'), height = 4, width = 6, dpi = 300)
}
