# Adjacent CPA vs nonCPA

# setup--------------------------------------------------------
rm(list = ls()) ; gc() ; gc()

easypackages::libraries(c("tidyverse",'edgeR','limma','readxl','biomaRt','sva',
                          'ComplexHeatmap','colorRamps', 'circlize','pals','RColorBrewer',
                          'msigdbr','fgsea',
                          'AnnotationDbi','org.Hs.eg.db','patchwork'))
wdir <- "set to your workdir"
setwd(wdir)

col_pal <- scales::hue_pal()(n=4) ; names(col_pal) <- c('Adjacent_CPA','SPN','CPA','Adjacent_nonCPA')
steroidogenic <- c("CYP11A1",'STAR',"CYP11B2",'HSD3B2','CYP11B1',"CYP17A1","CYP21A2","CYB5A","SULT2A1")

load(file = 'RDS/bulkRNAseq.rda')
all.equal(colnames(combat_count),rownames(metadata),colnames(d))

# significance DEG
lfc_thr <- 0.5 ; p_thr <- 0.1

# DE analysis--------------------------------------------------------------------------------
group1 <- 'Adjacent_CPA' ; group2 <- 'Adjacent_nonCPA'
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
  geom_vline(xintercept = lfc_thr, linewidth = 0.2, linetype = "dashed") +
  geom_vline(xintercept = -lfc_thr, linewidth = 0.2, linetype = "dashed") +
  geom_hline(yintercept = -log10(p_thr), linewidth = 0.2, linetype = "dashed") +
  geom_text(aes(x = max(abs(logFC)*1.2), y = max(abs(-log10(padj)))), 
            label = paste0('Adj.AC to CPA'," (", up,")"), hjust = 1, vjust = 0.5) + # group1
  geom_text(aes(x = -max(abs(logFC)*1.2), y = max(abs(-log10(padj)))), 
            label = paste0('Adj.AC to nonCPA'," (",down,")"), hjust = 1, vjust = 0.5) + # group2
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

mat <- select_lcpm[steroidogenic[steroidogenic%in%rownames(select_lcpm)],] %>% apply(.,MARGIN = 1,FUN = scale) %>% t()
colnames(mat) <- colnames(select_lcpm)

annotation_col <- data.frame(TissueType = ifelse(select_d$samples$TumorType=='Adjacent_CPA','Adj.AC to CPA','Adj.AC to nonCPA'),
                             Sample = select_d$samples$Label, row.names = 'Sample')

top_annotation <- HeatmapAnnotation(df = annotation_col, 
                                    col = list(TissueType = c('Adj.AC to CPA' = '#F8766D','Adj.AC to nonCPA' = '#C77CFF')),
                                    show_legend = T, show_annotation_name = F)

left_annot_df <- result %>% filter(SYMBOL %in% rownames(mat)) %>% 
  mutate(label = ifelse(abs(logFC) > lfc_thr & padj < p_thr,'*','')) %>% 
  dplyr::select(SYMBOL,label) %>% mutate(col = 'white') %>% column_to_rownames(var = 'SYMBOL')
left_annot_df <- left_annot_df[rownames(mat),]
left_annotation <- rowAnnotation(Star = anno_simple(left_annot_df$col, col = c('white' = 'white'), pch = left_annot_df$label),
                                 show_legend = F, show_annotation_name = F)

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
              top_annotation = top_annotation, left_annotation = left_annotation,
              
              heatmap_legend_param = list(direction =  'vertical', legend_height = unit(40, "mm"),
                                          title_position = 'leftcenter-rot', title = "Scaled Expression"))
draw(ht, merge_legend = T, heatmap_legend_side = "right")

# save plot
Fig <- grid.grabExpr(draw(ht, merge_legend = T, heatmap_legend_side = "right"))
ggsave(Fig, filename = paste0('Figure/bulkRNAseq/Heatmap',group1,'vs',group2,'.svg'), width = 7, height = 4, dpi = 300)
ggsave(Fig, filename = paste0('Figure/bulkRNAseq/Heatmap',group1,'vs',group2,'.pdf'), width = 7, height = 4, dpi = 300)

