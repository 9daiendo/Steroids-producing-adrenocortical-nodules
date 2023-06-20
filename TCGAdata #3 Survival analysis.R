# Survival analysis

# setup--------------------------------------------------------
rm(list = ls()) ; gc() ; gc()

easypackages::libraries(c("tidyverse",'recount3','edgeR','limma',
                          'AnnotationDbi','org.Hs.eg.db','patchwork',
                          'survival','survminer','msigdbr','fgsea'))
wdir <- "set to your workdir"
setwd(wdir)
load(file = 'RDS/TCGA_ACC_SPN_R_signature.rda')

# significance DEG
lfc_thr <- 0.5 ; p_thr <- 0.1

# survival analysis--------------------------------------------------------------------------------
fit <- survfit(formula = Surv(time = OS.Time, event = OS) ~ group, data = dat) 
p <- ggsurvplot(fit = fit, data = dat, 
                conf.int = F, conf.int.style = "ribbon", #step, ribbon
                conf.int.alpha = 0.15, 
                pval = T, pval.size = 5, 
                title = "", legend.title = "", xlab = "Time in days",
                risk.table = F, risk.table.y.text = F, risk.table.fontsize = 3, risk.table.title = "")

lab <- sapply(strsplit(names(fit$strata),'[=]'),'[[',2)
lab_col <- scales::hue_pal()(n=2) ; names(lab_col) <- names(fit$strata)
p$plot + scale_y_continuous(labels = scales::percent) + theme_classic() + 
  ggtitle('SPN-R signature') + scale_color_manual(values = lab_col, labels = lab) +
  theme(aspect.ratio = 1, legend.text = element_text(size = 10), legend.position = 'right',
        plot.title = element_text(size = 20, hjust = 0, face = 'bold'),
        axis.text = element_text(size = 10),axis.title = element_text(size = 15))
  
ggsave('Figure/TCGA/KaplanMeierPlot.pdf', height = 5, width = 6, dpi = 300)
ggsave('Figure/TCGA/KaplanMeierPlot.svg', height = 5, width = 6, dpi = 300)
  
# DE analysis--------------------------------------------------------------------------------
group1 <- 'High' ; group2 <- 'Low'
contrast <- paste0(group1,'-',group2)
group <- factor(dat$group,levels = c(group2,group1))
design <- model.matrix(~group)
d <- estimateDisp(d,design)
fit <- glmFit(d, design) ; lrt <- glmLRT(fit,coef = 2)
result <- topTags(lrt,  n = Inf, adjust.method="BH") %>% as.data.frame() %>% 
  arrange(desc(logFC)) %>% dplyr::select(logFC, padj = FDR,GENENAME,ENSEMBL,ENTREZID) %>% rownames_to_column(var = 'SYMBOL')

result %>% filter(abs(logFC) > lfc_thr & padj < p_thr) %>% 
  xlsx::write.xlsx(file = "Result/TCGA/DEG_SPN_R_signature.xlsx")

# violinplot--------------------------------------------------------------------------------
exp_df <- lcpm[gene.set[gene.set%in%rownames(lcpm)],] %>% t() %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "sample") %>%
  mutate(group = factor(group,levels = c('High','Low')))

lapply(gene.set[gene.set%in%rownames(lcpm)],function(gene){
  logFC <- result %>% filter(SYMBOL == gene) %>% .$logFC
  padj <- result %>% filter(SYMBOL == gene) %>% .$padj
  lab <- ifelse(abs(logFC) > lfc_thr & padj < p_thr, "*","")
  exp_df %>% dplyr::select(gene = dplyr::all_of(gene), group, sample) %>% 
    ggplot(aes(x = group, y = gene)) + 
    geom_violin(aes(fill = group), color = NA, alpha = 0.6) + 
    # scale_fill_manual(values = col_pal, name = "age") +
    geom_boxplot(width = 0.07, outlier.colour = NA, position = position_dodge(width = 0.9)) +
    ggsignif::geom_signif(comparisons = list(c("High","Low")), map_signif_level = F, 
                          annotations = lab, textsize = 8, fontface = "bold", vjust = 0.5) +
    theme_classic() + labs(title = gene, x = "", y = "Expression Level") +
    theme(plot.title = element_text(size = 20, face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15),
          axis.text.x = element_text(size = 15),
          legend.position = "none",
          aspect.ratio = 1)
}) %>% wrap_plots(ncol = 4) + plot_layout(guides = "collect")
ggsave('Figure/TCGA/Violin_SPN_R_signature.svg', height = 4*ceiling(length(gene.set)/4), width = 4*4, dpi = 300, limitsize = F)
ggsave('Figure/TCGA/Violin_SPN_R_signature.pdf', height = 4*ceiling(length(gene.set)/4), width = 4*4, dpi = 300, limitsize = F)

# Volcano plot--------------------------------------------------------------------------------
toptable <- result
up <- toptable %>% filter(logFC > lfc_thr & padj < p_thr) %>% nrow()
down <- toptable %>% filter(logFC < -lfc_thr & padj < p_thr) %>% nrow()
select_gene <- toptable %>% filter(abs(logFC) > lfc_thr & padj < p_thr) %>% group_by(logFC > 0) %>% top_n(n = 5, wt = abs(logFC)) %>% .$SYMBOL
highlight_df <- toptable %>% filter(SYMBOL %in% unique(c(gene.set,select_gene)))

vol_col_pal <- c('High' = "#F8766D", 'Low' = "#00BFC4", "n.p." = "lightgrey")
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

ggsave(p, filename = 'Figure/TCGA/Volcano_SPN_R_signature.pdf', height = 6, width = 6, dpi = 300)
ggsave(p, filename = 'Figure/TCGA/Volcano_SPN_R_signature.pdf', height = 6, width = 6, dpi = 300)

# GSEA analysis--------------------------------------------------------------------------------
human.genes <- msigdbr(species = "Homo sapiens")
table(human.genes$gs_subcat)
# genesets.interest <- human.genes %>% filter(gs_cat %in% c("H")) 
genesets.interest <- human.genes %>% filter(gs_subcat %in% c("CP:KEGG"))

pathways.interest <- genesets.interest %>% split(x = .$gene_symbol, f = .$gs_name)
pathways.interest$Androgen <- c('CYB5A','SULT2A1','CYP17A1')
names(pathways.interest) %>% length()

stats <- toptable %>% arrange(desc(logFC)) %>% dplyr::select(SYMBOL, logFC) %>% deframe()
gsea_res <- fgsea(pathways = pathways.interest, stats = stats) 
# gsea_res %>% filter(pval < 0.05 & padj < 0.05) %>% arrange(-NES) %>% view()

gsea_res %>% 
  filter(pval < 0.05 & padj < 0.05) %>% 
  # group_by(NES>0) %>% arrange(-abs(NES)) %>% dplyr::slice(1:20) %>% ungroup() %>%
  mutate(pathway = str_sub(pathway, start = 6)) %>%
  mutate(pathway = gsub(pathway,pattern = '_', replacement = ' ')) %>% 
  ggplot(aes(x = NES, y = fct_reorder(pathway,NES), fill = padj)) +
  geom_bar(stat = 'identity') +
  scale_fill_gradientn(
    colours = rev(pals::brewer.ylgn(n = 100)),
    # colours = rev(viridis::viridis(n = 100)),
    name = 'adjusted P', guide = guide_colorbar(reverse = T)) +
  scale_y_discrete(labels = scales::wrap_format(30)) +
  xlim(c(-max(abs(gsea_res$NES)),max(abs(gsea_res$NES)))) +
  theme_light() + 
  theme(
    aspect.ratio = 2, 
    plot.title = element_blank(), axis.title.y = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 15)
    )
ggsave(filename = 'Figure/TCGA/GSEA_SPN_R_signature.pdf', height = 10, width = 9, dpi = 300)
ggsave(filename = 'Figure/TCGA/GSEA_SPN_R_signature.svg', height = 6, width = 8, dpi = 300)

# Deconvolution analysis--------------------------------------------------------------------------------
tpm %>% write.csv('TIMERinput_tpm.csv')

# load TIMER2.0 result
deconv_res <- read.csv('Result/TCGA/TIMER_tpm.csv')


# summarise subpopulation
deconv_res$cell_type %>% grep('_QUANTISEQ$',.,value = T)
df <- deconv_res[grep('_QUANTISEQ$',deconv_res$cell_type),] %>% 
  mutate(cell_type = gsub(cell_type,pattern = '_QUANTISEQ', replace = '')) %>% 
  mutate(cell_type = ifelse(cell_type %in% grep('Macrophage',.data$cell_type,value = T), 'Macrophage', cell_type)) %>% # macrophageのサブタイプを統合
  mutate(cell_type = ifelse(cell_type %in% grep('T cell',.data$cell_type,value = T), 'T cell', cell_type)) %>% # Tcellのサブタイプを統合
  mutate(cell_type = ifelse(cell_type %in% grep('B cell',.data$cell_type,value = T), 'B cell', cell_type)) %>% # Bcellのサブタイプを統合
  mutate(cell_type = ifelse(cell_type %in% grep('NK cell',.data$cell_type,value = T), 'NK cell', cell_type)) %>% # NKcellのサブタイプを統合
  mutate(cell_type = ifelse(cell_type %in% grep('Myeloid',.data$cell_type,value = T), 'Myeloid dendritic cell', cell_type)) %>% # Myeloidcellのサブタイプを統合
  mutate(cell_type = ifelse(cell_type %in% grep('Mast',.data$cell_type,value = T), 'Mast cell', cell_type)) %>% # Mastcellのサブタイプを統合
  mutate(cell_type = ifelse(cell_type == 'uncharacterized cell', 'Tumor cell', cell_type))
colSums(df[,-1])

df <- df %>% pivot_longer(cols = -cell_type,names_to = 'sample') %>% 
  group_by(cell_type,sample) %>% summarise_each(sum) %>% 
  pivot_wider(names_from = sample, values_from = value) %>% as.data.frame()
rownames(df) <- df$cell_type
laborder <- df %>% dplyr::select(-cell_type) %>% 
  rowMeans() %>% sort(decreasing = T) %>% names()
col_pal <- pals::brewer.paired(n=length(laborder)) ; names(col_pal) <- sort(laborder)

plot.df <- df %>% gather(sample, fraction, -cell_type) %>% 
  left_join(dat %>% rownames_to_column(var = 'sample') %>% dplyr::select(sample,group), by = 'sample') %>% 
  mutate(cell_type = factor(cell_type, levels = laborder)) %>% 
  mutate(sample = factor(sample, levels = rev(colnames(df[,-1]))))
lapply(c('High','Low'),function(x){
  plot.df %>% filter(group == x) %>% 
    ggplot(aes(x = sample, y = fraction, fill = cell_type)) +
    geom_bar(stat = "identity") + scale_y_continuous(labels = scales::percent) +
    ggtitle(x) + ylab('Fraction') +
    scale_fill_manual(values = col_pal,breaks = sort(laborder)) + 
    coord_flip() + theme_light() + theme(aspect.ratio = 1, axis.title.y = element_blank())
}) %>% wrap_plots(ncol = 1) + plot_layout(guides = 'collect')

ggsave(filename = 'Figure/TCGA/deconv_fraction.pdf', height = 14, width = 10, dpi = 300)
ggsave(filename = 'Figure/TCGA/deconv_fraction.svg', height = 14, width = 10, dpi = 300)

# wilcox.test
test.df <- df %>% dplyr::select(-cell_type) %>% t() %>% as.data.frame() 
all.equal(rownames(test.df),rownames(dat))
test.df$group <- dat$group
test.df <- test.df %>% pivot_longer(cols = -group, names_to = 'cell_type', values_to = 'fraction')
test.df %>% group_by(group,cell_type) %>% summarise(n = mean(fraction)*100) %>% arrange(group,-n) %>% mutate(n = round(n,1))

pval.df <- lapply(df$cell_type,function(x){
  tmp <- test.df %>% filter(cell_type == x) %>% rstatix::wilcox_test(formula = fraction ~ group)  
  tmp$cell_type <- x
  return(tmp)
}) %>% do.call(rbind,.)

# volcano plot
logFC.df <- lapply(df$cell_type,function(x){
test.df %>% filter(cell_type == x) %>% 
    group_by(group,cell_type) %>% summarise_each(mean) %>% 
    pivot_wider(names_from = group, values_from = fraction) %>% 
    mutate(logFC = log2(High/Low))
}) %>% do.call(rbind,.)

col_pal <- pals::brewer.paired(n=length(laborder)) ; names(col_pal) <- sort(laborder)
plot.df <- pval.df %>% left_join(logFC.df,by = 'cell_type') %>% dplyr::select(logFC, p, cell_type)
plot.df %>% ggplot(aes(x = logFC, y = -log10(p))) +
  geom_vline(xintercept = 0.5, linetype = 'dashed', color = 'grey') + 
  geom_vline(xintercept = -0.5, linetype = 'dashed', color = 'grey') + 
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey') + 
  geom_point(aes(fill = cell_type),size = 5, pch = 21, color = 'black') +
  scale_fill_manual(values = col_pal) +
  geom_label(aes(x = max(abs(logFC)*1.2), y = 1.4), label = 'high SPN-R', hjust = 'inward', vjust = 'inward') +
  geom_label(aes(x = -max(abs(logFC)*1.2), y = 1.4), label = 'low SPN-R', hjust = 'inward', vjust = 'inward')  +
  ggrepel::geom_text_repel(aes(label = cell_type), box.padding = .5) +
  xlim(c(-max(abs(plot.df$logFC))*1.2,max(abs(plot.df$logFC))*1.2)) + ylim(c(0,NA)) +
  theme_bw() +  
  theme(legend.position = "none", aspect.ratio = 1,
        axis.title = element_text(size = 15), 
        plot.title = element_text(face = "bold", size = 20)) +
  xlab(bquote(~Log[2]~"FC")) + ylab(bquote(~-Log[10]~"P")) 

ggsave(filename = 'Figure/TCGA/deconv_volcano.svg', height = 5, width = 5, dpi = 300)
ggsave(filename = 'Figure/TCGA/deconv_volcano.pdf', height = 5, width = 5, dpi = 300)
