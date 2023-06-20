# Divided into two groups based on SPN-R signature

# setup--------------------------------------------------------
rm(list = ls()) ; gc() ; gc()

easypackages::libraries(c("tidyverse",'recount3','edgeR','limma',
                          'AnnotationDbi','org.Hs.eg.db','patchwork',
                          'survival','survminer','msigdbr','fgsea'))
wdir <- "set to your workdir"
setwd(wdir)
load(file = 'RDS/TCGA_ACC.rda')

# significance DEG
lfc_thr <- 0.5 ; p_thr <- 0.1


# prepare geneset --------------------------------------------------------------------------------
## logFC top10 genes in SPN-R
SPN_DEG <- xlsx::read.xlsx('Result/DEG_SPN.xlsx',sheetIndex = 1)
SPN_R_signature <- SPN_DEG %>% arrange(logFC) %>% 
  filter(abs(logFC) > 0.25 & padj < 0.05 & group == 'SPN-F') %>% 
  top_n(n = 10, wt = -logFC) %>% pull(feature)


# scale mean expression and divide into two group--------------------------------------------------------------------------------
gene.set <- SPN_R_signature

group_df <- lcpm[gene.set[gene.set%in%rownames(lcpm)],]
if(length(gene.set)>1){group_df <- colMeans(group_df)}
group_df <- group_df %>% scale() %>% as.data.frame() %>% 
  mutate(group = ifelse(V1 > 0, 'High', 'Low'))

all.equal(rownames(d$samples),rownames(group_df))

dat <- d$samples %>% 
  mutate(group = group_df$group,
         OS = ifelse(OS=='Alive',0,1))
dat %>% dplyr::select(-lib.size,-norm.factors) %>% 
  xlsx::write.xlsx('Result/TCGA/sampleinfo.xlsx')

df <- d$samples %>% mutate(group = group_df$group) %>%
  mutate(Stage1 = ifelse(Stage == 'Stage I',T,F),
         Stage2 = ifelse(Stage == 'Stage II',T,F),
         Stage3 = ifelse(Stage == 'Stage III',T,F),
         Stage4 = ifelse(Stage == 'Stage IV',T,F)) %>% 
  tableone::CreateTableOne(
    data = ., strata = "group",
    vars = c('OS','OS.Time','Stage1','Stage2','Stage3','Stage4','ACCsubtype','Sex','OncoSign','mRNA_K4'),
    )
df %>% print() %>% xlsx::write.xlsx('Result/TCGA/sample_table.xlsx')

all.equal(rownames(d$samples),rownames(dat))
d$samples$SPN_R_signature <- dat$group

save(d, cpm, lcpm, gene.set, dat, file = 'RDS/TCGA_ACC_SPN_R_signature.rda')
