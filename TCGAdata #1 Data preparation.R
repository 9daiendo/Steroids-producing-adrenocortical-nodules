# Data preparation

# setup--------------------------------------------------------
rm(list = ls()) ; gc() ; gc()

easypackages::libraries(c("tidyverse",'recount3','edgeR','limma',
                          'AnnotationDbi','org.Hs.eg.db','patchwork',
                          'survival','survminer','msigdbr','fgsea'))
wdir <- "set to your workdir"
setwd(wdir)

# significance DEG
lfc_thr <- 0.5 ; p_thr <- 0.1

# prepare dataset--------------------------------------------------------
ACC_df <- recount3::create_rse_manual(
  project = "ACC",
  project_home = "data_sources/tcga",
  organism = "human",
  annotation = "gencode_v29",
  type = "gene"
)

df <- ACC_df


## convert ENSEMBL to SYMBOL
total_count <- data.frame(gene_id = rownames(df), count = rowSums(df@assays@data$raw_counts)) %>% mutate(ENSEMBL = sapply(str_split(gene_id,"[.]"),"[[",1))
annot_df <- AnnotationDbi::select(x = org.Hs.eg.db, keys = unique(total_count$ENSEMBL), columns = c("ENSEMBL","SYMBOL"), keytype = "ENSEMBL")
total_count <- total_count %>% left_join(annot_df,by = "ENSEMBL") %>% 
  drop_na(SYMBOL) %>% arrange(-count) %>% 
  distinct(gene_id, .keep_all = T) %>% # remove duplicated gene symbols
  group_by(SYMBOL) %>% dplyr::slice(1:1)

df <- df[total_count$gene_id,]

## count data
countdata <- df@assays@data$raw_counts %>% as.data.frame()
rownames(countdata) <- total_count$SYMBOL

## feature data
featuredata <- AnnotationDbi::select(org.Hs.eg.db,keys = rownames(countdata), keytype = 'SYMBOL',columns = c('ENSEMBL','ENTREZID','GENENAME')) %>%
  distinct(SYMBOL,.keep_all = T) %>% column_to_rownames(var = 'SYMBOL')
featuredata <- featuredata[rownames(countdata),]
featuredata$bp_length <- rowRanges(df)$bp_length

## sample data(DOI: 10.1016/j.ccell.2016.04.002)
TCGAsubtype <- xlsx::read.xlsx('rawdata/TCGAsubtype.xlsx',sheetIndex = 1)

sample_df <- as.data.frame(df@colData)[colnames(countdata),] %>%   
  dplyr::select(tcga.tcga_barcode,tcga.gdc_cases.demographic.gender,
                tcga.xml_stage_event_pathologic_stage, 
                tcga.cgc_case_vital_status,
                tcga.xml_days_to_last_followup, 
                tcga.cgc_case_days_to_death) %>% 
  mutate(
    SAMPLE = str_sub(tcga.tcga_barcode,start = 1, end = 12),
    Sex = tcga.gdc_cases.demographic.gender,
    OS.Time = ifelse(tcga.cgc_case_vital_status=="Alive",tcga.xml_days_to_last_followup,tcga.cgc_case_days_to_death),
    OS = factor(tcga.cgc_case_vital_status,levels = c('Alive','Dead')),
    Stage = factor(tcga.xml_stage_event_pathologic_stage)) %>% 
  left_join(TCGAsubtype, by = 'SAMPLE') %>% 
  dplyr::select(all_of(colnames(TCGAsubtype)), Sex, OS.Time, OS, Stage, ACCsubtype) 
dim(sample_df)

## create edgeR object
colnames(countdata) <- paste0("ACC_",formatC(1:ncol(countdata),width = 3,flag = "0"))
rownames(sample_df) <- colnames(countdata)
d <- DGEList(counts = countdata, genes = featuredata, samples = sample_df)

keep <- filterByExpr(d) # remove low expressed genes
d <- d[keep, , keep.lib.sizes=FALSE]
d <- calcNormFactors(d, method = "TMM")
cpm <- edgeR::cpm(d, log=F) # CPM
lcpm <- edgeR::cpm(d, log=T) # log2

save(d,cpm,lcpm, file = 'RDS/TCGA_ACC.rda')

## calculate TPM
count_df <- as.data.frame(d$counts)
tmp <- count_df %>% rownames_to_column(var = 'SYMBOL') %>% 
  pivot_longer(cols = -SYMBOL, names_to = 'sample', values_to = 'count') %>% 
  left_join(as_tibble(d$genes,rownames = 'SYMBOL'),by = 'SYMBOL') %>% 
  mutate(normalize_genelength = (count / bp_length) * 1000) %>% 
  dplyr::select(SYMBOL,sample,normalize_genelength)
total_count <- tmp %>% group_by(sample) %>% summarise(total_count = sum(normalize_genelength))
tpm <- tmp %>% left_join(total_count, by = 'sample') %>% 
  mutate(tpm = (normalize_genelength / total_count) * 10^6) %>% 
  dplyr::select(SYMBOL,sample,tpm) %>% 
  pivot_wider(names_from = 'sample', values_from = 'tpm') %>% 
  column_to_rownames(var = 'SYMBOL')
write.csv(tpm,'Result/TCGA_tpm.csv')

## view sample info
tableone::CreateTableOne(data = d$samples,
                         vars = c('OS','OS.Time','Stage','ACCsubtype','Sex','OncoSign'), 
                         strata = "ACCsubtype") 