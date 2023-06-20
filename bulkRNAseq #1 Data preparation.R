# Bulk RNA-seq Data Processing

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

# Data load--------------------------------------------------------
countdata <- read.csv(file = "rawdata/count.csv", row.names = 1)
metadata <- read.csv(file = 'rawdata/metadata.csv', row.names = 1)
all.equal(colnames(countdata),rownames(metadata))

# make Featuredata--------------------------------------------------------
featuredata <- AnnotationDbi::select(org.Hs.eg.db,keys = rownames(countdata), keytype = 'SYMBOL',columns = c('ENSEMBL','ENTREZID','GENENAME')) %>%
  distinct(SYMBOL,.keep_all = T) %>% column_to_rownames(var = 'SYMBOL')
featuredata <- featuredata[rownames(countdata),]

# Batch correction--------------------------------------------------------
combat_count <- ComBat_seq(counts = as.matrix(countdata), batch = metadata$Batch)


# make DGElist--------------------------------------------------------
d <- DGEList(counts = combat_count, 
             genes = featuredata, 
             samples = metadata
             ) 
keep <- filterByExpr(d) # remove low expressed genes
d <- d[keep, , keep.lib.sizes=FALSE]
d <- calcNormFactors(d, method = "TMM")
cpm <- edgeR::cpm(d, log=F) # CPM
lcpm <- edgeR::cpm(d, log=T) # log2

# save--------------------------------------------------------
save(countdata,metadata,featuredata,combat_count,d,cpm,lcpm, file = 'RDS/bulkRNAseq.rda')