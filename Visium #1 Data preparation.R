#1 Data preparation, Quality control, Normalization

# setup------------------------------------------------------------------------------
rm(list = ls()) ; gc() ; gc()

easypackages::libraries(c("tidyverse","Seurat","sctransform","glmGamPoi","patchwork"))
wdir <- "set to your workdir"
setwd(wdir)

# load data in Seurat------------------------------------------------------------------------------
set.seed(1234)
se <- Load10X_Spatial(data.dir = "rawdata/2022_8_29_visium/SpaceRanger_Outputs/CPA_C20/outs", 
                      image = Read10X_Image(image.dir = "rawdata/2022_8_29_visium/SpaceRanger_Outputs/CPA_C20/outs/spatial"))
se$percent.mt <- PercentageFeatureSet(se, pattern = "^MT-")
# se <- CellCycleScoring(object = se, g2m.features = cc.genes.updated.2019$g2m.genes, s.features = cc.genes.updated.2019$s.genes)

# Annotation data of gross findings
annotation_df <- xlsx::read.xlsx('/Users/norifusa/R_workspace/Visium_CPCC/rawdata/2022_8_29_visium/SpaceRanger_Outputs/CPA_C20/HistologicalAnnotation.xlsx',sheetIndex = 1) %>% 
  column_to_rownames(var = 'Barcode')
glimpse(annotation_df) ; apply(annotation_df,2,table)
se$HistologicalAnnotation <- annotation_df[colnames(se),"HistologicalAnnotation"]
se$HistologicalAnnotation2 <- annotation_df[colnames(se),"HistologicalAnnotation2"]

# remove adipose tissue spot
se <- subset(se, HistologicalAnnotation!="Adipose")
print(paste0("number of spot: ",ncol(se))) ; print(paste0("number of genes: ",nrow(se)))

# Quality control------------------------------------------------------------------------------
VlnPlot(se, features = c("nCount_Spatial","nFeature_Spatial"), group.by = "orig.ident", pt.size = 0) +
  geom_hline(yintercept = 1000, color = "red") &
  theme_classic() & theme(legend.position = "none",axis.title = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank())
FeatureScatter(se,feature1 = "nCount_Spatial",feature2 = "nFeature_Spatial") + geom_hline(yintercept = 1000)
SpatialFeaturePlot(subset(se,nFeature_Spatial>1000),features = c("nCount_Spatial","nFeature_Spatial"))

preQC <- se
postQC <- subset(se,nFeature_Spatial>1000)
print(paste0(ncol(preQC)-ncol(postQC)," spots are removed"))

# Normalization------------------------------------------------------------------------------
set.seed(1234)
se <- postQC
se <- se %>% 
  SCTransform(assay = "Spatial", vst.flavor = "v2", verbose = FALSE) %>% 
  RunPCA(npcs = 50, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) 
DimPlot(se, group.by = 'HistologicalAnnotation')

# save object------------------------------------------------------------------------------
postQC <- se
save(preQC,postQC,file = "RDS/PostDataPreparation.rda")