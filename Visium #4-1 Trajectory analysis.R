# Trajectory analysis 1
# # # # # # # # # # 
# setup--------------------------------------------------------
rm(list = ls()) ; gc() ; gc()

reticulate::use_condaenv(condaenv = "/Users/norifusa/miniforge3-intel/envs/R_reticulate_intel")
easypackages::libraries(c("tidyverse","Seurat","sctransform","glmGamPoi","patchwork",
                          "SeuratWrappers","ElPiGraph.R","igraph","destiny","phateR",
                          'slingshot','plotly','presto',
                          "monocle3","CytoTRACE","reticulate"
                          ))
wdir <- "set to your workdir"
setwd(wdir)
load(file = "RDS/PostAnnotation.rda")

# Dimention reduction----------------------------------------------------
# Diffusion map
library(plotly)
se <- subset(se,annotation!='Capsule')
se <- se %>% RunPCA(npcs = 50, verbose = FALSE)

set.seed(1234)
dm <- DiffusionMap(Embeddings(se, "pca")[,1:30],n_eigs = 3)
se[["dc"]] <- CreateDimReducObject(embeddings = dm@eigenvectors,
                                   key = "DC_", assay = DefaultAssay(se))
dims <- c(1,3)
DimPlot(se, reduction = "dc", dims = dims,label=T,group.by = 'HistologicalAnnotation')
FeaturePlot(se,features = 'CYB5A', reduction = 'dc', dims = dims) + 
  xlim(c(min(se[['dc']][[,dims[1]]]),max(se[['dc']][[,dims[1]]]))) +
  ylim(c(min(se[['dc']][[,dims[2]]]),max(se[['dc']][[,dims[2]]])))

df <- as.data.frame(cbind(Embeddings(se,reduction = "dc"),se@meta.data))
plot_ly(data = df, x = ~DC_1, y = ~DC_2, z = ~DC_3, color = ~HistologicalAnnotation,
        type="scatter3d", mode="markers")

# PHATE
PHATE <- phate(Embeddings(se, "pca")[,1:30], ndim = 3)
se[["phate"]] <- CreateDimReducObject(embeddings = PHATE$embedding,
                                      key = "PHATE_", assay = DefaultAssay(se))
DimPlot(se,reduction = 'phate', dims = c(1,2), label=T)
df <- as.data.frame(cbind(Embeddings(se,reduction = "phate"),se@meta.data))
plot_ly(x = df$PHATE_1, y = df$PHATE_2, z = df$PHATE_3, color = df$annotation,
        type="scatter3d", mode="markers")
df <- as.data.frame(cbind(Embeddings(se,reduction = "pca"),se@meta.data))
plot_ly(x = df$PC_1, y = df$PC_2, z = df$PC_3, color = df$annotation,
        type="scatter3d", mode="markers")


# create trajectory(Elpigraph)----------------------------------------------------
ExpMat <- GetAssayData(se,assay = 'SCT', slot = 'data') %>% t() %>% as.matrix()
colnames(ExpMat) ; rownames(ExpMat)
tree_data <- as.matrix(Embeddings(se,reduction = "dc")[,c(1,2,3)])
grouplab <- se$annotation
boxplot(dist(tree_data)) # trimming parameter設定の参考に

# trajectory
set.seed(1234)
Tree <- computeElasticPrincipalTree(X = tree_data, 
                                    NumNodes = 30, nReps = 10, 
                                    ProbPoint = 1,TrimmingRadius = .04, DensityRadius = .04,
                                    ICOver = "Density", 
                                    ClusType = "FORK",
                                    ParallelRep = TRUE, n.cores = 8, 
                                    Do_PCA = F,
                                    drawAccuracyComplexity = T, drawEnergy = T, drawPCAView = T
                                    )
Extended <- ExtendLeaves(X = tree_data, TargetPG = Tree[[1]], Mode = "QuantCentroid")

Net <- ConstructGraph(PrintGraph = Extended)
Part <- PartitionData(X = tree_data, NodePositions = Extended$NodePositions)
PrjStruct <- project_point_onto_graph(X = tree_data,
                                      NodePositions = Extended$NodePositions,
                                      Edges = Extended$Edges$Edges,
                                      Partition = Part$Partition)

Paths <- GetSubGraph(Net, "end2end")

source('plot_2d.R') # custom function to plot trajectory in 2D
col_pal <- scales::hue_pal()(n = 4) ; names(col_pal) <- c('AdrenalCortex','SPN','CPA')
col_pal <- c("#F8766D", "#619CFF", "#F564E3", "#00BFC4") ; names(col_pal) <- c('AdrenalCortex','SPN-F','SPN-R','CPA')
plot_2d(PrjStruct = PrjStruct, SeuratObj = se, density_show = F, node_label_show = F, color_palette = col_pal,
        reduction = 'dc', dims = c(1,3), color_cells_by = 'annotation')
ggsave('Figure/DC_Traj.svg',dpi = 300, height = 6, width = 8)

plot_2d(PrjStruct = PrjStruct, SeuratObj = se, density_show = F, node_label_show = F, 
        reduction = 'dc', dims = c(1,3), color_cells_by = 'Pseudotime') +
  theme(axis.text = element_blank(),axis.ticks = element_blank())
ggsave('Figure/DC_Traj_PT.svg',dpi = 300, height = 6, width = 8)


## set start point of trajectory
start_node <- '34' # the node exists in AdrenalCortex
Paths_Norm <- lapply(Paths, function(x){
  if(x[length(x)] == start_node){
    return(rev(x))
  }
  if(x[1] == start_node){
    return(x)
  }
})
Paths_Norm <- Paths_Norm[sapply(Paths_Norm, function(x){!any(is.null(x))})]

Branches <- GetSubGraph(Net, "branches")
TB <- sapply(Paths_Norm, function(x){
  sapply(Branches, function(y) {
    all(as.integer(y) %in% as.integer(x))
  })
})
1*TB

BrPT <- which(degree(Net)>2)
Branches <- lapply(Branches, function(x){
  setdiff(as.integer(names(x)), BrPT)
})
CellsOnBrs <- lapply(Branches, function(x){
  which(Part$Partition %in% x)
})
Table <- sapply(CellsOnBrs, function(x){
  table(factor(grouplab[x], levels = unique(grouplab)))
})
Table ; Paths_Norm

## add path names
names(Paths_Norm) <- c("CPA_lineage", "SPN_lineage", "CPA_lineage")
names(Branches) <- c('CPA_lineage','CPA_lineage','CPA_lineage','SPN_lineage','Root')
names(CellsOnBrs) <- c('CPA_lineage','CPA_lineage','CPA_lineage','SPN_lineage','Root')

## seurat objectにbranchの情報を保存する
combinedf <- data.frame(cell_num = 1:ncol(se), Branch = NA)
for(i in 1:length(CellsOnBrs)){
  df <- reshape2::melt(CellsOnBrs[[i]])
  colnames(df) <- 'cell_num'
  df$Branch <- names(CellsOnBrs)[i]
  combinedf$Branch[combinedf$cell_num%in%df$cell_num] <- df$Branch
}
se$Trajectory <- combinedf %>% arrange(cell_num) %>% .$Branch

branch_col <- RColorBrewer::brewer.pal(name = 'Accent',n=3) ; names(branch_col) <- c('Root','SPN_lineage','CPA_lineage')
subse <- subset(se,annotation == 'SPN' & Trajectory != 'NA') 
embs <- Embeddings(subse, reduction = 'dc')
xlims <- c(min(embs[,1]), max(embs[,1]))
ylims <- c(min(embs[,3]), max(embs[,3]))
plot_2d(PrjStruct = PrjStruct, SeuratObj = subse, density_show = F, node_label_show = F,
          color_palette = branch_col, reduction = 'dc', dims = c(1,3), color_cells_by = 'Trajectory') +
  xlim(xlims) + ylim(ylims)
ggsave('Figure/Trajectory_Spatial.svg',dpi=300,height=6,width=8)


# calculate pseudotime
AllPt <- lapply(Paths_Norm, function(x){
  getPseudotime(ProjStruct = PrjStruct, NodeSeq = names(x))
})
PointsPT <- apply(sapply(AllPt, "[[", "Pt"), 1, function(x){unique(x[!is.na(x)])})
PlotPG(X = tree_data, TargetPG = Extended, GroupsLab = PointsPT, Do_PCA = F)
se$Pseudotime <- PointsPT
SpatialFeaturePlot(se,features = 'Pseudotime')


# save
save(se,PrjStruct,Extended,PointsPT, Paths_Norm,Branches,CellsOnBrs, file = "RDS/Trajectory.rda")


