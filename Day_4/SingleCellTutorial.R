# Suerat tutorial

library(dplyr)
library(Seurat)
library(patchwork)

# Load the prostate dataset
pbmc <- readRDS("~/Code/2022-MSTP-Bioinformatics-Bootcamp/Day_4/prostate.RDS")
pbmc

  #Standard pre-processing workflow
  #The steps below encompass the standard pre-processing workflow for scRNA-seq data in Seurat. These represent the selection and filtration of cells based on QC metrics, data normalization and scaling, and the detection of highly variable features.
  
  #QC and selecting cells for further analysis
  #Seurat allows you to easily explore QC metrics and filter cells based on any user-defined criteria. A few QC metrics commonly used by the community include
  
  #1. The number of unique genes detected in each cell.
  #2. Low-quality cells or empty droplets will often have very few genes
  #3. Cell doublets or multiplets may exhibit an aberrantly high gene count
  #4. Similarly, the total number of molecules detected within a cell (correlates strongly with unique genes)
  #5. The percentage of reads that map to the mitochondrial genome
  #6. Low-quality / dying cells often exhibit extensive mitochondrial contamination
  #    We calculate mitochondrial QC metrics with the PercentageFeatureSet() function, which calculates the percentage of counts originating from a set of features
  #    We use the set of all genes starting with MT- as a set of mitochondrial genes


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("n_genes", "n_counts", "percent.mt"), ncol = 3)

#normalize data
pbmc <- NormalizeData(pbmc)



#  Identification of highly variable features (feature selection)
# We next calculate a subset of features that exhibit high cell-to-cell variation 
# in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others). 
# We and others have found that focusing on these genes in downstream analysis helps to highlight 
# biological signal in single-cell datasets.

# Our procedure in Seurat is described in detail here, and improves on previous versions by 
# directly modeling the mean-variance relationship inherent in single-cell data, and is 
# implemented in the FindVariableFeatures() function. By default, we return 2,000 features per dataset. 
# These will be used in downstream analysis, like PCA.

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# scale data prior to PCA
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# perform PCA 
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

#### do you recognize any marker genes so far? ####

# visualize genes responsible for main principle components
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

# visualize pca
DimPlot(pbmc, reduction = "pca")

#### notice any clusters starting to form? ####

#examine first 5 PCs 

DimHeatmap(pbmc, dims = 1:5, cells = 500, balanced = TRUE)


# Next we determine the ‘dimensionality’ of the dataset
# To overcome the extensive technical noise in any single feature for scRNA-seq data, 
# Seurat clusters cells based on their PCA scores, with each PC essentially representing a 
# ‘metafeature’ that combines information across a correlated feature set. The 
# top principal components therefore represent a robust compression of the dataset. 
# However, how many components should we choose to include? 10? 20? 100?
  
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

#plot jackstraw results
ElbowPlot(pbmc)

#what does this plot tell us? 

#we probably want >20 clusters in our final grouping, but maybe not more than 30

pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.2) # this resolution number here is the most important.
                                              # see if you can alter it until you get ~20 communities

head(Idents(pbmc), 5)

#run umap analysis so you can effectively plot your different clusters
pbmc <- RunUMAP(pbmc, dims = 1:20)

#plot your umap!
DimPlot(pbmc, reduction = "umap")

#woah, that's cool. Look at how isolated some of these clusters are!


#### Finding differentially expressed features (cluster biomarkers) ####
# Seurat can help you find markers that define clusters via differential expression. 
# By default, it identifies positive and negative markers of a single cluster (specified in ident.1), 
# compared to all other cells. FindAllMarkers() automates this process for all clusters, 
# but you can also test groups of clusters vs. each other, or against all cells.

# The min.pct argument requires a feature to be detected at a minimum percentage 
# in either of the two groups of cells, and the thresh.test argument requires a 
# feature to be differentially expressed (on average) by some amount between the 
# two groups. You can set both of these to 0, but with a dramatic increase 
# in time - since this will test a large number of features that are unlikely 
# to be highly discriminatory. 

#find all upregulated cluster markers
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10<-pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

#look at expression of some example features
VlnPlot(pbmc, features = c("MS4A1", "CD79A")) # try with some other genes!
                                              # what type of cell do these markers belong to?

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"),reduction = "umap")
DotPlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))

# write out marker genes into table for examination

write.table(top10,"/Users/tjsears/Code/2022-MSTP-Bioinformatics-Bootcamp/markerGenes.txt",sep="\t",quote=F)

#################################
# AUTOMATIC BASELINE CLASSIFIER #
#################################

# let's start off with an automatic identification function
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

es.max = sctype_score(scRNAseqData = pbmc@assays$originalexp@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

cL_resutls = do.call("rbind", lapply(unique(pbmc@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(pbmc@meta.data[pbmc@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(pbmc@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

pbmc@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  pbmc@meta.data$customclassif[pbmc@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')  

###################################
# END OF AUTOMATIC CLASSIFIACTION #
###################################

#### HMMM WHICH clusters might be T cells? ####

new.cluster.ids <- c("Sperm!",              #0
                     "Cancer cells",        #1
                     "Memory CD8+ T cells", #2
                     "Cancer Cells2",       #3
                     "Plasmacytoid Dendritic cells",  #4
                     "Endothelial",         #5
                     "Non-classical monocytes", #6
                     "unknown1",            #7
                     "unknown2")            #8

names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()



# Now, given a list of true identities, can we accurately pinpoint which is which?

#B cell           Epithelial       CD4 Trm  CD4 naive-cm     CD8 Trm   CD8 cytotoxic             
#DC   Endothelial    Fibroblast        Luminal Epithelial          Mac1 
#Mac2       Mac-MT1   Mac-cycling     Mast cell          Mono    NK CD16neg    NK CD16pos 
#NK cell         Sperm        T cell          Treg 

# LEGEND:
# LE = luminal epithelial
# DC = dendritic cell
# 


