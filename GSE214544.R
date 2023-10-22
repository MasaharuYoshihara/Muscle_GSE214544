library(dplyr)
library(Seurat)
library(patchwork)

sample1.pre_count <- ReadMtx(mtx = "GSM6611295_P15306_5001_matrix.mtx.gz",
                       cells = "GSM6611295_P15306_5001_barcodes.tsv.gz",
                       features = "GSM6611295_P15306_5001_features.tsv.gz")

sample1.post_count <- ReadMtx(mtx = "GSM6611296_P15306_5002_matrix.mtx.gz",
                        cells = "GSM6611296_P15306_5002_barcodes.tsv.gz",
                        features = "GSM6611296_P15306_5002_features.tsv.gz")

sample2.pre_count <- ReadMtx(mtx = "GSM6611297_P14601_4004_matrix.mtx.gz",
                       cells = "GSM6611297_P14601_4004_barcodes.tsv.gz",
                       features = "GSM6611297_P14601_4004_features.tsv.gz")

sample2.post_count <- ReadMtx(mtx = "GSM6611298_P14601_4005_matrix.mtx.gz",
                        cells = "GSM6611298_P14601_4005_barcodes.tsv.gz",
                        features = "GSM6611298_P14601_4005_features.tsv.gz")

sample3.pre_count <- ReadMtx(mtx = "GSM6611299_P15306_5003_matrix.mtx.gz",
                       cells = "GSM6611299_P15306_5003_barcodes.tsv.gz",
                       features = "GSM6611299_P15306_5003_features.tsv.gz")

sample3.post_count <- ReadMtx(mtx = "GSM6611300_P15306_5004_matrix.mtx.gz",
                        cells = "GSM6611300_P15306_5004_barcodes.tsv.gz",
                        features = "GSM6611300_P15306_5004_features.tsv.gz")

sample1.pre <- CreateSeuratObject(counts = sample1.pre_count, project = "sample1.pre", min.cells = 3, min.features = 200)
sample1.post <- CreateSeuratObject(counts = sample1.post_count, project = "sample1.post", min.cells = 3, min.features = 200)
sample2.pre <- CreateSeuratObject(counts = sample2.pre_count, project = "sample2.pre", min.cells = 3, min.features = 200)
sample2.post <- CreateSeuratObject(counts = sample2.post_count, project = "sample2.post", min.cells = 3, min.features = 200)
sample3.pre <- CreateSeuratObject(counts = sample3.pre_count, project = "sample3.pre", min.cells = 3, min.features = 200)
sample3.post <- CreateSeuratObject(counts = sample3.post_count, project = "sample3.post", min.cells = 3, min.features = 200)

# sample1.pre 14760 features across 2425 samples
# sample1.post 14840 features across 3057 samples
# sample2.pre 16010 features across 4769 samples
# sample2.post 15595 features across 1435 samples
# sample3.pre 16814 features across 15516 samples
# sample3.post 15984 features across 9566 samples

sample.merge <- merge(sample1.pre, y = c(sample1.post, sample2.pre, sample2.post, sample3.pre, sample3.post)
                  , add.cell.ids = c("sample1.pre", "sample1.post", "sample2.pre", "sample2.post","sample3.pre", "sample3.post")
                  , project = "merge")

# sample.merge 18234 features across 36768 samples

unique(sapply(X = strsplit(colnames(sample.merge), split = "_"), FUN = "[", 1))
# table(sample.merge$orig.ident)

sample.merge[["percent.mt"]] <- PercentageFeatureSet(sample.merge, pattern = "^MT-")

VlnPlot(sample.merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

sample.merge <- subset(sample.merge, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 25)

sample.merge <- NormalizeData(sample.merge, normalization.method = "LogNormalize", scale.factor = 10000)

sample.merge <- FindVariableFeatures(sample.merge, selection.method = "vst", nfeatures = 2000)

plot1 <- VariableFeaturePlot(sample.merge)
plot1

all.genes <- rownames(sample.merge)
sample.merge <- ScaleData(sample.merge, features = all.genes)

sample.merge <- RunPCA(sample.merge, features = VariableFeatures(object = sample.merge))

sample.merge <- JackStraw(sample.merge, num.replicate = 100)
sample.merge <- ScoreJackStraw(sample.merge, dims = 1:20)
ElbowPlot(sample.merge)

sample.merge <- FindNeighbors(sample.merge, dims = 1:10)
sample.merge <- FindClusters(sample.merge, resolution = 1.0)

sample.merge <- RunUMAP(sample.merge, dims = 1:15)
DimPlot(sample.merge, reduction = "umap")

# saveRDS(sample.merge, file = "sample_merge_QC.rds")
# sample.merge <- readRDS("sample_merge_QC.rds")

DimPlot(sample.merge, reduction = "umap", label =  TRUE)

DimPlot(sample.merge, reduction = "umap", group.by = "orig.ident")

FeaturePlot(sample.merge, features = "NOTCH1")
FeaturePlot(sample.merge, features = "NOTCH2")

# VarID2 analysis
library(RaceID)

sample2.post_count <- ReadMtx(mtx = "GSM6611298_P14601_4005_matrix.mtx.gz",
                              cells = "GSM6611298_P14601_4005_barcodes.tsv.gz",
                              features = "GSM6611298_P14601_4005_features.tsv.gz")

sample2.post <- CreateSeuratObject(counts = sample2.post_count, project = "sample2.post", min.cells = 3, min.features = 200)

sample2.post[["percent.mt"]] <- PercentageFeatureSet(sample2.post, pattern = "^MT-")

# VlnPlot(sample2.post, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

sample2.post <- subset(sample2.post, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 25)

sample2.post <- NormalizeData(sample2.post, normalization.method = "LogNormalize", scale.factor = 10000)

sample2.post <- FindVariableFeatures(sample2.post, selection.method = "vst", nfeatures = 2000)

# plot1 <- VariableFeaturePlot(sample2.post)
# plot1

all.genes <- rownames(sample2.post)
sample2.post <- ScaleData(sample2.post, features = all.genes)

sample2.post <- RunPCA(sample2.post, features = VariableFeatures(object = sample2.post))

sample2.post <- JackStraw(sample2.post, num.replicate = 100)
sample2.post <- ScoreJackStraw(sample2.post, dims = 1:20)
# ElbowPlot(sample2.post)

sample2.post <- FindNeighbors(sample2.post, dims = 1:10)
sample2.post <- FindClusters(sample2.post, resolution = 1.0)

sample2.post.VarID2.res <- pruneKnn(sample2.post)
sample2.post.VarID2.sc  <- Seurat2SCseq(sample2.post)
sample2.post.VarID2.sc  <- comptsne(sample2.post.VarID2.sc)
plotmap(sample2.post.VarID2.sc)

noise <- compTBNoise(sample2.post.VarID2.res,getExpData(sample2.post.VarID2.sc),no_cores=1)
sample2.post.VarID2.sc <- updateSC(sample2.post.VarID2.sc,noise=noise)
plotexpmap(sample2.post.VarID2.sc,"NOTCH1")
plotexpmap(sample2.post.VarID2.sc,"NOTCH1",noise=TRUE, cex = 1)

sample2.post.VarID2.noise <- noise

# saveRDS(sample2.post.VarID2.res, "sample2_post_VarID2_res.RDS")
# saveRDS(sample2.post.VarID2.sc, "sample2_post_VarID2_sc.RDS")
# saveRDS(sample2.post.VarID2.noise, "sample2_post_VarID2_noise.RDS")

qn <- quantKnn(res = sample2.post.VarID2.res, 
               noise = sample2.post.VarID2.noise, 
               object = sample2.post.VarID2.sc)

# Error in checkForRemoteErrors(val) : 
# one node produced an error: supply both 'x' and 'y' or a matrix-like 'x'

# https://cran.r-project.org/web/packages/RaceID/vignettes/RaceID.html#integrating-seurat-with-varid2raceid3

# sessionInfo()
# R version 4.3.1 (2023-06-16)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Sonoma 14.0
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: Asia/Tokyo
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] BiocManager_1.30.22     RaceID_0.3.3            patchwork_1.1.3        
# [4] Seurat_4.9.9.9060       SeuratObject_4.9.9.9091 sp_2.1-1               
# [7] dplyr_1.1.3            
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3          rstudioapi_0.15.0          
# [3] jsonlite_1.8.7              umap_0.2.10.0              
# [5] magrittr_2.0.3              spatstat.utils_3.0-3       
# [7] farver_2.1.1                zlibbioc_1.46.0            
# [9] vctrs_0.6.4                 ROCR_1.0-11                
# [11] spatstat.explore_3.2-3      RCurl_1.98-1.12            
# [13] askpass_1.2.0               S4Arrays_1.0.6             
# [15] htmltools_0.5.6.1           sctransform_0.4.1          
# [17] parallelly_1.36.0           KernSmooth_2.23-22         
# [19] htmlwidgets_1.6.2           princurve_2.1.6            
# [21] ica_1.0-3                   plyr_1.8.9                 
# [23] plotly_4.10.3               zoo_1.8-12                 
# [25] igraph_1.5.1                mime_0.12                  
# [27] lifecycle_1.0.3             pkgconfig_2.0.3            
# [29] Matrix_1.6-1.1              R6_2.5.1                   
# [31] fastmap_1.1.1               GenomeInfoDbData_1.2.10    
# [33] MatrixGenerics_1.12.3       fitdistrplus_1.1-11        
# [35] future_1.33.0               shiny_1.7.5.1              
# [37] digest_0.6.33               colorspace_2.1-0           
# [39] S4Vectors_0.38.2            tensor_1.5                 
# [41] RSpectra_0.16-1             irlba_2.3.5.1              
# [43] GenomicRanges_1.52.1        vegan_2.6-4                
# [45] labeling_0.4.3              progressr_0.14.0           
# [47] randomForest_4.7-1.1        fansi_1.0.5                
# [49] spatstat.sparse_3.0-2       mgcv_1.9-0                 
# [51] runner_0.4.3                httr_1.4.7                 
# [53] polyclip_1.10-6             abind_1.4-5                
# [55] coop_0.6-3                  compiler_4.3.1             
# [57] withr_2.5.1                 fastDummies_1.7.3          
# [59] MASS_7.3-60                 openssl_2.1.1              
# [61] DelayedArray_0.26.7         permute_0.9-7              
# [63] tools_4.3.1                 lmtest_0.9-40              
# [65] httpuv_1.6.11               future.apply_1.11.0        
# [67] goftest_1.2-3               quadprog_1.5-8             
# [69] glue_1.6.2                  nlme_3.1-163               
# [71] promises_1.2.1              grid_4.3.1                 
# [73] Rtsne_0.16                  cluster_2.1.4              
# [75] reshape2_1.4.4              generics_0.1.3             
# [77] gtable_0.3.4                spatstat.data_3.0-1        
# [79] tidyr_1.3.0                 data.table_1.14.8          
# [81] XVector_0.40.0              utf8_1.2.3                 
# [83] BiocGenerics_0.46.0         spatstat.geom_3.2-7        
# [85] RcppAnnoy_0.0.21            ggrepel_0.9.4              
# [87] RANN_2.6.1                  pillar_1.9.0               
# [89] stringr_1.5.0               spam_2.9-1                 
# [91] RcppHNSW_0.5.0              later_1.3.1                
# [93] splines_4.3.1               lattice_0.21-9             
# [95] survival_3.5-7              FNN_1.1.3.2                
# [97] deldir_1.0-9                tidyselect_1.2.0           
# [99] SingleCellExperiment_1.22.0 locfit_1.5-9.8             
# [101] miniUI_0.1.1.1              pbapply_1.7-2              
# [103] gridExtra_2.3               IRanges_2.34.1             
# [105] SummarizedExperiment_1.30.2 scattermore_1.2            
# [107] stats4_4.3.1                Biobase_2.60.0             
# [109] matrixStats_1.0.0           pheatmap_1.0.12            
# [111] stringi_1.7.12              lazyeval_0.2.2             
# [113] codetools_0.2-19            som_0.3-5.1                
# [115] FateID_0.2.2                tibble_3.2.1               
# [117] cli_3.6.1                   uwot_0.1.16                
# [119] xtable_1.8-4                reticulate_1.34.0          
# [121] munsell_0.5.0               harmony_1.0.3              
# [123] GenomeInfoDb_1.36.4         Rcpp_1.0.11                
# [125] globals_0.16.2              spatstat.random_3.1-6      
# [127] png_0.1-8                   parallel_4.3.1             
# [129] ellipsis_0.3.2              ggplot2_3.4.4              
# [131] dotCall64_1.1-0             bitops_1.0-7               
# [133] listenv_0.9.0               viridisLite_0.4.2          
# [135] scales_1.2.1                ggridges_0.5.4             
# [137] crayon_1.5.2                leiden_0.4.3               
# [139] purrr_1.0.2                 rlang_1.1.1                
# [141] cowplot_1.1.1 

