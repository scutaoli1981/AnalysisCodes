# Analysis Codes
Custom data analysis codes for Integrated Multi-Omics Analysis Identifies Disrupted Branched-Chain Amino Acid Catabolism as Causal for Sarcopenia

**<font size='5'> 1. mRNA and Metabolite Codes.R. </font>** This file contains codes for differential expression analysis of mRNA data and metabolomics data.

**<font size='5'> 2. Machine Learning Analysis Codes.R. </font>** This file contains codes for random forest analysis.

#### 3. Preparation for local installation of R and packages
This tool is developed with R, so if you want to run it locally, you may do some preparatory work: 
**3.1. Install R.** You can download R from here: [https://www.r-project.org/](https://www.r-project.org/).  
**3.2. Install RStudio.** (Recommendatory but not necessary). You can download RStudio from here: [https://www.rstudio.com/](https://www.rstudio.com/).  
**3.3. Check packages.** After installing R and RStudio, you should check whether you have installed these packages (devtools, openxlsx, limma, qvalue, DESeq2, ggplot2, patchwork, pheatmap, impute, samr, randomForest, caret, pROC). You may run the codes below to check them:  

```r
if(!require(pacman)) install.packages("pacman")
pacman::p_load(devtools, openxlsx, limma, qvalue, DESeq2, ggplot2, patchwork, pheatmap, impute, samr, randomForest, caret, pROC)
```

#### 4. All software dependencies and operating systems (including version numbers)
```r
sessionInfo()
R version 4.4.1 (2024-06-14 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 10 x64 (build 19045)

Matrix products: default


locale:
[1] LC_COLLATE=Chinese (Simplified)_China.utf8  LC_CTYPE=Chinese (Simplified)_China.utf8   
[3] LC_MONETARY=Chinese (Simplified)_China.utf8 LC_NUMERIC=C                               
[5] LC_TIME=Chinese (Simplified)_China.utf8    

time zone: Asia/Shanghai
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] pROC_1.18.5                 caret_6.0-94                lattice_0.22-6             
 [4] randomForest_4.7-1.1        samr_3.0                    patchwork_1.2.0            
 [7] ggplot2_3.5.1               DESeq2_1.44.0               SummarizedExperiment_1.34.0
[10] Biobase_2.64.0              MatrixGenerics_1.16.0       matrixStats_1.3.0          
[13] GenomicRanges_1.56.1        GenomeInfoDb_1.40.1         IRanges_2.38.1             
[16] S4Vectors_0.42.1            BiocGenerics_0.50.0         qvalue_2.36.0              
[19] limma_3.60.4                openxlsx_4.2.6.1           

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3      rstudioapi_0.16.0       jsonlite_1.8.8          magrittr_2.0.3         
  [5] farver_2.1.2            GSA_1.03.3              fs_1.6.4                zlibbioc_1.50.0        
  [9] vctrs_0.6.5             memoise_2.0.1           ggtree_3.12.0           htmltools_0.5.8.1      
 [13] S4Arrays_1.4.1          SparseArray_1.4.8       gridGraphics_0.5-1      parallelly_1.38.0      
 [17] shinyFiles_0.9.3        plyr_1.8.9              httr2_1.0.3             lubridate_1.9.3        
 [21] impute_1.78.0           cachem_1.1.0            igraph_2.0.3            iterators_1.0.14       
 [25] mime_0.12               lifecycle_1.0.4         pkgconfig_2.0.3         Matrix_1.7-0           
 [29] R6_2.5.1                fastmap_1.2.0           future_1.34.0           GenomeInfoDbData_1.2.12
 [33] shiny_1.9.1             digest_0.6.37           aplot_0.2.3             enrichplot_1.24.2      
 [37] colorspace_2.1-1        AnnotationDbi_1.66.0    RSQLite_2.3.7           timechange_0.3.0       
 [41] fansi_1.0.6             httr_1.4.7              polyclip_1.10-7         abind_1.4-5            
 [45] compiler_4.4.1          bit64_4.0.5             withr_3.0.1             BiocParallel_1.38.0    
 [49] viridis_0.6.5           DBI_1.2.3               ggforce_0.4.2           R.utils_2.12.3         
 [53] lava_1.8.0              MASS_7.3-60.2           rappdirs_0.3.3          DelayedArray_0.30.1    
 [57] ModelMetrics_1.2.2.2    tools_4.4.1             ape_5.8                 scatterpie_0.2.3       
 [61] zip_2.3.1               httpuv_1.6.15           future.apply_1.11.2     nnet_7.3-19            
 [65] R.oo_1.26.0             glue_1.7.0              nlme_3.1-164            GOSemSim_2.30.2        
 [69] promises_1.3.0          grid_4.4.1              shadowtext_0.1.4        reshape2_1.4.4         
 [73] recipes_1.1.0           fgsea_1.30.0            generics_0.1.3          gtable_0.3.5           
 [77] class_7.3-22            R.methodsS3_1.8.2       tidyr_1.3.1             data.table_1.16.0      
 [81] tidygraph_1.3.1         utf8_1.2.4              XVector_0.44.0          foreach_1.5.2          
 [85] ggrepel_0.9.5           pillar_1.9.0            stringr_1.5.1           yulab.utils_0.1.7      
 [89] later_1.3.2             splines_4.4.1           dplyr_1.1.4             tweenr_2.0.3           
 [93] treeio_1.28.0           survival_3.6-4          bit_4.0.5               tidyselect_1.2.1       
 [97] GO.db_3.19.1            locfit_1.5-9.10         Biostrings_2.72.1       gridExtra_2.3          
[101] graphlayouts_1.1.1      hardhat_1.4.0           statmod_1.5.0           timeDate_4032.109      
[105] stringi_1.8.4           UCSC.utils_1.0.0        lazyeval_0.2.2          ggfun_0.1.5            
[109] codetools_0.2-20        ggraph_2.2.1            tibble_3.2.1            BiocManager_1.30.24    
[113] ggplotify_0.1.2         cli_3.6.3               rpart_4.1.23            xtable_1.8-4           
[117] munsell_0.5.1           Rcpp_1.0.13             globals_0.16.3          png_0.1-8              
[121] parallel_4.4.1          gower_1.0.1             blob_1.2.4              DOSE_3.30.4            
[125] listenv_0.9.1           viridisLite_0.4.2       tidytree_0.4.6          ipred_0.9-15           
[129] prodlim_2024.06.25      scales_1.3.0            purrr_1.0.2             crayon_1.5.3           
[133] rlang_1.1.4             cowplot_1.1.3           fastmatch_1.1-4         KEGGREST_1.44.1
```
