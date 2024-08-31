# Analysis Codes
Data analysis codes for Integrated Multi-Omics Analysis Identifies Disrupted Branched-Chain Amino Acid Catabolism as Causal for Sarcopenia

**<font size='5'> 1. mRNA and Metabolite Codes.R. </font>** This file contains codes for differential expression analysis of mRNA data and metabolomics data.

**<font size='5'> 2. Machine Learning Analysis Codes.R. </font>** This file contains codes for random forest analysis.

## All software dependencies and operating systems (including version numbers)
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
 [1] samr_3.0                    patchwork_1.2.0             ggplot2_3.5.1              
 [4] DESeq2_1.44.0               SummarizedExperiment_1.34.0 Biobase_2.64.0             
 [7] MatrixGenerics_1.16.0       matrixStats_1.3.0           GenomicRanges_1.56.1       
[10] GenomeInfoDb_1.40.1         IRanges_2.38.1              S4Vectors_0.42.1           
[13] BiocGenerics_0.50.0         qvalue_2.36.0               limma_3.60.4               
[16] openxlsx_4.2.6.1           

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3      rstudioapi_0.16.0       jsonlite_1.8.8          magrittr_2.0.3         
  [5] farver_2.1.2            GSA_1.03.3              fs_1.6.4                zlibbioc_1.50.0        
  [9] vctrs_0.6.5             memoise_2.0.1           ggtree_3.12.0           htmltools_0.5.8.1      
 [13] S4Arrays_1.4.1          SparseArray_1.4.8       gridGraphics_0.5-1      shinyFiles_0.9.3       
 [17] plyr_1.8.9              httr2_1.0.3             impute_1.78.0           cachem_1.1.0           
 [21] igraph_2.0.3            mime_0.12               lifecycle_1.0.4         pkgconfig_2.0.3        
 [25] Matrix_1.7-0            R6_2.5.1                fastmap_1.2.0           GenomeInfoDbData_1.2.12
 [29] shiny_1.9.1             digest_0.6.37           aplot_0.2.3             enrichplot_1.24.2      
 [33] colorspace_2.1-1        AnnotationDbi_1.66.0    RSQLite_2.3.7           fansi_1.0.6            
 [37] httr_1.4.7              polyclip_1.10-7         abind_1.4-5             compiler_4.4.1         
 [41] bit64_4.0.5             withr_3.0.1             BiocParallel_1.38.0     viridis_0.6.5          
 [45] DBI_1.2.3               ggforce_0.4.2           R.utils_2.12.3          MASS_7.3-60.2          
 [49] rappdirs_0.3.3          DelayedArray_0.30.1     tools_4.4.1             ape_5.8                
 [53] scatterpie_0.2.3        zip_2.3.1               httpuv_1.6.15           R.oo_1.26.0            
 [57] glue_1.7.0              nlme_3.1-164            GOSemSim_2.30.2         promises_1.3.0         
 [61] grid_4.4.1              shadowtext_0.1.4        reshape2_1.4.4          fgsea_1.30.0           
 [65] generics_0.1.3          gtable_0.3.5            R.methodsS3_1.8.2       tidyr_1.3.1            
 [69] data.table_1.16.0       tidygraph_1.3.1         utf8_1.2.4              XVector_0.44.0         
 [73] ggrepel_0.9.5           pillar_1.9.0            stringr_1.5.1           yulab.utils_0.1.7      
 [77] later_1.3.2             splines_4.4.1           dplyr_1.1.4             tweenr_2.0.3           
 [81] treeio_1.28.0           lattice_0.22-6          bit_4.0.5               tidyselect_1.2.1       
 [85] GO.db_3.19.1            locfit_1.5-9.10         Biostrings_2.72.1       gridExtra_2.3          
 [89] graphlayouts_1.1.1      statmod_1.5.0           stringi_1.8.4           UCSC.utils_1.0.0       
 [93] lazyeval_0.2.2          ggfun_0.1.5             codetools_0.2-20        ggraph_2.2.1           
 [97] tibble_3.2.1            BiocManager_1.30.24     ggplotify_0.1.2         cli_3.6.3              
[101] xtable_1.8-4            munsell_0.5.1           Rcpp_1.0.13             png_0.1-8              
[105] parallel_4.4.1          blob_1.2.4              DOSE_3.30.4             viridisLite_0.4.2      
[109] tidytree_0.4.6          scales_1.3.0            purrr_1.0.2             crayon_1.5.3           
[113] rlang_1.1.4             cowplot_1.1.3           fastmatch_1.1-4         KEGGREST_1.44.1
```
