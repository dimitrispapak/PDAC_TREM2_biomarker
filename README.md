# Molecular analysis highlights TREM2 as a discriminating biomarker for patients suffering from pancreatic ductal adenocarcinoma
This repository contains the scripts that were utilized in the above mentioned study of PDAC cancer.

## Installation
- The Session information for the R scripts are the following :
<details open>
<summary> <span title='Click to Expand'> Current session info </span> </summary>

```r

- Session info ---------------------------------------------------------------
 setting  value
 version  R version 4.1.3 (2022-03-10)
 os       CentOS Linux 7 (Core)
 system   x86_64, linux-gnu
 ui       X11
 language (EN)
 collate  C
 ctype    C
 tz       Europe/Paris
 date     2025-01-23
 pandoc   2.11 @ /mnt/beegfs/software/anaconda3/2020-11/bin/pandoc

- Packages -------------------------------------------------------------------
 package              * version    date (UTC) lib source
 abind                  1.4-8      2024-09-12 [2] CRAN (R 4.1.3)
 annotate               1.72.0     2021-10-26 [2] Bioconductor
 AnnotationDbi        * 1.56.2     2021-11-09 [2] Bioconductor
 ape                    5.8        2024-04-11 [2] CRAN (R 4.1.3)
 aplot                  0.2.3      2024-06-17 [2] CRAN (R 4.1.3)
 askpass                1.2.1      2024-10-04 [2] CRAN (R 4.1.3)
 base64enc              0.1-3      2015-07-28 [2] CRAN (R 4.1.3)
 basilisk               1.6.0      2021-10-26 [2] Bioconductor
 basilisk.utils         1.6.0      2021-10-26 [2] Bioconductor
 Biobase              * 2.54.0     2021-10-26 [2] Bioconductor
 BiocFileCache          2.2.1      2022-01-23 [2] Bioconductor
 BiocGenerics         * 0.40.0     2021-10-26 [2] Bioconductor
 BiocParallel         * 1.28.3     2021-12-09 [2] Bioconductor
 biomaRt              * 2.50.3     2022-02-03 [2] Bioconductor
 Biostrings             2.62.0     2021-10-26 [2] Bioconductor
 bit                    4.5.0      2024-09-20 [2] CRAN (R 4.1.3)
 bit64                  4.5.2      2024-09-22 [2] CRAN (R 4.1.3)
 bitops                 1.0-9      2024-10-03 [2] CRAN (R 4.1.3)
 blob                   1.2.4      2023-03-17 [2] CRAN (R 4.1.3)
 cachem                 1.1.0      2024-05-16 [2] CRAN (R 4.1.3)
 Cairo                  1.6-2      2023-11-28 [2] CRAN (R 4.1.3)
 caTools                1.18.3     2024-09-04 [2] CRAN (R 4.1.3)
 circlize             * 0.4.16     2024-02-20 [2] CRAN (R 4.1.3)
 class                  7.3-22     2023-05-03 [2] CRAN (R 4.1.3)
 cli                    3.6.3      2024-06-21 [2] CRAN (R 4.1.3)
 clipr                  0.8.0      2022-02-22 [2] CRAN (R 4.1.3)
 clue                   0.3-65     2023-09-23 [2] CRAN (R 4.1.3)
 cluster              * 2.1.6      2023-12-01 [2] CRAN (R 4.1.3)
 clusterProfiler      * 4.2.2      2022-01-13 [2] Bioconductor
 codetools              0.2-20     2024-03-31 [2] CRAN (R 4.1.3)
 colorspace             2.1-1      2024-07-26 [2] CRAN (R 4.1.3)
 ComplexHeatmap       * 2.15.4     2023-07-10 [2] Github (jokergoo/ComplexHeatmap@ae0ec42)
 cowplot                1.1.3      2024-01-22 [2] CRAN (R 4.1.3)
 crayon                 1.5.3      2024-06-20 [2] CRAN (R 4.1.3)
 curl                   5.2.3      2024-09-20 [2] CRAN (R 4.1.3)
 data.table             1.16.2     2024-10-10 [2] CRAN (R 4.1.3)
 DBI                    1.2.3      2024-06-02 [2] CRAN (R 4.1.3)
 dbplyr                 2.5.0      2024-03-19 [2] CRAN (R 4.1.3)
 DelayedArray           0.20.0     2021-10-26 [2] Bioconductor
 deldir                 2.0-4      2024-02-28 [2] CRAN (R 4.1.3)
 DEoptimR               1.1-3      2023-10-07 [2] CRAN (R 4.1.3)
 desc                   1.4.3      2023-12-10 [2] CRAN (R 4.1.3)
 DESeq2               * 1.34.0     2021-10-26 [2] Bioconductor
 details              * 0.3.0      2022-03-27 [1] CRAN (R 4.1.3)
 digest                 0.6.37     2024-08-19 [2] CRAN (R 4.1.3)
 diptest                0.77-1     2024-04-10 [2] CRAN (R 4.1.3)
 dir.expiry             1.2.0      2021-10-26 [2] Bioconductor
 DO.db                  2.9        2024-07-31 [2] Bioconductor
 doParallel             1.0.17     2022-02-07 [2] CRAN (R 4.1.3)
 DOSE                   3.20.1     2021-11-18 [2] Bioconductor
 downloader             0.4        2015-07-09 [2] CRAN (R 4.1.3)
 dplyr                * 1.1.4      2023-11-17 [2] CRAN (R 4.1.3)
 edgeR                * 3.36.0     2021-10-26 [1] Bioconductor
 enrichplot             1.14.2     2022-02-24 [2] Bioconductor
 evaluate               1.0.1      2024-10-10 [2] CRAN (R 4.1.3)
 fansi                  1.0.6      2023-12-08 [2] CRAN (R 4.1.3)
 farver                 2.1.2      2024-05-13 [2] CRAN (R 4.1.3)
 fastmap                1.2.0      2024-05-15 [2] CRAN (R 4.1.3)
 fastmatch              1.1-4      2023-08-18 [2] CRAN (R 4.1.3)
 fgsea                  1.20.0     2021-10-26 [2] Bioconductor
 filelock               1.0.3      2023-12-11 [2] CRAN (R 4.1.3)
 fitdistrplus           1.2-1      2024-07-12 [2] CRAN (R 4.1.3)
 flexmix                2.3-19     2023-03-16 [2] CRAN (R 4.1.3)
 forcats              * 1.0.0      2023-01-29 [2] CRAN (R 4.1.3)
 foreach                1.5.2      2022-02-02 [2] CRAN (R 4.1.3)
 fpc                  * 2.2-13     2024-09-24 [2] CRAN (R 4.1.3)
 fs                     1.6.4      2024-04-25 [2] CRAN (R 4.1.3)
 future                 1.34.0     2024-07-29 [2] CRAN (R 4.1.3)
 future.apply           1.11.2     2024-03-28 [2] CRAN (R 4.1.3)
 gdata                * 3.0.1      2024-10-22 [2] CRAN (R 4.1.3)
 genefilter           * 1.76.0     2021-10-26 [2] Bioconductor
 geneplotter            1.72.0     2021-10-26 [2] Bioconductor
 generics               0.1.3      2022-07-05 [2] CRAN (R 4.1.3)
 GenomeInfoDb         * 1.30.1     2022-01-30 [2] Bioconductor
 GenomeInfoDbData       1.2.7      2024-07-31 [2] Bioconductor
 GenomicRanges        * 1.46.1     2021-11-18 [2] Bioconductor
 GetoptLong             1.0.5      2020-12-15 [2] CRAN (R 4.1.3)
 ggforce                0.4.2      2024-02-19 [2] CRAN (R 4.1.3)
 ggfun                  0.1.7      2024-10-24 [2] CRAN (R 4.1.3)
 ggplot2              * 3.5.1      2024-04-23 [2] CRAN (R 4.1.3)
 ggplotify              0.1.2      2023-08-09 [2] CRAN (R 4.1.3)
 ggraph                 2.2.1      2024-03-07 [2] CRAN (R 4.1.3)
 ggrepel              * 0.9.6      2024-09-07 [2] CRAN (R 4.1.3)
 ggridges               0.5.6      2024-01-23 [2] CRAN (R 4.1.3)
 ggtext               * 0.1.2      2022-09-16 [2] CRAN (R 4.1.3)
 ggtree                 3.2.0      2021-10-26 [2] Bioconductor
 GlobalOptions          0.1.2      2020-06-10 [2] CRAN (R 4.1.3)
 globals                0.16.3     2024-03-08 [2] CRAN (R 4.1.3)
 glue                   1.8.0      2024-09-30 [2] CRAN (R 4.1.3)
 GO.db                  3.14.0     2024-07-31 [2] Bioconductor
 goftest                1.2-3      2021-10-07 [2] CRAN (R 4.1.3)
 GOSemSim               2.20.0     2021-10-26 [2] Bioconductor
 gplots               * 3.2.0      2024-10-05 [2] CRAN (R 4.1.3)
 graphlayouts           1.2.0      2024-09-24 [2] CRAN (R 4.1.3)
 gridExtra            * 2.3        2017-09-09 [2] CRAN (R 4.1.3)
 gridGraphics           0.5-1      2020-12-13 [2] CRAN (R 4.1.3)
 gridtext               0.1.5      2022-09-16 [2] CRAN (R 4.1.3)
 gtable                 0.3.5      2024-04-22 [2] CRAN (R 4.1.3)
 gtools                 3.9.5      2023-11-20 [2] CRAN (R 4.1.3)
 hdf5r                  1.3.11     2024-07-07 [2] CRAN (R 4.1.3)
 hms                    1.1.3      2023-03-21 [2] CRAN (R 4.1.3)
 htmltools              0.5.8.1    2024-04-04 [2] CRAN (R 4.1.3)
 htmlwidgets            1.6.4      2023-12-06 [2] CRAN (R 4.1.3)
 httpuv                 1.6.15     2024-03-26 [2] CRAN (R 4.1.3)
 httr                   1.4.7      2023-08-15 [2] CRAN (R 4.1.3)
 ica                    1.0-3      2022-07-08 [2] CRAN (R 4.1.3)
 igraph                 2.1.1      2024-10-19 [2] CRAN (R 4.1.3)
 IRanges              * 2.28.0     2021-10-26 [2] Bioconductor
 IRdisplay              1.1        2022-10-10 [2] local
 IRkernel               1.3.2      2023-01-20 [2] CRAN (R 4.1.3)
 irlba                  2.3.5.1    2022-10-03 [2] CRAN (R 4.1.3)
 iterators              1.0.14     2022-02-05 [2] CRAN (R 4.1.3)
 jsonlite               1.8.9      2024-09-20 [2] CRAN (R 4.1.3)
 KEGGREST               1.34.0     2021-10-26 [2] Bioconductor
 kernlab                0.9-33     2024-08-13 [2] CRAN (R 4.1.3)
 KernSmooth             2.23-24    2024-05-17 [2] CRAN (R 4.1.3)
 knitr                * 1.48       2024-07-07 [1] CRAN (R 4.1.3)
 later                  1.3.2      2023-12-06 [2] CRAN (R 4.1.3)
 lattice                0.22-6     2024-03-20 [2] CRAN (R 4.1.3)
 lazyeval               0.2.2      2019-03-15 [2] CRAN (R 4.1.3)
 leiden                 0.4.3.1    2023-11-17 [2] CRAN (R 4.1.3)
 lifecycle              1.0.4      2023-11-07 [2] CRAN (R 4.1.3)
 limma                * 3.50.3     2022-04-07 [2] Bioconductor
 listenv                0.9.1      2024-01-29 [2] CRAN (R 4.1.3)
 lmtest                 0.9-40     2022-03-21 [2] CRAN (R 4.1.3)
 locfit                 1.5-9.10   2024-06-24 [2] CRAN (R 4.1.3)
 lubridate            * 1.9.3      2023-09-27 [2] CRAN (R 4.1.3)
 maftools             * 2.10.05    2022-02-24 [1] Bioconductor
 magrittr               2.0.3      2022-03-30 [2] CRAN (R 4.1.3)
 MASS                   7.3-58.3   2023-03-07 [2] CRAN (R 4.1.3)
 Matrix                 1.5-3      2022-11-11 [2] CRAN (R 4.1.3)
 MatrixGenerics       * 1.6.0      2021-10-26 [2] Bioconductor
 matrixStats          * 1.4.1      2024-09-08 [1] CRAN (R 4.1.3)
 mclust                 6.1.1      2024-04-29 [2] CRAN (R 4.1.3)
 memoise                2.0.1      2021-11-26 [2] CRAN (R 4.1.3)
 mgcv                 * 1.9-1      2023-12-21 [2] CRAN (R 4.1.3)
 mime                   0.12       2021-09-28 [2] CRAN (R 4.1.3)
 miniUI                 0.1.1.1    2018-05-18 [2] CRAN (R 4.1.3)
 modeltools             0.2-23     2020-03-05 [2] CRAN (R 4.1.3)
 munsell                0.5.1      2024-04-01 [2] CRAN (R 4.1.3)
 nlme                 * 3.1-166    2024-08-14 [2] CRAN (R 4.1.3)
 nnet                   7.3-19     2023-05-03 [2] CRAN (R 4.1.3)
 openssl                2.2.2      2024-09-20 [2] CRAN (R 4.1.3)
 org.Hs.eg.db         * 3.14.0     2023-02-18 [2] Bioconductor
 parallelly             1.38.0     2024-07-27 [2] CRAN (R 4.1.3)
 patchwork            * 1.3.0      2024-09-16 [2] CRAN (R 4.1.3)
 pbapply                1.7-2      2023-06-27 [2] CRAN (R 4.1.3)
 pbdZMQ                 0.3-13     2024-09-17 [2] CRAN (R 4.1.3)
 pheatmap             * 1.0.12     2019-01-04 [2] CRAN (R 4.1.3)
 pillar                 1.9.0      2023-03-22 [2] CRAN (R 4.1.3)
 pkgconfig              2.0.3      2019-09-22 [2] CRAN (R 4.1.3)
 plotly                 4.10.4     2024-01-13 [2] CRAN (R 4.1.3)
 plyr                   1.8.9      2023-10-02 [2] CRAN (R 4.1.3)
 png                    0.1-8      2022-11-29 [2] CRAN (R 4.1.3)
 polyclip               1.10-4     2022-10-20 [2] CRAN (R 4.1.3)
 prabclus               2.3-4      2024-09-24 [2] CRAN (R 4.1.3)
 prettyunits            1.2.0      2023-09-24 [2] CRAN (R 4.1.3)
 progress               1.2.3      2023-12-06 [2] CRAN (R 4.1.3)
 progressr              0.14.0     2023-08-10 [2] CRAN (R 4.1.3)
 promises               1.3.0      2024-04-05 [2] CRAN (R 4.1.3)
 purrr                * 1.0.2      2023-08-10 [2] CRAN (R 4.1.3)
 qvalue                 2.26.0     2021-10-26 [2] Bioconductor
 R6                     2.5.1      2021-08-19 [2] CRAN (R 4.1.3)
 RANN                   2.6.2      2024-08-25 [2] CRAN (R 4.1.3)
 rappdirs               0.3.3      2021-01-31 [2] CRAN (R 4.1.3)
 RColorBrewer         * 1.1-3      2022-04-03 [2] CRAN (R 4.1.3)
 Rcpp                   1.0.13     2024-07-17 [2] CRAN (R 4.1.3)
 RcppAnnoy              0.0.22     2024-01-23 [2] CRAN (R 4.1.3)
 RCurl                  1.98-1.16  2024-07-11 [2] CRAN (R 4.1.3)
 readr                * 2.1.5      2024-01-10 [2] CRAN (R 4.1.3)
 repr                   1.1.7      2024-03-22 [2] CRAN (R 4.1.3)
 reshape2             * 1.4.4      2020-04-09 [2] CRAN (R 4.1.3)
 reticulate             1.39.0     2024-09-05 [2] CRAN (R 4.1.3)
 rjson                  0.2.23     2024-09-16 [2] CRAN (R 4.1.3)
 rlang                  1.1.4      2024-06-04 [2] CRAN (R 4.1.3)
 robustbase             0.99-4-1   2024-09-27 [2] CRAN (R 4.1.3)
 ROCR                   1.0-11     2020-05-02 [2] CRAN (R 4.1.3)
 RSpectra               0.16-2     2024-07-18 [2] CRAN (R 4.1.3)
 RSQLite                2.3.7      2024-05-27 [2] CRAN (R 4.1.3)
 Rtsne                  0.17       2023-12-07 [2] CRAN (R 4.1.3)
 S4Vectors            * 0.32.4     2022-03-24 [2] Bioconductor
 scales                 1.3.0      2023-11-28 [2] CRAN (R 4.1.3)
 scattermore            1.2        2023-06-12 [2] CRAN (R 4.1.3)
 scatterpie             0.2.4      2024-08-28 [2] CRAN (R 4.1.3)
 sctransform            0.4.1      2023-10-19 [2] CRAN (R 4.1.3)
 sessioninfo            1.2.2      2021-12-06 [2] CRAN (R 4.1.3)
 Seurat               * 4.3.0      2022-11-18 [2] CRAN (R 4.1.3)
 SeuratDisk           * 0.0.0.9020 2023-10-09 [2] Github (mojaveazure/seurat-disk@9b89970)
 SeuratObject         * 4.1.3      2022-11-07 [2] CRAN (R 4.1.3)
 shadowtext             0.1.4      2024-07-18 [2] CRAN (R 4.1.3)
 shape                  1.4.6.1    2024-02-23 [2] CRAN (R 4.1.3)
 shiny                  1.9.1      2024-08-01 [2] CRAN (R 4.1.3)
 SingleCellExperiment * 1.16.0     2021-10-26 [2] Bioconductor
 sp                     2.1-4      2024-04-30 [2] CRAN (R 4.1.3)
 spatstat.data          3.1-2      2024-06-21 [2] CRAN (R 4.1.3)
 spatstat.explore       3.3-3      2024-10-22 [2] CRAN (R 4.1.3)
 spatstat.geom          3.3-3      2024-09-18 [2] CRAN (R 4.1.3)
 spatstat.random        3.3-2      2024-09-18 [2] CRAN (R 4.1.3)
 spatstat.sparse        3.1-0      2024-06-21 [2] CRAN (R 4.1.3)
 spatstat.univar        3.0-1      2024-09-05 [2] CRAN (R 4.1.3)
 spatstat.utils         3.1-0      2024-08-17 [2] CRAN (R 4.1.3)
 stringi                1.8.4      2024-05-06 [2] CRAN (R 4.1.3)
 stringr              * 1.5.1      2023-11-14 [2] CRAN (R 4.1.3)
 SummarizedExperiment * 1.24.0     2021-10-26 [2] Bioconductor
 survival               3.7-0      2024-06-05 [2] CRAN (R 4.1.3)
 sva                  * 3.42.0     2021-10-26 [2] Bioconductor
 tensor                 1.5        2012-05-05 [2] CRAN (R 4.1.3)
 tibble               * 3.2.1      2023-03-20 [2] CRAN (R 4.1.3)
 tidygraph              1.3.1      2024-01-30 [2] CRAN (R 4.1.3)
 tidyr                * 1.3.1      2024-01-24 [2] CRAN (R 4.1.3)
 tidyselect             1.2.1      2024-03-11 [2] CRAN (R 4.1.3)
 tidytree               0.4.6      2023-12-12 [2] CRAN (R 4.1.3)
 tidyverse            * 2.0.0      2023-02-22 [2] CRAN (R 4.1.3)
 timechange             0.3.0      2024-01-18 [2] CRAN (R 4.1.3)
 treeio                 1.18.1     2021-11-14 [2] Bioconductor
 tweenr                 2.0.3      2024-02-26 [2] CRAN (R 4.1.3)
 tzdb                   0.4.0      2023-05-12 [2] CRAN (R 4.1.3)
 umap                 * 0.2.10.0   2023-02-01 [2] CRAN (R 4.1.3)
 utf8                   1.2.4      2023-10-22 [2] CRAN (R 4.1.3)
 uuid                   1.2-1      2024-07-29 [2] CRAN (R 4.1.3)
 uwot                   0.2.2      2024-04-21 [2] CRAN (R 4.1.3)
 vctrs                  0.6.5      2023-12-01 [2] CRAN (R 4.1.3)
 viridis              * 0.6.5      2024-01-29 [2] CRAN (R 4.1.3)
 viridisLite          * 0.4.2      2023-05-02 [2] CRAN (R 4.1.3)
 withr                  3.0.1      2024-07-31 [2] CRAN (R 4.1.3)
 xfun                   0.48       2024-10-03 [2] CRAN (R 4.1.3)
 XML                    3.99-0.13  2022-12-04 [2] CRAN (R 4.1.3)
 xml2                   1.3.6      2023-12-04 [2] CRAN (R 4.1.3)
 xtable                 1.8-4      2019-04-21 [2] CRAN (R 4.1.3)
 XVector                0.34.0     2021-10-26 [2] Bioconductor
 yulab.utils            0.1.7      2024-08-26 [2] CRAN (R 4.1.3)
 zellkonverter        * 1.4.0      2021-10-26 [2] Bioconductor
 zlibbioc               1.40.0     2021-10-26 [2] Bioconductor
 zoo                    1.8-12     2023-04-13 [2] CRAN (R 4.1.3)

 [1] /mnt/beegfs/userdata/d_papakonstantinou/R
 [2] /mnt/beegfs/userdata/d_papakonstantinou/R/lib/R/library

------------------------------------------------------------------------------

```

</details>
<br>	


- License

This project is licensed under the GNU General Public License v3.0
