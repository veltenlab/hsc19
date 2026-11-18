Package to decompose HSC gene expression data into 19 gene regulatory programs.

To install, in R run `devtools::install_local("~/Downloads/hsc19-main")` (assuming that the package was downloaded to `~/Downloads/hsc19-main`).

For usage info and to reproduce figures 3, 4 and S8 from Bowness et al., 2026, see package vignettes. We recommend starting with the Ageing vignette, which illustrates the main functionalities of the package. Pre-rendered vignettes are available:

* [Vignette 01: Ageing (Figure S8)] (https://htmlpreview.github.io/?https://github.com/veltenlab/hsc19/blob/main/vignettes/output/01_Ageing.html)
* [Vignette 02: LARRY/DARLIN (Figures 3B-G)] (https://htmlpreview.github.io/?https://github.com/veltenlab/hsc19/blob/main/vignettes/output/02_LARRY_DARLIN.html)
* [Vignette 03: AML survival (Figures 4A-D)] (https://htmlpreview.github.io/?https://github.com/veltenlab/hsc19/blob/main/vignettes/output/03_AML_survival.html)
* [Vignette 03: AML Venetoclax response (Figures 4E-G)] (https://htmlpreview.github.io/?https://github.com/veltenlab/hsc19/blob/main/vignettes/output/04_venetoclax.html)

To reproduce figures 1&2, see folder `manuscript_figures`.

To build the vignettes, you need additional data from figshare. This is downloaded automatically, but fails on some systems. download https://ndownloader.figshare.com/files/42587815 and save it e.g. as `~/Downloads/hsc19-main/vignettes/data_LARRY/seurat_LK_LSK.RDS`. Finally, build vignettes with `devtools::install_local("~/Downloads/hsc19-main/", build_vignettes=T)`. Building vignettes will take approx. 30 minutes on a system with 10 cores.


