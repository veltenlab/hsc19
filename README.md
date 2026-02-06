Package to decompose HSC gene expression data into 19 gene regulatory programs.

To install, in R run `devtools::install_local("~/Downloads/hsc19-main")` (assuming that the package was downloaded to `~/Downloads/hsc19-main`).

For usage info and to reproduce figures 3&4 from Bowness et al., 2026, see package vignettes.

To reproduce figures 1&2, see folder `manuscript_figures`.

To build the vignettes, you need additional data from figshare. This is downloaded automatically, but fails on some systems. download https://ndownloader.figshare.com/files/42587815 and save it e.g. as `~/Downloads/hsc19-main/vignettes/data_LARRY/seurat_LK_LSK.RDS`. Finally, build vignettes with `devtools::install_local("~/Downloads/hsc19-main/", build_vignettes=T)`. Building vignettes will take approx. 10 minutes.


