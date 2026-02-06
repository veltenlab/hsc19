Package to decompose HSC gene expression data into 19 gene regulatory programs.

To install, run `devtools::install_github("veltenlab/hsc19")`

For usage info and to reproduce figures 3&4 from Bowness et al., 2026, see package vignettes.

To reproduce figures 1&2, see folder `manuscript_figures`.

To build the vignettes, you need additional data from figshare. This is downloaded automatically, but fails on some systems. If this happens, download the package as a zip file, unpack to a directory of your choice such as `~/Downloads/hsc19-main/`. Then, download https://ndownloader.figshare.com/files/42587815 and save it as `~/Downloads/hsc19-main/vignettes/data_LARRY/seurat_LK_LSK.RDS`. Finally, build vignettes with `devtools::install_local("~/Downloads/hsc19-main/", build_vignettes=T)`. Building vignettes will take approx. 10 minutes.


