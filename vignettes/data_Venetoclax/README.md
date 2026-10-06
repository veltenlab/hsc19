These inputs accompany `../Venetoclax_Renders.Rmd`:

- `Count_matrix_LSC_MATURE.csv.gz`: supplied bulk RNA-seq count matrix for LSC
  and mature-cell samples. The vignette selects diagnostic LSC samples.
- `E-MTAB-14142.sdrf.txt`: study sample annotations for ArrayExpress accession
  E-MTAB-14142, including treatment-response labels. Blank fields denote missing
  annotations, not nonresponse.
- `Metadata_Lars_10022026.csv`: supplied supplementary sample identifiers,
  genetic risk classifications, transcriptomic scores, and PDX annotations.
  The vignette joins selected score columns to the original study annotations;
  response labels come from the SDRF table.

The input files are copied without modification. All filtering, score recoding,
normalization, and decomposition are shown in the vignette. The gzipped count
matrix uses the repository's existing Git LFS rule.
