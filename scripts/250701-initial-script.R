# 250701 initial script

#> 1. load matrix / seurat objects
#> 2. QC
#> 3. identify DGE
#> 4. generate UMAP
#> 5. label cells
#> 6. perform RNA velocity
#> 7. generate a table where for each cluster, identify the nearest
#>    developmental neighbors based on RNA velocity
#> 8. re-run DGE analysis for each cluster, restrict comparison to
#>    the nearest neighbors identified by RNA velocity
#> 9. filter the DGE-by-nearest-neighbors for surface proteins using
#>    verified annotation databases: UniProt, RCSB, GO terms database
#> 10. connect to antibody databases to identify available antibodies


# load libraries
library(Seurat)



#1 load seurat objects

