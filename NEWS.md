# `scLANE` 0.7.3

* Added a function named `embedGenes()` that takes a smoothed counts matrix as input & returns PCA & UMAP embeddings along with a graph-based clustering. 
* Updated the `clusterGenes()` function to be much more efficient as well as changing the distance metric used to be cosine dissimilarity. 
* Added `theme_scLANE()` for output plots. 
* Enhanced documentation. 
* Increased test coverage.

# `scLANE` 0.7.2 

* Added a function named `sortGenesHeatmap()` that aids in the creation of expression cascade heatmaps by sorting genes according to where in pseudotime their peak expression is. 
* Changed the parameter `approx.knot` in the `testDynamic()` function to use (stochasticity-controlled) subsampling instead of `seq()` to reduce candidate knot space. 
* Added `summarizeModels()` to sum up slopes across pseudotime intervals. 

# `scLANE` 0.7.1

* Changed input format of all functions to allow counts matrices formatted as `SingleCellExperiment` or `Seurat` objects, sparse matrices, or dense matrices.
* Updated visualization functions to reflect changes made in `ggplot2` v3.4 (mostly changing the `size` parameter in line-based geoms to be `linewidth` instead). 

# `scLANE` 0.6.3

* Added a `NEWS.md` file to track changes to the package.
