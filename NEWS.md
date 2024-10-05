# Changes in v0.8.4

+ Minor bug fixes. 
+ Improved matrix inversion through judicious usage of C++. 
+ Added S3 method `summary.scLANE()` to be used on output from `testDynamic()`. 

# Changes in version 0.8.3

+ Sped up GLMM mode. 
+ Switched GEE mode to use model-based variance and made computation of sandwich variance-covariance matrix conditional to speed things up and reduce memory usage. 

# Changes in version 0.8.2

+ Sped up the NB LASSSO implementation in `fitGLMM()`. 
+ Fixed some errors related to intercept-only `marge` models. 
+ Added DF and KC corrections to new function `biasCorrectGEE()`, usage of which is set to `FALSE` by default in `testDynamic()`. 
+ Added support for `monocle3` objects of class `cell_data_set` throughout the package, notably in `testDynamic()`. 

# Changes in version 0.8.1

+ Added small-sample bias correction method to GEE sandwich variance-covariance matrix, results in smaller Wald test statistics. 
+ Parallelized `getResultsDE()` using `future` backend. 
+ Added a function called `chooseCandidateGenes()` to identify good genes for trajectory DE testing based on mean / SD expression and sparsity. 

# Changes in version 0.8.0

+ Added implicit regularization of selected basis functions to the GLMM mode using a NB LASSO. 
+ Switched candidate knot subsampling to a uniform sequence of candidate knots across pseudotime's support. 

# Changes in version 0.7.9

+ Added `geneProgramDrivers()` function to compute & test correlations of expression with gene module scores.
+ Updated documentation & unit tests.
+ Added `geneProgramSignificance()` function to estimate associations between gene program module scores and pseudotime. 

# Changes in version 0.7.8

+ Added progress bar to `testDynamic()`. 
  + Changed parallel backend in `testDynamic()` from `doParallel` to `doSNOW` in order to make this possible. 
+ Updated documentation with more runnable examples. 
+ Passing `BiocCheck` with no errors. 
+ Reduced set of exported functions to just what's necessary for model fitting & downstream analysis. 
+ Reverted changes made to null model scoring internally that lead to worse knot selection i.e., knot selection is improved now. 

# Changes in version 0.7.7

+ Added DOI badge to README. 
+ Better compression of included datasets.
+ Added `geneProgramScoring()` for module scoring of dynamic gene clusters.
+ Added `plotModelCoefs()` to annotate gene dynamics plots with a table of model coefficients. 
+ Added citation file with link to Zenodo repository (until preprint is up).
+ Added runnable examples to most functions. 

# Changes in version 0.7.6

+ Added [Zenodo tracking](https://doi.org/10.5281/zenodo.10182497).
+ Added simulated dataset to `data/`.

# Changes in version 0.7.5

+ Preparing for BioConductor submission i.e., reformatting code, adding documentation, etc.
+ Added convolution function `npConvolve()` to be used for e.g., heatmap smoothing.

# Changes in version 0.7.4

+ Added the `getKnotDist()` function to pull the set of empirically-identified knots for a user-provided gene set.
+ Minor enhancements & documentation improvements. 

# Changes in version 0.7.3

+ Added a function named `embedGenes()` that takes a smoothed counts matrix as input & returns PCA & UMAP embeddings along with a graph-based clustering. 
+ Updated the `clusterGenes()` function to be much more efficient as well as changing the distance metric used to be cosine distance. 
+ Added `theme_scLANE()` for output plots. 
+ Enhanced documentation. 
+ Increased test coverage.

# Changes in version 0.7.2 

+ Added a function named `sortGenesHeatmap()` that aids in the creation of expression cascade heatmaps by sorting genes according to where in pseudotime their peak expression is. 
+ Changed the parameter `approx.knot` in the `testDynamic()` function to use (stochasticity-controlled) subsampling instead of `seq()` to reduce candidate knot space. 
+ Added `summarizeModels()` to sum up slopes across pseudotime intervals. 

# Changes in version 0.7.1

+ Changed input format of all functions to allow counts matrices formatted as `SingleCellExperiment` or `Seurat` objects, sparse matrices, or dense matrices.
+ Updated visualization functions to reflect changes made in `ggplot2` v3.4 (mostly changing the `size` parameter in line-based geoms to be `linewidth` instead). 

# Changes in version 0.6.3

+ Added a `NEWS.md` file to track changes to the package.
