
- [`scLANE`](#sclane)
  - [Installation](#installation)
  - [Model Structure](#model-structure)
- [Usage](#usage)
  - [Input Data](#input-data)
  - [Trajectory DE Testing](#trajectory-de-testing)
  - [Downstream Analysis &
    Visualization](#downstream-analysis--visualization)
- [Conclusions & Best Practices](#conclusions--best-practices)
- [Contact Information](#contact-information)
- [References](#references)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# `scLANE`

<!-- badges: start -->

[![codecov](https://codecov.io/gh/jr-leary7/scLANE/branch/main/graph/badge.svg?token=U2U5RTF2VW)](https://codecov.io/gh/jr-leary7/scLANE)
[![R-CMD-check](https://github.com/jr-leary7/scLANE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jr-leary7/scLANE/actions/workflows/R-CMD-check.yaml)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![release](https://img.shields.io/github/v/release/jr-leary7/scLANE?color=purple)
[![CodeFactor](https://www.codefactor.io/repository/github/jr-leary7/sclane/badge)](https://www.codefactor.io/repository/github/jr-leary7/sclane)
![last
commit](https://img.shields.io/github/last-commit/jr-leary7/scLANE/main?color=darkgreen)
<!-- badges: end -->

## Installation

You can install the most recent version of `scLANE` with:

``` r
remotes::install_github("jr-leary7/scLANE")
```

## Model Structure

The `scLANE` package enables users to accurately determine differential
expression of genes over pseudotime or latent time, and to characterize
gene’s dynamics using interpretable model coefficients. `scLANE` builds
upon the [`marge` modeling
framework](https://github.com/JakubStats/marge), allowing users to
characterize their trajectory’s effects on mRNA expression using
negative binomial
[GLMs](https://en.wikipedia.org/wiki/Generalized_linear_model),
[GEEs](https://en.wikipedia.org/wiki/Generalized_estimating_equation),
or
[GLMMs](https://en.wikipedia.org/wiki/Generalized_linear_mixed_model),
depending on the experimental design & biological questions of interest.
This modeling framework is an extension of the well-known [Multivariate
Adapative Regression Splines
(MARS)](https://en.wikipedia.org/wiki/Multivariate_adaptive_regression_spline)
method, which uses truncated power basis splines to build nonlinear
models using a [generalized cross validation (GCV)
criterion](https://doi.org/10.48550/arXiv.1706.02495).

A quickstart guide on how to use `scLANE` with simulated data continues
below, and a more detailed vignette showcasing its performance on real
data can be found
[here](https://jr-leary7.github.io/quarto-site/tutorials/scLANE_Trajectory_DE.html).

# Usage

Our method relies on a relatively simple test in order to define whether
a given gene is differentially expressed (or “dynamic”) over the
provided trajectory. While the exact structure of the test differs by
model backend, the concept is the same: the spline-based NB GLM / GEE /
GLMM is treated as the alternate model, and a null model is fit using
the corresponding model backend. If the GLM backend is used, then the
null model is simply an intercept-only NB GLM; the GEE backend fits an
intercept-only model with the same working correlation structure as the
alternate model, and if the GLMM backend is used then the null model is
an intercept-only model with random intercepts for each subject. The
alternate hypothesis is thus that at least one of the estimated
coefficients is significantly different from zero. We predict a given
gene to be dynamic if the adjusted *p*-value of the test is less than an
*a priori* threshold; the default threshold is 0.01, and the default
adjustment method is [the Holm
correction](https://en.wikipedia.org/wiki/Holm–Bonferroni_method).

## Input Data

`scLANE` has a few basic inputs with standardized formats, though
different modeling frameworks require different inputs; for example, the
GEE & GLMM backends require a vector of subject IDs corresponding to
each cell. We’ll start by simulating some scRNA-seq counts data with a
subset of genes being dynamic over a trajectory using the [`scaffold` R
package](https://github.com/rhondabacher/scaffold). As part of our
simulation studies I wrote the function sourced below to generate
multi-subject scRNA-seq datasets with trajectories partially conserved
between subjects.

``` r
source("https://raw.githubusercontent.com/jr-leary7/scLANE-Sims/main/R/functions_simulation.R")
```

We’ll also need to load a couple dependencies & resolve a function
conflict.

``` r
library(dplyr)
library(scran)
library(purrr)
library(scater)
library(scLANE)
library(ggplot2)
library(scaffold)
library(SingleCellExperiment)
select <- dplyr::select
filter <- dplyr::filter
```

We’ll use the `scRNAseq` R package to pull the human pancreas data from
[Baron *et al* 2017](https://doi.org/10.1016/j.cels.2016.08.011) to
serve as our reference dataset. We then simulate a dataset with 3
subjects, 1200 cells allocated equally between subjects, and 10% of the
genes dynamic over the simulated trajectory. We label a gene dynamic at
the population level if (in this case) it’s dynamic in 2/3 subjects.
Internally, the function ensures that 80% of the dynamic genes are
common across subject - a fairly homogeneous trajectory. The underlying
true pseudotime values are stored in the `colData` slot of the
`SingleCellExperiment` object under the name **cell_time_normed**.

``` r
set.seed(312)
panc <- scRNAseq::BaronPancreasData(which = "human")
sim_data <- simulate_multi_subject(ref.dataset = panc, 
                                   perc.dyn.genes = 0.1, 
                                   n.cells = 1200, 
                                   perc.allocation = rep(1/3, 3), 
                                   n.subjects = 3, 
                                   gene.dyn.threshold = 2)
```

The PCA embeddings show us a pretty simple trajectory that’s strongly
correlated with the first principal component.

``` r
plotPCA(sim_data, colour_by = "cell_time_normed") + 
  theme_scLANE()
```

<img src="man/figures/README-plot-sims-pt-1.png" width="100%" />

We also see that the data are not clustered by subject, which indicates
that gene dynamics are mostly homogeneous across subjects.

``` r
plotPCA(sim_data, colour_by = "subject") + 
  theme_scLANE()
```

<img src="man/figures/README-plot-sims-subj-1.png" width="100%" />

## Trajectory DE Testing

Since we have multi-subject data, we can use any of the three model
backends to run our DE testing. We’ll start with the simplest model, the
GLM, then work our way through the other options in order of increasing
complexity. We first prepare our inputs - a dataframe containing our
cell ordering, a set of genes to build models for, and a vector of
per-cell size factors to be used as offsets during estimation. In
reality, it’s usually unnecessary to fit a model for every single gene
in a dataset, as trajectories are usually estimated using a subset of
the entire set of genes (usually a few thousand most highly variable
genes). For the purpose of demonstration, we’ll select 50 genes each
from the dynamic and non-dynamic populations. **Note**: in this case
we’re working with a single pseudotime lineage, though in real datasets
several lineages often exist; in order to fit models for a subset of
lineages simply remove the corresponding columns from the cell ordering
dataframe passed as input to `testDynamic()`.

``` r
set.seed(312)
gene_sample <- c(sample(rownames(sim_data)[rowData(sim_data)$geneStatus_overall == "Dynamic"], size = 50), 
                 sample(rownames(sim_data)[rowData(sim_data)$geneStatus_overall == "NotDynamic"], size = 50))
order_df <- data.frame(X = sim_data$cell_time_normed)
cell_offset <- createCellOffset(sim_data)
```

### GLM Backend

Running `testDynamic()` provides us with a nested list containing model
output & DE test results for each gene over each pseudotime / latent
time lineage. In this case, since we have a true cell ordering we only
have one lineage. Parallel processing is turned on by default, and we
use 2 cores here to speed up runtime.

``` r
scLANE_models_glm <- testDynamic(sim_data, 
                                 pt = order_df, 
                                 genes = gene_sample, 
                                 size.factor.offset = cell_offset, 
                                 n.cores = 2, 
                                 track.time = TRUE)
#> [1] "testDynamic evaluated 100 genes across 1 lineage in 29.86 secs"
```

After the function finishes running, we use `getResultsDE()` to generate
a sorted table of DE test results, with one row for each gene & lineage.
The GLM backend uses a simple likelihood ratio test to compare the null
& alternate models, with the test statistic assumed to be
[asymptotically Chi-squared
distributed](https://en.wikipedia.org/wiki/Likelihood-ratio_test).

``` r
scLANE_res_glm <- getResultsDE(scLANE_models_glm)
select(scLANE_res_glm, Gene, Lineage, Test_Stat, P_Val, P_Val_Adj, Gene_Dynamic_Overall) %>% 
  slice_head(n = 5) %>% 
  knitr::kable(format = "pipe", 
               digits = 3, 
               col.names = c("Gene", "Lineage", "LRT Stat.", "P-value", "Adj. P-value", "Predicted Dynamic Status"))
```

| Gene   | Lineage | LRT Stat. | P-value | Adj. P-value | Predicted Dynamic Status |
|:-------|:--------|----------:|--------:|-------------:|-------------------------:|
| MPG    | A       |   430.567 |       0 |            0 |                        1 |
| WAPAL  | A       |   358.718 |       0 |            0 |                        1 |
| JARID2 | A       |   352.197 |       0 |            0 |                        1 |
| ERGIC3 | A       |   291.170 |       0 |            0 |                        1 |
| PFDN2  | A       |   259.313 |       0 |            0 |                        1 |

After creating a reference table of the ground truth status of each
gene - `1` denotes a dynamic gene and `0` a non-dynamic one - and adding
that binary indicator the the DE results table, we can generate some
classification metrics using [a confusion
matrix](https://en.wikipedia.org/wiki/Confusion_matrix).

``` r
gene_status_df <- data.frame(gene = gene_sample, 
                             True_Gene_Status = ifelse(rowData(sim_data)[gene_sample, ]$geneStatus_overall == "Dynamic", 1, 0))
scLANE_res_glm <- inner_join(scLANE_res_glm, 
                             gene_status_df, 
                             by = c("Gene" = "gene"))
caret::confusionMatrix(factor(scLANE_res_glm$Gene_Dynamic_Overall, levels = c(0, 1)), 
                       factor(scLANE_res_glm$True_Gene_Status, levels = c(0, 1)), 
                       positive = "1")
#> Confusion Matrix and Statistics
#> 
#>           Reference
#> Prediction  0  1
#>          0 50  4
#>          1  0 46
#>                                          
#>                Accuracy : 0.96           
#>                  95% CI : (0.9007, 0.989)
#>     No Information Rate : 0.5            
#>     P-Value [Acc > NIR] : <2e-16         
#>                                          
#>                   Kappa : 0.92           
#>                                          
#>  Mcnemar's Test P-Value : 0.1336         
#>                                          
#>             Sensitivity : 0.9200         
#>             Specificity : 1.0000         
#>          Pos Pred Value : 1.0000         
#>          Neg Pred Value : 0.9259         
#>              Prevalence : 0.5000         
#>          Detection Rate : 0.4600         
#>    Detection Prevalence : 0.4600         
#>       Balanced Accuracy : 0.9600         
#>                                          
#>        'Positive' Class : 1              
#> 
```

### GEE Backend

The function call is essentially the same when using the GLM backend,
with the exception of needing to provide a sorted vector of subject IDs
& a desired correlation structure. We also need to flip the `is.gee`
flag in order to indicate that we’d like to fit estimating equations
models (instead of mixed models). Since fitting GEEs is a fair bit more
computationally complex than fitting GLMs, DE testing with the GEE
backend takes a bit longer. Using more cores and / or running the tests
on an HPC cluster (as this document is being compiled on a 2020 MacBook
Pro) speeds things up considerably. Here we specify an [AR1 correlation
structure](https://rdrr.io/cran/nlme/man/corAR1.html), which is the
default for the GEE backend.

``` r
scLANE_models_gee <- testDynamic(sim_data, 
                                 pt = order_df, 
                                 genes = gene_sample, 
                                 size.factor.offset = cell_offset, 
                                 is.gee = TRUE, 
                                 id.vec = sim_data$subject, 
                                 cor.structure = "ar1", 
                                 n.cores = 2, 
                                 track.time = TRUE)
#> [1] "testDynamic evaluated 100 genes across 1 lineage in 3.671 mins"
```

We again generate the table of DE test results. The variance of the
estimated coefficients is determined using [the sandwich
estimator](https://online.stat.psu.edu/stat504/lesson/12/12.3), and a
Wald test is used to compare the null & alternate models.

``` r
scLANE_res_gee <- getResultsDE(scLANE_models_gee)
select(scLANE_res_gee, Gene, Lineage, Test_Stat, P_Val, P_Val_Adj, Gene_Dynamic_Overall) %>% 
  slice_head(n = 5) %>% 
  knitr::kable("pipe", 
               digits = 3, 
               col.names = c("Gene", "Lineage", "Wald Stat.", "P-value", "Adj. P-value", "Predicted Dynamic Status"))
```

| Gene    | Lineage | Wald Stat. | P-value | Adj. P-value | Predicted Dynamic Status |
|:--------|:--------|-----------:|--------:|-------------:|-------------------------:|
| KLHDC10 | A       |     79.087 |       0 |            0 |                        1 |
| SPCS3   | A       |     86.248 |       0 |            0 |                        1 |
| MRPL20  | A       |    101.487 |       0 |            0 |                        1 |
| TSPAN1  | A       |    111.245 |       0 |            0 |                        1 |
| GGNBP2  | A       |    117.029 |       0 |            0 |                        1 |

We create the same confusion matrix as before. Empirically speaking,
when the underlying dynamics don’t differ much between subjects, GEEs
tend to be more conservative (and thus perform slightly worse) than
GLMs. This is shown below, where the GEE backend has decent accuracy,
but the false negative rate is higher than that of the GLM backend.

``` r
scLANE_res_gee <- inner_join(scLANE_res_gee, 
                             gene_status_df, 
                             by = c("Gene" = "gene"))
caret::confusionMatrix(factor(scLANE_res_gee$Gene_Dynamic_Overall, levels = c(0, 1)), 
                       factor(scLANE_res_gee$True_Gene_Status, levels = c(0, 1)), 
                       positive = "1")
#> Confusion Matrix and Statistics
#> 
#>           Reference
#> Prediction  0  1
#>          0 46 15
#>          1  4 35
#>                                           
#>                Accuracy : 0.81            
#>                  95% CI : (0.7193, 0.8816)
#>     No Information Rate : 0.5             
#>     P-Value [Acc > NIR] : 1.351e-10       
#>                                           
#>                   Kappa : 0.62            
#>                                           
#>  Mcnemar's Test P-Value : 0.02178         
#>                                           
#>             Sensitivity : 0.7000          
#>             Specificity : 0.9200          
#>          Pos Pred Value : 0.8974          
#>          Neg Pred Value : 0.7541          
#>              Prevalence : 0.5000          
#>          Detection Rate : 0.3500          
#>    Detection Prevalence : 0.3900          
#>       Balanced Accuracy : 0.8100          
#>                                           
#>        'Positive' Class : 1               
#> 
```

### GLMM Backend

We re-run the DE tests a final time using the GLMM backend. This is the
most complex model architecture we support, and is the trickiest to
interpret. We recommend using it when you’re most interested in how a
trajectory differs between subjects e.g., if the subjects belong to
groups like Treatment & Control, and you expect the Treatment group to
experience a different progression through the biological process.
Executing the function with the GLMM backend differs only in that we
switch the `is.glmm` flag to `TRUE` and no longer need to specify a
working correlation structure. **Note**: the GLMM backend is still under
development, as we are working on further reducing runtime and
increasing the odds of the underlying optimization process converging
successfully. As such, updates will be frequent and functionality /
results may shift slightly.

``` r
scLANE_models_glmm <- testDynamic(sim_data, 
                                  pt = order_df, 
                                  genes = gene_sample, 
                                  size.factor.offset = cell_offset, 
                                  n.potential.basis.fns = 3, 
                                  is.glmm = TRUE, 
                                  glmm.adaptive = TRUE, 
                                  id.vec = sim_data$subject, 
                                  n.cores = 2, 
                                  track.time = TRUE)
#> [1] "testDynamic evaluated 100 genes across 1 lineage in 3.493 mins"
```

Like the GLM backend, the GLMM backend uses a likelihood ratio test to
compare the null & alternate models. We fit the two nested models using
maximum likelihood estimation instead of
[REML](https://en.wikipedia.org/wiki/Restricted_maximum_likelihood) in
order to perform this test; the null model in this case is a negative
binomial GLMM with a random intercept for each subject.

``` r
scLANE_res_glmm <- getResultsDE(scLANE_models_glmm)
select(scLANE_res_glmm, Gene, Lineage, Test_Stat, P_Val, P_Val_Adj, Gene_Dynamic_Overall) %>% 
  slice_head(n = 5) %>% 
  knitr::kable("pipe", 
               digits = 3, 
               col.names = c("Gene", "Lineage", "LRT Stat.", "P-value", "Adj. P-value", "Predicted Dynamic Status"))
```

| Gene   | Lineage | LRT Stat. | P-value | Adj. P-value | Predicted Dynamic Status |
|:-------|:--------|----------:|--------:|-------------:|-------------------------:|
| WAPAL  | A       |   366.154 |       0 |            0 |                        1 |
| JARID2 | A       |   328.033 |       0 |            0 |                        1 |
| CBX6   | A       |   255.565 |       0 |            0 |                        1 |
| ISOC2  | A       |   225.666 |       0 |            0 |                        1 |
| PFDN2  | A       |   212.453 |       0 |            0 |                        1 |

The GLMM backend performs about as well as the GEE backend. Like with
the GEE backend, it’s more appropriate to use these more complex models
if expression dynamics might differ between subjects, with the
difference being that you should use the GEE backend if you’re
interested in population-level trends & the GLMM backend if you’re
interested in per-subject trends. Since the dynamics in our simulated
data are strongly conserved across subjects, it follows that the simpler
GLMs perform the best.

``` r
scLANE_res_glmm <- inner_join(scLANE_res_glmm, 
                              gene_status_df, 
                              by = c("Gene" = "gene"))
caret::confusionMatrix(factor(scLANE_res_glmm$Gene_Dynamic_Overall, levels = c(0, 1)), 
                       factor(scLANE_res_glmm$True_Gene_Status, levels = c(0, 1)), 
                       positive = "1")
#> Confusion Matrix and Statistics
#> 
#>           Reference
#> Prediction  0  1
#>          0 49 17
#>          1  1 33
#>                                           
#>                Accuracy : 0.82            
#>                  95% CI : (0.7305, 0.8897)
#>     No Information Rate : 0.5             
#>     P-Value [Acc > NIR] : 3.074e-11       
#>                                           
#>                   Kappa : 0.64            
#>                                           
#>  Mcnemar's Test P-Value : 0.000407        
#>                                           
#>             Sensitivity : 0.6600          
#>             Specificity : 0.9800          
#>          Pos Pred Value : 0.9706          
#>          Neg Pred Value : 0.7424          
#>              Prevalence : 0.5000          
#>          Detection Rate : 0.3300          
#>    Detection Prevalence : 0.3400          
#>       Balanced Accuracy : 0.8200          
#>                                           
#>        'Positive' Class : 1               
#> 
```

## Downstream Analysis & Visualization

### Model Comparison

We can use the `plotModels()` to visually compare different types of
modeling backends. It takes as input the results from `testDynamic()`,
as well as a few specifications for which models & lineages should be
plotted. While more complex visualizations can be created from our model
output, this function gives us a good first glance at which models fit
the underlying trend the best. Here we show the output generated using
the GLM backend, split by model type. The intercept-only model shows the
null hypothesis against which the scLANE model is compared using the
likelihood ratio test and the GLM displays the inadequacy of monotonic
modeling architectures for nonlinear dynamics. A GAM shows essentially
the same trend as the `scLANE` model, though the fitted trend from
`scLANE` is smoother & of course more interpretable.

``` r
plotModels(scLANE_models_glm, 
           gene = "JARID2", 
           pt = order_df, 
           expr.mat = sim_data, 
           size.factor.offset = cell_offset)
```

<img src="man/figures/README-plot-models-glm-1.png" width="100%" />

When plotting the models generated using the GLMM backend, we split by
lineage & color the points by subject ID instead of by lineage.

``` r
plotModels(scLANE_models_glmm, 
           gene = "WAPAL", 
           pt = order_df, 
           expr.mat = sim_data, 
           size.factor.offset = cell_offset, 
           is.glmm = TRUE, 
           plot.null = FALSE, 
           id.vec = sim_data$subject)
```

<img src="man/figures/README-plot-models-glmm-1.png" width="100%" />

### Gene Clustering

After generating a suitable set of models, we can cluster the genes in a
semi-supervised fashion using `clusterGenes()`. This function uses the
Leiden algorithm (the default), hierarchical clustering, or *k*-means
and selects the best set of clustering hyperparameters using the
silhouette score. If desired PCA can be run prior to clustering, but the
default is to cluster on the fitted values themselves. We then use the
results as input to `plotClusteredGenes()`, which generates a table of
fitted values per-gene, per-lineage over pseudotime along with the
accompanying cluster labels.

``` r
gene_clusters <- clusterGenes(scLANE_models_glm, 
                              pt = order_df, 
                              size.factor.offset = cell_offset, 
                              clust.algo = "leiden")
gene_clust_table <- plotClusteredGenes(scLANE_models_glm, 
                                       gene.clusters = gene_clusters, 
                                       pt = order_df, 
                                       size.factor.offset = cell_offset, 
                                       n.cores = 2)
slice_sample(gene_clust_table, n = 5) %>% 
  knitr::kable("pipe", 
               digits = 3, 
               row.names = FALSE, 
               col.names = c("Gene", "Lineage", "Cell", "Fitted (link)", "Fitted (response)", "Pseudotime", "Cluster"))
```

| Gene   | Lineage | Cell | Fitted (link) | Fitted (response) | Pseudotime | Cluster |
|:-------|:--------|:-----|--------------:|------------------:|-----------:|:--------|
| MFSD2B | A       | 9    |        -2.823 |             0.302 |      0.022 | 1       |
| DLC1   | A       | 853  |        -2.091 |             0.142 |      0.132 | 1       |
| SMG1   | A       | 244  |        -0.975 |             0.200 |      0.610 | 1       |
| TPT1   | A       | 39   |         4.272 |            35.151 |      0.098 | 2       |
| PCF11  | A       | 441  |        -0.917 |             0.311 |      0.102 | 1       |

The results can then be plotted as desired using `ggplot2` or another
visualization package. Upon visual inspection, the genes seem to cluster
based on whether they are more dynamic or more static over pseudotime.

``` r
ggplot(gene_clust_table, aes(x = PT, y = FITTED, color = CLUSTER, group = GENE)) + 
  facet_wrap(~paste0("Cluster ", CLUSTER)) + 
  geom_line(alpha = 0.75, show.legend = FALSE) + 
  scale_y_continuous(labels = scales::label_number(accuracy = 1)) + 
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) + 
  labs(x = "Pseudotime", 
       y = "Fitted Values", 
       title = "Leiden clustering of gene dynamics") + 
  theme_scLANE()
```

<img src="man/figures/README-plot-clust-1.png" width="100%" />

Checking the true frequency of dynamic genes in each cluster seems to
somewhat confirm that hypothesis:

``` r
distinct(gene_clust_table, GENE, CLUSTER) %>% 
  inner_join(gene_status_df, by = c("GENE" = "gene")) %>% 
  with_groups(CLUSTER, 
              summarise, 
              P_DYNAMIC = mean(True_Gene_Status)) %>% 
  knitr::kable("pipe", 
               digits = 3, 
               col.names = c("Leiden Cluster", "Dynamic Gene Frequency"))
```

| Leiden Cluster | Dynamic Gene Frequency |
|:---------------|-----------------------:|
| 1              |                  0.476 |
| 2              |                  0.968 |

# Conclusions & Best Practices

In general, starting with the GLM backend is probably your best bet
unless you have a strong prior belief that expression trends will differ
significantly between subjects. If that is the case, you should use the
GEE backend if you’re interested in population-level estimates, but are
worried about wrongly predicting differential expression when
differences in expression are actually caused by inter-subject
variation. If you’re interested in generating subject-specific estimates
then the GLMM backend should be used; take care when interpreting the
fixed vs. random effects though, and consult a biostatistician if
necessary.

If you have a large dataset (10,000+ cells), you should start with the
GLM backend, since estimates don’t differ much between modeling methods
given high enough *n*. In addition, running the tests on an HPC cluster
with 4+ CPUs and 64+ GB of RAM will help your computations to complete
swiftly. Datasets with smaller numbers of cells or fewer genes of
interest may be easily analyzed in an R session on a local machine.

# Contact Information

This package is developed & maintained by Jack Leary. Feel free to reach
out by [opening an issue](https://github.com/jr-leary7/scLANE/issues) or
by email (<j.leary@ufl.edu>) if more detailed assistance is needed.

# References

1.  Bacher, R. *et al*. [Enhancing biological signals and detection
    rates in single-cell RNA-seq experiments with cDNA library
    equalization](https://doi.org/10.1093/nar/gkab1071). *Nucleic Acids
    Research* (2021).

2.  Warton, D. & J. Stoklosa. [A generalized estimating equation
    approach to multivariate adaptive regression
    splines](https://doi.org/10.1080/10618600.2017.1360780). *Journal of
    Computational and Graphical Statistics* (2018).

3.  Nelder, J. & R. Wedderburn. [Generalized linear
    models](https://doi.org/10.2307/2344614). *Journal of the Royal
    Statistical Society* (1972).

4.  Liang, K. & S. Zeger. [Longitudinal data analysis using generalized
    linear models](https://doi.org/10.1093/biomet/73.1.13). *Biometrika*
    (1986).

5.  Laird, N. & J. Ware. [Random-effects models for longitudinal
    data](https://doi.org/10.2307/2529876). *Biometrics* (1988).
