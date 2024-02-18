
- <a href="#sclane" id="toc-sclane"><code>scLANE</code></a>
  - <a href="#installation" id="toc-installation">Installation</a>
  - <a href="#model-structure" id="toc-model-structure">Model structure</a>
- <a href="#usage" id="toc-usage">Usage</a>
  - <a href="#libraries" id="toc-libraries">Libraries</a>
  - <a href="#input-data" id="toc-input-data">Input data</a>
  - <a href="#trajectory-de-testing"
    id="toc-trajectory-de-testing">Trajectory DE testing</a>
  - <a href="#downstream-analysis--visualization"
    id="toc-downstream-analysis--visualization">Downstream analysis &amp;
    visualization</a>
- <a href="#conclusions--best-practices"
  id="toc-conclusions--best-practices">Conclusions &amp; best
  practices</a>
- <a href="#contact-information" id="toc-contact-information">Contact
  information</a>
- <a href="#references" id="toc-references">References</a>

<!-- README.md is generated from README.Rmd. Please edit that file -->

# `scLANE`

<!-- badges: start -->

[![R-CMD-check](https://github.com/jr-leary7/scLANE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jr-leary7/scLANE/actions/workflows/R-CMD-check.yaml)
[![Bioc-check](https://github.com/jr-leary7/scLANE/actions/workflows/bioc-check.yaml/badge.svg)](https://github.com/jr-leary7/scLANE/actions/workflows/bioc-check.yaml)
![release](https://img.shields.io/github/v/release/jr-leary7/scLANE?color=purple)
![last
commit](https://img.shields.io/github/last-commit/jr-leary7/scLANE/main?color=darkgreen)
[![codecov](https://codecov.io/gh/jr-leary7/scLANE/branch/main/graph/badge.svg?token=U2U5RTF2VW)](https://codecov.io/gh/jr-leary7/scLANE)
[![CodeFactor](https://www.codefactor.io/repository/github/jr-leary7/sclane/badge)](https://www.codefactor.io/repository/github/jr-leary7/sclane)
[![DOI](https://img.shields.io/static/v1?label=DOI&message=10.1101/2023.12.19.572477&color=blue)](https://doi.org/10.1101/2023.12.19.572477)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

## Installation

You can install the most recent version of `scLANE` using:

``` r
remotes::install_github("jr-leary7/scLANE")
```

## Model structure

The `scLANE` package enables users to accurately determine differential
expression of genes over pseudotime or latent time, and to characterize
gene’s dynamics using interpretable model coefficients. `scLANE` builds
upon the `marge` modeling
framework([GitHub](https://github.com/JakubStats/marge),
[paper](https://doi.org/10.1080/10618600.2017.1360780)), allowing users
to characterize their trajectory’s effects on gene expression using
negative binomial GLMs, GEEs, or GLMMs depending on the experimental
design & biological questions of interest. This modeling framework is an
extension of the [Multivariate Adapative Regression Splines
(MARS)](https://en.wikipedia.org/wiki/Multivariate_adaptive_regression_spline)
method, which builds nonlinear models out of piecewise linear
components. `scLANE` is agnostic with respect to the ordering estimation
method used, and can be implemented downstream of any pseudotime or RNA
velocity method.

A quickstart guide on how to use `scLANE` with simulated data continues
below, and a more detailed vignette showcasing its performance on real
data can be found
[here](https://jr-leary7.github.io/quarto-site/tutorials/scLANE_Trajectory_DE.html).

# Usage

Our method relies on a relatively simple test in order to define whether
a given gene is differentially expressed (or “dynamic”) over the
provided trajectory. While the exact structure of the test differs by
model mode, the concept is the same: the spline-based NB GLM / GEE /
GLMM is treated as the alternate model, and a null model is fit using
the corresponding model mode. If the GLM mode is used, then the null
model is simply an intercept-only NB GLM; the GEE mode fits an
intercept-only model with the same working correlation structure as the
alternate model, and if the GLMM mode is used then the null model is an
intercept-only model with random intercepts for each subject. The
alternate hypothesis is that at least one of the estimated coefficients
is significantly different from zero. We predict a given gene to be
dynamic if the adjusted *p*-value of the test is less than the default
$\alpha = 0.01$ threshold, and classify it as static otherwise.

## Libraries

``` r
library(dplyr)
library(scLANE)
library(ggplot2)
```

## Input data

We read a previously-simulated dataset comprised of cells from 3
subjects exhibiting a homogeneous trajectory structure from [the Zenodo
repository](https://doi.org/10.5281/zenodo.8433077). The underlying true
pseudotime values are stored in the `colData` slot of the
`SingleCellExperiment` object under the name **cell_time_normed**.

``` r
sim_data <- readRDS(url("https://zenodo.org/records/8433077/files/scLANE_sim_data.Rds"))
```

The PCA embeddings show us a pretty simple trajectory that’s strongly
correlated with the first principal component.

``` r
data.frame(sim_data@int_colData$reducedDims@listData$PCA[, 1:2]) %>% 
  mutate(pseudotime = sim_data$cell_time_normed) %>% 
  ggplot(aes(x = PC1, y = PC2, color = pseudotime)) + 
  geom_point(size = 2, alpha = 0.75, stroke = 0) + 
  scale_color_gradientn(colors = viridisLite::plasma(n = 20)) + 
  labs(x = "PC 1", y = "PC 2", color = "Pseudotime") + 
  theme_scLANE(umap = TRUE)
```

<img src="man/figures/README-plot-sims-pt-1.png" width="100%" />

We also see that the data are not clustered by subject, which indicates
that gene dynamics are mostly homogeneous across subjects.

``` r
data.frame(sim_data@int_colData$reducedDims@listData$PCA[, 1:2]) %>% 
  mutate(subject = sim_data$subject) %>% 
  ggplot(aes(x = PC1, y = PC2, color = subject)) + 
  geom_point(size = 2, alpha = 0.75, stroke = 0) + 
  labs(x = "PC 1", y = "PC 2", color = "Subject ID") + 
  theme_scLANE(umap = TRUE)
```

<img src="man/figures/README-plot-sims-subj-1.png" width="100%" />

## Trajectory DE testing

Since we have multi-subject data, we can use any of the three model
modes to run our DE testing. We’ll start with the simplest model, the
GLM, then work our way through the other options in order of increasing
complexity. We first prepare our inputs - a dataframe containing our
cell ordering, a set of genes to build models for, and a vector of
per-cell size factors to be used as offsets during estimation. In
reality, it’s usually unnecessary to fit a model for every single gene
in a dataset, as trajectories are usually estimated using a subset of
the entire set of genes (usually a few thousand most highly variable
genes). For the purpose of demonstration, we’ll select 50 genes each
from the dynamic and non-dynamic populations.

**Note:** In this case we’re working with a single pseudotime lineage,
though in real datasets several lineages often exist; in order to fit
models for a subset of lineages simply remove the corresponding columns
from the cell ordering dataframe passed as input to `testDynamic()`.

``` r
set.seed(312)
gene_sample <- c(sample(rownames(sim_data)[rowData(sim_data)$geneStatus_overall == "Dynamic"], size = 50), 
                 sample(rownames(sim_data)[rowData(sim_data)$geneStatus_overall == "NotDynamic"], size = 50))
order_df <- data.frame(X = sim_data$cell_time_normed)
cell_offset <- createCellOffset(sim_data)
```

### GLM mode

Running `testDynamic()` provides us with a nested list containing model
output & DE test results for each gene over each pseudotime / latent
time lineage. In this case, since we have a true cell ordering we only
have one lineage. Parallel processing is turned on by default, and we
use 4 cores here to speed up runtime.

``` r
scLANE_models_glm <- testDynamic(sim_data, 
                                 pt = order_df, 
                                 genes = gene_sample, 
                                 size.factor.offset = cell_offset, 
                                 n.cores = 4L, 
                                 verbose = FALSE)
#> Registered S3 method overwritten by 'bit':
#>   method   from  
#>   print.ri gamlss
#> scLANE testing completed for 100 genes across 1 lineage in 35.37 secs
```

After the function finishes running, we use `getResultsDE()` to generate
a sorted table of DE test results, with one row for each gene & lineage.
The GLM mode uses a simple likelihood ratio test to compare the null &
alternate models, with the test statistic assumed to be asymptotically
Chi-squared distributed.

``` r
scLANE_res_glm <- getResultsDE(scLANE_models_glm)
select(scLANE_res_glm, Gene, Lineage, Test_Stat, P_Val, P_Val_Adj, Gene_Dynamic_Overall) %>% 
  slice_sample(n = 5) %>% 
  knitr::kable(format = "pipe", 
               digits = 3, 
               col.names = c("Gene", "Lineage", "LRT stat.", "P-value", "Adj. p-value", "Predicted dynamic status"))
```

| Gene       | Lineage | LRT stat. | P-value | Adj. p-value | Predicted dynamic status |
|:-----------|:--------|----------:|--------:|-------------:|-------------------------:|
| RAB1B      | A       |   219.950 |   0.000 |        0.000 |                        1 |
| LY75.CD302 | A       |     4.858 |   0.028 |        0.541 |                        0 |
| UAP1L1     | A       |     9.894 |   0.007 |        0.163 |                        0 |
| TMCO3      | A       |   167.311 |   0.000 |        0.000 |                        1 |
| GOLGA8EP   | A       |     4.201 |   0.040 |        0.567 |                        0 |

### GEE mode

The function call is essentially the same when using the GLM mode, with
the exception of needing to provide a sorted vector of subject IDs & a
desired correlation structure. We also need to flip the `is.gee` flag in
order to indicate that we’d like to fit estimating equations models
(instead of mixed models). Since fitting GEEs is more computationally
complex than fitting GLMs, DE testing with the GEE mode takes a bit
longer. Using more cores and / or running the tests on an HPC cluster
speeds things up considerably.

``` r
scLANE_models_gee <- testDynamic(sim_data, 
                                 pt = order_df, 
                                 genes = gene_sample, 
                                 size.factor.offset = cell_offset, 
                                 is.gee = TRUE, 
                                 id.vec = sim_data$subject, 
                                 cor.structure = "ar1", 
                                 n.cores = 4L, 
                                 verbose = FALSE)
#> scLANE testing completed for 100 genes across 1 lineage in 1.525 mins
```

We again generate the table of DE test results. The variance of the
estimated coefficients is determined using the sandwich estimator, and a
Wald test is used to compare the null & alternate models.

``` r
scLANE_res_gee <- getResultsDE(scLANE_models_gee)
select(scLANE_res_gee, Gene, Lineage, Test_Stat, P_Val, P_Val_Adj, Gene_Dynamic_Overall) %>% 
  slice_sample(n = 5) %>% 
  knitr::kable("pipe", 
               digits = 3, 
               col.names = c("Gene", "Lineage", "Wald stat.", "P-value", "Adj. p-value", "Predicted dynamic status"))
```

| Gene     | Lineage | Wald stat. | P-value | Adj. p-value | Predicted dynamic status |
|:---------|:--------|-----------:|--------:|-------------:|-------------------------:|
| DGUOK    | A       | 200675.460 |       0 |            0 |                        1 |
| TBCC     | A       |     40.399 |       0 |            0 |                        1 |
| GOLGA8EP | A       |         NA |      NA |           NA |                        0 |
| EMC3     | A       |   8397.337 |       0 |            0 |                        1 |
| DDX41    | A       |   3486.998 |       0 |            0 |                        1 |

### GLMM mode

We re-run the DE tests a final time using the GLMM mode. This is the
most complex model architecture we support, and is the trickiest to
interpret. We recommend using it when you’re most interested in how a
trajectory differs between subjects e.g., if the subjects belong to
groups like Treatment & Control, and you expect the Treatment group to
experience a different progression through the biological process.
Executing the function with the GLMM mode differs only in that we switch
the `is.glmm` flag to `TRUE` and no longer need to specify a working
correlation structure.

``` r
scLANE_models_glmm <- testDynamic(sim_data, 
                                  pt = order_df, 
                                  genes = gene_sample, 
                                  size.factor.offset = cell_offset, 
                                  n.potential.basis.fns = 3, 
                                  is.glmm = TRUE, 
                                  id.vec = sim_data$subject, 
                                  n.cores = 4L, 
                                  verbose = FALSE)
#> scLANE testing completed for 100 genes across 1 lineage in 3.133 mins
```

**Note:** The GLMM mode is still under development, as we are working on
further reducing runtime and increasing the odds of the underlying
optimization process converging successfully. As such, updates will be
frequent and functionality / results may shift slightly.

Like the GLM mode, the GLMM mode uses a likelihood ratio test to compare
the null & alternate models.

``` r
scLANE_res_glmm <- getResultsDE(scLANE_models_glmm)
select(scLANE_res_glmm, Gene, Lineage, Test_Stat, P_Val, P_Val_Adj, Gene_Dynamic_Overall) %>% 
  slice_sample(n = 5) %>% 
  knitr::kable("pipe", 
               digits = 3, 
               col.names = c("Gene", "Lineage", "LRT stat.", "P-value", "Adj. p-value", "Predicted dynamic status"))
```

| Gene   | Lineage | LRT stat. | P-value | Adj. p-value | Predicted dynamic status |
|:-------|:--------|----------:|--------:|-------------:|-------------------------:|
| VDAC1  | A       |   205.644 |   0.000 |            0 |                        1 |
| CKAP4  | A       |   139.246 |   0.000 |            0 |                        1 |
| MRTO4  | A       |     5.048 |   0.998 |            1 |                        0 |
| FOXD3  | A       |     2.586 |   1.000 |            1 |                        0 |
| UBE2Q1 | A       |    94.077 |   0.000 |            0 |                        1 |

## Downstream analysis & visualization

### Model comparison

We can use the `plotModels()` to visually compare different types of
models. It takes as input the results from `testDynamic()`, as well as a
few specifications for which models & lineages should be plotted. While
more complex visualizations can be created from our model output, this
function gives us a good first glance at which models fit the underlying
trend the best. Here we show the output generated using the GLM mode,
split by model type. The intercept-only model shows the null hypothesis
against which the scLANE model is compared using the likelihood ratio
test and the GLM displays the inadequacy of monotonic modeling
architectures for nonlinear dynamics. A GAM shows essentially the same
trend as the `scLANE` model, though the fitted trend from `scLANE` is
more interpretable.

``` r
plotModels(scLANE_models_glm, 
           gene = scLANE_res_glm$Gene[1], 
           pt = order_df, 
           expr.mat = sim_data, 
           size.factor.offset = cell_offset, 
           plot.null = TRUE, 
           plot.glm = TRUE, 
           plot.gam = TRUE, 
           plot.scLANE = TRUE)
```

<img src="man/figures/README-plot-models-glm-1.png" width="100%" />

When plotting the models generated using the GLMM mode, we split by
lineage & color the points by subject ID instead of by lineage. The gene
in question highlights the utility of the scLANE model, since the gene
dynamics differ significantly by subject.

``` r
plotModels(scLANE_models_glmm, 
           gene = scLANE_res_glmm$Gene[1], 
           pt = order_df, 
           expr.mat = sim_data, 
           size.factor.offset = cell_offset, 
           id.vec = sim_data$subject, 
           is.glmm = TRUE, 
           plot.glm = TRUE,
           plot.gam = TRUE, 
           plot.scLANE = TRUE)
```

<img src="man/figures/README-plot-models-glmm-1.png" width="100%" />

### Coefficient summaries

A key feature of `scLANE` is the ability to obtain a quantitative,
interpretable coefficient for the effect of pseudotime on gene
expression. This functionality is currently available for the GLM & GEE
frameworks, and each coefficient carries the [interpretation of a
generalized linear
model](https://stats.oarc.ucla.edu/r/dae/negative-binomial-regression/).

``` r
scLANE_models_glm[["JARID2"]]$Lineage_A$Gene_Dynamics %>% 
  knitr::kable("pipe", 
               digits = 2, 
               col.names = c("Gene", "Lineage", "Breakpoint", "First Slope", "Second Slope", "First Trend", "Second Trend"))
```

| Gene   | Lineage | Breakpoint | First Slope | Second Slope | First Trend | Second Trend |
|:-------|:--------|-----------:|------------:|-------------:|------------:|-------------:|
| JARID2 | A       |        0.1 |      -38.09 |         3.93 |          -1 |            1 |

Coefficients can also be plotted like so:

``` r
plotModelCoefs(scLANE_models_glm, 
               gene = "JARID2", 
               pt = order_df, 
               expr.mat = sim_data,
               size.factor.offset = cell_offset)
```

<img src="man/figures/README-plot-model-coefs-1.png" width="100%" />

### Knot distribution

Lastly, we can pull the locations in pseudotime of all the knots fitted
by `scLANE`. Visualizing this distribution gives us some idea of where
transcriptional switches are occurring in the set of genes classified as
dynamic.

``` r
dyn_genes <- filter(scLANE_res_glm, Gene_Dynamic_Overall == 1) %>% 
             pull(Gene)
knot_dist <- getKnotDist(scLANE_models_glm, dyn.genes = dyn_genes)
ggplot(knot_dist, aes(x = knot)) + 
  geom_histogram(aes(y = after_stat(density)), 
                 color = "black", 
                 fill = "white", 
                 linewidth = 0.5) + 
  geom_density(color = "forestgreen", 
               fill = "forestgreen", 
               alpha = 0.5, 
               linewidth = 0.75) + 
  labs(x = "Knot Location", y = "Density") + 
  theme_scLANE()
#> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

<img src="man/figures/README-plot-knot-dist-1.png" width="100%" />

### Smoothed dynamics matrix

We can extract matrix of the fitted values for each dynamic gene using
the `smoothedCountsMatrix()` function.

``` r
smoothed_dynamics <- smoothedCountsMatrix(scLANE_models_glm, 
                                          size.factor.offset = cell_offset, 
                                          pt = order_df, 
                                          genes = dyn_genes)
```

The smoothed dynamics can then be used to generate expression cascade
heatmaps, cluster genes, etc. For more information on downstream
analysis of gene dynamics, see [the corresponding
vignette](https://jr-leary7.github.io/quarto-site/tutorials/scLANE_Trajectory_DE.html#downstream-analysis).

# Conclusions & best practices

In general, starting with the GLM mode is probably your best bet unless
you have a strong prior belief that expression trends will differ
significantly between subjects. If that is the case, you should use the
GEE mode if you’re interested in population-level estimates, but are
worried about wrongly predicting differential expression when
differences in expression are actually caused by inter-subject
variation. If you’re interested in generating subject-specific estimates
then the GLMM mode should be used; take care when interpreting the fixed
vs. random effects though, and consult a biostatistician if necessary.

If you have a large dataset (10,000+ cells), you should start with the
GLM mode, since standard error estimates don’t differ much between
modeling methods given high enough *n*. In addition, running the tests
on an HPC cluster with 4+ CPUs and 64+ GB of RAM will help your
computations to complete swiftly. Datasets with smaller numbers of cells
or fewer genes of interest may be easily analyzed in an R session on a
local machine.

# Contact information

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
