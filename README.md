[![R-CMD-check](https://github.com/jr-leary7/scLANE/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main)](https://github.com/jr-leary7/scLANE/actions/workflows/R-CMD-check.yaml)

# scLANE

The `scLANE` package enables users to accurately determine differential expression of genes over pseudotime or latent time, and to characterize gene's dynamics using interpretable model coefficients. `scLANE` builds upon the [`marge` modeling framework](https://github.com/JakubStats/marge) developed by Stoklosa & Warton, and applies `marge`'s negative binomial GLM spline models to scRNA-seq data. 

# Installation 

You can install `scLANE` (once the repository becomes public) as follows:

```
remotes::install_github("jr-leary7/scLANE")
```

# Identifying Differentially Expressed Genes

Say we have a matrix of integer counts (non-normalized) named `expr_mat` and a dataframe of pseudotime values named `pt_df`. We can then use those as inputs to the main test function in `scLANE`, `testDynamic()`, in order to determine which genes are differentially expressed over pseudotime. If the number of cells and/or genes is high, we'll probably want to use parallel processing to speed up the function's execution. 

```
gene_stats <- testDynamic(expr.mat = expr_mat, 
                          pt = pt_df, 
                          parallel = TRUE, 
                          n.cores = 5)
```

After running `testDynamic()`, which returns a list containing model statistics, test results, & fitted values, we'll use `getResultsDE()` to format that list into a tidy table of results at the per-gene level. 

```
test_results <- getResultsDE(gene_stats)
```

Lastly, fitted values from `marge`, intercept-only, GLM, and GAM models can be plotted over gene expression & pseudotime using `plotModels()`. In this case we'd be plotting the results for the gene AURKA. 

```
plotModels(test.dyn.res = gene_stats, 
           gene = "AURKA", 
           pt = pt_df, 
           gene.counts = sim_counts)
```

# Contact Information 

This package is developed & maintained by Jack Leary. Feel free to reach out by opening an issue or by email (j.leary@ufl.edu) if more detailed assistance is needed. 
