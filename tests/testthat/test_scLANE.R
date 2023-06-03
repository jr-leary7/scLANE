# load data & prepare for testing
load(system.file("testdata/sim_test_data.RData", package = "scLANE"))
genes_to_test <- c(rownames(sim_data)[SummarizedExperiment::rowData(sim_data)$geneStatus_overall == "Dynamic"][1:5],
                   rownames(sim_data)[SummarizedExperiment::rowData(sim_data)$geneStatus_overall == "NotDynamic"][1:5])
counts_test <- t(as.matrix(SingleCellExperiment::counts(sim_data)[genes_to_test, ]))
pt_test <- data.frame(PT = sim_data$cell_time_normed)
sink(tempfile())
# run GLM, GEE, & GLMM tests
glm_gene_stats <- testDynamic(expr.mat = counts_test,
                              pt = pt_test,
                              n.potential.basis.fns = 3,
                              is.gee = FALSE,
                              parallel.exec = TRUE,
                              n.cores = 2,
                              approx.knot = TRUE,
                              track.time = TRUE)
gee_gene_stats <- testDynamic(expr.mat = counts_test,
                              pt = pt_test,
                              n.potential.basis.fns = 3,
                              is.gee = TRUE,
                              id.vec = sim_data$subject,
                              parallel.exec = TRUE,
                              n.cores = 2,
                              approx.knot = TRUE,
                              track.time = TRUE)
glmm_gene_stats <- testDynamic(expr.mat = counts_test,
                               pt = pt_test,
                               n.potential.basis.fns = 3,
                               is.glmm = TRUE,
                               id.vec = sim_data$subject,
                               parallel.exec = TRUE,
                               n.cores = 2,
                               approx.knot = TRUE,
                               track.time = TRUE)
# get results tables overall
glm_test_results <- getResultsDE(glm_gene_stats)
gee_test_results <- getResultsDE(gee_gene_stats)
glmm_test_results <- getResultsDE(glmm_gene_stats)
# get results tables by interval
glm_slope_test <- testSlope(glm_gene_stats)
gee_slope_test <- testSlope(gee_gene_stats)
glmm_slope_test <- testSlope(glmm_gene_stats)
# run NB GAM
gam_mod_ps <- nbGAM(expr = counts_test[, 1],
                    pt = pt_test,
                    use.b.spline = FALSE)
gam_mod_bs <- nbGAM(expr = counts_test[, 1],
                    pt = pt_test,
                    use.b.spline = TRUE)
# run original marge model
marge_mod <- marge2(X_pred = pt_test,
                    Y = counts_test[, 3],
                    M = 3,
                    return.basis = TRUE,
                    return.GCV = TRUE,
                    return.wic = TRUE)
# run NB GLMM
glmm_mod <- fitGLMM(X_pred = pt_test,
                    Y = counts_test[, 1],
                    id.vec = sim_data$subject,
                    adaptive = FALSE,
                    approx.knot = TRUE,
                    M.glm = 2,
                    return.basis = TRUE,
                    return.GCV = TRUE)
# generate plots
plot_glm <- plotModels(glm_gene_stats,
                       gene = "ABR",
                       pt = pt_test,
                       gene.counts = counts_test)
plot_gee <- plotModels(glm_gene_stats,
                       gene = "ABR",
                       pt = pt_test,
                       gene.counts = counts_test,
                       is.gee = TRUE,
                       id.vec = sim_data$subject,
                       cor.structure = "exchangeable")
# downstream analysis
set.seed(312)
gene_clusters <- clusterGenes(glm_gene_stats, clust.algo = "leiden")
gene_clust_table <- plotClusteredGenes(glm_gene_stats,
                                       gene.clusters = gene_clusters,
                                       pt = pt_test,
                                       n.cores = 2)
smoothed_counts <- smoothedCountsMatrix(test.dyn.res = glm_gene_stats,
                                        parallel.exec = TRUE,
                                        n.cores = 2)
fitted_values_table <- getFittedValues(glm_gene_stats,
                                       genes = names(glm_gene_stats),
                                       pt = pt_test,
                                       expr.mat = counts_test,
                                       cell.meta.data = as.data.frame(SummarizedExperiment::colData(sim_data)),
                                       id.vec = sim_data$subject)
gsea_res <- enrichDynamicGenes(glm_test_results, gene.set.cat = "C2")
sink()

test_that("testDynamic() output", {
  expect_s3_class(glm_gene_stats, "scLANE")
  expect_s3_class(gee_gene_stats, "scLANE")
  expect_s3_class(glmm_gene_stats, "scLANE")
  expect_equal(length(glm_gene_stats), 10)
  expect_equal(length(gee_gene_stats), 10)
  expect_equal(length(glmm_gene_stats), 10)
  expect_equal(nrow(glm_gene_stats$ABCF1$Lineage_A$MARGE_Summary), 2)
  expect_equal(ncol(glm_gene_stats$ABCF1$Lineage_A$MARGE_Summary), 5)
  expect_equal(nrow(gee_gene_stats$ABCF1$Lineage_A$MARGE_Summary), 2)
  expect_equal(ncol(gee_gene_stats$ABCF1$Lineage_A$MARGE_Summary), 5)
  expect_equal(nrow(glmm_gene_stats$ABCF1$Lineage_A$MARGE_Summary), 5)
  expect_equal(ncol(glmm_gene_stats$ABCF1$Lineage_A$MARGE_Summary), 5)
  expect_equal(sum(unname(sapply(glm_gene_stats, \(x) x$Lineage_A$Model_Status)) == "MARGE model OK, null model OK"), 10)
  expect_equal(sum(unname(sapply(gee_gene_stats, \(x) x$Lineage_A$Model_Status)) == "MARGE model OK, null model OK"), 10)
  expect_equal(sum(unname(sapply(glmm_gene_stats, \(x) x$Lineage_A$Model_Status)) == "MARGE model OK, null model OK"), 8)
  expect_equal(unique(sapply(glm_gene_stats, \(x) length(x$Lineage_A))), 17)
  expect_equal(unique(sapply(gee_gene_stats, \(x) length(x$Lineage_A))), 17)
  expect_equal(unique(sapply(glmm_gene_stats, \(x) length(x$Lineage_A))), 17)
})

test_that("getResultsDE() output", {
  expect_s3_class(glm_test_results, "data.frame")
  expect_s3_class(gee_test_results, "data.frame")
  expect_s3_class(glmm_test_results, "data.frame")
  expect_equal(nrow(glm_test_results), 10)
  expect_equal(nrow(gee_test_results), 10)
  expect_equal(nrow(glmm_test_results), 10)
  expect_equal(ncol(glm_test_results), 14)
  expect_equal(ncol(gee_test_results), 15)
  expect_equal(ncol(glmm_test_results), 15)
  expect_gt(sum(glm_test_results$Gene_Dynamic_Overall), 0)
  expect_gt(sum(gee_test_results$Gene_Dynamic_Overall), 0)
  expect_gt(sum(glmm_test_results$Gene_Dynamic_Overall), 0)
})

test_that("testSlope() output", {
  expect_s3_class(glm_slope_test, "data.frame")
  expect_s3_class(gee_slope_test, "data.frame")
  expect_s3_class(glmm_slope_test, "data.frame")
  expect_equal(nrow(glm_slope_test), 10)
  expect_equal(nrow(gee_slope_test), 10)
  expect_equal(nrow(glmm_slope_test), 34)
  expect_equal(ncol(glm_slope_test), 11)
  expect_equal(ncol(gee_slope_test), 11)
  expect_equal(ncol(glmm_slope_test), 11)
  expect_equal(sum(glm_slope_test$Notes == "MARGE model error, null model OK", na.rm = TRUE), 0)
  expect_equal(sum(gee_slope_test$Notes == "MARGE model error, null model OK", na.rm = TRUE), 0)
  expect_equal(sum(glmm_slope_test$Notes == "MARGE model error, null model OK", na.rm = TRUE), 2)
  expect_equal(sum(glm_slope_test$P_Val_Adj < 0.01, na.rm = TRUE), 5)
  expect_equal(sum(gee_slope_test$P_Val_Adj < 0.01, na.rm = TRUE), 4)
  expect_equal(sum(glmm_slope_test$P_Val_Adj < 0.01, na.rm = TRUE), 0)
})

test_that("nbGAM() output", {
  expect_s3_class(gam_mod_bs, "gamlss")
  expect_s3_class(gam_mod_ps, "gamlss")
  expect_equal(length(coef(gam_mod_bs)), 4)
  expect_equal(length(coef(gam_mod_ps)), 2)
  expect_true(gam_mod_bs$converged)
  expect_true(gam_mod_ps$converged)
})

test_that("marge2() output", {
  expect_s3_class(marge_mod, "marge")
  expect_s3_class(marge_mod$final_mod, "negbin")
  expect_equal(length(coef(marge_mod$final_mod)), 2)
  expect_equal(ncol(marge_mod$bx), 2)
  expect_equal(round(marge_mod$final_mod$theta, 1), 3.8)
  expect_equal(marge_mod$model_type, "GLM")
  expect_true(marge_mod$final_mod$converged)
})
test_that("fitGLMM() output", {
  expect_s3_class(glmm_mod$final_mod, "glmmTMB")
  expect_equal(nrow(coef(glmm_mod$final_mod)$cond$subject), 2)
  expect_equal(ncol(coef(glmm_mod$final_mod)$cond$subject), 5)
  expect_false(glmm_mod$final_mod$modelInfo$REML)
  expect_equal(length(fitted(glmm_mod$final_mod)), 300)
})

test_that("plotModels() output", {
  expect_s3_class(plot_glm, "ggplot")
  expect_s3_class(plot_gee, "ggplot")
  expect_equal(ncol(plot_glm$data), 11)
  expect_equal(ncol(plot_gee$data), 11)
  expect_equal(nrow(plot_glm$data), 1200)
  expect_equal(nrow(plot_gee$data), 1200)
})

test_that("clusterGenes() output", {
  expect_s3_class(gene_clusters, "data.frame")
  expect_equal(ncol(gene_clusters), 3)
  expect_equal(nrow(gene_clusters), 10)
})

test_that("plotClusteredGenes() output", {
  expect_s3_class(gene_clust_table, "data.frame")
  expect_equal(ncol(gene_clust_table), 7)
  expect_equal(nrow(gene_clust_table), 3000)
})

test_that("smoothedCountsMatrix() output", {
  expect_type(smoothed_counts, "list")
  expect_equal(length(smoothed_counts), 1)
  expect_type(smoothed_counts$Lineage_A, "double")
  expect_equal(ncol(smoothed_counts$Lineage_A), 10)
  expect_equal(nrow(smoothed_counts$Lineage_A), 300)
})

test_that("getFittedValues() output", {
  expect_s3_class(fitted_values_table, "data.frame")
  expect_equal(ncol(fitted_values_table), 17)
  expect_equal(nrow(fitted_values_table), 3000)
})

test_that("enrichDynamicGenes() output", {
  expect_s4_class(gsea_res, "enrichResult")
  expect_s3_class(gsea_res@result, "data.frame")
  expect_equal(ncol(gsea_res@result), 9)
})
