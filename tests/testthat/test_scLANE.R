# load data & prepare for testing
load(system.file("testdata/sim_test_data.RData", package = "scLANE"))
cell_offset <- createCellOffset(expr.mat = sim_data)
genes_to_test <- c(rownames(sim_data)[SummarizedExperiment::rowData(sim_data)$geneStatus_overall == "Dynamic"][1:5],
                   rownames(sim_data)[SummarizedExperiment::rowData(sim_data)$geneStatus_overall == "NotDynamic"][1:5])
counts_test <- t(as.matrix(SingleCellExperiment::counts(sim_data)[genes_to_test, ]))
pt_test <- data.frame(PT = sim_data$cell_time_normed)

# generate results w/ all functions of interest
withr::with_output_sink(tempfile(), {
  # run GLM, GEE, & GLMM tests
  glm_gene_stats <- testDynamic(expr.mat = counts_test,
                                pt = pt_test,
                                n.potential.basis.fns = 5,
                                size.factor.offset = cell_offset,
                                n.cores = 2,
                                track.time = TRUE)
  gee_gene_stats <- testDynamic(expr.mat = counts_test,
                                pt = pt_test,
                                n.potential.basis.fns = 5,
                                size.factor.offset = cell_offset,
                                is.gee = TRUE,
                                cor.structure = "ar1",
                                id.vec = sim_data$subject,
                                n.cores = 2,
                                track.time = TRUE)
  glmm_gene_stats <- testDynamic(expr.mat = counts_test,
                                 pt = pt_test,
                                 size.factor.offset = cell_offset,
                                 n.potential.basis.fns = 3,
                                 is.glmm = TRUE,
                                 glmm.adaptive = TRUE,
                                 id.vec = sim_data$subject,
                                 n.cores = 2,
                                 track.time = TRUE)
  # get results tables overall
  glm_test_results <- getResultsDE(glm_gene_stats)
  gee_test_results <- getResultsDE(gee_gene_stats)
  glmm_test_results <- getResultsDE(glmm_gene_stats)
  # get results tables by interval
  glm_slope_test <- testSlope(glm_gene_stats)
  gee_slope_test <- testSlope(gee_gene_stats)
  glmm_slope_test <- testSlope(glmm_gene_stats)
  # run NB GAMs of varying structure
  gam_mod_bs <- nbGAM(expr = counts_test[, 1],
                      pt = pt_test,
                      Y.offset = cell_offset)
  gam_mod_ps <- nbGAM(expr = counts_test[, 1],
                      pt = pt_test,
                      Y.offset = cell_offset,
                      penalize.spline = TRUE)
  gam_mod_ps_mix <- nbGAM(expr = counts_test[, 1],
                          pt = pt_test,
                          Y.offset = cell_offset,
                          id.vec = sim_data$subject,
                          penalize.spline = TRUE)
  # run GLM model -- no offset
  marge_mod <- marge2(X_pred = pt_test,
                      Y = counts_test[, 3],
                      M = 5,
                      return.basis = TRUE,
                      return.GCV = TRUE,
                      return.WIC = TRUE)
  # run GLM model -- with offset
  marge_mod_offset <- marge2(X_pred = pt_test,
                             Y = counts_test[, 3],
                             Y.offset = cell_offset,
                             M = 5,
                             return.basis = TRUE,
                             return.GCV = TRUE,
                             return.WIC = TRUE)
  # run GEE model -- no offset
  marge_mod_GEE <- marge2(X_pred = pt_test,
                          Y = counts_test[, 3],
                          M = 5,
                          is.gee = TRUE,
                          id.vec = sim_data$subject,
                          cor.structure = "ar1",
                          return.basis = TRUE,
                          return.GCV = TRUE,
                          return.WIC = TRUE)
  # run GEE model -- with offset
  marge_mod_GEE_offset <- marge2(X_pred = pt_test,
                                 Y = counts_test[, 3],
                                 Y.offset = cell_offset,
                                 M = 5,
                                 is.gee = TRUE,
                                 id.vec = sim_data$subject,
                                 cor.structure = "ar1",
                                 return.basis = TRUE,
                                 return.GCV = TRUE,
                                 return.WIC = TRUE)
  # run GLMM model -- no offset
  glmm_mod <- fitGLMM(X_pred = pt_test,
                      Y = counts_test[, 3],
                      id.vec = sim_data$subject,
                      adaptive = TRUE,
                      M.glm = 3,
                      return.basis = TRUE,
                      return.GCV = TRUE)
  # run GLMM model -- with offset
  glmm_mod_offset <- fitGLMM(X_pred = pt_test,
                             Y = counts_test[, 3],
                             Y.offset = cell_offset,
                             id.vec = sim_data$subject,
                             adaptive = TRUE,
                             M.glm = 3,
                             return.basis = TRUE,
                             return.GCV = TRUE)
  # generate plots
  plot_glm <- plotModels(test.dyn.res = glm_gene_stats,
                         size.factor.offset = cell_offset,
                         gene = "ABR",
                         pt = pt_test,
                         expr.mat = counts_test)
  plot_gee <- plotModels(test.dyn.res = gee_gene_stats,
                         gene = "ABR",
                         pt = pt_test,
                         expr.mat = counts_test,
                         size.factor.offset = cell_offset,
                         is.gee = TRUE,
                         id.vec = sim_data$subject,
                         cor.structure = "ar1")
  plot_glmm <- plotModels(test.dyn.res = glmm_gene_stats,
                          size.factor.offset = cell_offset,
                          gene = "ABR",
                          pt = pt_test,
                          expr.mat = counts_test,
                          is.glmm = TRUE,
                          id.vec = sim_data$subject)
  # downstream analysis
  set.seed(312)
  gene_clusters <- clusterGenes(glm_gene_stats,
                                pt = pt_test,
                                size.factor.offset = cell_offset,
                                clust.algo = "leiden")
  gene_clust_table <- plotClusteredGenes(glm_gene_stats,
                                         gene.clusters = gene_clusters,
                                         size.factor.offset = cell_offset,
                                         pt = pt_test,
                                         n.cores = 2)
  smoothed_counts <- smoothedCountsMatrix(test.dyn.res = glm_gene_stats,
                                          pt = pt_test,
                                          size.factor.offset = cell_offset,
                                          parallel.exec = TRUE,
                                          n.cores = 2)
  fitted_values_table <- getFittedValues(test.dyn.res = glm_gene_stats,
                                         genes = names(glm_gene_stats),
                                         pt = pt_test,
                                         expr.mat = counts_test,
                                         cell.meta.data = as.data.frame(SummarizedExperiment::colData(sim_data)),
                                         id.vec = sim_data$subject)
  gsea_res <- enrichDynamicGenes(glm_test_results, gene.set.cat = "C2")
})

# run tests
test_that("createCellOffset() output", {
  expect_type(cell_offset, "double")
  expect_length(cell_offset, 300)
  expect_equal(sum(is.na(cell_offset)), 0)
})

test_that("testDynamic() output", {
  expect_s3_class(glm_gene_stats, "scLANE")
  expect_s3_class(gee_gene_stats, "scLANE")
  expect_s3_class(glmm_gene_stats, "scLANE")
  expect_equal(length(glm_gene_stats), 10)
  expect_equal(length(gee_gene_stats), 10)
  expect_equal(length(glmm_gene_stats), 10)
  expect_equal(nrow(glm_gene_stats$ABCF1$Lineage_A$MARGE_Summary), 3)
  expect_equal(ncol(glm_gene_stats$ABCF1$Lineage_A$MARGE_Summary), 5)
  expect_equal(nrow(gee_gene_stats$ABCF1$Lineage_A$MARGE_Summary), 3)
  expect_equal(ncol(gee_gene_stats$ABCF1$Lineage_A$MARGE_Summary), 5)
  expect_equal(nrow(glmm_gene_stats$ABCF1$Lineage_A$MARGE_Summary), 3)
  expect_equal(ncol(glmm_gene_stats$ABCF1$Lineage_A$MARGE_Summary), 5)
  expect_equal(sum(purrr::map_lgl(glm_gene_stats, \(x) x$Lineage_A$Model_Status == "MARGE model OK, null model OK")), 7)
  expect_equal(sum(purrr::map_lgl(gee_gene_stats, \(x) x$Lineage_A$Model_Status == "MARGE model OK, null model OK")), 6)
  expect_equal(sum(purrr::map_lgl(glmm_gene_stats, \(x) x$Lineage_A$Model_Status == "MARGE model OK, null model OK")), 6)
  expect_equal(unique(purrr::map_dbl(glm_gene_stats, \(x) length(x$Lineage_A))), 17)
  expect_equal(unique(purrr::map_dbl(glm_gene_stats, \(x) length(x$Lineage_A))), 17)
  expect_equal(unique(purrr::map_dbl(glm_gene_stats, \(x) length(x$Lineage_A))), 17)
})

test_that("getResultsDE() output", {
  expect_s3_class(glm_test_results, "data.frame")
  expect_s3_class(gee_test_results, "data.frame")
  expect_s3_class(glmm_test_results, "data.frame")
  expect_equal(nrow(glm_test_results), 10)
  expect_equal(nrow(gee_test_results), 10)
  expect_equal(nrow(glmm_test_results), 10)
  expect_equal(ncol(glm_test_results), 15)
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
  expect_equal(nrow(glm_slope_test), 15)
  expect_equal(nrow(gee_slope_test), 15)
  expect_equal(nrow(glmm_slope_test), 16)
  expect_equal(ncol(glm_slope_test), 11)
  expect_equal(ncol(gee_slope_test), 11)
  expect_equal(ncol(glmm_slope_test), 11)
  expect_equal(sum(glm_slope_test$Notes == "MARGE model error, null model OK", na.rm = TRUE), 3)
  expect_equal(sum(gee_slope_test$Notes == "MARGE model error, null model OK", na.rm = TRUE), 4)
  expect_equal(sum(glmm_slope_test$Notes == "MARGE model error, null model OK", na.rm = TRUE), 4)
  expect_equal(sum(glm_slope_test$P_Val_Adj < 0.01, na.rm = TRUE), 10)
  expect_equal(sum(gee_slope_test$P_Val_Adj < 0.01, na.rm = TRUE), 11)
  expect_equal(sum(glmm_slope_test$P_Val_Adj < 0.01, na.rm = TRUE), 7)
})

test_that("nbGAM() output", {
  expect_s3_class(gam_mod_bs, "gamlss")
  expect_s3_class(gam_mod_ps, "gamlss")
  expect_equal(length(coef(gam_mod_bs)), 6)
  expect_equal(length(coef(gam_mod_ps)), 2)
  expect_true(gam_mod_bs$converged)
  expect_true(gam_mod_ps$converged)
})

test_that("marge2() output", {
  expect_s3_class(marge_mod, "marge")
  expect_s3_class(marge_mod_GEE_offset, "marge")
  expect_s3_class(marge_mod$final_mod, "negbin")
  expect_s3_class(marge_mod_GEE_offset$final_mod, "geem")
  expect_length(marge_mod$coef_names, 3)
  expect_equal(ncol(marge_mod$basis_mtx), 3)
  expect_equal(round(marge_mod$final_mod$theta, 1), 6.5)
  expect_equal(marge_mod$model_type, "GLM")
  expect_equal(marge_mod_GEE_offset$model_type, "GEE")
  expect_true(marge_mod$final_mod$converged)
  expect_true(marge_mod_GEE_offset$final_mod$converged)
})
test_that("fitGLMM() output", {
  expect_s3_class(glmm_mod$final_mod, "glmmTMB")
  expect_equal(nrow(coef(glmm_mod$final_mod)$cond$subject), 2)
  expect_equal(ncol(coef(glmm_mod$final_mod)$cond$subject), 3)
  expect_false(glmm_mod$final_mod$modelInfo$REML)
  expect_equal(length(fitted(glmm_mod$final_mod)), 300)
})

test_that("plotModels() output", {
  expect_s3_class(plot_glm, "ggplot")
  expect_s3_class(plot_gee, "ggplot")
  expect_s3_class(plot_glmm, "ggplot")
  expect_equal(ncol(plot_glm$data), 12)
  expect_equal(ncol(plot_gee$data), 12)
  expect_equal(ncol(plot_glmm$data), 12)
  expect_equal(nrow(plot_glm$data), 1200)
  expect_equal(nrow(plot_gee$data), 1200)
  expect_equal(nrow(plot_glmm$data), 1200)
})

test_that("clusterGenes() output", {
  expect_s3_class(gene_clusters, "data.frame")
  expect_equal(ncol(gene_clusters), 3)
  expect_equal(nrow(gene_clusters), 7)
})

test_that("plotClusteredGenes() output", {
  expect_s3_class(gene_clust_table, "data.frame")
  expect_equal(ncol(gene_clust_table), 7)
  expect_equal(nrow(gene_clust_table), 2100)
})

test_that("smoothedCountsMatrix() output", {
  expect_type(smoothed_counts, "list")
  expect_length(smoothed_counts, 1)
  expect_type(smoothed_counts$Lineage_A, "double")
  expect_equal(ncol(smoothed_counts$Lineage_A), 7)
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
