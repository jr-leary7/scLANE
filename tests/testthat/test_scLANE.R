# load data & prepare for testing
load(system.file("testdata/sim_test_data.rda", package = "scLANE"))
sim_data_seu <- Seurat::as.Seurat(sim_data)
cell_offset <- createCellOffset(sim_data)
cell_offset_seu <- createCellOffset(sim_data_seu)
genes_to_test <- c(rownames(sim_data)[SummarizedExperiment::rowData(sim_data)$geneStatus_overall == "Dynamic"][1:10],
                   rownames(sim_data)[SummarizedExperiment::rowData(sim_data)$geneStatus_overall == "NotDynamic"][1:10])
counts_test <- t(as.matrix(SingleCellExperiment::counts(sim_data)[genes_to_test, ]))
pt_test <- data.frame(PT = sim_data$cell_time_normed)

# test internal functions used in model fitting
max_span_res <- max_span(X_red = pt_test[, 1], q = 1)
min_span_res <- min_span(X_red = pt_test[, 1], q = 1)
Y_exp <- sim_data@assays@data$counts[genes_to_test[1], ]
stat_out_res <- stat_out(Y = Y_exp,
                         B1 = matrix(rnorm(ncol(sim_data)), nrow = ncol(sim_data), ncol = 1),
                         TSS = sum((Y_exp - mean(Y_exp))^2),
                         GCV.null = sum((Y_exp - mean(Y_exp))^2) / (ncol(sim_data) * (1 - (1 / ncol(sim_data)))^2))
null_stat_glm <- stat_out_score_glm_null(Y = Y_exp,
                                         B_null = rep(1, ncol(sim_data)))
null_stat_gee <- stat_out_score_gee_null(Y = Y_exp,
                                         B_null = matrix(1, ncol = 1, nrow = ncol(sim_data)),
                                         id.vec = sim_data$subject,
                                         cor.structure = "ar1",
                                         theta.hat = 1)
tp1_res <- tp1(x = rnorm(30), t = 0)
tp2_res <- tp2(x = rnorm(30), t = 0)

# generate scLANE results w/ all three modes
withr::with_output_sink(tempfile(), {
  # choose candidate genes
  candidate_genes <- chooseCandidateGenes(sim_data_seu,
                                          group.by.subject = TRUE,
                                          id.vec = sim_data_seu$subject,
                                          n.desired.genes = 50L)
  # run GLM, GEE, & GLMM tests
  glm_gene_stats <- testDynamic(sim_data,
                                pt = pt_test,
                                genes = genes_to_test,
                                n.potential.basis.fns = 5,
                                size.factor.offset = cell_offset,
                                n.cores = 1L)
  gee_gene_stats <- testDynamic(sim_data,
                                pt = pt_test,
                                genes = genes_to_test,
                                n.potential.basis.fns = 5,
                                size.factor.offset = cell_offset,
                                is.gee = TRUE,
                                cor.structure = "ar1",
                                id.vec = sim_data$subject,
                                n.cores = 2L)
  glmm_gene_stats <- testDynamic(sim_data,
                                 pt = pt_test,
                                 genes = genes_to_test,
                                 size.factor.offset = cell_offset,
                                 n.potential.basis.fns = 3,
                                 is.glmm = TRUE,
                                 glmm.adaptive = TRUE,
                                 id.vec = sim_data$subject,
                                 n.cores = 2L,
                                 verbose = FALSE)
  # S3 summary method
  scLANE_summary <- summary(glm_gene_stats)
  # get results tables overall
  glm_test_results <- getResultsDE(glm_gene_stats, n.cores = 1L)
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
  marge_mod_stripped <- stripGLM(marge_mod$final_mod)
  # run GLM model -- with offset
  marge_mod_offset <- marge2(X_pred = pt_test,
                             Y = counts_test[, 3],
                             Y.offset = cell_offset,
                             M = 5,
                             return.basis = TRUE,
                             return.GCV = TRUE,
                             return.WIC = TRUE)
  null_mod_offset <- MASS::glm.nb(counts_test[, 3] ~ pt_test$PT + offset(log(1 / cell_offset)),
                                  method = "glm.fit2",
                                  y = FALSE,
                                  model = FALSE,
                                  init.theta = 1,
                                  link = log)
  glm_lrt <- modelLRT(mod.1 = marge_mod_offset, mod.0 = null_mod_offset)
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
  # fit null GEE model
  null_mod_GEE <- geeM::geem(Y ~ 1 + offset(log(1 / cell_offset)),
                             family = MASS::negative.binomial(50, link = log),
                             data = data.frame(Y = counts_test[, 3]),
                             corstr = "ar1",
                             scale.fix = FALSE,
                             sandwich = FALSE)
  # run GEE Wald test
  wald_test <- waldTestGEE(marge_mod_GEE_offset, mod.0 = null_mod_GEE)
  # run GEE Score test
  score_test <- scoreTestGEE(marge_mod_GEE_offset, 
                             mod.0 = null_mod_GEE, 
                             alt.df = as.data.frame(marge_mod_GEE_offset$basis_mtx), 
                             null.df = data.frame(Y = counts_test[, 3]))
  # run GLMM model -- no offset
  glmm_mod <- fitGLMM(pt_test,
                      Y = counts_test[, 4],
                      id.vec = sim_data$subject,
                      adaptive = TRUE,
                      M.glm = 3,
                      return.basis = TRUE,
                      return.GCV = TRUE)
  # run GLMM model -- with offset
  glmm_mod_offset <- fitGLMM(pt_test,
                             Y = counts_test[, 4],
                             Y.offset = cell_offset,
                             id.vec = sim_data$subject,
                             adaptive = TRUE,
                             M.glm = 3,
                             return.basis = TRUE,
                             return.GCV = FALSE)
  # bootstrap GLMM random effects
  re_sumy <- bootstrapRandomEffects(glmm_mod_offset,
                                    id.vec = sim_data$subject,
                                    Y.offset = cell_offset,
                                    n.boot = 100L,
                                    n.cores = 2L)
  # generate plots
  plot_glm <- plotModels(test.dyn.res = glm_gene_stats,
                         size.factor.offset = cell_offset,
                         gene = "ABR",
                         pt = pt_test,
                         expr.mat = sim_data,
                         plot.null = TRUE,
                         plot.glm = TRUE,
                         plot.gam = TRUE) +
              theme_scLANE()
  plot_gee <- plotModels(test.dyn.res = gee_gene_stats,
                         gene = "ABR",
                         pt = pt_test,
                         expr.mat = sim_data,
                         size.factor.offset = cell_offset,
                         is.gee = TRUE,
                         id.vec = sim_data$subject,
                         cor.structure = "ar1",
                         plot.null = TRUE,
                         plot.glm = TRUE,
                         plot.gam = TRUE)
  plot_glmm <- plotModels(test.dyn.res = glmm_gene_stats,
                          size.factor.offset = cell_offset,
                          gene = "ABR",
                          pt = pt_test,
                          expr.mat = sim_data,
                          is.glmm = TRUE,
                          id.vec = sim_data$subject,
                          plot.null = TRUE,
                          plot.glm = TRUE,
                          plot.gam = TRUE)
  coef_plot_glm <- plotModelCoefs(test.dyn.res = glm_gene_stats,
                                  gene = glm_test_results$Gene[10],
                                  pt = pt_test,
                                  expr.mat = sim_data,
                                  size.factor.offset = cell_offset,
                                  lineage = "A",
                                  log1p.norm = TRUE)
  # gene clustering
  set.seed(312)
  gene_clusters_leiden <- clusterGenes(test.dyn.res = glm_gene_stats,
                                       pt = pt_test,
                                       size.factor.offset = cell_offset,
                                       clust.algo = "leiden")
  gene_clusters_kmeans <- clusterGenes(test.dyn.res = glm_gene_stats,
                                       pt = pt_test,
                                       size.factor.offset = cell_offset,
                                       clust.algo = "kmeans")
  gene_clusters_hclust <- clusterGenes(test.dyn.res = glm_gene_stats,
                                       pt = pt_test,
                                       size.factor.offset = cell_offset,
                                       clust.algo = "hclust")
  gene_clust_table <- plotClusteredGenes(test.dyn.res = glm_gene_stats,
                                         gene.clusters = gene_clusters_leiden,
                                         size.factor.offset = cell_offset,
                                         pt = pt_test,
                                         n.cores = 2L)
  # smoothed dynamics
  smoothed_counts <- smoothedCountsMatrix(test.dyn.res = glm_gene_stats,
                                          pt = pt_test,
                                          size.factor.offset = cell_offset,
                                          n.cores = 2L)
  sorted_genes <- sortGenesHeatmap(heatmap.mat = smoothed_counts$Lineage_A,
                                   pt.vec = pt_test$PT)
  fitted_values_table <- getFittedValues(test.dyn.res = glm_gene_stats,
                                         genes = names(glm_gene_stats),
                                         pt = pt_test,
                                         size.factor.offset = cell_offset,
                                         expr.mat = sim_data,
                                         cell.meta.data = as.data.frame(SummarizedExperiment::colData(sim_data)),
                                         id.vec = sim_data$subject)
  # gene embeddings
  gene_embedding_pca <- embedGenes(smoothed.counts = smoothed_counts$Lineage_A,
                                   pca.init = TRUE,
                                   pc.embed = 5,
                                   pc.return = 2,
                                   k.param = 5,
                                   random.seed = 312)
  gene_embedding <- embedGenes(smoothed.counts = smoothed_counts$Lineage_A,
                               pc.embed = 5,
                               pc.return = 2,
                               k.param = 5,
                               random.seed = 312)
  # gene program scoring
  sim_data <- geneProgramScoring(sim_data,
                                 genes = gene_embedding$gene,
                                 gene.clusters = gene_embedding$leiden)
  sim_data_seu <- geneProgramScoring(sim_data_seu,
                                     genes = gene_embedding$gene,
                                     gene.clusters = gene_embedding$leiden,
                                     n.cores = 1L)
  # gene program significance
  program_significance <- geneProgramSignificance(list(sim_data$cluster_0),
                                                  pt = pt_test$PT,
                                                  program.labels = c("Cluster0"))
  # gene program drivers
  program_drivers <- geneProgramDrivers(sim_data,
                                        genes = gene_embedding$gene,
                                        gene.program = sim_data$cluster_0)
  program_drivers_seu <- geneProgramDrivers(sim_data_seu,
                                            genes = gene_embedding$gene,
                                            gene.program = sim_data_seu$cluster_0)
  # enrichment analysis
  gsea_res <- enrichDynamicGenes(glm_test_results, species = "hsapiens")
  # coefficients
  coef_summary_glm <- summarizeModel(marge_mod_offset, pt = pt_test)
  coef_summary_gee <- summarizeModel(marge_mod_GEE_offset, pt = pt_test)
  coef_summary_glmm <- summarizeModel(glmm_mod_offset,
                                      pt = pt_test,
                                      is.glmm = TRUE)
  # cutpoints
  knot_df <- getKnotDist(glm_gene_stats)
  # convolution
  dyn_convolve <- npConvolve(counts_test[, 1], conv.kernel = rep(1/10, 10))
})

# run tests
test_that("internal marge functions", {
  expect_type(min_span_res, "double")
  expect_type(max_span_res, "double")
  expect_type(stat_out_res, "list")
  expect_type(null_stat_glm, "list")
  expect_type(null_stat_gee, "list")
  expect_length(stat_out_res, 4)
  expect_length(null_stat_glm, 5)
  expect_length(null_stat_gee, 8)
  expect_type(tp1_res, "double")
  expect_type(tp2_res, "double")
})

test_that("chooseCandidateGenes() output", {
  expect_type(candidate_genes, "character")
  expect_length(candidate_genes, 50)
})

test_that("createCellOffset() output", {
  expect_type(cell_offset, "double")
  expect_type(cell_offset_seu, "double")
  expect_length(cell_offset, 300)
  expect_length(cell_offset_seu, 300)
  expect_false(any(is.na(cell_offset)))
  expect_false(any(is.na(cell_offset_seu)))
})

test_that("testDynamic() output", {
  expect_s3_class(glm_gene_stats, "scLANE")
  expect_s3_class(gee_gene_stats, "scLANE")
  expect_s3_class(glmm_gene_stats, "scLANE")
  expect_length(glm_gene_stats, 20)
  expect_length(gee_gene_stats, 20)
  expect_length(glmm_gene_stats, 20)
  expect_s3_class(glm_gene_stats[[glm_test_results$Gene[1]]]$Lineage_A$MARGE_Summary, "data.frame")
  expect_s3_class(gee_gene_stats[[gee_test_results$Gene[1]]]$Lineage_A$MARGE_Summary, "data.frame")
  expect_s3_class(glmm_gene_stats[[glmm_test_results$Gene[1]]]$Lineage_A$MARGE_Summary, "data.frame")
  expect_gt(nrow(glm_gene_stats[[glm_test_results$Gene[1]]]$Lineage_A$MARGE_Summary), 0)
  expect_gt(nrow(gee_gene_stats[[gee_test_results$Gene[1]]]$Lineage_A$MARGE_Summary), 0)
  expect_gt(nrow(glmm_gene_stats[[glmm_test_results$Gene[1]]]$Lineage_A$MARGE_Summary), 0)
  expect_gt(sum(purrr::map_lgl(glm_gene_stats, \(x) x$Lineage_A$Model_Status == "MARGE model OK, null model OK")), 0)
  expect_gt(sum(purrr::map_lgl(gee_gene_stats, \(x) x$Lineage_A$Model_Status == "MARGE model OK, null model OK")), 0)
  expect_gt(sum(purrr::map_lgl(glmm_gene_stats, \(x) x$Lineage_A$Model_Status == "MARGE model OK, null model OK")), 0)
})

test_that("summary() output", {
  expect_s3_class(scLANE_summary, "summary.scLANE")
  expect_length(scLANE_summary, 5)
})

test_that("getResultsDE() output", {
  expect_s3_class(glm_test_results, "data.frame")
  expect_s3_class(gee_test_results, "data.frame")
  expect_s3_class(glmm_test_results, "data.frame")
  expect_gt(nrow(glm_test_results), 0)
  expect_gt(nrow(gee_test_results), 0)
  expect_gt(nrow(glmm_test_results), 0)
  expect_gt(sum(glm_test_results$Gene_Dynamic_Overall), 0)
  expect_gt(sum(gee_test_results$Gene_Dynamic_Overall), 0)
  expect_gt(sum(glmm_test_results$Gene_Dynamic_Overall), 0)
})

test_that("testSlope() output", {
  expect_s3_class(glm_slope_test, "data.frame")
  expect_s3_class(gee_slope_test, "data.frame")
  expect_s3_class(glmm_slope_test, "data.frame")
  expect_gt(nrow(glm_slope_test), 0)
  expect_gt(nrow(gee_slope_test), 0)
  expect_gt(nrow(glmm_slope_test), 0)
  expect_gt(sum(glm_slope_test$P_Val_Adj < 0.01, na.rm = TRUE), 0)
  expect_gt(sum(gee_slope_test$P_Val_Adj < 0.01, na.rm = TRUE), 0)
  expect_gt(sum(glmm_slope_test$P_Val_Adj < 0.01, na.rm = TRUE), 0)
})

test_that("nbGAM() output", {
  expect_s3_class(gam_mod_bs, "gamlss")
  expect_s3_class(gam_mod_ps, "gamlss")
  expect_s3_class(gam_mod_ps_mix, "gamlss")
  expect_true(gam_mod_bs$converged)
  expect_true(gam_mod_ps$converged)
  expect_true(gam_mod_ps_mix$converged)
})

test_that("marge2() output -- GLM backend", {
  expect_s3_class(marge_mod, "marge")
  expect_s3_class(marge_mod_offset, "marge")
  expect_s3_class(marge_mod$final_mod, "negbin")
  expect_s3_class(marge_mod_offset$final_mod, "negbin")
  expect_s3_class(marge_mod_stripped, "negbin")
  expect_equal(marge_mod$model_type, "GLM")
  expect_equal(marge_mod_offset$model_type, "GLM")
  expect_true(marge_mod$final_mod$converged)
  expect_true(marge_mod_offset$final_mod$converged)
})

test_that("marge2() output -- GEE backend", {
  expect_s3_class(marge_mod_GEE, "marge")
  expect_s3_class(marge_mod_GEE_offset, "marge")
  expect_s3_class(marge_mod_GEE$final_mod, "geem")
  expect_s3_class(marge_mod_GEE_offset$final_mod, "geem")
  expect_equal(marge_mod_GEE$model_type, "GEE")
  expect_equal(marge_mod_GEE_offset$model_type, "GEE")
  expect_true(marge_mod_GEE$final_mod$converged)
  expect_true(marge_mod_GEE_offset$final_mod$converged)
})

test_that("Statistical testing output", {
  expect_type(wald_test, "list")
  expect_length(wald_test, 4)
  expect_type(wald_test$P_Val, "double")
  expect_true(is.na(wald_test$Notes))
  expect_type(score_test, "list")
  expect_length(score_test, 4)
  expect_type(score_test$P_Val, "double")
  expect_true(is.na(score_test$Notes))
  expect_type(glm_lrt, "list")
  expect_length(glm_lrt, 7)
  expect_type(glm_lrt$P_Val, "double")
  expect_true(is.na(glm_lrt$Notes))
})

test_that("fitGLMM() output", {
  expect_s3_class(glmm_mod$final_mod, "glmmTMB")
  expect_s3_class(glmm_mod_offset$final_mod, "glmmTMB")
  expect_equal(nrow(coef(glmm_mod$final_mod)$cond$subject), 2)
  expect_equal(nrow(coef(glmm_mod_offset$final_mod)$cond$subject), 2)
  expect_false(glmm_mod$final_mod$modelInfo$REML)
  expect_false(glmm_mod_offset$final_mod$modelInfo$REML)
  expect_length(fitted(glmm_mod$final_mod), 300)
  expect_length(fitted(glmm_mod_offset$final_mod), 300)
  expect_equal(glmm_mod$model_type, "GLMM")
  expect_equal(glmm_mod_offset$model_type, "GLMM")
  expect_s3_class(re_sumy, "data.frame")
  expect_equal(nrow(re_sumy), 4)
  expect_equal(ncol(re_sumy), 4)
})

test_that("plotModels() output", {
  expect_s3_class(plot_glm, "ggplot")
  expect_s3_class(plot_gee, "ggplot")
  expect_s3_class(plot_glmm, "ggplot")
  expect_equal(ncol(plot_glm$data), 12)
  expect_equal(ncol(plot_gee$data), 12)
  expect_equal(ncol(plot_glmm$data), 12)
})

test_that("plotModelCoefs() output", {
  expect_s3_class(coef_plot_glm, "ggplot")
  expect_s3_class(coef_plot_glm, "ggarrange")
})

test_that("clusterGenes() output", {
  expect_s3_class(gene_clusters_leiden, "data.frame")
  expect_s3_class(gene_clusters_kmeans, "data.frame")
  expect_s3_class(gene_clusters_hclust, "data.frame")
  expect_equal(ncol(gene_clusters_leiden), 3)
  expect_equal(ncol(gene_clusters_kmeans), 3)
  expect_equal(ncol(gene_clusters_hclust), 3)
})

test_that("plotClusteredGenes() output", {
  expect_s3_class(gene_clust_table, "data.frame")
  expect_equal(ncol(gene_clust_table), 7)
})

test_that("smoothedCountsMatrix() output", {
  expect_type(smoothed_counts, "list")
  expect_length(smoothed_counts, 1)
  expect_type(smoothed_counts$Lineage_A, "double")
})

test_that("embedGenes() output", {
  expect_s3_class(gene_embedding, "data.frame")
  expect_s3_class(gene_embedding_pca, "data.frame")
  expect_equal(ncol(gene_embedding), 6)
  expect_equal(ncol(gene_embedding_pca), 6)
})

test_that("geneProgramScoring() output", {
  expect_equal(colnames(SummarizedExperiment::colData(sim_data))[7], "cluster_0")
  expect_equal(colnames(SummarizedExperiment::colData(sim_data))[8], "cluster_1")
  expect_equal(colnames(sim_data_seu@meta.data)[10], "cluster_0")
  expect_equal(colnames(sim_data_seu@meta.data)[11], "cluster_1")
})

test_that("geneProgramSignificance() output", {
  expect_s3_class(program_significance, "data.frame")
  expect_equal(ncol(program_significance), 5)
})

test_that("geneProgramDrivers() output", {
  expect_s3_class(program_drivers, "data.frame")
  expect_s3_class(program_drivers_seu, "data.frame")
  expect_equal(ncol(program_drivers), 4)
  expect_equal(ncol(program_drivers_seu), 4)
})

test_that("sortGenesHeatmap() output", {
  expect_type(sorted_genes, "character")
  expect_length(sorted_genes, ncol(smoothed_counts$Lineage_A))
})

test_that("getFittedValues() output", {
  expect_s3_class(fitted_values_table, "data.frame")
  expect_equal(ncol(fitted_values_table), 25)
})

test_that("enrichDynamicGenes() output", {
  expect_type(gsea_res, "list")
  expect_length(gsea_res, 2)
  expect_s3_class(gsea_res$result, "data.frame")
})

test_that("summarizeModels() output", {
  expect_type(coef_summary_glm, "list")
  expect_type(coef_summary_gee, "list")
  expect_length(coef_summary_glm, 3)
  expect_length(coef_summary_gee, 3)
  expect_type(coef_summary_glm$Slope.Segment, "double")
  expect_type(coef_summary_gee$Slope.Segment, "double")
})

test_that("getKnotDist() output", {
  expect_s3_class(knot_df, "data.frame")
  expect_equal(ncol(knot_df), 3)
})

test_that("npConvolve() output", {
  expect_type(dyn_convolve, "double")
  expect_length(dyn_convolve, nrow(counts_test))
})
