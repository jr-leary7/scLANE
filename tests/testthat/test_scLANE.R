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
glm_test_results <- getResultsDE(test.dyn.results = glm_gene_stats)
gee_test_results <- getResultsDE(test.dyn.results = gee_gene_stats)
glmm_test_results <- getResultsDE(test.dyn.results = glmm_gene_stats)
# get results tables by interval
glm_slope_test <- testSlope(test.dyn.results = glm_gene_stats)
gee_slope_test <- testSlope(test.dyn.results = gee_gene_stats)
glmm_slope_test <- testSlope(test.dyn.results = glmm_gene_stats)
# run NB GAM
gam_mod <- nbGAM(expr = counts_test[, 1], pt = pt_test)
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
  expect_equal(sum(glm_test_results$Gene_Dynamic_Overall), 5)
  expect_equal(sum(gee_test_results$Gene_Dynamic_Overall), 4)
  expect_equal(sum(glmm_test_results$Gene_Dynamic_Overall), 2)
  expect_equal(sum(is.na(glm_test_results$P_Val)), 0)
  expect_equal(sum(is.na(gee_test_results$P_Val)), 0)
  expect_equal(sum(is.na(glmm_test_results$P_Val)), 2)
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
  expect_equal(sum(glm_slope_test$Gene_Dynamic_Overall), 5)
  expect_equal(sum(gee_slope_test$Gene_Dynamic_Overall), 4)
  expect_equal(sum(glmm_slope_test$Gene_Dynamic_Overall), 0)
})

test_that("nbGAM() output", {
  expect_s3_class(gam_mod, "gamlss")
  expect_equal(length(coef(gam_mod)), 4)
  expect_equal(round(gam_mod$G.deviance), 585)
  expect_true(gam_mod$converged)
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
  expect_false(glmm_mod$final_mod$sdr$pdHess)
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
