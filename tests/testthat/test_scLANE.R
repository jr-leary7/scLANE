# load data & prepare for testing
load(system.file("testdata/sim_test_data.RData", package = "scLANE"))
genes_to_test <- c(rownames(sim_data)[SummarizedExperiment::rowData(sim_data)$geneStatus_overall == "Dynamic"][1:5],
                   rownames(sim_data)[SummarizedExperiment::rowData(sim_data)$geneStatus_overall == "NotDynamic"][1:5])
counts_test <- t(as.matrix(SingleCellExperiment::counts(sim_data)[genes_to_test, ]))
pt_test <- data.frame(PT = sim_data$cell_time_normed)
# run GLM & GEE tests
sink(tempfile())
glm_gene_stats <- testDynamic(expr.mat = counts_test,
                              pt = pt_test,
                              n.potential.basis.fns = 3,
                              is.gee = FALSE,
                              parallel.exec = FALSE,
                              approx.knot = TRUE,
                              track.time = FALSE)
gee_gene_stats <- testDynamic(expr.mat = counts_test,
                              pt = pt_test,
                              n.potential.basis.fns = 3,
                              is.gee = TRUE,
                              id.vec = sim_data$subject,
                              parallel.exec = FALSE,
                              approx.knot = TRUE,
                              track.time = FALSE)
sink()
# get results tables overall
glm_test_results <- getResultsDE(test.dyn.results = glm_gene_stats)
gee_test_results <- getResultsDE(test.dyn.results = gee_gene_stats)
# get results tables by interval
glm_slope_test <- testSlope(test.dyn.results = glm_gene_stats)
gee_slope_test <- testSlope(test.dyn.results = gee_gene_stats)

test_that("testDynamic() output", {
  expect_equal(class(glm_gene_stats), "scLANE")
  expect_equal(class(gee_gene_stats), "scLANE")
  expect_equal(length(glm_gene_stats), 10)
  expect_equal(length(gee_gene_stats), 10)
  expect_equal(nrow(glm_gene_stats$ABCF1$Lineage_A$MARGE_Summary), 2)
  expect_equal(ncol(glm_gene_stats$ABCF1$Lineage_A$MARGE_Summary), 5)
  expect_equal(nrow(gee_gene_stats$ABCF1$Lineage_A$MARGE_Summary), 2)
  expect_equal(ncol(gee_gene_stats$ABCF1$Lineage_A$MARGE_Summary), 5)
  expect_equal(sum(unname(sapply(glm_gene_stats, \(x) x$Lineage_A$Model_Status)) == "MARGE model OK, null model OK"), 10)
  expect_equal(sum(unname(sapply(gee_gene_stats, \(x) x$Lineage_A$Model_Status)) == "MARGE model OK, null model OK"), 10)
  expect_equal(unique(sapply(glm_gene_stats, \(x) length(x$Lineage_A))), 17)
  expect_equal(unique(sapply(gee_gene_stats, \(x) length(x$Lineage_A))), 17)
})

test_that("getResultsDE() output", {
  expect_equal(class(glm_test_results), "data.frame")
  expect_equal(class(gee_test_results), "data.frame")
  expect_equal(nrow(glm_test_results), 10)
  expect_equal(nrow(gee_test_results), 10)
  expect_equal(ncol(glm_test_results), 14)
  expect_equal(ncol(gee_test_results), 15)
  expect_equal(sum(glm_test_results$Gene_Dynamic_Overall), 5)
  expect_equal(sum(gee_test_results$Gene_Dynamic_Overall), 6)
  expect_equal(sum(is.na(glm_test_results$P_Val)), 0)
  expect_equal(sum(is.na(gee_test_results$P_Val)), 0)
})

test_that("testSlope() output", {
  expect_equal(class(glm_slope_test), "data.frame")
  expect_equal(class(gee_slope_test), "data.frame")
  expect_equal(nrow(glm_slope_test), 10)
  expect_equal(nrow(gee_slope_test), 10)
  expect_equal(ncol(glm_slope_test), 11)
  expect_equal(ncol(gee_slope_test), 11)
  expect_equal(sum(glm_slope_test$Notes == "MARGE model error, null model OK", na.rm = TRUE), 0)
  expect_equal(sum(gee_slope_test$Notes == "MARGE model error, null model OK", na.rm = TRUE), 0)
  expect_equal(sum(glm_slope_test$P_Val_Adj < 0.01, na.rm = TRUE), 5)
  expect_equal(sum(gee_slope_test$P_Val_Adj < 0.01, na.rm = TRUE), 6)
  expect_equal(sum(glm_slope_test$Gene_Dynamic_Overall), 5)
  expect_equal(sum(gee_slope_test$Gene_Dynamic_Overall), 6)
})
