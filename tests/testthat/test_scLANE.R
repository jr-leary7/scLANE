# load data & process
load(system.file("testdata/sim_test_data.RData", package = "scLANE"))
sink(tempfile())
gene_stats <- testDynamic(expr.mat = data_list$simdata,
                          pt = data_list$simpt,
                          parallel.exec = FALSE,
                          n.potential.basis.fns = 5,
                          track.time = FALSE)
sink()
test_results <- getResultsDE(gene_stats,
                             p.adj.method = "bonferroni",
                             fdr.cutoff = 0.01)
slope_test <- testSlope(dyn.results = gene_stats,
                        p.adj.method = "bonferroni",
                        fdr.cutoff = 0.01)

# testDynamic()
test_that("testDynamic() output", {
  expect_equal(class(gene_stats), "list")
  expect_equal(length(gene_stats), 10)
  expect_equal(nrow(gene_stats$PGAM1$MARGE_Summary), 3)
  expect_equal(ncol(gene_stats$PGAM1$MARGE_Summary), 5)
  expect_equal(sum(unname(sapply(gene_stats, function(x) x$Model_Status)) == "All clear"), 7)
  expect_equal(unique(sapply(gene_stats, length)), 13)
})

# getResultsDE()
test_that("getResultsDE() output", {
  expect_equal(class(test_results), "data.frame")
  expect_equal(nrow(test_results), 10)
  expect_equal(ncol(test_results), 10)
  expect_equal(sum(test_results$Gene_Dynamic), 2)
  expect_equal(sum(is.na(test_results$LogLik_Null)), 0)
})

# testSlope()
test_that("testSlope() output", {
  expect_equal(class(slope_test), "data.frame")
  expect_equal(nrow(slope_test), 15)
  expect_equal(ncol(slope_test), 9)
  expect_equal(sum(slope_test$Notes == "MARGE model error", na.rm = TRUE), 3)
  expect_equal(sum(slope_test$P_Val_Adj_Signif), 3)
  expect_equal(sum(slope_test$Gene_Dynamic), 4)
})
