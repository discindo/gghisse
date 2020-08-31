test_that("h_process_recon works", {
  ph <- h_process_recon(hisse_recon=diatoms$cid4_recon)
  
  expect_s3_class(ph$tip_rates, "data.frame")
  expect_s3_class(ph$node_rates, "data.frame")
  expect_s4_class(ph$tree_data, "treedata")
})

test_that("h_scatterplot works", {
  ph <- h_process_recon(hisse_recon=diatoms$cid4_recon)
  
  # parameters
  expect_s3_class(h_scatterplot(ph, parameter="turnover"), "ggplot")
  expect_s3_class(h_scatterplot(ph, parameter="extinct.frac"), "ggplot")
  expect_s3_class(h_scatterplot(ph, parameter="net.div"), "ggplot")
  expect_s3_class(h_scatterplot(ph, parameter="speciation"), "ggplot")
  expect_s3_class(h_scatterplot(ph, parameter="extinction"), "ggplot")
  
  # names
  expect_s3_class(h_scatterplot(ph, parameter="turnover", states_names = c("A", "B")), "ggplot")
  
  # colors
  expect_s3_class(h_scatterplot(ph, colors = c("red", "blue")), "ggplot")
  
  # waiting
  expect_s3_class(h_scatterplot(ph, plot_as_waiting_time = TRUE), "ggplot")
  
  # error
  expect_error(h_scatterplot(ph, parameter="turnoer", states_names = c("A", "B")))
  
})
