context("iNEXTbeta3D")
test_that("iNEXTbeta3D", {
  
  ## TD for abundance data
  # Coverage-based 
  data(Brazil_rainforests)
  output1c = iNEXTbeta3D(data = Brazil_rainforests, diversity = 'TD', q = 0,
                         datatype = 'abundance', base = "coverage", nboot = 10, level = 0.9)
  expect_is(output1c, "iNEXTbeta3D")
  expect_output(str(output1c), "List of 2")
  expect_equal(length(output1c), length(Brazil_rainforests))
  expect_equal(length(output1c[[1]]), 7)
  expect_equal(length(output1c[[1]]$gamma$SC), 1)
  
  
  # Size-based 
  data(Brazil_rainforests)
  output1s = iNEXTbeta3D(data = Brazil_rainforests, diversity = 'TD',  q = 1,
                         datatype = 'abundance', base = "size", nboot = 10, level = 200)
  expect_is(output1s, "iNEXTbeta3D")
  expect_output(str(output1s), "List of 2")
  expect_equal(length(output1s), length(Brazil_rainforests))
  expect_equal(length(output1s[[1]]), 2)
  expect_equal(length(output1s[[1]]$gamma$SC), 1)
  
  
  ## PD for abundance data
  # Coverage-based 
  data(Brazil_rainforests)
  data(Brazil_tree)
  output2c = iNEXTbeta3D(data = Brazil_rainforests[1:2], diversity = 'PD',  q = 2,
                         datatype = 'abundance', base = "coverage", nboot = 10, 
                         PDtree = Brazil_tree, PDreftime = NULL, PDtype = 'meanPD', level = 0.95)
  expect_is(output2c, "iNEXTbeta3D")
  expect_output(str(output2c), "List of 2")
  expect_equal(length(output2c), length(Brazil_rainforests[1:2]))
  expect_equal(length(output2c[[1]]), 7)
  expect_equal(length(output2c[[1]]$gamma$SC), 1)
  
  
  # Size-based 
  data(Brazil_rainforests)
  data(Brazil_tree)
  output2s = iNEXTbeta3D(data = Brazil_rainforests[1:2], diversity = 'PD',  q = 1,
                         datatype = 'abundance', base = "size", nboot = 10, 
                         PDtree = Brazil_tree, PDreftime = NULL, PDtype = 'meanPD', level = 200)
  expect_is(output2s, "iNEXTbeta3D")
  expect_output(str(output2s), "List of 2")
  expect_equal(length(output2s), length(Brazil_rainforests[1:2]))
  expect_equal(length(output2s[[1]]), 2)
  expect_equal(length(output2s[[1]]$gamma$SC), 1)
  
  
  ## FD for abundance data
  # Coverage-based 
  data(Brazil_rainforests)
  data(Brazil_distM)
  output3c = iNEXTbeta3D(data = Brazil_rainforests, diversity = 'FD',  q = 0,
                         datatype = 'abundance', base = "coverage", nboot = 10, 
                         FDdistM = Brazil_distM, FDtype = 'AUC', FDcut_number = 30, level = 0.99)
  expect_is(output3c, "iNEXTbeta3D")
  expect_output(str(output3c), "List of 2")
  expect_equal(length(output3c), length(Brazil_rainforests))
  expect_equal(length(output3c[[1]]), 7)
  expect_equal(length(output3c[[1]]$gamma$SC), 1)
  
  
  # Size-based 
  data(Brazil_rainforests)
  data(Brazil_distM)
  output3s = iNEXTbeta3D(data = Brazil_rainforests, diversity = 'FD',  q = 1,
                         datatype = 'abundance', base = "size", nboot = 10, 
                         FDdistM = Brazil_distM, FDtype = 'AUC', FDcut_number = 30, level = 200)
  expect_is(output3s, "iNEXTbeta3D")
  expect_output(str(output3s), "List of 2")
  expect_equal(length(output3s), length(Brazil_rainforests))
  expect_equal(length(output3s[[1]]), 2)
  expect_equal(length(output3s[[1]]$gamma$SC), 1)
  
  
  ## TD for incidence data
  # Coverage-based 
  data(Second_growth_forests)
  output4c = iNEXTbeta3D(data = Second_growth_forests, diversity = 'TD',  q = 1,
                         datatype = 'incidence_raw', base = "coverage", nboot = 10, level = 0.99)
  expect_is(output4c, "iNEXTbeta3D")
  expect_output(str(output4c), "List of 2")
  expect_equal(length(output4c), length(Second_growth_forests))
  expect_equal(length(output4c[[1]]), 7)
  expect_equal(length(output4c[[1]]$gamma$SC), 1)
  
  
  # Coverage-based (by pair)
  data(Second_growth_forests)
  output4c_pair = iNEXTbeta3D(data = Second_growth_forests, diversity = 'TD',  q = 2,
                              datatype = 'incidence_raw', base = "coverage", nboot = 10, level = 0.99, by_pair = TRUE)
  expect_is(output4c_pair, "iNEXTbeta3D")
  expect_output(str(output4c_pair), "List of 2")
  expect_equal(length(output4c_pair), length(Second_growth_forests))
  expect_equal(length(output4c_pair[[1]]), 7)
  expect_equal(length(unique(output4c_pair[[1]]$gamma$SC)), 1)
  expect_equal(length(unique(output4c_pair$CR_2005_2011_2017$gamma$Pair)), 3)
  
  
  # Size-based 
  data(Second_growth_forests)
  output4s = iNEXTbeta3D(data = Second_growth_forests, diversity = 'TD', q = 2,
                         datatype = 'incidence_raw', base = "size", nboot = 10, level = 120)
  expect_is(output4s, "iNEXTbeta3D")
  expect_output(str(output4s), "List of 2")
  expect_equal(length(output4s), length(Second_growth_forests))
  expect_equal(length(output4s[[1]]), 2)
  expect_equal(length(output4s[[1]]$gamma$SC), 1)
  
  
  # Size-based (by pair)
  data(Second_growth_forests)
  output4s_pair = iNEXTbeta3D(data = Second_growth_forests, diversity = 'TD', q = 1,
                              datatype = 'incidence_raw', base = "size", nboot = 10, level = 120, by_pair = TRUE)
  expect_is(output4s_pair, "iNEXTbeta3D")
  expect_output(str(output4s_pair), "List of 2")
  expect_equal(length(output4s_pair), length(Second_growth_forests))
  expect_equal(length(output4s_pair[[1]]), 2)
  expect_equal(length(unique(output4s_pair[[1]]$gamma$mT)), 1)
  expect_equal(length(unique(output4s_pair[[1]]$gamma$Pair)), 3)
  
  
  ## FD for abundance data under single threshold
  # Coverage-based 
  data(Brazil_rainforests)
  data(Brazil_distM)
  outputtauc = iNEXTbeta3D(data = Brazil_rainforests, diversity = 'FD',  q = 0,
                           datatype = 'abundance', base = "coverage", nboot = 10, 
                           FDdistM = Brazil_distM, FDtype = 'tau_value', FDtau = 0.2, level = 0.99)
  expect_is(outputtauc, "iNEXTbeta3D")
  expect_output(str(outputtauc), "List of 2")
  expect_equal(length(outputtauc), length(Brazil_rainforests))
  expect_equal(length(outputtauc[[1]]), 7)
  expect_equal(length(outputtauc[[1]]$gamma$SC), 1)
  
  
  # Size-based 
  data(Brazil_rainforests)
  data(Brazil_distM)
  outputtaus = iNEXTbeta3D(data = Brazil_rainforests, diversity = 'FD',  q = 1,
                           datatype = 'abundance', base = "size", nboot = 10, 
                           FDdistM = Brazil_distM, FDtype = 'tau_value', FDtau = NULL, level = 200)
  expect_is(outputtaus, "iNEXTbeta3D")
  expect_output(str(outputtaus), "List of 2")
  expect_equal(length(outputtaus), length(Brazil_rainforests))
  expect_equal(length(outputtaus[[1]]), 2)
  expect_equal(length(outputtaus[[1]]$gamma$SC), 1)
  
  
})

