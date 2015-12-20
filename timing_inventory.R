################################################################################
# Copyright (c) 2015 Farzad Noorian. Distributed under Apache License v2.0
################################################################################
# simulating and measuring the Times

# use this to accurately measure short running tests
REPEATS_FOR_ACCURACY_DJ_FULL = 10
REPEATS_FOR_ACCURACY_DJ_APPROX = 200 
REPEATS_FOR_ACCURACY_MSE = 200 

#######
# only do measurement for the first training range
train_range = train_ranges[[1]]$train
test_range = train_ranges[[1]]$test

#######
# Get the errors list
t1 = system.time({
  e_list = lapply(all_model_orders, function(model) {
    print(model)
    e = lapply(v_collection[1:num_total_tests], function(v) { 
      predict_inventory_error(model, v[train_range], v[test_range])
    })
  })
})

#######
## use Full matrix for ranking errors

t2 = system.time({
  for (count in 1:REPEATS_FOR_ACCURACY_DJ_FULL) {
    e_rhc_full_list = NULL
    model_counter = 0
    for (i in seq_along(all_model_orders)) {
      order_list = all_model_orders[[i]]
      
      e_rhc = 0
      for (test_number in 1:num_total_tests) {
        V = v_collection[[test_number]][test_range]
        e = e_list[[i]][[test_number]]

        e_rhc = e_rhc + (t(e) %*% theta %*% e + t(e) %*% omega %*% V) 
      }
      
      e_rhc_full_list = c(e_rhc_full_list, e_rhc / num_total_tests)
    }
  }
})

t2 = t2 / REPEATS_FOR_ACCURACY_DJ_FULL

#######
## use simulator to rank errors
N = length(horizons)

cost_U <- function(U, V) t(U) %*% JA(N) %*% U + 2 * t(U) %*% JB(N) %*% V + t(V) %*% JC(N) %*% V

# first collect the actual minimum costs
min_cost = NULL
for (test_number in 1:num_total_tests) {
  V = v_collection[[test_number]][test_range]

  U = RHC_simulator(V, e * 0, horizons)$U

  min_cost = c(min_cost, cost_U(U, V))
}

# now compute delta J using "cost - min_cost"
e_simulator_list = NULL

t4 = system.time({
  for (i in seq_along(all_model_orders)) {

    e_sim = 0
    for (test_number in 1:num_total_tests) {
      e = e_list[[i]][[test_number]]

      V = v_collection[[test_number]][test_range]
      U = RHC_simulator(V, e, horizons)$U 

      
      e_sim = e_sim + (cost_U(U, V) - min_cost[test_number])
    }
    
    e_simulator_list = c(e_simulator_list, e_sim / num_total_tests)
  }
})

#######
## use approximated matrix for ranking Errors
t3 = system.time({
  for (count in 1:REPEATS_FOR_ACCURACY_DJ_APPROX) {

    e_rhc_apprx_list = sapply(e_list, function(errs_for_model_order) {
      mean(sapply(errs_for_model_order, DeltaJ)) # get mean of DeltaJ for the certain model order
    })
  }
})
t3 = t3 / REPEATS_FOR_ACCURACY_DJ_APPROX

#######
## MSE error timing
t5 = system.time({
  for (count in 1: REPEATS_FOR_ACCURACY_MSE) {
    e_mse_full_list = sapply(e_list, function(errs_for_model_order) {
      mean(sapply(errs_for_model_order, MSE)) # get mean of MSEs for the certain model order
    })
  }
})

t5 = t5 / REPEATS_FOR_ACCURACY_MSE

#######
# create the table
table_timing = rbind(
  data.frame(names = "Step-by-step simulation", time = round(t4[1], 3), speedup = round(t4[1] / t4[1], 0) , computed_cost = mean(e_simulator_list)),
  data.frame(names = "Closed form $\\Delta J$", time = round(t2[1], 3), speedup = round(t4[1] / t2[1], 1) , computed_cost = mean(e_rhc_full_list)),
  data.frame(names = "Diagonal $\\Delta J$", time = round(t3[1], 3), speedup = round(t4[1] / t3[1], 1) , computed_cost = mean(e_rhc_full_list)),
  data.frame(names = "MSE", time = round(t5[1], 3), speedup = round(t4[1] / t5[1], 0) , computed_cost = mean(e_mse_full_list)))
colnames(table_timing) = c("Using", "Run-time (s)", "Speed-up", "Average Measured $\\Delta J$")

table_timing[4,4] = NA # set the MSE to NA, not to be confused with actual DJ

