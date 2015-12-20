################################################################################
# Copyright (c) 2015 Farzad Noorian. Distributed under Apache License v2.0
################################################################################
# simulating and measuring the Times

# use this to accurately measure short running tests
REPEATS_FOR_ACCURACY_DJ_FULL = 10
REPEATS_FOR_ACCURACY_DJ_APPROX = 100 
REPEATS_FOR_ACCURACY_MSE = 100 

#######
## use Full matrix for ranking errors

t2 = system.time({
  for (count in 1:REPEATS_FOR_ACCURACY_DJ_FULL) {
    e_rhc_full_list = NULL
    i = 1
    for (order_list in model_list) {
      
      e_rhc = 0
      for (test_number in 1:num_total_tests) {
        V1 = v_collection[[test_number]][test_range1]
        V2 = v_collection[[test_number]][test_range2]
        e1 = e1_list[[i]]
        e2 = e2_list[[i]]
        i = i + 1
        
        e_rhc = e_rhc + (t(e1) %*% theta %*% e1 + t(e2) %*% theta %*% e2 + 
                           t(e1) %*% omega %*% V1 + t(e2) %*% omega %*% V2)/2
      }
      
      e_rhc_full_list = c(e_rhc_full_list, e_rhc / num_total_tests)
    }
  }
})


t2 = t2 / REPEATS_FOR_ACCURACY_DJ_FULL

#######
## use simulator to rank errors
N = length(horizons)

cost_U <- function(U, V) t(U) %*% JA(N) %*% U + 2 * t(U) %*% JB(N) %*% V + t(V) %*% JC(N) %*% V1

# first collect the actual minimum costs
min_cost1 = NULL
min_cost2 = NULL
for (test_number in 1:num_total_tests) {
  V1 = v_collection[[test_number]][test_range1]
  V2 = v_collection[[test_number]][test_range2]
  
  U1 = RHC_simulator(V1, e1 * 0, horizons)$U
  U2 = RHC_simulator(V2, e2 * 0, horizons)$U
  
  min_cost1 = c(min_cost1, cost_U(U1, V1))
  min_cost2 = c(min_cost2, cost_U(U2, V2))
}

# now compute delta J using "cost - min_cost"
e_simulator_list = NULL
i = 1
t4 = system.time({
  for (order_list in model_list) {
      
    e_sim = 0
    for (test_number in 1:num_total_tests) {
      e1 = e1_list[[i]]
      e2 = e2_list[[i]]
      i = i + 1

      V1 = v_collection[[test_number]][test_range1]
      V2 = v_collection[[test_number]][test_range2]
      U1 = RHC_simulator(V1, e1, horizons)$U 
      U2 = RHC_simulator(V2, e2, horizons)$U 
          
      
      e_sim = e_sim + (cost_U(U1, V1) - min_cost1[test_number] + 
                         cost_U(U2, V2) - min_cost2[test_number])/2
    }
    
    e_simulator_list = c(e_simulator_list, e_sim / num_total_tests)
  }
})

#######
## use approximated matrix for ranking Errors
t3 = system.time({
  for (count in 1:REPEATS_FOR_ACCURACY_DJ_APPROX) {
    e_rhc_apprx_list = NULL
    i = 1
    
    for (order_list in model_list) {
      
      e_rhc = 0
      for (test_number in 1:num_total_tests) {
        e1 = e1_list[[i]]
        e2 = e2_list[[i]]
        i = i + 1
        
        e_rhc = e_rhc + (DeltaJ(e1, test_number, 1) + DeltaJ(e2, test_number, 2))/2
      }
      
      e_rhc_apprx_list = c(e_rhc_apprx_list, e_rhc / num_total_tests)
    }
  }
})
t3 = t3 / REPEATS_FOR_ACCURACY_DJ_APPROX

#######
## MSE error timing
t5 = system.time({
  for (count in 1: REPEATS_FOR_ACCURACY_MSE) {
    e_mse_full_list = NULL
    i = 1
    
    for (order_list in model_list) {
      e_mse = 0
      for (test_number in 1:num_total_tests) {
        e1 = e1_list[[i]]
        e2 = e2_list[[i]]
        i = i + 1
        
        e_mse = e_mse + (MSE(e1) + MSE(e2))/2
      }
      
      e_mse_full_list = c(e_mse_full_list, e_rhc / num_total_tests)
    }
  }
})

t5 = t5 / REPEATS_FOR_ACCURACY_MSE

#######
# create the table
table_timing = rbind(
  data.frame(names = "Step-by-step simulation", time = round(t4[1], 3), speedup = round(t4[1] / t4[1], 0) , computed_cost = mean(e_simulator_list)),
  data.frame(names = "Closed form $\\Delta J$", time = round(t2[1], 3), speedup = round(t4[1] / t2[1], 1) , computed_cost = mean(e_rhc_full_list)),
  data.frame(names = "Approximated $\\Delta J$", time = round(t3[1], 3), speedup = round(t4[1] / t3[1], 1) , computed_cost = mean(e_rhc_full_list)),
  data.frame(names = "MSE", time = round(t5[1], 3), speedup = round(t4[1] / t5[1], 0) , computed_cost = mean(e_mse_full_list)))
colnames(table_timing) = c("Using", "Run-time (s)", "Speed-up", "Measured $\\Delta J$")

table_timing[4,4] = NA # set the MSE to NA, not to be confused with actual DJ

