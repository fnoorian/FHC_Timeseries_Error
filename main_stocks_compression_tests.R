################################################################################
# (c) Copyright 2015 Farzad Noorian. 
#
# Licensed under the Apache License, Version 2.0 (the "License"); 
# you may not use this file except in compliance with the License. 
# You may obtain a copy of the License at 
# 
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software 
# distributed under the License is distributed on an "AS IS" BASIS, 
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
# See the License for the specific language governing permissions and 
# limitations under the License.
################################################################################
# Testing the effectiveness matrix compression in getting results
# Reuses stocks portfolio problem structure

rm(list=ls())

source("formatting/curve_plot.R")
source("formatting/print_table.R")
OUTPUT.DIR = "Output"

source("utils/matrix_tools.R")
source("fhc/mpc_matrices.R")
source("fhc/theta.R")
source("fhc/simulator.R")
source("prediction/linear_predictor.R")

######## Experiment repeatition settings
set.seed(1) # for reproducibility
num_total_tests = 20

######## the problem definition
x0 = 0

A_coeff = matrix(1, 1, 1)
B_coeff = matrix(1, 1, 1)
C_coeff = matrix(1, 1, 1)

# stage costs for a horizon of length n
P <- function(n) diag(n)
Q <- function(n) diag(n)

horizons = c(rep.int(5,5), 5:1) # the horizons are 5 receding and then 5 shrinking

####### Train and test range
len_total_data = 100

train_range1 = 1:60   # validation 1
test_range1 = 61:70

train_range2 = 1:70   # validation 2
test_range2 = 71:80

train_range3 = 1:80   # test 1
test_range3 = 81:90

train_range4 = 1:90   # test 2
test_range4 = 91:100

####### The AR model and its prediction errors
ar.alpha = c(2.75771015, -3.13187110,  1.78730851, -0.49586865,  0.05065899)

v_collection = list()
for (i in 1:num_total_tests) {
  v_collection[[i]] = arima.sim(n = len_total_data, list(ar=ar.alpha))
}

# precompute errors for only a single model
model_orders = c(4,4,4,4,4)

error_collection = list()
for (test_number in 1:num_total_tests) {
  v = v_collection[[test_number]]
  
  err1 = predict_horizons_error(model_orders, v[train_range1], v[test_range1]) 
  err2 = predict_horizons_error(model_orders, v[train_range2], v[test_range2])
  err3 = predict_horizons_error(model_orders, v[train_range3], v[test_range3]) 
  err4 = predict_horizons_error(model_orders, v[train_range4], v[test_range4])
  
  error_collection[[test_number]] = list()
  error_collection[[test_number]][[1]] = err1
  error_collection[[test_number]][[2]] = err2
  error_collection[[test_number]][[3]] = err3
  error_collection[[test_number]][[4]] = err4
}


#########################################################
## The Theta and its compressed version
theta = Theta_matrix(horizons)
omega = Omega_matrix(horizons)

eigv = eigen(theta)$value

lambda_performance = NULL
for (L in c(2:10, 20, 40)) {
  
  print(paste("Testing L = ", L))
  
  # Lambda is the energy
  lambda = sum(eigv[1:L])/sum(eigv)

  # precompute omegas
  omega_collection = list()
  for (i in 1:num_total_tests) {
    v = v_collection[[i]]
    
    omega_v = list()
    omega_v[[1]] = (omega %*% v[test_range1])
    omega_v[[2]] = (omega %*% v[test_range2])
    omega_v[[3]] = (omega %*% v[test_range3])
    omega_v[[4]] = (omega %*% v[test_range4])
    
    omega_collection[[i]] = omega_v
  }
  
  if (L < 40) {
    # W is the compressed Theta
    w = compress_Theta(theta, L)
    
    # The error function
    DeltaJ <- function(e, test_number, range_num) sum((e %*% w)^2) + e %*% omega_collection[[test_number]][[range_num]]
  } else {
    # 
    DeltaJ <- function(e, test_number, range_num) sum(t(e) %*% theta %*% e) + e %*% omega_collection[[test_number]][[range_num]]
  }
  
  exec_time = system.time({
    for (repeater in 1:10000) {
      err = 0
      for (test_number in 1:num_total_tests) {
        err3 = error_collection[[test_number]][[3]]
        err4 = error_collection[[test_number]][[4]]

        err = err + (DeltaJ(err3, test_number, 3) + DeltaJ(err4, test_number, 4))/2
      }
      err = err / num_total_tests
    }
  })

  # collect the performace
  lambda_performance = rbind(lambda_performance,
                             data.frame(L = L, lambda = lambda, Runtime = exec_time[1], DeltaJ = err))
}

#########################################################
## Create Error Table

# add error percent
lambda_performance = cbind(lambda_performance[,1:3], 
                           SpeedUp = tail(lambda_performance$Runtime, 1) / lambda_performance$Runtime,
                           lambda_performance$DeltaJ,
                           DeltaJ_Percent = (lambda_performance$DeltaJ/tail(lambda_performance$DeltaJ, 1)) * 100)

colnames(lambda_performance) = c("L", "$\\lambda$", "Run-time (s)", "Speed-up", "Measured $\\Delta J$", "$\\Delta J$ Accuracy \\%")

#########################################################
## Save results
print(lambda_performance)

# some formatting fixes
performance_table = round(lambda_performance, 4)
performance_table[,1] = as.character(performance_table[,1])
performance_table[performance_table == 1] = "1"
performance_table[performance_table == 100] = "100"

save_tex_table(performance_table, "table_theta_compression.tex", save.dir = OUTPUT.DIR,
               floating = FALSE,
               align="ccccccc")
