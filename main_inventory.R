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
# LQ Error Test for inventory management

rm(list=ls())

source("formatting/curve_plot.R")
source("formatting/print_table.R")
OUTPUT.DIR = "Output"

source("utils/matrix_tools.R")
source("fhc/mpc_matrices.R")
source("fhc/theta.R")
source("fhc/simulator.R")
source("prediction/linear_predictor_power.R", chdir=TRUE)

set.seed(1) # for reproducibility

######## the problem definition
N = 10 # length of simulation
h = 3

x0 = rep(0, h + 1)
A_coeff = upper.diagonal(h + 1)

B_coeff = diag(h)
B_coeff = rbind(0, B_coeff)

C_coeff = matrix(0, h+1, 1)
C_coeff[1, 1] = -1

gamma = 0.7
spot_price = 4

Q <- function(N) diag(rep(c(spot_price, 0, 0, 0), N)) # state weights
P <- function(N) diag(rep(c(gamma ^ 3 * spot_price, gamma ^ 2 * spot_price, gamma * spot_price), N)) # inputs weights

####### Train and test range
num_total_tests = 100
len_total_data = 100

train_range1 = 1:60   # validation 1
test_range1 = 61:70

train_range2 = 1:70   # validation 2
test_range2 = 71:80

train_range3 = 1:80   # test 1
test_range3 = 81:90

train_range4 = 1:90   # test 2
test_range4 = 91:100

####### The AR model
ar.alpha = c( 2.75771015, -3.13187110,  1.78730851, -0.49586865,  0.05065899)

v_collection = list()
for (i in 1:num_total_tests) {
  v_arima = arima.sim(n = len_total_data, list(ar=ar.alpha), rand.gen = function(n, ...) rnorm(n, sd=0.25))
  v_collection[[i]] = exp(v_arima/10)
}

#########################################################
## The Theta and its compressed version
K = h + 1 # must include current time as well
horizons = c(rep.int(K, N-K), K:1)

theta = Theta_matrix(horizons)
omega = Omega_matrix(horizons)

stopifnot(sum(abs(theta - diag(diag(theta)))) < 1e-10)
stopifnot(sum(abs(omega)) < 1e-10)

dtheta = diag(theta) # remove the obvious zeros (for current time-step)

# the first index per prediction (the obvious zeros)
i = cumsum(c(1, head(horizons, -1)))
stopifnot(sum(abs(theta[i])) < 1e-10)
#theta = theta[-i]

#########################################################
# The cost functions
MSE <- function(e) sum(e^2)/length(e) # by some unkown reason mean(x) is slower than sum(x)/length(x)
MAE <- function(e) mean(abs(e))
DeltaJ <- function(e) sum(e^2 * dtheta)

#########################################################
# The model orders
orders_range = 1:8

direct_orders = list()

for (l in orders_range) {
  for (m in orders_range) {
    for (n in orders_range) {
      i = length(direct_orders)
      direct_orders[[i + 1]] = c(l, m, n)
    }
  }
}

#########################################################
# get model errors
MAX_iterate = 16

#######
print("Exhaustive + Random (Hybrid) Models")

hybrid_dj_insample = 0
hybrid_dj_outsample = 0
hybrid_mae_insample = 0
hybrid_mae_outsample = 0
hybrid_mse_insample = 0
hybrid_mse_outsample = 0

for (test_num in 1:num_total_tests) {
  print(paste("Hybrid test ", test_num))
  v = v_collection[[test_num]]
  
  # the exhaustive part, find the best first order
  index_first_element = head(cumsum(c(1, head(horizons, -1))) + 1, -1)
  
  err_all = NULL
  for (l in orders_range) {
    model = c(l, min(orders_range), min(orders_range))
    
    pred1 = predict_horizons_error_power(model, v[train_range1], v[test_range1])
    pred2 = predict_horizons_error_power(model, v[train_range2], v[test_range2])
    
    e = MSE(pred1[index_first_element]) + MSE(pred2[index_first_element])
    err_all = c(err_all, e)
  }

  best_index1 = which.min(err_all)

  # the random search part, find the best second and third
  model_err = list()
  
  remaining_iterations = MAX_iterate - length(orders_range)
  for (i in 1:remaining_iterations) {
    # select random model
    model = c(best_index1, sample(orders_range, 1), min(orders_range))
    
    pred1 = predict_horizons_error_power(model, v[train_range1], v[test_range1])
    pred2 = predict_horizons_error_power(model, v[train_range2], v[test_range2])
    
    model_err[[i]] = list(err = DeltaJ(pred1) + DeltaJ(pred2),
                          order = model)
  }
  
  best_order = model_err[[which.min(sapply(model_err, function(e) e$err))]]$order
  
  # use the best order to predict all, measure in-sample & outsample
  pred_in_1 = predict_horizons_error_power(model, v[train_range1], v[test_range1])
  pred_in_2 = predict_horizons_error_power(model, v[train_range2], v[test_range2])
  pred_out_1 = predict_horizons_error_power(model, v[train_range3], v[test_range3])
  pred_out_2 = predict_horizons_error_power(model, v[train_range4], v[test_range4])
  
  hybrid_dj_insample = hybrid_dj_insample + DeltaJ(pred_in_1) + DeltaJ(pred_in_2)
  hybrid_dj_outsample = hybrid_dj_outsample + DeltaJ(pred_out_1) + DeltaJ(pred_out_2)
  hybrid_mae_insample = hybrid_mae_insample + MAE(pred_in_1) + MAE(pred_in_2)
  hybrid_mae_outsample = hybrid_mae_outsample + MAE(pred_out_1) + MAE(pred_out_2)
  hybrid_mse_insample = hybrid_mse_insample + MSE(pred_in_1) + MSE(pred_in_2)
  hybrid_mse_outsample = hybrid_mse_outsample + MSE(pred_out_1) + MSE(pred_out_2)
}

#############
print("Direct Models")

direct_dj_insample = 0
direct_dj_outsample = 0
direct_mae_insample = 0
direct_mae_outsample = 0
direct_mse_insample = 0
direct_mse_outsample = 0

for (test_num in 1:num_total_tests) {
  print(paste("direct test ", test_num))
  v = v_collection[[test_num]]
  
  # the random search part, find the best model
  model_err = list()
  
  for (i in 1:MAX_iterate) {
    # select random model
    model = sample(direct_orders, 1)[[1]]
    
    pred1 = predict_horizons_error_power(model, v[train_range1], v[test_range1])
    pred2 = predict_horizons_error_power(model, v[train_range2], v[test_range2])
    
    model_err[[i]] = list(err = DeltaJ(pred1) + DeltaJ(pred2),
                          order = model)
  }
  
  best_order = model_err[[which.min(sapply(model_err, function(e) e$err))]]$order
  
  # use the best order to predict all, measure in-sample & outsample
  pred_in_1 = predict_horizons_error_power(model, v[train_range1], v[test_range1])
  pred_in_2 = predict_horizons_error_power(model, v[train_range2], v[test_range2])
  pred_out_1 = predict_horizons_error_power(model, v[train_range3], v[test_range3])
  pred_out_2 = predict_horizons_error_power(model, v[train_range4], v[test_range4])
  
  direct_dj_insample = direct_dj_insample + DeltaJ(pred_in_1) + DeltaJ(pred_in_2)
  direct_dj_outsample = direct_dj_outsample + DeltaJ(pred_out_1) + DeltaJ(pred_out_2)
  direct_mae_insample = direct_mae_insample + MAE(pred_in_1) + MAE(pred_in_2)
  direct_mae_outsample = direct_mae_outsample + MAE(pred_out_1) + MAE(pred_out_2)
  direct_mse_insample = direct_mse_insample + MSE(pred_in_1) + MSE(pred_in_2)
  direct_mse_outsample = direct_mse_outsample + MSE(pred_out_1) + MSE(pred_out_2)
}

err_table = rbind(c(direct_mae_insample, direct_mse_insample, direct_dj_insample, direct_mae_outsample, direct_mse_outsample, direct_dj_outsample),
                  c(hybrid_mae_insample, hybrid_mse_insample, hybrid_dj_insample, hybrid_mae_outsample, hybrid_mse_outsample, hybrid_dj_outsample))
err_table = err_table / num_total_tests # the averaging
err_table = data.frame(err_table)
err_table = cbind(c("Random", "Hybrid"), err_table)

colnames(err_table) = c("Search Method", "In-sample MAE", "In-sample MSE", "In-sample $\\Delta J$",  "Out-sample MAE", "Out-sample MSE" , "Out-sample $\\Delta J$")

# remove useless columns
err_table = err_table[,-c(2,5)]

save_tex_table(err_table, "table_preorder_performance.tex", save.dir = OUTPUT.DIR,
               label = "tab:lqr_preorder_performance")

# print results
print(err_table)

###############################################
# timing simulator
# fill the errors list
e1_list = list()
e2_list = list()

t1 = system.time({
  for (model in direct_orders) {
    print(model)
    e1 = lapply(v_collection[1:num_total_tests], function(v) predict_horizons_error_power(model, v[train_range1], v[test_range1]))
    e2 = lapply(v_collection[1:num_total_tests], function(v) predict_horizons_error_power(model, v[train_range2], v[test_range2]))
    
    e1_list = c(e1_list, e1)
    e2_list = c(e2_list, e2)
  }
})

model_list = direct_orders
source("timing_inventory.R")

save_tex_table(table_timing, "table_preorder_timings.tex", save.dir = OUTPUT.DIR,
             label = "tab:lqr_preorder_timings")

print(table_timing)
