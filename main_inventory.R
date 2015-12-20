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

library("Mcomp")  # using data from the MComp

source("formatting/curve_plot.R")
source("formatting/print_table.R")
OUTPUT.DIR = "Output"

source("utils/matrix_tools.R")
source("fhc/mpc_matrices.R")
source("fhc/theta.R")
source("fhc/simulator.R")
source("prediction/linear_predictor_inventory.R", chdir=TRUE)

######## Experiment repeatition settings
set.seed(1) # for reproducibility
num_total_tests = 1020

######## the problem definition
N = 10 # length of simulation
n = 3  # max length of horizon

kappa = n # size of preorders selected as the max length of horizon
x0 = rep(0, kappa + 1) # 
A_coeff = upper.diagonal(kappa + 1)

B_coeff = diag(kappa)
B_coeff = rbind(0, B_coeff)

C_coeff = matrix(0, kappa+1, 1)
C_coeff[1, 1] = -1

gamma = 0.7
spot_price = 4

# P and Q for custom horizon length n
Q <- function(n) diag(rep(c(spot_price, 0, 0, 0), n)) # state weights
P <- function(n) diag(rep(c(gamma ^ 3 * spot_price, gamma ^ 2 * spot_price, gamma * spot_price), n)) # inputs weights

####### Train and test range
len_total_data = 100

train_ranges = list(list(train = 1:50, test = 51:60),
                    list(train = 1:55, test = 56:65),
                    list(train = 1:60, test = 61:70),
                    list(train = 1:65, test = 66:75),
                    list(train = 1:70, test = 71:80))

validation_ranges = list(list(train = 1:80, test = 81:90),
                         list(train = 1:85, test = 86:95),
                         list(train = 1:90, test = 91:100))

####### The AR model
v_collection = list()
counter = 1
for (i in 1:num_total_tests) {
  m3_data = M3[[counter]]$x
  
  while ((length(m3_data) < len_total_data) | (any(m3_data < 0))) {
    counter = counter + 1
    m3_data = M3[[counter]]$x
  }
  
  print(paste("Using: ", M3[[counter]]$st, counter))
  v_collection[[i]] = m3_data
  counter = counter + 1
}

#########################################################
## The Theta and its compressed version
horizons = c(rep.int(n + 1, N-(n+1)), (n + 1):1) # using n+1, as it must include the current time as well

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
# Returns a list of prediction errors for the ts_data using linear model orders specified
get_prediction_errors <- function(model, ts_data, ts_range) {
  lapply(ts_range, function(rng) {
    predict_inventory_error(model, ts_data[rng$train], ts_data[rng$test])
  })
}

# The in-sample and out-of-sample error function
get_in_out_sample_errs <- function(ts_data, model) {

  # measure in-sample errors
  pred_errs = get_prediction_errors(model, ts_data, train_ranges)
  dj_insample = mean(sapply(pred_errs, DeltaJ))
  mae_insample = mean(sapply(pred_errs, MAE))
  mse_insample = mean(sapply(pred_errs, MSE))
  
  # measure out-of-sample errors
  pred_errs = get_prediction_errors(model, ts_data, validation_ranges)
  dj_outsample = mean(sapply(pred_errs, DeltaJ))
  mae_outsample = mean(sapply(pred_errs, MAE))
  mse_outsample = mean(sapply(pred_errs, MSE))

  errs = c(mae_insample, mse_insample, dj_insample, mae_outsample, mse_outsample, dj_outsample)
  names(errs) = c("In-sample MAE", "In-sample MSE", "In-sample $\\Delta J$",  "Out-sample MAE", "Out-sample MSE" , "Out-sample $\\Delta J$")
  return(errs)
}

get_in_out_sample_naive_errs <- function(ts_data) {

  # measure in-sample errors
  pred_errs =   lapply(train_ranges, function(rng) {
    predict_inventory_naive_error(ts_data[rng$train], ts_data[rng$test], 3)
  })
    
  dj_insample = mean(sapply(pred_errs, DeltaJ))
  mae_insample = mean(sapply(pred_errs, MAE))
  mse_insample = mean(sapply(pred_errs, MSE))
  
  # measure out-of-sample errors
  pred_errs =   lapply(validation_ranges, function(rng) {
    predict_inventory_naive_error(ts_data[rng$train], ts_data[rng$test], 3)
  })

  dj_outsample = mean(sapply(pred_errs, DeltaJ))
  mae_outsample = mean(sapply(pred_errs, MAE))
  mse_outsample = mean(sapply(pred_errs, MSE))
  
  errs = c(mae_insample, mse_insample, dj_insample, mae_outsample, mse_outsample, dj_outsample)
  names(errs) = c("In-sample MAE", "In-sample MSE", "In-sample $\\Delta J$",  "Out-sample MAE", "Out-sample MSE" , "Out-sample $\\Delta J$")
  return(errs)
}

#########################################################
# The model orders
orders_range = 1:8

all_model_orders = list()

for (i1 in orders_range) {
  for (i2 in orders_range) {
    for (i3 in orders_range) {
      i = length(all_model_orders)
      all_model_orders[[i + 1]] = c(i1, i2, i3)
    }
  }
}

#########################################################
# index of the first, second, and third element in the prediction vector output (t, t+1, t+2 for each theta)
index_first_elements = head(cumsum(c(1, head(horizons, -1))) + 1, -1)
index_second_elements = head(index_first_elements + 1, -1)
index_third_elements = head(index_second_elements + 1, -1)

# the number of maximum iterations for random and hybrid search
MAX_iterate = 16

#######
# Model selection using AIC
# Notice that AIC is equal to leave-one-out cross-validation using MSE
# and use 3 * 8 evaluations (unlike random search which was limited to 16)
# so better results are expected

print("All AIC Models")

AIC_errs = 0

build_training_matrix <- function(ts, h, order) {
  ee = embed(ts, order + h)

  colnames(ee) = c(paste0("h", 1:h), paste0("x", 1:order))
  colnames(ee)[1] = "target"

  if (h >= 2) {
    ee = ee[,-c(2:h)]
  }

  return(ee)
}

AIC_errs = NULL
for (test_num in 1:num_total_tests) {
  print(paste("AIC test ", test_num))
  v = v_collection[[test_num]]
  v_train = v[1:80]

  # use MSE error to select model orders
  best_orders = sapply(1:3, function(h) {
    aic = sapply(orders_range, function(i) {
      ee = build_training_matrix(v_train, h, i)
      AIC(lm(target ~ 0+., as.data.frame(ee)))
    })
    orders_range[which.min(aic)]
  })

  # use the best order to predict all, measure in-sample & outsample
  AIC_errs = rbind(AIC_errs, get_in_out_sample_errs(v, best_orders))
}

#######
print("Hybrid (Exhaustive + Random) Models")

hybrid_errs = NULL

for (test_num in 1:num_total_tests) {
  print(paste("Hybrid test ", test_num))
  v = v_collection[[test_num]]
  
  # the exhaustive part, find the best first order
  err_all = NULL
  for (l in orders_range) {
    model = c(l, min(orders_range), min(orders_range))
    
    pred_errs = get_prediction_errors(model, v, train_ranges)
    err = sapply(pred_errs, function(e) MSE(e[index_first_elements]))

    err_all = c(err_all, mean(err))
  }

  best_index1 = which.min(err_all)

  # the random search part, find the best second and third
  model_err = list()
  
  remaining_iterations = MAX_iterate - length(orders_range)
  for (i in 1:remaining_iterations) {
    # select random model
    model = c(best_index1, sample(orders_range, 1), sample(orders_range, 1))
    
    pred_errs = get_prediction_errors(model, v, train_ranges)
    err = sapply(pred_errs, DeltaJ)

    model_err[[i]] = list(err = mean(err), order = model)
  }
  
  best_order = model_err[[which.min(sapply(model_err, function(e) e$err))]]$order

  # use the best order to predict all, measure in-sample & outsample
  hybrid_errs = rbind(hybrid_errs, get_in_out_sample_errs(v, best_order))
}

#############
print("Random Model with DeltaJ Selection")

random_dj_errs = NULL

for (test_num in 1:num_total_tests) {
  print(paste("Random DJ test ", test_num))
  v = v_collection[[test_num]]
  
  # the random search part, find the best model
  model_err = list()
  
  for (i in 1:MAX_iterate) {
    # select random model
    model = sample(all_model_orders, 1)[[1]]
    
    pred_errs = get_prediction_errors(model, v, train_ranges)
    err = sapply(pred_errs, DeltaJ)

    model_err[[i]] = list(err = mean(err), order = model)
  }
  
  best_order = model_err[[which.min(sapply(model_err, function(e) e$err))]]$order

  # use the best order to predict all, measure in-sample & outsample
  random_dj_errs = rbind(random_dj_errs, get_in_out_sample_errs(v, best_order))
}

#############
print("Random Model Selection with MSE")

random_mse_errs = NULL

for (test_num in 1:num_total_tests) {
  print(paste("Random MSE test ", test_num))
  v = v_collection[[test_num]]
  
  # the random search part, find the best model
  model_err = list()
  
  for (i in 1:MAX_iterate) {
    # select random model
    model = sample(all_model_orders, 1)[[1]]

    pred_errs = get_prediction_errors(model, v, train_ranges)
    err = sapply(pred_errs, MSE)
    
    model_err[[i]] = list(err = mean(err), order = model)
  }

  best_order = model_err[[which.min(sapply(model_err, function(e) e$err))]]$order

  # use the best order to predict all, measure in-sample & outsample
  random_mse_errs = rbind(random_mse_errs, get_in_out_sample_errs(v, best_order))
}

#############
print("Naive Models")

naive_errs = NULL

for (test_num in 1:num_total_tests) {
  print(paste("Naive test ", test_num))
  v = v_collection[[test_num]]
  
  # use the best order to predict all, measure in-sample & outsample
  naive_errs = rbind(naive_errs, get_in_out_sample_naive_errs(v))
}

#########################################################
# Save Raw Data
save(AIC_errs, random_mse_errs, random_dj_errs, 
     hybrid_errs, naive_errs, file=file.path(OUTPUT.DIR, "inventory_performance.RData"))

#########################################################
# Create Error tables
#ET <- function(errs) paste(signif(apply(errs, 2, mean), 4), "\\pm", signif(apply(errs, 2, sd), 4))
ET <- function(errs) signif(apply(errs, 2, mean), 4)
err_table_with_sd = rbind(
  data.frame(Selection.Method="AIC", t(ET(AIC_errs / naive_errs))),
  data.frame(Selection.Method="Random Search (using MSE)", t(ET(random_mse_errs / naive_errs))),
  data.frame(Selection.Method="Random Search (using $\\Delta J$)", t(ET(random_dj_errs / naive_errs))),
  data.frame(Selection.Method="Hybrid Search (using $\\Delta J$)", t(ET(hybrid_errs / naive_errs))))
colnames(err_table_with_sd) = c("Selection Method", colnames(random_mse_errs))

normalized_outsample_dj_for_random_mse = (random_mse_errs / naive_errs)[,6]
normalized_outsample_dj_for_hybrid= (hybrid_errs / naive_errs)[,6]
hybrid_random_ratio = (mean(normalized_outsample_dj_for_random_mse) - mean(normalized_outsample_dj_for_hybrid)) / mean(normalized_outsample_dj_for_random_mse)

#########################################################
# save results
bolded_err_table = err_table_with_sd
bolded_err_table[,2:7] = bold_min_table_cols(bolded_err_table[,2:7], digits=4)

sink(file.path(OUTPUT.DIR, "preorder_results.txt"))
options(width=180)
print(bolded_err_table)

cat("\n")

print("Cost Improvement of Hybrid over Random MSE (in %):")
print(round(hybrid_random_ratio * 100, 1))

cat("\n")

cat("Random MSE and Hybrid  out-of-sample results t-test: ")
print(t.test(normalized_outsample_dj_for_random_mse, normalized_outsample_dj_for_hybrid, paired=TRUE))

cat("\n")

cat("Correlation between errors: \n")
print(cor(cbind(Random.MSE=random_mse_errs[,6], Random.DeltaJ=random_dj_errs[,6], Hybrid.DeltaJ=hybrid_errs[,6])))

sink()

# remove useless columns
err_table = err_table_with_sd[-1,-c(2,5)]
err_table[,2:5] = bold_min_table_cols(err_table[,2:5], digits=4)

save_tex_table(err_table, "table_preorder_performance.tex", save.dir = OUTPUT.DIR,
               align = c("cccccc"),
               floating = FALSE)

###############################################
# timing simulator
print("Timing started...")

source("timing_inventory.R")

save_tex_table(table_timing, "table_preorder_timings.tex", save.dir = OUTPUT.DIR,
             label = "tab:lqr_preorder_timings")

print(table_timing)
