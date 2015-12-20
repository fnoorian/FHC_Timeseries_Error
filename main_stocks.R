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
# Synthetic test for LQ Error

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
num_total_tests = 10

######## the problem definition
x0 = 0

A_coeff = matrix(1, 1, 1)
B_coeff = matrix(1, 1, 1)
C_coeff = matrix(1, 1, 1)

# stage costs for a horizon of length n
P <- function(n) diag(n)
Q <- function(n) diag(n)

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

####### The AR model
ar.alpha = c(2.75771015, -3.13187110,  1.78730851, -0.49586865,  0.05065899)

v_collection = list()
for (i in 1:num_total_tests) {
  v_collection[[i]] = arima.sim(n = len_total_data, list(ar=ar.alpha))
}
#########################################################
## The Theta and its compressed version
horizons = c(rep.int(5,5), 5:1) # the horizons are 5 receding and then 5 shrinking

theta = Theta_matrix(horizons)

eigv = eigen(theta)$value
L = head(which(cumsum(eigv)/sum(eigv) > 0.99), 1) # only keep 99% of energy
w = compress_Theta(theta, L)

omega = Omega_matrix(horizons)

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

#########################################################
# The cost functions
MSE <- function(e) sum(e^2)/length(e) # by some unkown reason mean(x) is slower than sum(x)/length(x)
MAE <- function(e) mean(abs(e))
DeltaJ <- function(e, test_number, range_num) sum((e %*% w)^2) + e %*% omega_collection[[test_number]][[range_num]]

#########################################################
# model selection
model_range = 2:8

model_list = NULL
  
for (i1 in model_range) {
  for (i2 in model_range) {
    for (i3 in model_range) {
      for (i4 in model_range) {
        for (i5 in model_range) {
          model_list[[length(model_list) + 1]] = c(i1,i2,i3,i4,i5)
        }
      }
    }
  }
}

#########################################################
## Predict and store the prediction error vectors
e1_list = list()
e2_list = list()
e3_list = list()
e4_list = list()

for (order_list in model_list) {
  cat("Testing Order:", order_list, "\n")
  
  for (test_number in 1:num_total_tests) {
    v = v_collection[[test_number]]

    e1_list[[length(e1_list) + 1]] = predict_horizons_error(order_list, v[train_range1], v[test_range1])
    e2_list[[length(e2_list) + 1]] = predict_horizons_error(order_list, v[train_range2], v[test_range2])
    e3_list[[length(e3_list) + 1]] = predict_horizons_error(order_list, v[train_range3], v[test_range3])
    e4_list[[length(e4_list) + 1]] = predict_horizons_error(order_list, v[train_range4], v[test_range4])
  }
}

save(e1_list, e2_list, e3_list, e4_list, file=file.path(OUTPUT.DIR, "stock_errors.RData"))

#########################################################
models_errors = NULL

i = 1
for (order_list in model_list) {
  e_mse = 0
  e_mae = 0
  e_rhc = 0
  
  for (test_number in 1:num_total_tests) {
    e1 = e1_list[[i]]
    e2 = e2_list[[i]]
    i = i + 1
    
    e_mse = e_mse + (MSE(e1) + MSE(e2))/2
    e_mae = e_mae + (MAE(e1) + MAE(e2))/2
    e_rhc = e_rhc + (DeltaJ(e1, test_number, 1) + DeltaJ(e2, test_number, 2))/2
  }
  
  ec = c(e_mse / num_total_tests, e_rhc / num_total_tests, e_mae / num_total_tests, order_list)
  models_errors = rbind(models_errors, ec)
}

#########################################################
# find the best model
best_mse = which.min(models_errors[,1])
best_rhc = which.min(models_errors[,2])
best_mae = which.min(models_errors[,3])

pdf_create("MSE_VS_RHC.pdf", image.dir=OUTPUT.DIR)
plot(models_errors[,1], models_errors[,2], ylab=expression(Delta * J), xlab="MSE")
abline(v=models_errors[best_mse,1], lty=2)
abline(h=models_errors[best_rhc,2], lty=2)
dev.off()

pdf_create("MAE_VS_RHC.pdf", image.dir=OUTPUT.DIR)
plot(models_errors[,3], models_errors[,2], ylab=expression(Delta * J), xlab="MAE")
abline(v=models_errors[best_mae,3], lty=2)
abline(h=models_errors[best_rhc,2], lty=2)
dev.off()

# simplify plots by using raster graphics
png_create("MSE_VS_RHC.png", image.dir=OUTPUT.DIR)
plot(models_errors[,1], models_errors[,2], ylab=expression(Delta * J), xlab="MSE")
abline(v=models_errors[best_mse,1], lty=2)
abline(h=models_errors[best_rhc,2], lty=2)
dev.off()

#########################################################
# out-sample test

print("Testing Out-sample values")

orders_mae = models_errors[best_mae, 3 + 1:5]
orders_mse = models_errors[best_mse, 3 + 1:5]
orders_rhc = models_errors[best_rhc, 3 + 1:5]

cost_mse = NULL
cost_mae = NULL
cost_rhc = NULL
mse_mae = NULL
mse_mse = NULL
mse_rhc = NULL
mae_mse = NULL
mae_rhc = NULL
mae_mae = NULL

for (test_number in 1:num_total_tests) {
  v = v_collection[[test_number]]
  
  e1_mse = predict_horizons_error(orders_mse, v[train_range3], v[test_range3])
  e2_mse = predict_horizons_error(orders_mse, v[train_range4], v[test_range4])
  
  e1_rhc = predict_horizons_error(orders_rhc, v[train_range3], v[test_range3])
  e2_rhc = predict_horizons_error(orders_rhc, v[train_range4], v[test_range4])

  e1_mae = predict_horizons_error(orders_mae, v[train_range3], v[test_range3])
  e2_mae = predict_horizons_error(orders_mae, v[train_range4], v[test_range4])
  
  cost_mse = rbind(cost_mse, (DeltaJ(e1_mse, test_number, 3) + DeltaJ(e2_mse, test_number, 4))/2)
  cost_mae = rbind(cost_mae, (DeltaJ(e1_mae, test_number, 3) + DeltaJ(e2_mae, test_number, 4))/2)
  cost_rhc = rbind(cost_rhc, (DeltaJ(e1_rhc, test_number, 3) + DeltaJ(e2_rhc, test_number, 4))/2)
  
  mse_mse =  rbind(mse_mse, (MSE(e1_mse) + MSE(e2_mse))/2)
  mse_mae =  rbind(mse_mae, (MSE(e1_mae) + MSE(e2_mae))/2)
  mse_rhc =  rbind(mse_mse, (MSE(e1_rhc) + MSE(e2_rhc))/2)
  
  mae_mse =  rbind(mae_mse, (MAE(e1_mse) + MAE(e2_mse))/2)
  mae_mae =  rbind(mae_mae, (MAE(e1_mae) + MAE(e2_mae))/2)
  mae_rhc =  rbind(mae_rhc, (MAE(e1_rhc) + MAE(e2_rhc))/2)
}
#########################################################
# Find the best AIC model
#  If there is no overfitting in the MSE tests, 
#  AIC and MSE should give identical results
print("Testing AIC Models")

cost_aic_i = NULL
mse_aic_i = NULL
mae_aic_i = NULL
cost_aic_o = NULL
mse_aic_o = NULL
mae_aic_o = NULL

build_training_matrix <- function(ts, h, order) {
  ee = embed(ts, order + h)
  
  colnames(ee) = c(paste0("h", 1:h), paste0("x", 1:order))
  colnames(ee)[1] = "target"
  
  if (h >= 2) {
    ee = ee[,-c(2:h)]
  }
  
  return(ee)
}

for (test_num in 1:num_total_tests) {
  v = v_collection[[test_num]]
  v_train = v[1:80]
  
  # use MSE error to select model orders
  orders_aic = sapply(1:5, function(h) {
    aic = sapply(2:8, function(i) {
      ee = build_training_matrix(v_train, h, i)
      AIC(lm(target ~ 0+., as.data.frame(ee)))
    })
    model_range[which.min(aic)]
  })
  
  # use the best order to predict all
  e1_aic = predict_horizons_error(orders_aic, v[train_range1], v[test_range1])
  e2_aic = predict_horizons_error(orders_aic, v[train_range2], v[test_range2])
  e3_aic = predict_horizons_error(orders_aic, v[train_range3], v[test_range3])
  e4_aic = predict_horizons_error(orders_aic, v[train_range4], v[test_range4])
  
  # measure in-sample & outsample
  cost_aic_i = rbind(cost_aic_i, (DeltaJ(e1_aic, test_num, 1) + DeltaJ(e2_aic, test_num, 2))/2)
  mse_aic_i =  rbind(mse_aic_i, (MSE(e1_aic) + MSE(e2_aic))/2)
  mae_aic_i =  rbind(mae_aic_i, (MAE(e1_aic) + MAE(e2_aic))/2)

  cost_aic_o = rbind(cost_aic_o, (DeltaJ(e3_aic, test_num, 3) + DeltaJ(e4_aic, test_num, 4))/2)
  mse_aic_o =  rbind(mse_aic_o, (MSE(e3_aic) + MSE(e4_aic))/2)
  mae_aic_o =  rbind(mae_aic_o, (MAE(e3_aic) + MAE(e4_aic))/2)
}

#########################################################
# Save raw results
save(models_errors, best_mae, best_mse, best_rhc, 
     cost_mse, cost_mae, cost_rhc, mse_mae, mse_mse, mse_rhc, mae_mse, mae_rhc, mae_mae,
     file=file.path(OUTPUT.DIR, "stocks_performance.RData"))

#########################################################
# the covariance
cov_partition = models_errors[,1:3]
colnames(cov_partition)[1:3]=c("MSE", "DeltaJ", "MAE")

#########################################################
# Create Error tables
#ET <- function(err_vector) paste(signif(mean(err_vector), 4), "\\pm", signif(sd(err_vector), 4))
ET <- function(err_vector) round(mean(err_vector), 4)

# prepare latex output
perform_table = rbind(
  data.frame(A = "AIC",  
             B = ET(mse_aic_i), 
             C = ET(cost_aic_i), 
             D = ET(mse_aic_o),
             E = ET(cost_aic_o)),
  data.frame(A = "MSE",  
              B = models_errors[best_mse,1], 
              C = models_errors[best_mse,2], 
              D = ET(mse_mse),
              E = ET(cost_mse)),
   data.frame(A = "$\\Delta J$",  
              B = models_errors[best_rhc,1], 
              C = models_errors[best_rhc,2], 
              D = ET(mse_rhc),
              E = ET(cost_rhc)))
colnames(perform_table) = c("Measure for selection", "In-sample MSE", "In-sample $\\Delta J$", "Out-of-sample MSE" , "Out-of-sample $\\Delta J$")

ratio_dj_over_mse = (mean(cost_mse) - mean(cost_rhc)) / mean(cost_mse)
#########################################################
# save results

perform_table_full = perform_table
perform_table_full[,2:5] = bold_min_table_cols(perform_table_full[,2:5], digits=4)

sink(file.path(OUTPUT.DIR, "stocks_results.txt"))
options(width=150)
print(perform_table_full)

cat("\n")
print("Correlation ceofficient of errors:")
print(cor(cov_partition))

cat("\n")
print("Cost Improvement of DeltaJ over MSE (in %):")
print(signif(ratio_dj_over_mse * 100, 3))

cat("\n")
cat("MSE and DeltaJ out-of-sample results t-test:")
print(t.test(cost_mse, cost_rhc, paired=TRUE))

sink()

perform_table_paper = perform_table
perform_table_paper = perform_table_paper[-1,] # drop the AIC for paper
perform_table_paper[,2:5] = bold_min_table_cols(perform_table_paper[,2:5], digits=4)

save_tex_table(perform_table_paper, "table_stocks_performance.tex", save.dir = OUTPUT.DIR, 
               floating=FALSE, 
               digits = 4)

#########################################################
#The timing experiments
print("Timing started...")

source("timing_stocks.R")

save(t2, t3, t4, table_timing, file=file.path(OUTPUT.DIR, "stocks_timings.RData"))

save_tex_table(table_timing, "table_stocks_timings.tex", save.dir = OUTPUT.DIR,
             label = "tab:lqr_stocks_timings")
