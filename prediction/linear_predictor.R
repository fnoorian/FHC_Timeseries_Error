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
# linear prediction functions

fit_beta <- function(h, order, ts) {
  # fit AR model of "order" such that it predicts "h" step ahead of "ts
  # h: horizon to predict
  # order: order of AR model
  # ts: time-series to train on
  # returns: beta vector in linear model Y = X * BETA 
  ee = embed(ts, order + h)
  x = ee[,-(1:h)]
  y = ee[,1]
  
  solve(t(x) %*% x, t(x) %*% y)
}

fit_beta_matrix <- function(order_list, ts) {
  # fit multiple AR model to predict  many step ahead of "ts
  # order_list: order of AR model to predict horizons 1 to n
  # ts: time-series to train on
  # returns: list of beta vectors in linear model Y = X * BETA
  
  beta_list = list()
  for (i in seq_along(order_list)) {
    beta_list[[i]] = fit_beta(i, order_list[i], ts)
  }
  beta_list
}

predict_beta_matrix <- function(beta_list, ts) {
  # predicts "ts" using given AR models
  # beta_list: list of BETA coefficient to use as model
  # ts: time-series to predict from
  # returns: the vector of predictions

  pred = c()
  ts = rev(ts)
  for (i in seq_along(beta_list)) {
    b = beta_list[[i]]
    x = ts[1:length(b)]
    pred = c(pred, x %*% b)
  }
  pred
}

predict_horizons_error <- function(order_list, train, test) {
  # computes error of predicting "train" against "test"
  # train: observed values
  # test: target of prediction
  
  # returns: error for each horizon, concatenated
  
  pred = c()
  actual = c()
  max_horizon = length(test)
  
  for (i in 1:length(test)) {
    
    h = min(max_horizon, length(order_list))
    
    bl = fit_beta_matrix(order_list[1:h], train)
    p = predict_beta_matrix(bl, train)
    pred = c(pred, p)
    actual = c(actual, test[i:(i+h-1)])
    
    train = c(train, test[i])
    max_horizon = max_horizon - 1
  }
  
  e = actual - pred
  
  return (e)
}

