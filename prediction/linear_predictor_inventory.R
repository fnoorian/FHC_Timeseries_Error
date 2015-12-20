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
# linear prediction functions for the inventory problem

source("linear_predictor.R")

predict_inventory_error <- function(orders, train, test) {
  # computes error of predicting "train" against "test" using AR models for the inventory problem
  # train: observed values
  # test: target of prediction
  # orders: p for ar(p)
  # returns: error for each horizon, concatenated
  
  max.horizon = length(orders)
  
  pred = list()
  actual = list()
  
  for (i in 1:length(test)) {
    new_train = c(train, test[0:i])
    new_traing_mean = mean(new_train)
    new_train = new_train - new_traing_mean
    
    # fit the beta matrix
    beta = fit_beta_matrix(order = orders, ts = new_train)
    
    # get the horizon
    horizon = min(length(test) - i , max.horizon)
    
    # forecast 
    prediction = c()
    if (horizon > 0) {
      p = predict_beta_matrix(beta, new_train)
      prediction = p[1:horizon] + new_traing_mean
    }
    
    pred[[i]] = c(test[i], prediction)
    actual[[i]] = test[i + 0:(horizon)]
  }
  
  e = do.call(c, pred) - do.call(c, actual)
  
  return (e)
}

predict_inventory_naive_error <- function(train, test, max.horizon) {
  # computes error of predicting "train" against "test" using Naive technique for the inventory problem
  # train: observed values
  # test: target of prediction
  # max.horizon: max length of horizon
  # returns: error for each horizon, concatenated
  
  pred = list()
  actual = list()
  
  for (i in 1:length(test)) {
    new_train = c(train, test[0:i])

    # get the horizon
    horizon = min(length(test) - i , max.horizon)
    
    # forecast (i.e., repeat the last identified value)
    prediction = c()
    if (horizon > 0) {
      prediction = rep.int(tail(new_train, 1), horizon)
    }
    
    pred[[i]] = c(test[i], prediction)
    actual[[i]] = test[i + 0:(horizon)]
  }
  
  e = do.call(c, pred) - do.call(c, actual)
  
  return (e)
}
