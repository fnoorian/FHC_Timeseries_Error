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
# Unit test for linear predictor

rm(list=ls())

source("prediction/linear_predictor.R")
source("unittests/utils.R")

test_linear_predictor <- function() {
  # a simple test to validate predictor working
  train_ts = 1:20 # this is very predictable!
  test_ts = 21:25
  length_pred = length(test_ts)
  err = predict_horizons_error(rep.int(2, length_pred), train_ts, test_ts)
  stopifnot(length(err) == (length_pred * (length_pred+1)/2))   # there should be 5-4-3-2-1 predictions
  stopifnot(sum(abs(err)) < 1e-10)  
}

test_linear_predictor()
