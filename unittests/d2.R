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
# Unit test 2 for the LQ Err

rm(list=ls())

source("unittests/utils.R")
source("utils/matrix_tools.R")
source("fhc/mpc_matrices.R")
source("fhc/theta.R")
source("fhc/simulator.R")
source("prediction/linear_predictor.R")

test_simulator <- function(V, E, horizons) {
  
  N = length(horizons)
  
  r1 = RHC_simulator(V, E, horizons)
  r2 = RHC_simulator(V, E * 0, horizons)
  
  U_x = Theta_simulator(V, E, horizons)
  U_y = Theta_simulator(V, 0*E, horizons)
  
  assert_vector_equality(U_x, r1$U)
  assert_vector_equality(U_y, r2$U)
  
  cost_U <- function(U) cost_UV(U, V[1:N], N)
  cost_XU <- function(X, U) cost_XUV(X, U, V[1:N], N)
  
  assert_vector_equality(cost_U(r1$U), cost_XU(r1$X, r1$U))
  assert_vector_equality(cost_U(r2$U), cost_XU(r2$X, r2$U))

  theta = Theta_matrix(horizons)
  omega = Omega_matrix(horizons)
  Omega_V = omega %*% V
  
  delta_j1 = cost_U(r1$U) - cost_U(r2$U)
  delta_j2 = t(E) %*% theta %*% E + t(E) %*% Omega_V
  
  assert_vector_equality(delta_j1, delta_j2)
}

h = 3

x0 = rep(0, h + 1)
A_coeff = upper.diagonal(h)
A_coeff = rbind(cbind(A_coeff, 0), 0)
A_coeff[h+1,1] = 1

B_coeff = diag(h)
B_coeff = rbind(B_coeff, 0)

C_coeff = matrix(0, h+1, 1)
C_coeff[h+1, 1] = -1

Q <- function(N) diag(rep(c(0, 0, 0, 8), N)) # state weights
P <- function(N) diag(rep(c(1, 2, 4), N)) # inputs weights

horizons = c(4, 3, 3, 3)
N_run = length(horizons) + tail(horizons,1) - 1
N_pred = sum(horizons)
test_simulator(runif(N_run), runif(N_pred), horizons)

