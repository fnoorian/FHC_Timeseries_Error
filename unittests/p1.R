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
# Unit test 3 for the LQ Err (using step by step proof)

rm(list=ls())

source("unittests/utils.R")
source("utils/matrix_tools.R")
source("fhc/mpc_matrices.R")
source("fhc/theta.R")
source("fhc/simulator.R")

#################################
# System definitions

x0 = 0
A_coeff = matrix(2, 1, 1)
B_coeff = matrix(3, 1, 1)
C_coeff = matrix(4, 1, 1)

P <- function(n) diag(n) * 2.5
Q <- function(n) diag(n) * 3.5

horizons = c(5)#c(rep.int(5,5), 5:1)

################################################
# The problem setup
# TODO: remove 1, see what happens and why

horizons = c(5, 5, 4, 3, 2, 1) # must be 6 horizons (hard coded)
N_steps = length(horizons)
N_data = length(horizons) + tail(horizons,1) - 1
N_pred = sum(horizons)

V = runif(N_data)
E = runif(N_pred)

Y = V + pinv(SC(N_data)) %*% SA(N_data) %*% x0 
################################################
VL = lapply(1:N_steps, function(i) V[1:horizons[i] + (i-1)])
EL = lapply(1:N_steps, function(i) E[sum(horizons[0:(i-1)]) + 1:horizons[i]])
################################################
X0L = list()
YL = list()
CYL = list()
UL = list()

U_executed = NULL
X_executed = NULL

n_A = ncol(A_coeff)
n_B = ncol(B_coeff)
n_C = ncol(C_coeff); stopifnot(n_C == 1) # only 1 input is supported in the emulator
x = x0

for (i in 1:N_steps) {
  n = horizons[i]
  current_E = EL[[i]]
  current_V = VL[[i]]
  
  projected_Y = pinv(SC(n)) %*% SA(n) %*% x + current_V + current_E
  C_projected_Y = SA(n) %*% x + SC(n) %*% (current_V + current_E)
  
  X0L[[i]] = x
  YL[[i]] = projected_Y
  CYL[[i]] = C_projected_Y
  
  U_x = JAB(n) %*% projected_Y
  U =  JAB_c(n) %*% C_projected_Y
  UL[[i]] = U
  
  assert_vector_equality(U, U_x)
  
  U_current = U[1:n_B]
  x = A_coeff %*% x +  B_coeff %*% U_current + C_coeff %*% current_V[1]
  
  U_executed = c(U_executed, U_current)
  X_executed = c(X_executed, x)
}
#####################################
Y = V

SC_Y1 = SC_M_Y(N_data, horizons[1], t=1) %*% Y
SC_E1 = SC(horizons[1]) %*% EL[[1]]

stopifnot(all(SC(horizons[1]) %*% VL[[1]] == SC_Y1))
stopifnot(all(SC(horizons[1]) %*% EL[[1]] == SC_E1))

SC_UPS1 = SC_Y1 + SC_E1

U1 = JAB_c(horizons[1]) %*% SC_UPS1
assert_vector_equality(UL[[1]], U1)

SB_U1 = SA(horizons[2]) %*% B_coeff %*% M_U(horizons[1]) %*% U1

SC_Y2 = SC_M_Y(N_data, horizons[2], t=2) %*% Y
SC_E2 = SC(horizons[2]) %*% EL[[2]]

SC_UPS2 = SC_Y2 + SC_E2 + SB_U1
U2 = JAB_c(horizons[2]) %*% SC_UPS2
assert_vector_equality(UL[[2]], U2)

SB_U12 = SA(horizons[3]) %*% A_coeff %*% B_coeff %*% M_U(horizons[1]) %*% U1
SB_U2 = SA(horizons[3]) %*% B_coeff %*% M_U(horizons[2]) %*% U2

SC_Y3 = SC_M_Y(N_data, horizons[3], t=3) %*% Y
SC_E3 = SC(horizons[3]) %*% EL[[3]]

SC_UPS3 = SC_Y3 + SC_E3 + SB_U2 + SB_U12
U3 = JAB_c(horizons[3]) %*% SC_UPS3
assert_vector_equality(UL[[3]], U3)

SB_U13 = SA(horizons[4]) %*% A_coeff %*% A_coeff %*% B_coeff %*% M_U(horizons[1]) %*% U1
SB_U23 = SA(horizons[4]) %*% A_coeff %*% B_coeff %*% M_U(horizons[2]) %*% U2
SB_U3 = SA(horizons[4]) %*% B_coeff %*% M_U(horizons[3]) %*% U3

SC_Y4 = SC_M_Y(N_data, horizons[4], t=4) %*% Y
SC_E4 = SC(horizons[4]) %*% EL[[4]]

SC_UPS4 = SC_Y4 + SC_E4 + SB_U3 + SB_U23 + SB_U13
U4 = JAB_c(horizons[4]) %*% SC_UPS4
assert_vector_equality(UL[[4]], U4)
#####################################
Y = V

SC_Y1 = SC_M_Y(N_data, horizons[1], t=1) %*% Y
SC_E1 = SC(horizons[1]) %*% EL[[1]]
SC_YE_1 = SC_Y1 + SC_E1

SC_UPS1 = Gamma(horizons, 1, 1) %*% SC_YE_1

U1 = JAB_c(horizons[1]) %*% SC_UPS1
assert_vector_equality(UL[[1]], U1)

SB_U1 =  B_coeff %*% M_U(horizons[1]) %*% U1

SC_Y2 = SC_M_Y(N_data, horizons[2], t=2) %*% Y
SC_E2 = SC(horizons[2]) %*% EL[[2]]
SC_YE_2 = SC_Y2 + SC_E2

SC_UPS2 = SC_Y2 + SC_E2 + SA(horizons[2]) %*% SB_U1
SC_UPS2 = SC_Y2 + SC_E2 + Gamma(horizons, 1, 2) %*% SC_UPS1
SC_UPS2 = Gamma(horizons, 2, 2) %*% SC_YE_2 + Gamma(horizons, 1, 2)  %*% SC_YE_1
U2 = JAB_c(horizons[2]) %*% SC_UPS2
assert_vector_equality(UL[[2]], U2)

SB_U2 = B_coeff %*% M_U(horizons[2]) %*% U2

SC_Y3 = SC_M_Y(N_data, horizons[3], t=3) %*% Y
SC_E3 = SC(horizons[3]) %*% EL[[3]]
SC_YE_3 = SC_Y3 + SC_E3

SC_UPS3 = SC_Y3 + SC_E3 + SA(horizons[3]) %*% (A_coeff %*% SB_U1 + SB_U2)
SC_UPS3 = Gamma(horizons, 3,3) %*% SC_YE_3 + 
  (Gamma_internal(2, 3) %*% Gamma_internal(1, 2) %*% Gamma_internal(1, 1) + Gamma_internal(1, 3) %*% Gamma_internal(1, 1)) %*% SC_YE_1 + 
  Gamma_internal(2, 3) %*% Gamma_internal(2, 2) %*% SC_YE_2 
SC_UPS3 = Gamma(horizons, 3,3) %*% SC_YE_3 + Gamma(horizons, 1, 3) %*% SC_YE_1 + Gamma(horizons, 2, 3) %*% SC_YE_2
U3 = JAB_c(horizons[3]) %*% SC_UPS3
assert_vector_equality(UL[[3]], U3)

SB_U3 =  B_coeff %*% M_U(horizons[3]) %*% U3

SC_Y4 = SC_M_Y(N_data, horizons[4], t=4) %*% Y
SC_E4 = SC(horizons[4]) %*% EL[[4]]
SC_YE_4 = SC_Y4 + SC_E4

SC_UPS4 = SC_Y4 + SC_E4 + SA(horizons[4]) %*% (SB_U3 + A_coeff %*% SB_U2 + matrix_power(A_coeff, 2) %*% SB_U1)
SC_UPS4 = SC_Y4 + SC_E4 + Gamma_internal(1, 4) %*% SC_UPS1 + Gamma_internal(2, 4) %*% SC_UPS2 + Gamma_internal(3, 4) %*% SC_UPS3
SC_UPS4 = SC_Y4 + SC_E4 + Gamma(horizons, 1, 4) %*% SC_YE_1 + Gamma(horizons, 2, 4) %*% SC_YE_2 + Gamma(horizons, 3, 4) %*% SC_YE_3
U4 = JAB_c(horizons[4]) %*% SC_UPS4
assert_vector_equality(UL[[4]], U4)

Phi1 = M_phi(N_steps, horizons[1], 1) %*% JAB_c(horizons[1]) %*% Gamma(horizons, 1,1) + 
  M_phi(N_steps, horizons[2], 2) %*% JAB_c(horizons[2]) %*% Gamma(horizons, 1,2)+ 
  M_phi(N_steps, horizons[3], 3) %*% JAB_c(horizons[3]) %*% Gamma(horizons, 1,3)+ 
  M_phi(N_steps, horizons[4], 4) %*% JAB_c(horizons[4]) %*% Gamma(horizons, 1,4)
Phi2 = M_phi(N_steps, horizons[2], 2) %*% JAB_c(horizons[2]) %*% Gamma(horizons, 2,2)+ 
  M_phi(N_steps, horizons[3], 3) %*% JAB_c(horizons[3]) %*% Gamma(horizons, 2,3)+ 
  M_phi(N_steps, horizons[4], 4) %*% JAB_c(horizons[4]) %*% Gamma(horizons, 2,4)
Phi3 = M_phi(N_steps, horizons[3], 3) %*% JAB_c(horizons[3]) %*% Gamma(horizons, 3,3)+ 
  M_phi(N_steps, horizons[4], 4) %*% JAB_c(horizons[4]) %*% Gamma(horizons, 3,4)
Phi4 = M_phi(N_steps, horizons[4], 4) %*% JAB_c(horizons[4]) %*% Gamma(horizons, 4,4)

Ux = Phi1 %*% SC_YE_1 + Phi2 %*% SC_YE_2 + Phi3 %*% SC_YE_3 + Phi4 %*% SC_YE_4
n_h4 = 4 * n_B
assert_vector_equality(U_executed[1:n_h4],Ux[1:n_h4])

Phi_list_x = Phi_list(horizons)

SC_Y5 = SC_M_Y(N_data, horizons[5], t=5) %*% Y
SC_E5 = SC(horizons[5]) %*% EL[[5]]
SC_YE_5 = SC_Y5 + SC_E5

SC_Y6 = SC_M_Y(N_data, horizons[6], t=6) %*% Y
SC_E6 = SC(horizons[6]) %*% EL[[6]]
SC_YE_6 = SC_Y6 + SC_E6

Ux2 = Phi_list_x[[1]] %*% SC_YE_1 + 
  Phi_list_x[[2]] %*% SC_YE_2 +
  Phi_list_x[[3]] %*% SC_YE_3 +
  Phi_list_x[[4]] %*% SC_YE_4 +
  Phi_list_x[[5]] %*% SC_YE_5 +
  Phi_list_x[[6]] %*% SC_YE_6

assert_vector_equality(Ux2, U_executed)

Ux_y0 = Phi_list_x[[1]] %*% SC_Y1 + 
  Phi_list_x[[2]] %*% SC_Y2 +
  Phi_list_x[[3]] %*% SC_Y3 +
  Phi_list_x[[4]] %*% SC_Y4 +
  Phi_list_x[[5]] %*% SC_Y5 +
  Phi_list_x[[6]] %*% SC_Y6
Ux_y1 = Psi(horizons) %*% Y
assert_vector_equality(Ux_y0, Ux_y1)

Ux3 = Psi(horizons) %*% Y + Phi(horizons) %*% E
assert_vector_equality(Ux2, Ux3)

rhc1 = RHC_simulator(Y, E, horizons)
assert_vector_equality(Ux3, rhc1$U)
