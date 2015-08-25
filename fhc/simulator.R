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
# receding horizon simulator

RHC_simulator <- function(V, E, horizons) {
  # the actual RHC simulator
  E_pointer = 1
  
  U_executed = NULL
  X_executed = NULL
  
  n_B = ncol(B_coeff) # used for selecting the current action
  n_C = ncol(C_coeff); stopifnot(n_C == 1) # only 1 input is supported in the emulator
  x = x0
  
  for (i in seq_along(horizons)) {
    N = horizons[i]
    current_E = E[E_pointer:(E_pointer+N-1)]
    current_V = V[i:(i+N-1)]
    E_pointer = E_pointer + N
    
    # old style Y
    #projected_Y = pinv(SC(N)) %*% SA(N) %*% x + current_V + current_E
    #U = JAB(N) %*% projected_Y
    
    # avoiding the pseudo inverse
    C_projected_Y = SA(N) %*% x + SC(N) %*% (current_V + current_E)
    U =  JAB_c(N) %*% C_projected_Y
    
    U_current = U[1:n_B]
    x = A_coeff %*% x +  B_coeff %*% U_current + C_coeff %*% current_V[1]
    
    U_executed = c(U_executed, U_current)
    X_executed = c(X_executed, x)
  }
  
  return (list(U = U_executed, X = X_executed))
}

Theta_simulator <- function(V, E, horizons) {
  # a simulator using theta and omega matrices  
  N = length(horizons)
  
  psi_y = Psi(horizons)
  phi_e = Phi(horizons)
  U_sim = psi_y %*% V + phi_e %*% E
  
  U_sim
}
