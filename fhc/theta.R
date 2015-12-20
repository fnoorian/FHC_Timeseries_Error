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
# Definition of Matrices Theta and Omega

# M_U function, as described in the paper
# M_U_c <- function(N, n, t) {
#   # N is size of all Y
#   # n is the horizon
#   # t is the current step
#   ret = matrix(0, N * n_C, n * n_B)
#   
#   replacement = pinv(C_coeff) %*% matrix_power(A_coeff, t) %*% B_coeff
#   ret = matrix_block_replace(ret, 1, 1, replacement)
#   ret
# }

# M_U function, simplified
M_U <- function(n) {
  # n is the current horizon len
  n_B = ncol(B_coeff)
  
  ret = matrix(0, n_B, n * n_B)
  
  replacement = diag(n_B)
  ret = matrix_block_replace(ret, 1, 1, replacement)
  ret
}

M_phi <- function(N, n, t) {
  n_B = ncol(B_coeff)

  I = diag(n_B)
  
  ret = matrix(0, N * n_B, n * n_B) 
  ret = matrix_block_replace(ret, t, 1, I)
  ret
}

SC_M_Y <- function(N, n, t) {
  # returns the SC %*% M_Y = SI %*% (C * MY)
  # N : size of all Y
  # n : the current horizon length
  # t : the current time step
  n_A = ncol(A_coeff)
  n_C = ncol(C_coeff)
  
  ret = matrix(0, n * n_A, N * n_C)
  
  for (i in 1:n) {
    ret = matrix_block_replace(ret, i, (t-1)+i, C_coeff)
  }
  if (t > 1) {
    for (i in 1:(t-1)) {
      replacement = matrix_power(A_coeff, t - i) %*% C_coeff
      ret = matrix_block_replace(ret, 1, i, replacement)
    }
  }
  
  SI(n) %*% ret
}

SC_M_E <- function(N, n, start) {
  # returns the SC %*% M_E 
  # M_E is similar to M_Y but doesn't accumulate error
  # N : size of all E
  # n : the current horizon length
  # start : start of the current time step in the vector
  n_A = ncol(A_coeff)
  n_C = ncol(C_coeff)
  
  ret = matrix(0, n * n_A, N * n_C)
  
  for (i in 1:n) {
    ret = matrix_block_replace(ret,  i, (start - 1) + i, C_coeff)
  }
  
  
  SI(n) %*% ret
}

Gamma_internal <- function(instance, t) {
  
  N_data = length(horizons) + tail(horizons,1) - 1
  N = horizons[t]
  M = horizons[instance]
  
  if (instance > t) {
    ret = matrix(0, nrow= n_A, ncol=N_data)
  } else if (instance == t) {
    return (diag(n_A * N))
  } else {
    ret = matrix_power(A_coeff, t - instance - 1) %*% B_coeff %*% M_U(M) %*% JAB_c(M)
  }
  
  SA(N) %*% ret 
}

Gamma_nomemoise <- function(horizons, instance, t) {
  n_A = ncol(A_coeff)
  
  N_data = length(horizons) + tail(horizons,1) - 1
  N = horizons[t]
  M = horizons[instance]
  
  if (instance > t) {
    return (matrix(0, nrow= n_A * N, ncol=n_A * M))
  } else if (instance == t) {
    return (diag(n_A * N))
  } else {
    ret = 0
    for (inst in 1:(t-1)) {
      M = horizons[inst]
      g_inst = matrix_power(A_coeff, t - inst - 1) %*% B_coeff %*% M_U(M) %*% JAB_c(M)
      ret = ret + g_inst %*% Gamma(horizons, instance, inst)
    }
  }
  return (SA(N) %*% ret)
}

# Memoise gamma if available
if (require("memoise")) {
  Gamma = memoise(Gamma_nomemoise)
} else {
  Gamma = Gamma_nomemoise
}

Phi_list <- function(horizons) {
  #horizons: horizon for each step
  
  # check size of variable per
  n_A = ncol(A_coeff)
  n_B = ncol(B_coeff)
  n_C = ncol(C_coeff)
  stopifnot(nrow(A_coeff) == n_A)
  stopifnot(nrow(B_coeff) == n_A)
  stopifnot(nrow(C_coeff) == n_A)
  
  # compute Phi
  N_steps = length(horizons) # total steps
  Phi = list()
  for (i in 1:N_steps) {
    Phi[[i]] = 0
    for (j in i:N_steps) {
      Phi[[i]] = Phi[[i]] + M_phi(N_steps, horizons[j], j) %*% JAB_c(horizons[j]) %*% Gamma(horizons,i,j)
    }
  }
  
  return (Phi)
}


Phi_vector <- function(horizons) {
	Phi = Phi_list(horizons)

  # Convert Phi to Theta
  bd = do.call(cbind, Phi)
	
	return (bd)
}
  
Psi <- function(horizons) {
  # Convert Phi to Psi for Y
  phi = Phi_list(horizons)
  
  N = length(horizons) + tail(horizons,1) - 1 # N: total steps
  steps = length(horizons)
  
  psi_sum = 0
  for (i in 1:steps) {
    psi_sum = psi_sum + phi[[i]] %*% SC_M_Y(N, horizons[[i]], i)
  }
  
  return (psi_sum)
}

Phi <- function(horizons) {
  # Convert Phi to Phi for E
  phi = Phi_list(horizons)
  
  N = sum(horizons) # N: length of error vector
  steps = length(horizons)
  
  psi_sum = 0
  for (i in 1:steps) {
    prev_horizons = 0
    if (i > 1) {
      prev_horizons = sum(horizons[1:(i-1)])
    }

    # use SC_M_E to get JAB_C %*% SC %*% ME 
    psi_sum = psi_sum + phi[[i]] %*% SC_M_E(N, horizons[[i]], prev_horizons + 1)
  }
  
  return (psi_sum)
}

Theta_matrix <- function(horizons) {
  phi = Phi(horizons)
  N = length(horizons) # N: total steps
  
  return (t(phi) %*% JA(N) %*% phi)
}

compress_Theta <- function(Theta, L) {
  a = eigen(Theta)
  sigma = (a$values)
  v = (a$vectors)
  
  w = v[,1:L] %*% diag(sqrt(sigma[1:L]))
  w = Re(w) # for times when eigenvalue decomposition returns complex numbers
  
  return (w)
}

Omega_matrix <- function(horizons) {
  N = length(horizons) # N: total steps
  N_data = length(horizons) + tail(horizons,1) - 1
  
  Phi = Phi(horizons)
  psi = Psi(horizons)  
  
  JB_executed = JB(N)%*% diag(1,N,N_data)
  
  return (2 * t(Phi) %*% (JA(N) %*% psi + JB_executed))
}
