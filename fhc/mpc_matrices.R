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
# Definition of Matrices S and J (MPC basic definitions)

MPC_S <- function(N, A, B) {
  
  nr = nrow(B)
  nc = ncol(B)
  
  S = matrix(0, nrow=N*nr, ncol=N*nc)
  
  for (i in 1:N) {
    posx = (i-1) * nr + (1:nr)
    posy = (i-1) * nc + (1:nc)
    S[posx, posy] = B
    
    if (i > 1) {
      for (j in 1:(i-1)) {
        posy = (j-1) * nc + (1:nc)
        S[posx, posy] = matrix_power(A, (i-j)) %*% B
      }
    }
  }

  S
}

SA <- function(N) {
  ret = NULL
  for (i in 1:N) {
    ret = rbind(ret, matrix_power(A_coeff, i))
  }
  
  return (ret)
}

SB <- function(N) MPC_S(N, A_coeff, B_coeff)
SC <- function(N) MPC_S(N, A_coeff, C_coeff)
SI <- function(N) MPC_S(N, A_coeff, diag(1,nrow(A_coeff),ncol(A_coeff)))

#############
# The cost function of LQ control
# J = U * JA * U + 2 * U * JB * f + f * JC * f
JA <- function(N) P(N) + t(SB(N)) %*% Q(N) %*% SB(N)
JB <- function(N) t(SB(N)) %*% t(Q(N)) %*% SC(N)
JC <- function(N) t(SC(N)) %*% Q(N)  %*% SC(N)
JAB <- function(N) -solve(JA(N)) %*% JB(N)
JAB_c <- function(N) -solve(P(N) + t(SB(N)) %*% Q(N) %*% SB(N), t(SB(N)) %*% t(Q(N))) # this is JAB = JAB_c %*% SC(N)

#############
# Closed form cost function of LQ control
cost_UV <- function(U, V, N) t(U) %*% JA(N) %*% U + 2 * t(U) %*% JB(N) %*% V + t(V) %*% JC(N) %*% V
cost_XUV <- function(X, U, V, N) t(X) %*% Q(N) %*% X + t(U) %*% P(N) %*% U
