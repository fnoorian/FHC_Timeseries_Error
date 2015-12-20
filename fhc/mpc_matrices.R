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

MPC_S <- function(n, A, B) {
  
  nr = nrow(B)
  nc = ncol(B)
  
  S = matrix(0, nrow=n*nr, ncol=n*nc)
  
  for (i in 1:n) {
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

SA <- function(n) {
  ret = NULL
  for (i in 1:n) {
    ret = rbind(ret, matrix_power(A_coeff, i))
  }
  
  return (ret)
}

SB <- function(n) MPC_S(n, A_coeff, B_coeff)
SC <- function(n) MPC_S(n, A_coeff, C_coeff)
SI <- function(n) MPC_S(n, A_coeff, diag(1,nrow(A_coeff),ncol(A_coeff)))

#############
# The cost function of LQ control
# J = U * JA * U + 2 * U * JB * f + f * JC * f
JA <- function(n) P(n) + t(SB(n)) %*% Q(n) %*% SB(n)
JB <- function(n) t(SB(n)) %*% t(Q(n)) %*% SC(n)
JC <- function(n) t(SC(n)) %*% Q(n)  %*% SC(n)
JAB <- function(n) -solve(JA(n)) %*% JB(n)
JAB_c <- function(n) -solve(P(n) + t(SB(n)) %*% Q(n) %*% SB(n), t(SB(n)) %*% t(Q(n))) # this is JAB = JAB_c %*% SC(n)

#############
# Closed form cost function of LQ control
cost_UV <- function(U, V, n) t(U) %*% JA(n) %*% U + 2 * t(U) %*% JB(n) %*% V + t(V) %*% JC(n) %*% V
cost_XUV <- function(X, U, V, n) t(X) %*% Q(n) %*% X + t(U) %*% P(n) %*% U
