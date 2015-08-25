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
# usefull matrix function

pinv <- function(x) {
  #pseudo inverse
  solve(t(x) %*% x, t(x))
}

matrix_power <- function(A, N) {
  # A ^ N = A %*% A %*% ... %*% A
  ret = A
  if (N == 0) {
    ret = diag(1, nrow=nrow(A), ncol=ncol(A))
  } else if (N >= 2) {
    for (i in 2:N) {
      ret = ret %*% A
    }
  }
  
  ret
}

matrix_block_replace = function(mat, i, j, block) {
  # replaces a matrix as a block within another matrix
  # i & j are coordinates of the block (i * nrow(block), j * ncol(block))
  r = nrow(block)
  l = ncol(block)
  
  mat[(i-1)*r + 1:r, (j-1)*l + 1:l] = block
  
  mat
}

upper.diagonal = function(N) {
  # generate an upper diagonal matrix
  # i.e., only the single element above the diagonal is 1
  x = matrix(0, N, N)
  
  for (i in 1:(N-1)) {
    x[i, i+1] = 1
  }
  x
}

lower.diagonal <- function(N) {
  # generate an upper diagonal matrix
  # i.e., only the single element above the diagonal is 1
  x = matrix(0, N, N)
  for (i in 1:(N-1)) {
    x[i+1, i] = 1
  }
  x
}

