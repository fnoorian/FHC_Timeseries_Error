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
# Unit test utilities

assert_vector_equality <- function(a, b, THRESHOLD = 1e-8) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  stopifnot(length(a) == length(b))
  ratio = abs((a - b)/(abs(a) + abs(b)))
  ratio[is.nan(ratio)] = 0
  stopifnot(max(ratio) < THRESHOLD)
}

