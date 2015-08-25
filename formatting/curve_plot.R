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
# Functions for creating PDF fitted for Publications

pdf_create <- function(filename, width = 16, height = 12, no.individ.title=TRUE,
                       image.dir = ".") {

  dir.create(image.dir, showWarnings = FALSE, recursive = TRUE)
  
  filename = file.path(image.dir, filename)                       
  pdf(filename, width=width/2.54, height=height/2.54)
  
  # no margin
  par(oma=c(0,0,0,0))
  
  if (no.individ.title) {
    par(mar=c(2.3, 2.5, 0.7  , 0.5))
    par(mgp=c(1.3,0.4,0))
  }
  else {
    par(mar=c(2.3, 2.5, 1.5	, 0.5))
    par(mgp=c(1.3,0.4,0))
  }  
}

png_create <- function(filename, width = 16, height = 12, no.individ.title=TRUE,
                       image.dir = ".") {

  dir.create(image.dir, showWarnings = FALSE, recursive = TRUE)
  
  filename = file.path(image.dir, filename)                       
  png(filename, width=width/2.54, height=height/2.54, units="in", res=600)
  
  # no margin
  par(oma=c(0,0,0,0))
  
  if (no.individ.title) {
    par(mar=c(2.3, 2.5, 0.7  , 0.5))
    par(mgp=c(1.3,0.4,0))
  }
  else {
    par(mar=c(2.3, 2.5, 1.5  , 0.5))
    par(mgp=c(1.3,0.4,0))
  }  
}


