################################################################################
# (c) Copyright 2015 Farzad Noorian. 
#
# Licensed under the Apache License, Version 2.0 (the "License"); 
# you may not use this file except in compliance with the License. 
# You may obtain a copy of the License at 
# 
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software 
# distributed under the License is distributed on an "AS IS" BASIS, 
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
# See the License for the specific language governing permissions and 
# limitations under the License.
################################################################################
# Functions for saving a data.frame as a latex table

save_tex_table <- function(table.data, filename, save.dir = ".", 
                           digits = 4,
                           caption = NULL,
                           caption.placement = "top",
                           table.placement = "t",
                           include.rownames = FALSE,
                           ...) {
  dir.create(save.dir, showWarnings = FALSE, recursive = TRUE)
  
  sink(file.path(save.dir, filename))

  if (require("xtable")) {
      latex_table = xtable(table.data, digits = digits, caption = caption, ...)

  print(latex_table,
        sanitize.rownames.function = identity, # function(x) x,
        sanitize.colnames.function = identity, # function(x) x,
        sanitize.text.function = identity, # function(x) x,
        digits = digits,
        caption.placement = caption.placement,
        table.placement = table.placement,
        include.rownames = include.rownames,
        hline.after=NULL, 
        add.to.row=list(pos=list(-1,0, nrow(table.data)),
                        command=c('\\toprule\n',
                                  '\\midrule\n',
                                  '\\bottomrule\n')),
        ...)
  } else {
    print(table.data, digits = digits)
  }

  sink()
}

bold_min_table_rows <- function(x, digits = 3) {
  # Wraps a "\\textbf{}" around the minimum element of each row
  # useful to point out the minimum element after optimisations for LATEX
  # INPUTS:
  #   x: a table
  #   digits: number of digits to round
  # RETURNS:
  #   a data.frame of strings, with minimum of each row wrapped in latex bold

  x = as.data.frame(round(x, digits))
  
  # bold the minimum of each row
  for (i in 1:nrow(x)) {
    ind.min = which.min(x[i,])
    x[i,ind.min] = paste0("\\textbf{", x[i,ind.min], "}")
  }
  
  return (x)
}

bold_min_table_cols <- function(x, digits = 3) {
  # Wraps a "\\textbf{}" around the minimum element of each column
  # useful to point out the minimum element after optimisations for LATEX
  # INPUTS:
  #   x: a table
  #   digits: number of digits to round
  # RETURNS:
  #   a data.frame of strings, with minimum of each column wrapped in latex bold
  
  x = as.data.frame(round(x, digits))
  
  # bold the minimum of each column
  for (i in 1:ncol(x)) {
    ind.min = which.min(x[,i])
    x[ind.min, i] = paste0("\\textbf{", x[ind.min, i], "}")
  }
  
  return (x)
}

