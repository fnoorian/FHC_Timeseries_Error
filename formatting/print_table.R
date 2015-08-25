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

  print(xtable(table.data, digits = digits, caption = caption, ...),
        sanitize.rownames.function = function(x) x,
        sanitize.colnames.function = function(x) x,
        sanitize.text.function = function(x) x,
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
    print(table.data)
  }

  sink()
}
