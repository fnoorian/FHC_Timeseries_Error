On Time Series Forecasting Error Measures for Finite Horizon Control
====================================================================

The code contains implementation of the example from the paper: 
"On Time Series Forecasting Error Measures for Finite Horizon Control", 
by F. Noorian and P. H. W. Leong, 
submitted to IEEE Transactions on Control Systems Technology for review.

### Running the example
If not installed, obtain and install R from http://www.r-project.org/ .

In the command prompt, type:

 > source("main.R")

The code will generate the latex tables and PDF images used in the paper.
It will create an "Output" directory in the same folder to store its results.

Although not required, we suggest installing memoise and xtable packages.

 > install.packages(c("xtable", "memoise"))

### Files
- main.R: Runs all tests.
- main_inventory.R: The inventory problem simulation.
- timing_inventory.R: timing tests for the inventory problem.
- main_stocks.R: the stocks problem simulation.
- timing_stocks.R: timing tests for the stocks problem.
- fhc/mpc_matrices.R: Basic matrix definitions for Finite horizon control
- fhc/simulator.R: Step-by-step FHC simulator, to be used in comparison.
- fhc/theta.R: Generation of matrices used in formulas based on above parameters.
- prediction/linear_predictor.R: Function for linear AR model fitting and evaluation.
- prediction/linear_predictor_power.R: Function for linear AR model fitting and evaluation with power transform.
- utils/matrix_tools.R: Useful matrix manipulation functions
- formatting/*: Functions for formatting plots and tables
- unittests/*: Unit-tests


### Contact Information
 * Farzad Noorian <farzad.noorian@sydney.edu.au>
 * Philip H. W. Leong <philip.leong@sydney.edu.au>

### License
Copyright (c) 2015 Farzad Noorian <farzad.noorian@sydney.edu.au>.

All files in this package licensed under the Apache License, Version 2.0 (the "License");
you may not use these files except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
