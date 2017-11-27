# Bucky's Archive for Data Analysis in the Social Sciences

This is an R package that provides functions for various statistical techniques commonly used in the social sciences, including functions to compute clustered robust standard errors, combine results across multiply-imputed data sets, and simplify the addition of robust and clustered robust standard errors. The package was originally developed, in part, to assist porting of replication code from Stata and attempts to replicate default options from Stata where possible.

## Installation

Most users should use latest stable release of the packge, which can be installed from [CRAN](https://cran.r-project.org/) by running
```R
install.packages("bucky")
```

The development version can be installed directly from [GitHub](https://github.com/atahk/pscl) by running
```R
## load 'devtools', installing if not already installed
suppressWarnings(require("devtools")) || {install.packages("devtools"); library("devtools")}
## install development version of 'bucky'
install_github("atahk/bucky")
```
