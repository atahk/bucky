# NEWS for Bucky's Archive for Data Analysis in the Social Sciences

## Release 1.0.6 (2019-12-17)
   * Fixed small parentheses bug in mi.eval
   * Updated robustify to mirror Stata defaults with linear models
   * Updated mi.eval docs regarding texreg & stargazer
   * Fixed bug with cluster argument in mi.eval
   * Fixed bad link to mice docs

## Release 1.0.5 (2018-10-29)
   * Fixed bugs with predict method for robustified objects
   * Added vcovHC method for polr objects (from MASS)
   * Added predict method for mi.estimates objects based on lm or glm

## Release 1.0.4 (2017-11-30)
   * Moved methods, sandwich, lmtest from "Depends" to "Imports"
   * Fixed robustify examples ("coeftest" replaced by "summary")
   * Added option for lazy evaluation to "mi.eval"
   * Changed default for robust standard errors to HC0 to match Stata
   * Added predict method for robustified objects based on lm or glm
     models
   * Improved compatibility between "mi.eval" and other packages
     * Added compatibility with 'mice'
     * Fixed compatibility with 'mitools'
     * Added compatibility with 'stargazer' when using lm, glm, and a
       few other types of models

## Release 1.0.3 (2017-11-20)

   * Initial submission to CRAN
   * Initial check-in to Github
