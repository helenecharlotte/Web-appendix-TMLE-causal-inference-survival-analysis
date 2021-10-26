### load.R --- 
#----------------------------------------------------------------------

#-------------------------------------------------------------------------------------------#
## packages and generic functions
#-------------------------------------------------------------------------------------------#

library(data.table)
library(zoo)
library(glmnet)
library(stringr)
library(nleqslv)
library(prodlim)
library(ggplot2)
library(gridExtra)
library(survival)
library(riskRegression)
library(MASS)

#-------------------------------------------------------------------------------------------#
## source code
#-------------------------------------------------------------------------------------------#

source("./R/sim.data.continuous.R") 
source("./R/contmle.R")
source("./R/cox.loss.fun.R") 
source("./R/lebesgue.loss.fun.R")
source("./R/cv.fun.R")     
source("./R/basis.fun.R")
source("./R/hal.screening.R")
source("./R/fit.hal.R")   
source("./R/cox.sl.R")  
source("./R/fit.categorical.R")
source("./R/predict.catfit.R")

######################################################################
### load.R ends here




