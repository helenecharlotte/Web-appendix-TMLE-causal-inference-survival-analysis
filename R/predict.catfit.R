### predict.catfit.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 29 2021 (14:33) 
## Version: 
## Last-Updated: Jul 29 2021 (17:00) 
##           By: Thomas Alexander Gerds
##     Update #: 22
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' @export predict.catfit
#' @export
predict.catfit <- function(object, newdata,...){
    pseudo.newdata <- get.pseudodata(object$catfit.call$formula,data=newdata,atrisk=FALSE)
    pseudo.newdata <- data.table::data.table(pseudo.newdata)
    outcome.values <- object$y
    outcome <- attr(pseudo.newdata,"outcome")
    # predict pseudo hazard
    data.table::set(pseudo.newdata,j="predict.Y",value=stats::predict.glm(object, newdata=pseudo.newdata, type="response"))
    # atrisk probability
    pseudo.newdata[,predict.atrisk:=c(1,cumprod(1-predict.Y[-length(predict.Y)])),by=id]
    # map to outcome probabilities
    pseudo.newdata[,{
        # for the first value of Y we have set atrisk=1 
        prob:=predict.Y*predict.atrisk
    }]
    # flatten
    predicted.values <- dcast(id~pseudo.aa,data=pseudo.newdata,value.var="prob",fill=0)
    return(cbind(newdata, predicted.values))
}


#----------------------------------------------------------------------
### predict.catfit.R ends here
