##' bla
##'
##' blub
##' @title atitle
##' @param train.fit model fitted to training data. 
##' @param dt dataset. 
##' @param risk.set risk.set used for partial likelihood. 
##' @param test.set validation data. 
##' @param X design matrix if one of c("coxnet", "cv.glmnet", "glmnet") is used.
##' @param lambda.cv grid over which to choose penalization if one of c("coxnet", "cv.glmnet", "glmnet") is used.
##' @param delta.var name of event type variable. 
##' @param delta.value type of event of interest here-
##' @param change.point specified if there is a changepoint in the effect of treatment across time.
##' @return a
##' @seealso b
##' @examples c
##' @export 
##' @author Helene C. W. Rytgaard <hely@@biostat.ku.dk>
cox.loss.fun <- function(train.fit, dt, risk.set, test.set, X=NULL, lambda.cv=NULL,
                         delta.var="delta", delta.value=1, change.point=NULL) {

    tmp <- copy(dt)
   
    if (ifelse(length(change.point)>0, change.point>0, FALSE)) {
        tmp <- rbind(tmp, tmp)[order(id)]
        tmp[, time.period:=1:.N, by="id"]
    } 
    
    if (any(class(train.fit)%in%c("coxnet", "cv.glmnet", "glmnet"))) {
        tmp[, fit.lp:=predict(train.fit, type="link", newx=X, s=lambda.cv)]
    } else {
        tmp[, fit.lp:=predict(train.fit, type="lp", newdata=tmp)]
    }

    tmp[, risk:=0]
    tmp[id%in%risk.set, risk:=1]

    tmp <- tmp[rev(order(time))]
   
    if (ifelse(length(change.point)>0, change.point>0, FALSE)) {
        tmp[, term2:=cumsum(risk*exp(fit.lp)), by="time.period"]
    } else {
        tmp[, term2:=cumsum(risk*exp(fit.lp))]
    }
    
    tmp[term2==0, term2:=1]

    if (ifelse(length(change.point)>0, change.point>0, FALSE)) {
        return(-sum(tmp[id%in%test.set &
                        ((time.period==1 & time<=change.point) | (time.period==2 & time>change.point)),
        (get(delta.var)==delta.value)*(fit.lp - log(term2))]))
    } else {
        return(-sum(tmp[id%in%test.set, (get(delta.var)==delta.value)*(fit.lp - log(term2))]))
    }
    
}
