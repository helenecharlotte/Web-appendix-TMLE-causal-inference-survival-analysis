##' bla
##'
##' blub
##' @title atitle
##' @param train.fit model fitted to training data. 
##' @param dt dataset. 
##' @param risk.set risk.set used for partial likelihood; here NULL. 
##' @param test.set validation data. 
##' @param time.var name of time variable.
##' @param X HAL design matrix. 
##' @param lambda.cv grid over which to choose penalization.
##' @param delta.var name of event type variable. 
##' @param delta.value type of event of interest here-
##' @return a
##' @seealso b
##' @examples c
##' @export 
##' @author Helene C. W. Rytgaard <hely@@biostat.ku.dk>
lebesgue.loss.fun <- function(train.fit, dt, risk.set=NULL, test.set, time.var=NULL, X=NULL, lambda.cv=NULL,
                              delta.var="delta", delta.value=1) {

    tmp <- copy(dt)

    tmp[id %in% test.set, fit.lambda:=exp(predict(train.fit, X[tmp$id %in% test.set,],
                                                  newoffset=0, s=lambda.cv))]
    tmp[id %in% test.set, fit.dLambda:=fit.lambda*risk.time]
    tmp[id %in% test.set, fit.Lambda:=cumsum(fit.dLambda), by="id"]

    tmp[id %in% test.set, dN:=1*(time.obs==grid.time & get(delta.var)==delta.value)]

    return(tmp[id %in% test.set & time.obs==grid.time, -sum(log(fit.lambda)*dN - fit.Lambda)])
}
