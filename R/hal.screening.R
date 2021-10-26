##' bla
##'
##' blub
##' @title atitle
##' @param covars covariates to include in HAL model. 
##' @param dt dataset. 
##' @param V number of folds in cross-validation. 
##' @param cut.one.way cut-offs for main effect indicators. 
##' @param method.risk option to pick different cross-validation schemes for cox models; basically it picks the
##' risk set for the partial likelihood. Should be chosen as 'test'.
##' @param cv.glmnet if TRUE, default of glmnet is used to find penalization (do not want this). 
##' @param seed random seed :). 
##' @param grouped option to set cross-validation method in glmnet. Keep grouped=TRUE. 
##' @param use.min combind with cv.glmnet, uses the minimal optimal value of penalization.
##' @param order order of interactions screened. 
##' @param cut.time.treatment cut-off for time/treatment interaction.
##' @param delta.var name of event type variable. 
##' @param delta.value type of event of interest here-
##' @param treatment name of treatment variable.
##' @param time.var name of time varaible. 
##' @param cut.time cut-off for time variable. 
##' @param mat dataset with information for HAL.
##' @param browse browser() for helene :). 
##' @param cut.two.way cut-offs for two-way interactions.
##' @return a
##' @seealso b
##' @examples c
##' @export 
##' @author Helene C. W. Rytgaard <hely@@biostat.ku.dk>
hal.screening <- function(covars, dt, V=5, cut.one.way=18, method.risk="test", cv.glmnet=FALSE,
                          seed=13349, grouped=TRUE, use.min=TRUE, order=1,
                          cut.time.treatment=NULL,
                          delta.var="delta", delta.value=1,
                          treatment=NULL, time.var=NULL, cut.time=5,
                          mat=NULL, browse=FALSE,
                          cut.two.way=5) {

    if (order==1) {
        fit <- fit.hal(covars=covars, dt=dt, V=V, cut.one.way=cut.one.way, method.risk=method.risk,
                       cut.time.treatment=cut.time.treatment,
                       delta.var=delta.var, delta.value=delta.value,
                       mat=mat, browse=browse,
                       cv.glmnet=cv.glmnet, treatment=treatment, time.var=time.var, cut.time=cut.time,
                       seed=seed, grouped=grouped, use.min=use.min)
        picked.covars <- covars[sapply(covars, function(covar) length(grep(covar, coef(fit$fit)@Dimnames[[1]][coef(fit$fit)@i+1]))>0)]
        if (length(picked.covars)==0) {
            for (jj in 1:10) {
                fit <- fit.hal(covars=covars, dt=dt, V=V, cut.one.way=cut.one.way*2, method.risk=method.risk,
                               cut.time.treatment=cut.time.treatment,
                               delta.var=delta.var, delta.value=delta.value,
                               mat=mat, browse=browse, lambda.cvs=fit$cve$min$lambda.cv/3,
                               cv.glmnet=cv.glmnet, treatment=treatment, time.var=time.var, cut.time=cut.time,
                               seed=seed, grouped=grouped, use.min=use.min)
                picked.covars <- covars[sapply(covars, function(covar) length(grep(covar, coef(fit$fit)@Dimnames[[1]][coef(fit$fit)@i+1]))>0)]
                if (length(picked.covars)>0) break
            }
        }
        return(picked.covars)
    } else {
        two.way <- expand.grid(c(treatment, covars), c(treatment, covars), stringsAsFactors=FALSE)
        fit <- fit.hal(covars=covars, dt=dt, V=V, cut.one.way=cut.one.way, method.risk=method.risk, cv.glmnet=cv.glmnet,
                       seed=seed, grouped=grouped, use.min=use.min, treatment=treatment,
                       cut.time.treatment=cut.time.treatment,
                       delta.var=delta.var, delta.value=delta.value,
                       mat=mat, browse=browse,
                       time.var=time.var, cut.time=cut.time,
                       two.way=two.way, cut.two.way=cut.two.way)$fit
        two.way.out <- two.way[apply(two.way, 1, function(row2) {
            if (row2[1]!=row2[2]) {
                (length(grep(row2[1], grep(paste0(":", row2[2]), coef(fit)@Dimnames[[1]][coef(fit)@i+1], value=TRUE)))>0)
            } else return(FALSE)
        }),, drop=FALSE]
        if (is.na(two.way.out[1,1])) two.way.out <- cbind(var1="", var2="")
        return(two.way.out)
        }
    
}
