##' bla
##'
##' blub
##' @title atitle
##' @param covars covariates to include in HAL model. 
##' @param cut.one.way cut-offs for main effect indicators. 
##' @param dt dataset. 
##' @param treatment name of treatment variable. 
##' @param Y outcome object. 
##' @param delta.var name of event type variable. 
##' @param delta.value type of event of interest here. 
##' @param two.way two-way interactions of variables to include. 
##' @param cut.two.way cut-offs for two-way interactions.
##' @param cut.time.treatment cut-off for time/treatment interaction.
##' @param predict time-horizon to make prediction (temporary option). 
##' @param pseudo.dt dataset with information for HAL. 
##' @param time.var name of time-variable. 
##' @param cut.time cut-off for time variable. 
##' @return a
##' @seealso b
##' @examples c
##' @export 
##' @author Helene C. W. Rytgaard <hely@@biostat.ku.dk>
basis.fun <- function(covars, cut.one.way, dt, treatment=NULL,
                      Y=paste0("(", delta.var, "==", delta.value, ")"),
                      delta.var="delta", delta.value=1,
                      two.way=cbind(var1="", var2=""), cut.two.way=5,
                      cut.time.treatment=NULL, predict=FALSE, pseudo.dt=NULL, 
                      time.var=NULL, cut.time=5) {

    if (length(covars)>0) {
        dt[, (covars):=lapply(.SD, function(x) {
            if (is.character(x)) {
                return(as.numeric(as.factor(x)))
            } else return(x)
        }), .SDcols=covars]
    }

    if (length(time.var)>0) {
        if (length(pseudo.dt)>0) {
            grid.times <- c(0, indicator.basis.fun(dt, time.var, cut.time, return.grid=TRUE), Inf)
            if (predict) grid.times <- sort(c(predict, grid.times))
            pseudo.dt[, grid.period:=as.numeric(cut(get(time.var), grid.times, include.lowest=TRUE, right=FALSE))]
            pseudo.dt <- unique(pseudo.dt, by=c("id", "grid.period"))
            pseudo.dt[, grid.time:=grid.times[grid.period]]
            if (!predict) pseudo.dt <- pseudo.dt[grid.time<time.obs]
            pseudo.dt[, grid.time:=grid.times[grid.period+1]]
            if (!predict) pseudo.dt[time.obs<=grid.time, grid.time:=time.obs]
            pseudo.dt[, (time.var):=grid.time]
        } else {
            pseudo.dt <- copy(dt)        
            if (length(time.var)>0) {
                grid.times <- c(0, indicator.basis.fun(dt, time.var, cut.time, return.grid=TRUE), Inf)
                if (predict) grid.times <- sort(c(predict, grid.times))
                pseudo.dt <- do.call("rbind", lapply(1:length(grid.times), function(jj) {
                    tmp <- copy(pseudo.dt)[, grid.period:=jj]
                }))
            }
            pseudo.dt[, grid.time:=grid.times[grid.period]]
            if (!predict) pseudo.dt <- pseudo.dt[grid.time<get(time.var)]
            pseudo.dt <- pseudo.dt[order(id, get(time.var))]
            pseudo.dt[, grid.time:=grid.times[grid.period+1]]
            if (!predict) pseudo.dt[get(time.var)<=grid.time, grid.time:=get(time.var)]
            pseudo.dt[, time.obs:=get(time.var)]
            pseudo.dt[, (time.var):=grid.time]
        }
        pseudo.dt <- pseudo.dt[!is.na(time)]
    } else {
        pseudo.dt <- copy(dt)
    }
    X <- Matrix(
        model.matrix(formula(paste0(
            Y, "~-1",
            ifelse(length(treatment)>0, paste0("+", treatment), ""),
            ifelse(length(time.var)>0, paste0("+",
                                              paste0(indicator.basis.fun(pseudo.dt, "grid.period", cut.time), collapse="+")), ""),
            ifelse(length(time.var)>0 & cut.time.treatment>0,
                   paste0("+", paste0(paste0(indicator.basis.fun(pseudo.dt, "grid.period", cut.time.treatment), ":", treatment, "+"), collapse="")), ""),
            ifelse(length(covars)>0 & cut.one.way>2,
                   paste0("+", paste0(sapply(covars, function(covar) paste0(indicator.basis.fun(pseudo.dt, covar, cut.one.way), collapse="+")), collapse="+")),
                   ""),
            ifelse(two.way[1,1]!="" & cut.two.way>2,
                   paste0(apply(two.way, 1, function(row) {
                       if (row[1]!=row[2]) {
                           if (any(row==treatment)) {
                               paste0("+", paste0(apply(expand.grid(ifelse(row[1]==treatment, treatment, indicator.basis.fun(dt, row[1], cut.two.way)),
                                                                    ifelse(row[2]==treatment, treatment, indicator.basis.fun(dt, row[2], cut.two.way))), 1,
                                                        function(row2) paste0(row2, collapse=":")), collapse="+"))
                           } else {
                               paste0("+", paste0(apply(expand.grid(indicator.basis.fun(dt, row[1], cut.two.way),
                                                                    indicator.basis.fun(dt, row[2], cut.two.way)), 1,
                                                        function(row2) paste0(row2, collapse=":")), collapse="+"))
                           }
                       } else return("")
                   }), collapse=""), "")
        )), 
        data=pseudo.dt), sparse=FALSE)

    if (length(time.var)>0) {
        x.vector <- apply(X, 1, function(x) paste0(x, collapse=","))
        pseudo.dt[, x:=x.vector]
        pseudo.dt[, risk.time:=grid.time-c(0, grid.time[-.N]), by="id"]
        return(list(X=X, pseudo.dt=pseudo.dt))
    } else {
        return(X)
    }

}


indicator.basis.fun <- function(mat, xvar, xcut, type="obs", return.grid=FALSE, seed=220) {
    set.seed(seed)
    if (type=="obs") {
        xvar.values <- mat[, sort(unique(get(xvar)))]
        xvar.pick <- seq(1, length(xvar.values), length=min(xcut, length(xvar.values)))
        xgrid <- xvar.values[xvar.pick][-c(1,xcut)]
    } else {
        xgrid <- round(seq(mat[, min(get(xvar))],
                           mat[, max(get(xvar))],
                           length=xcut)[-c(1,xcut)], 2)
    }
    if (return.grid) return(c(unique(xgrid))) else return(paste0("(", xvar, ">=", unique(xgrid), ")"))
}
