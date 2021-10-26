##' bla
##'
##' blub
##' @title atitle
##' @param loss.fun should just be set to cox.loss.fun.
##' @param dt dataset.
##' @param V number of folds in cross-validation.
##' @param seed random seed :). 
##' @param delta.var name of event type variable. 
##' @param delta.value type of event of interest here-
##' @param treatment name of treatment variable. 
##' @param change.points specified if there is a changepoint in the effect of treatment across time.
##' @param cox.models a list of Cox models to be compared with cross-validation.
##' @return a
##' @seealso b
##' @examples c
##' @export 
##' @author Helene C. W. Rytgaard <hely@@biostat.ku.dk>
cox.sl <- function(loss.fun, dt, V=5, seed=19192,
                   delta.var=NULL, delta.value=NULL, treatment="A",
                   change.points=(0:12)/10, 
                   cox.models=list(mod1=list(Surv(time, delta==1)~A+L1+L2+L3),
                                   mod2=list(Surv(time, delta==1)~A+L1.squared+L2+L3),
                                   mod3=list(Surv(time, delta==1)~L2.squared+A+L1.squared+L2+L3),
                                   mod4=list(Surv(time, delta==1)~A+L1.squared),
                                   mod5=list(Surv(time, delta==1)~A*L1+L2+L3),
                                   mod6=list(Surv(time, delta==1)~A*L1.squared+L2+L3))) {

    if (length(delta.var)==0) {
        form <- as.character(cox.models[[1]][[1]])[2]
        delta <- str_split(form, ",", simplify=TRUE)[,2]
        if (length(grep("==", delta))>0) {
            delta.var <- gsub(" ", "", str_split(delta, "==", simplify=TRUE)[,1])
            if (length(delta.value)==0) delta.value <- as.numeric(gsub(")", "", str_split(delta, "==", simplify=TRUE)[,2]))
        } else {
            delta.var <- gsub(" ", "", gsub(")", "", delta))
            if (length(delta.value)==0) delta.value <- 1
        }
    }

    if (length(delta.value)>0) {
        if (delta.value!=as.numeric(
                             gsub(")", "", gsub("\\=\\=", "", str_split(as.character(cox.models[[1]][[1]])[2],
                                                                        delta.var, simplify=TRUE)[1,2])))) {
            cox.models <- lapply(cox.models, function(cox.model) {
                form <- as.character(cox.model[[1]])[2]
                delta.form <- paste0(delta.var, "==", delta.value, ")")
                return(list(as.formula(paste0(paste0(str_split(form, ",", simplify=TRUE)[1,1], ",", delta.form),
                                              "~", as.character(cox.model[[1]])[3]))))
            })
        }
    }

    cox.cve <- lapply(cox.models, function(cox.model) {
        if (length(grep(".squared", as.character(cox.model[[1]])[3]))>0) {
            squared.vars <- c(str_split(as.character(cox.model[[1]])[3], ".squared", simplify=TRUE))
            squared.vars <- sapply(squared.vars[-length(squared.vars)], function(squared.var) {
                var.out <- c(str_split(gsub(" ", "", squared.var), "\\+|\\*", simplify=TRUE))
                return(var.out[length(var.out)])
            })
            for (square.var in squared.vars)
                dt[, (paste0(square.var, ".squared")):=get(square.var)^2]
        }
        if (length(cox.model)>1) change.points <- cox.model[[2]] #else change.points <- NULL
        if (length(change.points)>1) {
            cve.tmp <- sapply(change.points, function(change.point)
                cv.fun(loss.fun=cox.loss.fun, dt=dt, cox.model=cox.model[[1]],
                       delta.var=delta.var, delta.value=delta.value,
                       change.point=change.point, treatment=treatment))
            names(cve.tmp) <- paste0("changepoint=", change.points)
            return(cve.tmp)
        } else {
            return(cv.fun(loss.fun=cox.loss.fun, dt=dt, cox.model=cox.model[[1]],
                          delta.var=delta.var, delta.value=delta.value,
                          change.point=change.points, treatment=treatment))
        }
    })

    picked.model <- unlist(cox.cve)[unlist(cox.cve)==min(unlist(cox.cve))]

    if (length(grep("\\.changepoint", names(picked.model)))>0) {
        picked.cox.model <- str_split(names(picked.model), "\\.changepoint", simplify=TRUE)[1,1]
        picked.change.point <- as.numeric(gsub("\\=", "", str_split(names(picked.model), "\\.changepoint", simplify=TRUE)[1,2]))
        picked.cox.model <- list(form=cox.models[[picked.cox.model]][[1]],
                                 change.point=picked.change.point, 
                                 cve=picked.model[[1]])
    } else {
        picked.cox.model <- list(form=cox.models[[gsub("\\.cve", "", names(picked.model))]][[1]],
                                 cve=picked.model[[1]])
    }

    return(list(picked.cox.model=picked.cox.model, 
                cox.cve.all=cox.cve))

}
