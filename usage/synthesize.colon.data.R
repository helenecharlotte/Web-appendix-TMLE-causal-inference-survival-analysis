fit.colon.fun <- function(formula.1=Surv(time, status==1)~rx+sex+nodes+age+obstruct+perfor+adhere+extent+surg,
                          formula.2=NULL,
                          formula.0=Surv(time, status==0)~rx+sex+nodes+age+obstruct+perfor+adhere+extent+surg,
                          formula.treat=trt~sex+age+nodes+obstruct+perfor+adhere+extent+surg,
                          initialize.scale.1=FALSE,
                          initialize.scale.2=FALSE,
                          initialize.scale.0=TRUE, #-- otherwise weibull for cens did not converge
                          d=colon) {

    covars <- unique(na.omit(gsub(" ", "", c(strsplit(as.character(formula.treat)[3], "\\+")[[1]],
                                             strsplit(as.character(formula.0)[3], "\\+")[[1]],
                                             strsplit(as.character(formula.1)[3], "\\+")[[1]],
                                             strsplit(as.character(formula.2)[3], "\\+")[[1]]))))

    covars.squared <- unlist(grep("\\.squared", covars))
    
    if (length(covars.squared)>0) {
        for (jj in covars.squared) {
            d[, (covars[jj]):=get(gsub("\\.squared", "", covars[jj]))^2]
        }
    }

    formula.treat.character <- as.character(formula.treat)
    d[, (paste0("num.", formula.treat.character[2])):=as.numeric(get(formula.treat.character[2]))]

    fit.treat <- fit.categorical(as.formula(paste0(paste0("num.", formula.treat.character[2]),
                                                   formula.treat.character[1],
                                                   formula.treat.character[3])), data=d, pseudo.out=TRUE)
    
    if (initialize.scale.1) {
        fit.weibull.1 <- survreg(as.formula(paste0(as.character(formula.1)[2], "~1")),
                                 data=na.omit(d), dist="weibull")
        fit.weibull.1 <- survreg(formula.1,
                                 data=na.omit(d), dist="weibull", scale=fit.weibull.1$scale)
    } else {
        fit.weibull.1 <- survreg(formula.1,
                                 data=na.omit(d), dist="weibull")
    }

    if (length(formula.2)>0) {
        if (initialize.scale.2) {
            fit.weibull.2 <- survreg(as.formula(paste0(as.character(formula.2)[2], "~1")),
                                     data=na.omit(d), dist="weibull")
            fit.weibull.2 <- survreg(formula.2,
                                     data=na.omit(d), dist="weibull", scale=fit.weibull.2$scale)
        } else {
            fit.weibull.2 <- survreg(formula.2,
                                     data=na.omit(d), dist="weibull")
        }
    }

    if (initialize.scale.0) {
        fit.weibull.0 <- survreg(as.formula(paste0(as.character(formula.0)[2], "~1")),
                                 data=na.omit(d), dist="weibull")
        fit.weibull.0 <- survreg(formula.0,
                                 data=na.omit(d), dist="weibull", scale=fit.weibull.0$scale)
    } else {
        fit.weibull.0 <- survreg(formula.0,
                                 data=na.omit(d), dist="weibull")
    }

    if (length(formula.2)>0) {
        return(list(fit.treat=fit.treat, fit.weibull.1=fit.weibull.1, fit.weibull.2=fit.weibull.2, fit.weibull.0=fit.weibull.0))
    } else {
        return(list(fit.treat=fit.treat, fit.weibull.1=fit.weibull.1, fit.weibull.0=fit.weibull.0))
    }

}

synthesize.colon.fun <- function(fit.colon, n=nrow(na.omit(d)), d=colon, get.true.value=NULL, tau=1000,
                                 name.treat="rx", event.name=NULL, time.name=NULL) {

    covars <- unique(unlist(lapply(fit.colon, function(fit) names(coef(fit))[-1])))
    covars <- covars[-unlist(lapply(c("pseudo", name.treat), function(x) grep(x, covars)))]

    d2 <- cbind(na.omit(d), fit.colon$fit.treat[[2]])
    
    covars2 <- names(d2)[unlist(sapply(names(d2), function(x) length(unique(grep(x, covars)>0))>0))]
    covars2 <- covars2[!covars2 %in% paste0(1:10)]
    
    sim.d <- setDT(d2)[sample(1:nrow(d2), n, replace=TRUE), c(covars2, colnames(fit.colon$fit.treat[[2]])), with=FALSE]
    
    if (length(get.true.value)>0) {
        sim.d[, (paste0(name.treat, ".num")):=get.true.value]
        sim.d[, (name.treat):=factor(get(paste0(name.treat, ".num")),
                                     levels=0:(length(levels(d2$rx))-1),
                                     labels=levels(d2$rx))]
    } else {
        sim.d[, (paste0(name.treat, ".num")):=as.numeric(sapply(1:nrow(sim.d), function(ii) sample(colnames(fit.colon$fit.treat[[2]]), 1, prob=sim.d[ii,names(sim.d)%in%colnames(fit.colon$fit.treat[[2]]), with=FALSE])))]
        sim.d[, (name.treat):=factor(get(paste0(name.treat, ".num")), labels=levels(d2$rx))]
    }

    shape.1 <- 1/fit.colon$fit.weibull.1$scale
    scale.1 <- exp(predict(fit.colon$fit.weibull.1, newdata=sim.d, type="lp"))

    shape.0 <- 1/fit.colon$fit.weibull.0$scale
    scale.0 <- exp(predict(fit.colon$fit.weibull.0, newdata=sim.d, type="lp"))

    time.1 <- rweibull(n, shape=shape.1, scale=scale.1)
    time.0 <- rweibull(n, shape=shape.0, scale=scale.0)

    if ("fit.weibull.2" %in% names(fit.colon)) {
        shape.2 <- 1/fit.colon$fit.weibull.2$scale
        scale.2 <- exp(predict(fit.colon$fit.weibull.2, newdata=sim.d, type="lp"))
        time.2 <- rweibull(n, shape=shape.2, scale=scale.2)
        sim.d[, time:=sapply(1:n, function(i) min(time.2[i], time.1[i], time.0[i]))]
        sim.d[, status:=1*(time.1<=time.0 & time.1<=time.2) +
                    2*(time.2<=time.0 & time.2<time.1)]
        if (length(get.true.value)>0) {
            status.uncensored <- 1*(time.1<=time.2) + 2*(time.2<time.1)
            return(c(F1=mean(time.1<=tau & status.uncensored==1),
                     F2=mean(time.1<=tau & status.uncensored==2)))
        } 
    } else {
        sim.d[, time:=sapply(1:n, function(i) min(time.1[i], time.0[i]))]
        sim.d[, status:=1*(time.1<=time.0)]
        if (length(get.true.value)>0) {
            return(mean(time.1<=tau))
        } 
    }

    if (length(time.name)>0) {
        names(sim.d)[names(sim.d)=="time"] <- time.name
    }
    if (length(event.name)>0) {
        names(sim.d)[names(sim.d)=="status"] <- event.name
    }
        
    return(sim.d[, -colnames(fit.colon$fit.treat[[2]]), with=FALSE])
}
