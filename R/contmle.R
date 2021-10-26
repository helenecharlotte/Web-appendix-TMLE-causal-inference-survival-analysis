contmle <- function(dt,
                    #-- outcome model;
                    estimation=list("outcome"=list(fit=c("sl", "sl", "hal", "km"),
                                                   model=Surv(time, delta==1)~A+L1+L2+L3,#A+L1.squared,
                                                   changepoint=NULL),
                                    "cens"=list(fit=c("sl", "sl", "hal", "km"),
                                                model=Surv(time, delta==0)~L2+L3+A*L1,
                                                changepoint=NULL)#,
                                    #"cr"=list(fit=c("cox", "sl", "hal", "km"),
                                    #          model=Surv(time, delta==2)~A+L1+L2+L3,
                                    #         changepoint=NULL)                                          
                                    ),
                    #-- when there are competing risks, what is the target?
                    target=1, # if =0, only survival is targeted
                    #-- use poisson-based or cox-based hal?
                    hal.screening=FALSE, 
                    #-- use iterative or one-step tmle; (for competing risks, one-step is default)
                    one.step=FALSE, deps.size=0.1, no.small.steps=200,
                    push.criterion=FALSE,
                    iterative=FALSE,
                    #-- option to not update in directions currently already solved;
                    no.update.if.solved=FALSE, 
                    #-- target cause-specific hazards separately or together in one-step tmle; 
                    separate.cr=FALSE,
                    #-- simulatenous inference;
                    simultaneous.ci=FALSE,
                    #-- treatment model;
                    treat.model=A~L1+L2+L3,
                    #-- cut-off for truncation of censoring weights;
                    cut.off.cens=1e-3,
                    #-- use weighted norm in multivariate one-step;
                    weighted.norm=c(FALSE, "sigma", "Sigma"),
                    #-- treatment effect of interest; 
                    treat.effect=c("1", "0", "ate", "stochastic"),
                    #-- specify stochastic intervention; 
                    pi.star.fun=function(L) L$L1*0.2,
                    #-- time-point(s) of interest;
                    tau=c(1.2),
                    #-- increasing grid? ;
                    use.observed.times=FALSE, length.times=length(tau), check.times.size=NULL,#100,
                    delta.min=0.01, check.min=min(tau)+0.1, check.max=max(tau)-0.1, 
                    #-- pick super learning loss (should not change this);
                    sl.method=3, #DO NOT CHANGE.
                    #-- number of folds in cross-validation;
                    V=5,
                    #-- specify penalization in hal? 
                    lambda.cv=NULL,
                    maxit=1e3, 
                    #-- specify grid over which to pick penalization in hal by cross-validation; 
                    #lambda.cvs=seq(0.00000000001, 0.1, length=50),
                    init.lambda.cvs=c(sapply(1:5, function(jjj) (9:1)/(10^jjj))),
                    lambda.cvs=seq(0.0000001, 0.01, length=50),#seq(0, 0.008, length=51)[-1],
                    lambda.grid.size=50,
                    #-- penalize time indicators in hal? 
                    penalize.time=FALSE,
                    #-- pick grid for indicators in hal;
                    cut.covars=8, cut.time=10,
                    cut.time.A=10,
                    cut.L.A=8, cut.two.way=5,
                    #-- maximum number of iterations in iterative tmle; 
                    maxIter=10,
                    verbose=FALSE, verbose.sl=FALSE, check.sup=FALSE,
                    #-- for comparison; output kaplan-meier and hr; 
                    output.km=FALSE, only.km=FALSE, only.cox.sl=FALSE, output.RR=FALSE,
                    #-- dont use;
                    output.mat=FALSE, save.sl.pick=FALSE, 
                    #-- models incorporated in super learner;
                    sl.change.points=0,#(0:12)/10, 
                    sl.models=NULL
                    ## list(mod1=list(Surv(time, delta==1)~A+L1+L2+L3),
                    ##      mod2=list(Surv(time, delta==1)~A+L1.squared+L2+L3),
                    ##      mod3=list(Surv(time, delta==1)~L2.squared+A+L1.squared+L2+L3),
                    ##      mod4=list(Surv(time, delta==1)~A+L1.squared),
                    ##      mod5=list(Surv(time, delta==1)~A*L1+L2+L3),
                    ##      mod6=list(Surv(time, delta==1)~A*L1.squared+L2+L3))
                    ) {

    dt[, id:=1:.N]

    if (verbose) verbose.sl <- TRUE
   
    not.fit.list <- list()

    #-- names of time variable and event (delta) variable; 
    time.var <- gsub("Surv\\(", "", unlist(strsplit(as.character(estimation[[1]][["model"]])[2], ","))[1])
    delta.var <- gsub(" ", "",
                      gsub("==", "",
                           gsub("[0-9]+\\)", "",
                                unlist(strsplit(as.character(estimation[[1]][["model"]])[2], ","))[2])))

    #-- add event value to list of estimation, and only include those observed; 
    estimation <- lapply(estimation, function(x) {
        event <- as.numeric(gsub("\\D", "", unlist(strsplit(as.character(x[["model"]])[2], ","))[2]))
        x[[length(x)+1]] <- event
        names(x)[length(x)] <- "event"
        if (!(event %in% dt[, unique(get(delta.var))]))
            message(paste0("model specified for ", delta.var, "=", event,
                           ", but there were no observations")) else return(x)
    })

    estimation <- estimation[sapply(estimation, function(each) length(each)>0)]

    outcome.index <- unlist(lapply(1:length(estimation), function(each) {
        event <- estimation[[each]][["event"]]
        each[event>0]
    }))

    if (all(target==0)) {
        target <- 1:length(outcome.index)
        target.S <- TRUE
    } else {
        target.S <- FALSE
    }

    #-- is model specified multiple times for one event type?
    events <- unlist(lapply(estimation, function(each) each[["event"]]))
    if (any(table(events)>1)) {
        stop(paste0("multiple models specified for ", delta.var, "=",
                    paste0(names(table(events))[table(events)>1], collapse=",")))
    }

    #-- a missing model for one of the deltas?
    delta.missing <- dt[, unique(get(delta.var))][!(dt[, unique(get(delta.var))] %in% unlist(lapply(estimation, function(each) each[["event"]])))]
    if (length(delta.missing)) {
        warning(paste0("No model specified for event type=", paste0(delta.missing, collapse=","),
                       "; will use the one specified for event type=",
                       estimation[[1]][["event"]]))
        for (delta in delta.missing) {
            estimation[[length(estimation)+1]] <- estimation[[1]]
            estimation[[length(estimation)]][["event"]] <- delta
            estimation[[length(estimation)]][["model"]] <-
                as.formula(paste0(gsub(estimation[[1]][["event"]], delta, estimation[[length(estimation)]][["model"]][2]),
                                  estimation[[length(estimation)]][["model"]][1],
                                  estimation[[length(estimation)]][["model"]][3]))
        }
    }

   
    #-- 0 -- some initializations:

    #-- are there competing risks?
    if (length(dt[get(delta.var)>0, unique(get(delta.var))])>1) cr <- TRUE else cr <- FALSE
    if (!cr) target <- 1 else if (length(dt[get(delta.var)>0, unique(get(delta.var))])<3)
                             target <- target[target<3]

    #-- get number of subjects:
    n <- length(dt[, unique(id)])
    
    #-- get treatment colname:
    A.name <- as.character(treat.model)[2]

    #-- list of covariates
    covars <- NULL
    for (mod in c(lapply(sl.models, function(x) x[[1]]),
                  unlist(lapply(estimation, function(x) x[["model"]])))) {
        mod3 <- as.character(mod)[3]
        covars <- unique(c(covars, unlist(strsplit(gsub("\\+", " ",
                                                        mod3), " "))))
        covars <- covars[!covars%in%c(A.name, "", "*")]
        if (length(grep(".squared", mod3))>0) {
            names.squared <- unique(gsub(".squared", "",
                                         grep(".squared", unlist(strsplit(gsub("\\+", " ",
                                                                               mod3), " ")),
                                              value=TRUE)))
            for (col in names.squared)
                dt[, (paste0(col, ".squared")):=get(col)^2]
        }
        if (length(grep(".log", mod3))>0) {
            names.log <- unique(gsub(".log", "",
                                     grep(".log", unlist(strsplit(gsub("\\+", " ",
                                                                       mod3), " ")),
                                          value=TRUE)))
            for (col in names.log)
                dt[, (paste0(col, ".log")):=log(get(col))]
        }
    }
    if (length(grep(".squared", covars))>0) covars <- covars[-grep(".squared", covars)]
    if (length(grep(".log", covars))>0) covars <- covars[-grep(".log", covars)]
    if (length(grep("cut.", covars))>0) covars <- covars[-grep("cut.", covars)]

    covars <- covars[covars!="1"]
    for (covar in covars) {
        if (dt[, length(unique(get(covar)))]==1) covars <- covars[covars!=covar]
    }

    #if (verbose) print(covars)
    
    #-- get unique times in dataset
    unique.times <- sort(unique(dt[, get(time.var)]))
    unique.times2 <- sort(unique(dt[get(delta.var)>0, get(time.var)]))

    #-- intervals in tau where no obs?

    if (use.observed.times) {
        n <- nrow(dt)
        unique.times3 <- unique.times[unique.times<=max(tau)]
        tau <- unique.times3[(1:length(unique.times3)) %% floor(length(unique.times3)/length.times)==1]
        #set.seed(1)
        sort(unique.times3[sample(length(unique.times3), length.times)])
    } else {    
        test.tau <- findInterval(tau, unique.times2)
        if (FALSE & !length(test.tau)==length(unique(test.tau))) {
            tau <- na.omit(tau[(1:length(tau))[unique(findInterval(unique.times2, tau))+1]])
            warning("no observations between tau as specified, truncating tau")
        }
    }

    #-- which parameters are we interested in?
    if (dt[, length(unique(get(A.name)))]>2) {
        a <- treat.effect
    } else {
        if (treat.effect[1]=="1") a <- 1 else if (treat.effect[1]=="0") a <- 0 else a <- c(1, 0)
    }

    #-- initialize dataset to be used later; 
    dt2 <- NULL
    bhaz.cox <- do.call("rbind", lapply(a, function(aa) data.table(time=c(0, unique.times), A=aa)))
    setnames(bhaz.cox, "A", A.name)

    #-- if there is any of the outcome models that uses coxnet
    sl.models.tmp <- sl.models
    sl.models <- list()
    #-- add separate sl models when specified with, e.g., multiple changepoints
    for (k1 in 1:length(sl.models.tmp)) {
        if (length(sl.models.tmp[[k1]])>1) {
            for (k2 in 2:length(sl.models.tmp[[k1]])) {
                sl.models[[length(sl.models)+1]] <- c(sl.models.tmp[[k1]][1],
                                                      sl.models.tmp[[k1]][k2])
                if (length(sl.models.tmp[[k1]])>=3) {
                    names(sl.models)[length(sl.models)] <- paste0(names(sl.models.tmp)[k1], k2)
                } else {
                    names(sl.models)[length(sl.models)] <- names(sl.models.tmp)[k1]
                }
            }
        } else {
            sl.models[[length(sl.models)+1]] <- c(sl.models.tmp[[k1]][1])
            names(sl.models)[length(sl.models)] <- names(sl.models.tmp)[k1]
        }
    }

    #-- 2 -- estimate treatment propensity:
    if (dt[, length(unique(get(A.name)))]>2) {
        dt[, (paste0(A.name, ".num")):=as.numeric(get(A.name))]
        prob.A <- fit.categorical(as.formula(paste0(A.name, ".num", "~", as.character(treat.model)[3])),
                                  data=dt, pseudo.out=TRUE)[[2]]
        prob.A <- sapply(1:nrow(dt), function(ii) prob.A[ii, dt[ii, get(paste0(A.name, ".num"))]])        
    } else {
        prob.A <- predict(glm(as.formula(treat.model), data=dt), type="response")
    }

    if (verbose & dt[, length(unique(get(A.name)))]==2) print(summary(glm(as.formula(treat.model), data=dt)))
    if (verbose & dt[, length(unique(get(A.name)))]==2) print(summary(prob.A))

    #-- 3 -- estimation -- loop over causes (including censoring):

    for (each in 1:length(estimation)) {
        
        fit <- estimation[[each]][["fit"]][1]
        fit.model <- estimation[[each]][["model"]]
        if (any(names(estimation[[each]])=="changepoint"))
            fit.changepoint <- estimation[[each]][["changepoint"]] else fit.changepoint <- NULL
        fit.delta <- estimation[[each]][["event"]]
        fit.name <- names(estimation)[each]

        if (fit[1]=="sl" | fit[1]=="cox.hal.sl") { #-- cox-sl
           
            if (verbose) print(paste0("use sl for ", fit.name))

            set.seed(1)
            sl.pick <- suppressWarnings(
                cox.sl(loss.fun=cox.loss.fun, dt=dt, delta.var=delta.var,
                       treatment=A.name, V=V, delta.value=fit.delta,
                       cox.models=sl.models, change.points=sl.change.points))

            if (save.sl.pick) return(sl.pick)
            
            cve.sl.pick <- sl.pick$picked.cox.model$cve
            fit.model <- sl.pick$picked.cox.model$form

            #rm(dt.tmp)
        
            if (verbose.sl) print(paste0("model picked for ", fit.name, ": "))
            if (verbose) print(fit.model)
            
            estimation[[each]]$model <- fit.model

            if (any(names(sl.pick$picked.cox.model)=="change.point")) {
                fit.changepoint <- sl.pick$picked.cox.model$change.point
                if (fit.changepoint>0) {
                    if (verbose.sl) print(paste0("changepoint picked: ", fit.changepoint))
                    estimation[[each]]$changepoint <- fit.changepoint
                } else {
                    fit.changepoint <- NULL
                    estimation[[each]]$changepoint <- NULL
                }
            } else {
                estimation[[each]]$changepoint <- NULL
            }

        } else {
            
            sl.pick <- ""
            cve.sl.pick <- ""
            
            if (fit[1] %in% c("km")) { #-- later for hal or if uses km
                tmp.model <- as.character(fit.model)
                if (fit[1]=="km") tmp.model[3] <- "strata(A)" else tmp.model[3] <- "1"
                estimation[[each]]$model <- fit.model <- formula(paste0(tmp.model[2], tmp.model[1], tmp.model[3]))
                estimation[[each]]$changepoint <- fit.changepoint <- NULL
            }
        }

        estimation[[each]]$sl.pick <- fit.model
        estimation[[each]]$cve.sl.pick <- cve.sl.pick

        #-- changepoint?
        if (length(fit.changepoint)>0 & length(dt2)==0) {
            dt2 <- rbind(dt, dt)[order(id)]
        }

        #   fit.cox.fun <- function(mod, changepoint, fit, dt, dt2, dd=1, sl.pick="") {
        if (length(fit.changepoint)>0) { #-- if there is a change-point:
            delta1 <- abs(fit.delta-1)
            dt2[, time.indicator:=(get(time.var)<=fit.changepoint)]
            dt2[, (paste0("period", fit.delta)):=1:.N, by="id"]
            dt2[get(paste0("period", fit.delta))==1, (paste0("tstart", fit.delta)):=0]
            dt2[get(paste0("period", fit.delta))==1, (paste0("tstop", fit.delta)):=(get(time.var)<=fit.changepoint)*get(time.var)+
                                                         (get(time.var)>fit.changepoint)*fit.changepoint]
            dt2[get(paste0("period", fit.delta))==1, (paste0("tstart", fit.delta)):=0]
            dt2[get(paste0("period", fit.delta))==1, (paste0("tstop", fit.delta)):=(get(time.var)<=fit.changepoint)*get(time.var)+
                                                         (get(time.var)>fit.changepoint)*fit.changepoint]
            dt2[get(paste0("period", fit.delta))==2, (paste0("tstart", fit.delta)):=fit.changepoint]
            dt2[get(paste0("period", fit.delta))==2, (paste0("tstop", fit.delta)):=get(time.var)]
            dt2[get(paste0("period", fit.delta))==1 & !time.indicator, (delta.var):=delta1]
            mod1 <- as.character(fit.model)
            mod2 <- paste0(gsub(substr(mod1[2], which(strsplit(mod1[2], "")[[1]]=="(")+1,
                                       which(strsplit(mod1[2], "")[[1]]==",")-1), paste0("tstart", fit.delta, ", tstop", fit.delta), mod1[2]),
                           "~", 
                           gsub(paste0("\\+", A.name, "\\+"), "", gsub(paste0("\\+", A.name, " "), "", gsub(" ", "", paste0("I((period", fit.delta,"==1)&(", A.name, "==1))",
                                                                                                                            " + I((period", fit.delta, "==2)&(", A.name, "==1))", " + ",
                                                                                                                            paste0("+",mod1[3]))))))
            fit.cox <- coxph(formula(mod2), data=dt2[!time.indicator | get(paste0("period", fit.delta))==1])
        } else { #-- if there is not a change-point:
            if (fit[1]=="sl" & length(grep("coxnet", sl.pick))>0) {
                X <- model.matrix(as.formula(deparse(fit.model)), data=dt)
                y <- dt[, Surv(get(time.var), get(delta.var)==fit.delta)]
                fit.cox <- glmnet(x=X, y=y, family="cox", maxit=1000,
                                  lambda=fit.penalty)
            } else {
                fit.cox <- coxph(as.formula(paste0(deparse(fit.model), collapse="")), data=dt)
            }
        }

        estimation[[each]]$fit.cox <- fit.cox
        
        if (verbose) print(fit.cox)

        #-- 6 -- get baseline hazard:

        if (fit[1]=="km") {
            tmp <- suppressWarnings(setDT(basehaz(fit.cox, centered=TRUE)))
            setnames(tmp, "strata", "A")
            tmp[, A:=as.numeric(gsub("A=", "", A))]
            bhaz.cox <- merge(bhaz.cox, 
                              rbind(do.call("rbind", lapply(a, function(aa) data.table(time=0, hazard=0, A=aa))),
                                    tmp),
                              by=c("time", "A"), all.x=TRUE)
            bhaz.cox[, hazard:=na.locf(hazard), by="A"]
            bhaz.cox[, (paste0("dhaz.", fit.delta)):=c(0, diff(hazard)), by="A"]
            setnames(bhaz.cox, "hazard", paste0("chaz", fit.delta))
        } else {
            if (fit[1]=="sl" & length(grep("coxnet", sl.pick))>0) { 
                basehaz <- glmnet_basesurv(dt[, get(time.var)],
                                           dt[, get(delta.var)==fit.delta], X, centered=TRUE)
                bhaz.cox <- merge(bhaz.cox, rbind(data.table(time=0, hazard=0),
                                                  data.table(time=basehaz$time,
                                                             hazard=basehaz$cumulative_base_hazard)),
                                  by="time", all.x=TRUE)
            } else {
                bhaz.cox <- merge(bhaz.cox, rbind(data.table(time=0, hazard=0),
                                                  suppressWarnings(setDT(basehaz(fit.cox, centered=TRUE)))),
                                  by="time", all.x=TRUE)
            }
            bhaz.cox[, hazard:=na.locf(hazard), by=A.name]
            bhaz.cox[, (paste0("dhaz", fit.delta)):=c(0, diff(hazard)), by=A.name]
            setnames(bhaz.cox, "hazard", paste0("chaz", fit.delta))
        }

    }

    
    #-- set names of bhaz.cox to match observed data
    setnames(bhaz.cox, c("time", A.name), c(time.var, A.name))
    
    #-- Xc -- get censoring survival one time-point back: 

    bhaz.cox[, chaz0.1:=c(0, chaz0[-.N])]
    
    #-- Y -- output Kaplan-Meier and/or crude HR?
    
    if (output.km) {
        if (cr & !target.S) {
            km.mod <- paste0(gsub(estimation[[1]][["event"]], "",
                                  gsub("\\=", "",
                                       gsub("Surv", "Hist", as.character(estimation[[1]][["model"]])[2]))),
                             "~", A.name)
            km.fit <- summary(prodlim(formula(km.mod), data=dt),
                              cause=target, times=tau, asMatrix=TRUE)$table
            if (length(a)==1) {
                km.est <- as.numeric(km.fit[km.fit[,2]==paste0(A.name, "=", a),,drop=FALSE][,"cuminc"])
                km.se <- as.numeric(km.fit[km.fit[,2]==paste0(A.name, "=", a),,drop=FALSE][,"se.cuminc"])
            } else {
                km.est <- as.numeric(km.fit[km.fit[,2]==paste0(A.name, "=", "1"),,drop=FALSE][,"cuminc"])-
                    (as.numeric(km.fit[km.fit[,2]==paste0(A.name, "=", "0"),,drop=FALSE][,"cuminc"]))
                km.se <- sqrt((as.numeric(km.fit[km.fit[,2]==paste0(A.name, "=", "1"),,drop=FALSE][,"se.cuminc"]))^2+
                              (as.numeric(km.fit[km.fit[,2]==paste0(A.name, "=", "0"),,drop=FALSE][,"se.cuminc"]))^2)
            }
        } else {
            if (target.S) {
                km.mod <- paste0("Hist(", time.var, ",", delta.var, "!=0)", "~", A.name)
            } else {
                km.mod <- paste0(gsub("Surv", "Hist", as.character(estimation[[1]][["model"]])[2]),
                                 "~", A.name)
            }
            km.fit <- summary(fit.km <- prodlim(formula(km.mod), data=dt),
                              times=tau, asMatrix=TRUE)$table
            if (length(a)==1) {
                if (cr) {
                    km.est <- as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", a),,drop=FALSE][,"surv"])
                } else {
                    km.est <- 1-as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", a),,drop=FALSE][,"surv"])
                }
                km.se <- as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", a),,drop=FALSE][,"se.surv"])
            } else {
                if (cr) {
                    km.est <- as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", "1"),,drop=FALSE][,"surv"])-
                        as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", "0"),,drop=FALSE][,"surv"])
                } else {
                    km.est <- 1-as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", "1"),,drop=FALSE][,"surv"])-
                        (1-as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", "0"),,drop=FALSE][,"surv"]))
                }
                km.se <- sqrt((as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", "1"),,drop=FALSE][,"se.surv"]))^2+
                              (as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", "0"),,drop=FALSE][,"se.surv"]))^2)
            }
        }
    }
   
    if (length(dt2)==0) {
        dt2 <- copy(dt)
    }

    if (only.km) return(lapply(1:length(target), function(each.index) {
        out <- rbind(km.est=km.est[((each.index-1)*length(tau)+1):(each.index*length(tau))],
                     km.se=km.se[((each.index-1)*length(tau)+1):(each.index*length(tau))])
        colnames(out) <- paste0("tau=", tau)
        return(out)
    }))

    #-- hal?
    any.hal <- unlist(lapply(estimation, function(each) length(grep("hal", each[["fit"]]))>0))
        
    if (!any(any.hal)) {
        bhaz.cox <- bhaz.cox[get(time.var)<=max(tau)]
    }

    #-- 8 -- add subject-specific information:

    dt2.a <- do.call("rbind", lapply(a, function(aa) {
        dt.tmp <- copy(dt2)
        dt.tmp[, (A.name):=aa]
    }))

    for (each in 1:length(estimation)) {

        if (estimation[[each]][["fit"]][1]=="sl" & length(grep("coxnet", estimation[[each]][["sl.pick"]]))>0) {
            X2.a <- model.matrix(as.formula(deparse(estimation[[each]][["model"]])), data=dt2.a)
            dt2.a[, (paste0("fit.cox", estimation[[each]][["event"]])):=
                        exp(predict(estimation[[each]][["fit.cox"]], newx=X2.a, type="link"))]
        } else {
            dt2.a[, (paste0("fit.cox", estimation[[each]][["event"]])):=
                        predict(estimation[[each]][["fit.cox"]], newdata=dt2.a, type="risk")]
        }

        dt2.a[get(paste0("fit.cox", estimation[[each]][["event"]]))>500, (paste0("fit.cox", estimation[[each]][["event"]])):=500]
            
    }
    
    mat <- do.call("rbind", lapply(1:n, function(i) {
        tmp <- cbind(id=i, bhaz.cox)
        tmp[, time.obs:=dt[id==i, get(time.var)]]
        tmp[, delta.obs:=dt[id==i, get(delta.var)]]
        tmp[, A.obs:=dt[id==i, get(A.name)]]
        tmp[, prob.A:=prob.A[i]]
        if (treat.effect[1]=="stochastic") {
            tmp <- merge(tmp, dt[, c("id", covars), with=FALSE], by="id")
            tmp[, pi.star:=pi.star.fun(.SD), .SDcols=covars]
            tmp[, pi.star:=pi.star*get(A.name)+(1-pi.star)*(1-get(A.name))]
            if (dt[, length(unique(get(A.name)))]>2) {
                tmp[, Ht:=-(pi.star / # treatment and censoring weights
                            ((prob.A)))]
            } else {
                tmp[, Ht:=-(pi.star / # treatment and censoring weights
                            ((prob.A^get(A.name) * (1-prob.A)^(1-get(A.name)))))]
            }
        } else {
            if (length(a)==2) {
                if (dt[, length(unique(get(A.name)))]>2) {
                    tmp[, Ht:=-((get(A.name)==1) - (get(A.name)==0)) / # treatment and censoring weights
                              ((prob.A))]
                } else {
                    tmp[, Ht:=-((get(A.name)==1) - (get(A.name)==0)) / # treatment and censoring weights
                              ((prob.A^get(A.name) * (1-prob.A)^(1-get(A.name))))]
                }
            } else {
                if (dt[, length(unique(get(A.name)))]>2) {
                    tmp[, Ht:=-((get(A.name)==a)) / # treatment and censoring weights
                              ((prob.A))] 
                } else {
                    tmp[, Ht:=-((get(A.name)==a)) / # treatment and censoring weights
                              ((prob.A^get(A.name) * (1-prob.A)^(1-get(A.name))))]
                }
            }            
        }
        for (each in 1:length(estimation)) {
            if (length(estimation[[each]][["changepoint"]])>0) {
                tmp[, (paste0("period", estimation[[each]][["event"]])):=
                          (get(time.var)<=estimation[[each]][["changepoint"]])*1+(get(time.var)>estimation[[each]][["changepoint"]])*2]
                tmp <- merge(tmp, dt2.a[id==i, c(paste0("period", estimation[[each]][["event"]]),
                                                 A.name,
                                                 paste0("fit.cox", estimation[[each]][["event"]])),
                                        with=FALSE], by=c(paste0("period", estimation[[each]][["event"]]), A.name))
            } else {
                tmp <- merge(tmp, unique(dt2.a[id==i, c(A.name,
                                                        paste0("fit.cox", estimation[[each]][["event"]])),
                                               with=FALSE]), by=c(A.name))
            }
        }

        if (any(unlist(lapply(estimation, function(x) x[["event"]]))==0)) {
            tmp[, surv.C1:=exp(-fit.cox0*chaz0.1)]
            if (tmp[, is.na(mean(is.na(surv.C1)))]) browser()
            if (tmp[, any(is.na(surv.C1))]) browser()
            if (tmp[, min(surv.C1)<cut.off.cens]) {
                tmp[surv.C1<cut.off.cens, surv.C1:=cut.off.cens]
                #warning(paste0("there seems to be positivity issues; truncated at level ",
                #               cut.off.cens))
            }
            tmp[, Ht:=Ht/surv.C1]
        } 
    }))

    if (verbose) paste0("min of censoring weights: ", mat[, min(surv.C1)])
   
    #-- 10 -- poisson-HAL used for initial:

    if (any(any.hal)) {

        set.seed(13444)

        count.hals <- 1

        mat <- merge(mat, dt[, c("id", c(covars[!covars%in%names(mat)])), with=FALSE], by="id")

        mat[, (covars):=lapply(.SD, function(x) {
            if (is.character(x)) {
                return(as.numeric(as.factor(x)))
            } else if (is.factor(x)) {
                return(as.numeric(x))
            } else return(x)
        }), .SDcols=covars]

        for (each in (1:length(estimation))[any.hal]) { 

            #-- initialize covariates
            covars1 <- covars
            
            fit.delta <- estimation[[each]][["event"]]
            fit.name <- names(estimation)[each]

            if (any(names(estimation[[each]])=="lambda.cvs")) {
                lambda.cvs.1 <- estimation[[each]][["lambda.cvs"]]
            } else {
                lambda.cvs.1 <- lambda.cvs
            }

            if (any(estimation[[each]]$fit=="cox.hal") | any(estimation[[each]]$fit=="cox.hal.sl")) { #--- only testing.
                
                if (hal.screening) { #--- first screening

                    (one.way.screening <- hal.screening(covars=covars, dt=dt, cut.one.way=5,# browse=TRUE,
                                                        mat=mat, delta.var=delta.var, delta.value=fit.delta,
                                                        treatment=A.name))

                    (two.way.screening <- hal.screening(covars=one.way.screening, dt=dt, cut.one.way=5,
                                                        delta.var=delta.var, delta.value=fit.delta,
                                                        treatment=A.name, order=2,
                                                        mat=mat))

                    if (verbose & length(one.way.screening)>0) print(paste0("variables picked by initial screening: ", paste0(one.way.screening, collapse=", ")))
                    
                    mat <- fit.hal(covars=one.way.screening, dt=dt, cut.one.way=15,
                                   mat=mat, delta.var=delta.var, V=V,
                                   intervention=a, 
                                   two.way=two.way.screening, cut.two.way=cut.two.way, 
                                   penalize.treatment=FALSE,
                                   delta.value=fit.delta,
                                   verbose=verbose, 
                                   predict=max(tau), treatment.prediction=A.name, 
                                   treatment=A.name)

                    if (FALSE) {
                        diff(mat[time<=tau, exp(-cumsum(dhaz1*fit.cox1))[.N], by=c("id", A.name)][, mean(V1), by=A.name][order(V1), V1])
                        diff(mat2[time<=tau, exp(-cumsum(dhaz1*fit.cox1))[.N], by=c("id", A.name)][, mean(V1), by=A.name][order(V1), V1])
                    }

                } else {

                    mat <- fit.hal(covars=covars, dt=dt, cut.one.way=15,
                                   mat=mat, delta.var=delta.var, V=V,
                                   intervention=a, 
                                   penalize.treatment=FALSE,
                                   verbose=verbose, delta.value=fit.delta,
                                   predict=max(tau), treatment.prediction=A.name, 
                                   treatment=A.name)
                    
                }
                
            } else {

                if (hal.screening) { #--- first screening
                  
                    print(paste0("EACH = ", each))

                    (one.way.screening <- hal.screening(covars=covars, dt=dt, cut.one.way=5, time.var="time",
                                                        cut.time=3, mat=mat, delta.var="delta.obs",
                                                        delta.value=fit.delta,
                                                        treatment="A.obs",
                                                        cut.time.treatment=cut.time.A))
                    
                    if (length(one.way.screening)>0) {
                        (two.way.screening <- hal.screening(covars=one.way.screening, dt=dt, cut.one.way=5,
                                                            time.var="time", cut.time=3, delta.var="delta.obs",
                                                            delta.value=fit.delta,
                                                            treatment="A.obs",
                                                            cut.time.treatment=cut.time.A, order=2,
                                                            mat=mat, cut.two.way=cut.two.way))                  
                        if (verbose & length(one.way.screening)>0) print(paste0("variables picked by initial screening: ", paste0(one.way.screening, collapse=", ")))
                        if (verbose & length(two.way.screening)>0) print(paste0("interactions picked by initial screening: "))
                        print(two.way.screening)
                    }
                   
                    if (FALSE) {
                        diff(mat2[time<=tau, exp(-fit.Lambda[.N]), by=c("id", A.name)][, mean(V1), by=A.name][order(V1), V1])
                        diff(mat[time<=tau, exp(-cumsum(dhaz1*fit.cox1))[.N], by=c("id", A.name)][, mean(V1), by=A.name][order(V1), V1])
                    }

                    if (FALSE) { ##-- 
                        cve1 <- lapply(list(NULL, 3, 4, 5), function(xx) fit.hal(covars=one.way.screening, dt=dt, cut.one.way=15,
                                                                                 time.var=time.var, cut.time=10,
                                                                                 mat=mat, delta.var="delta.obs", V=V,
                                                                                 intervention=a, 
                                                                                 verbose=FALSE, delta.value=fit.delta,
                                                                                 two.way=two.way.screening, cut.two.way=cut.two.way, 
                                                                                 penalize.treatment=FALSE, penalize.time=FALSE,
                                                                                 predict=max(tau), treatment.prediction=A.name, 
                                                                                 treatment="A.obs", cut.time.treatment=xx,
                                                                                 return.cve=TRUE)$cve$min$cve)
                        cut.time.A <- unlist(list(NULL, 3, 4, 5)[unlist(cve1)==min(unlist(cve1))])
                        cve2 <- lapply(list(NULL, 3, 4, 5, 7), function(xx) fit.hal(covars=one.way.screening, dt=dt, cut.one.way=15,
                                                                                    time.var=time.var, cut.time=10,
                                                                                    mat=mat, delta.var="delta.obs", V=V,
                                                                                    intervention=a, 
                                                                                    verbose=FALSE, delta.value=fit.delta,
                                                                                    two.way=two.way.screening, cut.two.way=xx, 
                                                                                    penalize.treatment=FALSE, penalize.time=FALSE,
                                                                                    predict=max(tau), treatment.prediction=A.name, 
                                                                                    treatment="A.obs", cut.time.treatment=cut.time.A,
                                                                                    return.cve=TRUE)$cve$min$cve)
                        cut.two.way <- unlist(list(NULL, 3, 4, 5, 7)[unlist(cve2)==min(unlist(cve2))])
                    }

                    mat <- fit.hal(covars=one.way.screening, dt=dt, cut.one.way=15,
                                   time.var=time.var, cut.time=10,
                                   mat=mat, delta.var="delta.obs", V=V,
                                   intervention=a, 
                                   verbose=verbose, delta.value=fit.delta,
                                   two.way=two.way.screening, cut.two.way=cut.two.way, 
                                   penalize.treatment=FALSE, penalize.time=FALSE,
                                   predict=max(tau), treatment.prediction=A.name, 
                                   treatment="A.obs", cut.time.treatment=cut.time.A)
                     
                    if (FALSE) {
                        mat2[time<=tau, exp(-cumsum(dhaz1*fit.cox1))[.N], by=c("id", A.name)][, mean(V1), by=A.name][order(V1), V1]
                        mat[time<=tau, exp(-cumsum(dhaz1*fit.cox1))[.N], by=c("id", A.name)][, mean(V1), by=A.name][order(V1), V1]
                    }
                        
                } else {

                    print(paste0("EACH = ", each))

                    mat <- fit.hal(covars=covars, dt=dt, cut.one.way=15, time.var=time.var, cut.time=10,
                                   mat=mat, delta.var="delta.obs", V=V,
                                   intervention=a, 
                                   verbose=verbose, delta.value=fit.delta,
                                   penalize.treatment=FALSE, penalize.time=FALSE,
                                   predict=max(tau), treatment.prediction=A.name, 
                                   treatment="A.obs", cut.time.treatment=cut.time.A)

                }                
            }

            if (FALSE) if (mat[["not.fit"]]) {
                warning(paste0("hal did not fit for ", fit.name, " (", delta.var, "=", fit.delta, ")"))
                not.fit.list[[length(not.fit.list)+1]] <-
                    paste0("hal did not fit for ", fit.name, " (", delta.var, "=", fit.delta, ")")
            }

            #if (!(count.hals<sum(any.hal)))
            #mat <- mat[["mat"]] #else mat[["not.fit"]] <- NULL

            count.hals <- count.hals+1
        }

        mat <- mat[get(time.var)<=max(tau)]

    }

    #-- report if censoring weights were truncated?
   
    if (mat[, min(surv.C1)<cut.off.cens]) {
        if (mat[get(time.var)<=time.obs, min(surv.C1)<cut.off.cens]) {
            not.fit.list[[length(not.fit.list)+1]] <-
                paste0(paste0("there seems to be positivity issues for censoring; truncated ",
                              sum(mat[get(time.var)<=time.obs, sum(unique(surv.C1<cut.off.cens)), by="id"][,2][[1]]), 
                              " observation(s), at level ",
                              cut.off.cens))
        }
        mat[surv.C1<=cut.off.cens, Ht:=Ht*surv.C1/cut.off.cens]
        mat[surv.C1<=cut.off.cens, surv.C1:=cut.off.cens]
    }

    if (mat[, min(prob.A)<cut.off.cens]) {
        if (mat[, min(prob.A)<cut.off.cens]) {
            not.fit.list[[length(not.fit.list)+1]] <-
                paste0(paste0("there seems to be positivity issues for treatment; truncated ",
                              sum(mat[, sum(unique(prob.A<cut.off.cens)), by="id"][,2][[1]]), 
                              " observation(s), at level ",
                              cut.off.cens))
        }
        if (dt[, length(unique(get(A.name)))]>2) {
            mat[prob.A<=cut.off.cens, Ht:=Ht*(prob.A)/cut.off.cens]
        } else {
            mat[prob.A<=cut.off.cens, Ht:=Ht*(prob.A^get(A.name)*(1-prob.A)^{1-get(A.name)})/cut.off.cens]
        }
        mat[prob.A<=cut.off.cens, prob.A:=cut.off.cens]
    }


    #-- 9 -- compute clever covariates:

    mat[, surv.t:=1]
    for (each in outcome.index) {
        fit.delta <- estimation[[each]][["event"]]
        mat[, surv.t:=surv.t*exp(-cumsum(get(paste0("dhaz", fit.delta))*
                                         get(paste0("fit.cox", fit.delta)))),
            by=c("id", A.name)]
    }

    mat[, surv.t1:=c(0, surv.t[-.N]), by=c("id", A.name)]

    for (kk in 1:length(tau)) {
        mat[, (paste0("surv.tau", kk)):=
                  surv.t[get(time.var)==max(get(time.var)[get(time.var)<=tau[kk]])],
            by=c("id", A.name)]
    }
   
    if (cr) {       
        for (each in outcome.index) {
            fit.delta <- estimation[[each]][["event"]]
            mat[, (paste0("F", fit.delta, ".t")):=cumsum(surv.t*get(paste0("dhaz", fit.delta))*
                                                         get(paste0("fit.cox", fit.delta))),
                by=c("id", A.name)]
            for (kk in 1:length(tau)) {
                mat[, (paste0("F", fit.delta, ".tau", kk)):=
                          get(paste0("F", fit.delta, ".t"))[get(time.var)==max(get(time.var)[get(time.var)<=tau[kk]])],
                    by=c("id", A.name)]
                for (each2 in outcome.index) {
                    fit.delta2 <- estimation[[each2]][["event"]]
                    if (fit.delta==fit.delta2) {
                        mat[surv.t>0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)):=
                                          -(1-(get(paste0("F", fit.delta, ".tau", kk)) - get(paste0("F", fit.delta, ".t"))) / surv.t)]
                        mat[surv.t==0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)):=-1]
                    } else {
                        mat[surv.t>0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)):=
                                          (get(paste0("F", fit.delta, ".tau", kk)) - get(paste0("F", fit.delta, ".t"))) / surv.t]
                        mat[surv.t==0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)):=0]
                    }
                }
            }
        }
        if (target.S) {
            for (each2 in outcome.index) {
                fit.delta2 <- estimation[[each2]][["event"]]
                for (kk in 1:length(tau)) {
                    mat[, (paste0("Ht", ".lambda", fit.delta2, ".", kk)):=0]   
                    for (each in outcome.index) {
                        fit.delta <- estimation[[each]][["event"]]
                        mat[, (paste0("Ht", ".lambda", fit.delta2, ".", kk)):=get(paste0("Ht", ".lambda", fit.delta2, ".", kk))-
                                  get(paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk))]
                    }
                }
            }
        }
    } else {
        for (kk in 1:length(tau)) {
            mat[surv.t>0, (paste0("Ht.lambda.", kk)):=get(paste0("surv.tau", kk)) / surv.t]
            mat[surv.t==0, (paste0("Ht.lambda.", kk)):=1]
        }
    }

    #-- 11 -- initial fit:
  
    if (cr) {
        if (treat.effect[1]=="stochastic") {
            init.fit <- lapply(target, function(each) {
                sapply(1:length(tau), function(kk) {
                    mean(rowSums(sapply(a, function(aa)
                    (mat[get(A.name)==aa, pi.star[1]*
                                          get(paste0("F", estimation[[outcome.index[each]]][["event"]],
                                                     ".tau", kk))[1],
                         by="id"][,2][[1]]))))
                })
            })
        } else {
            init.fit <- lapply(target, function(each) {
                sapply(1:length(tau), function(kk) {
                    mean(rowSums(sapply(a, function(aa)
                    (2*(aa==a[1])-1)*(mat[get(A.name)==aa, get(paste0("F", estimation[[outcome.index[each]]][["event"]],
                                                                      ".tau", kk))[1],
                                          by="id"][,2][[1]]))))
                })
            })
        }
        names(init.fit) <- paste0("F", sapply(outcome.index[target], function(each) estimation[[each]][["event"]]))
    } else {
        if (treat.effect[1]=="stochastic") {
            init.fit <- sapply(1:length(tau), function(kk) {
                mean(rowSums(sapply(a, function(aa)
                (mat[get(A.name)==aa, pi.star[1]*(1-get(paste0("surv.tau", kk))[1]), by="id"][,2][[1]]))))
            })
        } else {
            init.fit <- sapply(1:length(tau), function(kk) {
                mean(rowSums(sapply(a, function(aa)
                (2*(aa==a[1])-1)*(mat[get(A.name)==aa, 1-get(paste0("surv.tau", kk))[1], by="id"][,2][[1]]))))
            })
        }
    }
    
    if (cr) {
        eval.ic <- function(mat, fit, target.index=outcome.index, Sigma=FALSE, tau.values=tau,
                            survival=FALSE, tau.all=tau) {
            outer <- lapply(target.index, function(each) {
                fit.delta <- estimation[[each]][["event"]]
                each.index <- (1:length(target.index))[target.index==each]
                sapply(1:length(tau.values), function(kk) {
                    k2 <- (1:length(tau.all))[tau.all==max(tau.all[tau.all<=tau.values[kk]])]
                    out <- 0
                    for (each2 in outcome.index) {
                        fit.delta2 <- estimation[[each2]][["event"]]
                        mat[get(paste0("Ht", fit.delta, ".lambda", fit.delta2,".", k2))<=-500,
                        (paste0("Ht", fit.delta, ".lambda", fit.delta2,".", k2)):=-500]
                        out2 <- mat[, sum( (get(A.name)==A.obs) * (get(time.var)<=tau.values[kk]) *
                                           (get(time.var)<=time.obs) * Ht *
                                           get(paste0("Ht", fit.delta, ".lambda", fit.delta2,".", k2)) *
                                           ( (delta.obs==fit.delta2 & get(time.var)==time.obs) -
                                             get(paste0("dhaz", fit.delta2)) * get(paste0("fit.cox", fit.delta2)) ),
                                          na.rm=TRUE), by="id"]
                        out <- out + out2[, 2][[1]]
                    }
                    if (treat.effect[1]=="stochastic") {
                        ic.squared <- (out + rowSums(sapply(a, function(aa)
                        (mat[get(A.name)==aa, pi.star[1]*
                                              get(paste0("F", fit.delta,
                                                         ".tau", k2))[1],
                             by="id"][,2][[1]])))-
                        fit[[paste0("F", fit.delta)]][kk])
                    } else {
                        ic.squared <- (out + rowSums(sapply(a, function(aa)
                        (2*(aa==a[1])-1)*(mat[get(A.name)==aa, get(paste0("F", fit.delta,
                                                                          ".tau", k2))[1],
                                              by="id"][,2][[1]])))-
                        (length(a)==2)*fit[[paste0("F", fit.delta)]][kk]-
                        (length(a)==1)*fit[[paste0("F", fit.delta)]][kk])
                    }
                    if (survival) return(ic.squared)
                    if (Sigma) return(ic.squared) else return(sqrt(mean(ic.squared^2)/n))
                })
            })
            if (Sigma) {
                if (survival & length(outcome.index)==length(target)) {
                    outer2 <- Reduce("+", lapply(outer, function(xx) -xx))
                } else {
                    outer2 <- do.call("cbind", outer)
                }
                Sigma.list <- lapply(1:n, function(i) {
                    t(outer2[i,,drop=FALSE])%*%outer2[i,,drop=FALSE]
                })
                return(Reduce("+", Sigma.list) / length(Sigma.list))
            } else {
                names(outer) <- paste0("F", sapply(target.index, function(each) estimation[[each]]["event"]))
                return(outer)
            }
        }
    } else {
        eval.ic <- function(mat, fit, target.index=1, Sigma=FALSE, tau.values=tau, tau.all=tau) {
            outer <- sapply(1:length(tau.values), function(kk) {
                k2 <- (1:length(tau.all))[tau.all==max(tau.all[tau.all<=tau.values[kk]])]
                by.vars <- "id"
                out <- mat[, sum( (get(A.name)==A.obs) * (get(time.var)<=tau.values[kk]) *
                                  (get(time.var)<=time.obs) * Ht *
                                  get(paste0("Ht.lambda.", k2)) *
                                  ( (delta.obs==1 & get(time.var)==time.obs) -
                                    dhaz1 * fit.cox1 ), na.rm=TRUE), by=by.vars]
                if (treat.effect[1]=="stochastic") {
                    ic.squared <- (out[, 2][[1]] +
                                   rowSums(sapply(a, function(aa)
                                   (mat[get(A.name)==aa,
                                        pi.star[1]*(get(paste0("surv.tau", k2))[1]),
                                        by="id"][,2][[1]]))) -
                                   (1-fit[k2]))
                } else {
                    ic.squared <- (out[, 2][[1]] +
                                   rowSums(sapply(a, function(aa)
                                   (2*(aa==a[1])-1)*(mat[get(A.name)==aa,
                                                         get(paste0("surv.tau", k2))[1],
                                                         by="id"][,2][[1]]))) -
                                   (length(a)==2)*fit[kk] -
                                   (length(a)==1)*(1-fit[kk]))
                }
                if (Sigma) return(ic.squared) else return(sqrt(mean(ic.squared^2)/n))
            })
            if (Sigma) {
                Sigma.list <- lapply(1:nrow(outer), function(i) {
                    t(outer[i,,drop=FALSE])%*%outer[i,,drop=FALSE]
                })
                return(Reduce("+", Sigma.list) / nrow(outer))
            } else return(outer)
        }
    }

    init.ic <- eval.ic(mat, fit=init.fit, target.index=outcome.index[target])
    
    if (weighted.norm[1]=="Sigma") {
        if (target.S) {
            Sigma <- eval.ic(mat, init.fit, target.index=outcome.index[target], Sigma=TRUE, survival=TRUE)
        } else {
            Sigma <- eval.ic(mat, init.fit, target.index=outcome.index[target], Sigma=TRUE)
        }
        Sigma.inv <-  try(solve(Sigma))
        if (any(class(Sigma.inv)=="try-error")) {
            Sigma.inv <- solve(Sigma+diag(x=0.000001, nrow=nrow(Sigma)))
            warning("regularization of Sigma needed for inversion")
            #weighted.norm[1] <- "sigma"
            not.fit.list[[length(not.fit.list)+1]] <-
                paste0("Sigma was regularized")
        }
    }

    if (FALSE) {
        heatmap(Sigma)
        heatmap(Sigma.inv)
    }

    if (cr) {
        init.list <-  lapply(1:length(init.fit), function(each.index) {
            out <- rbind(init.est=init.fit[[each.index]],
                         init.se=init.ic[[each.index]])
            colnames(out) <- paste0("tau=", tau)
            return(out)
        })
        names(init.list) <- paste0("F", sapply(outcome.index[target], function(each) estimation[[each]]["event"]))
        if (length(target)==length(outcome.index)) {
            S.se.init <- sapply(tau, function(tt) sqrt(mean(rowSums(sapply(1:length(outcome.index), function(target11) {
                unlist(eval.ic(mat, unlist(lapply(init.list, function(xx) xx["init.est", paste0("tau=", tt)])),
                               target.index=outcome.index[target11],
                               tau.values=tt, survival=TRUE))
            }))^2)/n))
            S.fit.init <- sapply(tau, function(tt) (treat.effect!="ate")-sum(sapply(init.list, function(fl) fl["init.est",paste0("tau=", tt)])))
            init.list$S <- rbind(init.est=S.fit.init, init.se=S.se.init)
            colnames(init.list$S) <- paste0("tau=", tau)
        }
        if (target.S) {
            tmle.list <- list(init=init.list[names(init.list)=="S"])
        } else {
            tmle.list <- list(init=init.list)
        }
        if (output.km) {
            km.list <- lapply(1:min(length(target),length(tmle.list$init)), function(each.index) {
                out <- rbind(km.est=km.est[((each.index-1)*length(tau)+1):(each.index*length(tau))],
                             km.se=km.se[((each.index-1)*length(tau)+1):(each.index*length(tau))])
                colnames(out) <- paste0("tau=", tau)
                return(out)
            })
            if (target.S) {
                names(km.list) <- "S"
            } else {
                names(km.list) <- paste0("F", sapply(outcome.index[target], function(each) estimation[[each]]["event"]))
            }
            tmle.list$km <- km.list
        }
    }  else {
        init.list <- rbind(init.est=init.fit, init.se=init.ic)
        colnames(init.list) <- paste0("tau=", tau)
        tmle.list <- list(init=init.list)
        if (output.km) {
            km.list <- rbind(km.est=km.est, km.se=km.se)
            colnames(km.list) <- paste0("tau=", tau)
            tmle.list$km <- km.list
        }
    }

    if (output.mat=="init") {
        return(mat)
    }

    #---------------------------------------------------------------
    #-- 12 -- TMLE:

    if (cr) {
        eval.equation <- function(mat, eps=0, target.index=outcome.index, cr.index=outcome.index,
                                  tau.values=tau, tau.all=tau) {
            outer <- lapply(target.index, function(each) {
                fit.delta <- estimation[[each]][["event"]]
                sapply(1:length(tau.values), function(kk) {
                    k2 <- (1:length(tau.all))[tau.all==max(tau.all[tau.all<=tau.values[kk]])]
                    out <- 0
                    for (each2 in cr.index) {
                        fit.delta2 <- estimation[[each2]][["event"]]
                        out2 <- mat[(get(time.var)<=tau.values[kk]) & (get(time.var)<=time.obs),
                                    sum( (get(A.name)==A.obs) *
                                         (get(time.var)<=time.obs) * Ht *
                                         get(paste0("Ht", fit.delta, ".lambda", fit.delta2,".", k2)) *
                                         ( (delta.obs==fit.delta2 & get(time.var)==time.obs) -
                                           exp( eps * Ht *
                                                get(paste0("Ht", fit.delta, ".lambda", fit.delta2,".", k2)) ) *
                                           get(paste0("dhaz", fit.delta2)) * get(paste0("fit.cox", fit.delta2)) )), by="id"]
                        out <- out + out2[, 2][[1]]
                    }
                    return(mean(out))
                })
            })
            names(outer) <- paste0("F", sapply(target.index, function(each) estimation[[each]]["event"]))
            return(outer)
        }
    } else {
        eval.equation <- function(mat, eps=0, target.index=1, cr.index=1, tau.values=tau, tau.all=tau) {
            sapply(1:length(tau.values), function(kk) {
                k2 <- (1:length(tau.all))[tau.all==max(tau.all[tau.all<=tau.values[kk]])]
                out <- mat[(get(time.var)<=tau.values[kk]), sum( (get(A.name)==A.obs) *
                                                                 (get(time.var)<=time.obs) * Ht *
                                                                 get(paste0("Ht.lambda.", k2)) *
                                                                 ( (delta.obs==1 & get(time.var)==time.obs) -
                                                                   exp( eps * Ht *
                                                                        get(paste0("Ht.lambda.", k2)) ) *
                                                                   dhaz1 * fit.cox1 )), by="id"]
                return(mean(out[, 2][[1]]))
            })
        }
    }

    if (check.sup & weighted.norm[1]=="Sigma" & FALSE) {
        print("-----")
        print(c(matrix(eval.equation(mat, eps=0), nrow=1)%*%Sigma.inv%*%matrix(eval.equation(mat, eps=0), ncol=1)))
        print(c(matrix(eval.equation(mat, eps=0), nrow=1)%*%diag(1/(init.ic^2*n))%*%matrix(eval.equation(mat, eps=0), ncol=1)))
        print(nrow(Sigma.inv)*c(matrix(eval.equation(mat, eps=0), nrow=1)%*%Sigma.inv%*%matrix(eval.equation(mat, eps=0), ncol=1))>=
              c(matrix(eval.equation(mat, eps=0), nrow=1)%*%diag(1/(init.ic^2*n))%*%matrix(eval.equation(mat, eps=0), ncol=1)))
        print("-----")
    }

    #----------------------------------------
    #-- 12 (I) -- multivariate one-step tmle:

    if ((length(tau)>1 & !iterative) | one.step | (length(target)>1 & !iterative)) {

        second.round <- FALSE

        if (cr) {
            Pn.eic.fun <- function(mat) {
                eval <- eval.equation(mat, eps=0, target.index=outcome.index[target])
                if (target.S) { #-- if only want to target survival in competing risks setting
                    eval <- sapply(tau, function(tt) -sum(sapply(eval, function(xx) xx[tau==tt])))
                }
                return(eval)
            }
            if (separate.cr) {
                Pn.eic.fun.separate <- function(mat) {
                    eval <- lapply(outcome.index[target], function(jtarget) {
                        sapply(outcome.index, function(jhazard) {
                            eval.equation(mat, eps=0, target.index=jtarget, cr.index=jhazard)
                        })
                    })
                    return(eval)
                }
            }
        } else {
            Pn.eic.fun <- function(mat) {
                eval <- eval.equation(mat, eps=0)
                return(eval)
            }
        }
        
        Pn.eic <- Pn.eic.fun(mat)

        Pn.eic.norm.fun <- function(x2, x) {
            return(sqrt(sum(unlist(x2)*unlist(x))))
        }
 
        #-- which type of weighted norm are we using? 
        if (weighted.norm[1]=="sigma") {
            if (target.S) { #-- if only want to target survival in competing risks setting
                Pn.eic2.fun <- function(Pn.eic) Pn.eic/(ifelse(any(unlist(S.se.init)==0), S.se.init+0.00001, S.se.init)^2*n)
            } else {
                Pn.eic2.fun <- function(Pn.eic) {
                    if (is.list(Pn.eic)) {
                        lapply(1:length(Pn.eic), function(kk) Pn.eic[[kk]]/(ifelse(any(unlist(init.ic)==0), init.ic[[kk]]+0.00001, init.ic[[kk]])^2*n))
                    } else {
                        Pn.eic/(ifelse(any(unlist(init.ic)==0), init.ic+0.00001, init.ic)^2*n)
                    }
                }
            }
        } else if (weighted.norm[1]=="Sigma") {
            if (target.S) {
                Pn.eic2.fun <- function(Pn.eic) as.numeric(matrix(unlist(Pn.eic), 1, length(unlist(Pn.eic)))%*%Sigma.inv)
            } else {
                Pn.eic2.fun <- function(Pn.eic) {
                    Pn.eic.tmp <- as.numeric(matrix(unlist(Pn.eic), 1, length(unlist(Pn.eic)))%*%Sigma.inv)
                    if (length(target)>1) {
                        return(lapply(1:length(target), function(kk) Pn.eic.tmp[(1+(kk-1)*length(tau)):(kk*length(tau))]))
                    } else {
                        return(Pn.eic.tmp)
                    }
                }
            }
        } else {
            Pn.eic2.fun <- function(Pn.eic) Pn.eic
        }

        Pn.eic2 <- Pn.eic2.fun(Pn.eic)
      
        criterion <- 1/(sqrt(n)*log(n))
        
        Pn.eic.norm.prev <- Pn.eic.norm <- Pn.eic.norm.fun(Pn.eic2, Pn.eic)

        if (verbose) print(paste0("Pn.eic.norm=", Pn.eic.norm))

        for (step in 1:no.small.steps) {

            if (cr) {
                for (each in outcome.index) {
                    fit.delta <- estimation[[each]][["event"]]
                    mat[, (paste0("fit.cox", fit.delta, ".tmp")):=get(paste0("fit.cox", fit.delta))]
                    for (kk in 1:length(tau)) {
                        if (target.S) mat[, (paste0("Ht", ".lambda", fit.delta,".", kk, ".tmp")):=
                                                get(paste0("Ht", ".lambda", fit.delta,".", kk))]
                        for (each2 in outcome.index[target]) {
                            fit.delta2 <- estimation[[each2]][["event"]]
                            mat[, (paste0("Ht", fit.delta2, ".lambda", fit.delta,".", kk, ".tmp")):=
                                      get(paste0("Ht", fit.delta2, ".lambda", fit.delta,".", kk))]
                        }
                    }
                }
            } else {
                mat[, fit.cox1.tmp:=fit.cox1]
                for (kk in 1:length(tau)) {
                    mat[, (paste0("Ht.lambda.", kk, ".tmp")):=get((paste0("Ht.lambda.", kk)))]
                }
            }

            if (cr) {
                if (target.S) {
                    for (each in outcome.index) {
                        fit.delta <- estimation[[each]][["event"]]
                        mat[, (paste0("delta", fit.delta, ".dx")):=0]
                        for (kk in 1:length(tau)) {
                            mat[, (paste0("delta", fit.delta, ".dx")):=
                                      get(paste0("delta", fit.delta, ".dx"))+
                                      (get(time.var)<=tau[kk])*Ht*(
                                          (get(paste0("Ht",".lambda", fit.delta,".", kk)))*
                                          Pn.eic2[kk]
                                      )/Pn.eic.norm]
                        }
                    }
                } else {
                    if (!separate.cr) {
                        for (each in outcome.index) {
                            fit.delta <- estimation[[each]][["event"]]
                            mat[, (paste0("delta", fit.delta, ".dx")):=0]
                            for (each2 in outcome.index[target]) {
                                fit.delta2 <- estimation[[each2]][["event"]]
                                for (kk in 1:length(tau)) {
                                    mat[, (paste0("delta", fit.delta, ".dx")):=
                                              get(paste0("delta", fit.delta, ".dx"))+
                                              (get(time.var)<=tau[kk])*Ht*(
                                                  (get(paste0("Ht", fit.delta2,".lambda", fit.delta,".", kk)))*
                                                  Pn.eic2[each2==outcome.index[target]][[1]][kk]
                                              )/Pn.eic.norm]
                                }
                            }
                        }
                    } else {
                        Pn.eic.separate <- Pn.eic.fun.separate(mat)
                        for (each in outcome.index) {
                            fit.delta <- estimation[[each]][["event"]]
                            mat[, (paste0("delta", fit.delta, ".dx")):=0]
                            for (each2 in outcome.index[target]) {
                                fit.delta2 <- estimation[[each2]][["event"]]
                                for (kk in 1:length(tau)) {
                                    mat[, (paste0("delta", fit.delta, ".dx")):=
                                              get(paste0("delta", fit.delta, ".dx"))+
                                              (get(time.var)<=tau[kk])*Ht*(
                                                  (get(paste0("Ht", fit.delta2,".lambda", fit.delta,".", kk)))*
                                                  Pn.eic.separate[each2==outcome.index[target]][[1]][each==outcome.index[target]][[1]][kk]
                                              )/Pn.eic.norm]
                                }
                            }
                        }
                    }
                }
            } else {
                mat[, delta.dx:=0]
                for (kk in 1:length(tau)) {
                    mat[, delta.dx:=delta.dx+
                              (get(time.var)<=tau[kk])*
                              Ht*get(paste0("Ht.lambda.", kk))*Pn.eic2[kk]/Pn.eic.norm]
                }
            }

            deps <- deps.size

            if (cr) {
                mat[, surv.t:=1]
                for (each in outcome.index) {
                    fit.delta <- estimation[[each]][["event"]]
                    mat[, (paste0("fit.cox", fit.delta)):=
                              get(paste0("fit.cox", fit.delta))*
                              exp(deps*get(paste0("delta", fit.delta, ".dx")))]
                    mat[get(paste0("fit.cox", fit.delta))>500, (paste0("fit.cox", fit.delta)):=500]
                    mat[, surv.t:=surv.t*exp(-cumsum(get(paste0("dhaz", fit.delta))*
                                                     get(paste0("fit.cox", fit.delta)))),
                        by=c("id", A.name)]
                    mat[get(paste0("fit.cox", fit.delta))==Inf, surv.t:=0]
                }
                mat[, surv.t1:=c(0, surv.t[-.N]), by=c("id", A.name)]
                for (kk in 1:length(tau)) {
                    mat[, (paste0("surv.tau", kk)):=
                              surv.t[get(time.var)==max(get(time.var)[get(time.var)<=tau[kk]])],
                        by=c("id", A.name)]
                }
                for (each in outcome.index[target]) {
                    fit.delta <- estimation[[each]][["event"]]
                    mat[, (paste0("F", fit.delta, ".t")):=cumsum(surv.t*get(paste0("dhaz", fit.delta))*
                                                                 get(paste0("fit.cox", fit.delta))),
                        by=c("id", A.name)]
                    for (kk in 1:length(tau)) {
                        mat[, (paste0("F", fit.delta, ".tau", kk)):=
                                  get(paste0("F", fit.delta, ".t"))[get(time.var)==max(get(time.var)[get(time.var)<=tau[kk]])],
                            by=c("id", A.name)]
                        for (each2 in outcome.index) {
                            fit.delta2 <- estimation[[each2]][["event"]]
                            if (fit.delta==fit.delta2) {
                                mat[surv.t>0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)):=
                                                  -(1-(get(paste0("F", fit.delta, ".tau", kk)) - get(paste0("F", fit.delta, ".t"))) / surv.t)]
                                mat[round(surv.t,8)==0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)):=-1]
                            } else {
                                mat[surv.t>0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)):=
                                                  (get(paste0("F", fit.delta, ".tau", kk)) - get(paste0("F", fit.delta, ".t"))) / surv.t]
                                mat[round(surv.t,8)==0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)):=0]
                            }
                        }
                    }
                }
                if (target.S) {
                    for (each2 in outcome.index) {
                        fit.delta2 <- estimation[[each2]][["event"]]
                        for (kk in 1:length(tau)) {
                            mat[, (paste0("Ht", ".lambda", fit.delta2, ".", kk)):=0]   
                            for (each in outcome.index) {
                                fit.delta <- estimation[[each]][["event"]]
                                mat[, (paste0("Ht", ".lambda", fit.delta2, ".", kk)):=get(paste0("Ht", ".lambda", fit.delta2, ".", kk))-
                                          get(paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk))]
                            }
                        }
                    }
                }
            } else {
                mat[, fit.cox1:=fit.cox1*exp(deps*delta.dx)]
                mat[fit.cox1>500, fit.cox1:=500]
                mat[, surv.t:=exp(-cumsum(dhaz1*fit.cox1)), by=c("id", A.name)]
                for (kk in 1:length(tau)) {
                    mat[, (paste0("surv.tau", kk)):=
                              surv.t[get(time.var)==max(get(time.var)[get(time.var)<=tau[kk]])],
                        by=c("id", A.name)]
                    mat[surv.t>0, (paste0("Ht.lambda.", kk)):=get(paste0("surv.tau", kk)) / surv.t]
                    mat[surv.t==0, (paste0("Ht.lambda.", kk)):=1]
                }
            }

            Pn.eic <- Pn.eic.fun(mat)

            Pn.eic2 <- Pn.eic2.fun(Pn.eic)

            #print("CHECK 3")
            #print(Pn.eic2)
            
            Pn.eic.norm <- Pn.eic.norm.fun(Pn.eic2, Pn.eic)

            print(paste0("step = ", step))
            if (verbose) print(paste0("Pn.eic.norm=", Pn.eic.norm))
            
            if (Pn.eic.norm.prev<=Pn.eic.norm) { # reset all columns to '.tmp'
                
                if (verbose) {
                    print("----")
                    print(step)
                    print("----")
                }

                if (cr) {
                    for (each in outcome.index) {
                        fit.delta <- estimation[[each]][["event"]]
                        mat[, (paste0("fit.cox", fit.delta)):=get(paste0("fit.cox", fit.delta, ".tmp"))]
                        for (kk in 1:length(tau)) {
                            if (target.S) mat[, (paste0("Ht", ".lambda", fit.delta,".", kk)):=
                                                    get(paste0("Ht", ".lambda", fit.delta,".", kk, ".tmp"))]
                            for (each2 in outcome.index[target]) {
                                fit.delta2 <- estimation[[each2]][["event"]]
                                mat[, (paste0("Ht", fit.delta2, ".lambda", fit.delta,".", kk)):=
                                          get(paste0("Ht", fit.delta2, ".lambda", fit.delta,".", kk, ".tmp"))]
                            }
                        }
                    }
                } else {
                    mat[, fit.cox1.tmp:=fit.cox1]
                    for (kk in 1:length(tau)) {
                        mat[, (paste0("Ht.lambda.", kk)):=get((paste0("Ht.lambda.", kk, ".tmp")))]
                    }
                }

                Pn.eic <- Pn.eic.fun(mat)
                Pn.eic2 <- Pn.eic2.fun(Pn.eic)
                
                Pn.eic.norm <- Pn.eic.norm.fun(Pn.eic2, Pn.eic)
                deps.size <- 0.5*deps.size#0.1*deps.size

            } else {
                
                Pn.eic.norm.prev <- Pn.eic.norm

                if (target.S) {
                    Pn.eic3 <- sapply(1:length(Pn.eic), function(kk) Pn.eic[[kk]]/(ifelse(any(unlist(S.se.init)==0), S.se.init[kk]+0.001, S.se.init[kk])*sqrt(n)))
                } else {
                    Pn.eic3 <- lapply(1:length(Pn.eic), function(kk) Pn.eic[[kk]]/(ifelse(any(unlist(init.ic)==0), init.ic[[kk]]+0.001, init.ic[[kk]])*sqrt(n)))
                }
            }

            if (target.S) {
                Pn.eic3 <- sapply(1:length(Pn.eic), function(kk) Pn.eic[[kk]]/(ifelse(any(unlist(S.se.init)==0), S.se.init[kk]+0.001, S.se.init[kk])*sqrt(n)))
            } else {
                Pn.eic3 <- lapply(1:length(Pn.eic), function(kk) Pn.eic[[kk]]/(ifelse(any(unlist(init.ic)==0), init.ic[[kk]]+0.001, init.ic[[kk]])*sqrt(n)))
            }
            
            if (cr & length(target)==length(outcome.index) & !target.S) {
                #--- here in fact want to check that we solve survival eic well enough!
                #--- i.e., we add it to our vector of Pn.eic3; 
                Pn.eic3[[length(Pn.eic3)+1]] <- sapply(1:length(Pn.eic3[[1]]), function(jj) {
                    sum(sapply(Pn.eic, function(xx) xx[[jj]]))
                }) / (S.se.init*sqrt(n))
                check.sup.norm <- max(abs(unlist(Pn.eic3)))<=(criterion)
            } else {
                check.sup.norm <- max(abs(unlist(Pn.eic3)))<=(criterion)
            }

            if (verbose) print(paste0("target level = ", criterion))
            if (verbose) print(paste0("max eic = ", max(abs(unlist(Pn.eic3)))))

            if (no.update.if.solved) {
                
                if (length(target)>1) {
                    Pn.eic2 <- lapply(1:length(Pn.eic2), function(kk) {
                        tmp <- Pn.eic2[[kk]]
                        tmp[abs(Pn.eic3[[kk]])<=criterion] <- 0
                        return(tmp)
                    })
                } else {
                    Pn.eic2[abs(Pn.eic3)<=criterion] <- 0
                }

                if (all(unlist(Pn.eic2)==0)) check.sup.norm <- TRUE else Pn.eic.norm <- Pn.eic.norm.fun(Pn.eic2, Pn.eic)
                print(Pn.eic2)

            }
            
            if (check.sup.norm | step==no.small.steps) {#(left.criterion<=criterion) {
                if (step==no.small.steps) {
                    message("Warning: Algorithm did not converge")
                }
                
                if (verbose) print(paste0("converged", " at ", step, "th step"))
                if (verbose) print(paste0("eic = ", Pn.eic.fun(mat)))

                if (output.mat %in% paste0(1:10)) {
                    return(mat)
                }

                #-- 12c -- compute sd:
                if (cr) {
                    if (treat.effect[1]=="stochastic") {
                        final.fit <- lapply(target, function(each) {
                            sapply(1:length(tau), function(kk) {
                                mean(rowSums(sapply(a, function(aa)
                                (mat[get(A.name)==aa, pi.star[1]*
                                                      get(paste0("F", estimation[[outcome.index[each]]][["event"]],
                                                                 ".tau", kk))[1],
                                     by="id"][,2][[1]]))))
                            })
                        })
                    } else {
                        final.fit <- lapply(outcome.index[target], function(each) {
                            sapply(1:length(tau), function(kk) {
                                mean(rowSums(sapply(a, function(aa)
                                (2*(aa==a[1])-1)*(mat[get(A.name)==aa,
                                                      get(paste0("F", estimation[[each]][["event"]],
                                                                 ".tau", kk))[1],
                                                      by="id"][,2][[1]]))))
                            })
                        })
                    }
                    names(final.fit) <- paste0("F", sapply(outcome.index[target], function(each) estimation[[each]][["event"]]))
                    final.ic <- eval.ic(mat, final.fit, target.index=outcome.index[target])
                } else {
                    if (treat.effect[1]=="stochastic") {
                        final.fit <- list(sapply(1:length(tau), function(kk) {
                            mean(rowSums(sapply(a, function(aa)
                            (mat[get(A.name)==aa, pi.star[1]*(1-get(paste0("surv.tau", kk))[1]), by="id"][,2][[1]]))))
                        }))
                    } else {
                        final.fit <- list(sapply(1:length(tau), function(kk) {
                            mean(rowSums(sapply(a, function(aa)
                            (2*(aa==a[1])-1)*(mat[get(A.name)==aa, 1-get(paste0("surv.tau", kk))[1], by="id"][,2][[1]]))))
                        }))
                    }
                    final.ic <- list(eval.ic(mat, final.fit[[1]]))
                }
        
                final.list <- lapply(1:length(final.fit), function(each.index) {
                    out <- rbind(tmle.est=final.fit[[each.index]],
                                 tmle.se=final.ic[[each.index]])
                    colnames(out) <- paste0("tau=", tau)
                    return(out)
                })

                if (cr) names(final.list) <- paste0("F", sapply(outcome.index[target], function(each) estimation[[each]]["event"])) else final.list <- final.list[[1]]

                if (length(target)==length(outcome.index) & length(outcome.index)>1) {
                    S.se <- sapply(tau, function(tt) sqrt(mean(rowSums(sapply(1:length(outcome.index), function(target11) {
                        unlist(eval.ic(mat, unlist(lapply(final.list, function(xx) xx["tmle.est", paste0("tau=", tt)])),
                                       target.index=outcome.index[target11],
                                       tau.values=tt, survival=TRUE))
                    }))^2)/n))
                    S.fit <- sapply(tau, function(tt) (treat.effect!="ate")-sum(sapply(final.list, function(fl) fl["tmle.est",paste0("tau=", tt)])))
                    final.S <- rbind(S.fit, S.se)
                    colnames(final.S) <- paste0("tau=", tau)
                    rownames(final.S) <- c("tmle.est", "tmle.se")
                    final.list$S <- final.S
                    if (target.S) {
                        tmle.list$tmle <- final.list[names(final.list)=="S"]
                    } else {
                        tmle.list$tmle <- final.list
                    }
                } else {
                    tmle.list$tmle <- final.list
                }

                if (!second.round) {
                    if (check.sup) tmle.list$check.sup.norm <- list(
                                       check.sup.norm=check.sup.norm,
                                       lhs=max(abs(unlist(Pn.eic3))),
                                       rhs=criterion/sqrt(length(target)*length(tau)))
                    tmle.list$convergenced.at.step <- step
                   
                    if (!push.criterion)
                        break else criterion <- criterion/sqrt(n)
                    second.round <- TRUE
                } else {
                    if (target.S) {
                        tmle.list$tmle.second.round <- final.list[names(final.list)=="S"]
                    } else {
                        tmle.list$tmle.second.round <- final.list
                    }
                    if (check.sup) tmle.list$check.sup.norm.second.round <- list(
                                       check.sup.norm=
                                           max(abs(unlist(Pn.eic3)))<=(criterion),
                                       lhs=max(abs(unlist(Pn.eic3))),
                                       rhs=criterion)
                    tmle.list$convergenced.at.step.second.round <- step
                    break
                }               
            }
        }
        
    } else {

        #---------------------------------------------------------------
        # ITERATIVE TMLE

        mat1 <- copy(mat)

        if (length(target)==length(outcome.index) & length(outcome.index)>1)
            if (!target.S) target2 <- c(target, 0) else target2 <- 0 else target2 <- target

        update.fit <- lapply(target2, function(target1) {
            
            update.fit.inner <- lapply(tau, function(tt) {

                mat <- copy(mat1)
               # if (target1==0) browser()
                
                if (target1==0) {
                    target1 <- 1:length(outcome.index)
                }

                target1.index <- (1:length(target))[target==target1]

                #; target1 <- outcome.index

                for (iter in 1:maxIter) {

                    if (FALSE) {
                        
                        plot(sapply(seq(-1, 1, length=100), function(eps) unlist(eval.equation(mat, eps,
                                                                                               target.index=outcome.index[target1],
                                                                                               cr.index=index))))
                        lines(sapply(seq(-1, 1, length=100), function(eps) sum(unlist(
                                                                               eval.equation(mat, eps, target.index=outcome.index[target1]))[tau==tt])))
                    }

                    if (length(target1)==1) {
                        if (iter==1 & verbose) print(paste0("F", target1, ":"))
                        eps.hat <- sapply(outcome.index, function(index) { 
                            nleqslv(0.01, function(eps) unlist(eval.equation(mat, eps,
                                                                             target.index=outcome.index[target1],
                                                                             cr.index=index))[tau==tt])$x
                        })
                    } else {
                        if (iter==1 & verbose) print(paste0("S", ":"))
                        eps.hat <- sapply(outcome.index, function(index) {
                            nleqslv(0.01, function(eps)
                                sum(sapply(1:length(outcome.index), function(target1)
                                    unlist(eval.equation(mat, eps,
                                                         target.index=outcome.index[target1],
                                                         cr.index=index))[tau==tt])))$x
                        })
                    }
            
                    if (verbose) print(sapply(1:length(eps.hat),
                                              function(index) {
                                                  paste0("iter=", iter, ", eps",
                                                         ifelse(cr, paste0("(cause=",
                                                                           estimation[[outcome.index[index]]][["event"]],")"), ""),
                                                         "=", round(eps.hat[index], 4))}))

                    #-- 12b -- update hazard(s):
                    if (cr) {
                        mat[, surv.t:=1]
                        #for (target11 in target1) {
                        for (kk in (1:length(tau))[tau==tt]) {
                            for (each in outcome.index) {
                                fit.delta <- estimation[[each]][["event"]]
                                index <- (1:length(outcome.index))[outcome.index==each]
                                mat[, (paste0("fit.cox", fit.delta)):=
                                          get(paste0("fit.cox", fit.delta))*
                                          exp(eps.hat[index]*Ht*
                                              rowSums(sapply(target1, function(target11) get(paste0("Ht", target11,".lambda", fit.delta,".", kk)))))]
                                mat[get(paste0("fit.cox", fit.delta))>500, (paste0("fit.cox", fit.delta)):=500]
                                mat[, surv.t:=surv.t*exp(-cumsum(get(paste0("dhaz", fit.delta))*
                                                                 get(paste0("fit.cox", fit.delta)))),
                                    by=c("id", A.name)]
                                mat[get(paste0("fit.cox", fit.delta))==Inf, surv.t:=0]
                            }
                        }
                        #}
                        mat[, surv.t1:=c(0, surv.t[-.N]), by=c("id", A.name)]
                        for (each in outcome.index[target1]) {
                            fit.delta <- estimation[[each]][["event"]]
                            mat[, (paste0("F", fit.delta, ".t")):=cumsum(surv.t*get(paste0("dhaz", fit.delta))*
                                                                         get(paste0("fit.cox", fit.delta))),
                                by=c("id", A.name)]
                            for (kk in (1:length(tau))[tau==tt]) {
                                mat[, (paste0("F", fit.delta, ".tau", kk)):=
                                          get(paste0("F", fit.delta, ".t"))[get(time.var)==max(get(time.var)[get(time.var)<=tau[kk]])],
                                    by=c("id", A.name)]
                                for (each2 in outcome.index) {
                                    fit.delta2 <- estimation[[each2]][["event"]]
                                    if (fit.delta==fit.delta2) {
                                        mat[surv.t>0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)):=
                                                          -(1-(get(paste0("F", fit.delta, ".tau", kk)) - get(paste0("F", fit.delta, ".t"))) / surv.t)]
                                        mat[surv.t==0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)):=-1]
                                    } else {
                                        mat[surv.t>0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)):=
                                                          (get(paste0("F", fit.delta, ".tau", kk)) - get(paste0("F", fit.delta, ".t"))) / surv.t]
                                        mat[surv.t==0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)):=0]
                                    }
                                }
                            }
                        }
                    } else {
                        mat[, surv.t:=1]
                        for (kk in (1:length(tau))[tau==tt]) {
                            mat[, (paste0("fit.cox1")):=
                                      get(paste0("fit.cox1"))*
                                      exp(eps.hat*Ht*
                                          get(paste0("Ht.lambda.", kk)))]
                        }
                        mat[, surv.t:=surv.t*exp(-cumsum(get(paste0("dhaz1"))*
                                                         get(paste0("fit.cox1")))),
                            by=c("id", A.name)]
                        mat[get(paste0("fit.cox", target1))==Inf, surv.t:=0]
                        mat[, surv.t1:=c(0, surv.t[-.N]), by=c("id", A.name)]
                        for (kk in (1:length(tau))[tau==tt]) {
                            mat[, (paste0("surv.tau", kk)):=
                                      surv.t[get(time.var)==max(get(time.var)[get(time.var)<=tau[kk]])],
                                by=c("id", A.name)]
                            mat[surv.t>0, (paste0("Ht.lambda.", kk)):=get(paste0("surv.tau", kk)) / surv.t]
                            mat[surv.t==0, (paste0("Ht.lambda.", kk)):=1]
                        }
                    }

                    if (output.mat %in% paste0(1:10) & iter==as.numeric(output.mat)) {
                        return(mat)
                    }

                    if (length(target1)==1) {
                        eval.iter <- abs(sum(unlist(
                            eval.equation(mat, 0, target.index=outcome.index[target1]))[tau==tt]))
                        if (is.list(init.ic)) rhs <- (init.ic[[target1.index]][tau==tt]*
                                                      sqrt(n))/(sqrt(n)*log(n)) else rhs <- (init.ic[tau==tt]#[[target1.index]]
                                                          *sqrt(n))/(sqrt(n)*log(n))                        
                    } else {
                        eval.iter <- abs(sum(
                            sapply(1:length(outcome.index), function(target1) {
                                unlist(
                                    eval.equation(mat, 0, target.index=outcome.index[target1]))[tau==tt]})))
                        if (is.list(init.ic)) rhs <- sqrt(mean(rowSums(sapply(target1, function(target11) {
                            unlist(eval.ic(mat, init.fit, target.index=outcome.index[target11],
                                           tau.values=tt, survival=TRUE))
                        }))^2))/(sqrt(n)*log(n)) else rhs <- (init.ic[tau==tt]#[[target1.index]]
                            *sqrt(n))/(sqrt(n)*log(n))
                    }

                    if (push.criterion) rhs <- rhs/sqrt(n)

                    if (eval.iter<=rhs | iter==maxIter) {
                        if (iter==maxIter) {
                            message("Warning: Algorithm did not converge")
                        }

                        #eval.equation(mat, eps=0, target.index=outcome.index[target])

                        #if (length(target1)>1) browser()

                        eval.iter <- abs(sum(
                            sapply(1:length(outcome.index), function(target1) {
                                unlist(
                                    eval.equation(mat, 0, target.index=outcome.index[target1]))[tau==tt]})))
                        
                        if (cr) {
                            if (treat.effect[1]=="stochastic") {
                                update.est <- sapply((1:length(tau))[tau==tt], function(kk) {
                                    mean(rowSums(sapply(a, function(aa)
                                    (mat[get(A.name)==aa,
                                         pi.star[1]*
                                         get(paste0("F", estimation[[outcome.index[target1]]][["event"]],
                                                    ".tau", kk))[1],
                                         by="id"][,2][[1]]))))
                                })
                            } else {
                                update.est <- c(sapply((1:length(tau))[tau==tt], function(kk) {
                                    lapply(target1, function(target11) {
                                        mean(rowSums(sapply(a, function(aa)
                                        (2*(aa==a[1])-1)*(mat[get(A.name)==aa,
                                                              get(paste0("F", estimation[[outcome.index[target11]]][["event"]],
                                                                         ".tau", kk))[1],
                                                              by="id"][,2][[1]]))))
                                    })
                                }))
                            }
                        } else {
                            if (treat.effect[1]=="stochastic") {
                                update.est <- sapply((1:length(tau))[tau==tt], function(kk) {
                                    mean(rowSums(sapply(a, function(aa)
                                    (mat[get(A.name)==aa, pi.star[1]*(1-get(paste0("surv.tau", kk))[1]), by="id"][,2][[1]]))))
                                })
                            } else {
                                update.est <- sapply((1:length(tau))[tau==tt], function(kk) {
                                    mean(rowSums(sapply(a, function(aa)
                                    (2*(aa==a[1])-1)*(mat[get(A.name)==aa, 1-get(paste0("surv.tau", kk))[1], by="id"][,2][[1]]))))
                                })
                            }
                        }
                        #update.est <- list(update.est)
                        names(update.est) <- paste0("F", target1)
                        if (length(target1)>1 | target.S) {
                            # if (length(target1)>1) browser()
                            update.ic <- sqrt(mean(rowSums(sapply(target1, function(target11) {
                                unlist(eval.ic(mat, update.est, target.index=outcome.index[target11],
                                               tau.values=tt, survival=TRUE))
                            }))^2)/n)#[[1]]
                            if (treat.effect[1]=="ate") {
                                update.est <- list(S=-sum(unlist(update.est)))
                            } else {
                                update.est <- list(S=1-sum(unlist(update.est)))
                            }
                        } else {
                            update.ic <- unlist(eval.ic(mat, update.est, target.index=outcome.index[target1],
                                                        tau.values=tt))#[[1]]
                        }
                        out.list <- list(update.est=update.est, update.ic=update.ic)
                        if (cr) {
                            if (length(target1)==1) {
                                names(out.list)[length(out.list)] <- paste0("F", target1)
                            } else {
                                names(out.list)[length(out.list)] <- paste0("S")
                            }
                        }
                        #break
    
                        return(list(update.est=update.est, update.ic=update.ic))
                    } #else if (iter==maxIter) {
                    #  message("Warning: Algorithm did not converge")
                    #}
                }
            })

            #-- 12c -- evaluate target parameter:

            if (output.mat %in% paste0(1:10)) {
                return(update.fit.inner)
            }

            update.list.inner <- do.call("cbind", lapply(1:length(update.fit.inner), function(each.index) {
                out <- rbind(tmle.est=update.fit.inner[[each.index]]$update.est,
                             tmle.se=update.fit.inner[[each.index]]$update.ic)
                colnames(out) <- paste0("tau=", tau[each.index])
                return(out)
            }))

            return(update.list.inner)
        })

        if (output.mat %in% paste0(1:10)) {
            return(update.fit)
        }
        
        #-- 12c -- evaluate target parameter:

        if (FALSE) update.list <- lapply(1:length(update.fit), function(each.index) {
            out <- rbind(tmle.est=update.fit[[each.index]]$update.est,
                         tmle.se=update.fit[[each.index]]$update.ic)
            colnames(out) <- paste0("tau=", tt)
            return(out)
        }) else update.list <- update.fit

        if (cr) {
            if (length(target)==length(update.list)) {
                names(update.list) <- paste0("F", target)
            } else {
                names(update.list) <- c(paste0("F", target), "S")
            }
        } else {
            update.list <- update.list[[1]]
        }

        #-- 12e -- compute sd:
        
        tmle.list$tmle <- update.list
    }

    if (length(not.fit.list)>0) {
        tmle.list[[length(tmle.list)+1]] <- not.fit.list
        names(tmle.list)[length(tmle.list)] <- "messages"
    }
 
    if (simultaneous.ci) {
        if (target.S) {
            Sigma <- eval.ic(mat, init.fit, target.index=outcome.index[target], Sigma=TRUE, survival=TRUE)
        } else {
            Sigma <- eval.ic(mat, init.fit, target.index=outcome.index[target], Sigma=TRUE)
        }
        rho <- matrix(0, nrow=nrow(Sigma), ncol=ncol(Sigma))
        for (j1 in 1:nrow(rho)) {
            for (j2 in 1:nrow(rho)) {
                rho[j1,j2] <- Sigma[j1,j2] / sqrt(Sigma[j1,j1]*Sigma[j2,j2])
            }
        }
        generate.max <- list()
        repeat.no <- 50000
        for (mm in 1:repeat.no) {
            set.seed(10303+mm)
            Z <-  try(mvrnorm(n=1, rep(0, nrow(rho)), rho, tol=1e-6, empirical=FALSE))
            if (any(class(Z)=="try-error")) {
                generate.max[[mm]] <- NA
                warning("NA in Sigma; cannot estimate q95")
                break
            } else {
                generate.max[[mm]] <- max(Z)
            }
        }
        # plot(unlist(lapply(1:length(generate.max), function(jj) quantile(abs(unlist(generate.max[1:jj])), p=0.95))))
        if (length(generate.max)==repeat.no) {
            q.max.95 <- quantile(abs(unlist(generate.max)), p=0.95)
            tmle.list$q.max.95 <- q.max.95
        }
    }

    if (FALSE) {
        mean(mat[A==1, unique(surv.tau1), by="id"][,2][[1]])
        1-sum(unlist(lapply(tmle.list$tmle.second.round, function(xx) xx["tmle.est",])))
        sqrt(sum(unlist(lapply(tmle.list$tmle.second.round, function(xx) xx["tmle.se",]^2))))
        update.ic <- sqrt(mean(rowSums(sapply(1:length(outcome.index), function(target11) {
            unlist(eval.ic(mat, unlist(lapply(tmle.list$tmle.second.round, function(xx) xx["tmle.est",])), target.index=outcome.index[target11],
                           tau.values=tau, survival=TRUE))
        }))^2)/n)#[[1]]
    }
  
    if (length(check.times.size)>0) {
        unique.times3 <- unique.times[unique.times<=max(tau)]
        unique.times3 <- unique.times3[unique.times3>=check.min & unique.times3<=check.max]
        set.seed(2444231)
        unique.times2 <- sort(unique.times3[sample(length(unique.times3), check.times.size)])

        ## print(unique.times2)
        ## print(tau)

        for (kk in 1:length(unique.times2)) {
            mat[, (paste0("surv.tau", kk)):=
                      surv.t[get(time.var)==max(get(time.var)[get(time.var)<=unique.times2[kk]])],
                by=c("id", A.name)]
            mat[surv.t>0, (paste0("Ht.lambda.", kk)):=get(paste0("surv.tau", kk)) / surv.t]
            mat[surv.t==0, (paste0("Ht.lambda.", kk)):=1]
        }

        all.fit <- sapply(1:length(unique.times2), function(kk) {
            mean(rowSums(sapply(a, function(aa)
            (2*(aa==a[1])-1)*(mat[get(A.name)==aa, 1-get(paste0("surv.tau", kk))[1], by="id"][,2][[1]]))))
        })

        all.ic <- eval.ic(mat, all.fit, target.index=outcome.index[target], tau.values=unique.times2, tau=unique.times2)
        all.ic[all.ic<delta.min] <- delta.min

        Pn.eic.all <- eval.equation(mat, 0, target.index=outcome.index[target1], tau.values=unique.times2, tau=unique.times2)
        
        Pn.eic.all4 <- lapply(1:length(Pn.eic.all), function(kk) Pn.eic.all[[kk]]/(ifelse(any(unlist(all.ic)==0), all.ic[[kk]]+0.001, all.ic[[kk]])*sqrt(n)))
        check.sup.norm.all <- abs(unlist(Pn.eic.all4))<=1/(sqrt(n)*log(n))#criterion#/ifelse(push.criterion, 1, sqrt(length(target)*length(tau)))

        tmle.list$check.sup.norm.all <- list(
            check.sup.norm.all=mean(check.sup.norm.all),
            lhs=max(abs(unlist(Pn.eic.all4))),
            rhs=1/(sqrt(n)*log(n)),#criterion/ifelse(push.criterion, 1, sqrt(length(target)*length(tau))),
            Pn.eic=data.table(times=unique.times2,
                              Pn.eic=unlist(Pn.eic.all4),
                              tf=1*check.sup.norm.all))

    }

    return(tmle.list)    
}





