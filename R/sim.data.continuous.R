sim.data2 <- function(n, setting=1, competing.risk=FALSE, no.cr=2, cr.setting=1, 
                      censoring.informative=FALSE,
                      no.effect.A=FALSE,
                      randomize.A=TRUE,
                      print.forms=FALSE,
                      m=sample(373211, 1)) {

    if (!no.effect.A) betaA <- -0.15 else betaA <- 0
    betaL <- 1.1
    nu    <- 1.7
    eta   <- 0.7
    t0    <- 0.9
    
    if (setting==1) {
        square.effect1        <- TRUE
        square.effect2        <- FALSE
        interaction.Atime     <- TRUE
        reversed.setting      <- TRUE
    } else {
        square.effect1        <- FALSE
        square.effect2        <- TRUE
        interaction.Atime     <- FALSE
        reversed.setting      <- FALSE
    }
 
    no.censoring          <- FALSE

    if (interaction.Atime & !no.effect.A) betaA <- -0.7
    
    if (reversed.setting) {
        if (!no.effect.A) betaA <- 0.5
        t0 <- 0.7
    }

    return(sim.data(n, betaA=betaA, betaL=betaL, nu=nu, eta=eta, t0=t0,
                    seed=m+100, no.cr=no.cr,
                    competing.risk=competing.risk,
                    categorical=FALSE, randomize.A=randomize.A,
                    censoring.informative=censoring.informative,
                    censoring=!no.censoring,
                    square.effect2=square.effect2,
                    square.effect1=square.effect1,
                    print.forms=print.forms,
                    reversed.setting=reversed.setting,
                    cr.setting=cr.setting,
                    interaction.Atime=interaction.Atime))
    
}

sim.data <- function(n, loop.max=20, endoffollowup=30,
                     betaA=0.1, betaL=0.3,
                     nu=0.5, eta=4/sqrt(2)*(1/8),
                     firstevent=TRUE,
                     censoring=TRUE,
                     competing.risk=FALSE, cr.both=FALSE, no.cr=2, cr.setting=1,
                     seed=sample(4034244, 1),
                     interaction.AL=FALSE,
                     interaction.Atime=FALSE, t0=0.3,
                     browse=FALSE, verbose=FALSE,
                     randomize.A=FALSE, square.effect=FALSE,
                     square.effect2=FALSE,
                     square.effect1=FALSE,
                     new=FALSE,
                     print.forms=FALSE,
                     censoring.informative=TRUE, censoring.high=FALSE, 
                     categorical=TRUE, intervention.A=NULL, tau=2,
                     cvot.setting=FALSE, reversed.setting=FALSE
                     ) {
    
    set.seed(seed)

    if (square.effect1) square.effect2 <- TRUE

    censoring.alpha <- -2.1
    if (interaction.AL | square.effect) censoring.alpha <- -1.7
    if (cvot.setting) censoring.alpha <- -1.4
    if (reversed.setting) censoring.alpha <- -1.7
    if (square.effect2) censoring.alpha <- -2.5
    if (competing.risk & !interaction.Atime) censoring.alpha <- -0.6 else if (competing.risk & censoring.informative) censoring.alpha <- -2.3 else if (competing.risk) censoring.alpha <- -0.6
    if (censoring.high) censoring.alpha <- -1 

    if (firstevent) loop.max <- 1

    if (competing.risk & square.effect1) alphaT <- 0.1 else alphaT <- -0.6

    if (length(intervention.A)>0) censoring <- FALSE
    
    if (categorical) {
        Lcont <- runif(n, -3, 3)#rbinom(n, 1, 0.5)
        Lgrid <- seq(min(Lcont), max(Lcont), length=5)    
        L1 <- findInterval(Lcont, Lgrid)/10
        L2 <- rbinom(n, 1, 1/2)
        L3 <- rbinom(n, 1, 1/2)
    } else {
        #L <- runif(n, -3, 3)/3 #rbinom(n, 1, 0.5)
        if (square.effect) L1 <- rnorm(n, sd=0.5) else if (square.effect2) L1 <- runif(n, -1, 1) else L1 <- runif(n, 0, 1)
        L2 <- runif(n, 0, 1)
        L3 <- runif(n, 0, 1) 
    }

    if (interaction.AL | interaction.Atime) {
        if (FALSE) L1 <- runif(n, -1, 1) else  L1 <- runif(n, 0, 1)
        L2 <- runif(n, 0, 1)
        if (square.effect2 | !reversed.setting) {
            L2 <- runif(n, -1, 1)
            L2 <- sign(L2)*sqrt(abs(L2))
        }
        if (square.effect1) {
            L1 <- runif(n, -1, 1)
        }
        if ((interaction.AL & interaction.Atime) | new) L3 <- rbinom(n, 1, 0.35) else L3 <- runif(n, 0, 1)
        if (randomize.A) A <- rbinom(n, 1, plogis(qlogis(0.5))) else if (categorical)
                                                                    A <- rbinom(n, 1, plogis(0.4+0.3*L1)) else A <- rbinom(n, 1, plogis(0.4+0.3*L1-0.3*L2))
    } else {
        if (randomize.A) A <- rbinom(n, 1, plogis(qlogis(0.5))) else A <- rbinom(n, 1, plogis(-0.9-0.3*L1-1.5*(L2<0.2)))
        if (print.forms & !randomize.A) {
            print("treatment model:")
            print("(-0.9-0.3*L1-0.2*(L2<0.2))")
        }
    }
    
    if (length(intervention.A)>0 & is.numeric(intervention.A)) A <- intervention.A else if (length(intervention.A)>0 & is.function(intervention.A)) A <- rbinom(n, 1, intervention.A(cbind(L1,L2,L3)))

    if (interaction.AL & !interaction.Atime) {
        phiT <- function(t, A, L1, L2, L3, betaA, betaL) {
            return(exp(A*betaA+L1^2*betaL-0.15*A*L1^2+0.75*L2*L1-1.2*L3))
        }
    } else if (interaction.Atime) {
        if (square.effect2) {
            if (reversed.setting) {
                phiT <- function(t, A, L1, L2, L3, betaA, betaL) {
                    return(exp(
                    (t<=t0)*betaA*A+
                    (t>t0)*(-0.45)*betaA*A+
                    -1.2*L1^2+
                    alphaT ))
                }
            } else {
                phiT <- function(t, A, L1, L2, L3, betaA, betaL) {
                    return(exp(
                    (t<=t0)*betaA*A+
                    (t>t0)*(-0.45)*betaA*A-
                    L1*betaL-1.2*L2^2+0.8*L3+
                    + 0.1 ))
                }
            }
            if (verbose) print(paste0("time-varying HR, t0=", t0, ", betaA=", betaA, ", period2=", betaA*-0.45))
        } else if (interaction.AL) {
            phiT <- function(t, A, L1, L2, L3, betaA, betaL) {
                return(exp(
                (t<=t0)*betaA*A*(3.5*L3)+
                (t<=t0)*
                2.2*betaL*L3+
                (t>t0)*(0)*betaA*A-
                0.1*L1*1.0-0.1*0.6*L2+
                + 1.2))
            }
        } else if (new) {
            phiT <- function(t, A, L1, L2, L3, betaA, betaL) {
                return(exp(
                    betaA*A*(3.5*L3)+
                    4.6*betaL*L3+
                    (0)*betaA*A-
                    0.1*L1*1.0-0.1*0.6*L2+
                    + 0.8))
            }
        } else {
            phiT <- function(t, A, L1, L2, L3, betaA, betaL) {
                return(exp(
                (t<=t0)*betaA*A+
                (t>t0)*(-0.45)*betaA*A-
                L1*betaL-1.2*L2+0.8*L3+
                + 0.1 ))
            }
            if (verbose) print(paste0("time-varying HR, t0=", t0, ", betaA=", betaA, ", period2=", betaA*-0.45))
        }
    } else if (square.effect | square.effect2) {
        phiT <- function(t, A, L1, L2, L3, betaA, betaL) {
            return(exp(A*betaA+L1^2*1.2))
        }
        if (verbose) print("1.2*L1^2")
    } else {
        phiT <- function(t, A, L1, L2, L3, betaA, betaL) {
            return(exp(A*betaA+L1^2*betaL+0.75*sqrt(L2)*L1-1.2*sin(L3*6)))
        }
        if (verbose) print("sin(L3)")
    }
   
    lambdaT <- function(t, A, L1, L2, L3, betaA, betaL, eta, nu) {
        return(phiT(t, A, L1, L2, L3, betaA, betaL)*eta*nu*t^{nu-1})
    }

    if (!censoring.informative) {
        phiC <- function(t, A, L1, L2, L3) {
            return(exp(-0.1+censoring.alpha))
        }
    } else {
        if (verbose) print("informative censoring")
        if (interaction.Atime & (interaction.AL)) {
            phiC <- function(t, A, L1, L2, L3) {
                return(exp(-L3*0.1+2.6*L3 + 0.1 + censoring.alpha))
            }
        } else if (interaction.Atime & new) {
            phiC <- function(t, A, L1, L2, L3) {
                return(exp(-L3*0.8+2.2*L3 + 1.1 + censoring.alpha))
            }
        } else if (interaction.Atime & reversed.setting & square.effect1) {
            if (competing.risk) {
                phiC <- function(t, A, L1, L2, L3) {
                    return(exp(-L3*0.8+1.2*L1^2*A + 1.5 + censoring.alpha))
                }
            } else {
                phiC <- function(t, A, L1, L2, L3) {
                    return(exp(-L3*0.8+1.2*L1^2*(1-A) + 1.5 + censoring.alpha))
                }
            }
        }  else if (interaction.Atime & reversed.setting) {
            phiC <- function(t, A, L1, L2, L3) {
                return(exp(-L3*0.8+1.2*L1*(1-A) + 1.1 + censoring.alpha))
            }
        } else if (interaction.Atime & cvot.setting) {
            phiC <- function(t, A, L1, L2, L3) {
                return(exp(-L3*0.8-1.2*L1*A + 1.1 + censoring.alpha))
            }
        } else {
            phiC <- function(t, A, L1, L2, L3) {
                return(exp(-L3*0.8+1.2*L1*A + 1.1 + censoring.alpha))
            }
        }
    }
   
    lambdaC <- function(t, A, L1, L2, L3, eta, nu) {
        return(phiC(t, A, L1, L2, L3)*eta*nu*t^{nu-1})
    }

    if (square.effect2 & !interaction.Atime) {
        if (betaA==0) betaA2 <- 0.8 else betaA2 <- 0.4
        if (cr.setting==1) {
            phiT2 <- function(t, A, L1, L2, L3) {
                return(exp(+0.4+0.7*L1-betaA2*A))
            }
        } else {
            phiT2 <- function(t, A, L1, L2, L3) {
                return(exp(+0.4-0.7*L1+betaA2*A))
            }
        }
    } else if (square.effect1 & interaction.Atime) {
        if (betaA==0) betaA2 <- 0.8 else betaA2 <- 0.3
        if (cr.setting==1) {
            phiT2 <- function(t, A, L1, L2, L3) {
                return(exp(-0.2+0.7*L1-0.4*L3+betaA2*A))
            }
        } else {
            phiT2 <- function(t, A, L1, L2, L3) {
                return(exp(-0.2-0.7*L1-0.4*L3-betaA2*A))
            }
        }
    } else {
        if (betaA==0) betaA2 <- 0.8 else betaA2 <- 0.4
        if (cr.setting==1) {
            phiT2 <- function(t, A, L1, L2, L3) {
                return(exp(-1.4+0.7*L1-betaA2*A))
            }
        } else {
            phiT2 <- function(t, A, L1, L2, L3) {
                return(exp(-1.4-0.7*L1+betaA2*A))
            }
        }
    }
    
    lambdaT2 <- function(t, A, L1, L2, L3, eta, nu) {
        return(phiT2(t, A, L1, L2, L3)*eta*nu*t^{nu-1})
    }

    phiT3 <- function(t, A, L1, L2, L3) phiT2(t, A, L1, L2, L3)
    
    lambdaT3 <- function(t, A, L1, L2, L3, eta, nu) {
        return(phiT3(t, A, L1, L2, L3)*eta*nu*t^{nu-1})
    }

    if (print.forms) {
        print("Outcome hazard form: ")
        print(phiT)
        print(paste0("betaA=", betaA))
        if (length(t0)>0) print(paste0("t0=", t0))
        print("Censoring hazard form: ")
        print(phiC)
        if (competing.risk) {
            print("Cause-2 hazard form: ")
            print(phiT2)
            print(paste0("betaA2=", betaA2))
            if (no.cr>2) {
                print("Cause-3 hazard form: ")
                print(phiT3)
            }
        }
        print(paste0("shape parameter (eta) = ", eta))
        print(paste0("scale parameter (nu) = ", nu))
    }


    if (censoring) {
        if (competing.risk) {
            if (no.cr==3) {
                phi <- function(t, A, L1, L2, L3, betaA, betaL) phiT(t, A, L1, L2, L3, betaA, betaL) +
                                                                    phiC(t, A, L1, L2, L3) +
                                                                    phiT2(t, A, L1, L2, L3) +
                                                                    phiT3(t, A, L1, L2, L3)
            } else {
                phi <- function(t, A, L1, L2, L3, betaA, betaL) phiT(t, A, L1, L2, L3, betaA, betaL) +
                                                                    phiC(t, A, L1, L2, L3) +
                                                                    phiT2(t, A, L1, L2, L3)
            } 
        } else {
            phi <- function(t, A, L1, L2, L3, betaA, betaL) phiT(t, A, L1, L2, L3, betaA, betaL) +
                                                                phiC(t, A, L1, L2, L3)
        }
    } else {
        if (competing.risk) {
            if (no.cr==3) {
                phi <- function(t, A, L1, L2, L3, betaA, betaL) phiT(t, A, L1, L2, L3, betaA, betaL) +
                                                                    phiT2(t, A, L1, L2, L3) +
                                                                    phiT3(t, A, L1, L2, L3)
            } else {
                phi <- function(t, A, L1, L2, L3, betaA, betaL) phiT(t, A, L1, L2, L3, betaA, betaL) +
                                                                    phiT2(t, A, L1, L2, L3)
            } 
        } else {
            phi <- function(t, A, L1, L2, L3, betaA, betaL) phiT(t, A, L1, L2, L3, betaA, betaL)
        }
    }
        
    if (interaction.Atime) {
        if (interaction.AL) {
            Lambda.inv <- function(u, t, A, L1, L2, L3, nu, eta) {
                return( rowSums(cbind((u <= (eta*phi(t, A, L1, L2, L3, betaA, betaL))*t0^{nu}) *
                                      (( (u + eta*phi(t, A, L1, L2, L3, betaA, betaL*0)*t^{nu}) /
                                         (eta*phi(t, A, L1, L2, L3, betaA, betaL*0)) )^{1/nu} - t),
                (u > (eta*phi(t, A, L1, L2, L3, betaA, betaL))*t0^{nu}) *
                (( (u - (eta*phi(t, A, L1, L2, L3, betaA, betaL))*t0^{nunu} +
                    eta*phi(t, A, L1, L2, L3,-0*betaA, betaL)*t0^{nu}) /
                   (eta*phi(t, A, L1, L2, L3,-0*betaA, betaL)) )^{1/nu} - t)), na.rm=TRUE) )
            }
        } else if (new) {
            Lambda.inv <- function(u, t, A, L1, L2, L3, nu, eta) {
                return( rowSums(cbind((u <= (eta*phi(t, A, L1, L2, L3, betaA, betaL))*t0^{nu}) *
                                      (( (u + eta*phi(t, A, L1, L2, L3, betaA, betaL)*t^{nu}) /
                                         (eta*phi(t, A, L1, L2, L3, betaA, betaL)) )^{1/nu} - t),
                (u > (eta*phi(t, A, L1, L2, L3, betaA, betaL))*t0^{nu}) *
                (( (u - (eta*phi(t, A, L1, L2, L3, betaA, betaL))*t0^{nu} +
                    eta*phi(t, A, L1, L2, L3,betaA, betaL)*t0^{nu}) /
                   (eta*phi(t, A, L1, L2, L3,betaA, betaL)) )^{1/nu} - t)), na.rm=TRUE) )
            }
        } else  {
            Lambda.inv <- function(u, t, A, L1, L2, L3, nu, eta) {
                return( rowSums(cbind((u <= (eta*phi(t, A, L1, L2, L3, betaA, betaL))*t0^{nu}) *
                                      (( (u + eta*phi(t, A, L1, L2, L3, betaA, betaL)*t^{nu}) /
                                         (eta*phi(t, A, L1, L2, L3, betaA, betaL)) )^{1/nu} - t),
                (u > (eta*phi(t, A, L1, L2, L3, betaA, betaL))*t0^{nu}) *
                (( (u - (eta*phi(t, A, L1, L2, L3, betaA, betaL))*t0^{nu} +
                    eta*phi(t, A, L1, L2, L3,-0.45*betaA, betaL)*t0^{nu}) /
                   (eta*phi(t, A, L1, L2, L3,-0.45*betaA, betaL)) )^{1/nu} - t)), na.rm=TRUE) )
            }
        }
        Lambda.inv1 <- function(u, t, A, L1, L2, L3, nu, eta) {
            return(( (u + eta*phi(t, A, L1, L2, L3, betaA, betaL)*t^{nu}) /
                     (eta*phi(t, A, L1, L2, L3, betaA, betaL)) )^{1/nu} - t)
        }
    } else {
        Lambda.inv <- function(u, t, A, L1, L2, L3, nu, eta) {
            return(( (u + eta*phi(t, A, L1, L2, L3, betaA, betaL)*t^{nu}) /
                     (eta*phi(t, A, L1, L2, L3, betaA, betaL)) )^{1/nu} - t)
        }
    }
    
    #-- intialize monitoring times: 
    Tlist <- list(cbind(time=rep(0, n), delta=rep(0, n), id=1:n))

    #-- function to loop over for time-points: 
    loop.fun <- function(k, Tprev) {

        #-- simulate event time: 
        U <- -log(runif(n))
        Tout <- Lambda.inv(U, Tprev, A, L1, L2, L3, nu, eta) + Tprev

        #-- which event:
        if (censoring) {
            if (competing.risk) {
                denom <- (lambdaT(Tout, A, L1, L2, L3, betaA, betaL, eta, nu) +
                          lambdaC(Tout, A, L1, L2, L3, eta, nu) +
                          lambdaT2(Tout, A, L1, L2, L3, eta, nu))
                probT <- lambdaT(Tout, A, L1, L2, L3, betaA, betaL, eta, nu) / (denom)
                probC <- lambdaC(Tout, A, L1, L2, L3, eta, nu) / (denom)
                probT2 <- lambdaT2(Tout, A, L1, L2, L3, eta, nu) / (denom)
                if (no.cr==3) {
                    probT3 <- lambdaT3(Tout, A, L1, L2, L3, eta, nu) / (denom)
                    which <- apply(cbind(probC, probT, probT2, probT3), 1, function(p) sample(0:3, size=1, prob=p))
                } else {
                    which <- apply(cbind(probC, probT, probT2), 1, function(p) sample(0:2, size=1, prob=p))
                }
            } else {
                denom <- (lambdaT(Tout, A, L1, L2, L3, betaA, betaL, eta, nu) +
                          lambdaC(Tout, A, L1, L2, L3, eta, nu))
                probT <- lambdaT(Tout, A, L1, L2, L3, betaA, betaL, eta, nu) / (denom)
                probC <- lambdaC(Tout, A, L1, L2, L3, eta, nu) / (denom)
                which <- apply(cbind(probC, probT), 1, function(p) sample(0:1, size=1, prob=p))
            }
        } else {
            if (competing.risk) {
                if (no.cr==3) {
                    denom <- (lambdaT(Tout, A, L1, L2, L3, betaA, betaL, eta, nu) +
                              lambdaT2(Tout, A, L1, L2, L3, eta, nu)+
                              lambdaT3(Tout, A, L1, L2, L3, eta, nu))
                } else {
                    denom <- (lambdaT(Tout, A, L1, L2, L3, betaA, betaL, eta, nu) +
                              lambdaT2(Tout, A, L1, L2, L3, eta, nu))
                }
                probT <- lambdaT(Tout, A, L1, L2, L3, betaA, betaL, eta, nu) / (denom)
                probT2 <- lambdaT2(Tout, A, L1, L2, L3, eta, nu) / (denom)
                if (no.cr==3) {
                    probT3 <- lambdaT3(Tout, A, L1, L2, L3, eta, nu) / (denom)
                    which <- apply(cbind(probT, probT2, probT3), 1, function(p) sample(1:3, size=1, prob=p))
                } else {
                    which <- apply(cbind(probT, probT2), 1, function(p) sample(1:2, size=1, prob=p))
                }
            } else {
                which <- 1
            }
        }
        
        #-- return event time
        return(cbind(time=Tout, delta=which, id=1:n))
    }

    #-- run simulations: 
    for (k in 1:loop.max) {
        Tlist[[k+1]] <- loop.fun(k, Tlist[[k]][,1])
    }
    
    #-- collect data:
    dt <- data.table(do.call("rbind", Tlist))

    #-- order & throw away obs after end of followup: 
    setorder(dt, id, time, delta)
    dt <- dt[time<=endoffollowup]

    #-- only first event:
    if (firstevent) {
        dt[time>0, idN:=1:.N, by="id"]
        dt <- dt[idN==1][, -"idN", with=FALSE]
        
    }
    
    if (browse) browser()
    
    #-- merge with baseline information:
    baseline <- data.table(A=A, L1=L1, L2=L2, L3=L3, id=as.numeric(1:n))
    setkey(baseline, id); setkey(dt, id)
    dt <- merge(dt, baseline, by="id")

    if (verbose & !length(intervention.A)>0) {
        par(mfrow=c(1,2))
        dt[A==1 & delta==1, hist(time)]
        dt[A==0 & delta==1, hist(time)]
    } else if (verbose){
        dt[, hist(time)]
    }

    if (length(intervention.A)>0 & (competing.risk & cr.both==TRUE)) {
        return(do.call("rbind", lapply(tau, function(tau.kk) {
            out <- do.call("rbind", lapply(1:no.cr, function(each) {
                c(tau=tau.kk,
                  psi0=mean(dt[, time<=tau.kk & delta==each]))
            }))
            rownames(out) <- paste0("F", 1:no.cr)
            return(out)
        })))
    } else if (length(intervention.A)>0) {
        return(do.call("rbind", lapply(tau, function(tau.kk) {
            c(tau=tau.kk, S=mean(dt[, time<=tau.kk & delta==1]))
        })))
    } else return(dt)
}
