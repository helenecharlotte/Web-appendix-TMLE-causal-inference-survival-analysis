data(colon)

#-- create factor variable out of relevant covariates:
colon$differ <- factor(colon$differ)        
colon$extent <- factor(colon$extent)    

#-- create a dataset with one row per subject:
colon.dt <- merge(setDT(colon)[etype==2], setDT(colon)[etype==1, c("id", "time", "etype", "status"), with=FALSE], by="id")
colon.dt[, time.death:=time.x]  
colon.dt[, time.recurrence:=time.y]     
colon.dt[, status.death:=status.x]       
colon.dt[, status.recurrence:=status.y]
colon.dt[status.death==1 & status.recurrence==0, event:=2]
colon.dt[time.recurrence<=time.death & status.recurrence==1, event:=1]
colon.dt[status.death==0 & status.recurrence==0, event:=0] 
colon.dt <- colon.dt[, !(names(colon.dt) %in% c("status", grep("\\.x|\\.y", names(colon.dt), value=TRUE))), with=FALSE]

#-- create survival data set with event=1 if dead (with/without recurrence) and event=0 if censored: 
colon.surv <- na.omit(colon.dt)      
colon.surv[, event:=status.death]          
colon.surv[, time:=time.death]   
colon.surv <- colon.surv[, !(names(colon.surv) %in% grep("\\.death|\\.recurrence", names(colon.surv), value=TRUE)), with=FALSE]  

#-- create competing risks dataset with event=1 if cancer recurrence, event=2 if death without recurrence, event=0 if censored: 
colon.cr <- na.omit(colon.dt)                 
colon.cr[, time:=min(time.death, time.recurrence), by="id"]


