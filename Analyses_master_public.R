#!/usr/bin/env R

# 2022-11-17
# R-4.0.0

library("mice")
library("survival")
library("CPE")

source("Functions_for_manuscript.R")

mymid      <- readRDS("imputed.data.Rds")

#######################################################################################################################################################################################################
# Extract data for ER-negative and ER-positive subgroups

erstatus                                            <- 0 # set either 0 for negative or 1 for positive
m                                                   <- 50 # number of imputed data
brca12                                              <- c("BRCA1", "BRCA2")


discdata                                            <- vector("list", length = m)
for(i in 1:m) { 
               temp                                 <- complete(mymid, action=i)
               temp                                 <- subset(temp, ER==erstatus & M!=1)
               temp$bcdeathbc[which(temp$VTSbc==1)] <- temp$bcdeath[which(temp$VTSbc==1)]
               temp$GRD                             <- as.numeric(temp$GRD)
               discdata[[i]]                        <- temp
               }

#######################################################################################################################################################################################################
# Prognostic index based on PREDICT
# In ER+ analysis, a reduced score without tumor grade was calculated.

if(erstatus==0){
    discdata            <- lapply(discdata, ERNEGSCORE)
} else {
    discdata            <- lapply(discdata, ERPOSSCORE)
    discdata.grd        <- lapply(discdata, ERPOSGRD)
}

#######################################################################################################################################################################################################
# Adjust the analysis model

offset.results          <- lapply(discdata, function(dat) {
                                               coxph(Surv(TTE, FUTbc, as.numeric(bcdeathbc=="1")) ~ AgeDiagIndex + 
                                                                                                    YearDiagIndex + 
                                                                                                    GRD +
                                                                                                    PR +
                                                                                                    logsize  +
                                                                                                    nodecount +
                                                                                                    HER2 +
                                                                                                    strata(strat) +
                                                                                                    offset(PREDICT), 
                                                                                                    data = subset(dat, BRCA12==brca12[erstatus+1]))
                                                                                            })                                         
offset.results.pooled   <- pool.splinemodel(offset.results, m=50)

if(erstatus==1){
    offset.grd          <- lapply(discdata.grd, function(dat) {
                                                   coxph(Surv(TTE, FUTbc, as.numeric(bcdeathbc=="1")) ~ factor(GRD, levels=c(2,1,3)) + 
                                                                                                        YearDiagIndex + 
                                                                                                        strata(strat) +
                                                                                                        offset(PREDICT), 
                                                                                                        data = subset(dat, BRCA12==brca12[erstatus+1]))
                                                                                            })                                         
    offset.grd.pooled   <- pool.splinemodel(offset.grd, m=50)
}

#######################################################################################################################################################################################################
# Inspect the hazard associated with diagnosis age, when other factors are standardized
# Use offset analysis, predict hazard for dummy data, pool, and plot.

offsetage               <- lapply(discdata, function(dat) {
                                               coxph(Surv(TTE, FUTbc, as.numeric(bcdeathbc=="1")) ~ pspline(AgeDiagIndex) + 
                                                                                                    YearDiagIndex + 
                                                                                                    strata(strat) +
                                                                                                    offset(PREDnoage), 
                                                                                                    data = subset(dat, BRCA12==brca12[erstatus+1]))
                                                                                            })                                         

agepdata                <- expand.grid(AgeDiagIndex= 19:69, PREDnoage=0, strat="UK")

agerefdata              <- data.frame(AgeDiagIndex=50, PREDnoage=0, strat="UK")
                                         
agepred                 <- lapply(offsetage, function(mymod) {
                                                  temppred         <- predict(mymod, newdata=agepdata, se=TRUE)
                                                  tempyy           <- matrix(temppred$fit + outer(temppred$se, c(0, -1.96, 1.96),'*'), ncol=3)
                                                  colnames(tempyy) <- c("RR_age", "CI95L_age", "CI95U_age")
                                                  tempref          <- predict(mymod, newdata=agerefdata, se=FALSE)
                                                  tempyy           <- tempyy-tempref
                                                  return(list(age=tempyy))
                                                  })

agepooledpred           <- make.pred.splinemodel(agepred, ncat=1, m=50)

if(erstatus==0){
    agePREDICTpred          <- unlist(lapply(19:69, function(age) 0.0089*(age - 56.3))) - 0.0089*(50 - 56.3)
} else {
    agePREDICTpred          <- unlist(lapply(19:69, function(age) 34.53*((age/10)^-2 - 0.0287) + (-34.20)*((age/10)^-2 * log(age/10) - 0.0510))) - (34.53*((50/10)^-2 - 0.0287) + (-34.20)*((50/10)^-2 * log(50/10) - 0.0510))
}

#######################################################################################################################################################################################################
# Gönen & Heller unbiased concordance

reregres                <- list()

reregres[[1]]           <- lapply(discdata, function(mydat) coxph(Surv(TTE, FUTbc, bcdeathbc==1) ~ PREDICT + strata(strat), data=subset(mydat, BRCA12==brca12[erstatus+1])))
reregres[[2]]           <- lapply(discdata, function(mydat) coxph(Surv(TTE, FUTbc, bcdeathbc==1) ~ PREDICT + strata(strat), data=subset(mydat, BRCA12==brca12[-(erstatus+1)])))

GoHetab                 <- do.call(cbind, lapply(1:2, function(j) unlist(lapply(1:m, function(i) phcpe(reregres[[j]][[i]], CPE.SE=FALSE, out.ties=FALSE)))))

GoHepooled              <- cbind(apply(GoHetab, 2, median), apply(GoHetab, 2, IQR))
colnames(GoHepooled)    <- c("median", "IQR")

# In ER+ analysis, concordance is estimated also for the reduced model, excluding tumor grade
if(erstatus==1){
    reregres.grd                <- list()
    reregres.grd[[1]]           <- lapply(discdata.grd, function(mydat) coxph(Surv(TTE, FUTbc, bcdeathbc==1) ~ PREDICT + strata(strat), data=subset(mydat, BRCA12==brca12[erstatus+1])))
    GoHetab.grd                 <- do.call(cbind, lapply(1, function(j) unlist(lapply(1:m, function(i) phcpe(reregres.grd[[j]][[i]], CPE.SE=FALSE, out.ties=FALSE)))))
    GoHepooled.grd              <- cbind(apply(GoHetab.grd, 2, median), apply(GoHetab, 2, IQR))
    colnames(GoHepooled.grd)    <- c("median", "IQR")
}
#######################################################################################################################################################################################################
# CALIBRATION 

califunc                <- ifelse(erstatus==0, PREDcalineg, PREDcalipos)

cali_cats1              <- lapply(lapply(discdata, function(l) subset(l, BRCA12==brca12[erstatus+1])), califunc, ncats=1)
calipooled1_10y         <- poolcali(cali_cats1, "expected10y", "O10", "SE10")

N0_cali_cats1           <- lapply(lapply(discdata, function(l) subset(l, BRCA12==brca12[erstatus+1] & nodecount==0)), califunc, ncats=1)
N0_calipooled1_10y      <- poolcali(N0_cali_cats1, "expected10y", "O10", "SE10")

N1_cali_cats1           <- lapply(lapply(discdata, function(l) subset(l, BRCA12==brca12[erstatus+1] & nodecount>1)), califunc, ncats=1)
N1_calipooled1_10y      <- poolcali(N1_cali_cats1, "expected10y", "O10", "SE10")

GRD3_cali_cats1         <- lapply(lapply(discdata, function(l) subset(l, BRCA12==brca12[erstatus+1] & GRD==3)), califunc, ncats=1)
GRD3_calipooled1_10y    <- poolcali(GRD3_cali_cats1, "expected10y", "O10", "SE10")

GRD2_cali_cats1         <- lapply(lapply(discdata, function(l) subset(l, BRCA12==brca12[erstatus+1] & GRD==2)), califunc, ncats=1)
GRD2_calipooled1_10y    <- poolcali(GRDO_cali_cats1, "expected10y", "O10", "SE10")

T1_cali_cats1           <- lapply(lapply(discdata, function(l) subset(l, BRCA12==brca12[erstatus+1] & T==1)), califunc, ncats=1)
T1_calipooled1_10y      <- poolcali(T1_cali_cats1, "expected10y", "O10", "SE10")

TO_cali_cats1           <- lapply(lapply(discdata, function(l) subset(l, BRCA12==brca12[erstatus+1] & T!=1)), califunc, ncats=1)
TO_calipooled1_10y      <- poolcali(TO_cali_cats1, "expected10y", "O10", "SE10")

PR1_cali_cats1          <- lapply(lapply(discdata, function(l) subset(l, BRCA12==brca12[erstatus+1] & PR==1)), califunc, ncats=1)   
PR1_calipooled1_10y     <- poolcali(PR1_cali_cats1, "expected10y", "O10", "SE10")

PRO_cali_cats1          <- lapply(lapply(discdata, function(l) subset(l, BRCA12==brca12[erstatus+1] & PR==0)), califunc, ncats=1)
PRO_calipooled1_10y     <- poolcali(PRO_cali_cats1, "expected10y", "O10", "SE10")

HER21_cali_cats1        <- lapply(lapply(discdata, function(l) subset(l, BRCA12==brca12[erstatus+1] & HER2==1)), califunc, ncats=1)
HER21_calipooled1_10y   <- poolcali(HER21_cali_cats1, "expected10y", "O10", "SE10")

HER20_cali_cats1        <- lapply(lapply(discdata, function(l) subset(l, BRCA12==brca12[erstatus+1] & HER2==0)), califunc, ncats=1)
HER20_calipooled1_10y   <- poolcali(HER20_cali_cats1, "expected10y", "O10", "SE10")

YOUNG_cali_cats1        <- lapply(lapply(discdata, function(l) subset(l, BRCA12==brca12[erstatus+1] & AgeDiagIndex<35)), califunc, ncats=1)
YOUNG_calipooled1_10y   <- poolcali(YOUNG_cali_cats1, "expected10y", "O10", "SE10")

MIDD_cali_cats1         <- lapply(lapply(discdata, function(l) subset(l, BRCA12==brca12[erstatus+1] & AgeDiagIndex>=35 & AgeDiagIndex<45)), califunc, ncats=1)
MIDD_calipooled1_10y    <- poolcali(MIDD_cali_cats1, "expected10y", "O10", "SE10")

OLD_cali_cats1          <- lapply(lapply(discdata, function(l) subset(l, BRCA12==brca12[erstatus+1] & AgeDiagIndex>=45)), califunc, ncats=1)
OLD_calipooled1_10y     <- poolcali(OLD_cali_cats1, "expected10y", "O10", "SE10")

cali_cats1              <- lapply(lapply(discdata, function(l) subset(l, BRCA12==brca12[-(erstatus+1)])), califunc, ncats=1)
calipooled1_10y         <- poolcali(cali_cats1, "expected10y", "O10", "SE10")

cali_cats4              <- lapply(lapply(discdata, function(l) subset(l, BRCA12==brca12[erstatus+1])), califunc, ncats=4)

cal4_5y                 <- poolcali(cali_cats4, "expected5y", "O5", "SE5")
cal4_10y                <- poolcali(cali_cats4, "expected10y", "O10", "SE10")
cal4_15y                <- poolcali(cali_cats4, "expected15y", "O15", "SE15")

#######################################################################################################################################################################################################
# Prepare data for KM-plots 

bcdeathstudies          <- as.character(unique(mymid$data$study)[which(unlist(lapply(unique(mymid$data$study), function(s) {
                               length(which(mymid$data$study==s & mymid$data$bcdeath %in% 1:2))/length(which(mymid$data$study==s & mymid$data$VTS==1))>0.8})))])

obsvotedTher            <- VOTEStoOBS(discdata, as.character(unique(mymid$data$study)), brca12[erstatus+1], 4)

obs1voted               <- VOTEStoOBS(discdata, bcdeathstudies, brca12[erstatus+1], 4)
obs2voted               <- VOTEStoOBS(discdata, bcdeathstudies, brca12[-(erstatus+1)], 3)

obs1voted               <- subset(survSplit(Surv(TTE, FUTbc, bcdeathbc) ~ ., data= obs1voted, cut=c(10), episode= "tg", id="newid"), tg==1)
obs2voted               <- subset(survSplit(Surv(TTE, FUTbc, bcdeathbc) ~ ., data= obs2voted, cut=c(10), episode= "tg", id="newid"), tg==1)

obs1votedKM             <- survfit(Surv(TTE, FUTbc, bcdeathbc==1) ~ PREDvotes, data=obs1voted)
obs2votedKM             <- survfit(Surv(TTE, FUTbc, bcdeathbc==1) ~ PREDvotes, data=obs2voted)



