#!/usr/bin/env R

# 2022-11-17
# R-4.0.0

#######################################################################################################################################################################################################
# POOLING REGRESSION FUNCTION 
#######################################################################################################################################################################################################

                                                                
pool.splinemodel           <- function(reslist, m) {
                                       
                                       ########################################################################################
                                       # Pooling the coefficients according to Rubin's rules, and wrapping them into a table. #
                                       ########################################################################################
                                       
                                       ncoefs                 <- unique(unlist(lapply(1:m, function(i) length(reslist[[i]]$coefficients))))

                                       if(length(ncoefs)!=1){
                                         stop("number of coefficients differs")
                                       } 
  
                                       modcoefs               <- do.call(cbind, lapply(1:m, function(i) reslist[[i]]$coefficients))              
                                       colnames(modcoefs)     <- paste("imp", 1:m, sep="")

                                       modvar                 <- do.call(rbind, lapply(1:ncoefs, function(j) unlist(lapply(1:m, function(i) reslist[[i]]$var[j,j]))))
                                       row.names(modvar)      <- names(reslist[[1]]$coefficients)
                                       colnames(modvar)       <- paste("imp", 1:m, sep="")
                                       
                                       dummydfs               <- unlist(lapply(1:m, function(i) reslist[[i]]$nevent - nrow(summary(reslist[[i]])$coefficients)))
                                       
                                       modpooled              <- lapply(1:ncoefs, function(i) pool.scalar(Q=modcoefs[i,], U=modvar[i,], n=mean(dummydfs))) 
                                       names(modpooled)       <- row.names(modcoefs)

                                       modpooltab             <- as.data.frame(do.call(rbind, lapply(modpooled, function(mylist) mylist[c("m", "qbar", "ubar", "t", "df", "fmi")]
                                                                                                     )), stringsAsFactors=FALSE)                                       

                                       modpooltab$HR          <- exp(unlist(modpooltab$qbar))
                                       modpooltab$ci95.lower  <- exp(unlist(modpooltab$qbar)-1.96*sqrt(unlist(modpooltab$ubar)))
                                       modpooltab$ci95.upper  <- exp(unlist(modpooltab$qbar)+1.96*sqrt(unlist(modpooltab$ubar)))

                                       ###########################################################################                                                              
                                       # Pooling of the z-statistics per coefficient as in Marshall et al. 2009. #
                                       ###########################################################################
                                       
                                       modriv                 <- unlist(lapply(modpooled, function(mylist) mylist["r"]))
                                       
                                       # Degrees of freedom (df) for F-distribution
                                       # k = 1  --> numerator degrees of freedom
                                       # v = (m-1)*(1+r^-1)^2 --> denominator degrees of freedom for the F statistic
                                       modv2                 <- (m-1)*(1+modriv^-1)^2  

                                       # Pooled W
                                       # W = (Q0-qbar)^2/t   --> Q0=0
                                       modchisqpooledF       <- unlist(modpooltab$qbar)^2/unlist(modpooltab$t)

                                       # P-value per parameter
                                       modPspooled           <- 1-pf(modchisqpooledF, df1=1, df2=modv2)  
                                       # The calculated statistic follows F-distribution, with numerator and denominator dfs

                                       # combine values for parameter significance testing
                                       modpoolsum            <- data.frame(row.names = row.names(summary(reslist[[1]])$coefficients))
                                       modpoolsum$riv        <- modriv
                                       modpoolsum$v1         <- 1
                                       modpoolsum$v2         <- modv2
                                       modpoolsum$F          <- modchisqpooledF
                                       modpoolsum$P          <- modPspooled

                                       return(list("coefficients" = modpoolsum, "conf.int" = modpooltab))
                                       }

#######################################################################################################################################################################################################
# POOLING PREDICTIONS FUNCTION 
#######################################################################################################################################################################################################

make.pred.splinemodel      <- function(myimpred, ncats, m){
                                    lapply(1:ncats, function(mycat){
                                       avgest                <- rowMeans(do.call(cbind, lapply(1:m, function(i) myimpred[[i]][[mycat]][,1])))
                                       avgW                  <- rowMeans(do.call(cbind, lapply(1:m, function(i) myimpred[[i]][[mycat]][,1]-myimpred[[i]][[mycat]][,2]))) # W=1.92*sqrt(U) - CI-width
                                       myB                   <- (1+1/m)*(1/(m-1))*rowSums(do.call(cbind, lapply(1:m, function(i) (myimpred[[i]][[mycat]][,1]-avgest)^2)))
                                       newCIwidth            <- sqrt(avgW^2 + 1.96*myB)
                                       pooledpred            <- cbind(avgest, avgest-newCIwidth, avgest+newCIwidth)
                                       colnames(pooledpred)  <- colnames(myimpred[[1]][[mycat]])
                                       return(pooledpred)
                                       })
                                     }

#######################################################################################################################################################################################################
# CALCULATING PREDICT SCORES
#######################################################################################################################################################################################################

ERNEGSCORE     <- function(mydat){
                                          PREDICTrisks   <- c(AGE=0.0089, sqrtsize=2.09, nodesfun=0.626, GRD=1.13, HER2=1)

                                           mydat$PREDCDR  <- unlist(lapply(1:nrow(mydat), function(i){
                                                               testriskvalues   <- mydat[i, c("AgeDiagIndex", "logsize", "nodecount", "GRD", "HER2")]
                                                               testrisknumvals  <- c(testriskvalues[1, "AgeDiagIndex"]-56.3,
                                                                                      sqrt(exp(testriskvalues[1, "logsize"])/100)-0.509,
                                                                                      log(((testriskvalues[1, "nodecount"]+1)/10))+1.09,
                                                                                      as.numeric(testriskvalues[1,"GRD"]>1),                                                                    
                                                                                      HER2 = ifelse(testriskvalues[1,"HER2"]==1, 0.241, -0.076))
                                                               return(sum(testrisknumvals*PREDICTrisks))
                                                             }))
                                           
                                           chemoscores     <- with(mydat, ifelse(tax==1, -0.446, ifelse(anthrac==1, -0.248, 0)))
                                           
                                           mydat$PREDnother                          <- mydat$PREDCDR
                                           mydat$PREDnother[which(mydat$PR==0)]      <- mydat$PREDCDR[which(mydat$PR==0)]+log(1.03)             # from Grootes et al. Eur J Canc 2022
                                           mydat$PREDnother[which(mydat$PR==1)]      <- mydat$PREDCDR[which(mydat$PR==1)]+log(0.80)             # from Grootes et al. Eur J Canc 2022

                                           mydat$PREDICT                             <- mydat$PREDnother+chemoscores
                                           mydat$PREDICT[which(mydat$trastu==1 & mydat$HER2==1)]     <- mydat$PREDICT[which(mydat$trastu==1 & mydat$HER2==1)]-0.357

                                           mydat$PREDnoage                           <- mydat$PREDICT-(mydat$AgeDiagIndex-56.3)*PREDICTrisks["AGE"]
                                           
                                         return(mydat)
                                         }

#######################################################################################################################################################################################################
#######################################################################################################################################################################################################

ERPOSSCORE    <- function(mydat){
                                 PREDICTrisks   <- c(Age1=34.5, Age2=-34.20, logsize=0.75, nodesfun=0.71, GRD=0.75, HER2=1)  
                                 
                                 mydat$PREDCDR  <- unlist(lapply(1:nrow(mydat), function(i){
                                                       testriskvalues   <- mydat[i, c("AgeDiagIndex", "logsize", "nodecount", "GRD", "HER2")]
                                                       testrisknumvals  <- c(Age1 = (testriskvalues[1, "AgeDiagIndex"]/10)^-2 - 0.0288,                                                    
                                                                             Age2 = (testriskvalues[1, "AgeDiagIndex"]/10)^-2 * log(testriskvalues[1, "AgeDiagIndex"]/10) - 0.0510,
                                                                             logsize = testriskvalues[1, "logsize"] - log(100) + 1.55,                                                     
                                                                             nodesfun = log((testriskvalues[1,"nodecount"] + 1)/10) + 1.39,                          
                                                                             GRD = testriskvalues[1,"GRD"],                                                                    
                                                                             HER2 = ifelse(testriskvalues[1,"HER2"]==1, 0.241, -0.076))                                                    
                                                        return(sum(testrisknumvals*PREDICTrisks))
                                                        }))

                                 chemoscores     <- with(mydat, ifelse(tax==1, -0.446, ifelse(anthrac==1, -0.248, 0)))

                                 mydat$PREDnother                                      <- mydat$PREDCDR
                                 mydat$PREDnother[which(mydat$PR==0)]                  <- mydat$PREDCDR[which(mydat$PR==0)]+log(1.30)             # from Grootes et al. Eur J Canc 2022
                                 mydat$PREDnother[which(mydat$PR==1)]                  <- mydat$PREDCDR[which(mydat$PR==1)]+log(0.94)             # from Grootes et al. Eur J Canc 2022

                                 mydat$PREDICT                                         <- mydat$PREDnother+chemoscores
                                 mydat$PREDICT[which(mydat$endocr==1)]                 <- mydat$PREDICT[which(mydat$endocr==1)]-0.386 
                                 mydat$PREDICT[which(mydat$trastu==1 & mydat$HER2==1)] <- mydat$PREDICT[which(mydat$trastu==1 & mydat$HER2==1)]-0.357

                                 mydat$PREDnoage <- unlist(lapply(1:nrow(mydat), function(i){
                                                        testriskvalues   <- cbind(mydat[i, c("logsize", "nodecount", "GRD", "HER2")])
                                                        testrisknumvals  <- c(logsize = testriskvalues[1, "logsize"] - log(100) + 1.55,
                                                                              nodesfun = log((testriskvalues[1,"nodecount"] + 1)/10) + 1.39,
                                                                              GRD = testriskvalues[1,"GRD"],
                                                                              HER2 = ifelse(testriskvalues[1,"HER2"]==1, 0.241, -0.076))
                                                        return(sum(testrisknumvals*PREDICTrisks[c("logsize", "nodesfun", "GRD", "HER2")]))
                                                        }))

                                 mydat$PREDnoage[which(mydat$PR==0)]     <- mydat$PREDnoage[which(mydat$PR==0)]+log(1.30)             # from Grootes et al. Eur J Canc 2022
                                 mydat$PREDnoage[which(mydat$PR==1)]     <- mydat$PREDnoage[which(mydat$PR==1)]+log(0.94)             # from Grootes et al. Eur J Canc 2022
                                                                                                                  
                                 mydat$PREDnoage                         <- mydat$PREDnoage+chemoscores
                                 mydat$PREDnoage[which(mydat$endocr==1)] <- mydat$PREDnoage[which(mydat$endocr==1)]-0.386
                                 mydat$PREDnoage[which(mydat$trastu==1 & mydat$HER2==1)]     <- mydat$PREDnoage[which(mydat$trastu==1 & mydat$HER2==1)]-0.357   
                                                                                                                                
                                 return(mydat)
                                 }

#######################################################################################################################################################################################################
#######################################################################################################################################################################################################
# PARTIAL PREDICT SCORE FOR ER+ BC. GRADE EXCLUDED.
# Drop the factors, which have opposite effects in the full-offset model (GRD) by replacing the data with average value.

ERPOSGRD    <- function(mydat){
                                 meanGRD         <- round(mean(mydat$GRD),1)
                                           
                                 mydat$PREDCDR   <- unlist(lapply(1:nrow(mydat), function(i){
                                                        testriskvalues   <- cbind(mydat[i, c("AgeDiagIndex", "logsize", "nodecount", "HER2")], GRD=meanGRD)
                                                        testrisknumvals  <- c(Age1 = (testriskvalues[1, "AgeDiagIndex"]/10)^-2 - 0.0288,
                                                                              Age2 = (testriskvalues[1, "AgeDiagIndex"]/10)^-2 * log(testriskvalues[1, "AgeDiagIndex"]/10) - 0.0510,
                                                                              logsize = testriskvalues[1, "logsize"] - log(100) + 1.55,
                                                                              nodesfun = log((testriskvalues[1,"nodecount"] + 1)/10) + 1.39,
                                                                              GRD = as.numeric(testriskvalues[1,"GRD"]),                                                                    
                                                                              HER2 = ifelse(testriskvalues[1,"HER2"]==1, 0.241, -0.076))
                                                        return(sum(testrisknumvals*PREDICTrisks))
                                                        }))


                                 chemoscores     <- with(mydat, ifelse(tax==1, -0.446, ifelse(anthrac==1, -0.248, 0)))

                                 mydat$PREDnother                          <- mydat$PREDCDR
                                 mydat$PREDnother[which(mydat$PR==0)]      <- mydat$PREDCDR[which(mydat$PR==0)]+log(1.30)             # from Grootes et al. Eur J Canc 2022
                                 mydat$PREDnother[which(mydat$PR==1)]      <- mydat$PREDCDR[which(mydat$PR==1)]+log(0.94)             # from Grootes et al. Eur J Canc 2022

                                 mydat$PREDICT                             <- mydat$PREDnother+chemoscores
                                 mydat$PREDICT[which(mydat$endocr==1)]     <- mydat$PREDICT[which(mydat$endocr==1)]-0.386 
                                 mydat$PREDICT[which(mydat$trastu==1 & mydat$HER2==1)]     <- mydat$PREDICT[which(mydat$trastu==1 & mydat$HER2==1)]-0.357

                                 mydat$PREDnoage <- unlist(lapply(1:nrow(mydat), function(i){
                                                        testriskvalues   <- cbind(mydat[i, c("logsize", "nodecount", "HER2")], GRD=meanGRD)
                                                        testrisknumvals  <- c(logsize = testriskvalues[1, "logsize"] - log(100) + 1.55,
                                                                                        nodesfun = log((testriskvalues[1,"nodecount"] + 1)/10) + 1.39,
                                                                                        GRD = testriskvalues[1,"GRD"],
                                                                                        HER2 = ifelse(testriskvalues[1,"HER2"]==1, 0.241, -0.076))
                                                    return(sum(testrisknumvals*PREDICTrisks[c("logsize", "nodesfun", "GRD", "HER2")]))
                                                    })) 
                                                                                                                  
                                 mydat$PREDnoage[which(mydat$PR==0)]       <- mydat$PREDnoage[which(mydat$PR==0)]+log(1.30)             # from Grootes et al. Eur J Canc 2022
                                 mydat$PREDnoage[which(mydat$PR==1)]       <- mydat$PREDnoage[which(mydat$PR==1)]+log(0.94)             # from Grootes et al. Eur J Canc 2022
                                                                                                                  
                                 mydat$PREDnoage                           <- mydat$PREDnoage+chemoscores
                                 mydat$PREDnoage[which(mydat$endocr==1)]   <- mydat$PREDnoage[which(mydat$endocr==1)]-0.386
                                 mydat$PREDnoage[which(mydat$trastu==1 & mydat$HER2==1)]     <- mydat$PREDnoage[which(mydat$trastu==1 & mydat$HER2==1)]-0.357   
                                           
                                 return(mydat)
                                 }

#######################################################################################################################################################################################################
# BASELINE HAZARD FOR ER-NEGATIVE BREAST CANCER AND BASELINE HAZARD FOR ALL CAUSE MORTALITY
#######################################################################################################################################################################################################
ernegbase                       <- function(fut) {
                                      baseline   <- -1.156 + 0.4707/fut^2 -3.514/fut
                                      return(baseline)
                                      }
                                      
erposbase                       <- function(fut) {
                                      baseline   <- 0.7424 - 7.530/fut^(1/2) - 1.813*log(fut)/fut^(1/2)
                                      return(baseline)
                                      }
                                      
ocmbase                         <- function(fut) {
                                      baseline   <- (-6.053) + 1.080*log(fut) + 0.3255*fut^(1/2)
                                      return(baseline)
                                      }

#######################################################################################################################################################################################################
# CALIBRATION FUNCTIONS 
#######################################################################################################################################################################################################
                                      
                                      
PREDcalineg        <- function(mydat, ncats, split=TRUE){

                                if(ncats==4){
                                  mydat$patcats       <- cut(mydat$PREDICT, breaks=quantile(mydat$PREDICT, probs=c(c(0,0.2,0.5,0.8,1)), labels=c("low", "middle-low", "middle-high", "high")), include.lowest=TRUE) # PI quantiles
                                } else if(ncats==3){ 
                                  mydat$patcats       <- cut(mydat$PREDICT, breaks=quantile(mydat$PREDICT, probs=c(c(0,0.3,0.7,1)), labels=c("low", "middle", "high")), include.lowest=TRUE) # PI quantiles
                                } else {
                                  ncats               <- 1
                                  mydat$patcats       <- cut(mydat$PREDICT, breaks=quantile(mydat$PREDICT, probs=c(c(0,1)), labels=c("low")), include.lowest=TRUE)
                                }
# expected absolute risks
                                mydat$exp3y         <- (1-exp(-exp(ernegbase(3) + mydat$PREDICT))) * exp(-exp(0.0698*((mydat$AgeDiagIndex/10)^2 - 34.234) + ocmbase(3)))
                                mydat$exp5y         <- (1-exp(-exp(ernegbase(5) + mydat$PREDICT))) * exp(-exp(0.0698*((mydat$AgeDiagIndex/10)^2 - 34.234) + ocmbase(5)))
                                mydat$exp10y        <- (1-exp(-exp(ernegbase(10) + mydat$PREDICT))) * exp(-exp(0.0698*((mydat$AgeDiagIndex/10)^2 - 34.234) + ocmbase(10)))
                                mydat$exp15y        <- (1-exp(-exp(ernegbase(15) + mydat$PREDICT))) * exp(-exp(0.0698*((mydat$AgeDiagIndex/10)^2 - 34.234) + ocmbase(15)))

# return cumhaz in patcats
                                expected3y          <- unlist(lapply(1:ncats, function(i) with(subset(mydat, patcats==levels(patcats)[i]), mean(-log(1-exp3y)))))
                                expected5y          <- unlist(lapply(1:ncats, function(i) with(subset(mydat, patcats==levels(patcats)[i]), mean(-log(1-exp5y)))))
                                expected10y         <- unlist(lapply(1:ncats, function(i) with(subset(mydat, patcats==levels(patcats)[i]), mean(-log(1-exp10y)))))
                                expected15y         <- unlist(lapply(1:ncats, function(i) with(subset(mydat, patcats==levels(patcats)[i]), mean(-log(1-exp15y)))))

                                
                                expectedall         <- cbind(E3=as.data.frame(expected3y),
                                                             E5=as.data.frame(expected5y),
                                                             E10=as.data.frame(expected10y),
                                                             E15=as.data.frame(expected15y))

                                if(split){
                                  mydatsplit          <- survSplit(Surv(TTE, FUTbc, bcdeathbc) ~ ., 
                                                                   data= mydat, cut=c(3,5,10),
                                                                   episode= "tg", 
                                                                   id="newid")
                                  if(length(with(mydatsplit, which(FUTbc-TTE < 0.001)))>0){ 
                                    mydatsplit          <- mydatsplit[-with(mydatsplit, which(FUTbc-TTE < 0.001)),]
                                    }
                                    
                                  observedall         <- as.data.frame(do.call(cbind, lapply(1:4, function(j){
                                                                       tempsurvi           <- survfit(Surv(TTE, FUTbc, bcdeathbc==1) ~ patcats, data=subset(mydatsplit, tg %in% 1:j))
                                                                       if(ncats==1){
                                                                         return(cbind(tail(tempsurvi$cumhaz,1),
                                                                                      tail(tempsurvi$std.err,1)))
                                                                       } else {
                                                                         return(cbind(unlist(lapply(1:ncats, function(i) tempsurvi$cumhaz[sum(tempsurvi$strata[1:i])])),
                                                                                      unlist(lapply(1:ncats, function(i) tempsurvi$std.err[sum(tempsurvi$strata[1:i])]))))
                                                                       }
                                                                     })))
                                  colnames(observedall) <- paste(rep(c("O", "SE"),4), sort(rep(c(3,5,10,15),2)), sep="")
                                
                                  return(list(expected=expectedall, observed=observedall))
                                } else {
                                  mydatsplit          <- survSplit(Surv(TTE, FUTbc, bcdeathbc) ~ ., 
                                                                   data= mydat, cut=10,
                                                                   episode= "tg", 
                                                                   id="newid")
                                  if(length(with(mydatsplit, which(FUTbc-TTE < 0.001)))>0){ 
                                    mydatsplit          <- mydatsplit[-with(mydatsplit, which(FUTbc-TTE < 0.001)),]
                                    }
                                    
                                  observedall         <- as.data.frame(do.call(cbind, lapply(1:2, function(j){
                                                                       tempsurvi           <- survfit(Surv(TTE, FUTbc, bcdeathbc==1) ~ patcats, data=subset(mydatsplit, tg %in% 1:j))
                                                                       if(ncats==1){
                                                                         return(cbind(tail(tempsurvi$cumhaz,1),
                                                                                      tail(tempsurvi$std.err,1)))
                                                                       } else {
                                                                         return(cbind(unlist(lapply(1:ncats, function(i) tempsurvi$cumhaz[sum(tempsurvi$strata[1:i])])),
                                                                                      unlist(lapply(1:ncats, function(i) tempsurvi$std.err[sum(tempsurvi$strata[1:i])]))))
                                                                       }
                                                                     })))
                                  colnames(observedall) <- paste(rep(c("O", "SE"),2), sort(rep(c(10,15),2)), sep="")
                                
                                  return(list(expected=expectedall, observed=observedall))
                                
                                
                                }}

#######################################################################################################################################################################################################
#######################################################################################################################################################################################################

PREDcalipos        <- function(mydat, ncats){

                                if(ncats==4){
                                  mydat$patcats       <- cut(mydat$PREDICT, breaks=quantile(mydat$PREDICT, probs=c(c(0,0.2,0.5,0.8,1)), labels=c("low", "middle-low", "middle-high", "high")), include.lowest=TRUE) # PI quantiles
                                } else if(ncats==3){ 
                                  mydat$patcats       <- cut(mydat$PREDICT, breaks=quantile(mydat$PREDICT, probs=c(c(0,0.3,0.7,1)), labels=c("low", "middle", "high")), include.lowest=TRUE) # PI quantiles
                                } else {
                                  ncats               <- 1
                                  mydat$patcats       <- cut(mydat$PREDICT, breaks=quantile(mydat$PREDICT, probs=c(c(0,1)), labels=c("low")), include.lowest=TRUE)
                                }
# expected absolute risks
                                mydat$exp3y         <- (1-exp(-exp(erposbase(3) + mydat$PREDICT))) * exp(-exp(0.0698*((mydat$AgeDiagIndex/10)^2 - 34.234) + ocmbase(3)))
                                mydat$exp5y         <- (1-exp(-exp(erposbase(5) + mydat$PREDICT))) * exp(-exp(0.0698*((mydat$AgeDiagIndex/10)^2 - 34.234) + ocmbase(5)))
                                mydat$exp10y        <- (1-exp(-exp(erposbase(10) + mydat$PREDICT))) * exp(-exp(0.0698*((mydat$AgeDiagIndex/10)^2 - 34.234) + ocmbase(10)))
                                mydat$exp15y        <- (1-exp(-exp(erposbase(15) + mydat$PREDICT))) * exp(-exp(0.0698*((mydat$AgeDiagIndex/10)^2 - 34.234) + ocmbase(15)))

# return cumhaz in patcats
                                expected3y          <- unlist(lapply(1:ncats, function(i) with(subset(mydat, patcats==levels(patcats)[i]), mean(-log(1-exp3y)))))
                                expected5y          <- unlist(lapply(1:ncats, function(i) with(subset(mydat, patcats==levels(patcats)[i]), mean(-log(1-exp5y)))))
                                expected10y         <- unlist(lapply(1:ncats, function(i) with(subset(mydat, patcats==levels(patcats)[i]), mean(-log(1-exp10y)))))
                                expected15y         <- unlist(lapply(1:ncats, function(i) with(subset(mydat, patcats==levels(patcats)[i]), mean(-log(1-exp15y)))))

                                
                                expectedall         <- cbind(E3=as.data.frame(expected3y),
                                                             E5=as.data.frame(expected5y),
                                                             E10=as.data.frame(expected10y),
                                                             E15=as.data.frame(expected15y))

                                mydatsplit          <- survSplit(Surv(TTE, FUTbc, bcdeathbc) ~ ., 
                                                                 data= mydat, cut=c(3,5,10),
                                                                 episode= "tg", 
                                                                 id="newid")
                                if(length(with(mydatsplit, which(FUTbc-TTE < 0.001)))>0){ 
                                  mydatsplit          <- mydatsplit[-with(mydatsplit, which(FUTbc-TTE < 0.001)),]
                                  }
                                  
                                observedall         <- as.data.frame(do.call(cbind, lapply(1:4, function(j){
                                                                     tempsurvi           <- survfit(Surv(TTE, FUTbc, bcdeathbc==1) ~ patcats, data=subset(mydatsplit, tg %in% 1:j))
                                                                     if(ncats==1){
                                                                       return(cbind(tail(tempsurvi$cumhaz,1),
                                                                                    tail(tempsurvi$std.err,1)))
                                                                     } else {
                                                                       return(cbind(unlist(lapply(1:ncats, function(i) tempsurvi$cumhaz[sum(tempsurvi$strata[1:i])])),
                                                                                    unlist(lapply(1:ncats, function(i) tempsurvi$std.err[sum(tempsurvi$strata[1:i])]))))
                                                                     }
                                                                   })))
                                colnames(observedall) <- paste(rep(c("O", "SE"),4), sort(rep(c(3,5,10,15),2)), sep="")
                                
                                return(list(expected=expectedall, observed=observedall))
                                }




# pool expected and observed values
poolcali        <- function(caliobject, expname, obsname, sename) {
                     avgE            <- colMeans(do.call(rbind, lapply(1:m, function(i) caliobject[[i]]$expected[,expname])))

                     # pooled observed values
                     obsmatrix       <- do.call(rbind, lapply(1:m, function(i) caliobject[[i]]$observed[,obsname]))
                     avgQ            <- colMeans(obsmatrix)

                     # pooled variance
                     avgU            <- colMeans(do.call(rbind, lapply(1:m, function(i) caliobject[[i]]$observed[,sename]))^2)
                     if(ncol(obsmatrix)==1){
                       myB             <- sum(apply(obsmatrix, 1, function(r) (r-avgQ)^2))/(m-1)
                     } else {
                       myB             <- rowSums(apply(obsmatrix, 1, function(r) (r-avgQ)^2))/(m-1)
                     }
                     myT             <- avgU + (1+1/m)*myB
                     pooledSE        <- sqrt(myT)

                     poolupper       <- avgQ+1.96*pooledSE
                     poollower       <- avgQ-1.96*pooledSE
                     poollower[which(poollower<0)] <- 0

# This gives out pooled Pr(T<=t) = 1 - S(t)                     
                     myout           <- cbind(as.data.frame(1-exp(-avgE)),
                                              as.data.frame(1-exp(-avgQ)),
                                              as.data.frame(1-exp(-poolupper)),
                                              as.data.frame(1-exp(-poollower)))
                     colnames(myout) <- c("expected", "observed", "obsupper", "obslower")
                     
                     return(myout)

                     }

#######################################################################################################################################################################################################
# FUNCTIONS FOR RETURNING PREDICT CATEGORIES TO OBSERVED DATA 
#######################################################################################################################################################################################################
# VOTE the category
# Return to observed data
# Split to categories
# VOTING separates the categories better for survival

PREDICTtoOBS        <- function(mydat, mystudy, mutgene, ncats){

                                 totIDs       <- unique(unlist(lapply(mydat, function(l) subset(l, study %in% mystudy & BRCA12==mutgene, select="ID"))))
                                 countsIDs    <- unlist(lapply(totIDs, function(id) sum(unlist(lapply(mydat, function(l) id %in% l$ID)))))
                                 
                                 relIDs       <- totIDs[which(countsIDs>5)]
                                 relcounts    <- countsIDs[which(countsIDs>5)]
                                 
                                 PREDpooled   <- unlist(lapply(relIDs, function(id) mean(unlist(lapply(mydat, function(l) l$PREDICT[which(l$ID==id)])))))

                                 if(ncats==4){
                                   PREDcats     <- as.numeric(cut(PREDpooled, breaks=quantile(PREDpooled, probs=c(0,0.2,0.5, 0.8,1), labels=c("low", "middle-low", "middle-high", "high")), include.lowest=TRUE))
                                 } else if (ncats==3){
                                   PREDcats     <- as.numeric(cut(PREDpooled, breaks=quantile(PREDpooled, probs=c(0,0.3,0.7,1), labels=c("low", "middle", "high")), include.lowest=TRUE))
                                 }
                                 
                                 obsdata      <- merge(cbind(relIDs, PREDcats), mymid$data, by.x="relIDs", by.y="ID", all.x=TRUE, all.y=FALSE)
                                 
                                 return(obsdata)
                                 }


VOTEStoOBS          <- function(mydat, mystudy, mutgene, ncats){

                                 pdat        <- lapply(mydat, function(l){
                                                         l               <- subset(l, study %in% mystudy & BRCA12==mutgene)
                                                         if(ncats==4){
                                                           l$patcats     <- as.numeric(cut(l$PREDICT, breaks=quantile(l$PREDICT, probs=c(0,0.2,0.5,0.8,1), labels=c("low", "middle-low", "middle-high", "high")), include.lowest=TRUE)) 
                                                         } else if(ncats==3){
                                                           l$patcats     <- as.numeric(cut(l$PREDICT, breaks=quantile(l$PREDICT, probs=c(0,0.3,0.7,1), labels=c("low", "middle", "high")), include.lowest=TRUE)) 
                                                         }
                                                       return(l)
                                                  })

                                 totIDs       <- unique(unlist(lapply(pdat, function(l) l$ID)))
                                 countsIDs    <- unlist(lapply(totIDs, function(id) sum(unlist(lapply(pdat, function(l) id %in% l$ID)))))
                                 
                                 relIDs       <- totIDs[which(countsIDs>5)]
                                 relcounts    <- countsIDs[which(countsIDs>5)]
                                 
                                 PREDvotes    <- unlist(lapply(lapply(relIDs, function(id) unlist(lapply(pdat, function(l) as.numeric(l$patcat[which(l$ID==id)])))), function(votes) {
                                                                                                                          myvotes    <- sort(table(votes), decreasing=TRUE)
                                                                                                                            if(length(myvotes)==1){
                                                                                                                              return(names(myvotes)[1])
                                                                                                                            } else if(myvotes[1]==myvotes[2]){ # unidentified error with BCAC, here
                                                                                                                              return(sample(names(myvotes)[1:2], 1))
                                                                                                                            } else {
                                                                                                                              return(names(myvotes)[1])
                                                                                                                            }
                                                                                                                        }))

                                 obsdata      <- merge(cbind(relIDs, PREDvotes), mymid$data, by.x="relIDs", by.y="ID", all.x=TRUE, all.y=FALSE)
                                 
                                 return(obsdata)
                                 }

#######################################################################################################################################################################################################
#######################################################################################################################################################################################################
#######################################################################################################################################################################################################


