#!/usr/bin/env Rscript

# R-4.0.0

## Packages
library(survival)
library(mice)
library(dplyr)

load("sparse.data.RData")

# Input is a "sparsedata" data.frame with the following columns. The data is not public.
# ID               Patient ID                     character                              No missing values
# study            Study name                     character                              No missing values
# Country          Study country                  character                              No missing values
# famHist          Family history of BC/OC        dichotomous                            Missing values
# AgeDiag          Age of diagnosis               numeric                                No missing values
# YearDiag         Year of diagnosis              numeric                                No missing values       
# TTE              Years from dg to study entry   numeric                                Missing values (cases with missing data excluded after imputation before survival analyses)        
# FUT              Follow-up time in years        numeric                                Missing values (cases with missing data excluded after imputation before survival analyses)        
# VTS              Vital status                   dichotomous                            Missing values (cases with missing data excluded after imputation before survival analyses)        
# VTSbc            VTS censored for 2nd BC        factor(alive, BC-death, other-death)   Missing values (cases with missing data excluded after imputation before survival analyses)          
# FUTbc            FUT censored for 2nd BC        numeric                                Missing values (cases with missing data excluded after imputation before survival analyses)        
# BRCA12           BRCA1 or BRCA2 mutation        factor(BRCA1, BRCA2, not, filter)      No missing values
# ER               Estrogen receptor status       dichotomous                            Missing values
# PR               Progesterone receptor status   dichotomous                            Missing values
# GRD              Tumor grade                    factor(1,2,3)                          Missing values
# stage            TNM stage                      factor(0,1,2,3,4)                      Missing values
# morph            Tumor histological type        factor                                 Missing values
# logsize          Tumor size (mm) in logscale    numeric                                Missing values
# T                Tumor size category            factor(1,2,3)                          Missing values
# M                Metastasis at diagnosis        dichotomous                            Missing values
# nodecount        Number of affected lymph nodes numeric                                Missing values
# N                Lymph node status              dichotomous                            Missing values
# HER2             HER2 status                    dichotomous                            Missing values
# chemo            Adjuvant chemotherapy          dichotomous                            Missing values
# anthrac          Anthracycline-based chemo      dichotomous                            Missing values
# tax              Taxane-based chemo             dichotomous                            Missing values
# cmf              CMF chemo                      dichotomous                            Missing values       
# endocr           Adjuvant endocrine therapy     dichotomous                            Missing values
# tamox            Tamoxifen therapy              dichotomous                            Missing values
# ai               Aromatase inhibitor therapy    dichotomous                            Missing values
# trastu           Trastuzumab therapy            dichotomous                            Missing values
# surgery          Surgery type                   factor                                 Missing values
# radiation        Radiation therapy              dichotomous                            Missing values  
# chemoneo         Neoadjuvant chemotherapy       dichotomous                            Missing values
# anthneo          Anthracycline-based chemoneo   dichotomous                            Missing values
# taxneo           Taxane-based chemoneo          dichotomous                            Missing values
# cmfneo           CMF chemoneo                   dichotomous                            Missing values       
# bcdeath          BC-specific VTS                factor(alive, bc, other)               Missing values     
# bcdeathbc        BC-specific VTS cens for 2 BC  factor(alive, bc, other)               Missing values     

##############
# PREPROCESS #
##############
 
# Nelson-Aalen estimators of baseline hazard of death and bc-specific death
nelsonaalen_delayed       <- function (data, timevar_start, timevar_stop, statusvar, critvar, critoper, critval){
  if (!is.data.frame(data)) stop("Data must be a data frame")
  criterion               <- substitute(FUN(data[[column]], values), list(FUN=as.name(critoper), column=critvar, values=critval))
  timevar_start           <- as.character(substitute(timevar_start))
  timevar_stop            <- as.character(substitute(timevar_stop))
  statusvar               <- as.character(substitute(statusvar))
  survdat                 <- subset(data, eval(criterion))
  time_stop               <- survdat[, timevar_stop]
  time_start              <- survdat[, timevar_start]
  status                  <- as.numeric(survdat[, statusvar]==1)
  hazard                  <- basehaz(coxph(Surv(time_start,time_stop, status) ~ 1, data = survdat))
  idx                     <- match(data[,timevar_stop], hazard[, "time"])
  return(hazard[idx, "hazard"])
}

# bc-specific death was reported by a subset of studies
# baseline hazard of bc-specific death was based on these studies only
bcdeathstudies <- unique(sparsedata$study)[which(unlist(lapply(unique(sparsedata$study), 
                                                             function(s) {
                                                               length(which(sparsedata$study==s & sparsedata$bcdeath %in% 1:2))/length(which(sparsedata$study==s & sparsedata$VTS==1))>0.8
                                                             })))]

sparsedata$H.tilde_os       <- nelsonaalen_delayed(sparsedata, "TTE", "FUT", "VTS", "VTS", "%in%", 0:1)
sparsedata$H.tilde_bc       <- nelsonaalen_delayed(sparsedata, "TTE", "FUT", "bcdeath", "study", "%in%", bcdeathstudies)


######################
# PREPARE IMPUTATION # 
######################

# Define the variants to be used in imputation
predictors          <- c("study", "famHist", "AgeDiag", "YearDiag", "TTE", "VTS", "BRCA12", 
                         "ER", "PR", "GRD", "morph", "logsize", "T", "M", "nodecount", "N", "HER2",               
                         "chemo", "anthrac", "tax", "cmf", "endocr", "tamox", "ai", "trastu", "surgery", "radiation",          
                         "chemoneo", "anthneo", "taxneo", "cmfneo", "bcdeath", "H.tilde_os", "H.tilde_bc") 
notpred             <- colnames(sparsedata)[which(!(colnames(sparsedata) %in% predictors))]

# Create mice output structures that can be customized and used in the actual run  
ini                 <- mice(sparsedata, maxit=0, vis="monotone", print=FALSE)

# Set imputation methods
# N.B. M (metastasis) status was missing for most of the patients, only few had positive status, and proper predictors were not available in the data.
# PREDICT not allow M-positive patients. Sampling the M status was chosen as the conservative approach for patients with missing M status.
meth                <- ini$meth
meth[notpred]       <- ""
meth["H.tilde_os"]  <- "pmm"
meth["M"]           <- "sample"

# Create a prediction matrix
pred                <-quickpred(sparsedata, 
                                exclude=notpred, 
                                include = c("H.tilde_os", "H.tilde_bc", "study"), 
                                minpuc = 0.2, 
                                mincor = 0.125,
                                method="spearman")

# Exclude the unnecessary
pred[notpred,]                                                                       <- 0

# Edit dependencies
pred["logsize",]                                                                     <- 0
pred["logsize", "T"]                                                                 <- 1 
pred[, "logsize"]                                                                    <- 0
pred["ER", c("PR", "morph")]                                                         <- 1
pred[c("PR", "morph"), "ER"]                                                         <- 1

pred[c("anthrac","tax","cmf", "tamox", "ai", "anthneo", "taxneo", "cmfneo"),]        <- 0
pred[c("anthrac","tax","cmf"), c("chemo", "YearDiag", "AgeDiag", "study")]           <- 1
pred[c("tamox", "ai"), c("endocr", "YearDiag", "AgeDiag", "study")]                  <- 1
pred[c("anthneo","taxneo","cmfneo"), c("chemoneo", "YearDiag", "AgeDiag", "study")]  <- 1

pred[c("chemo"), c("anthrac","tax","cmf")]                                           <- 0
pred[c("endocr"), c("tamox", "ai")]                                                  <- 0
pred[c("chemoneo"), c("anthneo","taxneo","cmfneo")]                                  <- 0
pred["HER2", c("YearDiag", "TTE")]                                                   <- 0
pred["trastu", ]                                                                     <- 0
pred["trastu", c("YearDiag", "HER2", "study")]                                       <- 1

pred["bcdeath", "VTS"]                                                               <- 1   
pred[, "N"]                                                                          <- 0
pred["nodecount", "N"]                                                               <- 1


# Modify visit sequence, to exclude the unwanted variables
vis                 <- ini$visitSequence
vis                 <- vis[-which(vis %in% notpred)]

# Post-processing for treatment groups
post                <-ini$post # customize the structure created in the dry run
post["trastu"]      <- "imp[[j]][data$YearDiag[!r[,j]]<1998 | data$HER2[!r[,j]]==0,i] <- 0"
post["anthrac"]     <- "imp[[j]][data$chemo[!r[,j]]==0,i]           <- 0"
post["tax"]         <- "imp[[j]][data$chemo[!r[,j]]==0,i]           <- 0"
post["cmf"]         <- "imp[[j]][data$chemo[!r[,j]]==0,i]           <- 0"
post["ai"]          <- "imp[[j]][data$endocr[!r[,j]]==0,i]          <- 0"
post["tamox"]       <- "imp[[j]][data$endocr[!r[,j]]==0,i]          <- 0"
post["anthneo"]     <- "imp[[j]][data$chemoneo[!r[,j]]==0,i]        <- 0"
post["taxneo"]      <- "imp[[j]][data$chemoneo[!r[,j]]==0,i]        <- 0"
post["cmfneo"]      <- "imp[[j]][data$chemoneo[!r[,j]]==0,i]        <- 0"
post["nodecount"]   <- "imp[[j]][data$N[!r[,j]]==0,i]               <- 0"
post["N"]           <- "imp[[j]][data$nodecount[!r[,j]]==0,i]       <- 0"


##########
# IMPUTE #
##########

# Impute, first iteration
impute               <- mice(data=sparsedata, m=50, seed=seed, meth=meth, pred=pred, visitSequence=vis, post=post, maxit=1)

# Build on the first iteration and continue until the 30th
count                <- 1

repeat{

  if(count==30) break
  
  temp               <- mice.mids(impute, maxit=1)
  impute             <- temp
  count              <- count+1
  
  saveRDS(impute, file="imputed.data.Rds")
}


