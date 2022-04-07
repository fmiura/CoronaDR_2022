############################################################
#Purpose: Dose-response models for old coronavirus and SARS-CoV-2
#Final edit: 7 Apr 2022 (Netherlands) 9:00
#Editor: Fumi Miura
############################################################
###Procedure
#0. install packages 
#1. Data table
#2. Dose-Response model
#3. Likelihood 
#4. MLE
#5. Summary of model fits
#6. Estimated parameters
############################################################

###0. install packages -----
library("tidyverse")
library("fAsianOptions")
library("DescTools")
library("patchwork")

###1. Data table ----- 
##endemic coronaviruses ----
#raw data
data <- read_csv("Data_HC_oldcoronaviruses.csv")
#All 
data_all <- data 
##SARS-CoV-2 ----
data_covid <- data[1,]
data_covid$Dose <- 10
data_covid$Total <- 34
data_covid$Infected <- 18
data_covid$Ill <- 18
data_covid$Virus <- "SARS-CoV-2"
data_covid$Reference <- "[Killingley 2022]"
data_covid_binCI <- BinomCI(x = data_covid$Infected, n = data_covid$Total, method = "jeffreys")

covid_raw_data.frame <- data.frame(
  Dose = 10,
  Total = 34,
  Infected = 18,
  Ill = 18, #need to be checked
  Virus = "SARS-CoV-2",
  prob = data_covid_binCI[1],
  prob_low = data_covid_binCI[2],
  prob_upp = data_covid_binCI[3]
)
#---

######2. DR model -----
#delta
DR_delta <- function(qd, para){
  a <- para[1]
  return(1-exp(-a*qd))
  } 
#gamma 
DR_gamma <- function(qd, para){
  theta <- para[1]
  k <- para[2]
  return(1-(1+theta*qd)^(-k)) 
  } 
#tweedie
DR_tweedie <- function(qd, para){  
  lambda <- para[1] 
  gamma <- para[2] 
  alpha <- para[3] 
  return(1 - exp(lambda*((1+gamma*qd)^(-alpha) - 1)))
  } 
#two-level 
DR_twolevel <- function(qd, para){
  a1 <- para[1]
  a2 <- para[2]
  p1 <- para[3]
  return(1- (p1*exp(-a1*qd) + (1-p1)*exp(-a2*qd))) 
  } 
#Gamma+pointmass
DR_gammaPointmass <- function(qd, para){
  a1 <- para[1]
  theta <- para[2]
  k <- para[3]
  p1 <- para[4]
  return(1- (p1*exp(-a1*qd) + (1-p1)*(1+theta*qd)^(-k)))
  } 

######3. Likelihood -----
#delta
LogLikelihood_delta_e <- function(data){
  LL <- function(para){ #para[1]=a
    a <- exp(para[1])  
    sum( (-1)*(data[,2]-data[,3])*a*data[,1] + data[,3]*log( 1 - exp(-a*data[,1])))
    }
  return(LL)
} 
#gamma
LogLikelihood_gamma_e <- function(data){ 
  LL <- function(para){ #exp(para[1])=theta, exp(para[2])=k
    theta <- exp(para[1])
    k <- exp(para[2])
    sum( (-1)*(data[,2]-data[,3])*k*log(1 + theta*data[,1]) + data[,3]*log(1 - (1 + theta*data[,1])^(-k)) )
  }
  return(LL)
} 
#tweedie
LogLikelihood_tweedie_e <- function(data){
  LL <- function(para){
    lambda <- exp(para[1]) 
    gamma <- exp(para[2]) 
    alpha <- exp(para[3]) 
    sum( (data[,2]-data[,3])*(  lambda*((1+gamma*data[,1])^(-alpha) - 1) ) 
         + data[,3]*log(1 - exp(lambda*((1+gamma*data[,1])^(-alpha) - 1))) )
  }
  return(LL)
}
#two level
LogLikelihood_twolevel_logistic <- function(data){
  LL <- function(para){
    a1 <- exp(para[1])
    a2 <- exp(para[2])
    p1 <- 1/(1+exp(-para[3]))
    sum( (data[,2]-data[,3])*log((p1*exp(-a1*data[,1]) + (1-p1)*exp(-a2*data[,1])) ) + data[,3]*log(1- ( p1*exp(-a1*data[,1]) + (1-p1)*exp(-a2*data[,1])) ) )
  }
  return(LL)
}
#Gamma+pointmass
LogLikelihood_gammaPointmass_logistic_e <- function(data){
  LL <- function(para){ ## a1 = para[1], theta = para[2], k = para[3], p1= 0<para[4]<1 
    a1 <- exp(para[1])
    theta <- exp(para[2])
    k <- exp(para[3])
    p1 <- 1/(1+exp(-para[4]))
    sum( (data[,2]-data[,3])*log((p1*exp(-a1*data[,1]) + (1-p1)*(1+theta*data[,1])^(-k))) + data[,3]*log(1- (p1*exp(-a1*data[,1]) + (1-p1)*(1+theta*data[,1])^(-k))) )
  }
  return(LL)
}
#Gamma+pointmass_0immu
LogLikelihood_gammaPointmass_logistic_0immu_e <- function(data){
  LL <- function(para){ ## a1 = 0, theta = para[1], k = para[2], p1= 0<para[3]<1 
    theta <- exp(para[1])
    k <- exp(para[2])
    p1 <- 1/(1+exp(-para[3]))
    sum( (data[,2]-data[,3])*log((p1 + (1-p1)*(1+theta*data[,1])^(-k))) + data[,3]*log(1- (p1 + (1-p1)*(1+theta*data[,1])^(-k))) )
  }
  return(LL)
}
######4. MLE -----
#delta
ini_para_delta <- c(-8)
opt.delta_e <- list(
  optim( fn= LogLikelihood_delta_e(data_all), par=ini_para_delta, control = list(fnscale = -1), hessian = T)
)
#gamma
ini_para_gamma <- c(-1.5,-1.2) 
opt.gamma_e <- list(
  optim( fn= LogLikelihood_gamma_e(data_all), par=ini_para_gamma, control = list(fnscale = -1), hessian = T)
)
#tweedie
ini_para_twee <- c(log(10), log(5*10^(-5)), log(10))
opt.tweedie_e <- list(
  optim( fn= LogLikelihood_tweedie_e(data_all), par=ini_para_twee, control = list(fnscale = -1), hessian = T)
)
#two level
ini_para_twolevel <- c(log(10^(-8)),log(10^(-6)),-1)
opt.twolevel_logistic <- list(
  optim( fn= LogLikelihood_twolevel_logistic(data_all), par=ini_para_twolevel, control = list(fnscale = -1), hessian = T)
)
#Gamma+pointmass
ini_para_gamPM <- c(log(8*10^(-20)),log(5*10^(-6)), log(0.8), 0) #c(log(8*10^(-9)),log(5*10^(-6)), log(0.8), 0)
opt.gammaPointmass_logistic <- list(
  optim( fn= LogLikelihood_gammaPointmass_logistic_e(data_all), par=ini_para_gamPM, control = list(fnscale = -1), hessian = T)
)
#Gamma+pointmass_0immu
ini_para_gamPM0 <- c(log(5*10^(-6)), log(0.8), 0) #c(log(8*10^(-9)),log(5*10^(-6)), log(0.8), 0)
opt.gammaPointmass_logistic_0immu <- list(
  optim( fn= LogLikelihood_gammaPointmass_logistic_0immu_e(data_all), par=ini_para_gamPM0, control = list(fnscale = -1), hessian = T)
)

######5. Summary of model fits -----
###log-likelihood and AICs
#non vaccinated, general PDFs
loglike_est <- c(opt.delta_e[[1]]$value, opt.gamma_e[[1]]$value, opt.tweedie_e[[1]]$value, opt.twolevel_logistic[[1]]$value, opt.gammaPointmass_logistic[[1]]$value, opt.gammaPointmass_logistic_0immu[[1]]$value) 
num_para <- c(1,2,3,3,4,3)
AIC_est <- 2*num_para - 2*loglike_est
Sum_modelfit <- as.data.frame(rbind(loglike_est, num_para, AIC_est))
names(Sum_modelfit) <- c("delta", "gamma", "tweedie", "two-level", "gamma+point-mass", "gamma+point-mass_imm0")

#Table S2, summary of model fit
#write.csv(Sum_modelfit, "TableS2_DRmodel_OldCoronaSARSCoV2.csv")

######6. Estimated parameters -----
est_para_list <- list(
  #non vaccinated
  delta_nonvac   = list(est=exp(opt.delta_e[[1]]$par)),
  gamma_nonvac   = list(est=exp(opt.gamma_e[[1]]$par)),
  tweedie_nonvac = list(est=exp(opt.tweedie_e[[1]]$par)),
  twolev_nonvac  = list(est=c(exp(opt.twolevel_logistic[[1]]$par[1:2]),1/(1+exp(-opt.twolevel_logistic[[1]]$par[3])))),
  gamPM_nonvac   = list(est=c(exp(opt.gammaPointmass_logistic[[1]]$par[1:3]), 1/(1+exp(-opt.gammaPointmass_logistic[[1]]$par[4])))),
  gamPM0_nonvac  = list(est=c(exp(opt.gammaPointmass_logistic_0immu[[1]]$par[1:2]), 1/(1+exp(-opt.gammaPointmass_logistic_0immu[[1]]$par[3]))))
)

#write rds 
#write_rds(est_para_list, "est_para_list")