############################################################
#Purpose: Dose-response models for old coronavirus and SARS-CoV-2
#Final edit: 7 Apr 2022 (Netherlands) 9:00
#Editor: Fumi Miura
############################################################
###Procedure
#1. Run "1_DRmodel_CoronaDR.R"
#2. Bootstrap function
#3. Refitting 
#4. Draw trajectories of dose-response curves
#5. Estimated parameters with 95%CI (Table-S3)
#6. Bootstrapped dose-response curves
############################################################

###1. Run "1_DRmodel_CoronaDR.R" -----
#estimated parameters, likelihoods, dose-response models are required

###2. Bootstrap function -----
##Bootstrap function
Boot_func <- function(data, DRmodel, DRpara, N_boot=1000){
  Prob_est <- DRmodel(data$Dose, DRpara)
  Trial_num <- data$Total
  Boot_sample <- matrix(0,length(Trial_num),N_boot)
  for(i in 1:length(Trial_num)){
    Boot_sample[i,] <- rbinom(N_boot,Trial_num[i],Prob_est[i])
  }
  Boot_data <- matrix(0,length(Boot_sample[,1]),3)
  Boot_data[,1] <- data$Dose
  Boot_data[,2] <- data$Total
  
  return(list(Boot_sample = Boot_sample,
              Boot_data   = Boot_data))
}

###3. Refitting  -----
##Generate bootstrapped samples
Boot_delta_nonvac   <- Boot_func(data_all, DR_delta, est_para_list$delta_nonvac$est, 1000)
Boot_gamma_nonvac   <- Boot_func(data_all, DR_gamma, est_para_list$gamma_nonvac$est, 1000)
Boot_tweedie_nonvac   <- Boot_func(data_all, DR_tweedie, est_para_list$tweedie_nonvac$est, 1000)
Boot_twolev_nonvac   <- Boot_func(data_all, DR_twolevel, est_para_list$twolev_nonvac$est, 1000)
Boot_gamPM_nonvac   <- Boot_func(data_all, DR_gammaPointmass, est_para_list$gamPM_nonvac$est, 1000)
Boot_gamPM0_nonvac   <- Boot_func(data_all, DR_gammaPointmass, c(0, est_para_list$gamPM0_nonvac$est), 1000) #0 means "pointmass at exact 0"

##delta -----
delta_nonvac_refit <- list()
Boot_delta_nonvac_para <- matrix(0,1000,1)
for(i in 1:1000){
  Boot_delta_nonvac$Boot_data[,3] <- Boot_delta_nonvac$Boot_sample[,i]
  delta_nonvac_refit[i] <- optim( fn= LogLikelihood_delta_e(Boot_delta_nonvac$Boot_data), 
                                  par=ini_para_delta, control = list(fnscale = -1), hessian = T)
  Boot_delta_nonvac_para[i,] <- exp(delta_nonvac_refit[[i]])
}

##gamma -----
gamma_nonvac_refit <- list()
Boot_gamma_nonvac_para <- matrix(0,1000,2)
for(i in 1:1000){
  Boot_gamma_nonvac$Boot_data[,3] <- Boot_gamma_nonvac$Boot_sample[,i]
  gamma_nonvac_refit[i] <- optim( fn= LogLikelihood_gamma_e(Boot_gamma_nonvac$Boot_data), 
                                  par=ini_para_gamma, control = list(fnscale = -1), hessian = T)
  Boot_gamma_nonvac_para[i,] <- exp(gamma_nonvac_refit[[i]])
}

##tweedie -----
tweedie_nonvac_refit <- list()
Boot_tweedie_nonvac_para <- matrix(0,1000,3)
for(i in 1:1000){
  Boot_tweedie_nonvac$Boot_data[,3] <- Boot_tweedie_nonvac$Boot_sample[,i]
  tweedie_nonvac_refit[i] <- optim( fn= LogLikelihood_tweedie_e(Boot_tweedie_nonvac$Boot_data), 
                                  par=ini_para_twee, control = list(fnscale = -1), hessian = T)
  Boot_tweedie_nonvac_para[i,] <- exp(tweedie_nonvac_refit[[i]])
}

##two-level -----
twolev_nonvac_refit <- list()
Boot_twolev_nonvac_para <- matrix(0,1000,3)
for(i in 1:1000){
  Boot_twolev_nonvac$Boot_data[,3] <- Boot_twolev_nonvac$Boot_sample[,i]
  twolev_nonvac_refit[i] <- optim( fn= LogLikelihood_twolevel_logistic(Boot_twolev_nonvac$Boot_data), 
                                    par=ini_para_twolevel, control = list(fnscale = -1), hessian = T)
  Boot_twolev_nonvac_para[i,] <- c(exp(twolev_nonvac_refit[[i]][1:2]),1/(1+exp(-twolev_nonvac_refit[[i]][3])))
}

##gamma + PM -----
gamPM_nonvac_refit <- list()
Boot_gamPM_nonvac_para <- matrix(0,1000,4)
for(i in 1:1000){
  Boot_gamPM_nonvac$Boot_data[,3] <- Boot_gamPM_nonvac$Boot_sample[,i]
  gamPM_nonvac_refit[i] <- optim( fn= LogLikelihood_gammaPointmass_logistic_e(Boot_gamPM_nonvac$Boot_data), 
                                   par=ini_para_gamPM, control = list(fnscale = -1), hessian = T)
  Boot_gamPM_nonvac_para[i,] <- c(exp(gamPM_nonvac_refit[[i]][1:3]), 1/(1+exp(-gamPM_nonvac_refit[[i]][4])))
}

##gamma + PM0 -----
gamPM0_nonvac_refit <- list()
Boot_gamPM0_nonvac_para <- matrix(0,1000,4)
for(i in 1:1000){
  Boot_gamPM0_nonvac$Boot_data[,3] <- Boot_gamPM0_nonvac$Boot_sample[,i]
  gamPM0_nonvac_refit[i] <- optim( fn= LogLikelihood_gammaPointmass_logistic_0immu_e(Boot_gamPM0_nonvac$Boot_data), 
                                  par=ini_para_gamPM0, control = list(fnscale = -1), hessian = T)
  Boot_gamPM0_nonvac_para[i,] <- c(0, exp(gamPM0_nonvac_refit[[i]][1:2]), 1/(1+exp(-gamPM0_nonvac_refit[[i]][3])))
}

###4. Draw trajectories of dose-response curves -----
traject_func <- function(BootPara, DRmodel, x_plot){
  inter_bt <- matrix(0, length(BootPara[,1]), length(x_plot))
  for(i in 1:length(BootPara[,1])){
    inter_bt[i,] <- DRmodel(x_plot, BootPara[i,])
  }
  percentile_bt <- matrix(0, 3, length(x_plot))
  for(i in 1:length(x_plot)){
    percentile_bt[,i] <- quantile(inter_bt[,i], c(0.025, 0.5, 0.975))
  }
  return(list(inter_bt = inter_bt,
              percentile_bt = percentile_bt))
}
#draw trajectories
x_plot <- 10^seq(0,6,by=0.05)
delta_nonvac_traject <- traject_func(Boot_delta_nonvac_para, DR_delta, x_plot) #note: x_plot should be matched with x3 in "3_ggplot_CoronaDR.R". 
gamma_nonvac_traject <- traject_func(Boot_gamma_nonvac_para, DR_gamma, x_plot) 
tweedie_nonvac_traject <- traject_func(Boot_tweedie_nonvac_para, DR_tweedie, x_plot) 
twolev_nonvac_traject <- traject_func(Boot_twolev_nonvac_para, DR_twolevel, x_plot) 
gamPM_nonvac_traject <- traject_func(Boot_gamPM_nonvac_para, DR_gammaPointmass, x_plot) 
gamPM0_nonvac_traject <- traject_func(Boot_gamPM0_nonvac_para, DR_gammaPointmass, x_plot) 

###5. Estimated parameters with 95%CI (Table-S3) -----
#delta_nonvac
est_list_withCI <- data.frame(
delta_a = c(est_para_list$delta_nonvac$est, quantile(Boot_delta_nonvac_para[,1], c(0.025, 0.975))),
#gamma_nonvac
gamma_theta = c(est_para_list$gamma_nonvac$est[1], quantile(Boot_gamma_nonvac_para[,1], c(0.025, 0.975))),
gamma_k = c(est_para_list$gamma_nonvac$est[2], quantile(Boot_gamma_nonvac_para[,2], c(0.025, 0.975))),
#tweedie_nonvac
tweedie_lambda = c(est_para_list$tweedie_nonvac$est[1], quantile(Boot_tweedie_nonvac_para[,1], c(0.025, 0.975))),
tweedie_gamma = c(est_para_list$tweedie_nonvac$est[2], quantile(Boot_tweedie_nonvac_para[,2], c(0.025, 0.975))),
tweedie_alpha = c(est_para_list$tweedie_nonvac$est[3], quantile(Boot_tweedie_nonvac_para[,3], c(0.025, 0.975))),
#twolev_nonvac
twolev_a1 = c(est_para_list$twolev_nonvac$est[1], quantile(Boot_twolev_nonvac_para[,1], c(0.025, 0.975))),
twolev_a2 = c(est_para_list$twolev_nonvac$est[2], quantile(Boot_twolev_nonvac_para[,2], c(0.025, 0.975))),
twolev_p1 = c(est_para_list$twolev_nonvac$est[3], quantile(Boot_twolev_nonvac_para[,3], c(0.025, 0.975))),
#gamPM_nonvac
gamPM_a1 = c(est_para_list$gamPM_nonvac$est[1], quantile(Boot_gamPM_nonvac_para[,1], c(0.025, 0.975))),
gamPM_theta = c(est_para_list$gamPM_nonvac$est[2], quantile(Boot_gamPM_nonvac_para[,2], c(0.025, 0.975))),
gamPM_k = c(est_para_list$gamPM_nonvac$est[3], quantile(Boot_gamPM_nonvac_para[,3], c(0.025, 0.975))),
gamPM_p1 = c(est_para_list$gamPM_nonvac$est[4], quantile(Boot_gamPM_nonvac_para[,4], c(0.025, 0.975))),
#gamPM0_nonvac
gamPM0_a1 = c(0,                                  quantile(Boot_gamPM0_nonvac_para[,1], c(0.025, 0.975))), #fixed as 0
gamPM0_theta = c(est_para_list$gamPM0_nonvac$est[1], quantile(Boot_gamPM0_nonvac_para[,2], c(0.025, 0.975))),
gamPM0_k = c(est_para_list$gamPM0_nonvac$est[2], quantile(Boot_gamPM0_nonvac_para[,3], c(0.025, 0.975))),
gamPM0_p1 = c(est_para_list$gamPM0_nonvac$est[3], quantile(Boot_gamPM0_nonvac_para[,4], c(0.025, 0.975)))
)


###6. Bootstrapped dose-response curves -----
#intersection point (Gamma dose-reseponse model)
int_gamma <- function(prob, dose, theta){
  function(k){
    prob - DR_gamma(dose,c(theta,k))
  } 
}
#intersection point (Dirac delta dose-reseponse model)
int_delta <- function(prob, dose){
  function(a){
    prob - DR_delta(dose,a)
  } 
}
#generate bootstrap samples and determine the intersection point (Gamma DR model)
boot_int_gamma_func <- function(theta_fix=0, boot_N=100){
  bootSample <- rbinom(boot_N,data_covid$Total,data_covid$Infected/data_covid$Total)
  para_gam_boot <- c(0,0)
  for (n in bootSample) {
    f<-int_gamma(n/data_covid$Total,
                 data_covid$Dose,
                 theta_fix)
    para_gam_boot <- rbind(para_gam_boot, c(theta_fix, uniroot(f, c(0, 10^6))$root))
  }
  return(para_gam_boot[-1,])
}
#generate bootstrap samples and determine the intersection point (Delta DR model)
boot_int_delta_func <- function(boot_N=100){
  bootSample <- rbinom(boot_N,data_covid$Total,data_covid$Infected/data_covid$Total)
  para_delta_boot <- c(0)
  for (n in bootSample) {
    f<-int_delta(n/data_covid$Total,
                 data_covid$Dose)
    para_delta_boot <- rbind(para_delta_boot, uniroot(f, c(0, 10^6))$root)
  }
para_delta_boot <- para_delta_boot[-1,]
  return(para_delta_boot[-1])
}

###generate bootstrapped curves
#x-axis (dose)
x_plot <- 10^seq(0,6,by=0.05)
x_plot_lower <- 10^seq(-2,2,by=0.05)
##parameter sets, varying var(s) = theta
para_gam <- c(0,0)
for (theta in 10^(-6:2)) {
  f<-int_gamma(data_covid$Infected/data_covid$Total,
               data_covid$Dose,
               theta)
  para_gam <- rbind(para_gam, c(theta, uniroot(f, c(0, 10^6))$root))
}
para_gam <- para_gam[-1,]

##ggplot_data, varying var(s) = theta
#normal dose range
sim_dr_gamma_data <- data.frame(
  x_plot = rep(x_plot, length(para_gam[,1])),
  prob = as.numeric(apply(X=para_gam, FUN=function(x){DR_gamma(x_plot,c(x[1], x[2]))}, MARGIN = 1)),
  theta = as.numeric(apply(X=para_gam, FUN=function(x){rep(x[1],length(x_plot))}, MARGIN = 1))
)
sim_dr_gamma_data <- as_tibble(sim_dr_gamma_data)
#lower dose range
sim_dr_gamma_data_lower <- data.frame(
  x_plot = rep(x_plot_lower, length(para_gam[,1])),
  prob = as.numeric(apply(X=para_gam, FUN=function(x){DR_gamma(x_plot_lower,c(x[1], x[2]))}, MARGIN = 1)),
  theta = as.numeric(apply(X=para_gam, FUN=function(x){rep(x[1],length(x_plot_lower))}, MARGIN = 1))
)

##ggplot_data, animal dose-response model (SARS-CoV-1)
#normal dose range
animal_DR_delta_data <- data.frame(
  x_plot = x_plot,
  prob = DR_delta(x_plot,(1/410)*(1/0.7)) #Watanabe (2010) Risk Analysis + the ratio of PFU to TCID50 has been established as 0.7 (Covés-Datson et al. 2020)
)
#lower dose range
animal_DR_delta_data_lower <- data.frame(
  x_plot = x_plot_lower,
  prob = DR_delta(x_plot_lower,(1/410)*(1/0.7))
)

##ggplot_data, bootstrapped Gamma dose-response curves
para_boot_gam_var100 <- boot_int_gamma_func(theta_fix=100, boot_N=100)
boot_plot_gam_var100 <- apply(X=para_boot_gam_var100, FUN=function(x){DR_gamma(x_plot,c(x[1], x[2]))}, MARGIN = 1)
boot_gam_var100_data <- data.frame(
  x_plot = x_plot,
  prob = apply(X=boot_plot_gam_var100,function(x){quantile(x,0.5)},MARGIN = 1),
  prob_upp = apply(X=boot_plot_gam_var100,function(x){quantile(x,0.975)},MARGIN = 1),
  prob_low = apply(X=boot_plot_gam_var100,function(x){quantile(x,0.025)},MARGIN = 1)
)

para_boot_gam_var10 <- boot_int_gamma_func(theta_fix=10, boot_N=100)
boot_plot_gam_var10 <- apply(X=para_boot_gam_var10, FUN=function(x){DR_gamma(x_plot,c(x[1], x[2]))}, MARGIN = 1)
boot_gam_var10_data <- data.frame(
  x_plot = x_plot,
  prob = apply(X=boot_plot_gam_var10,function(x){quantile(x,0.5)},MARGIN = 1),
  prob_upp = apply(X=boot_plot_gam_var10,function(x){quantile(x,0.975)},MARGIN = 1),
  prob_low = apply(X=boot_plot_gam_var10,function(x){quantile(x,0.025)},MARGIN = 1)
)

para_boot_gam_var1 <- boot_int_gamma_func(theta_fix=1, boot_N=100)
boot_plot_gam_var1 <- apply(X=para_boot_gam_var1, FUN=function(x){DR_gamma(x_plot,c(x[1], x[2]))}, MARGIN = 1)
boot_gam_var1_data <- data.frame(
  x_plot = x_plot,
  prob = apply(X=boot_plot_gam_var1,function(x){quantile(x,0.5)},MARGIN = 1),
  prob_upp = apply(X=boot_plot_gam_var1,function(x){quantile(x,0.975)},MARGIN = 1),
  prob_low = apply(X=boot_plot_gam_var1,function(x){quantile(x,0.025)},MARGIN = 1)
)

para_boot_gam_var0.1 <- boot_int_gamma_func(theta_fix=0.1, boot_N=100)
boot_plot_gam_var0.1 <- apply(X=para_boot_gam_var0.1, FUN=function(x){DR_gamma(x_plot,c(x[1], x[2]))}, MARGIN = 1)
boot_gam_var0.1_data <- data.frame(
  x_plot = x_plot,
  prob = apply(X=boot_plot_gam_var0.1,function(x){quantile(x,0.5)},MARGIN = 1),
  prob_upp = apply(X=boot_plot_gam_var0.1,function(x){quantile(x,0.975)},MARGIN = 1),
  prob_low = apply(X=boot_plot_gam_var0.1,function(x){quantile(x,0.025)},MARGIN = 1)
)

para_boot_gam_var0.01 <- boot_int_gamma_func(theta_fix=0.01, boot_N=100)
boot_plot_gam_var0.01 <- apply(X=para_boot_gam_var0.01, FUN=function(x){DR_gamma(x_plot,c(x[1], x[2]))}, MARGIN = 1)
boot_gam_var0.01_data <- data.frame(
  x_plot = x_plot,
  prob = apply(X=boot_plot_gam_var0.01,function(x){quantile(x,0.5)},MARGIN = 1),
  prob_upp = apply(X=boot_plot_gam_var0.01,function(x){quantile(x,0.975)},MARGIN = 1),
  prob_low = apply(X=boot_plot_gam_var0.01,function(x){quantile(x,0.025)},MARGIN = 1)
)

para_boot_gam_var0.001 <- boot_int_gamma_func(theta_fix=0.001, boot_N=100)
boot_plot_gam_var0.001 <- apply(X=para_boot_gam_var0.001, FUN=function(x){DR_gamma(x_plot,c(x[1], x[2]))}, MARGIN = 1)
boot_gam_var0.001_data <- data.frame(
  x_plot = x_plot,
  prob = apply(X=boot_plot_gam_var0.001,function(x){quantile(x,0.5)},MARGIN = 1),
  prob_upp = apply(X=boot_plot_gam_var0.001,function(x){quantile(x,0.975)},MARGIN = 1),
  prob_low = apply(X=boot_plot_gam_var0.001,function(x){quantile(x,0.025)},MARGIN = 1)
)

##ggplot_data, bootstrapped Delta dose-response curves
boot_delta_plot <- sapply(X=boot_int_delta_func(boot_N=100), FUN=function(x){DR_delta(x_plot,x)})#boot_int_delta_func(boot_N=100)#
boot_DR_delta_data <- data.frame(
  x_plot = x_plot,
  prob = apply(X=boot_delta_plot,function(x){quantile(x,0.5)},MARGIN = 1),
  prob_upp = apply(X=boot_delta_plot,function(x){quantile(x,0.975)},MARGIN = 1),
  prob_low = apply(X=boot_delta_plot,function(x){quantile(x,0.025)},MARGIN = 1)
)
