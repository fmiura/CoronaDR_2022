############################################################
#Purpose: Dose-response models for old coronavirus and SARS-CoV-2
#Final edit: 7 Apr 2022 (Netherlands) 9:00
#Editor: Fumi Miura
############################################################
###Procedure
#0. Packages 
#1. SARS-CoV-2 dose-response models (Figure 1)
#2. Bootstrapped SARS-CoV-2 dose-reseponse curves (Figure S2)
#3. Estimated old coronavirus models (Figure S1)
############################################################

#Run "1_DRmodel_CoronaDR.R"
#Run "2_Bootstrap_CoronaDR.R"

###0. Packages -----
library(ggplot2)
library(patchwork)

###1.SARS-CoV-2 dose-response models (Figure 1) -----
#raw data
data_plot <- dplyr::mutate(data_all, Prob = Infected/Total, LogDose = log10(Dose))
data_plot <- as.data.frame(data_plot)
data_plot <- dplyr::mutate(data_plot,
                           Prob_low = as.numeric(BinomCI(data_plot$Infected, data_plot$Total)[,2]),
                           Prob_upp = as.numeric(BinomCI(data_plot$Infected, data_plot$Total)[,3]))
#ggplot
fig1A <- ggplot()+
  geom_line(data = sim_dr_gamma_data,
            aes(x=log10(x_plot), y=prob, group=theta, color=log10(theta)))+
  geom_point(data = covid_raw_data.frame,
             aes(x=log10(Dose), y=prob), color="red", size=2)+
  labs(color=expression(paste({Log[10]},"(var)", sep="")))+
  labs(x=expression(paste({Log[10]},"(dose)", sep="")), y = "Probability of infection")+
  theme_bw(base_size = 16)+
  geom_line(data = animal_DR_delta_data, 
            aes(x=log10(x_plot), y=prob), color=c("darkgrey"), linetype="dashed")

fig1B <- ggplot()+
  geom_line(data = sim_dr_gamma_data_lower,
            aes(x=log10(x_plot), y=prob, group=theta, color=log10(theta)))+
  geom_point(data = covid_raw_data.frame,
             aes(x=log10(Dose), y=prob), color="red", size=2)+
  labs(color=expression(paste({Log[10]},"(var)", sep="")))+
  labs(x=expression(paste({Log[10]},"(dose)", sep="")), y = "Probability of infection")+
  theme_bw(base_size = 16)+
  geom_line(data = animal_DR_delta_data_lower, 
            aes(x=log10(x_plot), y=prob), color=c("darkgrey"), linetype="dashed")

fig1C <- ggplot()+
  geom_line(data = boot_gam_var1_data,
            aes(x=log10(x_plot), y=prob))+
  geom_ribbon(data = boot_gam_var1_data,
              aes(x=log10(x_plot),ymin=prob_low,ymax=prob_upp), fill = "grey", alpha=0.6)+
  geom_point(data = covid_raw_data.frame,
             aes(x=log10(Dose), y=prob), color="red", size=2)+
  geom_linerange(data = covid_raw_data.frame,
                 aes(x=log10(Dose), ymin=prob_low, ymax=prob_upp, ), color="red")+
  labs(x=expression(paste({Log[10]},"(dose)", sep="")), y = "Probability of infection")+
  theme_bw(base_size = 16)+
  geom_line(data = animal_DR_delta_data, 
            aes(x=log10(x_plot), y=prob), color=c("darkgrey"), linetype="dashed")

fig1D <- ggplot()+
  geom_line(data = boot_DR_delta_data,
            aes(x=log10(x_plot), y=prob))+
  geom_ribbon(data = boot_DR_delta_data,
              aes(x=log10(x_plot),ymin=prob_low,ymax=prob_upp), fill = "grey", alpha=0.6)+
  geom_point(data = covid_raw_data.frame,
             aes(x=log10(Dose), y=prob), color="red", size=2)+
  geom_linerange(data = covid_raw_data.frame,
                 aes(x=log10(Dose), ymin=prob_low, ymax=prob_upp, ), color="red")+
  labs(x=expression(paste({Log[10]},"(dose)", sep="")), y = "Probability of infection")+
  theme_bw(base_size = 16)+
  geom_line(data = animal_DR_delta_data, 
            aes(x=log10(x_plot), y=prob), color=c("darkgrey"), linetype="dashed")

fig2 <- ggplot()+
  geom_line(data = sim_dr_gamma_data %>% dplyr::filter(theta==10^(-6)),
            aes(x=log10(x_plot), y=prob), color=RColorBrewer::brewer.pal(11, "RdBu")[11], size=1)+
  geom_line(data = sim_dr_gamma_data %>% dplyr::filter(theta==1),
            aes(x=log10(x_plot), y=prob), color=RColorBrewer::brewer.pal(11, "RdBu")[10], size=1)+
  #geom_line(data=gamma_plot, aes(x=log10dose, y=pred), color=RColorBrewer::brewer.pal(11, "RdBu")[1], size=1) +
  #geom_ribbon(data=gamma_plot, aes(x=log10dose,ymin=low,ymax=upp), fill = RColorBrewer::brewer.pal(11, "RdBu")[1], alpha=0.1)+
  labs(x=expression(paste({Log[10]},"(dose)", sep="")), y = "Probability of infection")+
  theme_bw(base_size = 20)+
  geom_point(data=data_plot, 
             aes(x=LogDose,y=Prob, size=Total), color=RColorBrewer::brewer.pal(11, "RdBu")[2], alpha=0.6)+
  geom_linerange(data = data_plot,
                 aes(x=LogDose, ymin=Prob_low, ymax=Prob_upp), color=RColorBrewer::brewer.pal(11, "RdBu")[2], alpha=0.6)

((fig1A/fig1C/fig1D) | fig2) + 
  plot_layout(widths = c(1, 2))+
  plot_annotation(tag_levels = "A") #1700x808


###2. Bootstrapped SARS-CoV-2 dose-reseponse curves (Figure S2) -----
###new bootstrapped curves 
traject_boot_gam_var100_data <- data.frame(
  x_plot = rep(x_plot,length(boot_plot_gam_var100[1,])),
  prob = as.numeric(boot_plot_gam_var100),
  loop = rep(1:length(boot_plot_gam_var100[1,]), each=length(boot_plot_gam_var100[,1]))
)
traject_boot_gam_var10_data <- data.frame(
  x_plot = rep(x_plot,length(boot_plot_gam_var10[1,])),
  prob = as.numeric(boot_plot_gam_var10),
  loop = rep(1:length(boot_plot_gam_var10[1,]), each=length(boot_plot_gam_var10[,1]))
)
traject_boot_gam_var1_data <- data.frame(
  x_plot = rep(x_plot,length(boot_plot_gam_var1[1,])),
  prob = as.numeric(boot_plot_gam_var1),
  loop = rep(1:length(boot_plot_gam_var1[1,]), each=length(boot_plot_gam_var1[,1]))
)
traject_boot_gam_var0.1_data <- data.frame(
  x_plot = rep(x_plot,length(boot_plot_gam_var0.1[1,])),
  prob = as.numeric(boot_plot_gam_var0.1),
  loop = rep(1:length(boot_plot_gam_var0.1[1,]), each=length(boot_plot_gam_var0.1[,1]))
)
traject_boot_gam_var0.01_data <- data.frame(
  x_plot = rep(x_plot,length(boot_plot_gam_var0.01[1,])),
  prob = as.numeric(boot_plot_gam_var0.01),
  loop = rep(1:length(boot_plot_gam_var0.01[1,]), each=length(boot_plot_gam_var0.01[,1]))
)
traject_boot_gam_var0.001_data <- data.frame(
  x_plot = rep(x_plot,length(boot_plot_gam_var0.001[1,])),
  prob = as.numeric(boot_plot_gam_var0.001),
  loop = rep(1:length(boot_plot_gam_var0.001[1,]), each=length(boot_plot_gam_var0.001[,1]))
)
##ggplot
fig_traj_var100 <- ggplot()+
  geom_line(data=traject_boot_gam_var100_data, aes(x=log10(x_plot), y=prob, group=loop), alpha=.1, color="blue")+
  geom_point(data = covid_raw_data.frame,
             aes(x=log10(Dose), y=prob), color="red", size=2)+
  geom_linerange(data = covid_raw_data.frame,
                 aes(x=log10(Dose), ymin=prob_low, ymax=prob_upp, ), color="red")+
  labs(x=expression(paste({Log[10]},"(dose)", sep="")), y = "Probability of infection", title = "Var=100")+
  theme_bw(base_size = 20)+
  geom_line(data = animal_DR_delta_data, 
            aes(x=log10(x_plot), y=prob), color=c("darkgrey"), linetype="dashed")
fig_traj_var10 <- ggplot()+
  geom_line(data=traject_boot_gam_var10_data, aes(x=log10(x_plot), y=prob, group=loop), alpha=.1, color="blue")+
  geom_point(data = covid_raw_data.frame,
             aes(x=log10(Dose), y=prob), color="red", size=2)+
  geom_linerange(data = covid_raw_data.frame,
                 aes(x=log10(Dose), ymin=prob_low, ymax=prob_upp, ), color="red")+
  labs(x=expression(paste({Log[10]},"(dose)", sep="")), y = "Probability of infection", title = "Var=10")+
  theme_bw(base_size = 20)+
  geom_line(data = animal_DR_delta_data, 
            aes(x=log10(x_plot), y=prob), color=c("darkgrey"), linetype="dashed")
fig_traj_var1 <- ggplot()+
  geom_line(data=traject_boot_gam_var1_data, aes(x=log10(x_plot), y=prob, group=loop), alpha=.1, color="blue")+
  geom_point(data = covid_raw_data.frame,
             aes(x=log10(Dose), y=prob), color="red", size=2)+
  geom_linerange(data = covid_raw_data.frame,
                 aes(x=log10(Dose), ymin=prob_low, ymax=prob_upp, ), color="red")+
  labs(x=expression(paste({Log[10]},"(dose)", sep="")), y = "Probability of infection", title = "Var=1")+
  theme_bw(base_size = 20)+
  geom_line(data = animal_DR_delta_data, 
            aes(x=log10(x_plot), y=prob), color=c("darkgrey"), linetype="dashed")
fig_traj_var0.1 <- ggplot()+
  geom_line(data=traject_boot_gam_var0.1_data, aes(x=log10(x_plot), y=prob, group=loop), alpha=.1, color="blue")+
  geom_point(data = covid_raw_data.frame,
             aes(x=log10(Dose), y=prob), color="red", size=2)+
  geom_linerange(data = covid_raw_data.frame,
                 aes(x=log10(Dose), ymin=prob_low, ymax=prob_upp, ), color="red")+
  labs(x=expression(paste({Log[10]},"(dose)", sep="")), y = "Probability of infection", title = "Var=0.1")+
  theme_bw(base_size = 20)+
  geom_line(data = animal_DR_delta_data, 
            aes(x=log10(x_plot), y=prob), color=c("darkgrey"), linetype="dashed")
fig_traj_var0.01 <- ggplot()+
  geom_line(data=traject_boot_gam_var0.01_data, aes(x=log10(x_plot), y=prob, group=loop), alpha=.1, color="blue")+
  geom_point(data = covid_raw_data.frame,
             aes(x=log10(Dose), y=prob), color="red", size=2)+
  geom_linerange(data = covid_raw_data.frame,
                 aes(x=log10(Dose), ymin=prob_low, ymax=prob_upp, ), color="red")+
  labs(x=expression(paste({Log[10]},"(dose)", sep="")), y = "Probability of infection", title = "Var=0.01")+
  theme_bw(base_size = 20)+
  geom_line(data = animal_DR_delta_data, 
            aes(x=log10(x_plot), y=prob), color=c("darkgrey"), linetype="dashed")
fig_traj_var0.001 <- ggplot()+
  geom_line(data=traject_boot_gam_var0.001_data, aes(x=log10(x_plot), y=prob, group=loop), alpha=.1, color="blue")+
  geom_point(data = covid_raw_data.frame,
             aes(x=log10(Dose), y=prob), color="red", size=2)+
  geom_linerange(data = covid_raw_data.frame,
                 aes(x=log10(Dose), ymin=prob_low, ymax=prob_upp, ), color="red")+
  labs(x=expression(paste({Log[10]},"(dose)", sep="")), y = "Probability of infection", title = "Var=0.001")+
  theme_bw(base_size = 20)+
  geom_line(data = animal_DR_delta_data, 
            aes(x=log10(x_plot), y=prob), color=c("darkgrey"), linetype="dashed")

((fig_traj_var100|fig_traj_var10|fig_traj_var1)/(fig_traj_var0.1|fig_traj_var0.01|fig_traj_var0.001))+
  plot_annotation(tag_levels = "A") #1700x808, Figure S2


###3. Estimated old coronavirus models (Figure S1) -----
#dose-response curves
delta_plot <- as.data.frame(list(dose=x_plot, log10dose=log10(x_plot), pred= DR_delta(x_plot,est_para_list$delta_nonvac$est), low = delta_nonvac_traject$percentile_bt[1,], upp = delta_nonvac_traject$percentile_bt[3,]))
gamma_plot <- as.data.frame(list(dose=x_plot, log10dose=log10(x_plot), pred= DR_gamma(x_plot,est_para_list$gamma_nonvac$est), low = gamma_nonvac_traject$percentile_bt[1,], upp = gamma_nonvac_traject$percentile_bt[3,]))
tweedie_plot <- as.data.frame(list(dose=x_plot, log10dose=log10(x_plot), pred= DR_tweedie(x_plot,est_para_list$tweedie_nonvac$est), low = tweedie_nonvac_traject$percentile_bt[1,], upp = tweedie_nonvac_traject$percentile_bt[3,]))
twolev_plot <- as.data.frame(list(dose=x_plot, log10dose=log10(x_plot), pred= DR_twolevel(x_plot,est_para_list$twolev_nonvac$est), low = twolev_nonvac_traject$percentile_bt[1,], upp = twolev_nonvac_traject$percentile_bt[3,]))
gamPM_plot <- as.data.frame(list(dose=x_plot, log10dose=log10(x_plot), pred= DR_gammaPointmass(x_plot,est_para_list$gamPM_nonvac$est), low = gamPM_nonvac_traject$percentile_bt[1,], upp = gamPM_nonvac_traject$percentile_bt[3,]))
gamPM0_plot <- as.data.frame(list(dose=x_plot, log10dose=log10(x_plot), pred= DR_gammaPointmass(x_plot,c(0,est_para_list$gamPM0_nonvac$est)), low = gamPM0_nonvac_traject$percentile_bt[1,], upp = gamPM0_nonvac_traject$percentile_bt[3,]))

#ggplot
fig_2Aa <- ggplot() + 
  xlim(0,5) +
  ylim(0,1.1) +
  geom_point(data=data_plot, aes(x=LogDose,y=Prob, size=Total),alpha=0.4) + 
  #scale_color_manual(values = c(RColorBrewer::brewer.pal(11, "RdBu")[11],RColorBrewer::brewer.pal(11, "RdBu")[1])) +
  geom_line(data=Nonvac_plot, aes(x=log10dose, y=pred), color=RColorBrewer::brewer.pal(11, "RdBu")[11], size=1) +
  geom_ribbon(data=Nonvac_plot, aes(x=log10dose,ymin=low,ymax=upp), fill = RColorBrewer::brewer.pal(11, "RdBu")[11], alpha=0.1) +
  theme_bw(base_size = 20) +
  labs(x=expression(paste({Log[10]},"(dose)", sep="")), y = "Probability of infection") +
  #annotate("text", x = c(Inf, Inf,Inf), y=c(1.1, 0.85,0.6),hjust=1,label = c("Unvaccinated","Leaky", "All-or-nothing"), color=c(RColorBrewer::brewer.pal(11, "RdBu")[11],RColorBrewer::brewer.pal(11, "RdBu")[3], RColorBrewer::brewer.pal(11, "RdBu")[4]), size=c(6,6,6))+
  scale_y_continuous(breaks=seq(0,1,0.25))

fig_delta <- ggplot() + 
  xlim(0,5) +
  ylim(0,1.1) +
  geom_point(data=data_plot, aes(x=LogDose,y=Prob, size=Total),alpha=0.4) + 
  geom_line(data=delta_plot, aes(x=log10dose, y=pred), color=RColorBrewer::brewer.pal(11, "RdBu")[11], size=1) +
  geom_ribbon(data=delta_plot, aes(x=log10dose,ymin=low,ymax=upp), fill = RColorBrewer::brewer.pal(11, "RdBu")[11], alpha=0.1) +
  theme_bw(base_size = 20) +
  labs(x=expression(paste({Log[10]},"(dose)", sep="")), y = "Probability of infection")

fig_gamma <- ggplot() + 
  xlim(0,5) +
  ylim(0,1.1) +
  geom_point(data=data_plot, aes(x=LogDose,y=Prob, size=Total),alpha=0.4) + 
  geom_line(data=gamma_plot, aes(x=log10dose, y=pred), color=RColorBrewer::brewer.pal(11, "RdBu")[11], size=1) +
  geom_ribbon(data=gamma_plot, aes(x=log10dose,ymin=low,ymax=upp), fill = RColorBrewer::brewer.pal(11, "RdBu")[11], alpha=0.1) +
  theme_bw(base_size = 20) +
  labs(x=expression(paste({Log[10]},"(dose)", sep="")), y = "Probability of infection")

fig_tweedie <- ggplot() + 
  xlim(0,5) +
  ylim(0,1.1) +
  geom_point(data=data_plot, aes(x=LogDose,y=Prob, size=Total),alpha=0.4) + 
  geom_line(data=tweedie_plot, aes(x=log10dose, y=pred), color=RColorBrewer::brewer.pal(11, "RdBu")[11], size=1) +
  geom_ribbon(data=tweedie_plot, aes(x=log10dose,ymin=low,ymax=upp), fill = RColorBrewer::brewer.pal(11, "RdBu")[11], alpha=0.1) +
  theme_bw(base_size = 20) +
  labs(x=expression(paste({Log[10]},"(dose)", sep="")), y = "Probability of infection")

fig_twolev <- ggplot() + 
  xlim(0,5) +
  ylim(0,1.1) +
  geom_point(data=data_plot, aes(x=LogDose,y=Prob, size=Total),alpha=0.4) + 
  geom_line(data=twolev_plot, aes(x=log10dose, y=pred), color=RColorBrewer::brewer.pal(11, "RdBu")[11], size=1) +
  geom_ribbon(data=twolev_plot, aes(x=log10dose,ymin=low,ymax=upp), fill = RColorBrewer::brewer.pal(11, "RdBu")[11], alpha=0.1) +
  theme_bw(base_size = 20) +
  labs(x=expression(paste({Log[10]},"(dose)", sep="")), y = "Probability of infection")

fig_gamPM <- ggplot() + 
  xlim(0,5) +
  ylim(0,1.1) +
  geom_point(data=data_plot, aes(x=LogDose,y=Prob, size=Total),alpha=0.4) + 
  geom_line(data=gamPM_plot, aes(x=log10dose, y=pred), color=RColorBrewer::brewer.pal(11, "RdBu")[11], size=1) +
  geom_ribbon(data=gamPM_plot, aes(x=log10dose,ymin=low,ymax=upp), fill = RColorBrewer::brewer.pal(11, "RdBu")[11], alpha=0.1) +
  theme_bw(base_size = 20) +
  labs(x=expression(paste({Log[10]},"(dose)", sep="")), y = "Probability of infection")

fig_gamPM0 <- ggplot() + 
  xlim(0,5) +
  ylim(0,1.1) +
  geom_point(data=data_plot, aes(x=LogDose,y=Prob, size=Total),alpha=0.4) + 
  geom_line(data=gamPM0_plot, aes(x=log10dose, y=pred), color=RColorBrewer::brewer.pal(11, "RdBu")[11], size=1) +
  geom_ribbon(data=gamPM0_plot, aes(x=log10dose,ymin=low,ymax=upp), fill = RColorBrewer::brewer.pal(11, "RdBu")[11], alpha=0.1) +
  theme_bw(base_size = 20) +
  labs(x=expression(paste({Log[10]},"(dose)", sep="")), y = "Probability of infection")

((fig_delta|fig_gamma)/(fig_twolev|fig_gamPM))+
  plot_annotation(tag_levels = "A") #1700x808, Figure S1
