source('markov_process_fit.R')

### Just the Treatment Group ###

gehan_times = c(6,6,6,6,7,9,10,10,11,13,16,17,19,20,22,23,25,32,32,34,35) 

gehan_cens = c(0,0,0,1,0,1,0,1,1,0,0,1,1,1,0,0,1,1,1,1,1)

gehan_treatment = as.matrix(rep(1,length(gehan_times)))

times = gehan_times; cens = gehan_cens; treatment = gehan_treatment

# Fit Gamma and Harmonic Models
gfit = gamma_fit(times,cens)

cbind(gfit$par,gfit$std_err) # Gamma process parameter estimates with standard errors

hfit = har_fit(times,cens)  # Not right yet

cbind(hfit$par,hfit$std_err) # Harmonic process parameter estimates with standard errors

# Conditional Plots
k_m <- function(t) {
  unique_times = unique(times)
  unique_times = c(0,unique_times[order(unique_times)])
  
  risk = unlist(lapply(unique_times,total_risk(0)))  
  
  deaths = death_unique = rep(0,0)
  for (i in 1:length(unique_times)) {
    t_i = unique_times[i]
    deaths[i] = length(which((times == t_i) & (cens == 0)))
    death_unique[i] = min(deaths[i],1)
  }
  
  atm_comp = prod(risk[(unique_times <= t) & (death_unique == TRUE)]/(risk[(unique_times <= t) & (death_unique == TRUE)]+deaths[(unique_times <= t) & (death_unique == TRUE)]))
  
  return(atm_comp)
}

rate = sum(cens == 0)/sum(times)
death_times = unique(times[cens==0])

obs_times = c(seq(0,40,0.1),death_times+0.001, death_times-0.001)
obs_times = obs_times[order(obs_times)]

gamma_surv = unlist(lapply(obs_times,gfit$cond_dist(1)))
har_surv = unlist(lapply(obs_times,hfit$cond_dist(1)))
Km_surv = unlist(lapply(obs_times,k_m))

exp_surv = 1-pexp(obs_times, rate = rate)

#png("cond_surv_plot.png", width = 6,height = 4, units = "in", res = 300)

par(mar= c(5,4,1,2)+0.1)

plot(obs_times,har_surv, type = "l", xlab = "Time (in Weeks)", ylab = "Conditional Survival Function", axes = FALSE, ylim = c(0,1))
axis(side = 1); axis(side = 2)

lines(obs_times, gamma_surv, lty = 2)

lines(obs_times, Km_surv, lty = 1, col = "red")

lines(obs_times, exp_surv, lty = 2, col = "red")

legend(2,0.4, c("Harmonic", "Gamma", "Kaplan-Meier", "Exponential"), lty = c(1,2,1,2), col = c("black","black","red","red"), cex = 0.75)
#dev.off()


### All Gehan Data ### 

treatment = as.matrix(c(rep(0, 21),rep(1,21)))

group_1 = c(1,1,2,2,3,4,4,5,5,8,8,8,8,11,11,12,12,15,17,22,23)
group_0 = c(6,6,6,6,7,9,10,10,11,13,16,17,19,20,22,23,25,32,32,34,35)

cens_1 = rep(0,21)
cens_0 = c(0,0,0,1,0,1,0,1,1,0,0,1,1,1,0,0,1,1,1,1,1)

times = c(group_1,group_0)
cens = c(cens_1,cens_0)

# Estimation
gfit_covariates <- gamma_fit(times,cens,weights_formula = ~treatment-1, parameterization = "ratio")
                      
cbind(gfit_covariates$par,gfit_covariates$std_err)

ilfit_covariates <- invlin_fit(times,cens,weights_formula = ~treatment-1)

cbind(ilfit_covariates$par,ilfit_covariates$std_err)

# Conditional Distribution 
treat_weight = min(gfit$weights); control_weight = max(gfit$weights)

death_times = unique(times[cens==0])

obs_times = c(seq(0,35,0.1),death_times+0.001, death_times-0.001)
obs_times = obs_times[order(obs_times)]

gamma_surv_treat = unlist(lapply(obs_times,gfit_covariates$cond_dist(treat_weight)))
gamma_surv_control = unlist(lapply(obs_times,gfit_covariates$cond_dist(control_weight)))

invlin_surv_treat = unlist(lapply(obs_times,ilfit_covariates$cond_dist(treat_weight)))
invlin_surv_control = unlist(lapply(obs_times,ilfit_covariates$cond_dist(control_weight)))


# png("cond_surv_plot.png", width = 6,height = 4, units = "in", res = 300)

par(mar= c(5,4,1,2)+0.1)

plot(obs_times,gamma_surv_treat, type = "l", xlab = "Time (in Weeks)", ylab = "Conditional Survival Function", axes = FALSE, ylim = c(0,1), lty = 2, col = "red")
axis(side = 1); axis(side = 2)

lines(obs_times, gamma_surv_control, lty = 1, col = "red")

lines(obs_times, invlin_surv_control, lty = 1)

lines(obs_times, invlin_surv_treat, lty = 2)

legend(25,1.0, c("Treatment", "Control"), lty = c(2,1), cex = 0.6)
legend(25,0.85, c("Harmonic", "Gamma"), lty = c(1,1), col = c("black", "red"), cex = 0.6)

# dev.off()

