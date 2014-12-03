har_fit <- function(times, cens, initial_params, weights_formula) {
  if(missing(weights_formula)) {missing_formula = TRUE} else {missing_formula = FALSE}
  
  if(missing_formula) {
    W = as.matrix(rep(1, length(times)))
  } else {
    W = model.matrix(weights_formula)
  }
  
  d_t = outer(times,times,function(x,y) {abs(x-y)})  #Time Differences
  delta = min(d_t[d_t >0])*10^{-4}
  
  har_lambda <- function(rho,R,d) {
    if(length(d) == 1) {
      return( digamma(sum(R) + sum(d)+rho) - digamma(rho+sum(R)))
    } else {
      R_new = c(R,d[1])
      d_new = d[2:length(d)]
      return( har_lambda(rho, R, d_new) - har_lambda(rho, R_new, d_new) )
    }
  }
  
  total_risk <- function(beta) {
    risk_t <- function(t) {
      weights = exp(W%*%beta)
      return(sum(weights[(times > t & cens == 0) | (times >= t & cens == 1)]))
    }
    return(risk_t)
  }
  
  har_cum_haz <- function(rho,beta) {
    k = length(unique(times[cens == 0]))
    unique_times = unique(times)
    unique_times = unique_times[order(unique_times)]
    diff_times = unique_times - c(0,unique_times[1:(length(unique_times)-1)])
    
    risk_times = c(0,unique_times[1:(length(unique_times)-1)])
    tot_risk = unlist(lapply(risk_times+delta,total_risk(beta)))
    return(sum(diff_times*(digamma(tot_risk+rho) - digamma(rho))))
  }
  
  har_log_lik <- function(params) {
    rho = params[1]
    nu = params[2]
    if(missing_formula) {beta = 0} else{beta = params[3:length(params)]}
    
    k = length(unique(times[cens == 0]))
    weights = exp(W%*%beta)
    
    cts_comp = -nu*har_cum_haz(rho,beta)+k*log(nu)
    atm_comp = 0
    
    for (t in unique(times[cens == 0])) {
      R = weights[  ((times > t) & (cens == 0)) | ((times >= t) & (cens == 1)) ]
      d = weights[  ((times == t) & (cens == 0))  ]	
      atm_comp = atm_comp+log(har_lambda(rho,R,d))
    }
    
    return(-atm_comp - cts_comp)	
  }
  
  if(missing_formula & missing(initial_params)) {
    initial_params = c(1,1)
  }
  
  if(!missing_formula & missing(initial_params)) {
    initial_params = c(1,1,rep(0,ncol(W)))
  }
  
  op_har_loglik = optim(initial_params, har_log_lik,hessian = TRUE, lower = c(0.0001,0.0001))
  
  
  har_par = op_har_loglik$par
  har_var = diag(solve(op_har_loglik$hessian))
  
  ### Cond Survival Functions ###
  har_rho = har_par[1];har_nu = har_par[2]; 
  
  if(missing_formula) {har_beta = 0} else {har_beta = har_par[3:length(har_par)]}
  
  har_cond_surv <- function(weight) {
    cond_surv_t <- function(t) {
      unique_times = unique(times)
      unique_times = unique_times[order(unique_times)]
      weights = exp(W%*%har_beta)
      
      risk_times = unique(c(0,unique_times))
      diff_times = c(risk_times[2:length(risk_times)] - risk_times[1:(length(risk_times)-1)],0)
      
      deaths = rep(0,0)
      for (i in 1:(length(unique_times))) {
        deaths[i] = min(length(which((times == unique_times[i]) & (cens == 0))),1)
      }
      
      R_risk = as.numeric(lapply(risk_times+delta,total_risk(har_beta)))
      L_risk = as.numeric(lapply(risk_times,total_risk(har_beta)))
      
      integrand = har_nu*(digamma(R_risk[which(unique_times <= t)]+weight+har_rho) - digamma(R_risk[which(unique_times <= t)]+har_rho))
      cts_comp = sum(diff_times[which(unique_times <= t)]*integrand)
      
      if(t >= max(unique_times)) {
        cts_comp = cts_comp + (t - max(unique_times))*har_nu*(digamma(weight+har_rho) - digamma(har_rho))			
      } else if (t < min(unique_times)){
        cts_comp = cts_comp + t*har_nu*(digamma(weight+total_risk(har_beta)(0)+har_rho) - digamma(total_risk(har_beta)(0)+har_rho))  		      
      } else {
        cts_comp = cts_comp + (t - max(unique_times[which(unique_times <= t)]))*har_nu*(digamma(weight+R_risk[max(which(unique_times <= t))]+har_rho) - digamma(R_risk[max(which(unique_times <= t))]+har_rho))			
      }	
      
      log_atm_comp = 0
      for(i in unique_times[(0 < unique_times) & (unique_times<= t) & (deaths == 1) ]) {
        death_weight = weights[(times == i) & (cens ==0)]
        log_atm_comp = log_atm_comp + log(har_lambda(har_rho,L_risk[i == risk_times]+weight,death_weight)) - log(har_lambda(har_rho,L_risk[i == risk_times],death_weight))
      }
      
      return(exp(-cts_comp + log_atm_comp))
    }
    return(cond_surv_t)
}
    
  return(list('par' = har_par , 
              'std_err' = sqrt(har_var), 
              'log_lik' = -op_har_loglik$value,
              'conv' = op_har_loglik$convergence,
              'cond_dist' = har_cond_surv
        ))
}

gamma_fit <- function(times, cens, weights_formula, 
                      initial_params, parameterization = "standard") {
  # Make sure there is no constant in the weights formula
  if(missing(weights_formula)) {missing_formula = TRUE} else {missing_formula = FALSE}
  
  if(missing_formula) {
    W = as.matrix(rep(1, length(times)))
  } else {
    W = model.matrix(weights_formula)
  }
  
  d_t = outer(times,times,function(x,y) {abs(x-y)})  #Time Differences
  delta = min(d_t[d_t >0])*10^{-2}
  
  lambda_calc <- function(rho, R, d) {
    if(length(d) == 1) {
      return(log(1+sum(d)/(sum(R)+rho)))
    } else {
      R_new = c(R,d[1])
      d_new = d[2:length(d)]
      return( lambda_calc(rho, R, d_new) - lambda_calc(rho, R_new, d_new) )
    }
  }
  
  total_risk <- function(beta) {
    risk_t <- function(t) {
      weights = exp(W%*%beta)
      return(sum(weights[(times > t & cens == 0) | (times >= t & cens == 1)]))
    }
    return(risk_t)
  }
  
  cum_haz <- function(rho,beta) {
    k = length(unique(times[cens == 0]))
    unique_times = unique(times)
    unique_times = unique_times[order(unique_times)]
    diff_times = unique_times - c(0,unique_times[1:(length(unique_times)-1)])
    
    risk_times = c(0,unique_times[1:(length(unique_times)-1)])
    tot_risk = unlist(lapply(risk_times+delta,total_risk(beta)))
    
    return(sum(diff_times*log(1+tot_risk/rho)))
  }
  
  log_lik <- function(params) {
    if(parameterization == "standard") {
      rho = params[1]
      nu = params[2]
      if(missing_formula) { beta = 0 } else {beta = params[3:length(params)]}
    }
    
    if(parameterization == "ratio") {
      rho = 1/params[1]
      kappa = params[2]
      if(missing_formula) { beta = 0 } else {beta = params[3:length(params)]}
    }
    
    
    k = length(unique(times[cens == 0]))
    weights = exp(W%*%beta)
    
    if (parameterization == "ratio") {
      cts_comp = -(kappa*rho)*cum_haz(rho,beta)+k*log(kappa*rho)
    } 
    if (parameterization == "standard") {
      cts_comp = -(nu)*cum_haz(rho,beta)+k*log(nu)
    }
    atm_comp = 0
    
    for (t in unique(times[cens == 0])) {
      R = weights[  ((times > t) & (cens == 0)) | ((times >= t) & (cens == 1)) ]
      d = weights[  ((times == t) & (cens == 0))  ]	
      atm_comp = atm_comp+log(lambda_calc(rho,R,d))
    }
    
    return(-atm_comp - cts_comp)	
  }
  
  if(missing(initial_params)) {
    if(missing_formula) { 
      initial_params = c(1,1)  
    } else {
      initial_params = c(1,1,rep(0,ncol(W)))  
    }      
  }
  op_loglik = optim(initial_params, log_lik, hessian = TRUE, lower = c(0.0001,0.0001,-Inf))
  
  par = op_loglik$par
  var = diag(solve(op_loglik$hessian))
    
  if(parameterization == "ratio") {
    gamma_rho = 1/par[1]
    gamma_nu = par[2]/par[1] 
  } 
  if (parameterization == "standard") {
    gamma_rho = par[1]
    gamma_nu = par[2]
  }
  if(missing_formula) {gamma_beta = 0} else {gamma_beta = par[3:length(par)]}
  
  gamma_cond_surv <- function(weight) {
    gamma_cond_surv_t <- function(t) {
      unique_times = unique(times)
      unique_times = unique_times[order(unique_times)]
      weights = exp(W%*%gamma_beta)
      
      risk_times = unique(c(0,unique_times))
      diff_times = c(risk_times[2:length(risk_times)] - risk_times[1:(length(risk_times)-1)],0)
      
      deaths = death_weight = rep(0,0)
      for (i in 1:(length(unique_times))) {
        death_weight[i] = sum(weights[which(times == unique_times[i])])
        deaths[i] = min(length(which((times == unique_times[i]) & (cens == 0))),1)
      }
      
      R_risk = as.numeric(lapply(risk_times+0.001,total_risk(gamma_beta)))	
      L_risk = as.numeric(lapply(risk_times,total_risk(gamma_beta)))  
      
      integrand = gamma_nu*log( 1 + weight/(R_risk[which(unique_times <= t)]+gamma_rho))
      cts_comp = sum(diff_times[which(unique_times <= t)]*integrand)
      
      if(t >= max(unique_times)) {
        cts_comp = cts_comp + (t - max(unique_times))*gamma_nu*(log(1+weight/(gamma_rho)))		
      } else if(t < min(unique_times)) {
        cts_comp = cts_comp + t*gamma_nu*(log(1+weight/(total_risk(gamma_beta)(0) + gamma_rho)))		
      } else {
        cts_comp = cts_comp + (t - max(unique_times[which(unique_times <= t)]))*gamma_nu*(log(1+weight/(R_risk[max(which(unique_times <= t))]+gamma_rho)))			
      }	
      
      log_atm_comp = 0
      for(i in unique_times[(0 < unique_times) & (unique_times <= t) & (deaths == 1) ]) {
        death_weight = weights[(times == i) & (cens ==0)]
        log_atm_comp = log_atm_comp + log(lambda_calc(gamma_rho,L_risk[i == risk_times]+weight,death_weight)) - log(lambda_calc(gamma_rho,L_risk[i == risk_times],death_weight))
      }
      
      return(exp(-cts_comp + log_atm_comp))
    }
    return(gamma_cond_surv_t)
  }
  
  return(list('par' = par,
              'std_err' = sqrt(var),
              'log_lik' = -op_loglik$value,
              'conv' = op_loglik$convergence,
              'cond_dist' = gamma_cond_surv, 
              'weights' = exp(W%*%gamma_beta)
        ))
}

invlin_fit <- function(times, cens, initial_params, weights_formula) {
  # Make sure there is no constant in the weights formula
  
  if(missing(weights_formula)) {
    W = rep(1, length(times))
  } else {
    W = model.matrix(weights_formula)
  }
  
  d_t = outer(times,times,function(x,y) {abs(x-y)})  #Time Differences
  delta = min(d_t[d_t >0])*10^{-2}
  
  limit_har_lambda <- function(ratio,R,d) {
    if(length(d) == 1) {
      numer = ratio*(sum(R) + sum(d))
      return( numer/(1+ numer) - ratio*sum(R)/(1+ratio*sum(R)) )
    } else {
      R_new = c(R,d[1])
      d_new = d[2:length(d)]
      return( limit_har_lambda(ratio, R, d_new) - limit_har_lambda(ratio, R_new, d_new) )
    }
  }
  
  limit_total_risk <- function(beta) {
    risk_t <- function(t) {
      weights = exp(W%*%beta)
      return(sum(weights[(times > t & cens == 0) | (times >= t & cens == 1)]))
    }
    return(risk_t)
  }
  
  limit_har_cum_haz <- function(gamma,beta) {
    k = length(unique(times[cens == 0]))
    unique_times = unique(times)
    unique_times = unique_times[order(unique_times)]
    diff_times = unique_times - c(0,unique_times[1:(length(unique_times)-1)])
    
    risk_times = c(0,unique_times[1:(length(unique_times)-1)])
    tot_risk = unlist(lapply(risk_times+delta,limit_total_risk(beta)))
    
    return(sum(diff_times*(gamma*tot_risk/(1+gamma*tot_risk))))
  }
  
  limit_har_loglik <- function(params) {
    theta = params[1]
    gamma = params[2]
    beta = params[3:length(params)]
    
    k = length(unique(times[cens == 0]))
    weights = exp(W%*%beta)
    
    cts_comp = -theta*limit_har_cum_haz(gamma,beta)+k*log(theta)
    atm_comp = 0
    
    for (t in unique(times[cens == 0])) {
      R = weights[  ((times > t) & (cens == 0)) | ((times >= t) & (cens == 1)) ]
      d = weights[  ((times == t) & (cens == 0))  ]  
      atm_comp = atm_comp+log(limit_har_lambda(gamma,R,d))
    }
    
    return(-atm_comp - cts_comp)	
    
  }
  
  if(missing(initial_params)) {
    initial_params = c(1,1,rep(0,ncol(W)))
  }
  
  limit_har_oplik = optim(initial_params, limit_har_loglik, lower = c(0.0001,0.0001,-Inf), hessian = TRUE)
  
  limit_har_par = limit_har_oplik$par
  limit_har_stderr = sqrt(diag(solve(limit_har_oplik$hessian)))
  
  
  ###  Conditional Survival
  theta = limit_har_par[1];gamma = limit_har_par[2]; beta = limit_har_par[3:length(limit_har_par)]
  
  limit_cond_surv <- function(weight) {
    cond_surv_t <- function(t) {
      unique_times = unique(times)
      unique_times = unique_times[order(unique_times)]
      weights = exp(W%*%beta)
      
      risk_times = unique(c(0,unique_times))
      diff_times = c(risk_times[2:length(risk_times)] - risk_times[1:(length(risk_times)-1)],0)
      
      deaths = rep(0,0)
      for (i in 1:(length(unique_times))) {
        deaths[i] = min(length(which((times == unique_times[i]) & (cens == 0))),1)
      }
      
      R_risk = as.numeric(lapply(risk_times+delta,limit_total_risk(beta)))
      L_risk = as.numeric(lapply(risk_times,limit_total_risk(beta)))
      
      integrand = theta*gamma*weight/((1+gamma*(R_risk[which(unique_times <= t)]+weight))*(1+gamma*R_risk[which(unique_times <= t)]))
      cts_comp = sum(diff_times[which(unique_times <= t)]*integrand)      
      
      if(t >= max(unique_times)) {
        cts_comp = cts_comp + (t - max(unique_times))*theta*(gamma*weight/(1+gamma*weight))			
      } else if (t < min(unique_times)){
        cts_comp = cts_comp + t*theta*gamma*weight/((1+(weight+limit_total_risk(beta)(0))*gamma)*(1+limit_total_risk(beta)(0)*gamma))
      } else {
        cts_comp = cts_comp + (t - max(unique_times[which(unique_times <= t)]))*theta*weight*gamma/ ((1+(weight+R_risk[max(which(unique_times <= t))])*gamma)*(1+(R_risk[max(which(unique_times <= t))])*gamma))
      }	
      
      log_atm_comp = 0
      for(i in unique_times[(0 < unique_times) & (unique_times<= t) & (deaths == 1) ]) {
        death_weight = weights[(times == i) & (cens ==0)]
        log_atm_comp = log_atm_comp + log(limit_har_lambda(gamma,L_risk[i == risk_times]+weight,death_weight)) - log(limit_har_lambda(gamma,L_risk[i == risk_times],death_weight))
      }
      
      return(exp(-cts_comp + log_atm_comp))
    }
    return(cond_surv_t)
  }
  
  return(list('par' = limit_har_par, 
              'std_err' = limit_har_stderr,
              'log_lik' = -limit_har_oplik$value,
              'conv' = limit_har_oplik$convergence,
              'cond_dist' = limit_cond_surv,
              'weights' = exp(W%*%beta)
        ))
}