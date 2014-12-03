Markov Survival Processes :  Fitting Software
===============

Two files, one which contains functions for estimation of parameters under three different markov survival processes (harmonic, gamma, and inverse linear) and another which applies these functions to a particular dataset (Gehan, (1965)).

1.  markov_process_fit.R :  Contains three functions (har_fit, gamma_fit, and invlin_fit)
  * Inputs : 
    * times : The observed survival or censoring times for each patient
    * cens : A vector indicating whether the observed time is a censoring or survival time for each patient
    * initial_params (optional) : Initialization of parameters.
    * weights_formula (optional) : The formula for the covariate matrix, W, associated with the weights = exp(W * beta).  If missing, weights are set to be identically one for each patient.
    * parameterization (optional) : For the gamma process, choose between either the 'ratio' or 'standard' parameterization.
  * Output : 
    * par : Maximum likelihood parameter estimates
    * std_err : Standard errors for parameter estimates
    * log_lik : Log-likelihood at the mle
    * conv : Convergence of the parameter estimates (useful when the mle is at a boundary case)
    * cond_dist : Function for the conditional survival distribution at t for a weight given observed times

2.  gehan_example.R : 
 * 6-MP subset of leukemia patients (Gehan, (1965))
   * Fit models under both harmonic and gamma process
    * Produce conditional survival distribution plots (including Kaplan-Meier product limit estimate and the exponential distribution estimate)
 * Complete leukemia dataset (Gehan, (1965)) : Treatment and Control Groups
   * Fit models under both inverese linear (limiting harmonic) and gamma process
    * Produce conditional survival distribution plot.
