Markov Survival Processes :  Fitting Software
===============

Two files, one which contains functions for estimation of parameters under three different markov survival processes (harmonic, gamma, and inverse linear) and another which applies these functions to the a particular dataset (Gehan, (1965)).

1.  markov_process_fit.R :  Contains three functions (har_fit, gamma_fit, and invlin_fit)
  * Inputs : 
    * times : The observed survival or censoring times for each patiet
    * cens : 
    * initial_params : 
    * weights_formula : 
2.  gehan_example.R : 
