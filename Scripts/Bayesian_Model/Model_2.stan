data {
  int<lower=0> n_obs_rec;						# Number of observations in secondary forests
  int<lower=0> n_obs_old;						# Number of observations in old-growth forests
  int<lower=0> n_site;  						# Number of chronosequences
  vector[n_obs_rec] EA_rec;						# Value of ecosystem attributes in secondary forests
  vector[n_obs_old] EA_old;						# Value of ecosystem attributes in old-growth forests
  vector[n_obs_rec] t;							# Age of secondary forests
  int<lower=1, upper=n_site> site_rec[n_obs_rec];			# Identity of the chronosequence of the observation in secondary forests
  int<lower=1, upper=n_site> site_old[n_obs_old];			# Identity of the chronosequence of the observation in old-growth forests
}

parameters {
  real <lower=0> theta_infty;						# Hyperparameter of the asymptotic value of the ecosystem attribute
  real <lower=0> theta_0;						# Hyperparameter of the initial value of the ecosystem attribute
  real <lower=0.001> lambda;						# Hyperparameter of the recovery rate of the ecosystem attribute
  real theta_over;							# Hyperparameter of the value of the over/below-shoot of the short-term process
  real <lower=5> theta_time;						# Hyperparameter of the time of the over/below-shoot of the short-term process
  vector <lower=0> [n_site] theta_infty_s;				# Per chronosequence asymptotic values of the ecosystem attribute
  vector <lower=0> [n_site] theta_0_s;					# Per chronosequence initial values of the ecosystem attribute
  vector <lower=0.001> [n_site] lambda_s;				# Per chronosequence recovery rate values of the ecosystem attribute
  vector [n_site] theta_over_s;						# Per chronosequence value of the over/below-shoot
  vector <lower=5> [n_site] theta_time_s;				# Per chronosequence value of the time of the over/below-shoot
  real<lower=0> sigma_rec;						# Dispersion parameter of the lognormal likelihood for secondary forests
  real<lower=0> sigma_old;						# Dispersion parameter of the lognormal likelihood for old-growth forests
  real<lower=0> sigma_lambda;						# Dispersion parameter of the hyperlaw on recovery rates lambda_s
  real<lower=0> sigma_theta_infty;					# Dispersion parameter of the hyperlaw on asymptotic values theta_infty_s
  real<lower=0> sigma_theta_0;						# Dispersion parameter of the hyperlaw on initial values theta_0_s
  real<lower=0> sigma_theta_over;					# Dispersion parameter of the hyperlaw on the value of the over/below-shoot of the short-term process
  real<lower=0> sigma_theta_time;					# Dispersion parameter of the hyperlaw on the time of the over/below-shoot
}

model {
  real mu_rec[n_obs_rec];						# Predicted values of ecosystem attributes for secondary plots
  real mu_old[n_obs_old];						# Predicted values of ecosystem attributes for old-growth plots

for (j in 1:n_obs_old)							# Loop on observations in old-growth forests to calculate the predicted values
  {
    mu_old[j]= theta_infty_s[site_old[j]];
  }

  for (i in 1:n_obs_rec)						# Loop on observations in secondary forests to calculate the predicted values
  {
    mu_rec[i]= 	theta_0_s[site_rec[i]] +
    		(theta_infty_s[site_rec[i]]-theta_0_s[site_rec[i]])*
    		(1-exp(-lambda_s[site_rec[i]]*t[i])
    		+
    		theta_over_s[site_rec[i]]* exp(-pow(log(t[i]/theta_time_s[site_rec[i]]),2))
    )
    ;
  }
  
  EA_rec ~ lognormal(log(mu_rec), sigma_rec);				# Model Likelihood for secondary forests
  EA_old ~ lognormal(log(mu_old), sigma_old);				# Model Likelihood for old-growth forests
  
  lambda_s~ lognormal(log(lambda), sigma_lambda);			# Hyperlaw on recovery rates
  theta_infty_s~ lognormal(log(theta_infty), sigma_theta_infty);	# Hyperlaw on asymptotic values
  theta_0_s~ lognormal(log(theta_0), sigma_theta_0);			# Hyperlaw on initial values
  theta_over_s~ normal(theta_over, sigma_theta_over); 			# Hyperlaw on the value of the over/below-shoot of the short-term process
  theta_time_s~ lognormal(log(theta_time), sigma_theta_time); 		# Hyperlaw on the time of the over/below-shoot
}