// ---------- Data block ----------
// All inputs required for the likelihood and covariates
data {
  int<lower=1> N;                                                               // number of fish
  int<lower=1> T;                                                               // number of timesteps
  int<lower=1> L;                                                               // number of locations
  vector[N] frac;                                                               // fraction of first interval experienced
  int<lower=1> t0[N];                                                           // release time
  int<lower=0> R[N];                                                            // recapture time (0 if censored)
  int<lower=1> loc[N];                                                          // location
  vector[N] Depth;                                                              // depth 
  vector[N] Temp;                                                               // temperature 
  vector[N] Selectivity;                                                        // size-selectivity
  matrix[T,L] Rec;                                                              // matrix of recreational effort at each time t and location l
  matrix[T,L] FHO;                                                              // matrix of for-hire fishing pressure at each time t and location l using observer data
}

// ---------- Parameters ----------
parameters {
  real gamma;                                                                   // intercept natural loss 
  vector[3] alpha;                                                              // intercept and slopes for encounter vs. Rec and FHO
  real<lower=0> Depth_thresh;                                                   // depth threshold determining asymptotic behavior
  real log_beta_D;                                                              // log depth effect
  real log_beta_0;                                                              // log baseline survival/reporting
}

// ---------- Transformed parameters ----------
transformed parameters {
  matrix[T,L] phi;                                                              // survival probability per interval (per t,l)
  matrix[T,L] h_recap;                                                          // encounter hazard per interval (per t,l)
  real beta_0 = -exp(log_beta_0);                                               // force beta < 0 (monotone decrease)
  real beta_D = -exp(log_beta_D);                                               // force beta < 0 (monotone decrease)
  vector[N] f_D = (log1p_exp((Depth- Depth_thresh))-log1p_exp(-Depth_thresh));  // depth effect: scaled saturating function with max = 0
  vector[N] psi = exp(beta_0 + beta_D*f_D);
  // Fill time × location arrays for survival and encounter
  for(l in 1:L){
    for(t in 1:T){
      // conditional survival between t and t+1
      phi[t,l] = exp(-exp(gamma)); 
      // probability of being re-encountered in step t
      h_recap[t,l] = exp(alpha[1] + alpha[2]*Rec[t,l] + alpha[3]*FHO[t,l]);
    }
  }
}

// ---------- Model block ----------
model {
  // Priors
  alpha[1] ~ normal(-2,1);                                                      // intercept for encounter
  alpha[2] ~ normal(0,2);                                                       // effect of Rec effort
  alpha[3] ~ normal(0,2);                                                       // effect of FHO effort
  gamma ~ normal(-2.67,1);                                                      // baseline natural loss (from stock assessment)
  log_beta_D ~ normal(0,1);                                                     // prior of log-scale
  log_beta_0 ~ normal(0,1);                                                     // prior of log-scale
  Depth_thresh ~ normal(20, 5);                                                 // depth threshold
  // ---- Likelihood ----
  for (i in 1:N) {
    int t_0i = t0[i];       // release step for fish i
    int end  = R[i];        // recapture time (0 if censored)
    int l    = loc[i];      // location
    vector[T] eta = 1-exp(-h_recap[,l]*Selectivity[i]);
    // survival/encounter for partial first interval
    real phi_first = pow(phi[t_0i,l], frac[i]);
    real eta_first = 1 - pow(1 - eta[t_0i], frac[i]);

    if (R[i] == 0) {
      // -------- Fish i never recaptured (right-censored) --------
      real log_p_recapture_sum = negative_infinity(); // start log-sum-exp at log(0)
      for (t in t_0i:T) {
        // log probability of first recapture occurring at time t
        real log_term = log(psi[i]) + log(phi_first) + log1m(eta_first);
        for (s in (t_0i + 1):(t - 1)){log_term += log(phi[s,l]) + log1m(eta[s]);}
        log_term += log(phi[t,l]) + log(eta[t]);
        log_p_recapture_sum = log_sum_exp(log_p_recapture_sum, log_term);
      }
      // add log probability of no recapture = log(1 - p(any recapture))
      target += log1m_exp(log_p_recapture_sum);

    } else {
      // -------- Fish i recaptured at 'end' --------
      real log_Recap_prob = log(psi[i]) + log(phi_first) + log1m(eta_first);
      for (t in (t_0i + 1):(end - 1)){log_Recap_prob += log(phi[t,l]) + log1m(eta[t]);}
      log_Recap_prob += log(phi[end,l]) + log(eta[end]);
      target += log_Recap_prob;
    }
  }
}

// ---------- Generated quantities block ----------
// Pointwise log-likelihood for model comparison (loo, waic)
generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    int t_0i = t0[i];
    int end  = R[i];
    int l    = loc[i];
    vector[T] eta = 1-exp(-h_recap[,l]*Selectivity[i]);
    real phi_first = pow(phi[t_0i,l], frac[i]);
    real eta_first = 1 - pow(1 - eta[t_0i], frac[i]);

    if (R[i] == 0) {
      real log_p_recapture_sum = negative_infinity();
      for (t in t_0i:T) {
        real log_term = log(psi[i]) + log(phi_first) + log1m(eta_first);
        for (s in (t_0i + 1):(t - 1)){log_term += log(phi[s,l]) + log1m(eta[s]);}
        log_term += log(phi[t,l]) + log(eta[t]);
        log_p_recapture_sum = log_sum_exp(log_p_recapture_sum, log_term);
      }
      log_lik[i] = log1m_exp(log_p_recapture_sum);
    } else {
      real log_Recap_prob = log(psi[i]) + log(phi_first) + log1m(eta_first);
      for (t in (t_0i + 1):(end - 1)){log_Recap_prob += log(phi[t,l]) + log1m(eta[t]);}
      log_Recap_prob += log(phi[end,l]) + log(eta[end]);
      log_lik[i] = log_Recap_prob;
    }
  }
}
