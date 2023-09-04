// peathamstr.stan with the following additions:
// the acrotelm-catotelm-boundary is modeled as a transition time (which is constant over time)
// the litter addition rate is allowed to vary over time
// acrotelm and catotelm have their own decomposition rates which are constant over time

// Hamstr with additional error due to bioturbation
// 11.10.2020 Andrew Dolman
// Bioturbation modelling
// Prior gamma distribution on L
// vector of n_ind  = no. individual particles in a 14C measurement
// latent variable approach - model bt_age,
data{
  // age control points
  int<lower=0> N;
  vector[N] depth;
  vector[N] obs_age;
  vector[N] obs_err;
  real min_age;

  // resolution of age-depth model
  int<lower = 1> n_lvls; // number of hierarchical levels, not including overall mean
  int<lower=0> K_fine;  // number of highest resolution sections
  int<lower=0> K_tot;  // total no of gamma parameters

  int parent[K_tot]; // index sections to their parent sections

  // modelled depths
  vector[K_fine] c_depth_bottom;
  vector[K_fine] c_depth_top;
  real<lower = 0> delta_c; // width of each highest resolution section


  // hyperparameters for the gamma innovations

  // prior for the oversall mean accumulation rate
  real<lower = 0> acc_mean_prior;

  // shape of the gamma distributions
  real<lower = 0> acc_shape;

  // scale the shape parameter to control the total variance of the alpha
  // innovations for the number of hierarchical levels
  int<lower=0, upper=1> scale_shape;

  // hyperparameters for prior distribution on memory strength (the AR1 coefficient)
  real<lower = 0> mem_mean;
  real<lower = 0> mem_strength;

  // observation error model parameters

  int<lower=0> nu; // degrees of freedom of t error distribution
  int which_c[N]; // index observations to their fine sections

  int<lower=0, upper=1> scale_R; // scale the AR1 coefficient or not
  int<lower=0, upper=1> inflate_errors; // use error inflation model or not


  // hyperparameters for the error inflation model
  real<lower = 0> infl_shape_shape;
  real<lower = 0> infl_shape_mean;

  real<lower = 0> infl_sigma_sd;


  // Additional data for bioturbation model

  // use bioturbation model or not
  int<lower=0, upper=1> model_bioturbation;

  int<lower=0> I;
  int smooth_i[I, N];

  real<lower = 0> L_prior_mean;
  real<lower = 0> L_prior_shape;

  vector[N*model_bioturbation] n_ind;

  // Additional data for modelling displacement
  int<lower=0, upper=1> model_displacement;
  real<lower = 0> D_prior_scale;

  int<lower=0, upper=1> smooth_s;

  // Model hiatuses

  int<lower=0, upper=1> model_hiatus;
  real H_top;
  real H_bottom;

  // Clymo model
  int<lower=0> N2; // number of sections for which there are measured mass values
  int<lower=0> N2_c; // number of simulated sections for mass data
  //vector[N2] cumulative_mass; // cumulative mass till depths where sample were taken
  real<lower = 0> p1_layer_mass_shape;
  real<lower = 0> p2_layer_mass_shape;
  real<lower = 0> p3_layer_mass_shape;
  real<lower = 0> p1_clymo_par;
  real<lower = 0> p2_clymo_par;
  real<lower = 0> p3_clymo_par;
  real<lower = 0> p1_clymo_alpha_1;
  real<lower = 0> p2_clymo_alpha_1;
  real<lower = 0> p3_clymo_alpha_1;
  real<lower = 0> cumulative_mass0;

  real<lower = 0> p1_clymo_alpha_2;
  real<lower = 0> p2_clymo_alpha_2;
  real<lower = 0> p3_clymo_alpha_2;

  real<lower = 0> p1_ac_age;
  real<lower = 0> p2_ac_age;
  real<lower = 0> p3_ac_age;

  int<lower = 1, upper = N2_c> index_has_mass_measurements[N2];
  int<lower = 1, upper = K_fine> which_c2_all[N2_c + 1]; // index for the sections in which all sample_depth_upper and the maximum sample_depth_lower fall
  vector<lower = 0>[N2_c + 1] depth2_all;
  vector[N2] layer_mass;

  real<lower = 0> clymo_par_memory_p1;
  real<lower = 0> clymo_par_memory_p2;

  int<lower = 1, upper = N2_c> index_clymo_par_constant;

  real<lower = 0> p1_age0;
  real<lower = 0> p2_age0;
  real<lower = 0> p3_age0;

}
transformed data{

  real min_depth = min(c_depth_top);
  real max_depth = max(c_depth_bottom);
  real data_age_range = max(obs_age) - min(obs_age);

  // inverse scale of the prior on L
  real L_rate;
  //real D_rate;

  int<lower = 0, upper = 1> sample_L;

  // transform mean and strength of memory beta distribution to alpha and beta
  // as used to parameterise beta dist in stan function
  real<lower=0> mem_alpha = mem_strength * mem_mean;
  real<lower=0> mem_beta = mem_strength * (1-mem_mean);

  // position of the first highest resolution innovation (alpha)
  int<lower = 1> first_K_fine = K_tot - K_fine+1;

  // scale shape
  real<lower = 1> acc_shape_adj;
  if (scale_shape == 1){
    acc_shape_adj = acc_shape * n_lvls;
  } else{
    acc_shape_adj = acc_shape;
  }

  L_rate = L_prior_shape / L_prior_mean;

  // bioturbation depth can be fixed rather than sampled
  if (L_prior_shape == 0) {
    sample_L = 0;
  } else {
    sample_L = 1;
  }


  //D_rate = 1/D_prior_scale;

}
parameters {
  // AR1 coeffiecient at 1 depth unit
  real<lower = 0, upper = 1> R;

  // the hierarchical gamma innovations in one long vector that will be indexed
  vector<lower = 0>[K_tot] alpha;

  // the age at the first modelled depth
  real<lower = min_age> age0;

  // the measurement error inflation factors
  // these have length 0 if inflate_errors == 0 meaning that the parameters are
  // in scope, so the model runs, but are zero length so nothing is sampled
  real<lower = 0> infl_mean[inflate_errors];
  real<lower = 0> infl_shape_1[inflate_errors];
  vector<lower = 0>[inflate_errors ? N : 0] infl;


  real<lower = 0> L[model_bioturbation * sample_L];

  vector<lower = 0>[model_bioturbation ? N : 0] bt_error;

  real<lower = 0> D[model_displacement];

  real<lower = H_top, upper = H_bottom> H_depth[model_hiatus];
  real<lower = 0, upper = data_age_range> H_length[model_hiatus];

  // Clymo model
  vector<lower = 0>[N2_c] clymo_par_raw; // peat addition rate at catotelm top
  real<lower=0>  clymo_alpha_1_raw; // peat decomposition rate at each modeled depth
  real<lower=0>  clymo_alpha_2_raw;
  real<lower=0>  ac_age_raw;
  real<lower = 0> layer_mass_shape;
  real<lower = 0, upper = 1> clymo_par_memory;

}

transformed parameters{

  // the AR1 coefficient scaled for the thickness of the modelled sediment sections
  real<lower = 0, upper = 1> w;

  // the highest resolution AR1 correlated innovations
  vector[K_fine] x;

  // the modelled ages
  vector[K_fine+1] c_ages;

  // the modelled ages interpolated to the positions of the data
  vector[N] Mod_age;

  // the inflated observation errors
  real<lower = 0> infl_shape[inflate_errors];

  vector[(inflate_errors == 1 || model_displacement == 1) ? N : 0] obs_err_infl;
  //vector[model_displacement ? N : 0] disp_err;

  // latent bioturbation corrected age
  //vector[model_bioturbation ? N : 0] bt_age;
  vector[model_bioturbation ? N : 0] bt_age;
  vector[(model_bioturbation == 1 || model_displacement == 1) ? N : 0] smooth_x;

  // age heterogeneity due to bioturbation at locations of observed ages
  vector[model_bioturbation ? N : 0] age_het;

  vector<lower = 0>[model_displacement ? N : 0] disp_yrs;

  if (scale_R == 1){
    w = R^(delta_c);
  } else {
    w = R;
  }

  // only the "fine" alphas
  // the first innovation
  x[1] = alpha[first_K_fine];

  // the remaining innovations with the AR1 parameter applied
  for(i in 2:K_fine){
    x[i] = w*x[i-1] + (1-w)*alpha[i + first_K_fine -1];
  }


  // the cumulative sum of the highest resolution innovations
  c_ages[1] = age0 * p3_age0;
  c_ages[2:(K_fine+1)] = c_ages[1] + cumulative_sum(x * delta_c);



  // hiatus vector
  if (model_hiatus == 1){
    for (i in 2:(K_fine+1)){
      if (H_depth[1] < c_depth_top[i-1]) c_ages[i] = c_ages[i] + H_length[1];
    }
  }


  // age model interpolated to the positions of the observations
  Mod_age = c_ages[which_c] + x[which_c] .* (depth - c_depth_top[which_c]);


  if (model_bioturbation == 1 || model_displacement == 1){
    if (smooth_s == 1){
      for (n in 1:N){
        smooth_x[n] =  mean(x[smooth_i[,n]]);
      }
    } else {
      for (n in 1:N){
        smooth_x[n] =  x[which_c[n]];
      }
    }
  }


  if (model_bioturbation == 1){

    if (sample_L == 1){
      age_het = L[1] * smooth_x;
      // the modelled (shifted) gamma distributed bioturbation error
      // subtract the age_het to centre the error around the obs_age
      bt_age = obs_age + bt_error - age_het;

    }
    if (sample_L == 0){
      age_het = L_prior_mean * smooth_x;
      // the modelled (shifted) gamma distributed bioturbation error
      // subtract the age_het to centre the error around the obs_age
      bt_age = obs_age + bt_error - age_het;


    }
  }

  if (inflate_errors == 1 && model_displacement == 1){
    disp_yrs = D[1] * smooth_x;
    obs_err_infl = sqrt((obs_err .* obs_err) + (infl .* infl) + (disp_yrs .* disp_yrs));
    infl_shape[1] = infl_shape_1[1] + 1;
  } else if (inflate_errors == 1 && model_displacement == 0){
    for (n in 1:N)
    //obs_err_infl[n] = obs_err[n] + infl_sigma[1] * infl[n];
    obs_err_infl[n] = sqrt((obs_err[n])^2 + (infl[n])^2);
    infl_shape[1] = infl_shape_1[1] + 1;
  } else if (inflate_errors == 0 && model_displacement == 1){
    disp_yrs = D[1] * smooth_x;
    obs_err_infl = sqrt((obs_err .* obs_err) + (disp_yrs .* disp_yrs));
  }

  // Clymo model
  vector<lower = 0>[N2_c] clymo_par;
  clymo_par[N2_c] = clymo_par_raw[N2_c];
  clymo_par[1:(N2_c - 1)] = clymo_par_raw[2:N2_c] * clymo_par_memory + clymo_par_raw[1:(N2_c - 1)] * (1 - clymo_par_memory);
  clymo_par[1:index_clymo_par_constant] = rep_vector(clymo_par[index_clymo_par_constant], index_clymo_par_constant);
  real<lower = 0> clymo_alpha_2 = clymo_alpha_2_raw * p3_clymo_alpha_2; // decomposition rate in catotelm
  real<lower = 0> clymo_alpha_1 = clymo_alpha_2 + clymo_alpha_1_raw * p3_clymo_alpha_1; // decomposition rate in acrotelm
  real<lower = 0> ac_age = ac_age_raw * p3_ac_age; // age at which acrotelm peat becomes catotelm peat
  vector<lower = 0>[N2_c] Mod_layer_mass; // mass at sampling depths
  vector<lower = 0>[N2_c] Mod_layer_mass_cumulative; // mass at sampling depths
  vector<lower = 0>[N2_c + 1] Mod_age2_all = c_ages[which_c2_all] + x[which_c2_all] .* (depth2_all - c_depth_top[which_c2_all]);
  vector<lower = 0>[N2_c] Mod_age2_duration;
  vector[N2_c] Mod_age2_duration_acrotelm_next;
  vector[N2_c] Mod_age2_duration_acrotelm_here;
  vector[N2_c] Mod_age2_duration_catotelm_next;
  vector[N2_c] Mod_age2_duration_catotelm_here;

  {
    vector[N2_c] Mod_age2_upper = Mod_age2_all[1:N2_c];
    vector[N2_c] Mod_age2_lower = Mod_age2_all[2:(N2_c + 1)];
    Mod_age2_duration = Mod_age2_lower - Mod_age2_upper;
    Mod_age2_duration_acrotelm_here = rep_vector(ac_age, N2_c);
    Mod_age2_duration_catotelm_here = Mod_age2_duration - ac_age;
    Mod_age2_duration_acrotelm_next = rep_vector(0.0, N2_c);
    Mod_age2_duration_catotelm_next = Mod_age2_upper;
    for(n in 1:N2_c) {
      if(ac_age > Mod_age2_duration[n]) { // AC boundary is outside the layer
        Mod_age2_duration_acrotelm_here[n] = Mod_age2_duration[n];
        Mod_age2_duration_catotelm_here[n] = 0.0;
        Mod_age2_duration_catotelm_next[n] = Mod_age2_upper[n] - (ac_age - Mod_age2_duration_acrotelm_here[n]);
        if(Mod_age2_duration_catotelm_next[n] < 0.0) {
          Mod_age2_duration_catotelm_next[n] = 0.0;
        }
        Mod_age2_duration_acrotelm_next[n] = Mod_age2_upper[n] - Mod_age2_duration_catotelm_next[n];
      }
    }

    // layer mass, under the assumption that the total mass was added to the top at the beginning of the time interval
    Mod_layer_mass =
      clymo_par * p3_clymo_par .* Mod_age2_duration .*
      exp(- clymo_alpha_1 * (Mod_age2_duration_acrotelm_here + Mod_age2_duration_acrotelm_next)) .*
      exp(- clymo_alpha_2 * (Mod_age2_duration_catotelm_here + Mod_age2_duration_catotelm_next));

    Mod_layer_mass_cumulative = cumulative_mass0 + cumulative_sum(Mod_layer_mass);

    //Mod_layer_mass =
    //  (clymo_par ./ clymo_alpha_1 .* (1 - exp(- clymo_alpha_1 * Mod_age2_duration_acrotelm_here)) + // part which is the acrotelm after the layer has been formed
    //  clymo_par ./ clymo_alpha_1 .* (1 - exp(- clymo_alpha_1 * Mod_age2_duration_catotelm_here)) .* exp(- clymo_alpha_2 * Mod_age2_duration_acrotelm_here)) .* // part of the layer which is in the catotelm after the layer has been formed
    //  exp(- clymo_alpha_1 * Mod_age2_duration_acrotelm_next) .* // decomposition in the acrotelm after the layer was formed
    //  exp(- clymo_alpha_2 * Mod_age2_duration_catotelm_next) ./ 1000; // decomposition in the catotelm until core collection
  }

}
model {

  // the overall mean accumulation rate
  // weak half normal prior
  alpha[1] ~ normal(0, 10*acc_mean_prior);

  // the gamma distributed innovations

  // prior parameterised by use set shape and the value of it's parent section
  alpha[2:K_tot] ~ gamma(acc_shape_adj, acc_shape_adj ./ alpha[parent[2:K_tot]]);

  // the memory parameters
  R ~ beta(mem_alpha, mem_beta);

  // the observation error inflation model

  if (inflate_errors){
    infl_shape ~ gamma(infl_shape_shape, infl_shape_shape / infl_shape_mean);
    infl_mean ~ normal(0, infl_sigma_sd);
    infl ~ gamma(infl_shape[1], infl_shape[1] / infl_mean[1]);
  }

  // bioturbation error model

  // parameters that are zero length do not get sampled
  L ~ gamma(L_prior_shape, L_rate);
  D ~ normal(0, D_prior_scale);

  // additional error in ages due to age-heterogeneity
  bt_error ~ gamma(n_ind, n_ind ./ age_het);

  if (model_bioturbation == 1){
    // bioturbation depth

    if (inflate_errors == 1 || model_displacement == 1){
      // the Likelihood of the data given the model
      bt_age ~ student_t(nu, Mod_age, obs_err_infl);
    } else {
      bt_age ~ student_t(nu, Mod_age, obs_err);
    }

  } else {
    if (inflate_errors == 1 || model_displacement == 1){
      // the Likelihood of the data given the model
      obs_age ~ student_t(nu, Mod_age, obs_err_infl);
    } else {
      obs_age ~ student_t(nu, Mod_age, obs_err);
    }
  }

  // Clymo model
  age0 ~ gamma(p1_age0, p2_age0 * p3_age0);
  clymo_par_raw ~ gamma(p1_clymo_par, p2_clymo_par * p3_clymo_par);
  clymo_par_memory ~ beta(clymo_par_memory_p1, clymo_par_memory_p2);
  clymo_alpha_1_raw ~ gamma(p1_clymo_alpha_1, p2_clymo_alpha_1 * p3_clymo_alpha_1);
  clymo_alpha_2_raw ~ gamma(p1_clymo_alpha_2, p2_clymo_alpha_2 * p3_clymo_alpha_2);
  layer_mass_shape ~ gamma(p1_layer_mass_shape, p2_layer_mass_shape * p3_layer_mass_shape);

  ac_age_raw ~ gamma(p1_ac_age, p2_ac_age * p3_ac_age);

  layer_mass ~ gamma(layer_mass_shape * p3_layer_mass_shape, layer_mass_shape * p3_layer_mass_shape ./ Mod_layer_mass_cumulative[index_has_mass_measurements]);

}

generated quantities {

  //vector[N2_c] layer_mass_rep;
  vector[N2_c] layer_cumulative_mass_rep;
  vector[N2_c] nmp_rep; // current (apparent) mass accumulation rate
  vector[N2_c] nmu_rep; // net mass uptake at the acrotelm-catotelm boundary (net mass uptake at the top of the acrotelm equals clymo_par)
  vector[N2_c] nmr_rep; // total C release by the entire peat core during a specified period

  // mass and cumulative mass
  for(n in 1:N2_c) {
    layer_cumulative_mass_rep[n] = gamma_rng(layer_mass_shape, layer_mass_shape ./ Mod_layer_mass_cumulative[n]);
  }
  //layer_cumulative_mass_rep = cumulative_mass0 + cumulative_sum(layer_mass_rep);

  {

    vector[N2_c] Mod_age2_upper_or = Mod_age2_all[1:N2_c];
    vector[N2_c] Mod_age2_lower_or = Mod_age2_all[2:(N2_c + 1)];
    vector[N2_c] fraction_mass_lost = (exp(- clymo_alpha_2 * (Mod_age2_duration_catotelm_next + Mod_age2_duration_catotelm_here)));

    // mass fluxes
    for(n in 1:N2_c) {
      if(Mod_age2_upper_or[n] < ac_age) { // need to decompose the layer down to the acrotelm-catotelm boundary
         fraction_mass_lost[n] = 1 / exp( - clymo_alpha_1 * (ac_age - Mod_age2_upper_or[n]));
      }
    }


    nmp_rep = Mod_layer_mass;
    nmu_rep = nmp_rep ./ fraction_mass_lost;

    // nmr_rep
    for(k in 1:N2_c) {

      // current time
      int N2_c_k = N2_c - k + 1;
      vector[N2_c_k] Mod_age2_upper_start_rep = Mod_age2_upper_or[k:N2_c] - Mod_age2_upper_or[k];
      vector[N2_c_k] Mod_age2_upper_end_rep = Mod_age2_upper_start_rep + Mod_age2_duration[k];

      // time information to reconstruct layer masses at the start of the current interval
      vector[N2_c_k] Mod_age2_duration_acrotelm_here_start_rep = rep_vector(ac_age, N2_c_k);
      vector[N2_c_k] Mod_age2_duration_catotelm_here_start_rep = rep_vector(Mod_age2_duration[k] - ac_age, N2_c_k);
      vector[N2_c_k] Mod_age2_duration_acrotelm_next_start_rep = rep_vector(0.0, N2_c_k);
      vector[N2_c_k] Mod_age2_duration_catotelm_next_start_rep = Mod_age2_upper_start_rep;
      if(ac_age > Mod_age2_duration[k]) { // AC boundary is outside the layer
        Mod_age2_duration_acrotelm_here_start_rep = rep_vector(Mod_age2_duration[k], N2_c_k);
        Mod_age2_duration_catotelm_here_start_rep = rep_vector(0.0, N2_c_k);
        Mod_age2_duration_catotelm_next_start_rep = Mod_age2_upper_start_rep - (ac_age - Mod_age2_duration_acrotelm_here_start_rep);
        for(m in 1:N2_c_k) {
          if(Mod_age2_duration_catotelm_next_start_rep[m] < 0.0) {
            Mod_age2_duration_catotelm_next_start_rep[m] = 0.0;
          }
        }
        Mod_age2_duration_acrotelm_next_start_rep = Mod_age2_upper_start_rep - Mod_age2_duration_catotelm_next_start_rep;
      }

      // time information to reconstruct layer masses at the end of the current interval
      vector[N2_c_k] Mod_age2_duration_acrotelm_here_end_rep = rep_vector(ac_age, N2_c_k);
      vector[N2_c_k] Mod_age2_duration_catotelm_here_end_rep = rep_vector(Mod_age2_duration[k] - ac_age, N2_c_k);
      vector[N2_c_k] Mod_age2_duration_acrotelm_next_end_rep = rep_vector(0.0, N2_c_k);
      vector[N2_c_k] Mod_age2_duration_catotelm_next_end_rep = Mod_age2_upper_end_rep;
      if(ac_age > Mod_age2_duration[k]) { // AC boundary is outside the layer
        Mod_age2_duration_acrotelm_here_end_rep = rep_vector(Mod_age2_duration[k], N2_c_k);
        Mod_age2_duration_catotelm_here_end_rep = rep_vector(0.0, N2_c_k);
        Mod_age2_duration_catotelm_next_end_rep = Mod_age2_upper_end_rep - (ac_age - Mod_age2_duration_acrotelm_here_end_rep);
        for(m in 1:N2_c_k) {
          if(Mod_age2_duration_catotelm_next_end_rep[m] < 0.0) {
            Mod_age2_duration_catotelm_next_end_rep[m] = 0.0;
          }
        }
        Mod_age2_duration_acrotelm_next_end_rep = Mod_age2_upper_end_rep - Mod_age2_duration_catotelm_next_end_rep;
      }

      // reconstruct layer masses
      vector[N2_c_k] ncu_1 = nmp_rep[k:N2_c] ./ (exp(- clymo_alpha_1 * Mod_age2_duration_acrotelm_next_end_rep - clymo_alpha_2 * Mod_age2_duration_catotelm_next_end_rep - clymo_alpha_1 * Mod_age2_duration_acrotelm_here_end_rep - clymo_alpha_2 * Mod_age2_duration_catotelm_here_end_rep));
      vector[N2_c_k] ncu_0 = nmp_rep[k:N2_c] ./ (exp(- clymo_alpha_1 * Mod_age2_duration_acrotelm_next_start_rep - clymo_alpha_2 * Mod_age2_duration_catotelm_next_start_rep - clymo_alpha_1 * Mod_age2_duration_acrotelm_here_start_rep - clymo_alpha_2 * Mod_age2_duration_catotelm_here_start_rep));
      nmr_rep[k] = sum(ncu_1 - ncu_0);
    }

  }

  // normalize all to their duration
  nmp_rep = nmp_rep ./ Mod_age2_duration;
  nmu_rep = nmu_rep ./ Mod_age2_duration;
  nmr_rep = nmr_rep ./ Mod_age2_duration;

  // new:start
  {

    vector[N2_c] Mod_age2_upper_or = Mod_age2_all[1:N2_c];
    vector[N2_c] Mod_age2_lower_or = Mod_age2_all[2:(N2_c + 1)];
    vector[N2_c] fraction_mass_lost = (exp(- clymo_alpha_2 * (Mod_age2_duration_catotelm_next + Mod_age2_duration_catotelm_here)));

    // mass fluxes
    for(n in 1:N2_c) {
      if(Mod_age2_upper_or[n] < ac_age) { // need to decompose the layer down to the acrotelm-catotelm boundary
         fraction_mass_lost[n] = 1 / exp( - clymo_alpha_1 * (ac_age - Mod_age2_upper_or[n]));
      }
    }

    nmp_rep = Mod_layer_mass;
    nmu_rep = nmp_rep ./ fraction_mass_lost;

    // duration each layer spent in the acrotelm and catotelm at the time of core collection
    vector[N2_c] Mod_age2_duration_acrotelm_tot_rep = Mod_age2_duration_acrotelm_next + Mod_age2_duration_acrotelm_here;
    vector[N2_c] Mod_age2_duration_catotelm_tot_rep = Mod_age2_duration_catotelm_next + Mod_age2_duration_catotelm_here;

    // nmr_rep
    for(k in 1:N2_c) {

      // at currently considered time step: age of the upper boundary at the start and end of the current interval
      int N2_c_k = N2_c - k + 1;
      vector[N2_c_k] Mod_age2_upper_start_rep = Mod_age2_upper_or[k:N2_c] - Mod_age2_upper_or[k];
      vector[N2_c_k] Mod_age2_upper_end_rep = Mod_age2_upper_start_rep + Mod_age2_duration[k];

      // time spent in acrotelm and catotelm while recomposing to the start of the current interval
      vector[N2_c_k] Mod_age2_duration_catotelm_start_rep = Mod_age2_duration_catotelm_tot_rep - Mod_age2_upper_start_rep;
      vector[N2_c_k] Mod_age2_duration_acrotelm_start_rep = rep_vector(0.0, N2_c_k);
      for(j in 1:N2_c_k) {
        if(Mod_age2_duration_catotelm_start_rep[j] < 0.0) { // layer also spent some time in the acrotelm
          Mod_age2_duration_catotelm_start_rep[j] = Mod_age2_duration_catotelm_tot_rep[j];
          Mod_age2_duration_acrotelm_start_rep[j] = Mod_age2_upper_start_rep[j] - Mod_age2_duration_catotelm_start_rep[j];
        }
      }

      // time spent in acrotelm and catotelm while recomposing to the end of the current interval
      vector[N2_c_k] Mod_age2_duration_catotelm_end_rep = Mod_age2_duration_catotelm_tot_rep - Mod_age2_upper_end_rep;
      vector[N2_c_k] Mod_age2_duration_acrotelm_end_rep = rep_vector(0.0, N2_c_k);
      for(j in 1:N2_c_k) {
        if(Mod_age2_duration_catotelm_end_rep[j] < 0.0) { // layer also spent some time in the acrotelm
          Mod_age2_duration_catotelm_end_rep[j] = Mod_age2_duration_catotelm_tot_rep[j];
          Mod_age2_duration_acrotelm_end_rep[j] = Mod_age2_upper_end_rep[j] - Mod_age2_duration_catotelm_end_rep[j];
        }
      }

      // reconstruct layer masses
      vector[N2_c_k] ncu_1 = nmp_rep[k:N2_c] ./ (exp(- clymo_alpha_1 * Mod_age2_duration_acrotelm_end_rep - clymo_alpha_2 * Mod_age2_duration_catotelm_end_rep));
      vector[N2_c_k] ncu_0 = nmp_rep[k:N2_c] ./ (exp(- clymo_alpha_1 * Mod_age2_duration_acrotelm_start_rep - clymo_alpha_2 * Mod_age2_duration_catotelm_start_rep));
      nmr_rep[k] = sum(ncu_1 - ncu_0);

    }

  }

  // normalize all to their duration
  nmp_rep = nmp_rep ./ Mod_age2_duration;
  nmu_rep = nmu_rep ./ Mod_age2_duration;
  nmr_rep = nmr_rep ./ Mod_age2_duration;

}

