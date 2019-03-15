data {
  int<lower=0> N;
  int<lower=0> N_y_obs;
  int<lower=0> N_y_mis;
  int<lower=0> N_y_cen;
  int<lower=1,upper=N> ii_y_obs[N_y_obs];
  int<lower=1,upper=N> ii_y_mis[N_y_mis];
  int<lower=1,upper=N> ii_y_cen[N_y_cen];
  vector<lower=0>[N_y_obs] y_obs;
  int<lower=0> N_a1c_obs;
  int<lower=0> N_a1c_mis;
  int<lower=1,upper=N> ii_a1c_obs[N_a1c_obs];
  int<lower=1,upper=N> ii_a1c_mis[N_a1c_mis];
  vector<lower=0>[N_a1c_obs] a1c_obs;
  int<lower=0> N_bmi_obs;
  int<lower=0> N_bmi_mis;
  int<lower=1,upper=N> ii_bmi_obs[N_bmi_obs];
  int<lower=1,upper=N> ii_bmi_mis[N_bmi_mis];
  vector<lower=0>[N_bmi_obs] bmi_obs;
  int<lower=1,upper=3> race[N];
  int<lower=1,upper=2> sex[N];
  int<lower=1,upper=2> diabetes[N];
  vector<lower=0>[N] age;
}
transformed data {
  real<lower=0> y_lod = min(y_obs);
}
parameters {
  vector<lower=0>[N_y_mis] y_mis;
  vector<lower=0,upper=y_lod>[N_y_cen] y_cen;
  vector<lower=0>[N_a1c_mis] a1c_mis;
  vector<lower=0>[N_bmi_mis] bmi_mis; 
  real bH;
  vector[3] bH_race;
  vector[2] bH_sex;
  vector[2] bH_diabetes;
  real bA;
  real bB;
  real a;
  vector[3] a_race;
  vector[2] a_sex;
  vector[2] a_diabetes;
  vector<lower=0>[3] mu_a1c;
  vector<lower=0>[3] sigma_a1c;
  vector<lower=0>[3] mu_bmi;
  vector<lower=0>[3] sigma_bmi;
  real<lower=0> sigma;
}
model {
  vector[N] y;
  vector[N] a1c;
  vector[N] bmi;
  vector[N] z_a1c;
  vector[N] z_bmi;
  vector[N] z_age;
  vector[N] A;
  vector[N] BH;
  
  y[ii_y_obs] = y_obs;
  y[ii_y_mis] = y_mis;
  y[ii_y_cen] = y_cen;
  
  a1c[ii_a1c_obs] = a1c_obs;
  a1c[ii_a1c_mis] = a1c_mis;
  
  bmi[ii_bmi_obs] = bmi_obs;
  bmi[ii_bmi_mis] = bmi_mis;
  
  bH ~ normal( 0, 0.5 ); 
  bH_race ~ normal( 0, 0.5 );
  bH_sex ~ normal( 0, 0.5 );
  bH_diabetes ~ normal( 0, 0.5 );
  bA ~ normal( 0, 0.5 );
  bB ~ normal( 0, 0.5 );
  a ~ normal( 0, 5 );
  a_race ~ normal( 0, 0.5 );
  a_sex ~ normal( 0, 0.5 );
  a_diabetes ~ normal( 0, 0.5 );
  mu_a1c ~ exponential( 0.025 );
  sigma_a1c ~ exponential( 0.2 );
  mu_bmi ~ exponential( 0.04 );
  sigma_bmi ~ exponential( 0.2 );
  sigma ~ exponential( 1 );
  
  a1c ~ normal( mu_a1c[race], sigma_a1c[race] );
  bmi ~ normal( mu_bmi[race], sigma_bmi[race] );
  
  z_a1c = ( a1c - 40 ) / 10;
  z_bmi = ( bmi - 30 ) / 5;
  z_age = ( age - 59 ) / 10; 
  
  A = a + a_race[race] + a_sex[sex] + a_diabetes[diabetes];
  BH = bH + bH_race[race] + bH_sex[sex] + bH_diabetes[diabetes];
  
  y ~ lognormal( A + rows_dot_product(BH, z_a1c) + bA * z_age + bB * z_bmi, sigma );
  
}
generated quantities {
  real BH[3,2,2];
  
  // slopes for ethnic groups
  for ( i in 1:3 ) {
    for ( j in 1:2 ) {
      for ( k in 1:2 ) {
        BH[i,j,k] = bH + bH_race[i] + bH_sex[j] + bH_diabetes[k];
      }
    }
  }
}
