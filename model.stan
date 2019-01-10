data {
  int<lower=0> N1;                               // # of y observed and x observed
  vector<lower=0>[N1] y1;                       
  vector[N1] x1;
  int race1[N1];
  vector[N1] age1;
  int sex1[N1];
  
  int<lower=0> N2;                               // # of y observed and x missing
  vector<lower=0>[N2] y2;
  int race2[N2];
  vector[N2] age2;
  int sex2[N2];
  
  int<lower=0> N3;                               // # of y missing and x observed
  vector<lower=0>[N3] x3;
  int race3[N3];
  vector[N3] age3;
  int sex3[N3];
  
  int<lower=0> N4;                               // # of y missing and x missing
  int race4[N4];
  vector[N4] age4;
  int sex4[N4];
  
  int<lower=0> N5;                               // # of y censered and x observed
  vector<lower=0>[N5] x5;                        
  int race5[N5];
  vector[N5] age5;
  int sex5[N5];
  
  int<lower=0> N6;                               // # of y censered and x missing
  int race6[N6];
  vector[N6] age6;
  int sex6[N6];
}
transformed data {
  vector[N1 + N2] y_obs;
  real<lower=0> y_lod;                           // lower limit of detection when y censered
  vector[N1 + N3 + N5] x_obs;
  vector[N1] zage1;
  vector[N2] zage2;
  vector[N3] zage3;
  vector[N4] zage4;
  vector[N5] zage5;
  vector[N6] zage6;
  
  x_obs = append_row(x1, append_row(x3, x5));
  
  y_obs = append_row(y1, y2);
  y_lod = min(y_obs);
  
  zage1 = (age1 - 59) / 10;
  zage2 = (age2 - 59) / 10;
  zage3 = (age3 - 59) / 10;
  zage4 = (age4 - 59) / 10;
  zage5 = (age5 - 59) / 10;
  zage6 = (age6 - 59) / 10;
}
parameters {
  vector<lower=0>[N2] x2;                          // x when y observed and x missing
  vector<lower=0>[N3] y3;                          // y when y missing and x observed
  vector<lower=0>[N4] x4;                          // x when y missing and x missing
  vector<lower=0>[N4] y4;                          // y when y missing and x missing
  vector<lower=0, upper=y_lod>[N5] y5;             // y when y censered and x observed
  vector<lower=0>[N6] x6;                          // x when y censered and x missing
  vector<lower=0, upper=y_lod>[N6] y6;             // y when y censored and x missing
  
  real bH;
  vector[3] bH_race;
  vector[2] bH_sex;
  real bA;
  real a;
  vector[3] a_race;
  vector[2] a_sex;
  vector<lower=0>[3] mu_x;
  vector<lower=0>[3] sigma_x;
  real<lower=0> sigma;
}
model {
  vector[N1] A1;
  vector[N2] A2;
  vector[N3] A3;
  vector[N4] A4;
  vector[N5] A5;
  vector[N6] A6;
  
  vector[N1] BH1;
  vector[N2] BH2;
  vector[N3] BH3;
  vector[N4] BH4;
  vector[N5] BH5;
  vector[N6] BH6;
  
  vector[N1] zx1;
  vector[N2] zx2;
  vector[N3] zx3;
  vector[N4] zx4;
  vector[N5] zx5;
  vector[N6] zx6;
  
  A1 = a + a_race[race1] + a_sex[sex1];
  A2 = a + a_race[race2] + a_sex[sex2];
  A3 = a + a_race[race3] + a_sex[sex3];
  A4 = a + a_race[race4] + a_sex[sex4];
  A5 = a + a_race[race5] + a_sex[sex5];
  A6 = a + a_race[race6] + a_sex[sex6];
  
  BH1 = bH + bH_race[race1] + bH_sex[sex1];
  BH2 = bH + bH_race[race2] + bH_sex[sex2];
  BH3 = bH + bH_race[race3] + bH_sex[sex3];
  BH4 = bH + bH_race[race4] + bH_sex[sex4];
  BH5 = bH + bH_race[race5] + bH_sex[sex5];
  BH6 = bH + bH_race[race6] + bH_sex[sex6];
  
  zx1 = (x1 - 40) / 10;
  zx2 = (x2 - 40) / 10;
  zx3 = (x3 - 40) / 10;
  zx4 = (x4 - 40) / 10;
  zx5 = (x5 - 40) / 10;
  zx6 = (x6 - 40) / 10;
  
  // prior
  bH ~ normal( 0, 10 ); 
  bH_race ~ normal( 0, 0.1 );
  bH_sex ~ normal( 0, 0.1 );
  bA ~ normal( 0, 10 );
  a ~ normal( 0, 10 );
  a_race ~ normal( 0, 0.1 );
  a_sex ~ normal( 0, 0.1 );
  mu_x ~ exponential( 0.025 );
  sigma_x ~ exponential( 0.2 ); 
  sigma ~ exponential( 1 );

  // imputation //
  
  // impute A1c
  x2 ~ normal( mu_x[race2], sigma_x[race2] );
  x4 ~ normal( mu_x[race4], sigma_x[race4] );
  x6 ~ normal( mu_x[race6], sigma_x[race6] );
  
  // impute y
  
  y3 ~ lognormal( A3 + rows_dot_product(BH3, zx3) + bA * zage3, sigma );
  
  y4 ~ lognormal( A4 + rows_dot_product(BH4, zx4) + bA * zage4, sigma );
  
  y5 ~ lognormal( A5 + rows_dot_product(BH5, zx5) + bA * zage5, sigma );
  
  y6 ~ lognormal( A6 + rows_dot_product(BH6, zx6) + bA * zage6, sigma );
  
  // likelihood //

  // observed A1c
  x1 ~ normal( mu_x[race1], sigma_x[race1] );
  x3 ~ normal( mu_x[race3], sigma_x[race3] );
  x5 ~ normal( mu_x[race5], sigma_x[race5] );
  
  // observed y
  y1 ~ lognormal( A1 + rows_dot_product(BH1, zx1) + bA * zage1, sigma );
  
  y2 ~ lognormal( A2 + rows_dot_product(BH2, zx2) + bA * zage2, sigma );
}
generated quantities {
  real b_EUR_m;
  real b_SUR_m;
  real b_GHA_m;
  real b_EUR_f;
  real b_SUR_f;
  real b_GHA_f;
  
  // slopes for ethnic groups
  // male
  b_EUR_m = bH + bH_race[1] + bH_sex[1];
  b_SUR_m = bH + bH_race[2] + bH_sex[1];
  b_GHA_m = bH + bH_race[3] + bH_sex[1];
  
  // female
  b_EUR_f = bH + bH_race[1] + bH_sex[2];
  b_SUR_f = bH + bH_race[2] + bH_sex[2];
  b_GHA_f = bH + bH_race[3] + bH_sex[2];
}
