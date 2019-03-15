library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(rstan)
library(multidplyr)
library(broom)
library(ggstance)
library(forcats)
library(stringr)
library(tidybayes)
library(gmodels)

# load metabolomics data
d <- read_rds("metabolomics_africa472_obesity300.rds")

d$diabetes <- ifelse(d$diabetes == 1, "Diabetes", "Prediabetes")
d$diabetes <- factor(d$diabetes, levels = c("Prediabetes", "Diabetes"))

# An adult cannot have BMI 12. 
d$BMI[which.min(d$BMI)] <- NA

## Table1
d1 <- d %>% 
  select(sampleid:HbA1c_mmolmol, diabetes)

test1 <- kruskal.test(age ~ ethnicity, data = d1)
test2 <- kruskal.test(BMI ~ ethnicity, data = d1)
test3 <- kruskal.test(waistcircum_cm ~ ethnicity, data = d1)
test4 <- kruskal.test(HbA1c_mmolmol ~ ethnicity, data = d1)

CrossTable(d1$ethnicity, d$gender, fisher = TRUE, chisq = TRUE, expected = TRUE, prop.c = FALSE, prop.t = FALSE, prop.chisq = FALSE, sresid = TRUE, format = "SPSS")
CrossTable(d1$ethnicity, d$diabetes, fisher = TRUE, chisq = TRUE, expected = TRUE, prop.c = FALSE, prop.t = FALSE, prop.chisq = FALSE, sresid = TRUE, format = "SPSS")

tab1 <- d1 %>% 
  mutate(`HbA1c_%` = (0.09148 * HbA1c_mmolmol) + 2.152) %>% 
  group_by(ethnicity) %>% 
  summarize(N = n(), 
            Female = mean(gender == "female"), 
            m_Age = mean(age), 
            sd_Age = sd(age), 
            m_BMI = mean(BMI, na.rm = TRUE), 
            sd_BMI = sd(BMI, na.rm = TRUE), 
            m_WC = mean(waistcircum_cm, na.rm = TRUE), 
            sd_WC = sd(waistcircum_cm, na.rm = TRUE), 
            m_a1c = mean(HbA1c_mmolmol, na.rm = TRUE), 
            sd_a1c = sd(HbA1c_mmolmol, na.rm = TRUE), 
            `m_a1c_%` = mean(`HbA1c_%`, na.rm = TRUE), 
            `sd_a1c_%` = sd(`HbA1c_%`, na.rm = TRUE), 
            N_DM = sum(diabetes == "Diabetes"), 
            P_DM = mean(diabetes == "Diabetes"))

a1c.mod <- stan_model("model.stan")

prep_data <- function(df, y_scale_factor) {
  
  N <- nrow(df)
  N_y_obs <- sum(df$y_observed)
  N_y_mis <- sum(df$y_missing)
  N_y_cen <- sum(df$y_censered)
  ii_y_obs <- which(df$y_observed == 1)
  ii_y_mis <- which(df$y_missing == 1)
  ii_y_cen <- which(df$y_censered == 1)
  y_obs <- df$y[ii_y_obs] * y_scale_factor
  N_a1c_obs <- sum(df$a1c_observed)
  N_a1c_mis <- sum(df$a1c_missing)
  ii_a1c_obs <- which(df$a1c_observed == 1)
  ii_a1c_mis <- which(df$a1c_missing == 1)
  a1c_obs <- df$HbA1c_mmolmol[ii_a1c_obs]
  N_bmi_obs <- sum(df$bmi_observed)
  N_bmi_mis <- sum(df$bmi_missing)
  ii_bmi_obs <- which(df$bmi_observed == 1)
  ii_bmi_mis <- which(df$bmi_missing == 1)
  bmi_obs <- df$BMI[ii_bmi_obs]
  race <- df$race
  sex <- df$sex
  diabetes <- df$diabetes
  age <- df$age
  
  list(
    N = N, 
    N_y_obs = N_y_obs, 
    N_y_mis = N_y_mis, 
    N_y_cen = N_y_cen, 
    ii_y_obs = ii_y_obs, 
    ii_y_mis = as.array(ii_y_mis), 
    ii_y_cen = as.array(ii_y_cen), 
    y_obs = y_obs, 
    N_a1c_obs = N_a1c_obs, 
    N_a1c_mis = N_a1c_mis, 
    ii_a1c_obs = ii_a1c_obs, 
    ii_a1c_mis = ii_a1c_mis, 
    a1c_obs = a1c_obs, 
    N_bmi_obs = N_bmi_obs, 
    N_bmi_mis = N_bmi_mis, 
    ii_bmi_obs = ii_bmi_obs, 
    ii_bmi_mis = ii_bmi_mis, 
    bmi_obs = bmi_obs, 
    race = race, 
    sex = sex, 
    diabetes = diabetes, 
    age = age
  )
}

plot_associations <- function(df) {
  ggplot(df) + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) + 
    geom_pointrangeh(aes(x = estimate, 
                         y = what,
                         xmin = conf.low, 
                         xmax = conf.high, 
                         shape = gender), 
                     position = position_dodgev(0.4)) + 
    scale_shape(solid = FALSE) + 
    facet_grid(diabetes~ethnicity) + 
    labs(x = "regression coefficient with 95% credible interval", 
         y = NULL) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(legend.position = "bottom", 
          legend.title = element_blank())
}

preprocess_dat <- function(df) {
  df %>% 
    mutate(y_missing = ifelse(is.na(y), 1L, 0L), 
           y_censered = ifelse(y_missing == 0 & y == 0, 1L, 0L), 
           y_observed = ifelse(y_missing == 0 & y_censered == 0, 1L, 0L), 
           a1c_observed = ifelse(!is.na(HbA1c_mmolmol), 1L, 0L), 
           a1c_missing = ifelse(is.na(HbA1c_mmolmol), 1L, 0L), 
           bmi_observed = ifelse(!is.na(BMI), 1L, 0L), 
           bmi_missing = ifelse(is.na(BMI), 1L, 0L), 
           race = as.integer(ethnicity), 
           sex = as.integer(gender), 
           diabetes = as.integer(diabetes))
}

AAs <- d %>% 
  select(sampleid, age, gender, ethnicity, HbA1c_mmolmol, diabetes, BMI, 
         Ala:Tyr) %>% 
  gather(what, y, Ala:Tyr) %>% 
  preprocess_dat() %>% 
  group_by(what) %>% 
  nest() %>% 
  mutate(stan_dat = map(data, prep_data, y_scale_factor = 1e3))

AAs_stan <- AAs %>% 
  mutate(nightingale_sim = map(stan_dat, ~sampling(a1c.mod, data = .x, chains = 4, cores = 4, 
                                                   control = list(max_treedepth = 11))))

#xiang_cluster <- create_cluster(8) %>% 
#  cluster_assign_value("nightingale_model", a1c.mod) %>% 
#  cluster_library("rstan") %>% 
#  cluster_library("purrr") %>% 
#  cluster_library("dplyr") %>% 
#  cluster_library("tibble")

#AAs_stan <- AAs %>% 
#  partition(what, cluster = xiang_cluster) %>% 
#  mutate(nightingale_sim = map(stan_dat, ~sampling(object = nightingale_model, 
#                                                   data = .x, 
#                                                   chains = 4, 
#                                                   iter = 4000))) %>% 
#  collect() %>% 
#  as_tibble()

write_rds(AAs_stan, "stan_out_AAs.rds")

AAs_res <- AAs_stan %>% 
  mutate(result = map(nightingale_sim, tidy, 
                      pars = "BH", 
                      conf.int = TRUE, 
                      conf.level = 0.95, 
                      conf.method = "HPDinterval")) %>% 
  select(what, result) %>% 
  ungroup() %>% 
  unnest(result) %>% 
  mutate(ethnicity = str_sub(term, 4, 4), 
         ethnicity = fct_recode(ethnicity, 
                                "European" = "1", 
                                "African Surinamese" = "2", 
                                "Ghanaian" = "3"), 
         gender = str_sub(term, 6, 6), 
         gender = fct_recode(gender, 
                             "Male" = "1", 
                             "Female" = "2"), 
         diabetes = str_sub(term, -2, -2),
         diabetes = fct_recode(diabetes, 
                               "Prediabetes" = "1", 
                               "Diabetes" = "2"),
  what = factor(what, 
                levels = c("Ile", "Leu", "Val", 
                           "Tyr", "Phe",
                           "Gln", "Ala", "His"),
                labels = c("Isoleucine", "Leucine", "Valine", 
                           "Tyrosine", "Phenylalanine", 
                           "Glutamine", "Alanine", "Histidine")), 
  what = fct_rev(what)) 

write_csv(AAs_res, "associations_amino_acids.csv")
plot_associations(AAs_res)
ggsave("Figure1.jpeg", width = 178, height = 250, units = "mm", dpi = 300)

## KB

KBs <- d %>% 
  select(sampleid, age, gender, ethnicity, HbA1c_mmolmol, diabetes, BMI, 
         AcAce:bOHBut) %>% 
  gather(what, y, AcAce:bOHBut) %>% 
  preprocess_dat() %>% 
  group_by(what) %>% 
  nest() %>% 
  mutate(stan_dat = map(data, prep_data, y_scale_factor = 1e3))

#xiang_cluster <- create_cluster(2) %>% 
#  cluster_assign_value("nightingale_model", a1c.mod) %>% 
#  cluster_library("rstan") %>% 
#  cluster_library("purrr") %>% 
#  cluster_library("dplyr") %>% 
#  cluster_library("tibble")

#KBs_stan <- KBs %>% 
#  partition(what, cluster = xiang_cluster) %>% 
#  mutate(nightingale_sim = map(stan_dat, ~sampling(object = nightingale_model, 
#                                                   data = .x, 
#                                                   chains = 4, 
#                                                   iter = 4000))) %>% 
#  collect() %>% 
#  as_tibble()


KBs_stan <- KBs %>% 
  mutate(nightingale_sim = map(stan_dat, ~sampling(a1c.mod, data = .x, chains = 4, cores = 4, 
                                                   control = list(max_treedepth = 11))))

write_rds(KBs_stan, "stan_out_KBs.rds")

KBs_res <- KBs_stan %>% 
  mutate(result = map(nightingale_sim, tidy, 
                      pars = "BH", 
                      conf.int = TRUE, 
                      conf.level = 0.95, 
                      conf.method = "HPDinterval")) %>% 
  select(what, result) %>% 
  ungroup() %>% 
  unnest(result) %>% 
  mutate(ethnicity = str_sub(term, 4, 4), 
         ethnicity = fct_recode(ethnicity, 
                                "European" = "1", 
                                "African Surinamese" = "2", 
                                "Ghanaian" = "3"), 
         gender = str_sub(term, 6, 6), 
         gender = fct_recode(gender, 
                             "Male" = "1", 
                             "Female" = "2"), 
         diabetes = str_sub(term, -2, -2),
         diabetes = fct_recode(diabetes, 
                           "Prediabetes" = "1", 
                           "Diabetes" = "2"),
         what = factor(what, levels = c("AcAce", 
                                        'bOHBut'), 
                       labels = c("Acetoacetate", 
                                  "3-Hydroxybutyrate")), 
         what = fct_rev(what)) 

write_csv(KBs_res, "associations_ketobodies.csv")
plot_associations(KBs_res)
ggsave("Figure2.jpeg", width = 178, height = 100, units = "mm", dpi = 300)

## apoA1

apoA1_lipoproteins <- d %>% 
  select(sampleid, age, gender, ethnicity, diabetes, BMI, 
         HbA1c_mmolmol,
         `XL-HDL-P`, 
         `L-HDL-P`, 
         `M-HDL-P`, 
         `S-HDL-P`) %>% 
  gather(what, y, `XL-HDL-P`:`S-HDL-P`) %>% 
  preprocess_dat() %>% 
  group_by(what) %>% 
  nest() %>% 
  mutate(stan_dat = map(data, prep_data, y_scale_factor = 1e6))

#xiang_cluster <- create_cluster(4) %>% 
#  cluster_assign_value("nightingale_model", a1c.mod) %>% 
#  cluster_library("rstan") %>% 
#  cluster_library("purrr") %>% 
#  cluster_library("dplyr") %>% 
#  cluster_library("tibble")

#apoA1_lipoproteins_stan <- apoA1_lipoproteins %>% 
#  partition(what, cluster = xiang_cluster) %>% 
#  mutate(nightingale_sim = map(stan_dat, ~sampling(object = nightingale_model, 
#                                                   data = .x, 
#                                                   chains = 4, 
#                                                   iter = 4000))) %>% 
#  collect() %>% 
#  as_tibble()

apoA1_lipoproteins_stan <- apoA1_lipoproteins %>% 
  mutate(nightingale_sim = map(stan_dat, ~sampling(a1c.mod, data = .x, chains = 4, cores = 4, 
                                                   control = list(max_treedepth = 11))))

write_rds(apoA1_lipoproteins_stan, "stan_out_apoA1.rds")

apoA1_lipoproteins_res <- apoA1_lipoproteins_stan %>% 
  mutate(result = map(nightingale_sim, tidy, 
                      pars = "BH", 
                      conf.int = TRUE, 
                      conf.level = 0.95, 
                      conf.method = "HPDinterval")) %>% 
  select(what, result) %>% 
  ungroup() %>% 
  unnest(result) %>% 
  mutate(ethnicity = str_sub(term, 4, 4), 
         ethnicity = fct_recode(ethnicity, 
                                "European" = "1", 
                                "African Surinamese" = "2", 
                                "Ghanaian" = "3"), 
         gender = str_sub(term, 6, 6), 
         gender = fct_recode(gender, 
                             "Male" = "1", 
                             "Female" = "2"), 
         diabetes = str_sub(term, -2, -2),
         diabetes = fct_recode(diabetes, 
                               "Prediabetes" = "1", 
                               "Diabetes" = "2"),
         what = factor(what, levels = c("XL-HDL-P", 
                                        'L-HDL-P', 
                                        "M-HDL-P", 
                                        "S-HDL-P"), 
                       labels = c("Very large HDL", 
                                  "Large HDL", 
                                  "Medium HDL", 
                                  "Small HDL")), 
         what = fct_rev(what)) 

write_csv(apoA1_lipoproteins_res, "associations_apoA1.csv")
plot_associations(apoA1_lipoproteins_res)
ggsave("Figure3.jpeg", width = 178, height = 160, units = "mm", dpi = 300)

## apoB

apoB_lipoproteins <- d %>% 
  select(sampleid, age, gender, ethnicity, diabetes, BMI, 
         HbA1c_mmolmol,
         `XXL-VLDL-P`, 
         `XL-VLDL-P`, 
         `L-VLDL-P`, 
         `M-VLDL-P`, 
         `S-VLDL-P`, 
         `XS-VLDL-P`, 
         `IDL-P`, 
         `L-LDL-P`, 
         `M-LDL-P`, 
         `S-LDL-P`) %>% 
  gather(what, y, `XXL-VLDL-P`:`S-LDL-P`) %>% 
  preprocess_dat() %>% 
  group_by(what) %>% 
  nest() %>% 
  mutate(stan_dat = map(data, prep_data, y_scale_factor = 1e9))

#xiang_cluster <- create_cluster(10) %>% 
#  cluster_assign_value("nightingale_model", a1c.mod) %>% 
#  cluster_library("rstan") %>% 
#  cluster_library("purrr") %>% 
#  cluster_library("dplyr") %>% 
#  cluster_library("tibble")

#apoB_lipoproteins_stan <- apoB_lipoproteins %>% 
#  partition(what, cluster = xiang_cluster) %>% 
#  mutate(nightingale_sim = map(stan_dat, ~sampling(object = nightingale_model, 
#                                                   data = .x, 
#                                                   chains = 4, 
#                                                   iter = 4000))) %>% 
#  collect() %>% 
#  as_tibble()

apoB_lipoproteins_stan <- apoB_lipoproteins %>% 
  mutate(nightingale_sim = map(stan_dat, ~sampling(a1c.mod, data = .x, chains = 4, cores = 4, 
                                                   control = list(max_treedepth = 11))))

write_rds(apoB_lipoproteins_stan, "stan_out_apoB.rds")

apoB_lipoproteins_res <- apoB_lipoproteins_stan %>% 
  mutate(result = map(nightingale_sim, tidy, 
                      pars = "BH", 
                      conf.int = TRUE, 
                      conf.level = 0.95, 
                      conf.method = "HPDinterval")) %>% 
  select(what, result) %>% 
  ungroup() %>% 
  unnest(result) %>% 
  mutate(ethnicity = str_sub(term, 4, 4), 
         ethnicity = fct_recode(ethnicity, 
                                "European" = "1", 
                                "African Surinamese" = "2", 
                                "Ghanaian" = "3"), 
         gender = str_sub(term, 6, 6), 
         gender = fct_recode(gender, 
                             "Male" = "1", 
                             "Female" = "2"), 
         diabetes = str_sub(term, -2, -2),
         diabetes = fct_recode(diabetes, 
                               "Prediabetes" = "1", 
                               "Diabetes" = "2"),
         what = factor(what, levels = c("XXL-VLDL-P", 
                                        'XL-VLDL-P', 
                                        "L-VLDL-P", 
                                        "M-VLDL-P", 
                                        "S-VLDL-P", 
                                        "XS-VLDL-P", 
                                        "IDL-P", 
                                        "L-LDL-P", 
                                        "M-LDL-P", 
                                        "S-LDL-P"), 
                       labels = c("Extreamly large VLDL", 
                                  "Very large VLDL", 
                                  "Large VLDL", 
                                  "Medium VLDL", 
                                  "Small VLDL", 
                                  "Very small VLDL", 
                                  "IDL", 
                                  "Large LDL", 
                                  "Medium LDL", 
                                  "Small LDL")), 
         what = fct_rev(what)) 

write_csv(apoB_lipoproteins_res, "associations_apoB.csv")

apoB_lipoproteins_res %>% 
  filter(what %in% c("Extreamly large VLDL", 
                     "Very large VLDL", 
                     "Large VLDL", 
                     "Medium VLDL", 
                     "Small VLDL", 
                     "Very small VLDL")) %>% 
  plot_associations()

apoB_lipoproteins_res %>% 
  filter(what %in% c("Large LDL", 
                     "Medium LDL", 
                     "Small LDL")) %>% 
  plot_associations()

ggsave("Figure4.jpeg", width = 178, height = 250, units = "mm", dpi = 300)
