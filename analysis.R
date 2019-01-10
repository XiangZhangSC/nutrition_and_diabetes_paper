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

# load metabolomics data
d <- read_rds("metabolomics_africa472_obesity300.rds")

# An adult cannot have BMI 12. 
d$BMI[which.min(d$BMI)] <- NA

## Table1
tab1 <- d %>% 
  select(sampleid:HbA1c_mmolmol) %>% 
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
            `sd_a1c_%` = sd(`HbA1c_%`, na.rm = TRUE))

a1c.mod <- stan_model("model.stan")

prep_data <- function(df, y_scale_factor) {
  
  # y observed and x observed
  df1 <- df %>% 
    filter(y_observed == 1, x_missing ==0)
  
  N1 = nrow(df1)
  y1 = as.array(df1$y)
  x1 = as.array(df1$HbA1c_mmolmol)
  race1 <- as.array(df1$race)
  age1 <- as.array(df1$age)
  sex1 <- as.array(df1$sex)
  
  # y observed and x missing
  df2 <- df %>% 
    filter(y_observed == 1, x_missing == 1)
  
  N2 = nrow(df2)
  y2 = as.array(df2$y)
  race2 <- as.array(df2$race)
  age2 <- as.array(df2$age)
  sex2 <- as.array(df2$sex)
  
  # y missing and x observed
  df3 <- df %>% 
    filter(y_missing == 1, x_missing == 0)
  
  N3 = nrow(df3)
  x3 = as.array(df3$HbA1c_mmolmol)
  race3 <- as.array(df3$race)
  age3 <- as.array(df3$age)
  sex3 <- as.array(df3$sex)
  
  # y missing and x missing
  df4 <- df %>% 
    filter(y_missing == 1, x_missing == 1)
  
  N4 = nrow(df4)
  race4 <- as.array(df4$race)
  age4 <- as.array(df4$age)
  sex4 <- as.array(df4$sex)
  
  # y censered and x observed
  df5 <- df %>% 
    filter(y_censered == 1, x_missing == 0)
  
  N5 <- nrow(df5)
  x5 <- as.array(df5$HbA1c_mmolmol)
  
  race5 <- as.array(df5$race)
  age5 <- as.array(df5$age)
  sex5 <- as.array(df5$sex)
  
  # y censered and x missing
  df6 <- df %>% 
    filter(y_censered == 1, x_missing == 1)
  
  N6 <- nrow(df6)
  
  race6 <- as.array(df6$race)
  age6 <- as.array(df6$age)
  sex6 <- as.array(df6$sex)
  
  list(
    N1 = N1, 
    y1 = y1 * y_scale_factor,
    x1 = x1, 
    race1 = race1, 
    age1 = age1,
    sex1 = sex1, 
    N2 = N2, 
    y2 = y2 * y_scale_factor,
    race2 = race2, 
    age2 = age2, 
    sex2 = sex2, 
    N3 = N3, 
    x3 = x3, 
    race3 = race3, 
    age3 = age3, 
    sex3 = sex3, 
    N4 = N4, 
    race4 = race4, 
    age4 = age4, 
    sex4 = sex4, 
    N5 = N5, 
    x5 = x5, 
    race5 = race5, 
    age5 = age5, 
    sex5 = sex5, 
    N6 = N6, 
    race6 = race6, 
    age6 = age6, 
    sex6 = sex6
  )
}

AAs <- d %>% 
  select(sampleid, age, gender, ethnicity, HbA1c_mmolmol, 
         Ala:Tyr) %>% 
  gather(what, y, Ala:Tyr) %>% 
  mutate(y_missing = ifelse(is.na(y), 1L, 0L), 
         y_censered = ifelse(y == 0, 1L, 0L), 
         y_observed = ifelse(y_missing == 0 & y_censered == 0, 1L, 0L), 
         x_missing = ifelse(is.na(HbA1c_mmolmol), 1L, 0L), 
         race = as.integer(ethnicity), 
         sex = as.integer(gender)) %>% 
  group_by(what) %>% 
  nest() %>% 
  mutate(stan_dat = map(data, prep_data, y_scale_factor = 1e3))

xiang_cluster <- create_cluster(8) %>% 
  cluster_assign_value("nightingale_model", a1c.mod) %>% 
  cluster_library("rstan") %>% 
  cluster_library("purrr") %>% 
  cluster_library("dplyr") %>% 
  cluster_library("tibble")

AAs_stan <- AAs %>% 
  partition(what, cluster = xiang_cluster) %>% 
  mutate(nightingale_sim = map(stan_dat, ~sampling(object = nightingale_model, 
                                                   data = .x, 
                                                   chains = 4, 
                                                   iter = 4000))) %>% 
  collect() %>% 
  as_tibble()


plot_associations <- function(df) {
  ggplot(df) + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
    geom_pointrangeh(aes(x = estimate, 
                         y = what,
                         xmin = conf.low, 
                         xmax = conf.high, 
                         shape = gender), 
                     position = position_dodgev(0.3)) + 
    scale_shape(solid = FALSE) + 
    facet_wrap(~term, ncol = 3) + 
    labs(x = "regression coefficient with 95% credible interval", 
         y = NULL) + 
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(legend.position = "bottom", 
          legend.title = element_blank())
}

AAs_res <- AAs_stan %>% 
  mutate(result = map(nightingale_sim, tidy, 
                      pars = c("b_EUR_m", "b_SUR_m", "b_GHA_m", 
                               "b_EUR_f", "b_SUR_f", "b_GHA_f"), 
                      conf.int = TRUE, 
                      conf.level = 0.95, 
                      conf.method = "HPDinterval")) %>% 
  select(what, result) %>% 
  ungroup() %>% 
  unnest(result) %>% 
  mutate(gender = str_sub(term, -1, -1), 
         gender = fct_recode(gender, 
                             "Female" = "f", 
                             "Male" = "m"), 
         term = str_sub(term, 3, -3),
         term = fct_recode(term, 
                           "European" = "EUR", 
                           "African Surinamese" = "SUR", 
                           "Ghanaian" = "GHA"
  ), 
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
ggsave("Figure1.jpeg", width = 178, height = 210, units = "mm", dpi = 300)

## KB

KBs <- d %>% 
  select(sampleid, age, gender, ethnicity, HbA1c_mmolmol, 
         AcAce:bOHBut) %>% 
  gather(what, y, AcAce:bOHBut) %>% 
  mutate(y_missing = ifelse(is.na(y), 1L, 0L), 
         y_censered = ifelse(y == 0, 1L, 0L), 
         y_observed = ifelse(y_missing == 0 & y_censered == 0, 1L, 0L), 
         x_missing = ifelse(is.na(HbA1c_mmolmol), 1L, 0L), 
         race = as.integer(ethnicity), 
         sex = as.integer(gender)) %>% 
  group_by(what) %>% 
  nest() %>% 
  mutate(stan_dat = map(data, prep_data, y_scale_factor = 1e3))

xiang_cluster <- create_cluster(2) %>% 
  cluster_assign_value("nightingale_model", a1c.mod) %>% 
  cluster_library("rstan") %>% 
  cluster_library("purrr") %>% 
  cluster_library("dplyr") %>% 
  cluster_library("tibble")

KBs_stan <- KBs %>% 
  partition(what, cluster = xiang_cluster) %>% 
  mutate(nightingale_sim = map(stan_dat, ~sampling(object = nightingale_model, 
                                                   data = .x, 
                                                   chains = 4, 
                                                   iter = 4000))) %>% 
  collect() %>% 
  as_tibble()

KBs_res <- KBs_stan %>% 
  mutate(result = map(nightingale_sim, tidy, 
                      pars = c("b_EUR_m", "b_SUR_m", "b_GHA_m", 
                               "b_EUR_f", "b_SUR_f", "b_GHA_f"), 
                      conf.int = TRUE, 
                      conf.level = 0.95, 
                      conf.method = "HPDinterval")) %>% 
  select(what, result) %>% 
  ungroup() %>% 
  unnest(result) %>% 
  mutate(gender = str_sub(term, -1, -1), 
         gender = fct_recode(gender, 
                             "Female" = "f", 
                             "Male" = "m"), 
         term = str_sub(term, 3, -3),
         term = fct_recode(term, 
                           "European" = "EUR", 
                           "African Surinamese" = "SUR", 
                           "Ghanaian" = "GHA"),
         what = factor(what, levels = c("AcAce", 
                                        'bOHBut'), 
                       labels = c("Acetoacetate", 
                                  "3-Hydroxybutyrate")), 
         what = fct_rev(what)) 

write_csv(KBs_res, "associations_ketobodies.csv")
plot_associations(KBs_res)
ggsave("Figure2.jpeg", width = 178, height = 80, units = "mm", dpi = 300)

## apoA1

apoA1_lipoproteins <- d %>% 
  select(sampleid, age, gender, ethnicity, 
         HbA1c_mmolmol,
         `XL-HDL-P`, 
         `L-HDL-P`, 
         `M-HDL-P`, 
         `S-HDL-P`) %>% 
  gather(what, y, `XL-HDL-P`:`S-HDL-P`) %>% 
  mutate(y_missing = ifelse(is.na(y), 1L, 0L), 
         y_censered = ifelse(y == 0, 1L, 0L), 
         y_observed = ifelse(y_missing == 0 & y_censered == 0, 1L, 0L), 
         x_missing = ifelse(is.na(HbA1c_mmolmol), 1L, 0L), 
         race = as.integer(ethnicity), 
         sex = as.integer(gender)) %>% 
  group_by(what) %>% 
  nest() %>% 
  mutate(stan_dat = map(data, prep_data, y_scale_factor = 1e6))

xiang_cluster <- create_cluster(4) %>% 
  cluster_assign_value("nightingale_model", a1c.mod) %>% 
  cluster_library("rstan") %>% 
  cluster_library("purrr") %>% 
  cluster_library("dplyr") %>% 
  cluster_library("tibble")

apoA1_lipoproteins_stan <- apoA1_lipoproteins %>% 
  partition(what, cluster = xiang_cluster) %>% 
  mutate(nightingale_sim = map(stan_dat, ~sampling(object = nightingale_model, 
                                                   data = .x, 
                                                   chains = 4, 
                                                   iter = 4000))) %>% 
  collect() %>% 
  as_tibble()

apoA1_lipoproteins_res <- apoA1_lipoproteins_stan %>% 
  mutate(result = map(nightingale_sim, tidy, 
                      pars = c("b_EUR_m", "b_SUR_m", "b_GHA_m", 
                               "b_EUR_f", "b_SUR_f", "b_GHA_f"), 
                      conf.int = TRUE, 
                      conf.level = 0.95, 
                      conf.method = "HPDinterval")) %>% 
  select(what, result) %>% 
  unnest(result) %>% 
  ungroup() %>% 
  mutate(gender = str_sub(term, -1, -1), 
         gender = fct_recode(gender, 
                             "Female" = "f", 
                             "Male" = "m"), 
         term = str_sub(term, 3, -3),
         term = fct_recode(term, 
                           "European" = "EUR", 
                           "African Surinamese" = "SUR", 
                           "Ghanaian" = "GHA"),
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
  select(sampleid, age, gender, ethnicity, 
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
  mutate(y_missing = ifelse(is.na(y), 1L, 0L), 
         y_censered = ifelse(y == 0, 1L, 0L), 
         y_observed = ifelse(y_missing == 0 & y_censered == 0, 1L, 0L), 
         x_missing = ifelse(is.na(HbA1c_mmolmol), 1L, 0L), 
         race = as.integer(ethnicity), 
         sex = as.integer(gender)) %>% 
  group_by(what) %>% 
  nest() %>% 
  mutate(stan_dat = map(data, prep_data, y_scale_factor = 1e9))

xiang_cluster <- create_cluster(10) %>% 
  cluster_assign_value("nightingale_model", a1c.mod) %>% 
  cluster_library("rstan") %>% 
  cluster_library("purrr") %>% 
  cluster_library("dplyr") %>% 
  cluster_library("tibble")

apoB_lipoproteins_stan <- apoB_lipoproteins %>% 
  partition(what, cluster = xiang_cluster) %>% 
  mutate(nightingale_sim = map(stan_dat, ~sampling(object = nightingale_model, 
                                                   data = .x, 
                                                   chains = 4, 
                                                   iter = 4000))) %>% 
  collect() %>% 
  as_tibble()

apoB_lipoproteins_res <- apoB_lipoproteins_stan %>% 
  mutate(result = map(nightingale_sim, tidy, 
                      pars = c("b_EUR_m", "b_SUR_m", "b_GHA_m", 
                               "b_EUR_f", "b_SUR_f", "b_GHA_f"), 
                      conf.int = TRUE, 
                      conf.level = 0.95, 
                      conf.method = "HPDinterval")) %>% 
  select(what, result) %>% 
  unnest(result) %>% 
  ungroup() %>% 
  mutate(gender = str_sub(term, -1, -1), 
         gender = fct_recode(gender, 
                             "Female" = "f", 
                             "Male" = "m"), 
         term = str_sub(term, 3, -3),
         term = fct_recode(term, 
                           "European" = "EUR", 
                           "African Surinamese" = "SUR", 
                           "Ghanaian" = "GHA"),
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
plot_associations(apoB_lipoproteins_res)



bEURm <- spread_draws(AAs_stan$nightingale_sim[[6]], b_EUR_m)
bEURf <- spread_draws(AAs_stan$nightingale_sim[[6]], b_EUR_f)
mean(bEURm$b_EUR_m > 0)
mean(bEURf$b_EUR_f > 0)

bGHAm <- spread_draws(AAs_stan$nightingale_sim[[6]], b_GHA_m)
bGHAf <- spread_draws(AAs_stan$nightingale_sim[[6]], b_GHA_f)
mean(bGHAm$b_GHA_m < 0)
mean(bGHAf$b_GHA_f < 0)

bSURm <- spread_draws(AAs_stan$nightingale_sim[[6]], b_SUR_m)
bSURf <- spread_draws(AAs_stan$nightingale_sim[[6]], b_SUR_f)
mean(bSURm$b_SUR_m < 0)
mean(bSURf$b_SUR_f < 0)

