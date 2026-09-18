# ---- Libraries ----
library(tidyverse)
library(conflicted)
library(ggpubr)
library(ggbeeswarm)
library(rstatix)
library(DHARMa)
library(ggeffects)
library(glmmTMB)
library(broom)
library(broom.mixed)
library(lmerTest)
library(MASS)
library(emmeans)
library(data.table)  # for fast operations
library(patchwork)
conflict_prefer_all("dplyr", quiet = TRUE)
conflict_prefer_all("lmerTest", quiet = TRUE)
set.seed(1444)
rm(list = ls())

abundances  <- read_tsv("~/PIMMSgit/data/phages/RPKM.tsv")
load("~/PIMMSgit/data/Rdata/PIMMs_vOTUs.Rdata")
load("~/PIMMSgit/data/Rdata/breadth_cutoff.Rdata")
load("~/PIMMSgit/data/Rdata/PIMMs_vOTUs_metadata.Rdata")

breadth_filtered <- breadth_cutoff %>% 
    filter(coverage >= 10) %>% #Genome-wide coverage
    pull(genome)

donor_shared <- 
    abundances %>%
    group_by(sample) %>%  
    mutate(tot_RPKM = sum(RPKM)) %>% ungroup() %>% 
    mutate(RPKM = RPKM / tot_RPKM * 1e6) %>% 
    filter(genome %in% PIMMs_vOTUs$representative) %>% 
    filter(genome %in% breadth_filtered) %>% 
    group_by(genome) %>% 
    mutate(donor_shared = if_else(any(timepoint == "D004"), "Donor-shared", "Recipient-unique")) %>% ungroup() %>% 
    filter(timepoint != "D004") %>% 
    mutate(timepoint = as.numeric(timepoint)) %>% 
    # mutate(timepoint = factor(timepoint, levels = c(-7, 0, 2, 4, 7, 14, 28))) %>% 
    group_by(genome, patient) %>% 
    mutate(method_shared = case_when(
        n_distinct(method) > 1 ~ "Method-shared",
        TRUE ~ "Method-unique"
        )) %>% ungroup()

donor_shared_metadata  <-   
    donor_shared %>% 
    mutate(genome = gsub("\\|\\|.*", "", genome)) %>% 
    left_join(PIMMs_vOTUs_metadata %>% select(-individual, -timepoint), by = c("genome" = "representative"))

induction  <- 
  donor_shared_metadata %>% 
  select(-method_shared) %>% #Redo the method_shared by section
    mutate(sections = case_when(
        timepoint < 2 ~ "Pre_transplant",
        timepoint < 14 ~ "First_week",
        TRUE ~ "Later"
    )) %>% 
    mutate(sections = factor(sections, levels = c("Pre_transplant", "First_week", "Later")))  %>% 
    filter(sections != "Later") %>% 
    group_by(genome, patient, sections) %>% 
    mutate(method_shared = if_else(n_distinct(method) > 1, "Method_shared", "Method_unique")) %>% ungroup() %>% 
    group_by(genome, patient, sections, method_shared, method, b_phylum, donor_shared, Pred, host_range, class) %>% 
    summarise(RPKM = mean(RPKM), .groups = "drop")  %>% 
    group_by(genome, patient) %>% 
    mutate(increase = case_when(
        mean(RPKM[sections == "First_week"]) > (mean(RPKM[sections == "Pre_transplant"]) + (mean(RPKM[sections == "Pre_transplant"])/10)) ~ "increase",
        mean(RPKM[sections == "Pre_transplant"]) > (mean(RPKM[sections == "First_week"]) + (mean(RPKM[sections == "First_week"])/10)) ~ "decrease",
        TRUE ~ "no_changes"
    )) %>% ungroup()
    
    
method_shared  <- 
    induction %>% 
    group_by(genome, patient) %>% 
    filter(any(method_shared == "Method_shared")) %>% ungroup() %>% 
    group_by(genome, patient) %>% 
    mutate(
            WGS_at_b = if_else(any(sections == "Pre_transplant" & method_shared == "Method_shared"), "Yes", "No"),
            VLP_after = if_else(any(sections == "First_week" & method_shared == "Method_shared"), "Yes", "No"),
            # collapse increase to ONE value per genome — pick the rule that's biologically right:
            presence = case_when(
                WGS_at_b == "Yes" & VLP_after == "Yes" ~ "Shared_both",
                WGS_at_b == "No" & VLP_after == "Yes" ~ "Shared_after",
                WGS_at_b == "Yes" & VLP_after == "No" ~ "Shared_pre",
                WGS_at_b == "No" & VLP_after == "No" ~ "Never_shared",
            )) %>% ungroup() %>% 
    group_by(genome, patient) %>% 
    mutate(sections_unique = case_when(
             presence == "Shared_pre" & n_distinct(sections) < 2 ~ "Yes",
             presence == "Shared_after" & n_distinct(sections) < 2 ~ "Yes",
             TRUE ~"No"
    )) %>% ungroup() %>% 
    group_by(genome, patient) %>% 
    mutate(porcodio = case_when(
        !str_detect(presence, "both") & sections_unique == "No" ~ method,
        TRUE ~ "porcodio"
    ))   %>% 
    mutate(porcodio = names(sort(table(method)))[2]) %>% ungroup() %>%
    mutate(sections_unique = case_when(
                presence == "Shared_after" & sections_unique == "No" & porcodio == "VLP" ~ "VLP_pre",
                presence == "Shared_after" & sections_unique == "No" & porcodio == "WGS" ~ "WGS_pre",
                presence == "Shared_pre" & sections_unique == "No" & porcodio == "VLP" ~ "VLP_after",
                presence == "Shared_pre" & sections_unique == "No" & porcodio == "WGS" ~ "WGS_after",
                presence == "Shared_after" & sections_unique == "Yes" ~ "None_pre",
                presence == "Shared_pre" & sections_unique == "Yes" ~ "None_after",
                TRUE ~ ""
    )) %>% 
    mutate(presence = paste(presence, sections_unique, sep = "_")) %>% select(-sections_unique, -porcodio)  %>% 
    mutate(method_increase = case_when(
            method == "VLP" & increase == "increase" ~ "VLP_increased",
            method == "WGS" & increase == "increase" ~ "WGS_increased",
            method == "VLP" & increase == "decrease" ~ "VLP_decreased",
            method == "WGS" & increase == "decrease" ~ "WGS_decreased")) %>% ungroup() %>%
    group_by(genome, patient) %>%
    mutate(method_increase = if (n() == 3) method_increase[sections == "Pre_transplant"][1] else method_increase) %>%
    ungroup()

long_data <- method_shared %>%
  filter(sections == "First_week", !is.na(method_increase), !is.na(b_phylum), !is.na(Pred)) %>%
  mutate(b_phylum = if_else(b_phylum %in% c("Bacillota", "Bacteroidota"), b_phylum, "Other")) %>% 
  distinct(genome, patient, method_increase, presence, donor_shared, b_phylum, class, host_range, Pred) %>%
  mutate(induced_binary = if_else(
    presence %in% c("Shared_after_WGS_pre"), 1L, 0L
  ))

table(long_data$induced_binary)
#######################################

model_full  <- glmmTMB(induced_binary ~ method_increase + donor_shared + b_phylum + class + Pred + (1|patient), 
family = binomial, data = long_data %>% 
  filter(str_detect(method_increase, "WGS")))

model_full$sdr$pdHess
any(is.nan(vcov(model_full)$cond))
diagnose(model_full)
df.residual(model_full)

testResiduals(simulateResiduals(model_full, n = 1000))

drop1(model_full, test = "Chisq")

emmeans(model_full, ~donor_shared, type = "response") %>%  pairs()
rbind(
  emmeans(model_full, ~ method_increase, type = "response") %>% pairs(),
  emmeans(model_full, ~ class,      type = "response") %>% pairs(),
  emmeans(model_full, ~ donor_shared,    type = "response") %>% pairs(),
  emmeans(model_full, ~ Pred,            type = "response") %>% pairs()
) %>% summary(adjust = "BH")


glimpse(long_data)


with(long_data, table(induced_binary, donor_shared, patient))
emmeans(model_full, ~ donor_shared, type = "response") %>% 
  pairs(infer = TRUE)    # adds the confidence interval
###########################################################################################

#####################################Method increase####################################

method_unique <- 
  induction %>% 
  group_by(genome, patient) %>% 
  filter(any(method_shared == "Method_unique")) %>% ungroup() %>% 
  mutate(location = case_when(
    sections == "Pre_transplant" & method == "WGS" ~ "WGS_pre",
    sections == "Pre_transplant" & method == "VLP" ~ "VLP_pre",
    sections == "First_week" & method == "WGS" ~ "WGS_after",
    sections == "First_week" & method == "VLP" ~ "VLP_after"
  )) %>% 
  mutate(method_increase = case_when(
            method == "VLP" & increase == "increase" ~ "VLP_increased",
            method == "WGS" & increase == "increase" ~ "WGS_increased",
            method == "VLP" & increase == "decrease" ~ "VLP_decreased",
            method == "WGS" & increase == "decrease" ~ "WGS_decreased"
            )) %>% ungroup() %>%
    group_by(genome, patient) %>%
    mutate(method_increase = if (n() == 3) method_increase[sections == "Pre_transplant"][1] else method_increase) %>%
    ungroup()  %>% 
    anti_join(method_shared, by = c("genome", "patient")) %>%
    group_by(patient, genome) %>% 
    mutate(presence = case_when(
        any(location == "WGS_pre") & any(location == "VLP_after") ~ "WGS_pre_VLP_after",
        any(location == "WGS_after") & any(location == "VLP_pre") ~ "VLP_pre_WGS_after",
        any(location == "WGS_pre") & any(location == "WGS_after") ~ "WGS_pre_WGS_after",
        any(location == "VLP_pre") & any(location == "VLP_after") ~ "VLP_pre_VLP_after"
    )) %>% ungroup() %>% select(-location) %>% 
    filter(!is.na(presence), sections == "First_week")

#####################################################################################################################################

model_data_unique  <- 
  method_unique %>% 
  filter(sections == "First_week", !is.na(Pred)) %>%
  mutate(b_phylum = if_else(b_phylum %in% c("Bacillota", "Bacteroidota"), b_phylum, "Other"))  %>% 
  distinct(genome, patient, presence, donor_shared, b_phylum, class, Pred) %>%
  mutate(induced_binary = if_else(
    presence %in% c("WGS_pre_VLP_after"), 1L, 0L
  ))

with(model_data_unique, table(presence))


#Add this to the long_data and do not include abundance increase###########################

long_no_abundance  <- long_data %>%  filter(str_detect(method_increase, "WGS")) %>% select(-method_increase)

combined_long  <- bind_rows(long_no_abundance, model_data_unique)

with(long_data, table(presence))

model_full  <- glmmTMB(induced_binary ~ donor_shared + b_phylum + Pred + class + host_range + (1 | patient), 
family = binomial, data = combined_long %>%
  filter(!presence %in% c("Shared_after_None_pre", 
  "Shared_after_VLP_pre", 
  "Shared_pre_VLP_after", 
  "VLP_pre_VLP_after", 
  "VLP_pre_WGS_after")) %>% #Remove the ones that could have not be induced
  drop_na(donor_shared, b_phylum, host_range, Pred, patient, induced_binary))

model_full$sdr$pdHess
any(is.nan(vcov(model_full)$cond))
diagnose(model_full)
df.residual(model_full)

testResiduals(simulateResiduals(model_full, n = 1000))
plot(simulateResiduals(model_full, n = 1000))

drop1(model_full, test = "Chisq")

emmeans(model_full, ~donor_shared, type = "response") %>% pairs()

model_clean  <- glmmTMB(induced_binary ~ donor_shared + Pred + b_phylum + class + (1 | patient), 
family = binomial, data = combined_long %>%
  filter(!presence %in% c("Shared_after_None_pre", 
  "Shared_after_VLP_pre", 
  "Shared_pre_VLP_after", 
  "VLP_pre_VLP_after", 
  "VLP_pre_WGS_after")) %>% 
  drop_na(donor_shared, b_phylum, host_range, Pred, patient, induced_binary))


model_clean$sdr$pdHess
any(is.nan(vcov(model_clean)$cond))
diagnose(model_clean)
df.residual(model_clean)

testDispersion(simulateResiduals(model_clean, n = 1000))
plot(simulateResiduals(model_clean, n = 1000))

drop1(model_clean, test = "Chisq")

AIC(model_full, model_clean)
BIC(model_full, model_clean)
anova(model_full, model_clean, test = "Chisq")

with(combined_long%>%
  filter(!presence %in% c("Shared_after_None_pre", 
  "Shared_after_VLP_pre", 
  "Shared_pre_VLP_after", 
  "VLP_pre_VLP_after", 
  "VLP_pre_WGS_after")) %>% 
  drop_na(donor_shared, b_phylum, host_range, Pred, patient, induced_binary), table(donor_shared, host_range))
##########################The differences we see die, for some reason after we get everything together################################
#Bayesan model#######################
library(brms)
conflict_prefer_all("brms", quiet = TRUE)
conflict_prefer_all("MASS", quiet = TRUE)
conflict_prefer_all("tidyr", quiet = TRUE)
conflict_prefer_all("DHARMa", quiet = TRUE)
conflict_prefer_all("data.table", quiet = TRUE)

d <- long_data %>%
  filter(!presence %in% c("Shared_after_None_pre","Shared_after_VLP_pre",
                          "Shared_pre_VLP_after","VLP_pre_VLP_after","VLP_pre_WGS_after")) %>%
  drop_na(donor_shared, b_phylum, host_range, Pred, class, patient, induced_binary)

m_bayes <- brm(
  induced_binary ~ donor_shared + b_phylum + Pred + class + host_range + (1 | patient),
  family = bernoulli(),
  data = d,
  prior = set_prior("normal(0, 2.5)", class = "b"),
  chains = 4, iter = 4000, cores = 4,
  control = list(adapt_delta = 0.95),   #you needed this last time for divergences
  seed = 1
)

summary(m_bayes)

m_with    <- m_bayes                              # host_range included
m_without <- update(m_bayes, . ~ . - host_range)  # host_range dropped

summary(m_without)
pp_check(m_bayes, ndraws = 100)
