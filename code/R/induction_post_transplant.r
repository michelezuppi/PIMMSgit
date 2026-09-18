source("~/PIMMSgit/data/Rcode/libraries_themes.r")
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


induction_later  <- 
  donor_shared_metadata %>% 
  select(-method_shared) %>% #Redo the method_shared by section
    mutate(sections = case_when(
        timepoint < 2 ~ "Pre_transplant",
        timepoint < 14 ~ "First_week",
        TRUE ~ "Later"
    )) %>% 
    mutate(sections = factor(sections, levels = c("Pre_transplant", "First_week", "Later")))  %>% 
    filter(sections != "Pre_transplant") %>% 
    group_by(genome, patient, sections) %>% 
    mutate(method_shared = if_else(n_distinct(method) > 1, "Method_shared", "Method_unique")) %>% ungroup() %>% 
    group_by(genome, patient, sections, method_shared, method, b_phylum, donor_shared, Pred, host_range, class) %>% 
    summarise(RPKM = mean(RPKM), .groups = "drop")  %>% 
    group_by(genome, patient) %>% 
    mutate(increase = case_when(
        mean(RPKM[sections == "Later"]) > (mean(RPKM[sections == "First_week"]) + (mean(RPKM[sections == "First_week"])/10)) ~ "increase",
        mean(RPKM[sections == "First_week"]) > (mean(RPKM[sections == "Later"]) + (mean(RPKM[sections == "Later"])/10)) ~ "decrease",
        TRUE ~ "no_changes"
    )) %>% ungroup()

method_shared_first_week  <- 
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
    ungroup() %>% 
    filter(sections == "First_week")


method_shared_later  <- 
    induction_later %>% 
    group_by(genome, patient) %>% 
    filter(any(method_shared == "Method_shared")) %>% ungroup() %>% 
    group_by(genome, patient) %>% 
    mutate(
            WGS_at_b = if_else(any(sections == "First_week" & method_shared == "Method_shared"), "Yes", "No"),
            VLP_after = if_else(any(sections == "Later" & method_shared == "Method_shared"), "Yes", "No"),
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
    filter(!is.na(method_increase)) %>% 
    group_by(genome, patient) %>%
            mutate(
            mode_mi = names(sort(table(method_increase), decreasing = TRUE))[1]
        ) %>% 
    mutate(method_increase = if_else(rep(n() == 3, n()), mode_mi, method_increase)) %>% ungroup() %>%
    select(-mode_mi) %>% 
    filter(sections == "Later")

method_unique_first_week <- 
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
    anti_join(method_shared_first_week, by = c("genome", "patient")) %>%
    group_by(patient, genome) %>% 
    mutate(presence = case_when(
        any(location == "WGS_pre") & any(location == "VLP_after") ~ "WGS_pre_VLP_after",
        any(location == "WGS_after") & any(location == "VLP_pre") ~ "VLP_pre_WGS_after",
        any(location == "WGS_pre") & any(location == "WGS_after") ~ "WGS_pre_WGS_after",
        any(location == "VLP_pre") & any(location == "VLP_after") ~ "VLP_pre_VLP_after"
    )) %>% ungroup() %>% select(-location) %>% 
    filter(!is.na(presence), sections == "First_week")

method_unique_later <- 
  induction_later %>% 
  group_by(genome, patient) %>% 
  filter(any(method_shared == "Method_unique")) %>% ungroup() %>% 
  mutate(location = case_when(
    sections == "First_week" & method == "WGS" ~ "WGS_pre",
    sections == "First_week" & method == "VLP" ~ "VLP_pre",
    sections == "Later" & method == "WGS" ~ "WGS_after",
    sections == "Later" & method == "VLP" ~ "VLP_after"
  )) %>% 
  mutate(method_increase = case_when(
            method == "VLP" & increase == "increase" ~ "VLP_increased",
            method == "WGS" & increase == "increase" ~ "WGS_increased",
            method == "VLP" & increase == "decrease" ~ "VLP_decreased",
            method == "WGS" & increase == "decrease" ~ "WGS_decreased"
            )) %>% ungroup() %>%
    anti_join(method_shared_later, by = c("genome", "patient")) %>%
    group_by(patient, genome) %>% 
    mutate(presence = case_when(
        any(location == "WGS_pre") & any(location == "VLP_after") ~ "WGS_pre_VLP_after",
        any(location == "WGS_after") & any(location == "VLP_pre") ~ "VLP_pre_WGS_after",
        any(location == "WGS_pre") & any(location == "WGS_after") ~ "WGS_pre_WGS_after",
        any(location == "VLP_pre") & any(location == "VLP_after") ~ "VLP_pre_VLP_after"
    )) %>% ungroup() %>% select(-location) %>% 
    filter(!is.na(presence), sections == "Later")

combined_method  <- bind_rows(method_shared_first_week, method_shared_later, method_unique_first_week, method_unique_later)

long_data <- combined_method %>%
  mutate(b_phylum = replace_na(b_phylum, "Unknown")) %>% 
  mutate(host_range = replace_na(host_range, "Unknown")) %>% 
  mutate(class = if_else(class == "" | is.na(class), "unplaced", class)) %>% ungroup() %>% 
  filter(sections != "Pre_transplant", !is.na(Pred)) %>%
  mutate(b_phylum = if_else(b_phylum %in% c("Bacillota", "Bacteroidota"), b_phylum, "Other")) %>% 
  distinct(genome, patient, sections, method_increase, presence, donor_shared, b_phylum, host_range, Pred, class) %>%
  mutate(induced_binary = if_else(
    presence %in% c("Shared_after_WGS_pre", "WGS_pre_VLP_after"), 1L, 0L
  )) 

model_full  <- glmmTMB(induced_binary ~ sections + b_phylum + donor_shared + class + Pred + (1 | patient), 
family = binomial, data = long_data %>% 
  filter(!presence %in% c("Shared_after_None_pre", "Shared_after_VLP_pre", "Shared_pre_VLP_after", "VLP_pre_VLP_after", "VLP_pre_WGS_after")))

model_full$sdr$pdHess
any(is.nan(vcov(model_full)$cond))
diagnose(model_full)
df.residual(model_full)

testDispersion(simulateResiduals(model_full, n = 1000))
plot(simulateResiduals(model_full, n = 1000))

drop1(model_full, test = "Chisq")

emmeans(model_full, ~ sections, type = "response")  %>% pairs()



with(long_data  %>% filter(str_detect(method_increase, "WGS")), table(induced_binary, donor_shared, sections))

model_clean <- glmmTMB(induced_binary ~ sections + (1 | patient), 
family = binomial, data = long_data %>% 
  filter(!presence %in% c("Shared_after_None_pre", "Shared_after_VLP_pre", "Shared_pre_VLP_after", "VLP_pre_VLP_after", "VLP_pre_WGS_after")))

AIC(model_full, model_clean)
anova(model_full, model_clean)
drop1(model_clean, test = "Chisq")

emmeans(model_clean, ~ sections | donor_shared, type = "response")  %>% pairs()

# induced_phages  <- long_data %>% 
#     mutate(induced = if_else(induced_binary == 1, "Induced", "Not_induced")) %>% 
#   filter(!presence %in% c("Shared_after_None_pre", "Shared_after_VLP_pre", "Shared_pre_VLP_after", "VLP_pre_VLP_after", "VLP_pre_WGS_after"))

# save(induced_phages, file = "~/PIMMSgit/data/Rdata/induced_phages.Rdata")

##################################################################################


#Select phages that were induced and than reverted back########################

UnD_first_shared  <- 
    method_shared_first_week %>% 
    semi_join(method_shared_later, by = c("genome", "patient")) %>% 
    distinct(genome, patient, sections, method_increase, presence, donor_shared, b_phylum, host_range, Pred, class) %>%
    mutate(induced_binary = if_else(
    presence %in% c("Shared_after_WGS_pre"), "Induced", "Not_induced"
         ))  %>% 
    select(genome, patient, induced_binary, method_increase, presence, b_phylum)


UnD_later_shared  <- 
    method_shared_later %>% 
    semi_join(method_shared_first_week, by = c("genome", "patient")) %>% 
    distinct(genome, patient, sections, method_increase, presence, donor_shared, b_phylum, host_range, Pred, class) %>%
    mutate(reduced_binary = if_else(
    presence %in% c("Shared_pre_WGS_after"), "Reduced", "Not_reduced"
         ))  %>% 
    select(genome, patient, reduced_binary, method_increase, presence, b_phylum)

UnD  <- left_join(UnD_first, UnD_later, by = c("genome", "patient")) %>%  
        mutate(UnD = case_when(
              induced_binary == "Induced" & reduced_binary == "Reduced" ~ "UnD",
              induced_binary == "Induced" & reduced_binary == "Not_reduced" ~ "Not_UnD", 
        )) %>% 
        filter(!is.na(UnD))


model_data  <- UnD %>% 
    filter(str_detect(method_increase.x, "WGS"), str_detect(method_increase.y, "WGS"))  %>% 
    mutate(UnD_binary = if_else(UnD == "UnD", 1L, 0L)) %>% 
    distinct(genome, patient, b_phylum.x, b_phylum.y, UnD_binary)


with(UnD, table(UnD, method_increase.y, presence.y)) %>% view()
