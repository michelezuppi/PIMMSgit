#Check if integrating multiple timepoints increases the proportion of shared SNVs

# ---- Libraries ----
library(tidyverse)
library(conflicted)
library(data.table)
library(patchwork)
conflict_prefer_all("dplyr", quiet = TRUE)

rm(list = ls())
load("~/PIMMSgit/data/Rdata/snvs.Rdata")
load("~/PIMMSgit/data/Rdata/PIMMs_vOTUs.Rdata")
load("~/PIMMSgit/data/Rdata/breadth_cutoff.Rdata")


scaffold_length <- PIMMs_vOTUs %>%
    select(representative, contig_length) %>%
    group_by(representative) %>%
    slice(1) %>%                                         # deduplicate upfront once
    ungroup()

snvs_breadth_cutoff <- 
    snvs %>% 
    filter(scaffold %in% breadth_cutoff$genome)

donor_scaffolds <- snvs_breadth_cutoff %>%
    filter(timepoint == "D004",
           scaffold %in% scaffold_length$representative) %>%
    distinct(scaffold)

# post_existing <- snvs_breadth_cutoff %>%
#     filter(scaffold %in% donor_scaffolds$scaffold,
#            timepoint %in% c(2, 4, 7, 14, 28)) %>%
#     distinct(scaffold, individual)

pre_fmt_donor_snvs_initial <-
    snvs_breadth_cutoff %>%
    # Filter to scaffolds with known length FIRST — reduces data size immediately
    filter(scaffold %in% donor_scaffolds$scaffold,       # only donor scaffolds
           timepoint %in% c(0, -7, "D004")) %>%         # only pre-FMT + donor rows
    # Keep only scaffold + individual combinations that include the donor timepoint
    # anti_join(post_existing,                              # same exclusion as post-FMT
    #           by = c("scaffold", "individual")) %>%
    group_by(scaffold, individual) %>%
    filter(any(timepoint == "D004")) %>%
    ungroup()

pre_fmt_donor_snvs  <- 
    pre_fmt_donor_snvs_initial %>% 
    # Filter to positions seen at 2+ timepoints
    group_by(scaffold, position, individual) %>%
    filter(n_distinct(timepoint) > 1 & any(timepoint == "D004")) %>%
    ungroup() %>%
    # Join contig length
    left_join(scaffold_length,
              by = c("scaffold" = "representative"))
# pre_existing <- snvs_breadth_cutoff %>%
#     filter(scaffold %in% donor_scaffolds$scaffold,
#            timepoint %in% c(-7, 0)) %>%
#     distinct(scaffold, individual)

post_fmt_donor_snvs_initial <- snvs_breadth_cutoff %>%
    filter(scaffold %in% donor_scaffolds$scaffold,
           !timepoint %in% c(0, -7)) %>%
    # anti_join(pre_existing,                              # same exclusion as pre-FMT
    #           by = c("scaffold", "individual")) %>%
    group_by(scaffold, individual) %>%
    filter(any(timepoint == "D004")) %>%
    ungroup()

# post_fmt_donor_snvs <- 
#     post_fmt_donor_snvs_initial %>% 
#     # Filter to positions seen at 2+ timepoints
#     group_by(scaffold, position) %>%
#     filter(n_distinct(timepoint) > 1 & any(timepoint == "D004")) %>% 
#     #This is so we don't get a biased amount of observation per contig. Now we have the same number of observation as the pre-FMT ones, 2 timepoints + donor
#     ungroup() %>%
#     # Join contig length
#     left_join(scaffold_length,
#               by = c("scaffold" = "representative"))

post_fmt_donor_snvs <- post_fmt_donor_snvs_initial %>%
    filter(timepoint != "D004") %>% 
    group_by(scaffold, position, individual) %>%
    filter(n_distinct(timepoint) > 1) %>%
    ungroup() %>%
    left_join(scaffold_length,
              by = c("scaffold" = "representative"))
              
post_fmt_donor_proportion <- 
    post_fmt_donor_snvs  %>% 
    group_by(scaffold, contig_length, individual) %>% 
    mutate(n_shared = n_distinct(position))  %>% ungroup() %>% 
    group_by(scaffold, position, individual) %>% 
    mutate(n_timepoints = n_distinct(timepoint)) %>% ungroup() %>% 
    mutate(prop_combined = n_shared/contig_length)

post_fmt_donor_proportion %>% 
    distinct(scaffold, n_timepoints, individual, prop_combined) %>% 
    ggplot(aes(x = as.factor(n_timepoints), y = log10(prop_combined)))+
    geom_boxplot()+
    geom_quasirandom(shape = 21)+
    stat_compare_means(method = "kruskal")+
    theme_classic2()+
    labs(x = "Number of timepoints the scaffold was identified at")


ggsave(file = "~/PIMMSgit/plots/correlation_between_number_of_timepoints_and_shared_SNVs_prop.pdf", plot = last_plot())
