# # ---- Libraries ----
# library(tidyverse)
# library(conflicted)
# library(data.table)
# library(patchwork)
# conflict_prefer_all("dplyr", quiet = TRUE)

# rm(list = ls())
# load("~/PIMMSgit/data/Rdata/snvs.Rdata")
# load("~/PIMMSgit/data/Rdata/PIMMs_vOTUs.Rdata")

scaffold_length <- PIMMs_vOTUs %>%
    select(representative, contig_length) %>%
    group_by(representative) %>%
    slice(1) %>%                                         # deduplicate upfront once
    ungroup()

donor_scaffolds <- snvs %>%
    filter(timepoint == "D004",
           scaffold %in% scaffold_length$representative) %>%
    distinct(scaffold)

pre_fmt_donor_snvs_initial <-
    snvs %>%
    # Filter to scaffolds with known length FIRST — reduces data size immediately
    filter(scaffold %in% donor_scaffolds$scaffold,       # only donor scaffolds
           timepoint %in% c(0, -7, "D004")) %>%         # only pre-FMT + donor rows
    # Keep only scaffold + individual combinations that include the donor timepoint
    group_by(scaffold, individual) %>%
    filter(any(timepoint == "D004")) %>%
    ungroup()

pre_fmt_donor_snvs  <- 
    pre_fmt_donor_snvs_initial %>% 
    # Filter to positions seen at 2+ timepoints
    group_by(scaffold, position) %>%
    filter(n_distinct(timepoint) > 1) %>%
    ungroup() %>%
    # Join contig length
    left_join(scaffold_length,
              by = c("scaffold" = "representative"))

post_fmt_donor_snvs_initial <- snvs %>%
    #Select the scaffolds found in donors
    filter(scaffold %in% donor_scaffolds$scaffold,       # only donor scaffolds
           !timepoint %in% c(0, -7)) %>%         # only pre-FMT + donor rows
    # Keep only scaffold+individual combinations that include the donor timepoint
    group_by(scaffold, individual) %>%
    filter(any(timepoint == "D004")) %>%
    ungroup()

post_fmt_donor_snvs <- 
    post_fmt_donor_snvs_initial %>% 
    # Filter to positions seen at 2+ timepoints
    group_by(scaffold, position) %>%
    filter(n_distinct(timepoint) > 1) %>%
    ungroup() %>%
    # Join contig length
    left_join(scaffold_length,
              by = c("scaffold" = "representative"))
#Sanity checks

#Check that all the scaffolds were found also in donors
pre_fmt_donor_snvs_initial  %>% 
    group_by(scaffold, individual) %>% 
    filter(!any(timepoint == "D004")) #0

pre_fmt_donor_snvs_initial  %>% 
    group_by(scaffold, individual) %>% 
    filter(!any(timepoint == "D004")) #0

#Check that all the positions are duplicated
pre_fmt_donor_snvs  %>% 
    group_by(scaffold, position) %>% 
    filter(n_distinct(timepoint) < 2) #0

post_fmt_donor_snvs  %>% 
    group_by(scaffold, position) %>% 
    filter(n_distinct(timepoint) < 2) #0
#Confirm the timepoints

table(pre_fmt_donor_snvs$timepoint)
table(post_fmt_donor_snvs$timepoint)
#Is everything 0? Good
#Good, now we have repeated positions within the same scaffold, between donor scaffolds found in an individual either pre or post transplant.

#Now we calculate the average amount of shared SNVs between scaffolds

#We compare all the contigs that have been found between the two groups, how many of these have a proportion of shared SNVs that is comparable
#To contigs that are surely the same


pre_fmt_snvs_proportion  <- 
    pre_fmt_donor_snvs %>% 
    group_by(scaffold, contig_length) %>% 
    summarise(n_shared = n_distinct(position), .groups = "drop") %>% 
    mutate(prop_combined = n_shared/contig_length) %>% 
    filter(prop_combined > 0) %>% 
    mutate(FMT = "pre_FMT")

post_fmt_donor_proportion <- 
    post_fmt_donor_snvs  %>% 
    group_by(scaffold, contig_length) %>% 
    summarise(n_shared = n_distinct(position), .groups = "drop") %>% 
    mutate(prop_combined = n_shared/contig_length) %>% 
    filter(prop_combined > 0) %>% 
    mutate(FMT = "post_FMT")

combined_prop <- bind_rows(post_fmt_donor_proportion, pre_fmt_snvs_proportion)

view(combined_prop)

combined_prop %>%
    mutate(FMT = factor(FMT, levels = c("pre_FMT", "post_FMT"))) %>%
    ggplot(aes(x = FMT, y = log10(prop_combined))) +
    geom_boxplot() +
    stat_compare_means(method = "wilcox", label.x = 1.3) +
    geom_quasirandom(shape = 21) +
    geom_hline(yintercept = -3.299, linetype = "dashed", colour = "red4") +
    annotate("text",
             x     = 1.6,
             y     = -4.1,
             label = "Same phage SNVs cutoff: \nlog10 = -3.299 \n(raw = 5e-04)",
             hjust = 1,
             vjust = 1,
             size  = 4,
             colour = "red4")+
    theme_classic2()+
    labs(x = NULL,
        y = "Average SNVs proportion adjusted by contig length",
        title = "Donor-shared phages")

ggsave(filename = "~/PIMMSgit/plots/SNVs_analysis/shared_snvs_frequency/Putative_engrafted_contigs_SNVs_proportion.pdf", plot = last_plot())

chisq_results <- combined_prop  %>% 
    mutate(cutoff = if_else(prop_combined >= 0.0005, "Yes", "No")) %>% 
    group_by(FMT, cutoff) %>% 
    count() %>% 
    ungroup() %>%
    pivot_wider(names_from = cutoff, values_from = n) %>% 
    select(-1) %>% 
    chisq.test()

p_label <- paste0("χ²   p = ", round(chisq_results$p.value, 5))

combined_prop  %>% 
    mutate(FMT = factor(FMT, 
                                    levels = c("pre_FMT", "post_FMT"))) %>%  
    mutate(cutoff = if_else(prop_combined >= 0.0005, "Yes", "No")) %>% 
    group_by(FMT, cutoff) %>% 
    count() %>% 
    ungroup()  %>% 
    group_by(FMT) %>% 
    mutate(tot = sum(n)) %>% 
    ungroup() %>% 
    mutate(prop = n / tot * 100) %>% 
    ggplot(aes(x = FMT, y = prop, fill = cutoff))+
    geom_col(color = "black") +
    scale_fill_viridis_d(option = "magma", begin = 0.3, end = 0.8)+
    theme_classic2()+
        annotate("text",
             x     = 1.5,                               # centred between the two bars
             y     = 105,                               # just above the bars
             label = p_label,
             size  = 5,
             fontface = "italic")+
    labs(x = NULL,
         y = "Percentage over total donor shared phages\n(%)",
        fill = "Same phage")

ggsave(filename = "~/PIMMSgit/plots/SNVs_analysis/shared_snvs_frequency/Percentage_of_same_donor_phages_before_and_after_transplant.pdf", plot = last_plot())

