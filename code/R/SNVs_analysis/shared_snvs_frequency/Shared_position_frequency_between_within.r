
# ---- Libraries ----
library(tidyverse)
library(conflicted)
library(data.table)  # for fast operations
conflict_prefer_all("dplyr", quiet = TRUE)

rm(list = ls())
load("~/PIMMSgit/data/Rdata/snvs.Rdata")
load("~/PIMMSgit/data/Rdata/PIMMs_vOTUs.Rdata")

scaffold_length <-
     PIMMs_vOTUs  %>% 
     select(representative, contig_length)

# ---- Compute total positions per scaffold per individual per timepoint BEFORE any filtering ----
total_positions <- snvs %>%
    filter(timepoint != -7) %>%                          # only remove -7, nothing else
    mutate(timepoint = case_when(
        timepoint == 0      ~ "Baseline",
        timepoint == 2      ~ "Day_2",
        timepoint == 4      ~ "Day_4",
        timepoint == 7      ~ "Day_7",
        timepoint == 14     ~ "Day_14",
        timepoint == 28     ~ "Day_28",
        timepoint == "D004" ~ "Donor",
        TRUE ~ as.character(timepoint)
    )) %>%
    group_by(scaffold) %>% #The total number of possible position found in that scaffold
    summarise(total_positions = n_distinct(position),    # true total before any filtering
              .groups = "drop")


# ---- Pre-filter for the main analysis ----
snvs_prefiltered <- snvs %>%
    filter(timepoint != -7) %>%
    group_by(scaffold) %>%
    filter(allele_count > 1) %>%
    filter(duplicated(position) | duplicated(position, fromLast = TRUE)) %>%
    ungroup() %>%
    group_by(position, scaffold) %>%
    filter(n_distinct(timepoint) > 1) %>% #Find the SNVs position that were present at more than 1 timepoint
    ungroup() %>%
    select(scaffold, position, timepoint, individual, var_base) %>%
    left_join(
        scaffold_length %>% distinct(representative, contig_length),
        by = c("scaffold" = "representative")
    ) %>%
    mutate(timepoint = case_when(
        timepoint == 0      ~ "Baseline",
        timepoint == 2      ~ "Day_2",
        timepoint == 4      ~ "Day_4",
        timepoint == 7      ~ "Day_7",
        timepoint == 14     ~ "Day_14",
        timepoint == 28     ~ "Day_28",
        timepoint == "D004" ~ "Donor",
        TRUE ~ as.character(timepoint)
    )) %>%
    # Join true total positions from unfiltered data
    left_join(total_positions,
              by = c("scaffold", "individual", "timepoint")) %>% 
    filter(!is.na(contig_length)) #Remove the ones without known length


# ---- Helper function: run analysis for a given anchor ----
run_analysis <- function(data, anchor) {
    anchor_data <- data %>%
        filter(timepoint == anchor) %>%
        select(scaffold, position,
               var_base_anchor    = var_base,
               individual_anchor  = individual,
               total_pos_anchor   = total_positions) %>%  # carry anchor total
        distinct() 
    #This first step selects the anchor timepoint. The anchor is used in the later analysis as the timepoint to compare the other ones to.

    data %>%
        filter(timepoint != anchor) %>%
        inner_join(anchor_data,
                   by = c("scaffold", "position"),
                   relationship = "many-to-many") %>%
        mutate(
            same_variant    = if_else(var_base == var_base_anchor, "Yes", "No"),
            same_individual = if_else(individual == individual_anchor, "Yes", "No"),
            timepoint_pair  = paste0(anchor, "_vs_", timepoint),
            min_total       = pmin(total_positions,        # minimum of the two timepoint totals
                                   total_pos_anchor)       # computed per row
        ) %>%
        distinct(scaffold, position, same_individual, same_variant,
                 timepoint_pair, contig_length, min_total)
}

# ---- Run analysis ----
anchor <- "Baseline"
results_0 <- run_analysis(snvs_prefiltered, anchor = anchor)

# ---- Summarise combined proportion per scaffold per group ----
combined_proportions <- results_0 %>%
    group_by(scaffold, same_individual, timepoint_pair, contig_length) %>%
    summarise(
        n_shared  = n_distinct(position),
        min_total = min(min_total),
        .groups   = "drop"
    ) %>%
    mutate(
        prop_combined   = n_shared / contig_length,  # now all values are scalar
        same_individual = if_else(same_individual == "Yes",
                                  "Within individual", "Between individuals")
    ) %>% 
    filter(!is.na(timepoint_pair))

# ---- Statistical test ----
fisher_combined <- combined_proportions %>%
    group_by(timepoint_pair, same_individual) %>%
    summarise(mean_prop = mean(prop_combined), .groups = "drop") %>%
    pivot_wider(names_from = same_individual, values_from = mean_prop)

# KS test per timepoint pair — comparing distributions of prop_combined
ks_results <- combined_proportions %>%
    group_by(timepoint_pair) %>%
    group_split() %>%
    map_dfr(function(df) {
        tp      <- unique(df$timepoint_pair)
        within  <- df %>% filter(same_individual == "Within individual")  %>% pull(prop_combined)
        between <- df %>% filter(same_individual == "Between individuals") %>% pull(prop_combined)

        test <- ks.test(within, between)

        tibble(
            timepoint_pair = tp,
            p_value        = test$p.value,
            d_statistic    = test$statistic
        )
    }) %>%
    mutate(
        p_adj     = p.adjust(p_value, method = "BH"),
        sig_label = case_when(
            p_adj < 0.001 ~ "***",
            p_adj < 0.01  ~ "**",
            p_adj < 0.05  ~ "*",
            TRUE          ~ "ns"
        )
    )

# ---- Plot ----
tp_order <- c("Baseline_vs_Day_2", "Baseline_vs_Day_4", "Baseline_vs_Day_7",
              "Baseline_vs_Day_14", "Baseline_vs_Day_28")

p <- combined_proportions %>%
    mutate(timepoint_pair = factor(timepoint_pair, levels = tp_order)) %>%
    left_join(ks_results %>% select(timepoint_pair, sig_label),
              by = "timepoint_pair") %>%              
    filter(!is.na(timepoint_pair)) %>% 
    ggplot(aes(x = same_individual, y = log10(prop_combined), fill = same_individual)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1, shape = 21) +
    # geom_text(
    #     data = ks_results %>%
    #         mutate(timepoint_pair = factor(timepoint_pair, levels = tp_order)),
    #     aes(x = 1.5, y = Inf, label = sig_label),
    #     inherit.aes = FALSE,
    #     vjust = 2, size = 5
    # ) +
    scale_fill_manual(values = c("Within individual"   = "#0072b2",
                                  "Between individuals" = "#d55e00")) +
    scale_y_continuous(labels = scales::scientific) +
    ggpubr::stat_compare_means(method = "wilcox", 
    label = "p.signif", size = 6)+
    facet_wrap(~ timepoint_pair, nrow = 1) +
    labs(
        x       = NULL,
        y       = "Shared positions / contig length",
        fill    = NULL,
        title   = "Combined SNV sharing proportion — within vs between individuals",
        caption = "KS test; * p.adj < 0.05, ** p.adj < 0.01, *** p.adj < 0.001"
    ) +
    theme_bw() +
    theme(
        legend.position  = "top",
        panel.grid.minor = element_blank(),
        axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank()
    )

p

ggsave(filename = "Shared_SNVs_position_freq_between_within.pdf", plot = p, height = 8, width = 6)
