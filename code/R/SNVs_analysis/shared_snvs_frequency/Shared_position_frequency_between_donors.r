
# # ---- Libraries ----
# library(tidyverse)
# library(conflicted)
# library(data.table)  # for fast operations
# conflict_prefer_all("dplyr", quiet = TRUE)

# rm(list = ls())
# load("~/PIMMSgit/data/Rdata/snvs.Rdata")
# load("~/PIMMSgit/data/Rdata/PIMMs_vOTUs.Rdata")

# scaffold_length <-
#      PIMMs_vOTUs  %>% 
#      select(representative, contig_length)

# ---- Compute total positions per scaffold per individual per timepoint BEFORE any filtering ----
total_positions <- snvs %>%
    filter(timepoint != -7) %>%
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
    group_by(scaffold, individual, timepoint) %>%
    summarise(total_positions = n_distinct(position),
              .groups = "drop")

# ---- Pre-filter for the main analysis ----
snvs_prefiltered <- snvs %>%
    filter(timepoint != -7) %>%
    group_by(scaffold) %>%
    filter(allele_count > 1) %>%
    filter(duplicated(position) | duplicated(position, fromLast = TRUE)) %>%
    ungroup() %>%
    group_by(position, scaffold) %>%
    filter(n_distinct(timepoint) > 1) %>%
    ungroup() %>%
    select(scaffold, position, timepoint, individual, var_base) %>%
    left_join(
        scaffold_length %>%
            group_by(representative) %>%
            slice(1) %>%
            ungroup() %>%
            select(representative, contig_length),
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
    left_join(total_positions,
              by = c("scaffold", "individual", "timepoint")) %>%
    filter(!is.na(contig_length))

# ---- Helper function: run analysis for a given anchor ----
run_analysis <- function(data, anchor) {
    anchor_data <- data %>%
        filter(timepoint == anchor) %>%
        select(scaffold, position,
               var_base_anchor   = var_base,
               individual_anchor = individual,
               total_pos_anchor  = total_positions) %>%
        distinct()

    data %>%
        filter(timepoint != anchor) %>%
        inner_join(anchor_data,
                   by = c("scaffold", "position"),
                   relationship = "many-to-many") %>%
        mutate(
            same_variant    = if_else(var_base == var_base_anchor, "Yes", "No"),
            same_individual = if_else(individual == individual_anchor, "Yes", "No"),
            timepoint_pair  = paste0(anchor, "_vs_", timepoint),
            min_total       = pmin(total_positions, total_pos_anchor)
        ) %>%
        distinct(scaffold, position, same_individual, same_variant,
                 timepoint_pair, contig_length, min_total)
}

# ---- Define anchor ----
anchor <- "Baseline"

# ---- Between-individual analysis ----
results_between <- run_analysis(snvs_prefiltered, anchor = anchor) %>%
    filter(same_individual == "No") %>%                  # keep only between-individual comparisons
    mutate(group = "Between individuals")

# ---- Donor vs recipient analysis ----
donor_data <- snvs_prefiltered %>%
    filter(timepoint == "Donor") %>%
    select(scaffold, position,
           var_base_donor  = var_base,
           total_pos_donor = total_positions) %>%        # carry donor total positions
    distinct()

results_donor <- snvs_prefiltered %>%
    filter(!timepoint %in% c("Donor", anchor)) %>%
    inner_join(donor_data,
               by = c("scaffold", "position"),
               relationship = "many-to-many") %>%
    mutate(
        same_variant   = if_else(var_base == var_base_donor, "Yes", "No"),
        timepoint_pair = paste0(anchor, "_vs_", timepoint),
        min_total      = pmin(total_positions, total_pos_donor)  # min of donor and recipient totals
    ) %>%
    distinct(scaffold, position, same_variant,
             timepoint_pair, contig_length, min_total) %>%
    mutate(group = "Donor")

# ---- Combine and summarise ----
tp_order <- c("Baseline_vs_Day_2", "Baseline_vs_Day_4", "Baseline_vs_Day_7",
              "Baseline_vs_Day_14", "Baseline_vs_Day_28")

combined_proportions <- bind_rows(
        results_between %>% select(scaffold, position, same_variant,
                                   timepoint_pair, contig_length, min_total, group),
        results_donor   %>% select(scaffold, position, same_variant,
                                   timepoint_pair, contig_length, min_total, group)
    ) %>%
    group_by(scaffold, group, timepoint_pair, contig_length) %>%
    summarise(
        n_shared  = n_distinct(position),
        min_total = min(min_total),
        .groups   = "drop"
    ) %>%
    mutate(
        prop_combined  = n_shared / contig_length,
        timepoint_pair = factor(timepoint_pair, levels = tp_order)
    ) %>%
    filter(!is.na(timepoint_pair))

# ---- KS test per timepoint pair ----
ks_results <- combined_proportions %>%
    group_by(timepoint_pair) %>%
    group_split() %>%
    map_dfr(function(df) {
        tp      <- unique(df$timepoint_pair)
        between <- df %>% filter(group == "Between individuals") %>% pull(prop_combined)
        donor   <- df %>% filter(group == "Donor")               %>% pull(prop_combined)

        test <- ks.test(between, donor)

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
p <- combined_proportions %>%
    left_join(ks_results %>%
                  mutate(timepoint_pair = factor(timepoint_pair, levels = tp_order)) %>%
                  select(timepoint_pair, sig_label),
              by = "timepoint_pair") %>%
    ggplot(aes(x = group, y = log10(prop_combined), fill = group)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    geom_text(
        data = ks_results %>%
            mutate(timepoint_pair = factor(timepoint_pair, levels = tp_order)),
        aes(x = 1.5, y = Inf, label = sig_label),
        inherit.aes = FALSE,
        vjust = 2, size = 5
    ) +
    # ggpubr::stat_compare_means(method = "wilcox", 
    # label = "p.signif", size = 6, 
    # label.x = 1.5)+
    scale_fill_manual(values = c(
        "Between individuals" = "#d55e00",
        "Donor"               = "#009e73"
    )) +
    facet_wrap(~ timepoint_pair, nrow = 1) +
    labs(
        x       = NULL,
        y       = "log10(Shared positions / contig length)",
        fill    = NULL,
        title   = "SNV sharing proportion — between individuals vs donor",
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

ggsave(filename = "Shared_SNVs_position_freq_between_donor.pdf", plot = p, height = 8, width = 6)