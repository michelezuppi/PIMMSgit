#Determine_snvs_shared_cutoff_value


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
    group_by(scaffold, individual, timepoint) %>% #The total number of possible position found in that scaffold
    summarise(total_positions = n_distinct(position),    # true total before any filtering
              .groups = "drop")

head(total_positions)

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

average <- combined_proportions %>%
    filter(prop_combined > 0) %>%
    group_by(scaffold, same_individual) %>%
    summarise(mean_prop_combined = mean(prop_combined), .groups = "drop")

view(average)
# ---- Extract log10 values per group ----
within_vals  <- average %>%
    filter(same_individual == "Within individual") %>%
    pull(mean_prop_combined) %>%
    log10()
quantile(within_vals, probs = c(0.05, 0.95))

# within_vals <- within_vals[within_vals >= -4.777852 & within_vals <= -2.412605]

between_vals <- average %>%
    filter(same_individual == "Between individuals") %>%
    pull(mean_prop_combined) %>%
    log10()
quantile(between_vals, probs = c(0.05, 0.95))
# between_vals <- between_vals[between_vals >= -4.777852 & between_vals <= -2.412605]

#Sanity check that we are only keeping one observation per group
average  %>% 
    group_by(scaffold, same_individual) %>% 
    filter(duplicated(scaffold))#Done


#Plot the distribution differences
ks_results  <- ks.test(within_vals, between_vals)

average  %>% 
    ggplot(aes(x = log10(mean_prop_combined), fill = same_individual))+
    geom_density(alpha = 0.5, color = "black")+
    theme_classic2()+
    scale_fill_manual(values = c("Within individual"   = "#0072b2",
                                  "Between individuals" = "#d55e00"))+
    labs(fill = NULL, x = "Average SNVs proportion adjusted by contig length")+
    annotate("text", 
            x = -2,
            y = 0.6,
            label = paste0("KS test = ", round(ks_results$p.value, 12)
                ))

ggsave(filename = "~/PIMMSgit/plots/SNVs_analysis/shared_snvs_frequency/Average_SNVs_proportion.pdf", plot = last_plot())
# Estimate densities on common grid
x_grid <- seq(
    min(c(within_vals, between_vals)),
    max(c(within_vals, between_vals)),
    length.out = 1000
)

dens_within  <- approxfun(density(within_vals))(x_grid)
dens_between <- approxfun(density(between_vals))(x_grid)

# Compute difference
dens_diff <- tibble(
    x    = x_grid,
    diff = dens_within - dens_between  # positive = within dominates, negative = between dominates
) %>%
    filter(!is.na(diff))

# Point of maximum divergence
max_diff <- dens_diff %>% slice_max(abs(diff), n = 1)
cat("Max divergence (log10):", round(max_diff$x, 4), "\n")
cat("Max divergence (raw):",   round(10^max_diff$x, 6), "\n")

max_diff
# Plot
dens_diff %>%
    ggplot(aes(x = x, y = diff)) +
    geom_line(colour = "#0072b2", linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_ribbon(aes(ymin = 0, ymax = diff),
                fill = "#0072b2", alpha = 0.2) +
     geom_vline(xintercept = max_diff$x,
               linetype = "dashed", colour = "firebrick") +
    annotate("text",
             x     = max_diff$x + 0.05,
             y     = Inf,
             label = paste0("log10 = ", round(max_diff$x, 3),
                            "\nraw = ",   round(10^max_diff$x, 5)),
             hjust = -0.1, vjust = 2,
             colour = "red4", size = 5) +
    labs(
        x       = "Average SNVs proportion adjusted by contig length",
        y       = "Density difference (within - between)",
        caption = "Positive = within individual dominates | Negative = between individuals dominates"
    ) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())

ggsave(filename = "~/PIMMSgit/plots/SNVs_analysis/shared_snvs_frequency/Highest_distribution_difference.pdf", plot = last_plot())
