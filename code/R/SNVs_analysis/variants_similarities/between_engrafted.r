
# # ---- Libraries ----
# library(tidyverse)
# library(conflicted)
# library(data.table)  # for fast operations
# conflict_prefer_all("dplyr", quiet = TRUE)
# rm(list = ls())

# load("~/PIMMSgit/data/Rdata/snvs.Rdata")
# ---- Pre-filter once before any analysis ----
snvs_prefiltered <- snvs %>%
    filter(timepoint != -7) %>%                          # Remove timepoint -7 from all analyses (numeric match)
    group_by(scaffold) %>%
    filter(allele_count > 1) %>%                         # Keep only polymorphic positions (true SNVs, not fixed sites)
    filter(duplicated(position) | duplicated(position,   # Keep positions appearing in 2+ samples within
           fromLast = TRUE)) %>%                         # the same scaffold (shared across samples)
    ungroup() %>%
    group_by(position, scaffold) %>%
    filter(n_distinct(timepoint) > 1) %>%                # Keep only positions seen at 2+ timepoints (persistent SNVs)
    ungroup() %>%
    select(scaffold, position, timepoint,                # Drop unused columns to reduce memory overhead
           individual, var_base) %>%
    mutate(timepoint = case_when(                        # Relabel timepoints AFTER filtering so -7 case is never needed
        timepoint == "0"    ~ "Baseline",
        timepoint == "2"    ~ "Day_2",
        timepoint == "4"    ~ "Day_4",
        timepoint == "7"    ~ "Day_7",
        timepoint == "14"   ~ "Day_14",
        timepoint == "28"   ~ "Day_28",
        timepoint == "D004" ~ "Donor",
        TRUE ~ timepoint                                 # Catch-all: keep anything unlisted as-is
    ))

# ---- Helper function: run analysis for a given anchor timepoint ----
run_analysis <- function(data, anchor) {
    anchor_data <- data %>%
        filter(timepoint == anchor) %>%                  # Extract all rows at the anchor timepoint
        select(scaffold, position,
               var_base_anchor   = var_base,             # Rename to avoid column conflicts after join
               individual_anchor = individual) %>%
        distinct()                                       # One row per unique scaffold+position+individual

    data %>%
        filter(timepoint != anchor) %>%                  # Left side of join: all non-anchor timepoint rows
        inner_join(anchor_data,                          # Join anchor rows onto non-anchor rows by scaffold+position;
                   by = c("scaffold", "position"),       # inner_join keeps ONLY positions present at BOTH the anchor
                   relationship = "many-to-many") %>%    # and the other timepoint (persistence filter)
        mutate(
            same_variant    = if_else(                   # Compare variant base at non-anchor timepoint
                var_base == var_base_anchor,             # to variant base at anchor timepoint;
                "Yes", "No"),                            # "Yes" = same variant base maintained over time
            same_individual = if_else(                   # Compare individual ID at non-anchor timepoint
                individual == individual_anchor,         # to individual ID at anchor timepoint;
                "Yes", "No"),                            # "Yes" = within-individual, "No" = between-individual
            timepoint_pair  = paste0(anchor, "_vs_",    # Create label for the timepoint pair using
                                     timepoint)          # the anchor variable — not hardcoded
        ) %>%
        distinct(scaffold, position, same_individual,    # Collapse to one row per unique combination
                 same_variant, timepoint_pair)
}

# ---- Define anchor ----
anchor <- "Baseline"                                     # Define once — used by both run_analysis and results_donor

# ---- Between-individual analysis ----
results_between <- run_analysis(snvs_prefiltered,
                                anchor = anchor) %>%
    filter(same_individual == "No") %>%                  # Keep only positions shared between DIFFERENT individuals
    mutate(group = "Between individuals")                # Label for plot

# ---- Donor vs recipient analysis ----
donor_data <- snvs_prefiltered %>%
    filter(timepoint == "Donor") %>%                     # Extract donor-only rows
    select(scaffold, position,
           var_base_donor = var_base) %>%                # Rename to avoid conflicts after join
    distinct()                                           # One row per unique scaffold+position

results_donor <- snvs_prefiltered %>%
    filter(!timepoint %in% c("Donor", anchor)) %>%       # Keep only recipient timepoints —
                                                         # excludes donor and anchor (-7 already removed upstream)
    inner_join(donor_data,                               # Join donor positions onto recipient rows;
               by = c("scaffold", "position"),           # keeps only scaffold+position combinations present
               relationship = "many-to-many") %>%        # in BOTH donor and recipient timepoints
    mutate(
        same_variant   = if_else(                        # Compare recipient variant base to donor variant base;
            var_base == var_base_donor,                  # "Yes" = recipient carries same variant as donor
            "Yes", "No"),                                # (potential signature of engraftment)
        timepoint_pair = paste0(anchor, "_vs_", timepoint) # Use anchor variable — not hardcoded string
    ) %>%
    distinct(scaffold, position,
             same_variant, timepoint_pair) %>%
    mutate(group = "Donor")                              # Label for plot

# ---- Sanity check: confirm timepoint pair labels match before proceeding ----
stopifnot(
    all(results_between %>% distinct(timepoint_pair) %>% pull() %in%
        c("Baseline_vs_Day_2", "Baseline_vs_Day_4", "Baseline_vs_Day_7",
          "Baseline_vs_Day_14", "Baseline_vs_Day_28")),
    all(results_donor %>% distinct(timepoint_pair) %>% pull() %in%
        c("Baseline_vs_Day_2", "Baseline_vs_Day_4", "Baseline_vs_Day_7",
          "Baseline_vs_Day_14", "Baseline_vs_Day_28"))
)

# ---- Statistical test: between-individual vs donor sharing rates ----
fisher_combined <- bind_rows(
        results_between %>%
            select(scaffold, position,
                   same_variant, timepoint_pair, group),
        results_donor %>%
            select(scaffold, position,
                   same_variant, timepoint_pair, group)
    ) %>%
    group_by(timepoint_pair, group, same_variant) %>%
    count() %>%                                          # Count positions in each group×same_variant combination
    ungroup() %>%
    pivot_wider(names_from = group,                      # Reshape to wide: groups become columns,
                values_from = n,                         # giving a 2×2 contingency table per timepoint pair
                values_fill = 0) %>%                     # Fill missing combinations with 0 not NA
    group_by(timepoint_pair) %>%
    group_split() %>%                                    # Split into list: one dataframe per timepoint pair
    map_dfr(function(df) {
        tp  <- unique(df$timepoint_pair)                 # Capture timepoint pair label for output
        mat <- df %>%
            select(-timepoint_pair) %>%
            column_to_rownames("same_variant") %>%       # Convert same_variant to row names (required by test)
            as.matrix()                                  # Convert to numeric matrix for statistical test

        test <- if (any(mat < 5))                        # Adaptive test: Fisher's exact for small counts,
                    fisher.test(mat)                     # chi-square for large — avoids invalid approximations
                else chisq.test(mat)

        tibble(
            timepoint_pair = tp,
            p_value        = test$p.value
        )
    }) %>%
    mutate(
        p_adj     = p.adjust(p_value, method = "BH"),   # FDR correction across all timepoint pairs
        sig_label = case_when(
            p_adj < 0.001 ~ "***",
            p_adj < 0.01  ~ "**",
            p_adj < 0.05  ~ "*",
            TRUE          ~ "ns"
        )
    )

# ---- Prepare dumbbell plot data ----
tp_order <- c("Baseline_vs_Day_2",  "Baseline_vs_Day_4",  "Baseline_vs_Day_7",
              "Baseline_vs_Day_14", "Baseline_vs_Day_28")  # y-axis order

dumbbell_data <- bind_rows(results_between, results_donor) %>%
    group_by(timepoint_pair, group, same_variant) %>%
    count() %>%                                          # Count positions per group×timepoint×variant
    group_by(timepoint_pair, group) %>%
    mutate(prop = n / sum(n)) %>%                        # Proportion sharing same variant base
    filter(same_variant == "Yes") %>%                    # One row per group per timepoint pair
    ungroup() %>%
    left_join(                                           # Add significance labels
        fisher_combined %>% select(timepoint_pair, sig_label),
        by = "timepoint_pair") %>%
    mutate(timepoint_pair = factor(timepoint_pair,       # Enforce y-axis order
                                   levels = tp_order)) %>% 
    filter(!is.na(timepoint_pair))

max_prop <- max(dumbbell_data$prop, na.rm = TRUE)        # For positioning significance labels

# ---- Dumbbell plot ----
p <- dumbbell_data %>%                                   # Assign to object for safe ggsave
    ggplot(aes(x = prop, y = timepoint_pair)) +
    geom_line(aes(group = timepoint_pair),               # Line connecting the two group points per timepoint pair
              colour = "grey70", linewidth = 1) +
    geom_point(aes(colour = group), size = 4) +          # One point per group per timepoint pair
    geom_text(
        data = dumbbell_data %>%
            distinct(timepoint_pair, sig_label),         # One significance label per timepoint pair
        aes(x = max_prop + 0.05,
            y = timepoint_pair,
            label = sig_label),
        inherit.aes = FALSE,
        hjust = 0, size = 5
    ) +
    scale_x_continuous(
        labels = scales::percent,                        # Format x-axis as percentages
        limits = c(NA, max_prop + 0.12)                  # Expand right margin to fit significance labels
    ) +
    scale_colour_manual(values = c(
        "Between individuals" = "#d55e00",               # Colourblind-friendly palette (Okabe-Ito)
        "Donor"               = "#009e73")) +
    labs(
        x       = "% positions sharing variant base",
        y       = "Timepoint (vs Baseline)",
        colour  = NULL,
        title   = "Between-individual vs donor variant base sharing",
        caption = "* p.adj < 0.05, ** p.adj < 0.01, *** p.adj < 0.001"
    ) +
    theme_bw() +
    theme(
        legend.position  = "top",
        panel.grid.minor = element_blank()
    )

p  # display plot

# ---- Save plot ----
ggsave(
    filename = "~/PIMMSgit/plots/Between_individuals_and_putative_engrafted.pdf",
    plot     = p,
    width    = 8,
    height   = 6
)