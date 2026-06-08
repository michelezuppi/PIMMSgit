# # ---- Libraries ----
# library(tidyverse)
# library(conflicted)
# library(data.table)
# conflict_prefer_all("dplyr", quiet = TRUE)
# rm(list = ls())

# load("~/PIMMSgit/data/Rdata/snvs.Rdata")

# ---- Pre-filter once before any analysis ----
snvs_prefiltered <- snvs %>%
    filter(timepoint != -7 &                             # Remove timepoint -7 (numeric match)
           timepoint != "D004") %>%                      # Remove donor timepoint upfront
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
    mutate(timepoint = case_when(                        # Relabel timepoints AFTER filtering
        timepoint == 0  ~ "Baseline",                    # Numeric comparisons — no quotes
        timepoint == 2  ~ "Day_2",
        timepoint == 4  ~ "Day_4",
        timepoint == 7  ~ "Day_7",
        timepoint == 14 ~ "Day_14",
        timepoint == 28 ~ "Day_28",
        TRUE ~ as.character(timepoint)                   # Safe catch-all
    ))

# ---- Sanity check: confirm no unexpected timepoint values ----
unexpected <- snvs_prefiltered %>%
    filter(!timepoint %in% c("Baseline", "Day_2", "Day_4",
                              "Day_7", "Day_14", "Day_28")) %>%
    distinct(timepoint)
if (nrow(unexpected) > 0) warning("Unexpected timepoint values: ",
                                   paste(unexpected$timepoint, collapse = ", "))

# ---- Helper function: run analysis for a given anchor ----
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
            timepoint_pair  = paste0(anchor, "_vs_",    # Create label using anchor variable — not hardcoded
                                     timepoint)
        ) %>%
        distinct(scaffold, position, same_individual,    # Collapse to one row per unique combination
                 same_variant, timepoint_pair)
}

# ---- Helper function: run statistical test ----
run_fisher <- function(data) {
    data %>%
        group_by(timepoint_pair, same_individual, same_variant) %>%
        count() %>%
        ungroup() %>%
        pivot_wider(names_from = same_individual,        # Reshape to wide: same_individual becomes columns
                    values_from = n,
                    values_fill = 0) %>%
        group_by(timepoint_pair) %>%
        group_split() %>%                                # Split into list: one dataframe per timepoint pair
        map_dfr(function(df) {
            tp  <- unique(df$timepoint_pair)
            mat <- df %>%
                select(-timepoint_pair) %>%
                column_to_rownames("same_variant") %>%   # Convert same_variant to row names (required by test)
                as.matrix()

            test <- if (any(mat < 5))                    # Adaptive test: Fisher's exact for small counts,
                        fisher.test(mat)                 # chi-square for large — avoids invalid approximations
                    else chisq.test(mat)

            tibble(
                timepoint_pair = tp,
                p_value        = test$p.value,
                n_within       = sum(df$Yes, na.rm = TRUE),  # Total within-individual positions
                n_between      = sum(df$No,  na.rm = TRUE)   # Total between-individual positions
            )
        }) %>%
        mutate(p_adj = p.adjust(p_value, method = "BH")) %>%  # FDR correction across all timepoint pairs
        arrange(p_adj)
}

# ---- Helper function: prepare dumbbell data ----
prep_dumbbell <- function(results, fisher_results, anchor_label) {
    results %>%
        group_by(timepoint_pair, same_individual, same_variant) %>%
        count() %>%
        group_by(timepoint_pair, same_individual) %>%
        mutate(prop = n / sum(n)) %>%                    # Proportion sharing same variant base
        filter(same_variant == "Yes") %>%                # One row per group per timepoint pair
        ungroup() %>%
        mutate(
            same_individual = if_else(                   # Convert Yes/No to readable labels
                same_individual == "Yes",
                "Within individual", "Between individuals"),
            anchor = anchor_label                        # Tag with anchor label for faceting
        ) %>%
        left_join(
            fisher_results %>%
                mutate(sig_label = case_when(            # Convert adjusted p-value to significance label
                    p_adj < 0.001 ~ "***",
                    p_adj < 0.01  ~ "**",
                    p_adj < 0.05  ~ "*",
                    TRUE          ~ "ns"
                )) %>%
                select(timepoint_pair, sig_label),
            by = "timepoint_pair"
        )
}

# ---- Define anchor ----
anchor <- "Baseline"                                     # Define once — used throughout

# ---- Run analysis ----
results_0  <- run_analysis(snvs_prefiltered, anchor = anchor)
fisher_0   <- run_fisher(results_0)
dumbbell_0 <- prep_dumbbell(results_0, fisher_0, anchor_label = "Timepoint 0")

# ---- Sanity check: confirm timepoint pair labels match tp_order ----
tp_order_0 <- c("Baseline_vs_Day_2", "Baseline_vs_Day_4", "Baseline_vs_Day_7",
                 "Baseline_vs_Day_14", "Baseline_vs_Day_28")

stopifnot(
    all(results_0 %>% distinct(timepoint_pair) %>% pull() %in% tp_order_0)
)

# ---- Apply factor order ----
dumbbell_0 <- dumbbell_0 %>%
    mutate(timepoint_pair = factor(timepoint_pair, levels = tp_order_0))

max_prop <- max(dumbbell_0$prop, na.rm = TRUE)           # For positioning significance labels

# ---- Dumbbell plot ----
p <- dumbbell_0 %>%                                      # Assign to object for safe ggsave
    ggplot(aes(x = prop, y = timepoint_pair)) +
    geom_line(aes(group = timepoint_pair),               # Line connecting two group points per timepoint pair
              colour = "grey70", linewidth = 1) +
    geom_point(aes(colour = same_individual), size = 4) + # One point per group per timepoint pair
    geom_text(
        data = dumbbell_0 %>%
            distinct(timepoint_pair, sig_label, anchor), # One significance label per timepoint pair
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
        "Within individual"   = "#0072b2",               # Colourblind-friendly palette (Okabe-Ito)
        "Between individuals" = "#d55e00")) +
    facet_wrap(~ anchor, ncol = 1, scales = "free_y") +  # Separate panels per anchor
    labs(
        x       = "% positions sharing variant base",
        y       = NULL,
        colour  = NULL,
        title   = "Within- vs between-individual variant base sharing",
        caption = "* p.adj < 0.05, ** p.adj < 0.01, *** p.adj < 0.001"
    ) +
    theme_bw() +
    theme(
        legend.position  = "top",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey90"),
        strip.text       = element_text(face = "bold")
    )

p  # display plot

# ---- Save plot ----
ggsave(
    filename = "~/PIMMSgit/plots/Within_vs_between_individual.pdf",
    plot     = p,
    width    = 8,
    height   = 6
)