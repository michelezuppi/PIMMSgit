#Quick analysis of tsv cluster files

library(tidyverse)
library(ggpubr)
setwd("~/PIMMSgit/data/results/gene_clusters/final_outputs/")
files <- list.files(pattern = "*.tsv", full.names = TRUE)
list.files()
view(files)
for (file in files) {
    
    # Extract file name for plot titles
    file_name <- tools::file_path_sans_ext(basename(file))
    
    # Read data
    data <- read_tsv(file, col_names = c("representative", "member")) %>%
        mutate(
            individual = str_extract(member, "M\\d+"),
            timepoint  = str_extract(member, "(?<=_)(D?-?\\d+)(?=\\|)")
        ) %>%
        filter(timepoint != "D004")
    
    # --- Plot 1 — Distribution of individuals per cluster ---
    p1 <- data %>% 
        group_by(representative) %>% 
        filter(!duplicated(individual)) %>%
        count() %>% 
        ggplot(aes(x = n)) +
        geom_histogram(binwidth = 1, fill = "black", color = "white") +
        geom_text(
            stat = "bin", 
            aes(label = after_stat(
                paste0(count, " (", round(count / sum(count) * 100, 1), "%)")
            )), 
            binwidth = 1, vjust = -0.5, size = 3
        ) +
        labs(x = "Number of individuals", y = "Count",
             title = paste0(file_name, " — individuals per cluster")) +
        theme_classic()
    
    # --- Define prevalent clusters ---
    prevalent <- data %>% 
        group_by(representative) %>% 
        filter(!duplicated(individual)) %>% 
        count() %>%
        filter(n == 4) %>% 
        ungroup() %>% 
        select(representative) %>% 
        distinct()
    
    # --- Plot 2 — Timepoint distribution for prevalent clusters ---
    p2 <- data %>% 
        filter(representative %in% prevalent$representative) %>% 
        group_by(representative, individual) %>%
        summarise(n_timepoints = n_distinct(timepoint), .groups = "drop") %>%
        group_by(representative) %>%
        summarise(mean_timepoints = mean(n_timepoints), .groups = "drop") %>% 
        ggplot(aes(x = mean_timepoints)) +
        geom_histogram(binwidth = 0.5, fill = "black", color = "white") +
        geom_text(
            stat = "bin",
            aes(label = after_stat(
                paste0(count, "\n(", round(count / sum(count) * 100, 1), "%)")
            )),
            binwidth = 0.5, vjust = -0.5, size = 3
        ) +
        labs(x = "Mean number of timepoints", y = "Count",
             title = paste0(file_name, " — timepoints per prevalent cluster")) +
        theme_classic()
    
    # --- Combine and save ---
    combined <- ggarrange(p1, p2, ncol = 2)
    ggsave(paste0("~/PIMMSgit/plots/", file_name, "_summary.pdf"), combined, width = 12, height = 5)
    
    cat("Done:", file_name, "| Total clusters:", n_distinct(data$representative),
        "| Prevalent:", nrow(prevalent), "\n")
}
