#Functional annotation of ppgcb_clusters from pharokka

rm(list = ls())
library(tidyverse)
library(ggpubr)

setwd("~/PIMMSgit/data/Rdata")
load("ppgcb_clusters.Rdata")

#Load the pharokka output

# view(ppgcb_clusters)
pharokka_output <- read_tsv("~/PIMMSgit/data/results/gene_clusters/pharokka_cds_final_merged_output.tsv")

pharokka_output$gene  <- gsub("CDS_0*", "", pharokka_output$gene) #Change the gene name to match the member name in the ppgcb_clusters dataframe. 

#For each cluster, select the most frequent annotation
clusters_annotations <- 
    ppgcb_clusters %>%
    left_join(pharokka_output, by = c("member" = "gene")) %>% 
    select(c(representative, member, annot, category)) %>% 
    group_by(representative, annot, category) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(representative) %>%
    mutate(percentage = n / sum(n) * 100) %>%
    slice_max(percentage, n = 1, with_ties = FALSE) %>%
    select(representative, annot, category) %>% 
    ungroup()

    view(clusters_annotations)

ppgcb_annotated <- 
    ppgcb_clusters %>%
    left_join(clusters_annotations, by = "representative")

view(ppgcb_annotated)

ppgcb_annotated %>%
    group_by(representative) %>%
    filter(!duplicated(representative)) %>%
    ungroup() %>%
    count(category) %>%
    mutate(percentage = n / sum(n) * 100) %>%
    ggplot(aes(y = category, x = percentage, fill = category)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = sprintf("%.1f%%", percentage)), hjust = -0.1, size = 3.5) +
    scale_fill_viridis_d() +
    theme_classic2() +
    theme(legend.position = "none") +
    labs(x = "Percentage (%)", y = NULL)

ggsave("~/PIMMSgit/plots/ppgcb_cluster_annotation_percentage.png", plot = last_plot(),width = 8, height = 5)


# view(pharokka_output)
#Check how many representatives were annotated in distinct ways

annotation  <- 
    ppgcb_annotated %>% 
    group_by(representative) %>% 
    mutate(n_cat = n_distinct(category)) %>% 
    ungroup() %>% 
    filter(n_cat > 1)

length(unique(annotation$representative)) # 622 genes have multiple annotations, out of 2965 total genes
view(annotation)

percentages  <-
    annotation %>%
    group_by(representative, category) %>%
    summarise(n = n(), mean_score = mean(score), .groups = "drop") %>%
    group_by(representative) %>%
    mutate(percentage = n / sum(n) * 100,
    score_proportion = min(mean_score) / max(mean_score)) %>% 
    group_by(representative) %>% 
    mutate(count_categories = n_distinct(category)) %>%
    ungroup()


view(percentages)


    

view(several_annotations)

length(unique(several_annotations$representative)) # 71
