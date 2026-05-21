#Functional annotation of ppgcb_clusters from pharokka

rm(list = ls())
library(tidyverse)
library(ggpubr)
library(patchwork)

setwd("~/PIMMSgit/data/Rdata")
load("ppgcb_phage_clusters.Rdata")

#Load the pharokka output

view(ppgcb_phage_clusters)
pharokka_output <- read_tsv("~/PIMMSgit/data/results/pharokka_results/pharokka_proteins_summary_output.tsv")

head(pharokka_output)
#For each cluster, select the most frequent annotation

clusters_annotations <- 
    ppgcb_phage_clusters %>%
    left_join(pharokka_output, by = c("member" = "ID")) %>%
    select(c(representative, member, annot, category)) %>%   
    group_by(representative, annot, category) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(representative) %>%
    mutate(percentage = n / sum(n) * 100) %>%
    slice_max(percentage, n = 1, with_ties = FALSE) %>% #Select the annotation with the highest percentage for each representative
    select(representative, annot, category) %>% 
    distinct(representative, annot, category) %>%
    ungroup()
    #We now have the annotations for each cluster

    view(clusters_annotations)

ppgcb_annotated <- 
    ppgcb_phage_clusters %>%
    left_join(clusters_annotations, by = "representative")

view(ppgcb_annotated)

p1  <- ppgcb_annotated %>%
    group_by(representative) %>%
    filter(!duplicated(representative)) %>%
    ungroup() %>%
    count(category) %>%
    mutate(percentage = n / sum(n) * 100) %>%
    ggplot(aes(y = category, x = percentage, fill = category)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = sprintf("%.1f%%", percentage)), hjust = -0.1, size = 3.5) +
    scale_fill_viridis_d(direction = -1) +
    theme_classic2() +
    theme(legend.position = "none") +
    labs(x = "Percentage (%)", y = NULL)

p2  <-    
   ppgcb_annotated %>%
        select(category, n_individuals, n_timepoint)  %>% 
        group_by(category, n_individuals) %>%
        slice_max(n_timepoint, n = 1, with_ties = FALSE) %>%
        ungroup() %>% 
        mutate(category = factor(category)) %>% 
        complete(category, n_individuals) %>%  # fills missing combos with NA
        ggplot(aes(x = category, y = n_timepoint, fill = category)) +
        geom_col()+
        scale_fill_viridis_d(direction = -1) +
        theme_classic2() +
        facet_wrap(~n_individuals, scales = "free", 
           labeller = labeller(n_individuals = function(x) paste("Number of individuals:", x)))+
        theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))+
        scale_y_continuous(breaks = seq(0, 6, by = 1))+
        labs(x = NULL, y = "Number of timepoints")



ggsave("~/PIMMSgit/plots/ppgcb_cluster_annotation_percentage.pdf", plot = last_plot())

#Check the annotations for the clusters that are present in 3 individuals


view(ppgcb_annotated)

  p3  <- ppgcb_annotated %>%    
    filter(category != "unknown function") %>% 
    group_by(representative, annot, n_individuals) %>%
    filter(!duplicated(representative)) %>%
    ungroup() %>%
    group_by(annot, category, n_individuals) %>%
    slice_max(n_timepoint, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    ggplot(aes(x = annot, y = n_timepoint, fill = category)) +
    geom_col() +
    coord_flip()+
    theme_classic2() +
    scale_fill_viridis_d(begin = 0.19, direction = -1)+
        facet_grid(~n_individuals, 
           labeller = labeller(n_individuals = function(x) paste("Number of individuals:", x)))

p1 / p2 / p3

ggsave("~/PIMMSgit/plots/ppgcb_cluster_annotation_percentage.pdf", plot = last_plot())
