#Check the ANI
library(tidyverse)
library(ggpubr)
library(ggbeeswarm)
library(rstatix)
library(twosamples)
library(conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
set.seed(1444)
setwd("~/PIMMSgit/data/phages/unified/vclust")
ANI <- read_tsv("_ani.tsv")
load("~/PIMMSgit/data/Rdata/vOTUs_clusters.Rdata")

view(vOTUs_clusters)
view(ANI_analysis)

vOTUs_clusters <- 
    vOTUs_clusters  %>% 
    mutate(object = str_remove(object, "\\|\\|.*")) %>% #Clean the sequence name
    mutate(cluster = str_remove(cluster, "\\|\\|.*")) #Clean the sequence name

   
ANI_analysis <- 
    ANI %>% 
    select(c(query, reference, tani)) %>% 
    filter(tani > 0.95) %>% 
    mutate(query = str_remove(query, "\\|\\|.*")) %>% #Clean the sequence name
    mutate(reference = str_remove(reference, "\\|\\|.*")) %>% #Clean the sequence name
    mutate(
        individual_q = str_extract(query, "(?<=_)[^_]+(?=_[^_]+$)"),
        timepoint_q  = str_extract(query, "[^_]+$")) %>% 
    mutate(
        individual_r = str_extract(reference, "(?<=_)[^_]+(?=_[^_]+$)"),
        timepoint_r  = str_extract(reference, "[^_]+$"))

    view(ANI_analysis)

vOTUs_ANI  <- ANI_analysis  %>% 
        mutate(pair_key = map2_chr(query, reference, \(q, r) paste(sort(c(q, r)), collapse = "|||"))) %>% # Create a sorted pair key so A-B and B-A get the same key
        distinct(pair_key, .keep_all = TRUE) %>% # Keep only one row per unique pair
        left_join(vOTUs_clusters, by = c("query" = "object")) %>% # Try joining vOTU on query
        left_join(vOTUs_clusters, by = c("reference" = "object"), suffix = c("_q", "_r")) %>% # Try joining vOTU on reference
        mutate(vOTU = coalesce(cluster_q, cluster_r)) %>% # Take whichever join succeeded
        group_by(vOTU, individual_q, individual_r, timepoint_q, timepoint_r) %>% # One tANI per vOTU per individual pair
        summarise(tani = mean(tani), .groups = "drop")

view(vOTUs_ANI)



FMT_analysis <- 
    vOTUs_ANI %>% 
    pivot_longer(cols = c("timepoint_q", "timepoint_r"), values_to = "timepoint", names_to = NULL) %>% 
    pivot_longer(cols = c("individual_q", "individual_r"), values_to = "individual", names_to = NULL) %>% 
    group_by(vOTU) %>% 
    filter(any(timepoint == "D004" & !all(timepoint == "D004"))) %>% #Select the vOTUs that had been found in donors and recipients
    select(c(vOTU, tani, timepoint, individual)) %>% 
    ungroup() %>% 
    group_by(vOTU, individual)  %>% 
    mutate(FMT_orientation = if_else(any(timepoint %in% c(-7, 0)), "Before", "After")) %>% #Divide on whether they were found within the individual before or after the transplant
    ungroup() %>% 
    distinct(vOTU, FMT_orientation, .keep_all = TRUE)

view(FMT_analysis)

head(FMT_analysis)

ANI_threshold  <- 
  FMT_analysis  %>% 
    distinct(vOTU, FMT_orientation, tani) %>% 
    mutate(threshold = if_else(tani >= 0.995, "Yes", "No"))  

view(ANI_threshold)

    

chisq_results  <- ANI_threshold  %>% 
    group_by(FMT_orientation, threshold) %>% 
    count() %>% 
    ungroup() %>% 
    pivot_wider(names_from = threshold, values_from = n, values_fill = 0) %>% 
    select(-1) %>% 
    chisq.test()

    chisq_results$p.value


p_label <- paste0("χ² p = ", round(chisq_results$p.value, 5))

ANI_threshold %>% 
    mutate(FMT_orientation = factor(FMT_orientation, 
                                    levels = c("Before", "After"))) %>%  
    group_by(FMT_orientation, threshold) %>% 
    count() %>% 
    group_by(FMT_orientation) %>% 
    mutate(tot = sum(n)) %>% 
    ungroup() %>% 
    mutate(prop = n / tot * 100) %>% 
    ggplot(aes(x = FMT_orientation, y = prop, fill = threshold)) +
    geom_col(color = "black") +
    scale_fill_viridis_d(option = "magma", begin = 0.3, end = 0.8)+
    theme_classic2()+
    labs(y = "Donor vOTUs in recipients", fill = "ANI >= 99.5%", x = "Before or after FMT")+
    annotate("text",
             x     = 1.5,                               # centred between the two bars
             y     = 105,                               # just above the bars
             label = p_label,
             size  = 4,
             fontface = "italic")

ggsave(filename = "~/PIMMSgit/plots/ANI_engraftment_analysis.pdf", plot = last_plot())


