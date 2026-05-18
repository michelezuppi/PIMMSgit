#Protein cluster analysis within and between individuals
rm(list = ls())
library(tidyverse)
library(ggpubr)
library(ggbeeswarm)
library(rstatix)
library(twosamples)
library(patchwork)
setwd("~/PIMMSgit/data/results/gene_clusters/nr_output/")

#What am I prevalence_analysising to do? See if there is a higher number of proteins identified within the same individual compared to between
blosum_clusters <- read_tsv("protein_clusters.tsv", col_names = c("representative", "member"))

# view(blosum_clusters)

blosum_clusters_no_donor  <- 
    blosum_clusters %>% 
    mutate(
            individual = str_extract(member, "M\\d+"),
            timepoint  = str_extract(member, "(?<=_)(D?-?\\d+)(?=\\|)")
        ) %>%
        filter(timepoint != "D004")
        
# view(blosum_clusters_no_donor)

# length(unique(blosum_clusters_no_donor$representative)) #106020 clusters

# singleton <- 
#     blosum_clusters_no_donor %>% 
#     group_by(representative) %>% 
#     filter(n_distinct(member) == 1) %>% #Select singleton
#     ungroup()
# # length(unique(singleton$representative)) #49575 clusters were singletons
# # view(singleton)

clusters <- 
    blosum_clusters_no_donor %>% 
    group_by(representative) %>% 
    filter(n_distinct(member) > 1) %>% #Remove singleton
    ungroup()
# length(unique(clusters$representative)) #56445 clusters were found in more than one sample
# view(clusters)

###Identify clusters that were found at at all timepoints after transplant, including baseline

    persistent_6 <- 
        clusters %>% 
        filter(timepoint != -7) %>% #We remove this timepoint as we cannot use it for comparison with donor clusters
        group_by(representative, individual) %>% 
        filter(!duplicated(timepoint)) %>% 
        group_by(representative, individual) %>% 
        summarise(n_timepoints = n_distinct(timepoint), .groups = "drop") %>% 
        filter(n_timepoints == 6) #Select only the ones identified at all 6 timepoints
        
# length(unique(persistent_6$representative)) #2965 clusters were found at all timepoints
# view(persistent_6)

#Define the number of individuals persistent clusters were found in and select only the ones found in more than one. 
    persistent_prevalent_gene_clusters_blosum62  <- 
        persistent_6 %>% 
        select(-n_timepoints) %>% 
        group_by(representative) %>% 
        summarise(n_individuals = n_distinct(individual), .groups = "drop") %>% 
        filter(n_individuals != 1) #Select only the one found in at least 2 individuals


# length(unique(persistent_prevalent_gene_clusters_blosum62$representative)) #1170 persistent (6) gene clusters were found in at least 2 individuals
# view(persistent_prevalent_gene_clusters_blosum62)


#Create a list of names of gene_cluster_representatives that were found across all 6 timepoints (0-28 days), excluding donor, in more than one individual.

setwd("~/PIMMSgit/data/Rdata")
save(persistent_prevalent_gene_clusters_blosum62, file = "persistent_prevalent_gene_clusters_blosum62.Rdata")

##### Plot the two distribution together
# Plot 1: n_timepoints per representative per individual
p1_data <- clusters %>%
    filter(timepoint != -7) %>%
    group_by(representative, individual) %>%
    filter(!duplicated(timepoint)) %>%
    group_by(representative, individual) %>%
    summarise(n_timepoints = n_distinct(timepoint), .groups = "drop") #The count of clusters present at 6 timepoints is different from what calculated previously 
    #(length(unique(persistent_6$representative))) because here we are counting the number of timepoints per representative per individual,
    # while previously we were counting the number of representatives that were present at 6 timepoints.
    #There can be multiple observation for the same unique clusters at 6 timepoints, if they are found in different individuals.


p1 <- p1_data %>%
    ggplot(aes(x = n_timepoints)) +
    geom_histogram(binwidth = 1, fill = "black") +
    stat_bin(
        binwidth = 1,
        geom = "text",
        aes(label = after_stat(ifelse(count > 0, count, ""))),
        vjust = -0.5,
        color = "black",
        size = 3.5
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    labs(
        title = "Timepoints per cluster per individual",
        x = "Number of timepoints",
        y = "Count"
    ) +
    theme_classic2()

# Plot 2: n_individuals per representative (persistent clusters only)
p2_data <- persistent_6 %>%
    select(-n_timepoints) %>%
    group_by(representative) %>%
    summarise(n_individuals = n_distinct(individual), .groups = "drop")

p2 <- p2_data %>%
    ggplot(aes(x = n_individuals)) +
    geom_histogram(binwidth = 1, fill = "black") +
    stat_bin(
        binwidth = 1,
        geom = "text",
        aes(label = after_stat(ifelse(count > 0, count, ""))),
        vjust = -0.5,
        color = "black",
        size = 3.5
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    labs(
        title = "Individuals per persistent cluster (6 timepoints)",
        x = "Number of individuals",
        y = "Count"
    ) +
    theme_classic2()

# Combine
p1 + p2

ggsave("~/PIMMSgit/plots/persistent_prevalent_gene_clusters_blosum62.png", plot = last_plot(), width = 10, height = 5)

#Create a dataframe containing the representatives and the clustering genes within the individuals and timepoints of interest. 

view(persistent_6)
# view(blosum_clusters_no_donor)

ppgcb_clusters <- 
    persistent_6 %>% 
    left_join(., blosum_clusters_no_donor  %>% 
        filter(timepoint != -7), #Remove the -7 ones
        by = c("representative", "individual")) %>% 
        bind_rows(filter(., !representative %in% member) %>% # For some clusters, the representative is not found in the member column
        mutate(member = representative) %>% 
        distinct(representative, .keep_all = TRUE)) %>% 
        #Redo the individual and timepoints, as we have added some new contigs
        select(representative, member) %>% 
        mutate(
            individual = str_extract(member, "M\\d+"),
            timepoint  = str_extract(member, "(?<=_)(D?-?\\d+)(?=\\|)")
        )

view(ppgcb_clusters)

setwd("~/PIMMSgit/data/Rdata")
save(ppgcb_clusters, file = "ppgcb_clusters.Rdata")
#Just a reminder, we now have some clusters in the pppgdb that are found at -7 and they results as if found in only 1 timepoints.
#We are keeping them to keep the representative gene of this cluster. These are marked as found at the timepoint -7. 

#Create a list of the genes id clustering within our clusters of interest

setwd("~/PIMMSgit/data/Rdata")

ppgcb_clusters %>% 
    pull(member) %>% 
    unique() %>% 
    writeLines("member_genes.txt")


#############################################################################################################################################################


#Check the median timepoints in which donor genes were found#####

blosum_clusters_donor  <- 
    blosum_clusters %>% 
    mutate(
            individual = str_extract(member, "M\\d+"),
            timepoint  = str_extract(member, "(?<=_)(D?-?\\d+)(?=\\|)")
        ) %>%
    group_by(representative, individual) %>% 
    filter(any(timepoint == "D004") & !any(timepoint == 0) & !any(timepoint == -7)) %>% #Select only the clusters supposedely engrafted
    filter(any(timepoint != "D004"))  %>% #Remove the clusters within donor samples
    ungroup()

view(blosum_clusters_donor)

length(unique(blosum_clusters_donor$representative)) #5198 clusters were found in donor after the transplant

blosum_clusters_donor %>% 
    group_by(representative) %>% 
    filter(n_distinct(member) > 1) %>% #Remove singleton
    ungroup() %>% 
    mutate(timepoint = if_else(timepoint == "D004", 0, as.numeric(timepoint)))  %>%  #Convert the donor timepoint to numeric
    group_by(representative, individual) %>% 
    summarise(n_timepoints = n_distinct(timepoint), .groups = "drop") %>%
    pull(n_timepoints) %>% 
    #quantile()    #   0%  25%  50%  75% 100%  #Median value = 2
                   #   2    2    2    3    6 
    {mean(. == 6) * 100} #4.74% of clusters were found at 6 timepoints

blosum_clusters_donor %>% 
    group_by(representative) %>% 
    filter(n_distinct(member) > 1) %>% #Remove singleton
    ungroup() %>% 
    select(representative) %>% 
    distinct() %>% 
    nrow() #5198 donor clusters were found in recipients 

blosum_clusters_donor %>% 
    group_by(representative) %>% 
    filter(n_distinct(member) > 1) %>% #Remove singleton
    ungroup() %>% 
    mutate(timepoint = if_else(timepoint == "D004", 0, as.numeric(timepoint)))  %>%  #Convert the donor timepoint to numeric
    group_by(representative, individual) %>% 
    summarise(n_timepoints = n_distinct(timepoint), .groups = "drop") %>%
    filter(n_timepoints == 6) %>% 
    select(representative) %>% 
    distinct() %>% 
    nrow() #295 donor clusters were found in all timepoints. 