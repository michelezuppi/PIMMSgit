

reads_count_WGS  <- read_tsv("~/PIMMSgit/data/combined_trimmed_WGS.tsv")
read_count_VLP  <- read_tsv("~/PIMMSgit/data/combined_trimmed_VLP.tsv")

reads_count_WGS  <- reads_count_WGS %>% 
    mutate(method = "WGS")

read_count_VLP  <- read_count_VLP %>% 
    mutate(method = "VLP")

read_count_combined <- bind_rows(read_count_VLP, reads_count_WGS)

reads_count_patient <- 
    read_count_combined %>% 
    mutate(individual = case_when(
            str_detect(file, "M049") ~ "M049",
            str_detect(file, "M053") ~ "M053",
            str_detect(file, "M060") ~ "M060",
            str_detect(file, "M064") ~ "M064"
    )) %>% 
    filter(!str_detect(file, "D004")) %>% 
    mutate(timepoint = str_extract(file, "(?<=_)[^_]+(?=_trimmed)")) %>% 
    mutate(timepoint = factor(timepoint, levels = c(-7, 0, 2, 4, 7, 14, 28))) %>% 
    group_by(individual, timepoint, method) %>% 
    summarise(depth = sum(num_seqs), .groups = "drop") %>% 
    view()

reads_count_patient  %>% 
    ggplot(aes(x = individual, y = depth, fill = individual))+    
    geom_boxplot(alpha = 0.5, outlier.alpha = 0)+
    geom_quasirandom(shape = 21, size = 3)+
    ggrepel::geom_text_repel(aes(label = timepoint),
                    position = position_quasirandom(),
                    size = 6, max.overlaps = Inf) +
    theme_classic2()+
    scale_fill_viridis_d(option = "turbo")+
    theme(legend.position = "none")+
    facet_wrap(~method)+
    labs(y = "Sequencing depth \n (trimmed reads)", x = "Individual")

ggsave(filename = "~/PIMMSgit/plots/selective_pressure_analysis/patients_sequencing_depth.pdf", plot = last_plot())


#Combine sequencing depth with SNVs count
load("~/PIMMSgit/data/Rdata/PIMMs_vOTUs.Rdata")
load("~/PIMMSgit/data/Rdata/breadth_cutoff.Rdata")
gene_info <- read_tsv("~/PIMMSgit/data/phages/instrain_combined_output/combined_gene_info.tsv")

clean_coverage <- quantile(gene_info$coverage, na.rm = TRUE, probs = c(0.95))

gene_info_prefiltered <- 
    gene_info %>% 
    filter(coverage >= 10, coverage < clean_coverage) %>% 
    filter(scaffold %in% breadth_cutoff$genome) %>% 
    filter(scaffold %in% PIMMs_vOTUs$representative) %>% 
    filter(!is.na(pNpS_variants))


SNVs_count_sample <- 
    gene_info_prefiltered %>% 
    group_by(method, individual, timepoint) %>% 
    summarise(SNV_count = sum(SNV_count), .groups = "drop") %>% 
    filter(timepoint != "D004") %>% 
    view()

SNVs_read  <- 
    reads_count_patient %>% 
    left_join(SNVs_count_sample)

SNVs_read$depth_s <- scale(SNVs_read$depth)

model <- glmmTMB(SNV_count ~ depth_s + (1 | individual),
                 family = nbinom2, data = SNVs_read)

model$sdr$pdHess
any(is.nan(vcov(model)$cond))
diagnose(model)

simulateResiduals(model) %>% 
    plot()

pvalue  <- tidy(model) %>% 
    filter(term == "depth_s") %>% 
    pull(p.value)




SNVs_read %>% 
    ggplot(aes(x = log10(SNV_count), y = log10(depth)))+
    geom_point()+
    geom_smooth(method = "lm")+
    theme_classic2()+
    annotate("text", x = 4, y = 8.5,
         label = paste0("pvalue = ", signif(pvalue, 2)))+
    labs(x = "SNVs count \n(log10)", y = "Sequencing depth \n (log10)", caption = "glmmTMB(SNV_count ~ depth_s + (1 | individual), family = nbinom2)")

 ggsave(filename = "~/PIMMSgit/plots/selective_pressure_analysis/depth_SNVs.pdf")


