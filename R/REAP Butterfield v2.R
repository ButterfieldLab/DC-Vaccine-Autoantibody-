
#Libraries

{
  
  # Some issues for solving rlang problems
  # I went to the my library 
  # C:\Users\pmunson\AppData\Local\Programs\R\R-4.2.2\library
  # and just deleted the old rlang folder 
  #install.packages("umap")
  #packageVersion("rlang")
  #install.packages("rlang")
  
  #remove.packages(rlang)
  #devtools::install_github("r-lib/rlang", build_vignettes = TRUE)
  
  ### this is how i can detach all packages
  #file.edit(file.path("~", ".Rprofile"))
  #install.packages("rlang")
  #install.packages("rlang", version="1.1.0")
  
  #install.packages("writexl")
  #install.packages("corrr")
  #install.packages("ggthemes")
  #install.packages("extrafont")
  #  install.packages("gridExtra")
  #  install.packages("VennDiagram")
  #  install.packages("gprofiler2")
  #  install.packages("ggtext")
  #  install.packages("tm")
  #  install.packages("hrbrthemes")
  #  install.packages("RColorBrewer")
  #  install.packages("ggpubr")
  #  install.packages("stringi") #I am trying to use this for multiple gsub arguments
  #  install.packages("purrr")
  #  install.packages("ggbreak")
  #  install.packages("stringr")
  #  install.packages("ggprism")
  #install.packages("umap")
  #install.packages("circlize")
  
  
  
  
  library(rlang)  
  library(dplyr)
  library(tidyr)
  library(tidyverse)
  library(readxl)
  library(stringr)
  library(writexl)
  library(corrr)
  library(Rcpp)
  library(ggrepel)
  library(data.table)
  library(readr)
  library(ggthemes)
  library(extrafont)
  library(gridExtra)
  library(VennDiagram)
  library(gprofiler2)
  library(writexl)
  library(ggtext)
  library(tm)
  #library(hrbrthemes)
  library(RColorBrewer)
  library(ggpubr)
  library(stringi) #I am trying to use this for multiple gsub arguments
  library(purrr)
  library(ggbreak)
  library(stringr)
  library(ggprism)
  
  ## I have had issues with annotating my volcano plots and it has to do with 
  ## Some packages overwriting functions within ggplot
  
  
  #detach("package:ggvenn", unload=TRUE)
  #detach("package:proustr", unload=TRUE)
  #detach("package:tm", unload=TRUE)
  #detach("package:ggpubr", unload=TRUE)
  #detach("package:ggrepel", unload=TRUE)
  #detach("package:ggplot2", unload=TRUE)
  
  
  library(ggvenn)
  library(ggplot2)
  library(ggprism)
  library(umap)
  ### for trying to add correct multiple comparison statistics
  library(PMCMRplus)
  
  
  ## heat map libraries
  library(cluster)
  library(digest)
  library(circlize)
  #install.packages("ComplexHeatmap") 
  library(ComplexHeatmap)
  library(RColorBrewer)
  #install.packages("pheatmap")
  library(rstatix)
  
  #  if (!require("BiocManager", quietly = TRUE))
  # install.packages("BiocManager")
  
  #BiocManager::install("ComplexHeatmap")
  #library(ComplexHeatmap)
}

### Butterfield REAP data
### Loading in the metadata and raw data
{
  metadata_df <- read_excel("C:/Users/pmunson/OneDrive - Nurix Therapeutics, Inc/Desktop/Autoantibodies/cancer vaccine Clinical Metadata.xlsx")
  raw_df <- read_csv("C:/Users/pmunson/OneDrive - Nurix Therapeutics, Inc/Desktop/Autoantibodies/REAP_cancer_vaccine_butterfield.csv")
  annotation_df <- read_csv("C:/Users/pmunson/OneDrive - Nurix Therapeutics, Inc/Desktop/Autoantibodies/library_annotation.csv")
}


####
## Cleaning metadata
{
  metadata_df_v2 <- metadata_df %>%
    ### Renaming the columns and patient information to be consistent with the data
    rename(Subject_ID = patient_numbers) %>%
    mutate(Subject_ID = gsub("pt_", "PT_", Subject_ID)) %>%
    mutate(Subject_ID = gsub("PT_H", "H", Subject_ID)) %>%
    ### Removing the HC, these are subjects that only exist in the PhiP seq data
    filter(!grepl("^HC", Subject_ID)) %>%
    ### Renaming cancer to HD and Melanoma
    mutate(Cancer = ifelse(Cancer == "None", "HD", ifelse(Cancer == "Yes", "Melanoma", Cancer)))
  
  
  
  
  ## It appears in the data it is KIR2DL2, but in the annotation table it is KI2L2
  ## There is also something called LOC348174 that isn't found in the annotation table, everything else appears accounted for
  
  
  ## Here I am showing what is present in the data, yet absent in the annoataion data frame
  setdiff(raw_df %>% select("...1") %>% unique() %>% pull(),annotation_df$protein)
  
  ## Changing the name of the protein to reflect its call in the data frame
  ## LOC348174 looks to be some kind of vector, potentially an artifact
  annotation_df_v2 <- annotation_df %>%
    mutate(protein = ifelse(protein == "KI2L2", "KIR2DL2", protein)) %>%
    rename(Protein = protein) %>%
    unique()
  
  ## The annotation table doesn't appear to be uniquely identified 
  problematic_annotations <- annotation_df_v2 %>%
    group_by(Protein) %>%
    count() %>%
    arrange(desc(n)) %>%
    filter(n > 1) %>% pull(Protein)
  
  
  raw_meta_df <- raw_df %>%
    rename(Protein = "...1") %>%
    pivot_longer(., cols = -c(Protein)) %>%
    mutate(Subject_ID = gsub("^([^_]+_[^_]+).*", "\\1", name)) %>%
    mutate(Visit_Name_time = gsub("^([^_]+_[^_]+_)(.*)", "\\2", name)) %>%
    
    ### Making the healthy donor timepoints equal to Day 0 
    mutate(Visit_Name_time = ifelse(str_starts(Subject_ID, "HD"), "DAY_0", Visit_Name_time)) %>% 
    left_join(.,annotation_df_v2) %>%
    rename(Subject_Visit = name) %>%
    left_join(., metadata_df_v2)
  
  
  raw_statistics_df <- raw_meta_df %>%
    filter(Visit_Name_time == "DAY_0") %>%
    select(Cancer, Protein, value) %>%
    group_by(Protein) %>%
    #### Calculating statistics
    summarise(
      lfc = log2(mean(value[Cancer == "Melanoma"]) / mean(value[Cancer == "HD"])),
      median_diff = median(value[Cancer == "Melanoma"])-median(value[Cancer == "HD"]),
      mean_diff = mean(value[Cancer == "Melanoma"])-mean(value[Cancer == "HD"]),
      max_value_cancer = max(value[Cancer == "Melanoma"]),
      max_value_hd = max(value[Cancer == "HD"]),
      p_value = t.test(value[Cancer == "Melanoma"], value[Cancer == "HD"], na.action = na.omit)$p.value
    ) %>%
    mutate(max_diff = max_value_cancer - max_value_hd)
  
  ### The infinite log fold change means it's detected in cancer patients, but in none of the healthy donors
  unique_cancer_patients <- raw_statistics_df %>%
    filter(lfc == Inf) %>% pull(Protein)
}

## Reapfigures
{





#### Volcano plot - FIG1A
### adapting volcano plots to focus on proteins with sufficiently high signals
{
  

  

volcano_plot <-  raw_statistics_df %>%
  ## removing the Epitope
  filter(!grepl("_Epitope", Protein)) %>%
  left_join(., annotation_df_v2) %>%
  #filter(lfc != Inf & lfc != -Inf) %>%
  ggplot(., aes(x = lfc, y = -log(p_value), color = max_diff))+
  geom_point(size =2) +
  theme_prism() +
  #xlim(-20,20) +
  #ylim(0,10.9) +
  ylab("-Log(p value)") +
  xlab("Log2 fold change") +
  #ggtitle("Baseline") +
    scale_color_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0) +
    #geom_text(x = Inf, y = Inf, label = "Up in Melanoma", hjust = 1.5, vjust = 1, color = "red", size =5) +
    #geom_text(x = -Inf, y = Inf, label = "Up in Healthy Donors", hjust = -0.75, vjust = 1, color = "blue", size = 5) +
    theme(legend.position = "none") 
      #geom_text_repel(data = subset(raw_statistics_df %>%
      #                                ## removing the Epitope
      #                                filter(!grepl("_Epitope", Protein)), 
      #                              #`p_value` < 0.25  & lfc >1 &
      #                              `p_value` < 0.25  &
      #                                lfc > 1 &
      #                                max_diff > 1), 
      #                aes(label = Protein),
      #                size = 3,
      #                direction = "both",
      #                max.overlaps = 20,min.segment.length = 0.01, box.padding = 0.25, color = "black") +
  
#setwd("C:/Users/pmunson/OneDrive - Nurix Therapeutics, Inc/Desktop/Autoantibodies/Graphs") 
#ggsave(filename = "volcano_plot.png", plot = volcano_plot, width = 4, height = 4, dpi = 600)





}

### Looking at individual proteins



# IL1F9
#IFNA6
#CCR5_Epitope-3-N
#CD164L2
#IFNA13
#CTLA4
#PDCD1
#IL6
#TNFAIP6
#SLC2A2
#CDON
#MMP2
#IL18RAP
#IFNA6
raw_meta_df %>%
  mutate(Responder = case_when(
    Best_Clinical_outcome %in% c("PR", "SD", "NED1") ~ "R",
    Best_Clinical_outcome %in% c("NED2", "PD") ~ "NR",
    TRUE ~ NA_character_  # Handle other cases if any
  )) %>%
  #filter(Visit_Name_time == "DAY_0") %>%
  filter(Protein == "PGLYRP1") %>%
  #filter(!is.na(Responder)) %>%
  ggplot(., aes(x = Visit_Name_time, y = value, color = Best_Clinical_outcome, group = Subject_ID)) +
  geom_point() +
  geom_line() +
  facet_wrap(~Responder + Protein) 


raw_meta_df  %>%
  filter(Visit_Name_time == "DAY_0") %>%
  #filter(Subject_ID == "PT_10") %>%
  arrange(desc(value)) %>%
  filter(Subject_ID == "PT_10") 
  

## if they have an infinite fold change that menas its in the melanomna but not
## in the HD
{
unique_to_melanoma <- raw_statistics_df %>%
  filter(!grepl("_Epitope", Protein)) %>%
  filter(lfc == Inf) %>% pull(Protein)

raw_statistics_df %>%
  filter(!grepl("_Epitope", Protein)) %>%
  filter(lfc == Inf) %>%
  arrange(desc(max_diff)) %>%
  print(n = 30)
}


### I want to filter my unique in melanoma proteins, 
## to ones that are present in at least 3 subjects

# Heatmap - FIG1B
{


  
### Finding DE bound proteins, excluding ones with unknown functions
  ## ACRV1 testis specific 
  ## AMEYL (y linked, tooth enamel)
  ## APOO apolipoprotein family, metabolism #Metabolism
  ## CDH19, a cadherin "loss of caherins may be associated with cancer formation"
  ## LRP11 (low density lipoportein receptor), #Metabolism
  ## OBP2B , olfactory, nose met? 
  ## PRAP1, binds lipids (#Metabolism), mainly expressed in epithelial cells..
  ## PTPRR ,,, #WNT
  ## SSPN #ECM
  # TMEM119 bone, bone marrow metastasis?
  
  
  volcano_proteins <- subset(raw_statistics_df %>%
                               ## removing the Epitope
                               filter(!grepl("_Epitope", Protein)), 
                             #`p_value` < 0.25  & lfc >1 &
                             `p_value` < 0.25  &
                               lfc > 1 &
                               max_diff > 1) %>% select(Protein)
  
  ### Manually assigning functions to potential proteins of interest
  unknown_proteins <- volcano_proteins %>%
    left_join(., annotation_df_v2) %>%
    #filter(reactome_group != "Other / Unknown") %>%
    mutate(annotation_and_reactome = ifelse(is.na(annotation_group), reactome_group, annotation_group)) %>%
    select(Protein, annotation_and_reactome) %>%
    filter(annotation_and_reactome == "Other / Unknown") %>%
    pull(Protein)
  

  ### select desired order of proteins based on pathway similarites on the heat map
  
pathway_levels <- volcano_proteins %>%
    left_join(., annotation_df_v2) %>%
    #filter(reactome_group != "Other / Unknown") %>%
    mutate(annotation_and_reactome = ifelse(is.na(annotation_group), reactome_group, annotation_group)) %>%
    select(Protein, annotation_and_reactome) %>%
    #filter(annotation_and_reactome != "Other / Unknown") %>%
    arrange(annotation_and_reactome) %>%
    select(annotation_and_reactome) %>% mutate(value = 0) %>%
    unique() %>%
    pivot_wider(., names_from = annotation_and_reactome, values_from = value) %>%
    ## define desired order here....
    select(c(Checkpoints.Drugs, Chemokines, Cytokines,"Interferons", "Innate Immune System", 
             "Signal Transduction", "Signaling by GPCR", "SLC-mediated transmembrane transport",
             "Transporters", "Ion channel transport", "GPCRs", "Metabolism", "Metabolism of proteins", "Proteases/Hydrolases",
             "WNT.related", everything())) %>%
    colnames()
  
heatmap_proteins <-volcano_proteins %>%
    left_join(., annotation_df_v2) %>%
    #filter(reactome_group != "Other / Unknown") %>%
    mutate(annotation_and_reactome = ifelse(is.na(annotation_group), reactome_group, annotation_group)) %>%
    select(Protein, annotation_and_reactome) %>%
    #filter(annotation_and_reactome != "Other / Unknown") %>%
    mutate(annotation_and_reactome = factor(annotation_and_reactome, levels = c(pathway_levels))) %>%
    arrange(annotation_and_reactome) %>%
    #filter(Protein != "229E-RBD") %>%
    pull(Protein)
  
heat_map_data <- raw_meta_df %>%
  filter(Visit_Name_time == "DAY_0") %>%
  select(Subject_ID, Protein, value) %>%
  filter(Protein %in% heatmap_proteins) %>%
  pivot_wider(., names_from = Protein, values_from = value) %>%
  select(-Subject_ID) %>%
  #scale() %>%
  t() 

heat_map_data <-heat_map_data[heatmap_proteins,]

colnames_heat_map <- raw_meta_df %>%
  filter(Visit_Name_time == "DAY_0") %>%
  select(Subject_ID, Protein, value) %>%
  filter(Protein %in% heatmap_proteins) %>%
  pivot_wider(., names_from = Protein, values_from = value) %>%
  select(Subject_ID) %>% pull()

colnames(heat_map_data) <- colnames_heat_map

### Annotation data

annotation_heatmap <-raw_meta_df%>%
  #select(Subject_ID, Best_Clinical_outcome, Overall_Survival,PFS, Metastasis_Site, Anti_CTLA_4, Anti_PD_1, Prior_IFN, Prior_IL2, Cancer) 
  select(Subject_ID,Cancer,  Anti_CTLA_4,Anti_PD_1, ) %>%
  mutate(Anti_CTLA_4 = factor(Anti_CTLA_4, levels = c("HD", "0", "1")))


# Reorder the annotation data to match the matrix column names
#annotation_data <- 

annotation_heatmap <- annotation_heatmap[match(colnames_heat_map, annotation_heatmap$Subject_ID), ]


## If this is zero, it means the 
## annotation rows are in the same order as the data columns 
sum(annotation_heatmap$Subject_ID!= colnames_heat_map)

## checking that the heatmap and the annotation table are compatible in dimensions
dim(heat_map_data)[2] == dim(annotation_heatmap)[1]

## I need to make the subject_visit a rowname
annotation_heatmap_df <- as.data.frame(annotation_heatmap)
rownames(annotation_heatmap_df) <-  annotation_heatmap_df[,1] 
annotation_heatmap_df_select <- annotation_heatmap_df %>% select(-Subject_ID)



##This should be zero, means the rows of the annotation perfectly match
# the columns of the data 
sum(rownames(annotation_heatmap_df_select) != colnames(heat_map_data))


colours <- list(
  Cancer = c('HD'= '#ADD8E6', 'Melanoma' = '#FFB6C1'),
  Anti_PD_1 = c('0' = 'black', '1' = '#C3B1E1', 'HD' = 'white' ),
  Anti_CTLA_4 = c('0' = 'black', '1' = '#C3B1E1', 'HD' = 'white')
  )


col_heat_map_annotation <- HeatmapAnnotation(
  df = annotation_heatmap_df_select,
  col = colours, 
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  gap = unit(1, 'mm'))



### color of the cells
col_fun = colorRamp2(c(0,1.5,3,4.5,6), c("black",  "red", "orange", "yellow", "white"))


### Row annotations
### this is copied from the earlier chunk of code to ensure the same order as the data 
color_levels <- c("#fc8c84", "#6454ac", "#5bc0de", "#1abc9c", "#fc6404", "#fcd444", "gray", "black")

Annot <- rownames(heat_map_data) %>%
  data.frame() %>%
  rename(Protein = ".") %>%
  left_join(., annotation_df_v2) %>%
  #filter(reactome_group != "Other / Unknown") %>%
  mutate(annotation_and_reactome = ifelse(is.na(annotation_group), reactome_group, annotation_group)) %>%
  select(Protein, annotation_and_reactome) %>%
  #filter(annotation_and_reactome != "Other / Unknown") %>%
  mutate(annotation_and_reactome = factor(annotation_and_reactome, 
                                          levels = c(pathway_levels))) %>%
  mutate(colors = case_when(
    annotation_and_reactome == "Checkpoints.Drugs" ~ "#fc8c84",
    annotation_and_reactome %in% c("Chemokines", "Cytokines", "Interferons") ~ "#6454ac",
    annotation_and_reactome %in% c("Innate Immune System") ~ "#5bc0de",
    annotation_and_reactome %in% c("Developmental Biology", "ECM organization", "Extracellular matrix organization", "Growth.Factors", "Hemostasis") ~ "#1abc9c",
    annotation_and_reactome %in% c("Ion channel transport",  "SLC-mediated transmembrane transport", "GPCRs", "Signal Transduction", "Signaling by GPCR", "Transporters") ~ "#fc6404",
    annotation_and_reactome %in% c("Metabolism", "Metabolism of proteins", "Proteases/Hydrolases", "WNT.related") ~ "#fcd444",
    annotation_and_reactome == "Viral" ~ "black",
    annotation_and_reactome %in% c("Other / Unknown") ~ "gray",
    TRUE ~ NA_character_ # default case
  )) %>% 
  mutate(colors = factor(colors, levels = color_levels)) 


# Here is defining of Annotation fo CATEGORIES
CATEGORIES = rowAnnotation(annotation_and_reactome = Annot$annotation_and_reactome,
                           col = list(annotation_and_reactome = c(`Checkpoints.Drugs` = "#fc8c84",
                                                   `Chemokines` = "#6454ac",
                                                   `Cytokines` = "#6454ac",
                                                   `Interferons` = "#6454ac",
                                                   `Innate Immune System` = "#5bc0de",
                                                   `Developmental Biology` = "#1abc9c",
                                                   `ECM organization` = "#1abc9c",
                                                   `Extracellular matrix organization` = "#1abc9c",
                                                   `Growth.Factors` = "#1abc9c",
                                                   `Hemostasis` = "#1abc9c",
                                                   `Ion channel transport` = "#fc6404",
                                                   `Metabolism` = "#fcd444",
                                                   `Metabolism of proteins` = "#fcd444",
                                                   `Proteases/Hydrolases` = "#fcd444",
                                                   `SLC-mediated transmembrane transport` = "#fc6404",
                                                   `GPCRs` ="#fc6404",
                                                   `Signal Transduction` = "#fc6404",
                                                   `Signaling by GPCR` = "#fc6404",
                                                   `Transporters` = "#fc6404",
                                                   `Other / Unknown` = "gray",
                                                   `Viral` = "black",
                                                   `WNT.related` = "#fcd444")),
                           gp = gpar(col = "white"),simple_anno_size = unit(0.4, "cm"),
                           annotation_name_gp= gpar(fontsize = 0))


#setwd("C:/Users/pmunson/OneDrive - Nurix Therapeutics, Inc/Desktop/Autoantibodies/Graphs") 
#png("heatmap_baseline_melanoma_vs_hd.png",width=8,height=8,units="in",res=600)

Heatmap(heat_map_data, 
        col = col_fun,
       heatmap_legend_param = list(
          title = "REAP score",
          direction = "horizontal"),
      # annotation_legend_param = list(
      #   title = "Cancer",
      #   direction = "horizontal"),
        cluster_rows = FALSE,
        show_column_names = FALSE,
        cluster_columns = FALSE,
        clustering_distance_columns = "euclidean",
        row_names_gp = gpar(fontsize = 8), 
        column_names_gp = gpar(fontsize = 0),
        #rect_gp = gpar(col = "white", lwd = 0.01),
        top_annotation = col_heat_map_annotation,
        row_split = Annot$color,
        row_title_gp = gpar(fontsize = 0),
        column_title_gp = gpar(fontsize = 0),
        #annotation_row = annotation_row,
        column_title = "", 
        column_split = annotation_heatmap_df_select$Anti_CTLA_4,
        right_annotation = CATEGORIES,
        # column_split = column_split_info
      )

#dev.off()





}



frequently_hit_proteins <- raw_meta_df %>%
  filter(Cancer == "Melanoma") %>%
  filter(Visit_Name_time == "DAY_0") %>%
  select(Subject_ID, Protein, value) %>%
  filter(Protein %in% heatmap_proteins) %>%
  left_join(., annotation_df_v2) %>%
  filter(value >1) %>%
  group_by(Protein) %>%
  count() %>%
  arrange(desc(n)) %>%
  filter(n >= 3) %>% pull(Protein)

#### Creating lists of proteins 
### per each patient with a value above a given threshold
raw_meta_df %>%
  filter(Visit_Name_time == "DAY_0") %>%
  group_by(Subject_ID) %>%
  filter(value > 1) %>%
  group_split()


raw_meta_df %>%
  filter(Cancer == "Melanoma") %>%
  filter(Visit_Name_time == "DAY_89/101") %>%
  mutate(Responder = case_when(
    Best_Clinical_outcome %in% c("PR", "SD", "NED1") ~ "R",
    Best_Clinical_outcome %in% c("NED2", "PD") ~ "NR",
    TRUE ~ NA_character_  # Handle other cases if any
  )) %>%
  filter(Protein %in% frequently_hit_proteins) %>%
  ggplot(., aes(x = Responder, y = value)) +
  geom_point() +
  stat_compare_means(method = "t.test", label = "p.format",
                     label.y = 6,
                     size= 6,
                     comparisons = list(c("R", "NR"))) + 
  facet_wrap(~Protein)


raw_meta_df %>%
  #filter(Visit_Name_time == "DAY_0") %>%
  mutate(Responder = case_when(
    Best_Clinical_outcome %in% c("PR", "SD", "NED1") ~ "R",
    Best_Clinical_outcome %in% c("NED2", "PD") ~ "NR",
    TRUE ~ NA_character_  # Handle other cases if any
  )) %>%
  filter(Protein == "IL15RA") %>%
  ggplot(., aes(x = Visit_Name_time, y = value, group = Subject_ID)) +
  geom_point() +
  geom_line() +
  facet_wrap(~Responder)



### High concentration hits figure
{
high_concentration_ab_hits_graph <- raw_meta_df %>%
  filter(Visit_Name_time == "DAY_0") %>%
  #filter(reactome_group != "Other / Unknown") %>%
  #filter(!grepl("_Epitope", Protein)) %>%
  group_by(Subject_ID) %>%
  summarise(hits = sum(value > 7))%>%
  ungroup() %>%
  #pivot_wider(., names_from = reactome_group, values_from = n, values_fill = 0) %>%
  #pivot_longer(cols = -c(Subject_ID), names_to = "Annotation", values_to = "n") %>%
  left_join(., metadata_df_v2)  %>%
  ggplot(., aes(x = Cancer, y = hits, color = Cancer)) +
  #ggplot(., aes(x = n, y = Cancer, color = Cancer)) +
  #geom_hline(yintercept = 0, color = "white")+
  #geom_boxplot() +
  geom_violin(linewidth = 1, color = "black") +
  geom_jitter(width = 0.15, height = 0, size =4) +
  #coord_flip() +
  #facet_wrap(~Annotation) +
  ylab("No. of Reactivities (score >= 7)") +
  xlab("") +
  scale_color_manual(values = c('#ADD8E6','#FFB6C1')) +
  stat_compare_means(method = "t.test", label = "p.format",
                     label.y = 6,
                     size= 6,
                     comparisons = list(c("HD", "Melanoma"))) +
  theme_prism() + 
  theme(legend.position = "none") 

setwd("C:/Users/pmunson/OneDrive - Nurix Therapeutics, Inc/Desktop/Autoantibodies/Graphs") 
#ggsave(filename = "high_concentration_ab_hits_graph.png", plot = high_concentration_ab_hits_graph, width = 4, height = 6, dpi = 600)
  }

### Pathway hits graph FIG1C
{
  
  heatmap_proteins_pathway <-volcano_proteins %>%
    left_join(., annotation_df_v2) %>%
    #filter(reactome_group != "Other / Unknown") %>%
    mutate(annotation_and_reactome = ifelse(is.na(annotation_group), reactome_group, annotation_group)) %>%
    select(Protein, annotation_and_reactome) %>%
    #filter(annotation_and_reactome != "Other / Unknown") %>%
    mutate(annotation_and_reactome = factor(annotation_and_reactome, levels = c(pathway_levels))) %>%
    arrange(annotation_and_reactome) %>%
    #filter(Protein != "229E-RBD") %>%
    pull(Protein)
  
  pathway_levels <- volcano_proteins %>%
    left_join(., annotation_df_v2) %>%
    #filter(reactome_group != "Other / Unknown") %>%
    mutate(annotation_and_reactome = ifelse(is.na(annotation_group), reactome_group, annotation_group)) %>%
    select(Protein, annotation_and_reactome) %>%
    #filter(annotation_and_reactome != "Other / Unknown") %>%
    arrange(annotation_and_reactome) %>%
    select(annotation_and_reactome) %>% mutate(value = 0) %>%
    unique() %>%
    pivot_wider(., names_from = annotation_and_reactome, values_from = value) %>%
    ## define desired order here....
    select(c(Checkpoints.Drugs, Chemokines, Cytokines,"Interferons", "Innate Immune System", 
             "Signal Transduction", "Signaling by GPCR", "SLC-mediated transmembrane transport",
             "Transporters", "Ion channel transport", "GPCRs", "Metabolism", "Metabolism of proteins", "Proteases/Hydrolases",
             "WNT.related", everything())) %>%
    colnames()
  
pathway_hits_graph <- raw_meta_df %>%
  filter(Protein %in% heatmap_proteins_pathway) %>%
  #filter(reactome_group != "Other / Unknown") %>%
  mutate(annotation_and_reactome = ifelse(is.na(annotation_group), reactome_group, annotation_group)) %>%
  mutate(annotation_and_reactome = factor(annotation_and_reactome, levels = c(pathway_levels))) %>%
  mutate(annotation_and_reactome = case_when(
    annotation_and_reactome == "Checkpoints.Drugs" ~ "Checkpoint Drugs",
    annotation_and_reactome %in% c("Chemokines", "Cytokines", "Interferons") ~ "Chemokines/Cytokines/Interferons",
    annotation_and_reactome %in% c("Innate Immune System") ~ "Innate Immunity",
    annotation_and_reactome %in% c("Developmental Biology", "ECM organization", "Extracellular matrix organization", "Growth.Factors", "Hemostasis") ~ "Growth Factors/ECM/Hemostasis",
    annotation_and_reactome %in% c("Ion channel transport",  "SLC-mediated transmembrane transport", "GPCRs", "Signal Transduction", "Signaling by GPCR", "Transporters") ~ "Signal Transduction",
    annotation_and_reactome %in% c("Metabolism", "Metabolism of proteins", "Proteases/Hydrolases", "WNT.related") ~ "Metabolism/WNT",
    annotation_and_reactome == "Viral" ~ "Viral",
    annotation_and_reactome %in% c("Other / Unknown") ~ "Other / Unknown",
    TRUE ~ NA_character_ # default case
  )) %>%
  filter(Visit_Name_time == "DAY_0") %>%
  filter(!grepl("_Epitope", Protein)) %>%
  group_by(Subject_ID, annotation_and_reactome) %>%
  summarise(hits = sum(value > 1))%>%
  ungroup()  %>%
  pivot_wider(., names_from = annotation_and_reactome, values_from = hits, values_fill = 0) %>%
  pivot_longer(cols = -c(Subject_ID), names_to = "annotation_and_reactome", values_to = "n") %>%
  left_join(., metadata_df_v2) %>%
  #mutate(annotation_and_reactome = factor(annotation_and_reactome, levels = c(pathway_levels))) %>%
  ggplot(., aes(x = Cancer, y = n, color = Cancer)) +
  #ggplot(., aes(x = n, y = Cancer, color = Cancer)) +
  geom_hline(yintercept = 8, color = "white")+
  #geom_boxplot() +
  geom_violin(linewidth = 1, color = "black") +
  geom_jitter(width = 0.15, height = 0, size =4) +
  #coord_flip() +
  facet_wrap(~annotation_and_reactome, nrow = 4, scales = "free") +
  ylab("No. of Reactivities (score >= 1)") +
  xlab("") +
  scale_color_manual(values = c('#ADD8E6','#FFB6C1')) +
  stat_compare_means(method = "t.test", label = "p.format",
                     label.y = 6,
                     size= 6,
                     comparisons = list(c("HD", "Melanoma"))) +
  theme_prism() + 
  theme(legend.position = "none")

setwd("C:/Users/pmunson/OneDrive - Nurix Therapeutics, Inc/Desktop/Autoantibodies/Graphs") 
#ggsave(filename = "pathway_hits_graph.png", plot = pathway_hits_graph, width = 6, height = 16, dpi = 600)
}


##### Checkpoint AND OTHER INDIVIDUAL PROTEINS OF INTEREST, FIG1D
{
checkpoint_blockade_graph <- raw_meta_df %>%
  mutate(Anti_PD_1 = factor(Anti_PD_1, levels = c("HD", "0", "1"))) %>%
  #filter(Cancer == "Melanoma") %>%
  filter(Visit_Name_time == "DAY_0") %>% 
  filter(Protein == "PDCD1") %>%
  #filter(value > 0.5) %>%
  mutate(Responder = case_when(
    Best_Clinical_outcome %in% c("PR", "SD", "NED1") ~ "R",
    Best_Clinical_outcome %in% c("NED2", "PD") ~ "NR",
    TRUE ~ NA_character_  # Handle other cases if any
  )) %>%
  #filter(value > 0.5) %>%
  ggplot(., aes(x = Anti_PD_1, y = value, color = Anti_PD_1)) +
  geom_hline(yintercept = 7, color = "white") +
  geom_jitter(height =0, size = 4, width = 0.15, alpha = 0.75) +
  ylab("Anti-PD-1 (REAP Score)") + 
  xlab("") +
  theme_prism() +
  theme(legend.position = "none") +
  scale_color_manual(values = c('#ADD8E6','#FFB6C1', "#FFB6C1")) +
  
  raw_meta_df %>%
  mutate(Anti_CTLA_4 = factor(Anti_CTLA_4, levels = c("HD", "0", "1"))) %>%
  filter(Visit_Name_time == "DAY_0") %>% 
  filter(Protein == "CTLA4") %>%
  #filter(value > 0.5) %>%
  mutate(Responder = case_when(
    Best_Clinical_outcome %in% c("PR", "SD", "NED1") ~ "R",
    Best_Clinical_outcome %in% c("NED2", "PD") ~ "NR",
    TRUE ~ NA_character_  # Handle other cases if any
  )) %>%
  #filter(value > 0.5) %>%
  ggplot(., aes(x = Anti_CTLA_4, y = value, color = Anti_CTLA_4)) +
  geom_hline(yintercept = 7, color = "white") +
  geom_jitter(height =0, size = 4, width = 0.15, alpha = 0.75) +
  ylab("Anti-CTLA-4 (REAP Score)") + 
  xlab("") +
  theme_prism() +
  theme(legend.position = "none") +
  scale_color_manual(values = c('#ADD8E6','#FFB6C1', "#FFB6C1")) 
  

  setwd("C:/Users/pmunson/OneDrive - Nurix Therapeutics, Inc/Desktop/Autoantibodies/Graphs") 
  #ggsave(filename = "checkpoint_blockade_graph.png", plot = checkpoint_blockade_graph, width = 5, height = 4, dpi = 600)
  
  

  
individual_protein_hd_vs_melanoma_graph <- raw_meta_df %>%
    mutate(Anti_PD_1 = factor(Anti_PD_1, levels = c("HD", "0", "1"))) %>%
    #filter(Cancer == "Melanoma") %>%
    filter(Visit_Name_time == "DAY_0") %>% 
    filter(Protein == "IL15RA") %>%
    #filter(value > 0.5) %>%
    mutate(Responder = case_when(
      Best_Clinical_outcome %in% c("PR", "SD", "NED1") ~ "R",
      Best_Clinical_outcome %in% c("NED2", "PD") ~ "NR",
      TRUE ~ NA_character_  # Handle other cases if any
    )) %>%
    #filter(value > 0.5) %>%
    ggplot(., aes(x = Cancer, y = value, color = Cancer)) +
    geom_hline(yintercept = 6, color = "white") +
    geom_hline(yintercept = 1, color = "gray", linetype = "dashed", size = 1) +
    geom_jitter(height =0, size = 4, width = 0.15, alpha = 0.75) +
    ylab("") + 
    xlab("") +
    theme_prism() +
    ggtitle("IL15RA") +
    theme(legend.position = "none") +
    scale_color_manual(values = c('#ADD8E6','#FFB6C1', "red")) +
    stat_compare_means(method = "t.test", label = "p.format",
                       label.y = 5,
                       size= 6,
                       comparisons = list(c("HD", "Melanoma"))) +
    
    raw_meta_df %>%
    mutate(Anti_PD_1 = factor(Anti_PD_1, levels = c("HD", "0", "1"))) %>%
    #filter(Cancer == "Melanoma") %>%
    filter(Visit_Name_time == "DAY_0") %>% 
    filter(Protein == "DEFB118-N") %>%
    #filter(value > 0.5) %>%
    mutate(Responder = case_when(
      Best_Clinical_outcome %in% c("PR", "SD", "NED1") ~ "R",
      Best_Clinical_outcome %in% c("NED2", "PD") ~ "NR",
      TRUE ~ NA_character_  # Handle other cases if any
    )) %>%
    #filter(value > 0.5) %>%
    ggplot(., aes(x = Cancer, y = value, color = Cancer)) +
    geom_hline(yintercept = 6, color = "white") +
    geom_hline(yintercept = 1, color = "gray", linetype = "dashed", size = 1) +
    geom_jitter(height =0, size = 4, width = 0.15, alpha = 0.75) +
    ylab("") + 
    xlab("") +
    theme_prism() +
    ggtitle("DEFB118-N") +
    theme(legend.position = "none") +
    scale_color_manual(values = c('#ADD8E6','#FFB6C1', "red")) +
    stat_compare_means(method = "t.test", label = "p.format",
                       label.y = 5,
                       size= 6,
                       comparisons = list(c("HD", "Melanoma"))) +
  
  raw_meta_df %>%
    mutate(Anti_PD_1 = factor(Anti_PD_1, levels = c("HD", "0", "1"))) %>%
    #filter(Cancer == "Melanoma") %>%
    filter(Visit_Name_time == "DAY_0") %>% 
    filter(Protein == "APOH") %>%
    #filter(value > 0.5) %>%
    mutate(Responder = case_when(
      Best_Clinical_outcome %in% c("PR", "SD", "NED1") ~ "R",
      Best_Clinical_outcome %in% c("NED2", "PD") ~ "NR",
      TRUE ~ NA_character_  # Handle other cases if any
    )) %>%
    #filter(value > 0.5) %>%
    ggplot(., aes(x = Cancer, y = value, color = Cancer)) +
    geom_hline(yintercept = 6, color = "white") +
    geom_hline(yintercept = 1, color = "gray", linetype = "dashed", size = 1) +
    geom_jitter(height =0, size = 4, width = 0.15, alpha = 0.75) +
    ylab("") + 
    xlab("") +
    theme_prism() +
    ggtitle("APOH") +
    theme(legend.position = "none") +
    scale_color_manual(values = c('#ADD8E6','#FFB6C1', "red")) +
    stat_compare_means(method = "t.test", label = "p.format",
                       label.y = 5,
                       size= 6,
                       comparisons = list(c("HD", "Melanoma"))) +
  
  raw_meta_df %>%
    mutate(Anti_PD_1 = factor(Anti_PD_1, levels = c("HD", "0", "1"))) %>%
    #filter(Cancer == "Melanoma") %>%
    filter(Visit_Name_time == "DAY_0") %>% 
    filter(Protein == "FGFBP3") %>%
    #filter(value > 0.5) %>%
    mutate(Responder = case_when(
      Best_Clinical_outcome %in% c("PR", "SD", "NED1") ~ "R",
      Best_Clinical_outcome %in% c("NED2", "PD") ~ "NR",
      TRUE ~ NA_character_  # Handle other cases if any
    )) %>%
    #filter(value > 0.5) %>%
    ggplot(., aes(x = Cancer, y = value, color = Cancer)) +
    geom_hline(yintercept = 6, color = "white") +
    geom_hline(yintercept = 1, color = "gray", linetype = "dashed", size = 1) +
    geom_jitter(height =0, size = 4, width = 0.15, alpha = 0.75) +
    ylab("") + 
    xlab("") +
    theme_prism() +
    ggtitle("FGFBP3") +
    theme(legend.position = "none") +
    scale_color_manual(values = c('#ADD8E6','#FFB6C1', "red")) +
    stat_compare_means(method = "t.test", label = "p.format",
                       label.y = 5,
                       size= 6,
                       comparisons = list(c("HD", "Melanoma"))) +
  
    raw_meta_df %>%
    mutate(Anti_PD_1 = factor(Anti_PD_1, levels = c("HD", "0", "1"))) %>%
    #filter(Cancer == "Melanoma") %>%
    filter(Visit_Name_time == "DAY_0") %>% 
    filter(Protein == "ASIP") %>%
    #filter(value > 0.5) %>%
    mutate(Responder = case_when(
      Best_Clinical_outcome %in% c("PR", "SD", "NED1") ~ "R",
      Best_Clinical_outcome %in% c("NED2", "PD") ~ "NR",
      TRUE ~ NA_character_  # Handle other cases if any
    )) %>%
    #filter(value > 0.5) %>%
    ggplot(., aes(x = Cancer, y = value, color = Cancer)) +
    geom_hline(yintercept = 6, color = "white") +
    geom_hline(yintercept = 1, color = "gray", linetype = "dashed", size = 1) +
    geom_jitter(height =0, size = 4, width = 0.15, alpha = 0.75) +
    ylab("") + 
    xlab("") +
    theme_prism() +
    ggtitle("ASIP") +
    theme(legend.position = "none") +
    scale_color_manual(values = c('#ADD8E6','#FFB6C1', "red")) +
    stat_compare_means(method = "t.test", label = "p.format",
                       label.y = 5,
                       size= 6,
                       comparisons = list(c("HD", "Melanoma"))) +
  
  raw_meta_df %>%
    mutate(Anti_PD_1 = factor(Anti_PD_1, levels = c("HD", "0", "1"))) %>%
    #filter(Cancer == "Melanoma") %>%
    filter(Visit_Name_time == "DAY_0") %>% 
    filter(Protein == "APOC3") %>%
    #filter(value > 0.5) %>%
    mutate(Responder = case_when(
      Best_Clinical_outcome %in% c("PR", "SD", "NED1") ~ "R",
      Best_Clinical_outcome %in% c("NED2", "PD") ~ "NR",
      TRUE ~ NA_character_  # Handle other cases if any
    )) %>%
    #filter(value > 0.5) %>%
    ggplot(., aes(x = Cancer, y = value, color = Cancer)) +
    geom_hline(yintercept = 6, color = "white") +
    geom_hline(yintercept = 1, color = "gray", linetype = "dashed", size = 1) +
    geom_jitter(height =0, size = 4, width = 0.15, alpha = 0.75) +
    ylab("") + 
    xlab("") +
    theme_prism() +
    ggtitle("APOC3") +
    theme(legend.position = "none") +
    scale_color_manual(values = c('#ADD8E6','#FFB6C1', "red")) +
    stat_compare_means(method = "t.test", label = "p.format",
                       label.y = 5,
                       size= 6,
                       comparisons = list(c("HD", "Melanoma"))) 
  

  
setwd("C:/Users/pmunson/OneDrive - Nurix Therapeutics, Inc/Desktop/Autoantibodies/Graphs") 
#ggsave(filename = "individual_protein_hd_vs_melanoma_graph.png", plot = individual_protein_hd_vs_melanoma_graph, width = 12, height = 8, dpi = 600)
raw_meta_df$Overall_Survival


 

  
  
  
}

### Trying to find key autoantibodies that distinguish
## HD and Melanoma

{
raw_meta_df %>%
  filter(Protein %in% heatmap_proteins) %>%
  #filter(reactome_group != "Other / Unknown") %>%
  mutate(annotation_and_reactome = ifelse(is.na(annotation_group), reactome_group, annotation_group)) %>%
  mutate(annotation_and_reactome = factor(annotation_and_reactome, levels = c(pathway_levels))) %>%
  mutate(annotation_and_reactome = case_when(
    annotation_and_reactome == "Checkpoints.Drugs" ~ "Checkpoint Drugs",
    annotation_and_reactome %in% c("Chemokines", "Cytokines", "Interferons") ~ "Chemokines/Cytokines/Interferons",
    annotation_and_reactome %in% c("Innate Immune System") ~ "Innate Immunity",
    annotation_and_reactome %in% c("Developmental Biology", "ECM organization", "Extracellular matrix organization", "Growth.Factors", "Hemostasis") ~ "Growth Factors/ECM/Hemostasis",
    annotation_and_reactome %in% c("Ion channel transport",  "SLC-mediated transmembrane transport", "GPCRs", "Signal Transduction", "Signaling by GPCR", "Transporters") ~ "Signal Transduction",
    annotation_and_reactome %in% c("Metabolism", "Metabolism of proteins", "Proteases/Hydrolases", "WNT.related") ~ "Metabolism/WNT",
    annotation_and_reactome == "Viral" ~ "Viral",
    annotation_and_reactome %in% c("Other / Unknown") ~ "Other / Unknown",
    TRUE ~ NA_character_ # default case
  )) %>%
  filter(Visit_Name_time == "DAY_0") %>%
  filter(!grepl("_Epitope", Protein)) %>%
  left_join(., metadata_df_v2) %>%
  filter(annotation_and_reactome == "Chemokines/Cytokines/Interferons") %>%
  #mutate(annotation_and_reactome = factor(annotation_and_reactome, levels = c(pathway_levels))) %>%
  ggplot(., aes(x = Cancer, y = value, color = Cancer)) +
  #ggplot(., aes(x = n, y = Cancer, color = Cancer)) +
  #geom_hline(yintercept = 8, color = "white")+
  #geom_boxplot() +
  geom_violin(linewidth = 1, color = "black") +
  geom_jitter(width = 0.15, height = 0, size =4) +
  facet_wrap(~Protein, scales = "free") +
  ylab("REAP Score") +
  xlab("") +
  scale_color_manual(values = c('#ADD8E6','#FFB6C1')) +
  stat_compare_means(method = "t.test", label = "p.format",
                     label.y = 0,
                     size= 6,
                     comparisons = list(c("HD", "Melanoma"))) +
  theme_prism() + 
  theme(legend.position = "none")
}


### I want to do the volcano looking at responders vs. non-responders
## Responder # PR, SD, NED1
## Non-responder # PD, NED2


## Does vaccination change responses to proteins
### I think I should for all proteins per patient
## look at fold change and magnitude change for each protein




### Responder vs. Non resonder volcano plot
{
  
  
  raw_statistics_responder_df <- raw_meta_df %>%
    filter(Cancer == "Melanoma") %>%
    mutate(Responder = case_when(
      Best_Clinical_outcome %in% c("PR", "SD", "NED1") ~ "R",
      Best_Clinical_outcome %in% c("NED2", "PD") ~ "NR",
      TRUE ~ NA_character_  # Handle other cases if any
    )) %>%
    filter(Visit_Name_time == "DAY_0") %>%
    select(Responder, Protein, value) %>%
    group_by(Protein) %>%
    #### Calculating statistics
    summarise(
      lfc = log2(mean(value[Responder == "R"]) / mean(value[Responder == "NR"])),
      median_diff = median(value[Responder == "R"])-median(value[Responder == "NR"]),
      mean_diff = mean(value[Responder == "R"])-mean(value[Responder == "NR"]),
      max_diff = max(value[Responder == "R"])-max(value[Responder == "NR"]),
      max_value = max(value),
      p_value = t.test(value[Responder == "R"], value[Responder == "NR"], na.action = na.omit)$p.value
    ) 
  
  
  raw_statistics_responder_df %>%
    ## removing the Epitope
    filter(!grepl("_Epitope", Protein)) %>%
    left_join(., annotation_df_v2) %>%
    #filter(lfc != Inf & lfc != -Inf) %>%
    ggplot(., aes(x = lfc, y = -log(p_value), color = max_diff))+
    geom_point() +
    theme_prism() +
    #xlim(-20,20) +
    #ylim(0,10.9) +
    ylab("-log(p value)") +
    xlab("log 2 fold change") +
    ggtitle("Baseline") +
    geom_text_repel(data = subset(raw_statistics_responder_df %>%
                                    filter(!grepl("_Epitope", Protein)), 
                                  p_value < 0.25 & 
                                    (max_diff > 1 | max_diff < -1)), 
                    aes(label = Protein),
                    size = 3,
                    direction = "both",
                    max.overlaps = 20,min.segment.length = 0.01, box.padding = 0.25, color = "black") +
    scale_color_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0) +
    geom_text(x = Inf, y = Inf, label = "Up in Responders", hjust = 1.5, vjust = 1, color = "#3cb371", size =5) +
    geom_text(x = -Inf, y = Inf, label = "Up in Non-Responders", hjust = -0.75, vjust = 1, color = "orange", size = 5) +
    theme(legend.position = "none") 
  
  
  

  
  }


### Vaccine response volcano plot
{
  raw_statistics_Vaccine_Response_df <- raw_meta_df %>%
    filter(Cancer == "Yes") %>%
    filter(Visit_Name_time == "DAY_0") %>%
    filter(!is.na(Vaccine_Response)) %>%
    filter(Vaccine_Response != "NA") %>%
    select(Vaccine_Response, Protein, value) %>%
    group_by(Protein) %>%
    #### Calculating statistics
    summarise(
      lfc = log2(mean(value[Vaccine_Response == "Yes"]) / mean(value[Vaccine_Response == "No"])),
      median_diff = median(value[Vaccine_Response == "Yes"])-median(value[Vaccine_Response == "No"]),
      mean_diff = mean(value[Vaccine_Response == "Yes"])-mean(value[Vaccine_Response == "No"]),
      max_diff = max(value[Vaccine_Response == "Yes"])-max(value[Vaccine_Response == "No"]),
      max_value = max(value),
      p_value = t.test(value[Vaccine_Response == "Yes"], value[Vaccine_Response == "No"], na.action = na.omit)$p.value
    ) 
  
  
  raw_statistics_Vaccine_Response_df %>%
    ## removing the Epitope
    filter(!grepl("_Epitope", Protein)) %>%
    left_join(., annotation_df_v2) %>%
    #filter(lfc != Inf & lfc != -Inf) %>%
    ggplot(., aes(x = lfc, y = -log(p_value), color = max_diff))+
    geom_point() +
    theme_prism() +
    #xlim(-20,20) +
    #ylim(0,10.9) +
    ylab("-log(p value)") +
    xlab("log 2 fold change") +
    ggtitle("Baseline") +
    geom_text_repel(data = subset(raw_statistics_Vaccine_Response_df %>%
                                    filter(!grepl("_Epitope", Protein)), 
                                  p_value < 0.25 & 
                                    (max_diff > 1 | max_diff < -1)), 
                    aes(label = Protein),
                    size = 3,
                    direction = "both",
                    max.overlaps = 20,min.segment.length = 0.01, box.padding = 0.25, color = "black") +
    scale_color_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0) +
    geom_text(x = Inf, y = Inf, label = "Up in Vaccine_Responses", hjust = 1.5, vjust = 1, color = "green", size =5) +
    geom_text(x = -Inf, y = Inf, label = "Up in Non-Vaccine_Responses", hjust = -0.75, vjust = 1, color = "orange", size = 5) +
    theme(legend.position = "none")
  
  
  raw_meta_df %>%
    filter(!is.na(Vaccine_Response)) %>%
    filter(Vaccine_Response != "NA") %>%
    filter(Visit_Name_time == "DAY_0") %>% 
    filter(Protein == "IFNA8") %>%
    ggplot(., aes(x = Vaccine_Response, y = value, color = Best_Clinical_outcome)) +
    geom_jitter(width = 0.15)
  
}

### Random forest
{
# Install and load the randomForest package
install.packages("randomForest")
library(randomForest)
install.packages("Rtsne")
library(Rtsne)

### making this for my data
mach_learn_df <- raw_meta_df %>%
  filter(Cancer == "Melanoma") %>%
  filter(Protein %in% heatmap_proteins) %>%
  mutate(Protein = gsub("-", "_", Protein)) %>%
  mutate(Protein = gsub("_", "", Protein)) %>%
  mutate(Protein = gsub("229ES1", "viral_229ES1", Protein)) %>%
  mutate(Protein = gsub("229ERBD", "viral_229ERBD", Protein)) %>%
  mutate(Protein = gsub("UNQ6190/PRO20217", "UNQ6190_PRO20217", Protein)) %>%
  filter(Visit_Name_time == "DAY_0") %>%
  mutate(Responder = case_when(
    Best_Clinical_outcome %in% c("PR", "SD", "NED1") ~ "R",
    Best_Clinical_outcome %in% c("NED2", "PD") ~ "NR",
    TRUE ~ NA_character_  # Handle other cases if any
  )) %>%
  select(Subject_ID, Responder, Protein, value) %>%
  pivot_wider(., names_from = Protein, values_from = value)  %>% 
  mutate(Responder = factor(Responder, levels = c("NR", "R"))) %>%
  select(Subject_ID, everything())


rf_baseline_model <-randomForest(Responder ~ ., data = random_forest_df)

data.frame(Gini=sort(importance(rf_baseline_model, type=2)[,], decreasing=T))

#### doing tSNE

tsne_coord <- Rtsne(mach_learn_df[,-c(1)], dims = 2, perplexity=10, verbose=TRUE, max_iter = 500)$Y

mach_learn_df %>%
  select(Subject_ID) %>%
  cbind(., tsne_coord) %>%
  left_join(., metadata_df_v2, by = "Subject_ID") %>%
  rename("tSNE1" = "1") %>%
  rename("tSNE2" = "2") %>%
  mutate(Responder = case_when(
    Best_Clinical_outcome %in% c("PR", "SD", "NED1") ~ "R",
    Best_Clinical_outcome %in% c("NED2", "PD") ~ "NR",
    TRUE ~ NA_character_  # Handle other cases if any
  )) %>%
  ggplot(., aes(x = tSNE1, y = tSNE2, color = Responder)) +
  geom_point(size = 4)


iris_features <- iris[, -5] 
# Compute t-SNE transformation
tsne_result <- Rtsne(iris_features)
Rtsne(random_forest_df[,-1])
}


### Figure 2A
### DC vaccination with TA has a limited impact on exoproteome-specific antibodies 

{
  Low_baseline_values <- raw_meta_df %>%
    filter(Cancer == "Melanoma") %>%
    filter(Visit_Name_time == "DAY_0") %>%
    filter(value < 1) %>%
    mutate(Subject_Protein = str_c(Subject_ID, Protein, sep = "_")) %>%
    pull(Subject_Protein)
  
  
  detectable_day43_values <- raw_meta_df %>%
    filter(Cancer == "Melanoma") %>%
    filter(Visit_Name_time == "DAY_43") %>%
    filter(value > 1) %>%
    mutate(Subject_Protein = str_c(Subject_ID, Protein, sep = "_")) %>%
    pull(Subject_Protein)
  
  detectable_day89_values <- raw_meta_df %>%
    filter(Cancer == "Melanoma") %>%
    filter(Visit_Name_time == "DAY_89/101") %>%
    filter(value > 1) %>%
    mutate(Subject_Protein = str_c(Subject_ID, Protein, sep = "_")) %>%
    pull(Subject_Protein)
  
  intersect(detectable_day43_values,detectable_day89_values) %>%
    intersect(.,Low_baseline_values)

  
 # raw_statistics_vaccine_postdose <- 
### We need to be able to distinguish between
### samples that were run, and samples that were not
samples_ran <- raw_meta_df %>%
  select(Subject_ID, Visit_Name_time) %>%
  mutate(Subject_Time = str_c(Subject_ID, Visit_Name_time, sep = "_")) %>% 
  pull(Subject_Time)

all_three_timepoints_subjects <- raw_meta_df %>%
  select(Subject_ID, Visit_Name_time) %>%
  unique() %>%
  group_by(Subject_ID) %>% count() %>% filter(n >= 3) %>% pull(Subject_ID)

raw_statistics_vaccine_postdose <-raw_meta_df %>%
  mutate(Subject_Time = str_c(Subject_ID, Visit_Name_time, sep = "_")) %>% 
 # filter(Subject_Time %in% samples_ran) %>%
  filter(!grepl("_Epitope", Protein)) %>%
  filter(Cancer == "Melanoma") %>%
  select( Visit_Name_time, Subject_ID, Protein, value) %>%
  mutate(Protein_timmepoint = str_c(Protein, Visit_Name_time, sep = "_")) %>%
  select(  Subject_ID, Protein, Visit_Name_time, value) %>%
  pivot_wider(., names_from = Visit_Name_time, values_from = value) %>%
  rename(day_0_value = "DAY_0") %>%
  rename(day_43_value = "DAY_43") %>%
  rename(day_89_value = "DAY_89/101") %>%
  mutate(day_43_minus_0 = day_43_value - day_0_value) %>%
  mutate(day_89_minus_0 = day_89_value - day_0_value) %>%
   mutate(
     Subject_ID = factor(Subject_ID, 
                         levels = unique(Subject_ID[order(as.numeric(gsub(".*_(\\d+)", "\\1", Subject_ID)))]))
   ) 
 
  

  
  raw_statistics_vaccine_postdose %>%
    mutate(custom_color = ifelse(day_43_minus_0 > 1 & day_89_minus_0 > 1, "red", 
                                 ifelse(day_43_minus_0 < -1 & day_89_minus_0 < -1, "blue", "gray"))) %>%
    #filter(abs(day_43_minus_0) > 0.1 & abs(day_89_minus_0) > 0.1) %>%
    #filter(Subject_ID == "PT_8") %>%
    ggplot(., aes(x = day_43_minus_0, y = day_89_minus_0, color = custom_color)) +
    geom_hline(yintercept = 1, color = "black", linetype = "dashed") +
    geom_vline(xintercept = 1, color = "black", linetype = "dashed") +
    geom_hline(yintercept = -6, color = "white") +
    geom_hline(yintercept = 7, color = "white") +
    geom_vline(xintercept = -5, color = "white") +
    geom_vline(xintercept = 5, color = "white") +
    geom_point(alpha = 0.75, size = 4) +
    ylab("Day 89 - Day 0 (REAP Score)") +
    xlab("Day 43 - Day 0 (REAP Score)") +
    geom_text_repel(data = raw_statistics_vaccine_postdose %>%
                      mutate(custom_color = ifelse(day_43_minus_0 > 1 & day_89_minus_0 > 1, "red", 
                                                   ifelse(day_43_minus_0 < -1 & day_89_minus_0 < -1, "blue", "gray"))) %>%
                      #filter(Subject_ID == "PT_8") %>%
                      filter(abs(day_43_minus_0) > 0.1 | abs(day_89_minus_0) > 0.1) %>%
                      filter((day_43_minus_0 > 1 & day_89_minus_0 > 1) | (day_43_minus_0 < -1 & day_89_minus_0 < -1)), 
                    aes(label = Protein), 
                    size = 3,
                    box.padding = 0.5, max.overlaps = 40) +
    #facet_wrap(~Subject_ID, scale = "free") +
    theme_prism() +
    scale_color_identity()
    theme(legend.position = "none")


#

differential_graph <- raw_statistics_vaccine_postdose %>% 
  filter(Subject_ID %in% all_three_timepoints_subjects) %>%
  mutate(custom_color = ifelse(day_43_minus_0 > 1, "red", 
                               ifelse(day_43_minus_0 < -1, "blue", "gray"))) %>%
  ggplot(., aes(x = Subject_ID, y = day_43_minus_0, color = custom_color)) +
  geom_hline(yintercept = 8, linetype = "dashed", linewidth = 1, color = "white") +
  geom_hline(yintercept = -8, linetype = "dashed", linewidth = 1, color = "white") +
  geom_jitter(width = 0.15,alpha = 0.5, height = 0, size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = -1, linetype = "dashed", linewidth = 1) +
  theme_prism() +
  xlab("") +
  ylab("Day 43 - Day 0 (REAP Score)") + 
  scale_color_identity() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))  +
  geom_text_repel(data = raw_statistics_vaccine_postdose %>% 
                    filter(Subject_ID %in% all_three_timepoints_subjects) %>%
                    mutate(custom_color = ifelse(day_43_minus_0 > 1, "red", 
                                                 ifelse(day_43_minus_0 < -1, "blue", "gray"))) %>%
                    filter(custom_color == "red" | custom_color == "blue"), 
                  aes(label = Protein), 
                  size = 3,min.segment.length = 0,
                  box.padding = 0.5, max.overlaps = 15) +

raw_statistics_vaccine_postdose %>% 
  filter(Subject_ID %in% all_three_timepoints_subjects) %>%
  mutate(custom_color = ifelse(day_89_minus_0 > 1, "red", 
                               ifelse(day_89_minus_0 < -1, "blue", "gray"))) %>%
  ggplot(., aes(x = Subject_ID, y = day_89_minus_0, color = custom_color)) +
  geom_hline(yintercept = 8, linetype = "dashed", linewidth = 1, color = "white") +
  geom_hline(yintercept = -8, linetype = "dashed", linewidth = 1, color = "white") +
  geom_jitter(width = 0.15,alpha = 0.5, height = 0, size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = -1, linetype = "dashed", linewidth = 1) +
  theme_prism() +
  xlab("") +
  ylab("Day 89 - Day 0 (REAP Score)") + 
  scale_color_identity() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  geom_text_repel(data = raw_statistics_vaccine_postdose %>% 
                    filter(Subject_ID %in% all_three_timepoints_subjects) %>%
                    mutate(custom_color = ifelse(day_89_minus_0 > 1, "red", 
                                                 ifelse(day_89_minus_0 < -1, "blue", "gray"))) %>%
                    filter(custom_color == "red" | custom_color == "blue"), 
                  aes(label = Protein), 
                  size = 3, min.segment.length = 0,
                  box.padding = 0.5, max.overlaps = 15) 

 
setwd("C:/Users/pmunson/OneDrive - Nurix Therapeutics, Inc/Desktop/Autoantibodies/Graphs") 
#ggsave(filename = "differential_graph.png", plot = differential_graph, width = 18, height = 6, dpi = 600)
}


### Venn diagram

Differential_ab_d43 <- raw_statistics_vaccine_postdose %>%
  filter(abs(day_43_minus_0) >1) %>%
  select(Subject_ID, Protein, day_43_minus_0) %>%
  mutate(Up_or_down_day_43 = ifelse(day_43_minus_0 >1, "Up", "Down")) %>%
  mutate(ID_Protein_direction = str_c(Subject_ID, Protein, Up_or_down_day_43, sep = "_")) %>%
  select(ID_Protein_direction)

Differential_ab_d89 <- raw_statistics_vaccine_postdose %>%
  filter(abs(day_89_minus_0) >1) %>%
  select(Subject_ID, Protein, day_89_minus_0) %>%
  mutate(Up_or_down_day_43 = ifelse(day_89_minus_0 >1, "Up", "Down")) %>%
  mutate(ID_Protein_direction = str_c(Subject_ID, Protein, Up_or_down_day_43, sep = "_")) %>%
  select(ID_Protein_direction)

ggvenn(data = list('D43' = Differential_ab_d43$ID_Protein_direction, 'D89' = Differential_ab_d89$ID_Protein_direction), 
       show_elements = F, 
       label_sep = "\n", 
       show_percentage = T,
       fill_color = c("black", "red"),
       text_size = 7.5,
       set_name_size = 7.5,
       fill_alpha =0.2)

intersect_result <- intersect(Differential_ab_d43,Differential_ab_d89) %>% data.frame()

glimpse(intersect_result)
intersect_result$ID_Protein_direction
#extracted_info <- 

consistently_altered_autoantibody <- intersect_result %>%
  ### This gives me what I originally asked for
  mutate(extracted_info = gsub(".*_([^_]+)_.*", "\\1", ID_Protein_direction)) %>%
  #mutate(extracted_info = gsub(".*_([^_]+_.*)", "\\1", ID_Protein_direction)) %>%
  group_by(extracted_info) %>%
  count() %>%
  arrange(desc(n)) %>% 
  pull(extracted_info)

### create a heatmap of proteins that are modulated in the same patient in both timepoints
### Figure 2 b/c
{
  Differential_ab_d43 <- raw_statistics_vaccine_postdose %>%
    filter(abs(day_43_minus_0) >1) %>%
    select(Subject_ID, Protein, day_43_minus_0) %>%
    mutate(Up_or_down_day_43 = ifelse(day_43_minus_0 >1, "Up", "Down")) %>%
    mutate(ID_Protein_direction = str_c(Subject_ID, Protein, Up_or_down_day_43, sep = "_")) %>%
    select(ID_Protein_direction)
  
  Differential_ab_d89 <- raw_statistics_vaccine_postdose %>%
    filter(abs(day_89_minus_0) >1) %>%
    select(Subject_ID, Protein, day_89_minus_0) %>%
    mutate(Up_or_down_day_43 = ifelse(day_89_minus_0 >1, "Up", "Down")) %>%
    mutate(ID_Protein_direction = str_c(Subject_ID, Protein, Up_or_down_day_43, sep = "_")) %>%
    select(ID_Protein_direction)
  
  ggvenn(data = list('D43' = Differential_ab_d43$ID_Protein_direction, 'D89' = Differential_ab_d89$ID_Protein_direction), 
         show_elements = F, 
         label_sep = "\n", 
         show_percentage = T,
         fill_color = c("black", "red"),
         text_size = 7.5,
         set_name_size = 7.5,
         fill_alpha =0.2)
  
  intersect_result <- intersect(Differential_ab_d43,Differential_ab_d89) %>% data.frame()
  
  glimpse(intersect_result)
  intersect_result$ID_Protein_direction
  #extracted_info <- 
  
  consistently_altered_autoantibody <- intersect_result %>%
    ### This gives me what I originally asked for
    mutate(extracted_info = gsub(".*_([^_]+)_.*", "\\1", ID_Protein_direction)) %>%
    #mutate(extracted_info = gsub(".*_([^_]+_.*)", "\\1", ID_Protein_direction)) %>%
    group_by(extracted_info) %>%
    count() %>%
    arrange(desc(n)) %>% 
    pull(extracted_info)
  
  #### Getting pathway information for these proteins
  
  pathway_levels <- consistently_altered_autoantibody %>%
    as.data.frame() %>% rename(Protein = ".") %>%
    left_join(., annotation_df_v2) %>%
    #filter(reactome_group != "Other / Unknown") %>%
    mutate(annotation_and_reactome = ifelse(is.na(annotation_group), reactome_group, annotation_group)) %>%
    select(Protein, annotation_and_reactome) %>%
    #filter(annotation_and_reactome != "Other / Unknown") %>%
    arrange(annotation_and_reactome) %>%
    select(annotation_and_reactome) %>% mutate(value = 0) %>%
    unique() %>%
    pivot_wider(., names_from = annotation_and_reactome, values_from = value) %>%
    ## define desired order here....
    select(c(Checkpoints.Drugs, "Checkpoints.Inhibitory",Chemokines, Cytokines,"Interferons", "Innate Immune System",
             "NKG2D.ligands","Adaptive Immune System",
             "Signal Transduction", "Signaling by GPCR", "Ion channel transport", "Sensory Receptors",
             "Metabolism", "Metabolism of proteins", "Proteases/Hydrolases", everything())) %>%
    colnames()
  
  heatmap_proteins <-consistently_altered_autoantibody %>%
    as.data.frame() %>% rename(Protein = ".") %>%
    left_join(., annotation_df_v2) %>%
    #filter(reactome_group != "Other / Unknown") %>%
    mutate(annotation_and_reactome = ifelse(is.na(annotation_group), reactome_group, annotation_group)) %>%
    select(Protein, annotation_and_reactome) %>%
    #filter(annotation_and_reactome != "Other / Unknown") %>%
    mutate(annotation_and_reactome = factor(annotation_and_reactome, levels = c(pathway_levels))) %>%
    arrange(annotation_and_reactome) %>%
    #filter(Protein != "229E-RBD") %>%
    pull(Protein)
  
  
  all_three_timepoints_subjects <- raw_meta_df %>%
    select(Subject_ID, Visit_Name_time) %>%
    unique() %>%
    group_by(Subject_ID) %>% count() %>% filter(n >= 3) %>% pull(Subject_ID)
  
  unique(raw_meta_df$Visit_Name_time)
  
  heat_map_data <- raw_statistics_vaccine_postdose %>%
    filter(Subject_ID %in% all_three_timepoints_subjects) %>%
    filter(Protein %in% consistently_altered_autoantibody) %>%
    #select(Subject_ID, Protein, day_43_minus_0) %>%
    #pivot_wider(., names_from = Protein, values_from = day_43_minus_0) %>%
    select(Subject_ID, Protein, day_89_minus_0) %>%
    pivot_wider(., names_from = Protein, values_from = day_89_minus_0) %>%
    select(-Subject_ID) %>% 
    #scale() %>%
    t()
  
  heat_map_data <-heat_map_data[heatmap_proteins,]
  
  colnames_heat_map <- raw_statistics_vaccine_postdose %>%
    filter(Subject_ID %in% all_three_timepoints_subjects) %>%
    filter(Protein %in% consistently_altered_autoantibody) %>%
    select(Subject_ID, Protein, day_43_minus_0) %>%
    pivot_wider(., names_from = Protein, values_from = day_43_minus_0) %>%
    mutate(Subject_ID = as.character(Subject_ID)) %>%
    select(Subject_ID) %>% pull()
  
  colnames(heat_map_data)  <- colnames_heat_map
  
  ### annotation data
  raw_meta_df$Prior_Checkpoint
  annotation_heatmap <-raw_meta_df%>%
    mutate(Responder = case_when(
      Best_Clinical_outcome %in% c("PR", "SD", "NED1") ~ "R",
      Best_Clinical_outcome %in% c("NED2", "PD") ~ "NR",
      TRUE ~ NA_character_  # Handle other cases if any
    )) %>%
    #select(Subject_ID, Best_Clinical_outcome, Overall_Survival,PFS, Metastasis_Site, Anti_CTLA_4, Anti_PD_1, Prior_IFN, Prior_IL2, Cancer) 
    select(Subject_ID,Responder,Prior_Checkpoint)
  
  # Reorder the annotation data to match the matrix column names
  #annotation_data <- 
  
  annotation_heatmap <- annotation_heatmap[match(colnames_heat_map, annotation_heatmap$Subject_ID), ]
  
  ## If this is zero, it means the 
  ## annotation rows are in the same order as the data columns 
  sum(annotation_heatmap$Subject_ID!= colnames_heat_map)
  
  ## checking that the heatmap and the annotation table are compatible in dimensions
  dim(heat_map_data)[2] == dim(annotation_heatmap)[1]
  
  ## I need to make the subject_visit a rowname
  annotation_heatmap_df <- as.data.frame(annotation_heatmap)
  rownames(annotation_heatmap_df) <-  annotation_heatmap_df[,1] 
  annotation_heatmap_df_select <- annotation_heatmap_df %>% select(-Subject_ID)
  
  ##This should be zero, means the rows of the annotation perfectly match
  # the columns of the data 
  sum(rownames(annotation_heatmap_df_select) != colnames(heat_map_data))
  
  #### Annotations
  
  ### Row annotations
  ### this is copied from the earlier chunk of code to ensure the same order as the data 
  color_levels <- c("#fc8c84", "#6454ac", "#5bc0de", "#1abc9c", "#fc6404", "#fcd444", "gray", "black")
  
  
  
  Annot <- rownames(heat_map_data) %>%
    data.frame() %>%
    rename(Protein = ".") %>%
    left_join(., annotation_df_v2) %>%
    #filter(reactome_group != "Other / Unknown") %>%
    mutate(annotation_and_reactome = ifelse(is.na(annotation_group), reactome_group, annotation_group)) %>%
    select(Protein, annotation_and_reactome) %>%
    #filter(annotation_and_reactome != "Other / Unknown") %>%
    mutate(annotation_and_reactome = factor(annotation_and_reactome, 
                                            levels = c(pathway_levels))) %>%
    mutate(colors = case_when(
      annotation_and_reactome %in% c("Checkpoints.Drugs","Checkpoints.Inhibitory") ~ "#fc8c84",
      annotation_and_reactome %in% c("Chemokines", "Cytokines", "Interferons") ~ "#6454ac",
      annotation_and_reactome %in% c("Innate Immune System","NKG2D.ligands","Adaptive Immune System") ~ "#5bc0de",
      annotation_and_reactome %in% c("Developmental Biology", "ECM organization", "Extracellular matrix organization", "Growth.Factors", "Hemostasis") ~ "#1abc9c",
      annotation_and_reactome %in% c("Ion channel transport",  "SLC-mediated transmembrane transport", "GPCRs", "Signal Transduction", "Signaling by GPCR", "Transporters",
                                     "Sensory Receptors") ~ "#fc6404",
      annotation_and_reactome %in% c("Metabolism", "Metabolism of proteins", "Proteases/Hydrolases", "WNT.related") ~ "#fcd444",
      annotation_and_reactome %in% "Viral" ~ "black",
      annotation_and_reactome %in% c("Other / Unknown") ~ "gray",
      TRUE ~ NA_character_ # default case
    )) %>% 
    mutate(colors = factor(colors, levels = color_levels)) 
  
  

  
  # Here is defining of Annotation fo CATEGORIES
  CATEGORIES = rowAnnotation(annotation_and_reactome = Annot$annotation_and_reactome,
                             col = list(annotation_and_reactome = c(`Checkpoints.Drugs` = "#fc8c84",
                                                                    `Checkpoints.Inhibitory` = "#fc8c84",
                                                                    `Chemokines` = "#6454ac",
                                                                    `Cytokines` = "#6454ac",
                                                                    `Interferons` = "#6454ac",
                                                                    `Innate Immune System` = "#5bc0de",
                                                                    `NKG2D.ligands` = "#5bc0de",
                                                                    `Adaptive Immune System` = "#5bc0de",
                                                                    `Developmental Biology` = "#1abc9c",
                                                                    `ECM organization` = "#1abc9c",
                                                                    `Extracellular matrix organization` = "#1abc9c",
                                                                    `Growth.Factors` = "#1abc9c",
                                                                    `Hemostasis` = "#1abc9c",
                                                                    `Ion channel transport` = "#fc6404",
                                                                    `Sensory Receptors` = "#fc6404",
                                                                    `Metabolism` = "#fcd444",
                                                                    `Metabolism of proteins` = "#fcd444",
                                                                    `Proteases/Hydrolases` = "#fcd444",
                                                                    `SLC-mediated transmembrane transport` = "#fc6404",
                                                                    `GPCRs` ="#fc6404",
                                                                    `Signal Transduction` = "#fc6404",
                                                                    `Signaling by GPCR` = "#fc6404",
                                                                    `Transporters` = "#fc6404",
                                                                    `Other / Unknown` = "gray",
                                                                    `Viral` = "black",
                                                                    `WNT.related` = "#fcd444")),
                             gp = gpar(col = "white"),simple_anno_size = unit(0.4, "cm"),
                             annotation_name_gp= gpar(fontsize = 0))
  
  #### annotation end
  
  colours <- list(
    Cancer = c('HD'= '#ADD8E6', 'Melanoma' = '#FFB6C1'),
    Anti_PD_1 = c('0' = 'black', '1' = '#fc8c84', 'HD' = 'white' ),
    Anti_CTLA_4 = c('0' = 'black', '1' = '#fc8c84', 'HD' = 'white'),
    Prior_IFN = c('0' = 'black', '1' = '#38B489', 'NA' = 'white'),
    Determinant_Spreading = c('No' = 'black', 'Yes' = '#fc8c84', 'NA' = 'white'),
    Responder = c('NR'= '#F6B32A', 'R' = '#3EB489')
  )
  
  col_heat_map_annotation <- HeatmapAnnotation(
    df = annotation_heatmap_df_select,
    col = colours, 
    annotation_height = 0.6,
    annotation_width = unit(1, 'cm'),
    gap = unit(1, 'mm'))
  
  
  ### color of the cells
  #col_fun = colorRamp2(c(-4,-1,0,1,4), c("#0041C2","#9eb9d4","black","#FFC1c", "red"))
  col_fun = colorRamp2(c(-1,0,1), c("#0041C2","black", "red"))
  
  ### Row annotations
  ### this is copied from the earlier chunk of code to ensure the same order as the data 
  color_levels <- c("#fc8c84", "#6454ac", "#5bc0de", "#1abc9c", "#fc6404", "#fcd444", "gray", "black")
  
  #setwd("C:/Users/pmunson/OneDrive - Nurix Therapeutics, Inc/Desktop/Autoantibodies/Graphs") 
  #png("day_89_differential_R_vs_NR.png",width=8,height=8,units="in",res=600)
  
  Heatmap(heat_map_data, 
          col = col_fun,
          heatmap_legend_param = list(
            title = "REAP score",
            direction = "horizontal"),
          # annotation_legend_param = list(
          #   title = "Cancer",
          #   direction = "horizontal"),
          cluster_rows = FALSE,
          show_column_names = FALSE,
          cluster_columns = FALSE,
          clustering_distance_columns = "manhattan",
          row_names_gp = gpar(fontsize = 8), 
          column_names_gp = gpar(fontsize = 0),
          rect_gp = gpar(col = "black", lwd = 0.01),
          top_annotation = col_heat_map_annotation,
          row_split = Annot$color,
          row_title_gp = gpar(fontsize = 0),
          column_title_gp = gpar(fontsize = 0),
          #annotation_row = annotation_row,
          column_title = "", 
          row_gap = unit(1, "mm"), 
          column_gap = unit(3, "mm"),
          column_split = annotation_heatmap_df_select$Responder,
          right_annotation = CATEGORIES,
          # column_split = column_split_info
  )
  
  #dev.off()
  
  
  ### Quantifying hits based on pathways/reactome
  
  
  
  
  
  responder_pathway_counts <- raw_statistics_vaccine_postdose %>%
    filter(Protein %in% consistently_altered_autoantibody)  %>%
    select(Subject_ID, Protein, day_89_minus_0) %>%
    left_join(., annotation_df_v2) %>%
    #filter(reactome_group != "Other / Unknown") %>%
    mutate(annotation_and_reactome = ifelse(is.na(annotation_group), reactome_group, annotation_group)) %>%
    mutate(annotation_and_reactome = factor(annotation_and_reactome, levels = c(pathway_levels))) %>%
    mutate(aggregated_paths = case_when(
      annotation_and_reactome %in% c("Checkpoints.Drugs","Checkpoints.Inhibitory") ~ "Checkpoint",
      annotation_and_reactome %in% c("Chemokines", "Cytokines", "Interferons") ~ "Chemo/Cyto/IFNs",
      annotation_and_reactome %in% c("Innate Immune System","NKG2D.ligands","Adaptive Immune System") ~ "Innate and Adaptive",
      annotation_and_reactome %in% c("Developmental Biology", "ECM organization", "Extracellular matrix organization", "Growth.Factors", "Hemostasis") ~ "Dev. Bio",
      annotation_and_reactome %in% c("Ion channel transport",  "SLC-mediated transmembrane transport", "GPCRs", "Signal Transduction", "Signaling by GPCR", "Transporters",
                                     "Sensory Receptors") ~ "Transporters",
      annotation_and_reactome %in% c("Metabolism", "Metabolism of proteins", "Proteases/Hydrolases", "WNT.related") ~ "Metabolism",
      annotation_and_reactome %in% "Viral" ~ "Viral",
      annotation_and_reactome %in% c("Other / Unknown") ~ "Other",
      TRUE ~ NA_character_ # default case
    )) %>%
    filter(!grepl("_Epitope", Protein)) %>%
    filter(!is.na(day_89_minus_0)) %>%
    left_join(., metadata_df_v2) %>%
    #filter(Subject_ID == "PT_21" & annotation_and_reactome == "Innate Immune System")
    mutate(Responder = case_when(
      Best_Clinical_outcome %in% c("PR", "SD", "NED1") ~ "R",
      Best_Clinical_outcome %in% c("NED2", "PD") ~ "NR",
      TRUE ~ NA_character_  # Handle other cases if any
    )) %>%
    group_by(Responder, annotation_and_reactome) %>%
    summarize(
      count_positive = sum(day_89_minus_0 > 1, na.rm = TRUE),
      count_negative = sum(day_89_minus_0 < 1, na.rm = TRUE)
    ) %>%
    ungroup()  %>%
    select(-(count_negative)) %>%
    pivot_wider(., names_from = annotation_and_reactome, values_from = count_positive, values_fill = 0) %>%
    pivot_longer(cols = -c(Responder), names_to = "annotation_and_reactome", values_to = "n") %>%
    filter(!is.na(Responder)) %>%
    #filter(!is.na(Determinant_Spreading)) %>%
    #filter(Determinant_Spreading != "NA") %>%
    #mutate(annotation_and_reactome = factor(annotation_and_reactome, levels = c(pathway_levels))) %>%
    filter(annotation_and_reactome != "Other / Unknown" &
             annotation_and_reactome != "Checkpoint.Drugs" &
             annotation_and_reactome != "Metabolism" &
             annotation_and_reactome != "NKG2D.ligands" &
             annotation_and_reactome != "Sensory Receptors" &
             annotation_and_reactome != "Signaling by GPCR") %>%
    mutate(n = ifelse(Responder == "NR", n*-1,n)) 
  
  graphed_levels <- responder_pathway_counts %>%
    filter(Responder == "R") %>%
    arrange(n) %>% pull(annotation_and_reactome)
  
pathway_differential_colomn_graph <- responder_pathway_counts %>%
    mutate(annotation_and_reactome = factor(annotation_and_reactome, levels = c(graphed_levels))) %>%
    ggplot(., aes(x = annotation_and_reactome, y = n, fill = Responder, color = "black")) +
    geom_hline(yintercept = 0, color = "white") +
    geom_hline(yintercept = 20, color = "white") +
    geom_hline(yintercept = -20, color = "white") +
    #ggplot(., aes(x = n, y = Cancer, color = Cancer)) +
    geom_hline(yintercept = 8, color = "white")+
    #geom_boxplot() +
    geom_col() +
    #coord_flip() +
    #facet_wrap(~Responder, scales = "free_x") +
    ylab("No. of Reactivies (REAP >0 over baseline) R vs. NR") +
    xlab("") +
    coord_flip() +
    scale_color_manual(values = c('black')) +
    scale_fill_manual(values = c('red', 'blue')) + # Specify fill colors
    theme_prism() + 
    theme(legend.position = "none")
  
  
  
  #setwd("C:/Users/pmunson/OneDrive - Nurix Therapeutics, Inc/Desktop/Autoantibodies/Graphs") 
  #ggsave(filename = "pathway_differential_colomn_graph.png", plot = pathway_differential_colomn_graph, width =8, height = 8, dpi = 600)
  }

### identify unique things of PRs

### volcano of responders and non-responders at baseline:

{
  
  {volcano_proteins_responders <- subset(raw_statistics_responders_day0_df %>%
                                           ## removing the Epitope
                                           filter(!grepl("_Epitope", Protein)), 
                                         #`p_value` < 0.25  & lfc >1 &
                                         ### this p value should detect if its in at least 2 subjects
                                         `p_value` < 0.4&
                                           lfc > 1 &
                                           abs(max_diff) > 1)
  ### trying to remove proteins that might be noisey
  
  
  detectable_in_healthy_donors <- raw_meta_df %>%
    filter(Cancer == "HD") %>%
    filter(value > 1) %>%
    pull(Protein) %>% unique()
  
  setdiff(volcano_proteins_responders %>% pull(Protein),detectable_in_healthy_donors)
  
  volcano_proteins_responders <-volcano_proteins_responders %>%
    filter(!Protein %in% detectable_in_healthy_donors)
  
  ### Manually assigning functions to potential proteins of interest
  unknown_proteins <- volcano_proteins_responders %>%
    left_join(., annotation_df_v2) %>%
    #filter(reactome_group != "Other / Unknown") %>%
    mutate(annotation_and_reactome = ifelse(is.na(annotation_group), reactome_group, annotation_group)) %>%
    select(Protein, annotation_and_reactome) %>%
    filter(annotation_and_reactome == "Other / Unknown") %>%
    pull(Protein)
  
  
  ### select desired order of proteins based on pathway similarites on the heat map
  
  
  pathway_levels <-  volcano_proteins_responders%>%
    left_join(., annotation_df_v2) %>%
    mutate(annotation_and_reactome = ifelse(is.na(annotation_group), reactome_group, annotation_group)) %>%
    select(Protein, annotation_and_reactome) %>%
    mutate(annotation_and_reactome = ifelse(Protein == "IL12-sc", "Cytokines", annotation_and_reactome)) %>%
    filter(annotation_and_reactome!= "Other / Unknown") %>%
    unique() %>%
    #filter(annotation_and_reactome != "Other / Unknown") %>%
    arrange(annotation_and_reactome) %>%
    select(annotation_and_reactome) %>%
    unique() %>%
    mutate(value = 0) %>%
    #print(n = Inf)
    pivot_wider(., names_from = annotation_and_reactome, values_from = value) %>%
    ## define desired order here....
    select(any_of(c("Checkpoints.Drugs", "Checkpoints.Inhibitory","Chemokines", "Cytokines","Interferons", "Innate Immune System",
                    "NKG2D.ligands","Adaptive Immune System", "Immune (other)",
                    everything())))%>%
    colnames() 
  
  heatmap_proteins <-volcano_proteins_responders %>%
    left_join(., annotation_df_v2) %>%
    mutate(annotation_and_reactome = ifelse(is.na(annotation_group), reactome_group, annotation_group)) %>%
    select(Protein, annotation_and_reactome) %>%
    mutate(annotation_and_reactome = ifelse(Protein == "IL12-sc", "Cytokines", annotation_and_reactome)) %>%
    filter(annotation_and_reactome != "Other / Unknown") %>% 
    #filter(annotation_and_reactome != "Other / Unknown") %>%
    mutate(annotation_and_reactome = factor(annotation_and_reactome, levels = c(pathway_levels))) %>%
    arrange(annotation_and_reactome) %>%
    filter(!is.na(annotation_and_reactome)) %>%
    #filter(Protein != "229E-RBD") %>%
    pull(Protein)
  
  heat_map_data <- raw_meta_df %>%
    filter(Visit_Name_time == "DAY_0") %>%
    select(Subject_ID, Protein, value) %>%
    filter(Protein %in% heatmap_proteins) %>%
    pivot_wider(., names_from = Protein, values_from = value) %>%
    select(-Subject_ID) %>%
    #scale() %>%
    t() 
  
  heat_map_data <-heat_map_data[heatmap_proteins,]
  
  colnames_heat_map <- raw_meta_df %>%
    filter(Visit_Name_time == "DAY_0") %>%
    select(Subject_ID, Protein, value) %>%
    filter(Protein %in% heatmap_proteins) %>%
    pivot_wider(., names_from = Protein, values_from = value) %>%
    select(Subject_ID) %>% pull()
  
  colnames(heat_map_data) <- colnames_heat_map
  
  ### Annotation data
  
  unique(raw_meta_df$Best_Clinical_outcome)
  
  annotation_heatmap <-raw_meta_df%>%
    mutate(Responder = case_when(
      Best_Clinical_outcome %in% c("PR", "SD", "NED1") ~ "R",
      Best_Clinical_outcome %in% c("NED2", "PD") ~ "NR",
      Best_Clinical_outcome %in% c("HD") ~ "HD",
      TRUE ~ NA_character_  # Handle other cases if any
    )) %>%
    mutate(Responder = factor(Responder, levels = c("HD", "R", "NR"))) %>%
    #select(Subject_ID, Best_Clinical_outcome, Overall_Survival,PFS, Metastasis_Site, Anti_CTLA_4, Anti_PD_1, Prior_IFN, Prior_IL2, Cancer) 
    select(Subject_ID,Cancer,  Anti_CTLA_4,Anti_PD_1, Responder ) %>%
    mutate(Anti_CTLA_4 = factor(Anti_CTLA_4, levels = c("HD", "0", "1"))) %>%
    select(Subject_ID, Responder)
  
  
  # Reorder the annotation data to match the matrix column names
  #annotation_data <- 
  
  annotation_heatmap <- annotation_heatmap[match(colnames_heat_map, annotation_heatmap$Subject_ID), ]
  
  
  ## If this is zero, it means the 
  ## annotation rows are in the same order as the data columns 
  sum(annotation_heatmap$Subject_ID!= colnames_heat_map)
  
  ## checking that the heatmap and the annotation table are compatible in dimensions
  dim(heat_map_data)[2] == dim(annotation_heatmap)[1]
  
  ## I need to make the subject_visit a rowname
  annotation_heatmap_df <- as.data.frame(annotation_heatmap)
  rownames(annotation_heatmap_df) <-  annotation_heatmap_df[,1] 
  annotation_heatmap_df_select <- annotation_heatmap_df %>% select(-Subject_ID)
  
  
  
  ##This should be zero, means the rows of the annotation perfectly match
  # the columns of the data 
  sum(rownames(annotation_heatmap_df_select) != colnames(heat_map_data))
  
  
  colours <- list(
    Cancer = c('HD'= '#ADD8E6', 'Melanoma' = '#FFB6C1'),
    Anti_PD_1 = c('0' = 'black', '1' = '#C3B1E1', 'HD' = 'white' ),
    Anti_CTLA_4 = c('0' = 'black', '1' = '#C3B1E1', 'HD' = 'white'),
    Responder = c('NR' = '#b03608', 'R' = '#738678','HD' = 'gray' )
  )
  
  
  col_heat_map_annotation <- HeatmapAnnotation(
    df = annotation_heatmap_df_select,
    col = colours, 
    annotation_height = 0.6,
    annotation_width = unit(1, 'cm'),
    gap = unit(1, 'mm'))
  
  
  
  ### color of the cells
  col_fun = colorRamp2(c(0,1.5,3,4.5,6), c("black",  "red", "orange", "yellow", "white"))
  
  
  ### Row annotations
  ### this is copied from the earlier chunk of code to ensure the same order as the data 
  color_levels <- c("#fc8c84", "#6454ac", "#5bc0de", "#1abc9c", "#fc6404", "#fcd444", "gray", "black")
  
  Annot <- rownames(heat_map_data) %>%
    data.frame() %>%
    rename(Protein = ".") %>%
    left_join(., annotation_df_v2) %>%
    #filter(reactome_group != "Other / Unknown") %>%
    mutate(annotation_and_reactome = ifelse(is.na(annotation_group), reactome_group, annotation_group)) %>%
    select(Protein, annotation_and_reactome) %>%
    mutate(annotation_and_reactome = ifelse(Protein == "IL12-sc", "Cytokines", annotation_and_reactome)) %>%
    #filter(annotation_and_reactome != "Other / Unknown") %>%
    mutate(annotation_and_reactome = factor(annotation_and_reactome, 
                                            levels = c(pathway_levels))) %>%
    mutate(colors = case_when(
      annotation_and_reactome %in% c("Checkpoints.Drugs", "Checkpoints.Inhibitory") ~ "#fc8c84",
      annotation_and_reactome %in% c("Chemokines", "Cytokines", "Interferons") ~ "#6454ac",
      annotation_and_reactome %in% c("Innate Immune System", "NKG2D.ligands", "Adaptive Immune System", "Immune (other)") ~ "#5bc0de",
      annotation_and_reactome %in% c("Developmental Biology", "ECM organization", "Extracellular matrix organization") ~ "#1abc9c",
      annotation_and_reactome %in% c("Ion channel transport","Cell-Cell communication",  "SLC-mediated transmembrane transport", "GPCRs", "Signal Transduction","Sensory Receptors", "Signaling by GPCR", "Transporters") ~ "#fc6404",
      annotation_and_reactome %in% c("Metabolism","Melanoma.Expressed", "Neuronal System","Growth factors", "Metabolism of proteins", "Proteases/Hydrolases","Hemostasis", "WNT.related") ~ "#fcd444",
      annotation_and_reactome == "Viral" ~ "black",
      annotation_and_reactome %in% c("Other / Unknown") ~ "gray",
      TRUE ~ NA_character_ # default case
    )) %>% 
    mutate(colors = factor(colors, levels = color_levels)) 
  
  
  # Here is defining of Annotation fo CATEGORIES
  CATEGORIES = rowAnnotation(annotation_and_reactome = Annot$annotation_and_reactome,
                             col = list(annotation_and_reactome = c(`Checkpoints.Drugs` = "#fc8c84",
                                                                    `Checkpoints.Inhibitory` = "#fc8c84",
                                                                    `Chemokines` = "#6454ac",
                                                                    `Cytokines` = "#6454ac",
                                                                    `Interferons` = "#6454ac",
                                                                    `Innate Immune System` = "#5bc0de",
                                                                    `Adaptive Immune System` = "#5bc0de",
                                                                    `Immune (other)` = "#5bc0de",
                                                                    `NKG2D.ligands` = "#5bc0de",
                                                                    `Developmental Biology` = "#1abc9c",
                                                                    `Cell-Cell communication` = "#fc6404",
                                                                    `ECM organization` = "#1abc9c",
                                                                    `Extracellular matrix organization` = "#1abc9c",
                                                                    `Growth.Factors` = "#fcd444",
                                                                    `Neuronal System` = "#fcd444",
                                                                    `Hemostasis` = "#fcd444",
                                                                    `Growth factors`= "#fcd444",
                                                                    `Ion channel transport` = "#fc6404",
                                                                    `Metabolism` = "#fcd444",
                                                                    `Metabolism of proteins` = "#fcd444",
                                                                    `Proteases/Hydrolases` = "#fcd444",
                                                                    `Melanoma.Expressed` = "#fcd444",
                                                                    `SLC-mediated transmembrane transport` = "#fc6404",
                                                                    `GPCRs` ="#fc6404",
                                                                    `Signal Transduction` = "#fc6404",
                                                                    `Signaling by GPCR` = "#fc6404",
                                                                    `Sensory Receptors` =  "#fc6404",
                                                                    `Transporters` = "#fc6404",
                                                                    `Other / Unknown` = "gray",
                                                                    `Viral` = "black",
                                                                    `WNT.related` = "#fcd444")),
                             gp = gpar(col = "white"),simple_anno_size = unit(0.4, "cm"),
                             annotation_name_gp= gpar(fontsize = 0))
  
  
  #setwd("C:/Users/pmunson/OneDrive - Nurix Therapeutics, Inc/Desktop/Autoantibodies/Graphs") 
  #png("heatmap_baseline_melanoma_vs_hd.png",width=8,height=8,units="in",res=600)
  
  heat_map_up_in_NRs <- Heatmap(heat_map_data, 
                                col = col_fun,
                                heatmap_legend_param = list(
                                  title = "REAP score",
                                  direction = "horizontal"),
                                # annotation_legend_param = list(
                                #   title = "Cancer",
                                #   direction = "horizontal"),
                                cluster_rows = FALSE,
                                show_column_names = FALSE,
                                cluster_columns = FALSE,
                                clustering_distance_columns = "euclidean",
                                row_names_gp = gpar(fontsize = 8), 
                                column_names_gp = gpar(fontsize = 0),
                                #rect_gp = gpar(col = "white", lwd = 0.01),
                                top_annotation = col_heat_map_annotation,
                                row_split = Annot$color,
                                row_title_gp = gpar(fontsize = 0),
                                column_title_gp = gpar(fontsize = 0),
                                #annotation_row = annotation_row,
                                column_title = "", 
                                column_split = annotation_heatmap_df_select$Responder,
                                right_annotation = CATEGORIES)
  }
  
  {volcano_proteins_responders <- subset(raw_statistics_responders_day0_df %>%
                                           ## removing the Epitope
                                           filter(!grepl("_Epitope", Protein)), 
                                         #`p_value` < 0.25  & lfc >1 &
                                         ### this p value should detect if its in at least 2 subjects
                                         `p_value` < 0.4&
                                           lfc < 1 &
                                           abs(max_diff) > 1)
    ### trying to remove proteins that might be noisey
    
    
    detectable_in_healthy_donors <- raw_meta_df %>%
      filter(Cancer == "HD") %>%
      filter(value > 1) %>%
      pull(Protein) %>% unique()
    
    setdiff(volcano_proteins_responders %>% pull(Protein),detectable_in_healthy_donors)
    
    volcano_proteins_responders <-volcano_proteins_responders %>%
      filter(!Protein %in% detectable_in_healthy_donors)
    
    ### Manually assigning functions to potential proteins of interest
    unknown_proteins <- volcano_proteins_responders %>%
      left_join(., annotation_df_v2) %>%
      #filter(reactome_group != "Other / Unknown") %>%
      mutate(annotation_and_reactome = ifelse(is.na(annotation_group), reactome_group, annotation_group)) %>%
      select(Protein, annotation_and_reactome) %>%
      filter(annotation_and_reactome == "Other / Unknown") %>%
      pull(Protein)
    
    
    ### select desired order of proteins based on pathway similarites on the heat map
    
    
    pathway_levels <-  volcano_proteins_responders%>%
      left_join(., annotation_df_v2) %>%
      mutate(annotation_and_reactome = ifelse(is.na(annotation_group), reactome_group, annotation_group)) %>%
      select(Protein, annotation_and_reactome) %>%
      mutate(annotation_and_reactome = ifelse(Protein == "IL12-sc", "Cytokines", annotation_and_reactome)) %>%
      filter(annotation_and_reactome!= "Other / Unknown") %>%
      unique() %>%
      #filter(annotation_and_reactome != "Other / Unknown") %>%
      arrange(annotation_and_reactome) %>%
      select(annotation_and_reactome) %>%
      unique() %>%
      mutate(value = 0) %>%
      #print(n = Inf)
      pivot_wider(., names_from = annotation_and_reactome, values_from = value) %>%
      ## define desired order here....
      select(any_of(c("Checkpoints.Drugs", "Checkpoints.Inhibitory","Chemokines", "Cytokines","Interferons", "Innate Immune System",
                      "NKG2D.ligands","Adaptive Immune System", "Immune (other)",
                      everything())))%>%
      colnames() 
    
    heatmap_proteins <-volcano_proteins_responders %>%
      left_join(., annotation_df_v2) %>%
      mutate(annotation_and_reactome = ifelse(is.na(annotation_group), reactome_group, annotation_group)) %>%
      select(Protein, annotation_and_reactome) %>%
      mutate(annotation_and_reactome = ifelse(Protein == "IL12-sc", "Cytokines", annotation_and_reactome)) %>%
      filter(annotation_and_reactome != "Other / Unknown") %>% 
      #filter(annotation_and_reactome != "Other / Unknown") %>%
      mutate(annotation_and_reactome = factor(annotation_and_reactome, levels = c(pathway_levels))) %>%
      arrange(annotation_and_reactome) %>%
      filter(!is.na(annotation_and_reactome)) %>%
      #filter(Protein != "229E-RBD") %>%
      pull(Protein)
    
    heat_map_data <- raw_meta_df %>%
      filter(Visit_Name_time == "DAY_0") %>%
      select(Subject_ID, Protein, value) %>%
      filter(Protein %in% heatmap_proteins) %>%
      pivot_wider(., names_from = Protein, values_from = value) %>%
      select(-Subject_ID) %>%
      #scale() %>%
      t() 
    
    heat_map_data <-heat_map_data[heatmap_proteins,]
    
    colnames_heat_map <- raw_meta_df %>%
      filter(Visit_Name_time == "DAY_0") %>%
      select(Subject_ID, Protein, value) %>%
      filter(Protein %in% heatmap_proteins) %>%
      pivot_wider(., names_from = Protein, values_from = value) %>%
      select(Subject_ID) %>% pull()
    
    colnames(heat_map_data) <- colnames_heat_map
    
    ### Annotation data
    
    unique(raw_meta_df$Best_Clinical_outcome)
    
    annotation_heatmap <-raw_meta_df%>%
      mutate(Responder = case_when(
        Best_Clinical_outcome %in% c("PR", "SD", "NED1") ~ "R",
        Best_Clinical_outcome %in% c("NED2", "PD") ~ "NR",
        Best_Clinical_outcome %in% c("HD") ~ "HD",
        TRUE ~ NA_character_  # Handle other cases if any
      )) %>%
      mutate(Responder = factor(Responder, levels = c("HD", "R", "NR"))) %>%
      #select(Subject_ID, Best_Clinical_outcome, Overall_Survival,PFS, Metastasis_Site, Anti_CTLA_4, Anti_PD_1, Prior_IFN, Prior_IL2, Cancer) 
      select(Subject_ID,Cancer,  Anti_CTLA_4,Anti_PD_1, Responder ) %>%
      mutate(Anti_CTLA_4 = factor(Anti_CTLA_4, levels = c("HD", "0", "1"))) %>%
      select(Subject_ID, Responder)
    
    
    # Reorder the annotation data to match the matrix column names
    #annotation_data <- 
    
    annotation_heatmap <- annotation_heatmap[match(colnames_heat_map, annotation_heatmap$Subject_ID), ]
    
    
    ## If this is zero, it means the 
    ## annotation rows are in the same order as the data columns 
    sum(annotation_heatmap$Subject_ID!= colnames_heat_map)
    
    ## checking that the heatmap and the annotation table are compatible in dimensions
    dim(heat_map_data)[2] == dim(annotation_heatmap)[1]
    
    ## I need to make the subject_visit a rowname
    annotation_heatmap_df <- as.data.frame(annotation_heatmap)
    rownames(annotation_heatmap_df) <-  annotation_heatmap_df[,1] 
    annotation_heatmap_df_select <- annotation_heatmap_df %>% select(-Subject_ID)
    
    
    
    ##This should be zero, means the rows of the annotation perfectly match
    # the columns of the data 
    sum(rownames(annotation_heatmap_df_select) != colnames(heat_map_data))
    
    
    colours <- list(
      Cancer = c('HD'= '#ADD8E6', 'Melanoma' = '#FFB6C1'),
      Anti_PD_1 = c('0' = 'black', '1' = '#C3B1E1', 'HD' = 'white' ),
      Anti_CTLA_4 = c('0' = 'black', '1' = '#C3B1E1', 'HD' = 'white'),
      Responder = c('NR' = '#b03608', 'R' = '#738678','HD' = 'gray' )
    )
    
    
    col_heat_map_annotation <- HeatmapAnnotation(
      df = annotation_heatmap_df_select,
      col = colours, 
      annotation_height = 0.6,
      annotation_width = unit(1, 'cm'),
      gap = unit(1, 'mm'))
    
    
    
    ### color of the cells
    col_fun = colorRamp2(c(0,1.5,3,4.5,6), c("black",  "red", "orange", "yellow", "white"))
    
    
    ### Row annotations
    ### this is copied from the earlier chunk of code to ensure the same order as the data 
    color_levels <- c("#fc8c84", "#6454ac", "#5bc0de", "#1abc9c", "#fc6404", "#fcd444", "gray", "black")
    
    Annot <- rownames(heat_map_data) %>%
      data.frame() %>%
      rename(Protein = ".") %>%
      left_join(., annotation_df_v2) %>%
      #filter(reactome_group != "Other / Unknown") %>%
      mutate(annotation_and_reactome = ifelse(is.na(annotation_group), reactome_group, annotation_group)) %>%
      select(Protein, annotation_and_reactome) %>%
      mutate(annotation_and_reactome = ifelse(Protein == "IL12-sc", "Cytokines", annotation_and_reactome)) %>%
      #filter(annotation_and_reactome != "Other / Unknown") %>%
      mutate(annotation_and_reactome = factor(annotation_and_reactome, 
                                              levels = c(pathway_levels))) %>%
      mutate(colors = case_when(
        annotation_and_reactome %in% c("Checkpoints.Drugs", "Checkpoints.Inhibitory") ~ "#fc8c84",
        annotation_and_reactome %in% c("Chemokines", "Cytokines", "Interferons") ~ "#6454ac",
        annotation_and_reactome %in% c("Innate Immune System", "NKG2D.ligands", "Adaptive Immune System", "Immune (other)") ~ "#5bc0de",
        annotation_and_reactome %in% c("Developmental Biology", "ECM organization", "Extracellular matrix organization") ~ "#1abc9c",
        annotation_and_reactome %in% c("Ion channel transport","Cell-Cell communication",  "SLC-mediated transmembrane transport", "GPCRs", "Signal Transduction","Sensory Receptors", "Signaling by GPCR", "Transporters") ~ "#fc6404",
        annotation_and_reactome %in% c("Metabolism","Melanoma.Expressed", "Neuronal System","Growth factors", "Metabolism of proteins", "Proteases/Hydrolases","Hemostasis", "WNT.related") ~ "#fcd444",
        annotation_and_reactome == "Viral" ~ "black",
        annotation_and_reactome %in% c("Other / Unknown") ~ "gray",
        TRUE ~ NA_character_ # default case
      )) %>% 
      mutate(colors = factor(colors, levels = color_levels)) 
    
    
    # Here is defining of Annotation fo CATEGORIES
    CATEGORIES = rowAnnotation(annotation_and_reactome = Annot$annotation_and_reactome,
                               col = list(annotation_and_reactome = c(`Checkpoints.Drugs` = "#fc8c84",
                                                                      `Checkpoints.Inhibitory` = "#fc8c84",
                                                                      `Chemokines` = "#6454ac",
                                                                      `Cytokines` = "#6454ac",
                                                                      `Interferons` = "#6454ac",
                                                                      `Innate Immune System` = "#5bc0de",
                                                                      `Adaptive Immune System` = "#5bc0de",
                                                                      `Immune (other)` = "#5bc0de",
                                                                      `NKG2D.ligands` = "#5bc0de",
                                                                      `Developmental Biology` = "#1abc9c",
                                                                      `Cell-Cell communication` = "#fc6404",
                                                                      `ECM organization` = "#1abc9c",
                                                                      `Extracellular matrix organization` = "#1abc9c",
                                                                      `Growth.Factors` = "#fcd444",
                                                                      `Neuronal System` = "#fcd444",
                                                                      `Hemostasis` = "#fcd444",
                                                                      `Growth factors`= "#fcd444",
                                                                      `Ion channel transport` = "#fc6404",
                                                                      `Metabolism` = "#fcd444",
                                                                      `Metabolism of proteins` = "#fcd444",
                                                                      `Proteases/Hydrolases` = "#fcd444",
                                                                      `Melanoma.Expressed` = "#fcd444",
                                                                      `SLC-mediated transmembrane transport` = "#fc6404",
                                                                      `GPCRs` ="#fc6404",
                                                                      `Signal Transduction` = "#fc6404",
                                                                      `Signaling by GPCR` = "#fc6404",
                                                                      `Sensory Receptors` =  "#fc6404",
                                                                      `Transporters` = "#fc6404",
                                                                      `Other / Unknown` = "gray",
                                                                      `Viral` = "black",
                                                                      `WNT.related` = "#fcd444")),
                               gp = gpar(col = "white"),simple_anno_size = unit(0.4, "cm"),
                               annotation_name_gp= gpar(fontsize = 0))
    
    
    #setwd("C:/Users/pmunson/OneDrive - Nurix Therapeutics, Inc/Desktop/Autoantibodies/Graphs") 
    #png("heatmap_baseline_melanoma_vs_hd.png",width=8,height=8,units="in",res=600)
    
    heat_map_up_in_Rs <- Heatmap(heat_map_data, 
                                 col = col_fun,
                                 heatmap_legend_param = list(
                                   title = "REAP score",
                                   direction = "horizontal"),
                                 # annotation_legend_param = list(
                                 #   title = "Cancer",
                                 #   direction = "horizontal"),
                                 cluster_rows = FALSE,
                                 show_column_names = FALSE,
                                 cluster_columns = FALSE,
                                 clustering_distance_columns = "euclidean",
                                 row_names_gp = gpar(fontsize = 8), 
                                 column_names_gp = gpar(fontsize = 0),
                                 #rect_gp = gpar(col = "white", lwd = 0.01),
                                 top_annotation = col_heat_map_annotation,
                                 row_split = Annot$color,
                                 row_title_gp = gpar(fontsize = 0),
                                 column_title_gp = gpar(fontsize = 0),
                                 #annotation_row = annotation_row,
                                 column_title = "", 
                                 column_split = annotation_heatmap_df_select$Responder,
                                 right_annotation = CATEGORIES)
  }
  setwd("C:/Users/pmunson/OneDrive - Nurix Therapeutics, Inc/Desktop/Autoantibodies/Graphs") 
  png("heatmap_baseline_melanoma_R_vs_NR.png",width=8,height=8,units="in",res=600)
  draw(heat_map_up_in_NRs %v% heat_map_up_in_Rs, merge_legend = TRUE)
  dev.off()
  
}

{
raw_statistics_responders_day0_df <- raw_meta_df %>%
  filter(Visit_Name_time == "DAY_0") %>%
  filter(Cancer == "Melanoma") %>%
  mutate(Responder = case_when(
    Best_Clinical_outcome %in% c("PR", "SD", "NED1") ~ "R",
    Best_Clinical_outcome %in% c("NED2", "PD") ~ "NR",
    TRUE ~ NA_character_  # Handle other cases if any
  )) %>%
  select(Responder, Protein, value) %>%
  group_by(Protein) %>%
  #### Calculating statistics
  summarise(
    lfc = log2(mean(value[Responder == "NR"]) / mean(value[Responder == "R"])),
    median_diff = median(value[Responder == "NR"])-median(value[Responder == "R"]),
    mean_diff = mean(value[Responder == "NR"])-mean(value[Responder == "R"]),
    max_value_Responder = max(value[Responder == "NR"]),
    max_value_R = max(value[Responder == "R"]),
    p_value = t.test(value[Responder == "NR"], value[Responder == "R"], na.action = na.omit)$p.value
  ) %>%
  mutate(max_diff = max_value_Responder - max_value_R)

  ####
  
  detectable_in_healthy_donors <- raw_meta_df %>%
    filter(Cancer == "HD") %>%
    filter(value > 1) %>%
    pull(Protein) %>% unique()
  

  
  volcano_proteins_responders <-volcano_proteins_responders %>%
    filter(!Protein %in% detectable_in_healthy_donors)
  

  
baseline_responder_volcano <- raw_statistics_responders_day0_df %>%
  filter(!grepl("_Epitope", Protein)) %>%
  left_join(., annotation_df_v2) %>%
  #filter(lfc != Inf & lfc != -Inf) %>%
  ggplot(., aes(x = lfc, y = -log(p_value), color = max_diff))+
  geom_point(size =2) +
  theme_prism() +
  #xlim(-20,20) +
  #ylim(0,10.9) +
  ylab("-Log(p value)") +
  xlab("Log2 fold change") +
  #ggtitle("Baseline") +
    #Responder = c('NR' = '#b03608', 'R' = '#738678','HD' = 'gray' )
  scale_color_gradient2(low = "green", mid = "gray", high = "#b03608", midpoint = 0) +
  #geom_text(x = Inf, y = Inf, label = "Up in Melanoma", hjust = 1.5, vjust = 1, color = "red", size =5) +
  #geom_text(x = -Inf, y = Inf, label = "Up in Healthy Donors", hjust = -0.75, vjust = 1, color = "blue", size = 5) +
  theme(legend.position = "none") +
  geom_text_repel(data = subset(raw_statistics_responders_day0_df %>%
                                  ## removing the Epitope
                                  filter(!grepl("_Epitope", Protein)) %>%
                                  filter(!Protein %in% detectable_in_healthy_donors), 
                                #`p_value` < 0.25  & lfc >1 &
                                ### this p value should detect if its in at least 2 subjects
                                `p_value` < 0.3 &
                                  abs(lfc) > 1 &
                                  abs(max_diff) > 1), 
                  aes(label = Protein),
                  size = 4,
                  direction = "both",
                  max.overlaps = 20,min.segment.length = 0.01, box.padding = 0.25, color = "black")
  
 setwd("C:/Users/pmunson/OneDrive - Nurix Therapeutics, Inc/Desktop/Autoantibodies/Graphs") 
  ggsave(filename = "baseline_responder_volcano.png", plot = baseline_responder_volcano, width = 8, height = 6, dpi = 600)

}

responder_checkpoint_graph <- raw_meta_df %>%
  mutate(Responder = case_when(
    Best_Clinical_outcome %in% c("PR", "SD", "NED1") ~ "R",
    Best_Clinical_outcome %in% c("NED2", "PD") ~ "NR",
    TRUE ~ NA_character_  # Handle other cases if any
  )) %>%
  filter(!is.na(Responder)) %>%
  filter(Visit_Name_time == "DAY_0") %>%
  filter(Protein == "PDCD1" | Protein == "CTLA4") %>%
  #filter(!is.na(Responder)) %>%
  ggplot(., aes(x = Responder, y = value, color = Responder)) +
  geom_hline(yintercept = 7, color = "white") +
  geom_jitter(width = 0.15, size =3) +
  theme_prism() +
  ylab("REAP score") +
  theme(legend.position = "none") +
  scale_color_manual(values = c('#b03608','#738678')) +
  facet_wrap(~Protein, scale = "free")  +
  xlab("")


raw_meta_df %>%
  mutate(Responder = case_when(
    Best_Clinical_outcome %in% c("PR", "SD", "NED1") ~ "R",
    Best_Clinical_outcome %in% c("NED2", "PD") ~ "NR",
    TRUE ~ NA_character_  # Handle other cases if any
  )) %>%
  filter(!is.na(Responder)) %>%
  filter(Visit_Name_time == "DAY_0") %>%
  filter(Protein == "PDCD1" | Protein == "CTLA4") %>%
  #filter(!is.na(Responder)) %>%
  ggplot(., aes(x = Responder, y = value, color = Responder)) +
  geom_hline(yintercept = 7, color = "white") +
  geom_jitter(width = 0.15, size =3) +
  theme_prism() +
  ylab("REAP score") +
  theme(legend.position = "none") +
  scale_color_manual(values = c('#b03608','#738678')) +
  facet_wrap(~Protein, scale = "free")  +
  xlab("") +
  geom_text_repel(data = raw_meta_df %>%
                    mutate(Responder = case_when(
                      Best_Clinical_outcome %in% c("PR", "SD", "NED1") ~ "R",
                      Best_Clinical_outcome %in% c("NED2", "PD") ~ "NR",
                      TRUE ~ NA_character_  # Handle other cases if any
                    )) %>%
                    filter(!is.na(Responder)) %>%
                    filter(Visit_Name_time == "DAY_0") %>%
                    filter(Protein == "PDCD1" | Protein == "CTLA4") %>%
                    filter(value >0.5), aes(label = Subject_ID), size = 4)

#setwd("C:/Users/pmunson/OneDrive - Nurix Therapeutics, Inc/Desktop/Autoantibodies/Graphs") 
#ggsave(filename = "responder_checkpoint_graph.png", plot = responder_checkpoint_graph, width = 4, height = 4, dpi = 600)
}

## PhIP seq


### Reading in the raw phip seq data
Phip_seq_raw <- read_csv("C:/Users/pmunson/OneDrive - Nurix Therapeutics, Inc/Desktop/Autoantibodies/downloads_study310_mapped_rpk_gene.csv")




### They did two experiments, one with 2 rounds of amplification and one with 3 rounds of amplification
Phip_seq_intermediate_v1 <- Phip_seq_raw %>%
  pivot_longer(cols = -c(gene)) %>%
  mutate(
    Subject_ID = gsub(".*?(PT_\\d+|HD\\d+|HC\\d+)_.*", "\\1", name),
    Subject_ID = ifelse(startsWith(Subject_ID, "PT_"), gsub("^(PT_\\d+)_.*", "\\1", Subject_ID), Subject_ID)
  ) %>%
  mutate(Subject_ID = toupper(Subject_ID)) %>%
  mutate(
    Subject_ID = ifelse(startsWith(Subject_ID, "PT_"), 
                        gsub("^(PT_[^_]+)_[^_]+.*", "\\1", Subject_ID),
                        Subject_ID)
  ) %>%
  mutate(
    Visit_Name_time = ifelse(grepl("day", name, ignore.case = TRUE), 
                             gsub(".*day_(\\d+)_.*", "DAY_\\1", name),
                             "DAY_0"),
    Lysate = ifelse(grepl("Round2", name), "Round2",
                    ifelse(grepl("Round3", name), "Round3", NA))
  ) 



Phip_seq_intermediate_v2 <- Phip_seq_intermediate_v1 %>%
  select(-c(name)) %>%
  #pivot_wider(., names_from = gene, values_from = value) %>%
  mutate(Visit_Name_time = gsub("DAY_89", "DAY_89/101", Visit_Name_time)) %>%
  mutate(Subject_Visit = str_c(Subject_ID, Visit_Name_time, sep = "_")) %>%
  rename(Protein = gene) %>%
  select(Subject_ID, Visit_Name_time, Subject_Visit, Lysate, Protein, value) %>%
  #### Filtering out the AG1 and GFAPs, I think these are some kind of controls
  filter(!grepl("^AG|^GFAP", Subject_ID))


## adding in the metdata
Phip_seq_final <- Phip_seq_intermediate_v2 %>%
  left_join(., metadata_df_v2) %>% 
  #### The HC columns are extra healthy controls that the DeRisi lab added in additon to our HD1-5 controls
  mutate(Cancer = ifelse(grepl("^HC", Subject_ID), "HD", Cancer)) %>%
  ### There are some proteins that look like dates? It must be some kind of artifact
  filter(!grepl("^\\d", Protein))


## calculating statistics 
### Baseline

Phip_statistics_df <- Phip_seq_final %>%
  filter(Lysate == "Round3") %>%
  filter(Visit_Name_time == "DAY_0") %>%
  select(Cancer, Protein, value) %>%
  group_by(Protein) %>%
  #### Calculating statistics
  summarise(
    lfc = log2(mean(value[Cancer == "Melanoma"]) / mean(value[Cancer == "HD"])),
    median_diff = median(value[Cancer == "Melanoma"])-median(value[Cancer == "HD"]),
    mean_diff = mean(value[Cancer == "Melanoma"])-mean(value[Cancer == "HD"]),
    max_value_cancer = max(value[Cancer == "Melanoma"]),
    max_value_hd = max(value[Cancer == "HD"]),
    p_value = t.test(value[Cancer == "Melanoma"], value[Cancer == "HD"], na.action = na.omit)$p.value
  ) %>%
  mutate(max_diff = max_value_cancer - max_value_hd) 

### volcano

## removing the Epitope



Phip_statistics_df %>%
  filter(max_value_cancer >50 | max_value_hd > 50) %>% 
  filter(!grepl("L0C", Protein)) %>%
  ggplot(., aes(x = lfc, y = -log(p_value), color = max_diff))+
  geom_point(size =2) +
  theme_prism() +
  #xlim(-20,20) +
  #ylim(0,10.9) +
  ylab("-Log(p value)") +
  xlab("Log2 fold change") +
  scale_color_gradient2(low = "gray",mid = "gray",  high = "red", midpoint = 0) +
  geom_text_repel(data = subset(Phip_statistics_df, 
                              `p_value` < 0.05  & lfc >1  & max_value_cancer >50 & max_value_hd < 50
                              | is.infinite(lfc) & lfc > 0 & max_value_cancer >50 & max_value_hd < 50|
                              #p_value < 0.001 |
                              max_diff >1000 & max_value_hd < 50),
                              # max_value_cancer >50), 
                aes(label = Protein),
                size = 3,
                direction = "both",
                max.overlaps = 40,min.segment.length = 0.01, box.padding = 0.25, color = "black") +
  theme(legend.position = "none") 

  

up_in_melanoma_phip <- Phip_statistics_df %>%
    filter(p_value < 0.05) %>%
    filter(lfc >1) %>%
    arrange(desc(max_value_cancer)) %>%
    filter(max_value_cancer >50) %>% 
    pull(Protein)

### MAGED1 is in melanoma
### CALR3 is a cancer testis antigen

### Now that I have a curated list of signal versus healthy people at baseline
## I'm going to have to manually go over this list, to find patterns
Phip_statistics_df %>%
  filter(`p_value` < 0.05  & lfc >1  & max_value_cancer >50 & max_value_hd < 50
         | is.infinite(lfc) & lfc > 0 & max_value_cancer >50 & max_value_hd < 50|
           #p_value < 0.001 |
           max_diff >1000 & max_value_hd < 50) %>%
  print(n = Inf)

  Phip_seq_final %>%
    filter(Lysate == "Round3") %>%
    filter(Visit_Name_time == "DAY_0") %>%
    #filter(Protein %in% up_in_melanoma_phip) %>%
    filter(Protein == "MGMT") %>%
    ggplot(., aes(x = Cancer, y = value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = Best_Clinical_outcome), width = 0.015) +
    facet_wrap(~Protein, scales = "free") +
    theme_prism()
  

  Phip_seq_final %>%
    filter(Lysate == "Round3") %>%
    #filter(Visit_Name_time == "DAY_0") %>%
    #filter(Protein %in% up_in_melanoma_phip) %>%
    filter(Protein == "TAF7L") %>%
    ggplot(., aes(x = Visit_Name_time, y = value, group = Subject_ID)) +
    geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
    geom_line() +
    geom_point() +
    facet_wrap(~Protein + Cancer) +
    theme_prism()

  ### Potentially vaccine induced 
cta_list <- c("ADAM2","ARMC3", "ATAD2", "CABYR",
  "CALR3","CASC5", "CTAG2","IGSF11",
  "LUZP4","OIP5", "PAGE3", "PIWIL2","POTEA", "POTED",
  "RGS22", "SAGE1", "SPAG17", "SPAG8","SSX4","TDRD1",
  "THEG", "TMEFF1", "TSGA10", "XAGE2", "IGSF11",
  ### Detectable at baseline
  "CEP290", "CPXCR1","CTCFL","CTNNA2", "MAEL","MORC1",
  "NOL4", "ODF2","OTOA", "PRSS54","SPAG4", "SPAG9",
  "SSX1","SSX5", "TEKT5", "TULP2", "ZNF645", "LYPD6B")
  
  ### Cancer testis antigens
  # http://www.cta.lncc.br/



cta_day_0 <- Phip_seq_final %>%
  filter(Protein %in% cta_list) %>%
  group_by(Subject_Visit) %>%
  summarise(CTA_hits = sum(`value` > 50, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    Subject_ID = gsub(".*?(PT_\\d+|HD\\d+|HC\\d+)_.*", "\\1", Subject_Visit),
    Subject_ID = ifelse(startsWith(Subject_ID, "PT_"), gsub("^(PT_\\d+)_.*", "\\1", Subject_ID), Subject_ID)
  ) %>% # Extracting Visit_Name_time
  mutate(Visit_Name_time = sub("^.*?(DAY.*)", "\\1", Subject_Visit)) %>%
  left_join(., metadata_df_v2) %>%
  mutate(Cancer = ifelse(is.na(Cancer), "HD", Cancer)) %>%
  filter(Visit_Name_time == "DAY_0") 
 
cta_day_0 %>% 
ggplot(., aes(x = Cancer, y = CTA_hits, color = Cancer)) +
  geom_boxplot() +
  stat_compare_means(method = "t.test", label = "p.format",
                     label.y = 6,
                     size= 6,
                     comparisons = list(c("HD", "Melanoma")))
 
 #### CTA hits over time 

cta_post_vaccination <- Phip_seq_final %>%
  filter(Protein %in% cta_list) %>%
  group_by(Subject_Visit) %>%
  summarise(CTA_hits = sum(`value` > 50, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    Subject_ID = gsub(".*?(PT_\\d+|HD\\d+|HC\\d+)_.*", "\\1", Subject_Visit),
    Subject_ID = ifelse(startsWith(Subject_ID, "PT_"), gsub("^(PT_\\d+)_.*", "\\1", Subject_ID), Subject_ID)
  ) %>% # Extracting Visit_Name_time
  mutate(Visit_Name_time = sub("^.*?(DAY.*)", "\\1", Subject_Visit)) %>%
  left_join(., metadata_df_v2) %>%
  mutate(Cancer = ifelse(is.na(Cancer), "HD", Cancer)) %>%
  filter(Cancer == "Melanoma") %>%
  mutate(Visit_Name_time = factor(Visit_Name_time, levels = c("DAY_0", "DAY_43", "DAY_89/101")))







cta_post_vaccination %>%
  filter(Prior_Checkpoint == "Yes") %>%
  pairwise_t_test(
    CTA_hits ~ Visit_Name_time, paired = FALSE,
    p.adjust.method = "bonferroni"
  )

cta_post_vaccination %>%
  select(Subject_ID, Visit_Name_time, Prior_Checkpoint, Checkpoint, CTA_hits) %>%
  pivot_wider(., names_from = Visit_Name_time, values_from = CTA_hits) %>%
  arrange(Prior_Checkpoint)


subjects_both_timepoints <- cta_post_vaccination %>%
  filter(Visit_Name_time %in% c("DAY_0", "DAY_89/101")) %>%
  group_by(Subject_ID) %>%
  filter(all(c("DAY_0", "DAY_89/101") %in% Visit_Name_time)) %>%
  pull(Subject_ID)

cta_post_vaccination %>% 
  filter(Subject_ID %in% subjects_both_timepoints) %>%
  ggplot(., aes(x = Visit_Name_time, y = CTA_hits, color = Checkpoint)) +
  geom_boxplot(aes(group = Visit_Name_time), color = "black", outlier.shape = NA)+ 
  geom_hline(yintercept = 7, color = "white") +
  geom_jitter(width = 0.15, size = 3) +
  facet_wrap(~Prior_Checkpoint,scales = "free") +
  xlab("") +
  theme_prism() +
  stat_compare_means(comparisons = list(c("DAY_89/101", "DAY_0")), method = "t.test", paired = TRUE)  # Use t-test for comparisons

### volcano plot of day 0 vs. day 89/101  


Phip_statistics_df_timepoint <- Phip_seq_final %>%
  filter(Lysate == "Round3") %>%
  filter(Visit_Name_time == "DAY_0") %>%
  select(Cancer, Protein, value) %>%
  group_by(Protein) %>%
  #### Calculating statistics
  summarise(
    lfc = log2(mean(value[Cancer == "Melanoma"]) / mean(value[Cancer == "HD"])),
    median_diff = median(value[Cancer == "Melanoma"])-median(value[Cancer == "HD"]),
    mean_diff = mean(value[Cancer == "Melanoma"])-mean(value[Cancer == "HD"]),
    max_value_cancer = max(value[Cancer == "Melanoma"]),
    max_value_hd = max(value[Cancer == "HD"]),
    p_value = t.test(value[Cancer == "Melanoma"], value[Cancer == "HD"], na.action = na.omit)$p.value
  ) %>%
  mutate(max_diff = max_value_cancer - max_value_hd) 


  Phip_statistics_df_timepoint <- Phip_seq_final %>%
    filter(Lysate == "Round3") %>%
    filter(Cancer == "Melanoma") %>%
    filter(Visit_Name_time == "DAY_0" |Visit_Name_time == "DAY_89/101" ) %>%
    select(Visit_Name_time, Protein, value) %>%
    group_by(Protein) %>%
    #### Calculating statistics
    summarise(
      lfc = log2(mean(value[Visit_Name_time == "DAY_89/101"]) / mean(value[Visit_Name_time == "DAY_0"])),
      median_diff = median(value[Visit_Name_time == "DAY_89/101"])-median(value[Visit_Name_time == "DAY_0"]),
      mean_diff = mean(value[Visit_Name_time == "DAY_89/101"])-mean(value[Visit_Name_time == "DAY_0"]),
      max_value_post_vax = max(value[Visit_Name_time == "DAY_89/101"]),
      max_value_baseline = max(value[Visit_Name_time == "DAY_0"]),
      p_value = t.test(value[Visit_Name_time == "DAY_89/101"], value[Visit_Name_time == "DAY_0"], na.action = na.omit)$p.value
    ) %>%
    mutate(max_diff = max_value_post_vax - max_value_baseline) 
  
  Phip_statistics_df_timepoint  %>%
    filter(max_value_post_vax >50 | max_value_baseline > 50) %>% 
    filter(!grepl("L0C", Protein)) %>%
    ggplot(., aes(x = lfc, y = -log(p_value), color = max_diff))+
    geom_point(size =2) +
    theme_prism() +
    #xlim(-20,20) +
    #ylim(0,10.9) +
    ylab("-Log(p value)") +
    xlab("Log2 fold change") +
    scale_color_gradient2(low = "gray",mid = "gray",  high = "red", midpoint = 0) +
    geom_text_repel(data = subset(Phip_statistics_df_timepoint, 
                                  `p_value` < 0.05  & lfc > 1  & max_value_post_vax >50 & max_value_baseline < 50
                                  | is.infinite(lfc) & lfc > 0 & max_value_post_vax >50 & max_value_baseline < 50|
                                    #p_value < 0.001 |
                                    max_diff >300 & max_value_baseline < 50),
                    # max_value_post_vax >50), 
                    aes(label = Protein),
                    size = 3,
                    direction = "both",
                    max.overlaps = 100,min.segment.length = 0.01, box.padding = 0.25, color = "black") +
    theme(legend.position = "none")   
  
  
  
  detected_in_hd <- Phip_seq_final %>%
    filter(Lysate == "Round3") %>%
    filter(Cancer != "Melanoma") %>%
    filter(value >50) %>%
    select(Protein) %>% unique() %>% pull()
  
  
  
  induced_proteins <- Phip_statistics_df_timepoint %>%
    filter(`p_value` < 0.05 & lfc > 1 & max_value_post_vax > 50 & max_value_baseline < 50 |
             is.infinite(lfc) & lfc > 0 & max_value_post_vax > 50 & max_value_baseline < 50 |
             max_diff > 100 & max_value_baseline < 50) %>%
    filter(!(Protein %in% detected_in_hd)) %>%
    select(Protein)
  
  induced_proteins %>% pull(Protein) 
  
  
  annotation_df_v2 %>%
    filter(Protein %in% induced_proteins)
  
  setwd("C:/Users/pmunson/OneDrive - Nurix Therapeutics, Inc/Desktop/Autoantibodies/")
  #write_xlsx(induced_proteins, "C:/Users/pmunson/OneDrive - Nurix Therapeutics, Inc/Desktop/Autoantibodies/induced_proteins.xlsx")
  
  
subjects_that_have_d43 <- Phip_seq_final %>%
    select(Subject_ID, Visit_Name_time) %>%
    filter(Visit_Name_time == "DAY_0" | Visit_Name_time == "DAY_43" | Visit_Name_time == "DAY_89/101" ) %>%
    unique() %>% 
    group_by(Subject_ID)  %>%
    count() %>%
    filter(n >= 2) %>%
    pull(Subject_ID)
  


  Phip_seq_final %>%
    filter(Lysate == "Round3") %>%
    #filter(Visit_Name_time == "DAY_0") %>%
    filter(Cancer == "Melanoma") %>%
    filter(Subject_ID %in% subjects_that_have_d43) %>%
    #filter(Visit_Name_time == "DAY_0" | Visit_Name_time == "DAY_43") %>%
    filter(Subject_ID == "PT_7") %>%
    #filter(Protein %in% up_in_melanoma_phip) %>%
    filter(Protein %in% induced_proteins$Protein) %>%
   # filter(Protein == "MAGEA10") %>%
    mutate(Subject_Protein = str_c(Subject_ID, Protein, sep = "_")) %>%
    mutate(Responder = case_when(
      Best_Clinical_outcome %in% c("PR", "SD", "NED1") ~ "R",
      Best_Clinical_outcome %in% c("NED2", "PD") ~ "NR",
      Best_Clinical_outcome %in% c("HD") ~ "HD",
      TRUE ~ NA_character_  # Handle other cases if any
    )) %>%
    mutate(Responder = factor(Responder, levels = c("HD", "R", "NR"))) %>%
    #filter(Checkpoint == "anti-CTLA-4 + anti-PD-1") %>%
    ggplot(., aes(x = Visit_Name_time, y = value, group = Subject_Protein, color = Checkpoint )) +
    geom_hline(yintercept = 250, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 50, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "white") +
    geom_hline(yintercept = 300, linetype = "dashed", color = "white") +
    geom_line() +
    geom_point() +
    facet_wrap(~Subject_ID) +
    theme_prism() +
    ylab("RPK") +
    xlab("") +
    geom_text_repel(data = subset(Phip_seq_final %>%
                                    filter(value >500) %>%
                                    filter(Lysate == "Round3") %>%
                                    #filter(Visit_Name_time == "DAY_0") %>%
                                    filter(Cancer == "Melanoma") %>%
                                    filter(Subject_ID %in% subjects_that_have_d43) %>%
                                    filter(Visit_Name_time == "DAY_89/101") %>%
                                    filter(Subject_ID == "PT_7") %>%
                                    #filter(Protein %in% up_in_melanoma_phip) %>%
                                    filter(Protein %in% induced_proteins$Protein)%>%
                                  mutate(Subject_Protein = str_c(Subject_ID, Protein, sep = "_"))), 
                    aes(label = Protein),
                    size = 3,
                    direction = "both",
                    max.overlaps = 100, min.segment.length = 0.01, box.padding = 0.25, color = "black")




### vaccine responses 
# Tyr = tyrosinase
# MLANA = MART-1
vaccine_proteins = c("TYR", "MLANA", "MAGEA6")

vaccine_responses <- Phip_seq_final %>%
  filter(Lysate == "Round3") %>%
  filter(Protein %in% vaccine_proteins) %>%
  mutate(Subject_Protein = str_c(Subject_ID, Protein, sep = "_")) %>%
  mutate(Best_Clinical_outcome = ifelse(is.na(Best_Clinical_outcome), "HD",Best_Clinical_outcome )) %>%
  mutate(Responder = case_when(
    Best_Clinical_outcome %in% c("PR", "SD", "NED1") ~ "R",
    Best_Clinical_outcome %in% c("NED2", "PD") ~ "NR",
    Best_Clinical_outcome %in% c("HD") ~ "HD",
    TRUE ~ NA_character_  # Handle other cases if any
  )) %>%
  mutate(Responder = factor(Responder, levels = c("HD", "R", "NR"))) %>%
  mutate(Timepoint = ifelse(Responder == "HD", "HD",
                            ifelse(Visit_Name_time == "DAY_0", 0,
                                   ifelse(Visit_Name_time == "DAY_43", 43,
                                          ifelse(Visit_Name_time == "DAY_89/101", "89/101", NA))))) %>%
  mutate(Timepoint = factor(Timepoint, levels = c("HD", "0", "43", "89/101"))) %>%
  ggplot(., aes(x = Timepoint, y = value, group = Subject_Protein, color = Responder)) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "white") +
  geom_hline(yintercept = 300, linetype = "dashed", color = "white") +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c('black','#738678','#b03608')) +
  geom_point(size =2) +
  facet_wrap(~Protein, scales = "free") +
  theme_prism() +
  ylab("RPK") +
  xlab("")

#setwd("C:/Users/pmunson/OneDrive - Nurix Therapeutics, Inc/Desktop/Autoantibodies/Graphs") 
#ggsave(filename = "vaccine_responses.png", plot = vaccine_responses, width = 9, height = 2.5, dpi = 600)

other_mages <- c("MAGEA4", "MAGEA10", "MAGEB1", "MAGEB3", "MAGEC1")

filtered_data <- Phip_seq_final %>%
  filter(Lysate == "Round3") %>%
  filter(Protein %in% other_mages) %>%
  mutate(Subject_Protein = str_c(Subject_ID, Protein, sep = "_")) %>%
  mutate(Best_Clinical_outcome = ifelse(is.na(Best_Clinical_outcome), "HD",Best_Clinical_outcome )) %>%
  mutate(Responder = case_when(
    Best_Clinical_outcome %in% c("PR", "SD", "NED1") ~ "R",
    Best_Clinical_outcome %in% c("NED2", "PD") ~ "NR",
    Best_Clinical_outcome %in% c("HD") ~ "HD",
    TRUE ~ NA_character_  # Handle other cases if any
  )) %>%
  mutate(Responder = factor(Responder, levels = c("HD", "R", "NR"))) %>%
  mutate(Timepoint = ifelse(Responder == "HD", "HD",
                            ifelse(Visit_Name_time == "DAY_0", 0,
                                   ifelse(Visit_Name_time == "DAY_43", 43,
                                          ifelse(Visit_Name_time == "DAY_89/101", "89/101", NA))))) %>%
  filter(value > 50) %>%
  mutate(Label = ifelse(Responder == "R" & Timepoint == "89/101", as.character(Subject_ID), 
                        ifelse(Responder == "NR" & Timepoint == "43", as.character(Subject_ID), "")))



graphing_df <- Phip_seq_final %>%
  filter(Lysate == "Round3") %>%
  filter(Protein %in% other_mages) %>%
  mutate(Subject_Protein = str_c(Subject_ID, Protein, sep = "_")) %>%
  mutate(Best_Clinical_outcome = ifelse(is.na(Best_Clinical_outcome), "HD",Best_Clinical_outcome )) %>%
  mutate(Responder = case_when(
    Best_Clinical_outcome %in% c("PR", "SD", "NED1") ~ "R",
    Best_Clinical_outcome %in% c("NED2", "PD") ~ "NR",
    Best_Clinical_outcome %in% c("HD") ~ "HD",
    TRUE ~ NA_character_  # Handle other cases if any
  )) %>%
  mutate(Responder = factor(Responder, levels = c("HD", "R", "NR"))) %>%
  mutate(Timepoint = ifelse(Responder == "HD", "HD",
                           ifelse(Visit_Name_time == "DAY_0", 0,
                                  ifelse(Visit_Name_time == "DAY_43", 43,
                                         ifelse(Visit_Name_time == "DAY_89/101", "89/101", NA))))) 
  


mage_graph <- graphing_df %>%
    mutate(Timepoint = factor(Timepoint, levels = c("HD", "0", "43", "89/101"))) %>%
  ggplot(., aes(x = Timepoint, y = value, group = Subject_Protein, color = Responder)) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "white") +
  geom_hline(yintercept = 700, linetype = "dashed", color = "white") +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c('black','#738678','#b03608')) +
  geom_point(size =2) +
  facet_wrap(~Protein, scales = "free", nrow = 1) +
  theme_prism() +
  ylab("RPK") +
    theme(legend.position = "none") +
  xlab("") +
  geom_text_repel(data = graphing_df %>%
                    filter(value >50 & Responder == "R" & Timepoint == "89/101" |
                             value >50 & Responder == "NR" & Timepoint == "43" ), 
                  max.overlaps = 20,min.segment.length = 0.01, box.padding = 2.5, color = "black",
                  aes(label = Subject_ID), size = 5)

#setwd("C:/Users/pmunson/OneDrive - Nurix Therapeutics, Inc/Desktop/Autoantibodies/Graphs") 
#ggsave(filename = "mage_graphs.png", plot = mage_graph, width = 13, height = 3, dpi = 600)

#### pt 7 deep dive

pt_7_graphing <- Phip_seq_final %>%
  filter(Subject_ID == "PT_7") %>%
  filter(Lysate == "Round3") %>%
  #filter(Visit_Name_time == "DAY_0") %>%
  filter(Cancer == "Melanoma") %>%
  filter(Subject_ID %in% subjects_that_have_d43) %>%
  #filter(Visit_Name_time == "DAY_0" | Visit_Name_time == "DAY_43") %>%
  #filter(Protein %in% up_in_melanoma_phip) %>%
  filter(Protein %in% induced_proteins$Protein | Protein %in% cta_list ) %>%
  # filter(Protein == "MAGEA10") %>%
  mutate(Subject_Protein = str_c(Subject_ID, Protein, sep = "_")) %>%
  mutate(Responder = case_when(
    Best_Clinical_outcome %in% c("PR", "SD", "NED1") ~ "R",
    Best_Clinical_outcome %in% c("NED2", "PD") ~ "NR",
    Best_Clinical_outcome %in% c("HD") ~ "HD",
    TRUE ~ NA_character_  # Handle other cases if any
  )) %>%
  mutate(Responder = factor(Responder, levels = c("HD", "R", "NR"))) %>%
  #filter(Checkpoint == "anti-CTLA-4 + anti-PD-1") %>%
  mutate(Timepoint = ifelse(Responder == "HD", "HD",
                            ifelse(Visit_Name_time == "DAY_0", 0,
                                   ifelse(Visit_Name_time == "DAY_43", 43,
                                          ifelse(Visit_Name_time == "DAY_89/101", "89/101", NA))))) 

pt_7_graph <- pt_7_graphing %>%
  ggplot(., aes(x = Timepoint, y = value, group = Subject_Protein, color = Responder)) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "white") +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~Subject_ID) +
  theme_prism() +
  ylab("RPK") +
  xlab("") +
  theme(legend.position = "none") +
  scale_color_manual(values = c('#738678','#b03608'))  +
  geom_text_repel(data = subset(pt_7_graphing %>% filter(value >500 & Timepoint == "89/101")), 
                  aes(label = Protein),
                  direction = "x",
                  max.overlaps = 100, min.segment.length = 0.01, box.padding = 2, color = "black", size =4)

#setwd("C:/Users/pmunson/OneDrive - Nurix Therapeutics, Inc/Desktop/Autoantibodies/Graphs") 
#ggsave(filename = "pt_7_graph.png", plot = pt_7_graph, width = 4, height = 4, dpi = 600)


### pt 10 case study 

pt_10_graphing <- Phip_seq_final %>%
  filter(Subject_ID == "PT_10") %>%
  filter(Lysate == "Round3") %>%
  #filter(Visit_Name_time == "DAY_0") %>%
  filter(Cancer == "Melanoma") %>%
  filter(Subject_ID %in% subjects_that_have_d43) %>%
  #filter(Visit_Name_time == "DAY_0" | Visit_Name_time == "DAY_43") %>%
  #filter(Protein %in% up_in_melanoma_phip) %>%
  filter(Protein %in% induced_proteins$Protein | Protein %in% cta_list ) %>%
  # filter(Protein == "MAGEA10") %>%
  mutate(Subject_Protein = str_c(Subject_ID, Protein, sep = "_")) %>%
  mutate(Responder = case_when(
    Best_Clinical_outcome %in% c("PR", "SD", "NED1") ~ "R",
    Best_Clinical_outcome %in% c("NED2", "PD") ~ "NR",
    Best_Clinical_outcome %in% c("HD") ~ "HD",
    TRUE ~ NA_character_  # Handle other cases if any
  )) %>%
  mutate(Responder = factor(Responder, levels = c("HD", "R", "NR"))) %>%
  #filter(Checkpoint == "anti-CTLA-4 + anti-PD-1") %>%
  mutate(Timepoint = ifelse(Responder == "HD", "HD",
                            ifelse(Visit_Name_time == "DAY_0", 0,
                                   ifelse(Visit_Name_time == "DAY_43", 43,
                                          ifelse(Visit_Name_time == "DAY_89/101", "89/101", NA))))) 

setwd("C:/Users/pmunson/OneDrive - Nurix Therapeutics, Inc/Desktop/Autoantibodies/Graphs") 
ggsave(filename = "pt_10_graph.png", plot = pt_10_graph, width = 4, height = 4, dpi = 600)

pt_10_graph <- pt_10_graphing %>%
  ggplot(., aes(x = Timepoint, y = value, group = Subject_Protein, color = Responder)) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "white") +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~Subject_ID) +
  theme_prism() +
  ylab("RPK") +
  xlab("") +
  theme(legend.position = "none") +
  scale_color_manual(values = c('#738678','#b03608'))  +
  geom_text_repel(data = subset(pt_10_graphing %>% filter(value >500 & Timepoint == "89/101")), 
                  aes(label = Protein),
                  direction = "x",
                  max.overlaps = 100, min.segment.length = 0.01, box.padding = 2, color = "black", size =4)

Phip_seq_final$Determinant_Spreading

Phip_seq_final %>%
  filter(!is.na(Determinant_Spreading)) %>%
  filter(Lysate == "Round3") %>%
  #filter(Visit_Name_time == "DAY_0") %>%
  filter(Cancer == "Melanoma") %>%
  filter(Subject_ID %in% subjects_that_have_d43) %>%
  #filter(Visit_Name_time == "DAY_0" | Visit_Name_time == "DAY_43") %>%
  #filter(Protein %in% up_in_melanoma_phip) %>%
  filter(Protein %in% cta_list | Protein %in% induced_proteins$Protein) %>%
  # filter(Protein == "MAGEA10") %>%
  mutate(Subject_Protein = str_c(Subject_ID, Protein, sep = "_")) %>%
  mutate(Responder = case_when(
    Best_Clinical_outcome %in% c("PR", "SD", "NED1") ~ "R",
    Best_Clinical_outcome %in% c("NED2", "PD") ~ "NR",
    Best_Clinical_outcome %in% c("HD") ~ "HD",
    TRUE ~ NA_character_  # Handle other cases if any
  )) %>%
  mutate(Responder = factor(Responder, levels = c("HD", "R", "NR"))) %>%
  #filter(Checkpoint == "anti-CTLA-4 + anti-PD-1") %>%
  ggplot(., aes(x = Visit_Name_time, y = value, group = Subject_Protein, color =Responder)) +
  #geom_hline(yintercept = 250, linetype = "dashed", color = "red") +
  #geom_hline(yintercept = 50, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "white") +
  #geom_hline(yintercept = 300, linetype = "dashed", color = "white") +
  geom_line() +
  geom_point() +
  facet_wrap(~Subject_ID) +
  theme_prism() +
  ylab("RPK") +
  xlab("") 
  geom_text_repel(data = subset(Phip_seq_final %>%
                                  filter(value >500) %>%
                                  filter(Lysate == "Round3") %>%
                                  #filter(Visit_Name_time == "DAY_0") %>%
                                  filter(Cancer == "Melanoma") %>%
                                  filter(Subject_ID %in% subjects_that_have_d43) %>%
                                  filter(Visit_Name_time == "DAY_89/101") %>%
                                  filter(Subject_ID == "PT_7") %>%
                                  #filter(Protein %in% up_in_melanoma_phip) %>%
                                  filter(Protein %in% induced_proteins$Protein)%>%
                                  mutate(Subject_Protein = str_c(Subject_ID, Protein, sep = "_"))), 
                  aes(label = Protein),
                  size = 3,
                  direction = "both",
                  max.overlaps = 100, min.segment.length = 0.01, box.padding = 0.25, color = "black")
  
  
  #### Identifying induced proteins
  
detected_in_hd <- Phip_seq_final %>%
    filter(Lysate == "Round3") %>%
    filter(Cancer != "Melanoma") %>%
    filter(value >100) %>%
    select(Protein) %>% unique() %>% pull()
  
  
detected_in_melanoma_at_baseline <- Phip_seq_final %>%
    filter(Lysate == "Round3") %>%
    filter(Cancer == "Melanoma") %>%
    filter(Visit_Name_time == "DAY_0") %>%
    mutate(Subject_Protein = str_c(Subject_ID, Protein, sep = "_")) %>%
    filter(value >50) %>%
    select(Subject_Protein) %>% unique() %>% pull()

detected_in_melanoma_at_post_vax <- Phip_seq_final %>%
  filter(Lysate == "Round3") %>%
  filter(Cancer == "Melanoma") %>%
  #filter(!(Protein %in%  detected_in_hd)) %>%
  filter(Visit_Name_time == "DAY_43" | Visit_Name_time == "DAY_89/101") %>%
  mutate(Subject_Protein = str_c(Subject_ID, Protein, sep = "_")) %>%
  filter(value >250) %>%
  select(Subject_Protein) %>% unique() %>% pull()

induced_protein_v2 <- setdiff(detected_in_melanoma_at_post_vax,detected_in_melanoma_at_baseline) 


Phip_seq_final %>%
  mutate(Subject_Protein = str_c(Subject_ID, Protein, sep = "_")) %>%
  filter(Subject_Protein %in% induced_protein_v2) %>%
  filter(Lysate == "Round3") %>%
  #filter(Visit_Name_time == "DAY_0") %>%
  filter(Cancer == "Melanoma") %>%
  filter(Subject_ID %in% subjects_that_have_d43) %>%
  #filter(Visit_Name_time == "DAY_0" | Visit_Name_time == "DAY_43") %>%
  #filter(Protein %in% up_in_melanoma_phip) %>%
  # filter(Protein == "MAGEA10") %>%
  mutate(Responder = case_when(
    Best_Clinical_outcome %in% c("PR", "SD", "NED1") ~ "R",
    Best_Clinical_outcome %in% c("NED2", "PD") ~ "NR",
    Best_Clinical_outcome %in% c("HD") ~ "HD",
    TRUE ~ NA_character_  # Handle other cases if any
  )) %>%
  mutate(Responder = factor(Responder, levels = c("HD", "R", "NR"))) %>%
  #filter(Checkpoint == "anti-CTLA-4 + anti-PD-1") %>%
  ggplot(., aes(x = Visit_Name_time, y = value, group = Subject_Protein, color =Responder)) +
  #geom_hline(yintercept = 250, linetype = "dashed", color = "red") +
  #geom_hline(yintercept = 50, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "white") +
  #geom_hline(yintercept = 300, linetype = "dashed", color = "white") +
  geom_line() +
  geom_point() +
  facet_wrap(~Subject_ID, scales = "free") +
  theme_prism() +
  ylab("RPK") +
  xlab("")
