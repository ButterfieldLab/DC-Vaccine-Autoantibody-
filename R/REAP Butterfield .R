## Analyzing the REAP data sent by Yile for the melanoma vaccine patients

#Libraries 
{
required_packages <- c(
  "readr", "dplyr", "ggplot2", "forcats", "data.table",
  "gridExtra", "gridtext", "egg", "ggrepel",
  "corrplot", "pROC", "randomForest", "tidyr", "readxl",
  "umap", "Rtsne", "ggthemes", "ggpubr", "patchwork",
  "cowplot", "Polychrome", "ggalt", "writexl",
  "pheatmap", "heatmap3"
)

# Function to install and load packages
install_and_load <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, dependencies = TRUE)
  }
  suppressMessages(library(package, character.only = TRUE))
}


# Apply the function to install and load each package
invisible(sapply(required_packages, install_and_load))
}


#Importing the data

Paul_score <- read.csv("C:/Users/pmunson/OneDrive - Nurix Therapeutics, Inc/Desktop/Autoantibodies/REAP_cancer_vaccine_butterfield.csv" , row.names = 1)
Melanoma_Vaccine_Clinical_Data_Butterfield_Lab <- read_excel("C:/Users/pmunson/OneDrive - Nurix Therapeutics, Inc/Desktop/Autoantibodies/cancer vaccine Clinical Metadata.xlsx")

# I wanted to separate the patient number from the time point for easier grouping 
# This is in the data frame Paul_patient_info


Paul_transpose <- as.data.frame(t(as.matrix(Paul_score)))
Patient_numbers <- gsub("^([^_]*_[^_]*)_.*$", "\\1", rownames(Paul_transpose))
Patient_timepoints <- gsub("^[^_]*_[^_]*_", "\\1", rownames(Paul_transpose))

#This is to make the healthy donors time point be called a "baseline" time point

Patient_numbers[87:91] <- c("HD1", "HD2", "HD3", "HD4", "HD5")
Patient_timepoints[87:91] <- c("DAY_0","DAY_0","DAY_0","DAY_0","DAY_0")

Paul_patient_info <- cbind(Patient_numbers, Patient_timepoints, Paul_transpose)


## Perform a tSNE on the data and looking at different parameters
set.seed(1)
tsne <- Rtsne(Paul_patient_info[,-c(1:2)], dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)
tsne_patient_data <- as.data.frame(cbind(Patient_numbers, Patient_timepoints, tsne$Y))

tsne_patient_data[,3] <- as.numeric(tsne_patient_data[,"V3"])
tsne_patient_data[,4] <- as.numeric(tsne_patient_data[,"V4"])

## Converting clinical data so patient # is consistent

Melanoma_Vaccine_Clinical_Data_Butterfield_Lab<- Melanoma_Vaccine_Clinical_Data_Butterfield_Lab %>%
  mutate(patient_numbers =  toupper(patient_numbers)) %>%
  rename(Patient_numbers = patient_numbers)

tsne_patient_data <- left_join(tsne_patient_data, Melanoma_Vaccine_Clinical_Data_Butterfield_Lab)


### I want to specify the order of the clinical outcomes

tsne_patient_data$Best_Clinical_outcome <- factor(tsne_patient_data$Best_Clinical_outcome, levels = c("HD", "PR", "SD", "NED1", "NED2", "PD"))


## By only doing rows 1-86 I can exclude the healthy donors
## Note, I am considering all of the healthy donors to be "DAY_0"

### Overview of patients - OS vs. PFS

ggplot(Melanoma_Vaccine_Clinical_Data_Butterfield_Lab[1:35,], aes(x = log2(PFS), y = log2(Overall_Survival), color = Best_Clinical_outcome)) + 
ggalt::geom_encircle(aes(group=Best_Clinical_outcome, spread = 0.0002)) +
  geom_point(size = 4)+
  scale_color_manual(values = c("#00AC8C",  "#FF6C2F", "#F65058" , "#685BC7", "#0083C3")) +
  xlab("PFS Log2(Months)") + 
  ylab("OS Log2(Months)") + 
  xlim(0.5,5.8) +
  ylim(1,5.8) +
  theme(legend.position = "none") +
  theme_par() +
  theme(axis.line = element_line(colour = 'black', size = 0.2)) +
  ggtitle("Clinical Outcome") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = c(0.9,0.5)) +
  guides(shape = FALSE) +
  guides(size = FALSE) +
  guides(col=guide_legend("RECIST Criteria")) + 
  theme(text=element_text(size=12,  family="Arial")) + 
  theme(plot.margin=unit(c(0,0,0,0), "cm"))

##########
##### Reorganizing tSNES's
{

fig1a <- ggplot(tsne_patient_data[1:86,], aes(x = V3, y = V4, color = Gender)) + 
  ggalt::geom_encircle(aes(group=Patient_numbers, spread = 0.0002)) +
  geom_point(alpha = 1, aes(size =4, shape = as.factor(Patient_timepoints)))+
  scale_color_manual(values = c("#F65058","#0083C3")) +
  xlab("") + 
  ylab("") + 
  xlim(-7.5,13) +
  ylim(-8,8.5) +
  theme(legend.position = "none") +
  theme_par() +
  theme(axis.line = element_line(colour = 'black', size = 0.2)) +
  ggtitle("Gender") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = c(0.9,0.5)) +
  guides(shape = FALSE) +
  guides(size = FALSE) +
  guides(col=guide_legend("Gender")) +
  theme(text=element_text(size=12,  family="Arial")) + 
  theme(plot.margin=unit(c(0,0,0,0), "cm")) 

fig1b <- ggplot(tsne_patient_data[1:86,], aes(x = V3, y = V4, color = Age)) + 
  geom_point(alpha = 1, aes(size =4, shape = as.factor(Patient_timepoints)))+
  ggalt::geom_encircle(aes(group=Patient_numbers, spread = 0.0002)) +
  scale_colour_gradient(low = "black", high = "#F65058") +
  xlab("") + 
  ylab("") + 
  xlim(-7.5,13) +
  ylim(-8,8.5) +
  theme_par() +
  theme(axis.line = element_line(colour = 'black', size = 0.2)) +
  ggtitle("Age")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = c(0.9,0.5)) +
  guides(shape = FALSE) +
  guides(size = FALSE) +
  guides(col=guide_legend("Age (Years)")) +
  theme(text=element_text(size=12,  family="Arial")) + 
  theme(plot.margin=unit(c(0,0,0,0), "cm"))

### Basic parameters gender & age

grid.arrange(fig1a, fig1b)

## Tumor characteristics, location, mutation 

fig1c <- ggplot(tsne_patient_data[1:86,], aes(x = V3, y = V4, color =  Metastasis_Site)) + 
  ggalt::geom_encircle(aes(group=Patient_numbers, spread = 0.0002)) +
  geom_point(alpha = 1, aes(size =4, shape = as.factor(Patient_timepoints)))+
  xlab("") + 
  ylab("") + 
  xlim(-7.5,13) +
  ylim(-8,8.5) +
  theme_par() +
  theme(axis.line = element_line(colour = 'black', size = 0.2)) +
  ggtitle("Metastasis Site") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = c(0.89,0.5)) +
  guides(shape = FALSE) +
  guides(size = FALSE) +
  guides(col=guide_legend(ncol=1)) +
  scale_colour_brewer(palette = "Paired", name = "Metastasis Site")+
  theme(text=element_text(size=12,  family="Arial")) + 
  theme(plot.margin=unit(c(0,0,0,0), "cm"))

fig1d <- ggplot(tsne_patient_data[1:86,], aes(x = V3, y = V4, color = as.factor(BRAF_Mutation))) + 
  geom_point(alpha = 1, aes(size =4, shape = as.factor(Patient_timepoints)))+
  ggalt::geom_encircle(aes(group=Patient_numbers, spread = 0.0002)) +
  scale_color_manual(values = c("#0083C3","#F65058","black"), labels = c("No", "Yes", "NA")) +
  xlab("") + 
  ylab("") + 
  xlim(-7.5,13) +
  ylim(-8,8.5) +
  theme_par() +
  theme(axis.line = element_line(colour = 'black', size = 0.2)) +
  ggtitle("BRAF Mutation")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = c(0.9,0.5)) +
  guides(shape = FALSE) +
  guides(size = FALSE) +
  guides(col=guide_legend("BRAF Mutation")) +
  theme(text=element_text(size=12,  family="Arial")) + 
  theme(plot.margin=unit(c(0,0,0,0), "cm"))

fig1e <- ggplot(tsne_patient_data[1:86,], aes(x = V3, y = V4, color = as.factor(NRAS_Mutation))) + 
  geom_point(alpha = 1, aes(size =4, shape = as.factor(Patient_timepoints)))+
  ggalt::geom_encircle(aes(group=Patient_numbers, spread = 0.0002)) +
  scale_color_manual(values = c("#0083C3","#F65058","black"), labels = c("No", "Yes", "NA")) +
  xlab("") + 
  ylab("") + 
  xlim(-7.5,13) +
  ylim(-8,8.5) +
  theme_par() +
  theme(axis.line = element_line(colour = 'black', size = 0.2)) +
  ggtitle("NRAS Mutation")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = c(0.9,0.5)) +
  guides(shape = FALSE) +
  guides(size = FALSE) +
  guides(col=guide_legend("NRAS Mutation")) +
  theme(text=element_text(size=12,  family="Arial")) + 
  theme(plot.margin=unit(c(0,0,0,0), "cm"))

## Metastasis, BRAF, NRAS

grid.arrange(fig1c, fig1d, fig1e)

## Prior Immunotherapy 

fig1f <- ggplot(tsne_patient_data[1:86,], aes(x = V3, y = V4, color = as.factor(Prior_IL2))) + 
  geom_point(alpha = 1, aes(size =4, shape = as.factor(Patient_timepoints)))+
  ggalt::geom_encircle(aes(group=Patient_numbers, spread = 0.0002)) +
  scale_color_manual(values = c("black", "#F65058"), labels = c("No", "Yes")) +
  xlab("") + 
  ylab("") + 
  xlim(-7.5,13) +
  ylim(-8,8.5) +
  theme_par() +
  theme(axis.line = element_line(colour = 'black', size = 0.2)) +
  ggtitle("Prior IL-2")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = c(0.9,0.5)) +
  guides(shape = FALSE) +
  guides(size = FALSE) +
  guides(col=guide_legend("IL-2")) +
  theme(text=element_text(size=12,  family="Arial")) + 
  theme(plot.margin=unit(c(0,0,0,0), "cm")) 

fig1g <- ggplot(tsne_patient_data[1:86,], aes(x = V3, y = V4, color = as.factor(Prior_IFN))) + 
  geom_point(alpha = 1, aes(size =4, shape = as.factor(Patient_timepoints)))+
  ggalt::geom_encircle(aes(group=Patient_numbers, spread = 0.0002)) +
  scale_color_manual(values = c("black", "#F65058"), labels = c("No", "Yes")) +
  xlab("") + 
  ylab("") + 
  xlim(-7.5,13) +
  ylim(-8,8.5) +
  theme_par() +
  theme(axis.line = element_line(colour = 'black', size = 0.2)) +
  ggtitle("Prior IFN")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = c(0.9,0.5)) +
  guides(shape = FALSE) +
  guides(size = FALSE) +
  guides(col=guide_legend("IFN")) +
  theme(text=element_text(size=12,  family="Arial")) + 
  theme(plot.margin=unit(c(0,0,0,0), "cm")) 

fig1h <- ggplot(tsne_patient_data[1:86,], aes(x = V3, y = V4, color = as.factor(Anti_PD_1))) + 
  geom_point(alpha = 1, aes(size =4, shape = as.factor(Patient_timepoints)))+
  ggalt::geom_encircle(aes(group=Patient_numbers, spread = 0.0002)) +
  scale_color_manual(values = c("black", "#F65058"), labels = c("No", "Yes")) +
  xlab("") + 
  ylab("") + 
  xlim(-7.5,13) +
  ylim(-8,8.5) +
  theme_par() +
  theme(axis.line = element_line(colour = 'black', size = 0.2)) +
  ggtitle("Prior Anti-PD-1")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = c(0.9,0.5)) +
  guides(shape = FALSE) +
  guides(size = FALSE) +
  guides(col=guide_legend("Anti-PD-1")) +
  theme(text=element_text(size=12,  family="Arial")) + 
  theme(plot.margin=unit(c(0,0,0,0), "cm"))

fig1i <- ggplot(tsne_patient_data[1:86,], aes(x = V3, y = V4, color = as.factor(Anti_CTLA_4))) + 
  geom_point(alpha = 1, aes(size =4, shape = as.factor(Patient_timepoints)))+
  ggalt::geom_encircle(aes(group=Patient_numbers, spread = 0.0002)) +
  scale_color_manual(values = c("black", "#F65058"), labels = c("No", "Yes")) +
  xlab("") + 
  ylab("") + 
  xlim(-7.5,13) +
  ylim(-8,8.5) +
  theme(legend.position = "none") +
  theme_par() +
  theme(axis.line = element_line(colour = 'black', size = 0.2)) +
  ggtitle("Prior Anti-CTLA-4")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = c(0.9,0.5)) +
  guides(shape = FALSE) +
  guides(size = FALSE) +
  guides(col=guide_legend("Anti-CTLA-4")) +
  theme(text=element_text(size=12,  family="Arial")) + 
  theme(plot.margin=unit(c(0,0,0,0), "cm"))

fig1j <- ggplot(tsne_patient_data[1:86,], aes(x = V3, y = V4, color = Checkpoint)) + 
  geom_point(alpha = 1, aes(size =4, shape = as.factor(Patient_timepoints)))+
  ggalt::geom_encircle(aes(group=Patient_numbers, spread = 0.0002)) +
  scale_color_manual(values = c("#0083C3","#F65058","black")) +
  xlab("") + 
  ylab("") + 
  xlim(-7.5,13) +
  ylim(-8,8.5) +
  theme_par() +
  theme(axis.line = element_line(colour = 'black', size = 0.2)) +
  ggtitle("Prior Checkpoint")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = c(0.9,0.5)) +
  guides(shape = FALSE) +
  guides(size = FALSE) +
  guides(col=guide_legend("Checkpoint")) +
  theme(text=element_text(size=12,  family="Arial")) + 
  theme(plot.margin=unit(c(0,0,0,0), "cm"))

grid.arrange(fig1f, fig1g, fig1h, fig1i, fig1j)

### Now I want to do the clinical outcome, PFS and OS

fig1k <- ggplot(tsne_patient_data, aes(x = V3, y = V4, color = Best_Clinical_outcome)) + 
  ggalt::geom_encircle(aes(group=Patient_numbers, spread = 0.0002)) +
  geom_point(alpha = 1, aes(size =4, shape = as.factor(Patient_timepoints)))+
  scale_color_manual(values = c("black", "#685BC7",  "#0083C3", "#00AC8C" , "#FF6C2F", "#F65058")) +
  xlab("") + 
  ylab("") + 
  xlim(-7.5,13) +
  ylim(-8,8.5) +
  theme(legend.position = "none") +
  theme_par() +
  theme(axis.line = element_line(colour = 'black', size = 0.2)) +
  ggtitle("Clinical Outcome") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = c(0.9,0.5)) +
  guides(shape = FALSE) +
  guides(size = FALSE) +
  guides(col=guide_legend("Clinical outcome")) + 
  theme(text=element_text(size=12,  family="Arial")) + 
  theme(plot.margin=unit(c(0,0,0,0), "cm"))

fig1l <- ggplot(tsne_patient_data[1:86,], aes(x = V3, y = V4, color = log2(PFS))) + 
  geom_point(alpha = 1, aes(size =4, shape = as.factor(Patient_timepoints)))+
  ggalt::geom_encircle(aes(group=Patient_numbers, spread = 0.0002)) +
  scale_colour_gradient(low = "black", high = "#F65058") +
  xlab("") + 
  ylab("") + 
  xlim(-7.5,13) +
  ylim(-8,8.5) +
  theme_par() +
  theme(axis.line = element_line(colour = 'black', size = 0.2)) +
  ggtitle("Progression Free Survival")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = c(0.9,0.5)) +
  guides(shape = FALSE) +
  guides(size = FALSE) +
  guides(col=guide_legend("PFS Log2(Months)")) +
  theme(text=element_text(size=12,  family="Arial")) + 
  theme(plot.margin=unit(c(0,0,0,0), "cm"))

fig1m <- ggplot(tsne_patient_data[1:86,], aes(x = V3, y = V4, color = Overall_Survival)) + 
  geom_point(alpha = 1, aes(size =4, shape = as.factor(Patient_timepoints)))+
  ggalt::geom_encircle(aes(group=Patient_numbers, spread = 0.0002)) +
  scale_colour_gradient(low = "black", high = "#F65058") +
  xlab("") + 
  ylab("") + 
  xlim(-7.5,13) +
  ylim(-8,8.5) +
  theme_par() +
  theme(axis.line = element_line(colour = 'black', size = 0.2)) +
  ggtitle("Overall Survival")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = c(0.9,0.5)) +
  guides(shape = FALSE) +
  guides(size = FALSE) +
  guides(col=guide_legend("OS (Months)")) +
  theme(text=element_text(size=12,  family="Arial")) + 
  theme(plot.margin=unit(c(0,0,0,0), "cm"))

tsne_patient_data %>%
  filter(Determinant_Spreading != "NA" & Cancer == "Yes")%>%
ggplot(., aes(x = V3, y = V4, color = as.factor(Determinant_Spreading))) + 
  geom_point(alpha = 1, aes(size =4, shape = as.factor(Patient_timepoints)))+
  ggalt::geom_encircle(aes(group=Patient_numbers, spread = 0.0002)) +
  #scale_colour_gradient(low = "black", high = "#F65058") +
  xlab("") + 
  ylab("") + 
  xlim(-7.5,13) +
  ylim(-8,8.5) +
  theme_par() +
  theme(axis.line = element_line(colour = 'black', size = 0.2)) +
  ggtitle("Determinant Spreading")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = c(0.9,0.5)) +
  guides(shape = FALSE) +
  guides(size = FALSE) +
  guides(col=guide_legend("Determinant Spreading)")) +
  theme(text=element_text(size=12,  family="Arial")) + 
  theme(plot.margin=unit(c(0,0,0,0), "cm"))

tsne_patient_data %>%
  filter(Vaccine_Response != "NA" & Cancer == "Yes")%>%
  ggplot(., aes(x = V3, y = V4, color = as.factor(Vaccine_Response))) + 
  geom_point(alpha = 1, aes(size =4, shape = as.factor(Patient_timepoints)))+
  ggalt::geom_encircle(aes(group=Patient_numbers, spread = 0.0002)) +
  #scale_colour_gradient(low = "black", high = "#F65058") +
  xlab("") + 
  ylab("") + 
  xlim(-7.5,13) +
  ylim(-8,8.5) +
  theme_par() +
  theme(axis.line = element_line(colour = 'black', size = 0.2)) +
  ggtitle("Vaccine Response")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = c(0.9,0.5)) +
  guides(shape = FALSE) +
  guides(size = FALSE) +
  guides(col=guide_legend("Vaccine Response)")) +
  theme(text=element_text(size=12,  family="Arial")) + 
  theme(plot.margin=unit(c(0,0,0,0), "cm"))


grid.arrange(fig1k, fig1l, fig1m)

#### Looking at Correlation with tSNE1
### exploratory plot, doesn't seem interesting
ggplot(tsne_patient_data[1:86,], aes(x = V3, y = Overall_Survival)) + 
  geom_point(alpha = 1, aes(size =4, shape = as.factor(Patient_timepoints)))+
  xlab("tSNE1") + 
  ylab("Overall Survival (Months)") + 
  #xlim(-7.5,13) +
  #ylim(-8,8.5) +
  theme_par() +
  theme(axis.line = element_line(colour = 'black', size = 0.2)) +
  ggtitle("Checkpoint vs. tSNE1")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = c(0.9,0.3)) +
  guides(shape = FALSE) +
  guides(size = FALSE) +
  guides(col=guide_legend("OS (Months)")) +
  theme(text=element_text(size=12,  family="Arial")) + 
  theme(plot.margin=unit(c(0,0,0,0), "cm"))

}

#############
### Now I want to start diving into the raw data 
## Initially, I wanted to calculate the total REAP scores over time by summing ALL of the REAP values

## I made a tidy version of the data to generate the individual value data 
Paul_patient_tidy <- Paul_patient_info %>%
                        gather(key = "Protein", value = "REAP", 3:6186)



Total_REAP <- Paul_patient_tidy %>%
  group_by(Patient_numbers, Patient_timepoints) %>%
  summarise(Total_REAP = sum(REAP))%>%
  as.data.frame()



## Adding clinical data to Total REAP
Total_REAP <-left_join(Total_REAP, Melanoma_Vaccine_Clinical_Data_Butterfield_Lab)

## exploratory, just seeing if the total reap signal differes among patients 
Total_REAP %>%
  ggplot(aes(x = Best_Clinical_outcome, y = log2(Total_REAP), color = Best_Clinical_outcome, group = Patient_numbers))+
  geom_point(size = 3) + 
  scale_color_manual(values = c("black", "#00AC8C",  "#FF6C2F", "#F65058" , "#685BC7", "#0083C3")) +
  theme_par() +
  xlab("Clinical Outcome") +
  ylab("Log2 (Sum of REAP)") +
  ylim(8,9)

## Creating a tidy version of the data 

### Experimenting on a single patient
Paul_patient_tidy %>%
  filter(Patient_numbers == "PT_20")%>%
    ggplot(., aes(x = Patient_timepoints, y = log2(REAP))) +
      geom_point(alpha = 0.4) +
      theme_par() +
      geom_jitter()+
      xlab("") 

## Defining how I want the patients to be ordered
Paul_patient_tidy$Patient_numbers <- factor(Paul_patient_tidy$Patient_numbers, levels = c("PT_1", "PT_2", "PT_3", "PT_4", "PT_5", "PT_6", "PT_7", "PT_8", "PT_9", "PT_10", "PT_11", "PT_12", "PT_13", "PT_14", "PT_15", "PT_16", "PT_17", "PT_18", "PT_19", "PT_20", "PT_21", "PT_22", "PT_23", "PT_24", "PT_25", "PT_26", "PT_27", "PT_28", "PT_29", "PT_30", "PT_31", "PT_32", "PT_33", "PT_34", "PT_35", "HD1", "HD2", "HD3", "HD4","HD5"))
                                            
                                                                                          
## Showing the distribution of REAP values over a threshold    
## REAP signal per protein among all patients
Paul_patient_tidy %>%
  filter(REAP >= 1.5)%>%
  ggplot(., aes(x = Patient_timepoints, y = REAP)) +
  geom_point(alpha = 0.4) +
  theme_par() +
  geom_jitter()+
  xlab("") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", color="#0083C3", width=0.2, size = 1) +
  stat_summary(fun=mean, geom="point", color="#0083C3", size = 2)+
  facet_wrap(~Patient_numbers) 


## I want to try and count the number of hits above a certain threshold
## This is where I exported the total number of hits

Number_of_hits <- Paul_patient_tidy %>%
                      group_by(Patient_numbers, Patient_timepoints) %>%
                          dplyr::summarize(n_any_hits = sum(REAP>0),
                                     n_moderate_hits = sum(REAP >1),     
                                     n_high_hits = sum(REAP>=6.5)) %>%
                          as.data.frame()

Total_REAP_hits <- left_join(Total_REAP, Number_of_hits)

Total_REAP_hits$Best_Clinical_outcome <- factor(Total_REAP_hits$Best_Clinical_outcome, levels = c("HD", "PR", "SD", "NED1", "NED2", "PD"))


### Code showing the statisttical difference between men and women 
### Also adding in difference between healthy and cancer patients
### I think I will also want to specifiy specifically "immunne proteins" but I will have to add this in later 

set.seed(1)

Total_REAP_hits %>%
  filter(Patient_timepoints == "DAY_0", Cancer == "Yes")%>%
  ggplot(aes(x = Gender, y = n_high_hits))+
  geom_point(size = 2.5, position = position_jitter(w =0.15,h= 0.01)) + 
  #scale_color_manual(values = c("#685BC7",  "#0083C3", "#00AC8C" , "#FF6C2F", "#F65058")) +
  theme_par() +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", color="#F65058", width=0.2, size = 1) +
  stat_summary(fun=mean, geom="point", color="#F65058", size = 3)+
  xlab("") +
  ylab("Number of REAP >6.5") + 
  ylim(-0.5,10) +
  ggtitle("Sex")+
  stat_compare_means(method = "t.test", label = "p.format", label.y = 9, label.x = 1.5, comparisons = list(c("F","M")))

set.seed(1)

Total_REAP_hits %>%
  mutate(Cancer = ifelse(is.na(Cancer), "None", Cancer)) %>%
  filter(Patient_timepoints == "DAY_0")%>%
  ggplot(aes(x = Cancer, y = n_high_hits))+
  geom_point(size = 2.5, position = position_jitter(w =0.15,h= 0.01)) + 
  #scale_color_manual(values = c("#685BC7",  "#0083C3", "#00AC8C" , "#FF6C2F", "#F65058")) +
  theme_par() +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", color="#F65058", width=0.2, size = 1) +
  stat_summary(fun=mean, geom="point", color="#F65058", size = 3)+
  xlab("") +
  ylab("Number of REAP >6.5") + 
  ylim(-0.5,10) +
  ggtitle("Cancer")+
  stat_compare_means(method = "t.test", label = "p.format", label.y = 9, label.x = 1.5, comparisons = list(c("None","Yes")))

set.seed(1)
Total_REAP_hits %>%
  mutate(Cancer = ifelse(is.na(Cancer), "None", Cancer)) %>%
  filter(Patient_timepoints == "DAY_0")%>%
  ggplot(aes(x = Cancer, y = n_moderate_hits))+
  geom_point(size = 2.5, position = position_jitter(w =0.15,h= 0.01)) + 
  #scale_color_manual(values = c("#685BC7",  "#0083C3", "#00AC8C" , "#FF6C2F", "#F65058")) +
  theme_par() +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", color="#F65058", width=0.2, size = 1) +
  stat_summary(fun=mean, geom="point", color="#F65058", size = 3)+
  xlab("") +
  ylab("Number of REAP >1") + 
  ylim(0,120) +
  ggtitle("Cancer")+
  stat_compare_means(method = "t.test", label = "p.format", label.y = 9, label.x = 1.5, comparisons = list(c("None","Yes")))
 
## Adding statistics to compare different groups based on number of hits

t1 <- Total_REAP_hits %>%
  filter(Patient_timepoints == "DAY_0", Cancer == "Yes", Gender == "F") %>%
    select(n_high_hits)
t2 <- Total_REAP_hits %>%
  filter(Patient_timepoints == "DAY_0", Cancer == "Yes", Gender == "M") %>%
  select(n_high_hits)

t_stat <- t.test(t1, t2 , paired=FALSE)
t_stat$p.value


Number_of_hits_day0 <- Number_of_hits %>% filter(Patient_timepoints == "DAY_0")

getwd()
#write_xlsx(Number_of_hits_day0,"Number_of_hits_day0_moderate.xlsx")

## This is a way to explore the number of hits above a given threshold
hits <- Paul_patient_tidy %>%
  filter(REAP == 0 & Patient_timepoints == "DAY_0") %>%
    group_by(Protein) %>%
    count()

Paul_patient_tidy_number_blank <- left_join(Paul_patient_tidy, as.data.frame(hits))
hist(Paul_patient_tidy_number_blank$n)

Paul_patient_tidy_number_blank %>%
  filter(n >= 32 & n < 40 & REAP > 0.5)

Proteins_of_interest <- Paul_patient_tidy_number_blank %>%
                          filter(n >= 30 & n < 40 & REAP > 0.5) 

unique(Proteins_of_interest$Protein)

Paul_patient_tidy %>%
  filter(Protein == "CTLA4") %>%
  ggplot(.,aes(x = Patient_numbers, y = REAP)) + 
  geom_point() +
  ylim(0,9) +
  theme_par() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))

## Current work
### Creating Paul_patient_info_REAP
Paul_patient_info_REAP <- left_join(Melanoma_Vaccine_Clinical_Data_Butterfield_Lab, Paul_patient_info)

Paul_patient_info_REAP <-Paul_patient_info_REAP %>%
  filter (Patient_timepoints == "DAY_0")

ggplot(Paul_patient_info_REAP, aes(x = PDCD1, y = CTLA4, color = Checkpoint, shape = as.factor(Patient_timepoints))) + 
  geom_point(position = position_jitter(w =0.1,h= 0.1)) + 
  scale_color_manual(values = c("#0083C3","#F65058","black")) +
  xlab("") + 
  ylab("") + 
  xlim(-1,7) +
  ylim(-1,6) +
  theme_par() +
  theme(axis.line = element_line(colour = 'black', size = 0.2)) +
  ggtitle("Prior Checkpoint")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = c(0.9,0.5)) +
  guides(shape = FALSE) +
  guides(size = FALSE) +
  guides(col=guide_legend("Checkpoint")) +
  theme(text=element_text(size=12)) + 
  theme(plot.margin=unit(c(0,0,0,0), "cm"))

### Cool protein list 
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

Tidy_REAP_Clinical <- left_join(Paul_patient_tidy, Melanoma_Vaccine_Clinical_Data_Butterfield_Lab)

fig3a <- Tidy_REAP_Clinical %>%
  #filter(Patient_timepoints == "DAY_0") %>%
filter(Protein == "PDCD1" & Cancer == "Yes")%>%
  ggplot(aes(x = as.factor(Anti_PD_1), y = REAP, color = Checkpoint))+

  #scale_color_manual(values = c("black", "#00AC8C",  "#FF6C2F", "#F65058" , "#685BC7", "#0083C3")) +
  theme_par() +
  xlab("") +
  ylab("Anti-PD-1 (REAP)") +
  ylim(-2,8) +
  ggtitle("Prior Nivolumab or Pembrolizumab") + 
    geom_point(size = 3, position = position_jitter(w =0.2,h= 0)) + 
  facet_wrap(~Patient_timepoints)

fig3b <- Tidy_REAP_Clinical %>%
  #filter(Patient_timepoints == "DAY_0") %>%
  filter(Protein == "CTLA4" & Cancer == "Yes")%>%
  ggplot(aes(x = as.factor(Anti_CTLA_4), y = REAP, color = Checkpoint))+
  geom_point(size = 3, position = position_jitter(w =0.2,h= 0)) + 
  #scale_color_manual(values = c("black", "#00AC8C",  "#FF6C2F", "#F65058" , "#685BC7", "#0083C3")) +
  theme_par() +
  xlab("") +
  ylab("Anti-CTLA-4 (REAP)") +
  ylim(-2,8) +
  ggtitle("Prior Ipilimumab") + 
  facet_wrap(~Patient_timepoints)

grid.arrange(fig3a, fig3b, nrow = 1)

########

## Graph showing patients with previous checkpoint have a stastically lower tSNE1 value

######


###########################################

## Correlation of data with Overall Survival

##

##########################################

Paul_patient_info_REAP <- left_join(Melanoma_Vaccine_Clinical_Data_Butterfield_Lab, Paul_patient_info)

Paul_patient_info_REAP <-Paul_patient_info_REAP %>%
                              filter (Patient_timepoints == "DAY_0")

REAP_correlations <- correlation_table(data=Paul_patient_info_REAP , target="Overall_Survival")
REAP_correlations

Paul_patient_tidy %>%
  filter(Protein == "TMEM149") %>%
  ggplot(.,aes(x = Patient_numbers, y = REAP)) + 
  geom_point() +
  ylim(0,9) +
  theme_par() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))

###########
#CXCL16


### Filtering on a defined threshold prior to tSNE

###








#### Volcano plots
### I needed to reformat the table so that each rown is a protein and each column is a patient
### So I don't over-represent patients with mulitple timepoints I decided to focus on baseline


a <-  Tidy_REAP_Clinical %>% 
  filter(Patient_timepoints == "DAY_0") %>%
  group_by(Prior_Checkpoint) %>%
  nest()



### Formatting the data for volcano plots

Patients_day0_df <- Paul_score[ , grepl( "DAY_0" , names(Paul_score) ) ]


## Identifying the number of categories for a give variable
unique(Melanoma_Vaccine_Clinical_Data_Butterfield_Lab$Prior_Checkpoint)

## Finding the patiends that either have, or do not have, a given clinical status
No_Checkpoint_DF <- Patients_day0_df[,Melanoma_Vaccine_Clinical_Data_Butterfield_Lab$Prior_Checkpoint[1:35] == unique(Melanoma_Vaccine_Clinical_Data_Butterfield_Lab$Prior_Checkpoint)[1]]
Checkpoint_DF    <- Patients_day0_df[,Melanoma_Vaccine_Clinical_Data_Butterfield_Lab$Prior_Checkpoint[1:35] == unique(Melanoma_Vaccine_Clinical_Data_Butterfield_Lab$Prior_Checkpoint)[2]]

#setwd("/Users/pmunson/Box/Paul Munson/Autoantibody Project")
#getwd()

## Exporting these data frames to explore some stats in Prism

#write_xlsx(No_Checkpoint_DF,"No_Checkpoint_DF.xlsx")
#write_xlsx(Checkpoint_DF,"Checkpoint_DF.xlsx")

## The volcano plot doesn't appear to yiled many significant hits - I think just because there are SO many analytes


################
##
##
##  Z score, SD, range of values
##
##
################
dim(Paul_patient_info)

scale(Paul_patient_info[,-c(1:2)])

means_df <- colMeans(Paul_patient_info[,-c(1:2)])
sd_df <- as.data.frame(apply(Paul_patient_info[,-c(1:2)],2,sd))
colnames(sd_df) <- "sd"
sd_mean_df <- sd_df/means_df 

### Finding outliers for each column 
FindOutliers <- function(data) {
  lowerq = quantile(data)[2]
  upperq = quantile(data)[4]
  iqr = upperq - lowerq 
  extreme.threshold.upper = (iqr * 3) + upperq
  extreme.threshold.lower = lowerq - (iqr * 3)
  result <- which(data > extreme.threshold.upper | data < extreme.threshold.lower)
}

## Calculating the location of outliers for each column 
pos <- lapply(Paul_patient_info[,-c(1:2)], FindOutliers)
outliers_df <- as.data.frame(lengths(pos))

colnames(sd_mean_df) <- "sd_mean"
colnames(outliers_df) <- "outliers"

df_sd_mean_outliers <- cbind(sd_mean_df, outliers_df)


ggplot(df_sd_mean_outliers, aes(x = sd_mean, y = outliers)) +
  geom_point() +
  theme_par() +
  xlab("SD / mean per Protein") +
  ylab("No. of Outliers per Protein")

df_sd_mean_outliers <- cbind(df_sd_mean_outliers, sd_df)
setDT(df_sd_mean_outliers, keep.rownames = TRUE)[]
colnames(df_sd_mean_outliers)[1] <- "protein"
colnames(df_sd_mean_outliers)


### Cool protein list 
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

reap_label_subset <- df_sd_mean_outliers %>%
                       filter(protein == "CTLA4" | protein == "PDCD1" | protein == "VEGFB" | protein == "IFNA6" | protein == "CCR5_Epitope-3-N" | protein == "CD164L2" | protein == "IFNA13" | protein == "IL6" | protein == "TNFAIP6" | protein == "SLC2A2" | protein == "CDON" | protein == "MMP2" | protein == "IL18RAP")

reap_label_subset_Yile_not_Paul <- df_sd_mean_outliers %>%
                        filter(protein == "IL15RA"| protein == "CCL11")

ggplot(df_sd_mean_outliers, aes(outliers)) +
  geom_density(adjust = 1, fill = "black", color = "black", alpha = 0.2) +
  theme_par() +
  xlab("Outliers") +
  ylab("Density") +
  annotate("rect", xmin=c(3), xmax=c(10), ymin=c(-0.001) , ymax=c(0.06), alpha=0.5, color="red", fill="pink") 


### I want to label where RCN2 lands  

artifact_proteins <- df_sd_mean_outliers %>%
                      filter(protein == "RCN2")

df_sd_mean_outliers <- as.data.frame(df_sd_mean_outliers)
glimpse(df_sd_mean_outliers)





proteins_of_interest <- df_sd_mean_outliers %>%
                             filter(outliers >= 1) %>%
                             filter(log10(sd) >= -1 & log10(sd) <= 0.4 & log10(sd_mean) >= 0.57 & log10(sd_mean) <= 1)


ggplot(df_sd_mean_outliers, aes(x = sd_mean, y = outliers)) +
  geom_point() +
  theme_par() +
  xlab("SD / mean per Protein") +
  ylab("No. of Outliers per Protein")+
  geom_text_repel(data = reap_label_subset, aes(label = reap_label_subset$protein, fontface = 'bold', color = "purple"), size = 6,  min.segment.length = unit(0, 'lines'), force = 90)
  theme(legend.position = "none") 



sum(proteins_of_interest$protein == "VEGFB")
dim(proteins_of_interest)
Paul_score_proteins_of_interest <- subset(Paul_score, rownames(Paul_score) %in% proteins_of_interest$protein)
#write.csv(Paul_score_proteins_of_interest,"Paul_score_proteins_of_interest.csv")

##########################
##
##
## I want to get my proteins of interest on a per patient basis
##
#############################
colnames(Tidy_REAP_Clinical)
proteins_of_interest$protein 

pt_all_proteins <- Tidy_REAP_Clinical %>%
                       filter(Patient_numbers == "PT_1") %>%
                       filter(Protein %in% proteins_of_interest$protein) %>%
                        filter(REAP >0) 

unique(pt_all_proteins$Protein)  




getwd()
setwd("/Users/pmunson/Box/Paul Munson/Autoantibody Project")

############################################################
##
## Creating a heat map of my curated immune REAP hits
##
##
############################################################

## Importing the list of anti-immune autoantibodies
Immune_autoantibodies <- read_excel("C:/Users/pmunson/OneDrive - Nurix Therapeutics, Inc/Desktop/Autoantibodies/Immune_autoantibodies.xlsx")





Tidy_Immune_autoantibodies  <- Tidy_REAP_Clinical %>%
                                  filter(Protein %in% Immune_autoantibodies$REAP_name)  

REAP_Immune_autoantibodies <- subset(Paul_score, rownames(Paul_score) %in% Immune_autoantibodies$REAP_name)[c(Immune_autoantibodies$REAP_name),]
#write.csv(REAP_Immune_autoantibodies,"REAP_Immune_autoantibodies.csv")

immune_hits <- Paul_patient_tidy %>%
  filter(REAP >= 6.5 & Patient_timepoints == "DAY_0") %>%
  filter(Protein %in% Immune_autoantibodies$REAP_name) %>%
  group_by(Protein) %>%
  count()



Tidy_Immune_autoantibodies %>% 
  mutate(Cancer = ifelse(is.na(Cancer), "None", Cancer)) %>%
  filter(Protein %in% Immune_autoantibodies$REAP_name & Patient_timepoints == "DAY_0") %>%
  ggplot(aes(x = Cancer, y = REAP))+
  geom_point(size = 2.5, position = position_jitter(w =0.15,h= 0.01)) + 
  #scale_color_manual(values = c("#685BC7",  "#0083C3", "#00AC8C" , "#FF6C2F", "#F65058")) +
  theme_par() +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", color="#F65058", width=0.2, size = 1) +
  stat_summary(fun=mean, geom="point", color="#F65058", size = 3)+
  xlab("") +
  ylab("Number of REAP >6.5") + 
  ylim(-0.5,10) +
  ggtitle("Cancer")+
  stat_compare_means(method = "t.test", label = "p.format", label.y = 9, label.x = 1.5, comparisons = list(c("None","Yes")))
  


set.seed(1)
Total_REAP_hits %>%
  mutate(Cancer = ifelse(is.na(Cancer), "None", Cancer)) %>%
  filter(Patient_timepoints == "DAY_0")%>%
  ggplot(aes(x = Cancer, y = n_high_hits))+
  geom_point(size = 2.5, position = position_jitter(w =0.15,h= 0.01)) + 
  #scale_color_manual(values = c("#685BC7",  "#0083C3", "#00AC8C" , "#FF6C2F", "#F65058")) +
  theme_par() +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", color="#F65058", width=0.2, size = 1) +
  stat_summary(fun=mean, geom="point", color="#F65058", size = 3)+
  xlab("") +
  ylab("Number of REAP >6.5") + 
  ylim(-0.5,10) +
  ggtitle("Cancer")+
  stat_compare_means(method = "t.test", label = "p.format", label.y = 9, label.x = 1.5, comparisons = list(c("None","Yes")))




############################################################
##
## Making a heat map of my immune hits 
##
##
############################################################

### Playing with heatmap3

############################################################
##
## Heatmap 3 tutorial
##
##
############################################################

#gererate data
set.seed(123456789)
rnormData<-matrix(rnorm(1000), 40, 25)
rnormData[1:15, seq(6, 25, 2)] = rnormData[1:15, seq(6, 25, 2)] + 2
rnormData[16:40, seq(7, 25, 2)] = rnormData[16:40, seq(7, 25, 2)] + 4
colnames(rnormData)<-c(paste("Control", 1:5, sep = ""), 
                       paste(c("TrtA", "TrtB"), rep(1:10,each=2), sep = ""))

class(rnormData)

rownames(rnormData)<-paste("Probe", 1:40, sep = "")
ColSideColors<-cbind(Group1=c(rep("steelblue2",5), rep(c("brown1", 
                                                         "mediumpurple2"),10)),Group2=sample(c("steelblue2","brown1", 
                                                                                               "mediumpurple2"),25,replace=TRUE))
colorCell<-data.frame(row=c(1,3,5),col=c(2,4,6),color=c("green4",
                                                        "black","orange2"),stringsAsFactors=FALSE)
highlightCell<-data.frame(row=c(2,4,6),col=c(1,3,5),color=c("black",
                                                            "green4","orange2"),lwd=1:3,stringsAsFactors=FALSE)
#A simple example
heatmap3(rnormData,ColSideColors=ColSideColors,showRowDendro=FALSE,
         colorCell=colorCell,highlightCell=highlightCell)

#A more detail example
ColSideAnn<-data.frame(Information=rnorm(25),Group=c(rep("Control",5),
                                                     rep(c("TrtA","TrtB"),10)),stringsAsFactors=TRUE)
row.names(ColSideAnn)<-colnames(rnormData)
RowSideColors<-colorRampPalette(c("chartreuse4", "white", 
                                  "firebrick"))(40)
result<-heatmap3(rnormData,ColSideCut=1.2,ColSideAnn=ColSideAnn,
                 ColSideFun=function(x) showAnn(x),ColSideWidth=0.8,
                 RowSideColors=RowSideColors,col=colorRampPalette(c("green","black"
                                                                    , "red"))(1024),RowAxisColors=1,legendfun=function() 
                                                                      showLegend(legend=c("Low","High"),col=c("chartreuse4","firebrick"))
                 ,verbose=TRUE)
#annotations distribution in different clusters and the result 
#of statistic tests
result$cutTable

############################################################
##
## Making a heat map of REAP
##
##
############################################################

### This is the excel file where I custom annotated each protein by its color and functionality
Immune_autoantibodies$color
REAP_Immune_autoantibodies_reordered <- as.data.frame(read_excel("C:/Users/pmunson/OneDrive - Nurix Therapeutics, Inc/Desktop/Autoantibodies/REAP_Immune_autoantibodies_reordered.xlsx"))
rownames(REAP_Immune_autoantibodies_reordered) <-REAP_Immune_autoantibodies_reordered[,1]
REAP_Immune_autoantibodies_reordered <- REAP_Immune_autoantibodies_reordered[,-c(1)]

############################
 install.packages("pheatmap")
library(pheatmap)

Patient_reordered_timepoints <- as.data.frame(read_excel("C:/Users/pmunson/OneDrive - Nurix Therapeutics, Inc/Desktop/Autoantibodies/Patient_reordered_timepoints.xlsx"))
Patient_reordered_timepoints[,2]
t1 <- gsub("DAY_0", "black", Patient_reordered_timepoints[,2])
t2 <- gsub("DAY_43", "blue", t1)
Timepoints <- gsub("DAY_89.101", "red", t2)

Clinical_Outcome <- rep(c("#006400","#4d934d","#99c199","#F1959B","#DC1C13","white"),c(6,24,18,14,24,5))
pheat_map_row_colors <- as.data.frame(cbind(rownames(REAP_Immune_autoantibodies_reordered), Immune_autoantibodies$color))
  
rownames(REAP_Immune_autoantibodies_reordered)


pheatmap(REAP_Immune_autoantibodies_reordered,RowSideColors =  Immune_autoantibodies$color, Colv = NA, Rowv = NA,  
         col=  colorRampPalette(c("black", "#483248", "#685BC7", "purple","#F65058","#FF6C2F","orange", "yellow"))(350),
         ColSideColors = cbind(Clinical_Outcome = Clinical_Outcome, Timepoints = Timepoints))
         #,margins = c(30,10),legendfun=function()
           #plot(0,0,bty="n",xaxt="n",yaxt="n",type="n"), cexRow = 0.2 + 1/log10(nrow(REAP_Immune_autoantibodies_reordered)))

pheatmap(as.matrix(REAP_Immune_autoantibodies_reordered), color =colorRampPalette(c("black", "#483248", "#685BC7", "purple","#F65058","#FF6C2F","orange", "yellow"))(350),
         cluster_row = FALSE, cluster_cols = FALSE,
         annotation_row = pheat_map_row_colors)

pheatmap(REAP_Immune_autoantibodies_reordered,RowSideColors =  Immune_autoantibodies$color,   
         col=  colorRampPalette(c("black", "#483248", "#685BC7", "purple","#F65058","#FF6C2F","orange", "yellow"))(350),
         ColSideColors = cbind(Clinical_Outcome = Clinical_Outcome, Timepoints = Timepoints),
         margins = c(20,10))

set.seed(2)
heatmap3(REAP_Immune_autoantibodies_reordered,RowSideColors =  Immune_autoantibodies$color, Rowv = NA,  
         col=  colorRampPalette(c("black", "#483248", "#685BC7", "purple","#F65058","#FF6C2F","orange", "yellow"))(350),
         ColSideColors = cbind(Clinical_Outcome = Clinical_Outcome, Timepoints = Timepoints),
         margins = c(30,10))

set.seed(2)
heatmap3(REAP_Immune_autoantibodies_reordered,RowSideColors =  Immune_autoantibodies$color, Rowv = NA,  
         col=  colorRampPalette(c("white","#685BC7", "purple","#F65058","#FF6C2F","orange", "yellow"))(350),
         ColSideColors = cbind(Clinical_Outcome = Clinical_Outcome, Timepoints = Timepoints),
         margins = c(30,10))














dim(Paul_score)

CTLA4_Baseline_REAP <- Tidy_REAP_Clinical %>%
                          filter(Checkpoint != "None" & Protein == "CTLA4")
 

unique(CTLA4_Baseline_REAP$Checkpoint)
 
setwd("/Users/pmunson/Box/Paul Munson/Autoantibody Project")
getwd()
#write.csv(CTLA4_Baseline_REAP,"CTLA4_Baseline_REAP.csv")

PD1_Baseline_REAP <- Tidy_REAP_Clinical %>%
  filter(Checkpoint == "anti-CTLA-4 + anti-PD-1" & Protein == "PDCD1")

setwd("/Users/pmunson/Box/Paul Munson/Autoantibody Project")
getwd()
#write.csv(PD1_Baseline_REAP,"PD1_Baseline_REAP.csv")


############

#

#  Number of immune hits and clinical outcomes

#

############

Melanoma_Vaccine_Clinical_Data_Butterfield_Lab <- read_excel("Melanoma Vaccine Clinical Data Butterfield Lab.xlsx")
colnames(Melanoma_Vaccine_Clinical_Data_Butterfield_Lab )

Melanoma_Vaccine_Clinical_Data_Butterfield_Lab %>%
  ggplot(., aes(x = PFS, y = Number_of_immune_hits, color = Best_Clinical_outcome)) +
  geom_jitter(size = 3) + 
  scale_color_manual(values = c("black", "#00AC8C",  "#FF6C2F", "#F65058" , "#685BC7", "#0083C3")) +
  theme_par() +
  xlab("PFS (Months)") +
  ylab("Immune Hits (Value >1))") 

cor()














