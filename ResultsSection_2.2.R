### Results - Section 2.2 ### ----

### Load Packages ### ----
library(dplyr)
library(psych)
library(ggplot2)
library(ggpubr)
library(dgof)
library(LSD)
library(readxl)
library(naniar)

############################################################################
########################## Scenario Calculations ###########################
############################################################################

### Read in Data ### ----
load("scen1_icc_s0_t1t2_ela01_revision.RData")
icc.comb.scen1 <- dplyr::filter(icc.comb.scen1, icc.comb.scen1$covariates == "Adjusted for DNAm_Monocytes, Methylation_Batch")
icc.comb.scen1 <- dplyr::filter(icc.comb.scen1, icc.comb.scen1$type == "ICC(2,1)")[,c(1,5)]
colnames(icc.comb.scen1)[2] <- "ICC_Scen1"
load("scen2_icc_s0_t1t3_ela01_revision.RData")
icc.comb.scen2 <- dplyr::filter(icc.comb.scen2, icc.comb.scen2$covariates == "Adjusted for DNAm_Monocytes, Methylation_Batch")
icc.comb.scen2 <- dplyr::filter(icc.comb.scen2, icc.comb.scen2$type == "ICC(2,1)")[,c(1,5)]
colnames(icc.comb.scen2)[2] <- "ICC_Scen2"
load("scen3_icc_s0_t1t4_ela01_revision.RData")
icc.comb.scen3 <- dplyr::filter(icc.comb.scen3, icc.comb.scen3$covariates == "Adjusted for DNAm_Monocytes, Methylation_Batch")
icc.comb.scen3 <- dplyr::filter(icc.comb.scen3, icc.comb.scen3$type == "ICC(2,1)")[,c(1,5)]
colnames(icc.comb.scen3)[2] <- "ICC_Scen3"
load("scen4_icc_s0_t3t4_ela01_revision.RData")
icc.comb.scen4 <- dplyr::filter(icc.comb.scen4, icc.comb.scen4$covariates == "Adjusted for DNAm_Monocytes, Methylation_Batch")
icc.comb.scen4 <- dplyr::filter(icc.comb.scen4, icc.comb.scen4$type == "ICC(2,1)")[,c(1,5)]
colnames(icc.comb.scen4)[2] <- "ICC_Scen4"
load("scen5_icc_s0_t1t2t3t4_ela01_revision.RData")
icc.comb.scen5 <- dplyr::filter(icc.comb.scen5, icc.comb.scen5$covariates == "Adjusted for DNAm_Monocytes, Methylation_Batch")
icc.comb.scen5 <- dplyr::filter(icc.comb.scen5, icc.comb.scen5$type == "ICC(2,1)")[,c(1,5)]
colnames(icc.comb.scen5)[2] <- "ICC_Scen5"

load("scen6sampled_icc_s1_t1t2_ela01_revision.RData")
icc.comb.scen6sampled <- dplyr::filter(icc.comb.scen6sampled, icc.comb.scen6sampled$covariates == "Adjusted for DNAm_Monocytes, Methylation_Batch")
icc.comb.scen6sampled <- dplyr::filter(icc.comb.scen6sampled, icc.comb.scen6sampled$type == "ICC(2,k)")[,c(1,5)]
colnames(icc.comb.scen6sampled)[2] <- "ICC_Scen6samp"
load("scen7sampled_icc_s1_t1t3_ela01_revision.RData")
icc.comb.scen7sampled <- dplyr::filter(icc.comb.scen7sampled, icc.comb.scen7sampled$covariates == "Adjusted for DNAm_Monocytes, Methylation_Batch")
icc.comb.scen7sampled <- dplyr::filter(icc.comb.scen7sampled, icc.comb.scen7sampled$type == "ICC(2,k)")[,c(1,5)]
colnames(icc.comb.scen7sampled)[2] <- "ICC_Scen7samp"
load("scen8sampled_icc_s1_t1t4_ela01_revision.RData")
icc.comb.scen8sampled <- dplyr::filter(icc.comb.scen8sampled, icc.comb.scen8sampled$covariates == "Adjusted for DNAm_Monocytes, Methylation_Batch")
icc.comb.scen8sampled <- dplyr::filter(icc.comb.scen8sampled, icc.comb.scen8sampled$type == "ICC(2,k)")[,c(1,5)]
colnames(icc.comb.scen8sampled)[2] <- "ICC_Scen8samp"
load("scen9sampled_icc_s1_t3t4_ela01_revision.RData")
icc.comb.scen9sampled <- dplyr::filter(icc.comb.scen9sampled, icc.comb.scen9sampled$covariates == "Adjusted for DNAm_Monocytes, Methylation_Batch")
icc.comb.scen9sampled <- dplyr::filter(icc.comb.scen9sampled, icc.comb.scen9sampled$type == "ICC(2,k)")[,c(1,5)]
colnames(icc.comb.scen9sampled)[2] <- "ICC_Scen9samp"
load("scen10_icc_s0s1_t1_ela01_revision.RData")
icc.comb.scen10 <- dplyr::filter(icc.comb.scen10, icc.comb.scen10$covariates == "Adjusted for DNAm_Monocytes, Methylation_Batch")
icc.comb.scen10 <- dplyr::filter(icc.comb.scen10, icc.comb.scen10$type == "ICC(2,1)")[,c(1,5)]
colnames(icc.comb.scen10)[2] <- "ICC_Scen10"

### Merge all Data ### ----
data <- merge(icc.comb.scen1, icc.comb.scen2, by = "DNAm.probe")
data <- merge(data, icc.comb.scen3, by = "DNAm.probe")
data <- merge(data, icc.comb.scen4, by = "DNAm.probe")
data <- merge(data, icc.comb.scen5, by = "DNAm.probe")
data <- merge(data, icc.comb.scen6sampled, by = "DNAm.probe")
data <- merge(data, icc.comb.scen7sampled, by = "DNAm.probe")
data <- merge(data, icc.comb.scen8sampled, by = "DNAm.probe")
data <- merge(data, icc.comb.scen9sampled, by = "DNAm.probe")
data <- merge(data, icc.comb.scen10, by = "DNAm.probe")
Exclude_Probes_Comb <- readRDS("Exclude_Probes_Comb.Rds")
keep <- !(data$DNAm.probe %in% Exclude_Probes_Comb$DNAm.probe)
data <- data[keep,]

# Delete old dataframes #
icc.comb.scen1 <- NULL
icc.comb.scen2 <- NULL
icc.comb.scen3 <- NULL
icc.comb.scen4 <- NULL
icc.comb.scen5 <- NULL
icc.comb.scen6sampled <- NULL
icc.comb.scen7sampled <- NULL
icc.comb.scen8sampled <- NULL
icc.comb.scen9sampled <- NULL
icc.comb.scen10 <- NULL

### Remove ICC < 0.01 values ### ----
data_noZeros <- data
data_noZeros[,2:11][data_noZeros[,2:11] < 0.01] <- NA
colnames(data_noZeros)[2:11] <- c("NoStressT1-2",
                                  "NoStressT1-3",
                                  "NoStressT1-4",
                                  "NoStressT3-4",
                                  "NoStressT1-2-3-4",
                                  "StressT1-2",
                                  "StressT1-3",
                                  "StressT1-4",
                                  "StressT3-4",
                                  "CrossSessionT1")

### Descriptive Statistics ### ----
# Descriptive statistics with no zero ICC values #
describe(data_noZeros)

### Counting the Number of Probes Stratified by ICC Category (with Zeros) ### ----
length(data$ICC_Scen1[data$ICC_Scen1 <= 0.01])
length(data$ICC_Scen1[data$ICC_Scen1 > 0.01 & data$ICC_Scen1 < 0.50])
length(data$ICC_Scen1[data$ICC_Scen1 > 0.50 & data$ICC_Scen1 < 0.75])
length(data$ICC_Scen1[data$ICC_Scen1 > 0.75 & data$ICC_Scen1 < 0.90])
length(data$ICC_Scen1[data$ICC_Scen1 > 0.90])

length(data$ICC_Scen2[data$ICC_Scen2 <= 0.01])
length(data$ICC_Scen2[data$ICC_Scen2 > 0.01 & data$ICC_Scen2 < 0.50])
length(data$ICC_Scen2[data$ICC_Scen2 > 0.50 & data$ICC_Scen2 < 0.75])
length(data$ICC_Scen2[data$ICC_Scen2 > 0.75 & data$ICC_Scen2 < 0.90])
length(data$ICC_Scen2[data$ICC_Scen2 > 0.90])

length(data$ICC_Scen3[data$ICC_Scen3 <= 0.01])
length(data$ICC_Scen3[data$ICC_Scen3 > 0.01 & data$ICC_Scen3 < 0.50])
length(data$ICC_Scen3[data$ICC_Scen3 > 0.50 & data$ICC_Scen3 < 0.75])
length(data$ICC_Scen3[data$ICC_Scen3 > 0.75 & data$ICC_Scen3 < 0.90])
length(data$ICC_Scen3[data$ICC_Scen3 > 0.90])

length(data$ICC_Scen4[data$ICC_Scen4 <= 0.01])
length(data$ICC_Scen4[data$ICC_Scen4 > 0.01 & data$ICC_Scen4 < 0.50])
length(data$ICC_Scen4[data$ICC_Scen4 > 0.50 & data$ICC_Scen4 < 0.75])
length(data$ICC_Scen4[data$ICC_Scen4 > 0.75 & data$ICC_Scen4 < 0.90])
length(data$ICC_Scen4[data$ICC_Scen4 > 0.90])

length(data$ICC_Scen5[data$ICC_Scen5 <= 0.01])
length(data$ICC_Scen5[data$ICC_Scen5 > 0.01 & data$ICC_Scen5 < 0.50])
length(data$ICC_Scen5[data$ICC_Scen5 > 0.50 & data$ICC_Scen5 < 0.75])
length(data$ICC_Scen5[data$ICC_Scen5 > 0.75 & data$ICC_Scen5 < 0.90])
length(data$ICC_Scen5[data$ICC_Scen5 > 0.90])

length(data$ICC_Scen6[data$ICC_Scen6 <= 0.01])
length(data$ICC_Scen6[data$ICC_Scen6 > 0.01 & data$ICC_Scen6 < 0.50])
length(data$ICC_Scen6[data$ICC_Scen6 > 0.50 & data$ICC_Scen6 < 0.75])
length(data$ICC_Scen6[data$ICC_Scen6 > 0.75 & data$ICC_Scen6 < 0.90])
length(data$ICC_Scen6[data$ICC_Scen6 > 0.90])

length(data$ICC_Scen7[data$ICC_Scen7 <= 0.01])
length(data$ICC_Scen7[data$ICC_Scen7 > 0.01 & data$ICC_Scen7 < 0.50])
length(data$ICC_Scen7[data$ICC_Scen7 > 0.50 & data$ICC_Scen7 < 0.75])
length(data$ICC_Scen7[data$ICC_Scen7 > 0.75 & data$ICC_Scen7 < 0.90])
length(data$ICC_Scen7[data$ICC_Scen7 > 0.90])

length(data$ICC_Scen8[data$ICC_Scen8 <= 0.01])
length(data$ICC_Scen8[data$ICC_Scen8 > 0.01 & data$ICC_Scen8 < 0.50])
length(data$ICC_Scen8[data$ICC_Scen8 > 0.50 & data$ICC_Scen8 < 0.75])
length(data$ICC_Scen8[data$ICC_Scen8 > 0.75 & data$ICC_Scen8 < 0.90])
length(data$ICC_Scen8[data$ICC_Scen8 > 0.90])

length(data$ICC_Scen9[data$ICC_Scen9 <= 0.01])
length(data$ICC_Scen9[data$ICC_Scen9 > 0.01 & data$ICC_Scen9 < 0.50])
length(data$ICC_Scen9[data$ICC_Scen9 > 0.50 & data$ICC_Scen9 < 0.75])
length(data$ICC_Scen9[data$ICC_Scen9 > 0.75 & data$ICC_Scen9 < 0.90])
length(data$ICC_Scen9[data$ICC_Scen9 > 0.90])

length(data$ICC_Scen10[data$ICC_Scen10 <= 0.01])
length(data$ICC_Scen10[data$ICC_Scen10 > 0.01 & data$ICC_Scen10 < 0.50])
length(data$ICC_Scen10[data$ICC_Scen10 > 0.50 & data$ICC_Scen10 < 0.75])
length(data$ICC_Scen10[data$ICC_Scen10 > 0.75 & data$ICC_Scen10 < 0.90])
length(data$ICC_Scen10[data$ICC_Scen10 > 0.90])

### Probe Correlations Across Scenarios ### ----
data_cor <- cor(data_noZeros[,2:11], use = "pairwise.complete.obs")

### Performing Hierarchical Clustering ### ----
clusters <- hclust(dist(t(data_noZeros[,2:11])))
plot(clusters, xlab = "")

### Plot All Scenarios Histograms ### ----
data_noZeros <- data
data_noZeros[,2:11][data_noZeros[,2:11] < 0.01] <- NA

hist_31 <- ggplot(data_noZeros, aes(x=ICC_Scen1)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.47),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,60000) +
  labs(x = NULL, y = NULL, title = ("NoStressT1-2")) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())  #remove y axis ticks

hist_32 <- ggplot(data_noZeros, aes(x=ICC_Scen2)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.45),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,60000) +
  labs(x = NULL, y = NULL, title = ("NoStressT1-3")) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        plot.title = element_text(hjust = 0.5))  

hist_33 <- ggplot(data_noZeros, aes(x=ICC_Scen3)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.31),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,60000) +
  labs(x = NULL, y = NULL, title = ("NoStressT1-4")) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        plot.title = element_text(hjust = 0.5))  

hist_34 <- ggplot(data_noZeros, aes(x=ICC_Scen4)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.50),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,60000) +
  labs(x = NULL, y = NULL, title = ("NoStressT3-4")) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        plot.title = element_text(hjust = 0.5))  

hist_35 <- ggplot(data_noZeros, aes(x=ICC_Scen5)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.44),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,60000) +
  labs(x = NULL, y = NULL, title = ("NoStressT1-2-3-4")) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        plot.title = element_text(hjust = 0.5))  

hist_36 <- ggplot(data_noZeros, aes(x=ICC_Scen6samp)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.64),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,60000) +
  labs(x = NULL, y = NULL, title = ("StressT1-2")) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        plot.title = element_text(hjust = 0.5))  

hist_37 <- ggplot(data_noZeros, aes(x=ICC_Scen7samp)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.64),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,60000) +
  labs(x = NULL, y = NULL, title = ("StressT1-3")) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        plot.title = element_text(hjust = 0.5))  

hist_38 <- ggplot(data_noZeros, aes(x=ICC_Scen8samp)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.65),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,60000) +
  labs(x = NULL, y = NULL, title = ("StressT1-4")) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        plot.title = element_text(hjust = 0.5))  

hist_39 <- ggplot(data_noZeros, aes(x=ICC_Scen9samp)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.62),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,60000) +
  labs(x = NULL, y = NULL, title = ("StressT3-4")) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        plot.title = element_text(hjust = 0.5))  

hist_310 <- ggplot(data_noZeros, aes(x=ICC_Scen10)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.28),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,60000) +
  labs(x = NULL, y = NULL, title = ("CrossSessionT1")) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        plot.title = element_text(hjust = 0.5)) 

ggarrange(hist_31, hist_32, hist_33, hist_34, hist_35,
          hist_36, hist_37, hist_38, hist_39, hist_310,
          nrow = 2, ncol = 5)


