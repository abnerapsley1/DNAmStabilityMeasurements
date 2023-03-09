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
load("scen1_icc_s0_t1t2.RData")
scen1.icc.comb <- dplyr::filter(scen1.icc.comb, scen1.icc.comb$covariates == "Adjusted for DNAm_Monocytes")
scen1.icc.comb <- dplyr::filter(scen1.icc.comb, scen1.icc.comb$type == "ICC(2,1)")[,c(1,5)]
colnames(scen1.icc.comb)[2] <- "ICC_Scen1"
load("scen2_icc_s0_t1t3.RData")
scen2.icc.comb <- dplyr::filter(scen3.icc.comb, scen3.icc.comb$covariates == "Adjusted for DNAm_Monocytes")
scen2.icc.comb <- dplyr::filter(scen2.icc.comb, scen2.icc.comb$type == "ICC(2,1)")[,c(1,5)]
colnames(scen2.icc.comb)[2] <- "ICC_Scen2"
load("scen3_icc_s0_t1t4.RData")
scen3.icc.comb <- dplyr::filter(scen11.icc.comb, scen11.icc.comb$covariates == "Adjusted for DNAm_Monocytes")
scen3.icc.comb <- dplyr::filter(scen3.icc.comb, scen3.icc.comb$type == "ICC(2,1)")[,c(1,5)]
colnames(scen3.icc.comb)[2] <- "ICC_Scen3"
load("scen4_icc_s0_t3t4.RData")
scen4.icc.comb <- dplyr::filter(scen4.icc.comb, scen4.icc.comb$covariates == "Adjusted for DNAm_Monocytes")
scen4.icc.comb <- dplyr::filter(scen4.icc.comb, scen4.icc.comb$type == "ICC(2,1)")[,c(1,5)]
colnames(scen4.icc.comb)[2] <- "ICC_Scen4"
load("scen5_icc_s0_t1t2t3t4.RData")
scen5.icc.comb <- dplyr::filter(scen5.icc.comb, scen5.icc.comb$covariates == "Adjusted for DNAm_Monocytes")
scen5.icc.comb <- dplyr::filter(scen5.icc.comb, scen5.icc.comb$type == "ICC(2,1)")[,c(1,5)]
colnames(scen5.icc.comb)[2] <- "ICC_Scen5"
load("scen6sampled_icc_s1_t1t2.RData")
scen6sampled.icc.comb <- dplyr::filter(scen6sampled.icc.comb, scen6sampled.icc.comb$covariates == "Adjusted for DNAm_Monocytes")
scen6sampled.icc.comb <- dplyr::filter(scen6sampled.icc.comb, scen6sampled.icc.comb$type == "ICC(2,1)")[,c(1,5)]
colnames(scen6sampled.icc.comb)[2] <- "ICC_Scen6"
load("scen7sampled_icc_s1_t1t3.RData")
scen7sampled.icc.comb <- dplyr::filter(scen7sampled.icc.comb, scen7sampled.icc.comb$covariates == "Adjusted for DNAm_Monocytes")
scen7sampled.icc.comb <- dplyr::filter(scen7sampled.icc.comb, scen7sampled.icc.comb$type == "ICC(2,1)")[,c(1,5)]
colnames(scen7sampled.icc.comb)[2] <- "ICC_Scen7"
load("scen8sampled_icc_s1_t1t4.RData")
scen8sampled.icc.comb <- dplyr::filter(scen8sampled.icc.comb, scen8sampled.icc.comb$covariates == "Adjusted for DNAm_Monocytes")
scen8sampled.icc.comb <- dplyr::filter(scen8sampled.icc.comb, scen8sampled.icc.comb$type == "ICC(2,1)")[,c(1,5)]
colnames(scen8sampled.icc.comb)[2] <- "ICC_Scen8"
load("scen9sampled_icc_s1_t3t4.RData")
scen9sampled.icc.comb <- dplyr::filter(scen9sampled.icc.comb, scen9sampled.icc.comb$covariates == "Adjusted for DNAm_Monocytes")
scen9sampled.icc.comb <- dplyr::filter(scen9sampled.icc.comb, scen9sampled.icc.comb$type == "ICC(2,1)")[,c(1,5)]
colnames(scen9sampled.icc.comb)[2] <- "ICC_Scen9"
load("scen10_icc_s1s2_t1.RData")
scen10.icc.comb <- dplyr::filter(scen10.icc.comb, scen10.icc.comb$covariates == "Adjusted for DNAm_Monocytes")
scen10.icc.comb <- dplyr::filter(scen10.icc.comb, scen10.icc.comb$type == "ICC(2,1)")[,c(1,5)]
colnames(scen10.icc.comb)[2] <- "ICC_Scen10"

### Merge all Data ### ----
data <- merge(scen1.icc.comb, scen2.icc.comb, by = "DNAm.probe")
data <- merge(data, scen3.icc.comb, by = "DNAm.probe")
data <- merge(data, scen4.icc.comb, by = "DNAm.probe")
data <- merge(data, scen5.icc.comb, by = "DNAm.probe")
data <- merge(data, scen6sampled.icc.comb, by = "DNAm.probe")
data <- merge(data, scen7sampled.icc.comb, by = "DNAm.probe")
data <- merge(data, scen8sampled.icc.comb, by = "DNAm.probe")
data <- merge(data, scen9sampled.icc.comb, by = "DNAm.probe")
data <- merge(data, scen10.icc.comb, by = "DNAm.probe")
Exclude_Probes_Comb <- readRDS("Exclude_Probes_Comb.Rds")
keep <- !(data$DNAm.probe %in% Exclude_Probes_Comb$DNAm.probe)
data <- data[keep,]

# Delete old dataframes #
scen1.icc.comb <- NULL
scen2.icc.comb <- NULL
scen3.icc.comb <- NULL
scen4.icc.comb <- NULL
scen5.icc.comb <- NULL
scen6sampled.icc.comb <- NULL
scen7sampled.icc.comb <- NULL
scen8sampled.icc.comb <- NULL
scen9sampled.icc.comb <- NULL
scen10.icc.comb <- NULL
scen11.icc.comb <- NULL

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
  geom_vline(aes(xintercept=0.44),
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
  geom_vline(aes(xintercept=0.30),
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
  geom_vline(aes(xintercept=0.42),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,60000) +
  labs(x = NULL, y = NULL, title = ("NoStressT1-2-3-4")) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        plot.title = element_text(hjust = 0.5))  

hist_36 <- ggplot(data_noZeros, aes(x=ICC_Scen6)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.50),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,60000) +
  labs(x = NULL, y = NULL, title = ("StressT1-2")) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        plot.title = element_text(hjust = 0.5))  

hist_37 <- ggplot(data_noZeros, aes(x=ICC_Scen7)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.50),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,60000) +
  labs(x = NULL, y = NULL, title = ("StressT1-3")) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        plot.title = element_text(hjust = 0.5))  

hist_38 <- ggplot(data_noZeros, aes(x=ICC_Scen8)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.52),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,60000) +
  labs(x = NULL, y = NULL, title = ("StressT1-4")) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        plot.title = element_text(hjust = 0.5))  

hist_39 <- ggplot(data_noZeros, aes(x=ICC_Scen9)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.48),
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
  geom_vline(aes(xintercept=0.29),
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


