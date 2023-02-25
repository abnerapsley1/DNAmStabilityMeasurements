### Testing Aims 5-7 ### ----

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
################################### Aim 5 ##################################
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

### Remove ICC = 0 values ### ----
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

### Statistical Testing ### ----
t.test(data_noZeros$ICC_Scen2, data_noZeros$ICC_Scen1, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen3, data_noZeros$ICC_Scen1, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen4, data_noZeros$ICC_Scen1, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen5, data_noZeros$ICC_Scen1, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen6, data_noZeros$ICC_Scen1, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen7, data_noZeros$ICC_Scen1, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen8, data_noZeros$ICC_Scen1, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen9, data_noZeros$ICC_Scen1, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen10, data_noZeros$ICC_Scen1, paired = "TRUE", na.action = na.omit)

t.test(data_noZeros$ICC_Scen3, data_noZeros$ICC_Scen2, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen4, data_noZeros$ICC_Scen2, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen5, data_noZeros$ICC_Scen2, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen6, data_noZeros$ICC_Scen2, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen7, data_noZeros$ICC_Scen2, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen8, data_noZeros$ICC_Scen2, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen9, data_noZeros$ICC_Scen2, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen10, data_noZeros$ICC_Scen2, paired = "TRUE", na.action = na.omit)

t.test(data_noZeros$ICC_Scen4, data_noZeros$ICC_Scen3, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen5, data_noZeros$ICC_Scen3, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen6, data_noZeros$ICC_Scen3, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen7, data_noZeros$ICC_Scen3, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen8, data_noZeros$ICC_Scen3, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen9, data_noZeros$ICC_Scen3, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen10, data_noZeros$ICC_Scen3, paired = "TRUE", na.action = na.omit)

t.test(data_noZeros$ICC_Scen5, data_noZeros$ICC_Scen4, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen6, data_noZeros$ICC_Scen4, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen7, data_noZeros$ICC_Scen4, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen8, data_noZeros$ICC_Scen4, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen9, data_noZeros$ICC_Scen4, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen10, data_noZeros$ICC_Scen4, paired = "TRUE", na.action = na.omit)

t.test(data_noZeros$ICC_Scen6, data_noZeros$ICC_Scen5, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen7, data_noZeros$ICC_Scen5, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen8, data_noZeros$ICC_Scen5, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen9, data_noZeros$ICC_Scen5, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen10, data_noZeros$ICC_Scen5, paired = "TRUE", na.action = na.omit)

t.test(data_noZeros$ICC_Scen7, data_noZeros$ICC_Scen6, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen8, data_noZeros$ICC_Scen6, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen9, data_noZeros$ICC_Scen6, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen10, data_noZeros$ICC_Scen6, paired = "TRUE", na.action = na.omit)

t.test(data_noZeros$ICC_Scen8, data_noZeros$ICC_Scen7, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen9, data_noZeros$ICC_Scen7, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen10, data_noZeros$ICC_Scen7, paired = "TRUE", na.action = na.omit)

t.test(data_noZeros$ICC_Scen9, data_noZeros$ICC_Scen8, paired = "TRUE", na.action = na.omit)
t.test(data_noZeros$ICC_Scen10, data_noZeros$ICC_Scen8, paired = "TRUE", na.action = na.omit)

t.test(data_noZeros$ICC_Scen10, data_noZeros$ICC_Scen9, paired = "TRUE", na.action = na.omit)

### Counting the Number of Zero Probes ### ----
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
# Descriptive statistics with no zero ICC values #
describe(data_noZeros)

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

############################################################################
################################### Aim 6 ##################################
############################################################################
### Read in Data ### ----
load("scen5a_icc_s0_t1t2t3t4_ela0.RData")
scen5a.ela0.icc.comb <- dplyr::filter(scen5a.ela0.icc.comb, scen5a.ela0.icc.comb$covariates == "Adjusted for DNAm_Monocytes")
scen5a.ela0.icc.comb <- dplyr::filter(scen5a.ela0.icc.comb, scen5a.ela0.icc.comb$type == "ICC(2,1)")[,c(1,5)]
colnames(scen5a.ela0.icc.comb)[2] <- "ICC_Scen5_ELA0"
load("scen5b_icc_s0_t1t2t3t4_ela1.RData")
scen5b.ela1.icc.comb <- dplyr::filter(scen5b.ela1.icc.comb, scen5b.ela1.icc.comb$covariates == "Adjusted for DNAm_Monocytes")
scen5b.ela1.icc.comb <- dplyr::filter(scen5b.ela1.icc.comb, scen5b.ela1.icc.comb$type == "ICC(2,1)")[,c(1,5)]
colnames(scen5b.ela1.icc.comb)[2] <- "ICC_Scen5_ELA1"

load("scen6a_icc_s1_t1t2_ela0.RData")
scen6a.ela0.icc.comb <- dplyr::filter(scen6a.ela0.icc.comb, scen6a.ela0.icc.comb$covariates == "Adjusted for DNAm_Monocytes")
scen6a.ela0.icc.comb <- dplyr::filter(scen6a.ela0.icc.comb, scen6a.ela0.icc.comb$type == "ICC(2,1)")[,c(1,5)]
colnames(scen6a.ela0.icc.comb)[2] <- "ICC_Scen6_ELA0"
load("scen6b_icc_s1_t1t2_ela1.RData")
scen6b.ela1.icc.comb <- dplyr::filter(scen6b.ela1.icc.comb, scen6b.ela1.icc.comb$covariates == "Adjusted for DNAm_Monocytes")
scen6b.ela1.icc.comb <- dplyr::filter(scen6b.ela1.icc.comb, scen6b.ela1.icc.comb$type == "ICC(2,1)")[,c(1,5)]
colnames(scen6b.ela1.icc.comb)[2] <- "ICC_Scen6_ELA1"

load("scen7a_icc_s1_t1t3_ela0.RData")
scen7a.ela0.icc.comb <- dplyr::filter(scen7a.ela0.icc.comb, scen7a.ela0.icc.comb$covariates == "Adjusted for DNAm_Monocytes")
scen7a.ela0.icc.comb <- dplyr::filter(scen7a.ela0.icc.comb, scen7a.ela0.icc.comb$type == "ICC(2,1)")[,c(1,5)]
colnames(scen7a.ela0.icc.comb)[2] <- "ICC_Scen7_ELA0"
load("scen7b_icc_s1_t1t3_ela1.RData")
scen7b.ela1.icc.comb <- dplyr::filter(scen7b.ela1.icc.comb, scen7b.ela1.icc.comb$covariates == "Adjusted for DNAm_Monocytes")
scen7b.ela1.icc.comb <- dplyr::filter(scen7b.ela1.icc.comb, scen7b.ela1.icc.comb$type == "ICC(2,1)")[,c(1,5)]
colnames(scen7b.ela1.icc.comb)[2] <- "ICC_Scen7_ELA1"

load("scen8a_icc_s1_t1t4_ela0.RData")
scen8a.ela0.icc.comb <- dplyr::filter(scen8a.ela0.icc.comb, scen8a.ela0.icc.comb$covariates == "Adjusted for DNAm_Monocytes")
scen8a.ela0.icc.comb <- dplyr::filter(scen8a.ela0.icc.comb, scen8a.ela0.icc.comb$type == "ICC(2,1)")[,c(1,5)]
colnames(scen8a.ela0.icc.comb)[2] <- "ICC_Scen8_ELA0"
load("scen8b_icc_s1_t1t4_ela1.RData")
scen8b.ela1.icc.comb <- dplyr::filter(scen8b.ela1.icc.comb, scen8b.ela1.icc.comb$covariates == "Adjusted for DNAm_Monocytes")
scen8b.ela1.icc.comb <- dplyr::filter(scen8b.ela1.icc.comb, scen8b.ela1.icc.comb$type == "ICC(2,1)")[,c(1,5)]
colnames(scen8b.ela1.icc.comb)[2] <- "ICC_Scen8_ELA1"

load("scen9a_icc_s1_t3t4_ela0.RData")
scen9a.ela0.icc.comb <- dplyr::filter(scen9a.ela0.icc.comb, scen9a.ela0.icc.comb$covariates == "Adjusted for DNAm_Monocytes")
scen9a.ela0.icc.comb <- dplyr::filter(scen9a.ela0.icc.comb, scen9a.ela0.icc.comb$type == "ICC(2,1)")[,c(1,5)]
colnames(scen9a.ela0.icc.comb)[2] <- "ICC_Scen9_ELA0"
load("scen9b_icc_s1_t3t4_ela1.RData")
scen9b.ela1.icc.comb <- dplyr::filter(scen9b.ela1.icc.comb, scen9b.ela1.icc.comb$covariates == "Adjusted for DNAm_Monocytes")
scen9b.ela1.icc.comb <- dplyr::filter(scen9b.ela1.icc.comb, scen9b.ela1.icc.comb$type == "ICC(2,1)")[,c(1,5)]
colnames(scen9b.ela1.icc.comb)[2] <- "ICC_Scen9_ELA1"

### Merge all Data ### ----
data <- merge(scen5a.ela0.icc.comb, scen5b.ela1.icc.comb, by = "DNAm.probe")
data <- merge(data, scen6a.ela0.icc.comb, by = "DNAm.probe")
data <- merge(data, scen6b.ela1.icc.comb, by = "DNAm.probe")
data <- merge(data, scen7a.ela0.icc.comb, by = "DNAm.probe")
data <- merge(data, scen7b.ela1.icc.comb, by = "DNAm.probe")
data <- merge(data, scen8a.ela0.icc.comb, by = "DNAm.probe")
data <- merge(data, scen8b.ela1.icc.comb, by = "DNAm.probe")
data <- merge(data, scen9a.ela0.icc.comb, by = "DNAm.probe")
data <- merge(data, scen9b.ela1.icc.comb, by = "DNAm.probe")
Exclude_Probes_Comb <- readRDS("Exclude_Probes_Comb.Rds")
keep <- !(data$DNAm.probe %in% Exclude_Probes_Comb$DNAm.probe)
data <- data[keep,]

# Delete old dataframes #
scen5a.ela0.icc.comb <- NULL
scen5b.ela1.icc.comb <- NULL
scen6a.ela0.icc.comb <- NULL
scen6b.ela1.icc.comb <- NULL
scen7a.ela0.icc.comb <- NULL
scen7b.ela1.icc.comb <- NULL
scen8a.ela0.icc.comb <- NULL
scen8b.ela1.icc.comb <- NULL
scen9a.ela0.icc.comb <- NULL
scen9b.ela1.icc.comb <- NULL

### Generate Summary Statistics ### ----
describe(data[,2:11])

### Generate Combined Histogram Plot ### ----
hist_5a <- ggplot(data, aes(x=ICC_Scen5_ELA0)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.41),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,100000) +
  labs(title = "Scenario 5", x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())  #remove y axis ticks

hist_5b <- ggplot(data, aes(x=ICC_Scen5_ELA1)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.53),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,100000) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        plot.title = element_text(hjust = 0.5))  

hist_6a <- ggplot(data, aes(x=ICC_Scen6_ELA0)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.53),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,100000) +
  labs(title = "Scenario 6", x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())  #remove y axis ticks

hist_6b <- ggplot(data, aes(x=ICC_Scen6_ELA1)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.39),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,100000) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        plot.title = element_text(hjust = 0.5))  

hist_7a <- ggplot(data, aes(x=ICC_Scen7_ELA0)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.49),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,100000) +
  labs(title = "Scenario 7", x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        plot.title = element_text(hjust = 0.5))  

hist_7b <- ggplot(data, aes(x=ICC_Scen7_ELA1)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.45),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,100000) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        plot.title = element_text(hjust = 0.5))  

hist_8a <- ggplot(data, aes(x=ICC_Scen8_ELA0)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.49),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,100000) +
  labs(title = "Scenario 8", x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        plot.title = element_text(hjust = 0.5))  

hist_8b <- ggplot(data, aes(x=ICC_Scen8_ELA1)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.44),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,100000) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        plot.title = element_text(hjust = 0.5))  

hist_9a <- ggplot(data, aes(x=ICC_Scen9_ELA0)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.45),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,100000) +
  labs(title = "Scenario 9", x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        plot.title = element_text(hjust = 0.5))  

hist_9b <- ggplot(data, aes(x=ICC_Scen9_ELA1)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.43),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,100000) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        plot.title = element_text(hjust = 0.5))  

### Combine All Histogram Plots ### ----
ggarrange(hist_5a, hist_6a, hist_7a, hist_8a, hist_9a,
          hist_5b, hist_6b, hist_7b, hist_8b, hist_9b,
          ncol = 5, nrow = 2)

### Perform Statistical Testing ### ----
## Scenario 5 ##
# Paired T-test #
t.test(data$ICC_Scen5_ELA0, data$ICC_Scen5_ELA1, paired = TRUE)
# Proportional bias test #
data$diff <- data$ICC_Scen5_ELA0 - data$ICC_Scen5_ELA1
data$mean <- (data$ICC_Scen5_ELA0 + data$ICC_Scen5_ELA1)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$ICC_Scen5_ELA0, data$ICC_Scen5_ELA1)
# Chi-square test #
data_chisq <- data.frame(ICC_Scen5_ELA0 = c(length(data$ICC_Scen5_ELA0[data$ICC_Scen5_ELA0 <= 0.01]), 
                                       length(data$ICC_Scen5_ELA0) - length(data$ICC_Scen5_ELA0[data$ICC_Scen5_ELA0 <= 0.01])),
                         ICC_Scen5_ELA1 = c(length(data$ICC_Scen5_ELA1[data$ICC_Scen5_ELA1 <= 0.01]), 
                                       length(data$ICC_Scen5_ELA1) - length(data$ICC_Scen5_ELA1[data$ICC_Scen5_ELA1 <= 0.01])))
chisq.test(data_chisq)

heatscatter(data$ICC_Scen5_ELA0,data$ICC_Scen5_ELA1,
            colpal="bl2gr2rd",
            main="NoStressT1-2-3-4",
            xlab = "Control Group",
            ylab = "ELA Group",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

## Scenario 6 ##
# Paired T-test #
t.test(data$ICC_Scen6_ELA0, data$ICC_Scen6_ELA1, paired = TRUE)
# Proportional bias test #
data$diff <- data$ICC_Scen6_ELA0 - data$ICC_Scen6_ELA1
data$mean <- (data$ICC_Scen6_ELA0 + data$ICC_Scen6_ELA1)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$ICC_Scen6_ELA0, data$ICC_Scen6_ELA1)
# Chi-square test #
data_chisq <- data.frame(ICC_Scen6_ELA0 = c(length(data$ICC_Scen6_ELA0[data$ICC_Scen6_ELA0 <= 0.01]), 
                                            length(data$ICC_Scen6_ELA0) - length(data$ICC_Scen6_ELA0[data$ICC_Scen6_ELA0 <= 0.01])),
                         ICC_Scen6_ELA1 = c(length(data$ICC_Scen6_ELA1[data$ICC_Scen6_ELA1 <= 0.01]), 
                                            length(data$ICC_Scen6_ELA1) - length(data$ICC_Scen6_ELA1[data$ICC_Scen6_ELA1 <= 0.01])))
chisq.test(data_chisq)

heatscatter(data$ICC_Scen6_ELA0,data$ICC_Scen6_ELA1,
            colpal="bl2gr2rd",
            main="StressT1-2",
            xlab = "Control Group",
            ylab = "ELA Group",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

## Scenario 7 ##
# Paired T-test #
t.test(data$ICC_Scen7_ELA0, data$ICC_Scen7_ELA1, paired = TRUE)
# Proportional bias test #
data$diff <- data$ICC_Scen7_ELA0 - data$ICC_Scen7_ELA1
data$mean <- (data$ICC_Scen7_ELA0 + data$ICC_Scen7_ELA1)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$ICC_Scen7_ELA0, data$ICC_Scen7_ELA1)
# Chi-square test #
data_chisq <- data.frame(ICC_Scen7_ELA0 = c(length(data$ICC_Scen7_ELA0[data$ICC_Scen7_ELA0 <= 0.01]), 
                                            length(data$ICC_Scen7_ELA0) - length(data$ICC_Scen7_ELA0[data$ICC_Scen7_ELA0 <= 0.01])),
                         ICC_Scen7_ELA1 = c(length(data$ICC_Scen7_ELA1[data$ICC_Scen7_ELA1 <= 0.01]), 
                                            length(data$ICC_Scen7_ELA1) - length(data$ICC_Scen7_ELA1[data$ICC_Scen7_ELA1 <= 0.01])))
chisq.test(data_chisq)

heatscatter(data$ICC_Scen7_ELA0,data$ICC_Scen7_ELA1,
            colpal="bl2gr2rd",
            main="StressT1-3",
            xlab = "Control Group",
            ylab = "ELA Group",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

## Scenario 8 ##
# Paired T-test #
t.test(data$ICC_Scen8_ELA0, data$ICC_Scen8_ELA1, paired = TRUE)
# Proportional bias test #
data$diff <- data$ICC_Scen8_ELA0 - data$ICC_Scen8_ELA1
data$mean <- (data$ICC_Scen8_ELA0 + data$ICC_Scen8_ELA1)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$ICC_Scen8_ELA0, data$ICC_Scen8_ELA1)
# Chi-square test #
data_chisq <- data.frame(ICC_Scen8_ELA0 = c(length(data$ICC_Scen8_ELA0[data$ICC_Scen8_ELA0 <= 0.01]), 
                                            length(data$ICC_Scen8_ELA0) - length(data$ICC_Scen8_ELA0[data$ICC_Scen8_ELA0 <= 0.01])),
                         ICC_Scen8_ELA1 = c(length(data$ICC_Scen8_ELA1[data$ICC_Scen8_ELA1 <= 0.01]), 
                                            length(data$ICC_Scen8_ELA1) - length(data$ICC_Scen8_ELA1[data$ICC_Scen8_ELA1 <= 0.01])))
chisq.test(data_chisq)

heatscatter(data$ICC_Scen8_ELA0,data$ICC_Scen8_ELA1,
            colpal="bl2gr2rd",
            main="StressT1-4",
            xlab = "Control Group",
            ylab = "ELA Group",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

## Scenario 9 ##
# Paired T-test #
t.test(data$ICC_Scen9_ELA0, data$ICC_Scen9_ELA1, paired = TRUE)
# Proportional bias test #
data$diff <- data$ICC_Scen9_ELA0 - data$ICC_Scen9_ELA1
data$mean <- (data$ICC_Scen9_ELA0 + data$ICC_Scen9_ELA1)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$ICC_Scen9_ELA0, data$ICC_Scen9_ELA1)
# Chi-square test #
data_chisq <- data.frame(ICC_Scen9_ELA0 = c(length(data$ICC_Scen9_ELA0[data$ICC_Scen9_ELA0 <= 0.01]), 
                                            length(data$ICC_Scen9_ELA0) - length(data$ICC_Scen9_ELA0[data$ICC_Scen9_ELA0 <= 0.01])),
                         ICC_Scen9_ELA1 = c(length(data$ICC_Scen9_ELA1[data$ICC_Scen9_ELA1 <= 0.01]), 
                                            length(data$ICC_Scen9_ELA1) - length(data$ICC_Scen9_ELA1[data$ICC_Scen9_ELA1 <= 0.01])))
chisq.test(data_chisq)

heatscatter(data$ICC_Scen9_ELA0,data$ICC_Scen9_ELA1,
            colpal="bl2gr2rd",
            main="StressT3-4",
            xlab = "Control Group",
            ylab = "ELA Group",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)


############################################################################
################################### Aim 7 ##################################
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

### Compute Correlations for Each Scenario ### ----
data_cor <- cor(data[,2:11])

### HeatScatter Plots of Our Probe Reliabilities vs Sugden Probe Reliabilities ### ----
# Read in Sugden Probes #
Sugden_Probes <- read_xlsx("mmc2.xlsx")
colnames(Sugden_Probes) <- c("DNAm.probe","ICC_Sugden")
#Sugden_Probes <- dplyr::filter(Sugden_Probes, Sugden_Probes$ICC_Sugden >= 0.8)
# Create combined dataframe #
data_comb <- merge(data, Sugden_Probes, by = "DNAm.probe")
# Create difference and mean score variables #
data_comb$Mean_S1_Sug <- (data_comb$ICC_Scen1 + data_comb$ICC_Sugden) / 2 
data_comb$Mean_S2_Sug <- (data_comb$ICC_Scen2 + data_comb$ICC_Sugden) / 2 
data_comb$Mean_S3_Sug <- (data_comb$ICC_Scen3 + data_comb$ICC_Sugden) / 2 
data_comb$Mean_S4_Sug <- (data_comb$ICC_Scen4 + data_comb$ICC_Sugden) / 2 
data_comb$Mean_S5_Sug <- (data_comb$ICC_Scen5 + data_comb$ICC_Sugden) / 2 
data_comb$Mean_S6_Sug <- (data_comb$ICC_Scen6 + data_comb$ICC_Sugden) / 2 
data_comb$Mean_S7_Sug <- (data_comb$ICC_Scen7 + data_comb$ICC_Sugden) / 2 
data_comb$Mean_S8_Sug <- (data_comb$ICC_Scen8 + data_comb$ICC_Sugden) / 2 
data_comb$Mean_S9_Sug <- (data_comb$ICC_Scen9 + data_comb$ICC_Sugden) / 2 
data_comb$Mean_S10_Sug <- (data_comb$ICC_Scen10 + data_comb$ICC_Sugden) / 2 

data_comb$Diff_S1_Sug <- (data_comb$ICC_Scen1 - data_comb$ICC_Sugden)
data_comb$Diff_S2_Sug <- (data_comb$ICC_Scen2 - data_comb$ICC_Sugden)
data_comb$Diff_S3_Sug <- (data_comb$ICC_Scen3 - data_comb$ICC_Sugden)
data_comb$Diff_S4_Sug <- (data_comb$ICC_Scen4 - data_comb$ICC_Sugden)
data_comb$Diff_S5_Sug <- (data_comb$ICC_Scen5 - data_comb$ICC_Sugden)
data_comb$Diff_S6_Sug <- (data_comb$ICC_Scen6 - data_comb$ICC_Sugden)
data_comb$Diff_S7_Sug <- (data_comb$ICC_Scen7 - data_comb$ICC_Sugden)
data_comb$Diff_S8_Sug <- (data_comb$ICC_Scen8 - data_comb$ICC_Sugden)
data_comb$Diff_S9_Sug <- (data_comb$ICC_Scen9 - data_comb$ICC_Sugden)
data_comb$Diff_S10_Sug <- (data_comb$ICC_Scen10 - data_comb$ICC_Sugden)

# Create HeatScatter plots #
# Scenario 1 # 
heatscatter(data_comb$ICC_Scen1,data_comb$ICC_Sugden,
            colpal="bl2gr2rd",
            main="Scenario 1",
            xlab = "Scenario 1 Reliability",
            ylab = "Sugden Reliability",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

# Scenario 2 # 
heatscatter(data_comb$ICC_Scen2,data_comb$ICC_Sugden,
            colpal="bl2gr2rd",
            main="Scenario 2",
            xlab = "Scenario 2 Reliability",
            ylab = "Sugden Reliability",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

# Scenario 3 # 
heatscatter(data_comb$ICC_Scen3,data_comb$ICC_Sugden,
            colpal="bl2gr2rd",
            main="Scenario 3",
            xlab = "Scenario 3 Reliability",
            ylab = "Sugden Reliability",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

# Scenario 4 # 
heatscatter(data_comb$ICC_Scen4,data_comb$ICC_Sugden,
            colpal="bl2gr2rd",
            main="Scenario 4",
            xlab = "Scenario 4 Reliability",
            ylab = "Sugden Reliability",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

# Scenario 5 # 
heatscatter(data_comb$ICC_Scen5,data_comb$ICC_Sugden,
            colpal="bl2gr2rd",
            main="Scenario 5",
            xlab = "Scenario 5 Reliability",
            ylab = "Sugden Reliability",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

# Scenario 6 # 
heatscatter(data_comb$ICC_Scen6,data_comb$ICC_Sugden,
            colpal="bl2gr2rd",
            main="Scenario 6",
            xlab = "Scenario 6 Reliability",
            ylab = "Sugden Reliability",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

# Scenario 7 # 
heatscatter(data_comb$ICC_Scen7,data_comb$ICC_Sugden,
            colpal="bl2gr2rd",
            main="Scenario 7",
            xlab = "Scenario 7 Reliability",
            ylab = "Sugden Reliability",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

# Scenario 8 # 
heatscatter(data_comb$ICC_Scen8,data_comb$ICC_Sugden,
            colpal="bl2gr2rd",
            main="Scenario 8",
            xlab = "Scenario 8 Reliability",
            ylab = "Sugden Reliability",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

# Scenario 9 # 
heatscatter(data_comb$ICC_Scen9,data_comb$ICC_Sugden,
            colpal="bl2gr2rd",
            main="Scenario 9",
            xlab = "Scenario 9 Reliability",
            ylab = "Sugden Reliability",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

# Scenario 10 # 
heatscatter(data_comb$ICC_Scen10,data_comb$ICC_Sugden,
            colpal="bl2gr2rd",
            main="Scenario 10",
            xlab = "Scenario 10 Reliability",
            ylab = "Sugden Reliability",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

### Compute ICC Stratified Correlation Values with Sugden Probes ### ----
data_comb_LowLowICC <- dplyr::filter(data_comb, data_comb$ICC_Sugden < 0.01)
data_cor_LowLowICC <- cor(data_comb_LowLowICC[,2:12])[,11]

data_comb_LowICC <- dplyr::filter(data_comb, data_comb$ICC_Sugden > 0.01 & data_comb$ICC_Sugden < 0.50)
data_cor_LowICC <- cor(data_comb_LowICC[,2:12])[,11]

data_comb_MedICC <- dplyr::filter(data_comb, data_comb$ICC_Sugden > 0.50 & data_comb$ICC_Sugden < 0.75)
data_cor_MedICC <- cor(data_comb_MedICC[,2:12])[,11]

data_comb_HighICC <- dplyr::filter(data_comb, data_comb$ICC_Sugden > 0.75 & data_comb$ICC_Sugden < 0.90)
data_cor_HighICC <- cor(data_comb_HighICC[,2:12])[,11]

data_comb_HighHighICC <- dplyr::filter(data_comb, data_comb$ICC_Sugden > 0.90)
data_cor_HighHighICC <- cor(data_comb_HighHighICC[,2:12])[,11]

SugdenCorrAll <- data.frame(LowLow = data_cor_LowLowICC,
                            Low = data_cor_LowICC,
                            Med = data_cor_MedICC,
                            High = data_cor_HighICC,
                            HighHigh = data_cor_HighHighICC)

write.csv(SugdenCorrAll, "Aim7_SugdenCorrAll.csv")
