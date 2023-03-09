### Results - Section 2.5 ### ----

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
################## Early-Life Adversity Calculations #######################
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
