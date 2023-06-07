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
### Read in Data ### ----
load("scen6a_icc_s1_t1t2_ela0_covar2_revision.RData")
icc.comb.scen6a.covar2 <- dplyr::filter(icc.comb.scen6a.covar2, icc.comb.scen6a.covar2$covariates == "Adjusted for DNAm_Monocytes, Methylation_Batch")
icc.comb.scen6a.covar2 <- dplyr::filter(icc.comb.scen6a.covar2, icc.comb.scen6a.covar2$type == "ICC(2,1)")[,c(1,5)]
colnames(icc.comb.scen6a.covar2)[2] <- "ICC_Scen6_ELA0"
load("scen6b_icc_s1_t1t2_ela1_covar2_revision.RData")
icc.comb.scen6b.covar2 <- dplyr::filter(icc.comb.scen6b.covar2, icc.comb.scen6b.covar2$covariates == "Adjusted for DNAm_Monocytes, Methylation_Batch")
icc.comb.scen6b.covar2 <- dplyr::filter(icc.comb.scen6b.covar2, icc.comb.scen6b.covar2$type == "ICC(2,1)")[,c(1,5)]
colnames(icc.comb.scen6b.covar2)[2] <- "ICC_Scen6_ELA1"

load("scen7a_icc_s1_t1t3_ela0_covar2_revision.RData")
icc.comb.scen7a.covar2 <- dplyr::filter(icc.comb.scen7a.covar2, icc.comb.scen7a.covar2$covariates == "Adjusted for DNAm_Monocytes, Methylation_Batch")
icc.comb.scen7a.covar2 <- dplyr::filter(icc.comb.scen7a.covar2, icc.comb.scen7a.covar2$type == "ICC(2,1)")[,c(1,5)]
colnames(icc.comb.scen7a.covar2)[2] <- "ICC_Scen7_ELA0"
load("scen7b_icc_s1_t1t3_ela1_covar2_revision.RData")
icc.comb.scen7b.covar2 <- dplyr::filter(icc.comb.scen7b.covar2, icc.comb.scen7b.covar2$covariates == "Adjusted for DNAm_Monocytes, Methylation_Batch")
icc.comb.scen7b.covar2 <- dplyr::filter(icc.comb.scen7b.covar2, icc.comb.scen7b.covar2$type == "ICC(2,1)")[,c(1,5)]
colnames(icc.comb.scen7b.covar2)[2] <- "ICC_Scen7_ELA1"

load("scen8a_icc_s1_t1t4_ela0_covar2_revision.RData")
icc.comb.scen8a.covar2 <- dplyr::filter(icc.comb.scen8a.covar2, icc.comb.scen8a.covar2$covariates == "Adjusted for DNAm_Monocytes, Methylation_Batch")
icc.comb.scen8a.covar2 <- dplyr::filter(icc.comb.scen8a.covar2, icc.comb.scen8a.covar2$type == "ICC(2,1)")[,c(1,5)]
colnames(icc.comb.scen8a.covar2)[2] <- "ICC_Scen8_ELA0"
load("scen8b_icc_s1_t1t4_ela1_covar2_revision.RData")
icc.comb.scen8b.covar2 <- dplyr::filter(icc.comb.scen8b.covar2, icc.comb.scen8b.covar2$covariates == "Adjusted for DNAm_Monocytes, Methylation_Batch")
icc.comb.scen8b.covar2 <- dplyr::filter(icc.comb.scen8b.covar2, icc.comb.scen8b.covar2$type == "ICC(2,1)")[,c(1,5)]
colnames(icc.comb.scen8b.covar2)[2] <- "ICC_Scen8_ELA1"

load("scen9a_icc_s1_t3t4_ela0_covar2_revision.RData")
icc.comb.scen9a.covar2 <- dplyr::filter(icc.comb.scen9a.covar2, icc.comb.scen9a.covar2$covariates == "Adjusted for DNAm_Monocytes, Methylation_Batch")
icc.comb.scen9a.covar2 <- dplyr::filter(icc.comb.scen9a.covar2, icc.comb.scen9a.covar2$type == "ICC(2,1)")[,c(1,5)]
colnames(icc.comb.scen9a.covar2)[2] <- "ICC_Scen9_ELA0"
load("scen9b_icc_s1_t3t4_ela1_covar2_revision.RData")
icc.comb.scen9b.covar2 <- dplyr::filter(icc.comb.scen9b.covar2, icc.comb.scen9b.covar2$covariates == "Adjusted for DNAm_Monocytes, Methylation_Batch")
icc.comb.scen9b.covar2 <- dplyr::filter(icc.comb.scen9b.covar2, icc.comb.scen9b.covar2$type == "ICC(2,1)")[,c(1,5)]
colnames(icc.comb.scen9b.covar2)[2] <- "ICC_Scen9_ELA1"

### Merge all Data ### ----
data <- merge(icc.comb.scen6a.covar2, icc.comb.scen6b.covar2, by = "DNAm.probe")
data <- merge(data, icc.comb.scen7a.covar2, by = "DNAm.probe")
data <- merge(data, icc.comb.scen7b.covar2, by = "DNAm.probe")
data <- merge(data, icc.comb.scen8a.covar2, by = "DNAm.probe")
data <- merge(data, icc.comb.scen8b.covar2, by = "DNAm.probe")
data <- merge(data, icc.comb.scen9a.covar2, by = "DNAm.probe")
data <- merge(data, icc.comb.scen9b.covar2, by = "DNAm.probe")
Exclude_Probes_Comb <- readRDS("Exclude_Probes_Comb.Rds")
keep <- !(data$DNAm.probe %in% Exclude_Probes_Comb$DNAm.probe)
data <- data[keep,]

# Delete old dataframes #
icc.comb.scen6a.covar2 <- NULL
icc.comb.scen6b.covar2 <- NULL
icc.comb.scen7a.covar2 <- NULL
icc.comb.scen7b.covar2 <- NULL
icc.comb.scen8a.covar2 <- NULL
icc.comb.scen8b.covar2 <- NULL
icc.comb.scen9a.covar2 <- NULL
icc.comb.scen9b.covar2 <- NULL


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
