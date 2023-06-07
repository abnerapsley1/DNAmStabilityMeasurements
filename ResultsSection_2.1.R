### Results - Section 2.1 ### ----

### Load Packages ### ----
library(dplyr)
library(psych)
library(ggplot2)
library(ggpubr)
library(dgof)
library(LSD)

############################################################################
###################### Type of ICC Value Calculations ######################
############################################################################

##################
### Scenario 1 ###
##################

### Read in Data ### ----
load("scen1_icc_s0_t1t2_ela01_revision.RData")
icc.comb.scen1 <- dplyr::filter(icc.comb.scen1, icc.comb.scen1$covariates == "No covariates")[,c(1,4:5)]
icc.comb.scen1_21 <- icc.comb.scen1[1:865859,]
icc.comb.scen1_2k <- icc.comb.scen1[865860:1731718,]
colnames(icc.comb.scen1_21)[3] <- "ICC(2,1)"
colnames(icc.comb.scen1_2k)[3] <- "ICC(2,k)"
data <- merge(icc.comb.scen1_21, icc.comb.scen1_2k, by = "DNAm.probe")
data <- data[,c(1,3,5)]
Exclude_Probes_Comb <- readRDS("Exclude_Probes_Comb.Rds")
keep <- !(data$DNAm.probe %in% Exclude_Probes_Comb$DNAm.probe)
data <- data[keep,]

### Perform Statistical Testing ### ----
# Paired T-test #
t.test(data$`ICC(2,1)`, data$`ICC(2,k)`, paired = TRUE)
# Proportional bias test #
data$diff <- data$`ICC(2,1)` - data$`ICC(2,k)`
data$mean <- (data$`ICC(2,1)` + data$`ICC(2,k)`)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$`ICC(2,1)`, data$`ICC(2,k)`)
# Chi-square test #
length(data$`ICC(2,1)`[data$`ICC(2,1)` <= 0.01])
length(data$`ICC(2,k)`[data$`ICC(2,k)` <= 0.01])
data_chisq <- data.frame(ICC_21 = c(63725, length(data$`ICC(2,1)`) - 63725),
                         ICC_2k = c(61287, length(data$`ICC(2,k)`) - 61287))
chisq.test(data_chisq)

### Plot Bland-Altman Comparison ### ----
heatscatter(data$`ICC(2,1)`,data$`ICC(2,k)`,
            colpal="bl2gr2rd",
            main = "",
            xlab = "ICC(2,1)",
            ylab = "ICC(2,k)",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

##################
### SCENARIO 2 ###
##################

### Read in Data ### ----
load("scen2_icc_s0_t1t3_ela01_revision.RData")
icc.comb.scen2 <- dplyr::filter(icc.comb.scen2, icc.comb.scen2$covariates == "No covariates")[,c(1,4:5)]
icc.comb.scen2_21 <- icc.comb.scen2[1:865859,]
icc.comb.scen2_2k <- icc.comb.scen2[865860:1731718,]
colnames(icc.comb.scen2_21)[3] <- "ICC(2,1)"
colnames(icc.comb.scen2_2k)[3] <- "ICC(2,k)"
data <- merge(icc.comb.scen2_21, icc.comb.scen2_2k, by = "DNAm.probe")
data <- data[,c(1,3,5)]
Exclude_Probes_Comb <- readRDS("Exclude_Probes_Comb.Rds")
keep <- !(data$DNAm.probe %in% Exclude_Probes_Comb$DNAm.probe)
data <- data[keep,]

### Perform Statistical Testing ### ----
# Paired T-test #
t.test(data$`ICC(2,1)`, data$`ICC(2,k)`, paired = TRUE)
# Proportional bias test #
data$diff <- data$`ICC(2,1)` - data$`ICC(2,k)`
data$mean <- (data$`ICC(2,1)` + data$`ICC(2,k)`)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$`ICC(2,1)`, data$`ICC(2,k)`)
# Chi-square test #
length(data$`ICC(2,1)`[data$`ICC(2,1)` <= 0.01])
length(data$`ICC(2,k)`[data$`ICC(2,k)` <= 0.01])
data_chisq <- data.frame(ICC_21 = c(63725, length(data$`ICC(2,1)`) - 63725),
                         ICC_2k = c(61287, length(data$`ICC(2,k)`) - 61287))
chisq.test(data_chisq)

### Plot Bland-Altman Comparison ### ----
heatscatter(data$mean,data$diff,
            colpal="bl2gr2rd",
            main="Scenario 2",
            xlab = "Mean Reliability",
            ylab = "ICC(2,1) - ICC(2,k)",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

##################
### SCENARIO 3 ###
##################

### Read in Data ### ----
load("scen3_icc_s0_t1t4_ela01_revision.RData")
icc.comb.scen3 <- dplyr::filter(icc.comb.scen3, icc.comb.scen3$covariates == "No covariates")[,c(1,4:5)]
icc.comb.scen3_21 <- icc.comb.scen3[1:865859,]
icc.comb.scen3_2k <- icc.comb.scen3[865860:1731718,]
colnames(icc.comb.scen3_21)[3] <- "ICC(2,1)"
colnames(icc.comb.scen3_2k)[3] <- "ICC(2,k)"
data <- merge(icc.comb.scen3_21, icc.comb.scen3_2k, by = "DNAm.probe")
data <- data[,c(1,3,5)]
Exclude_Probes_Comb <- readRDS("Exclude_Probes_Comb.Rds")
keep <- !(data$DNAm.probe %in% Exclude_Probes_Comb$DNAm.probe)
data <- data[keep,]

### Perform Statistical Testing ### ----
# Paired T-test #
t.test(data$`ICC(2,1)`, data$`ICC(2,k)`, paired = TRUE)
# Proportional bias test #
data$diff <- data$`ICC(2,1)` - data$`ICC(2,k)`
data$mean <- (data$`ICC(2,1)` + data$`ICC(2,k)`)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$`ICC(2,1)`, data$`ICC(2,k)`)
# Chi-square test #
length(data$`ICC(2,1)`[data$`ICC(2,1)` <= 0.01])
length(data$`ICC(2,k)`[data$`ICC(2,k)` <= 0.01])
data_chisq <- data.frame(ICC_21 = c(63725, length(data$`ICC(2,1)`) - 63725),
                         ICC_2k = c(61287, length(data$`ICC(2,k)`) - 61287))
chisq.test(data_chisq)

### Plot Bland-Altman Comparison ### ----
heatscatter(data$mean,data$diff,
            colpal="bl2gr2rd",
            main="Scenario 3",
            xlab = "Mean Reliability",
            ylab = "ICC(2,1) - ICC(2,k)",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

##################
### SCENARIO 4 ###
##################

### Read in Data ### ----
load("scen4_icc_s0_t3t4_ela01_revision.RData")
icc.comb.scen4 <- dplyr::filter(icc.comb.scen4, icc.comb.scen4$covariates == "No covariates")[,c(1,4:5)]
icc.comb.scen4_21 <- icc.comb.scen4[1:865859,]
icc.comb.scen4_2k <- icc.comb.scen4[865860:1731718,]
colnames(icc.comb.scen4_21)[3] <- "ICC(2,1)"
colnames(icc.comb.scen4_2k)[3] <- "ICC(2,k)"
data <- merge(icc.comb.scen4_21, icc.comb.scen4_2k, by = "DNAm.probe")
data <- data[,c(1,3,5)]
Exclude_Probes_Comb <- readRDS("Exclude_Probes_Comb.Rds")
keep <- !(data$DNAm.probe %in% Exclude_Probes_Comb$DNAm.probe)
data <- data[keep,]

### Perform Statistical Testing ### ----
# Paired T-test #
t.test(data$`ICC(2,1)`, data$`ICC(2,k)`, paired = TRUE)
# Proportional bias test #
data$diff <- data$`ICC(2,1)` - data$`ICC(2,k)`
data$mean <- (data$`ICC(2,1)` + data$`ICC(2,k)`)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$`ICC(2,1)`, data$`ICC(2,k)`)
# Chi-square test #
length(data$`ICC(2,1)`[data$`ICC(2,1)` <= 0.01])
length(data$`ICC(2,k)`[data$`ICC(2,k)` <= 0.01])
data_chisq <- data.frame(ICC_21 = c(63725, length(data$`ICC(2,1)`) - 63725),
                         ICC_2k = c(61287, length(data$`ICC(2,k)`) - 61287))
chisq.test(data_chisq)

### Plot Bland-Altman Comparison ### ----
heatscatter(data$mean,data$diff,
            colpal="bl2gr2rd",
            main="Scenario 4",
            xlab = "Mean Reliability",
            ylab = "ICC(2,1) - ICC(2,k)",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

##################
### SCENARIO 5 ###
##################

### Read in Data ### ----
load("scen5_icc_s0_t1t2t3t4_ela01_revision.RData")
icc.comb.scen5 <- dplyr::filter(icc.comb.scen5, icc.comb.scen5$covariates == "No covariates")[,c(1,4:5)]
icc.comb.scen5_21 <- icc.comb.scen5[1:865859,]
icc.comb.scen5_2k <- icc.comb.scen5[865860:1731718,]
colnames(icc.comb.scen5_21)[3] <- "ICC(2,1)"
colnames(icc.comb.scen5_2k)[3] <- "ICC(2,k)"
data <- merge(icc.comb.scen5_21, icc.comb.scen5_2k, by = "DNAm.probe")
data <- data[,c(1,3,5)]
Exclude_Probes_Comb <- readRDS("Exclude_Probes_Comb.Rds")
keep <- !(data$DNAm.probe %in% Exclude_Probes_Comb$DNAm.probe)
data <- data[keep,]

### Perform Statistical Testing ### ----
# Paired T-test #
t.test(data$`ICC(2,1)`, data$`ICC(2,k)`, paired = TRUE)
# Proportional bias test #
data$diff <- data$`ICC(2,1)` - data$`ICC(2,k)`
data$mean <- (data$`ICC(2,1)` + data$`ICC(2,k)`)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$`ICC(2,1)`, data$`ICC(2,k)`)
# Chi-square test #
length(data$`ICC(2,1)`[data$`ICC(2,1)` <= 0.01])
length(data$`ICC(2,k)`[data$`ICC(2,k)` <= 0.01])
data_chisq <- data.frame(ICC_21 = c(63725, length(data$`ICC(2,1)`) - 63725),
                         ICC_2k = c(61287, length(data$`ICC(2,k)`) - 61287))
chisq.test(data_chisq)

### Plot Bland-Altman Comparison ### ----
heatscatter(data$mean,data$diff,
            colpal="bl2gr2rd",
            main="Scenario 5",
            xlab = "Mean Reliability",
            ylab = "ICC(2,1) - ICC(2,k)",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

##################
### SCENARIO 6 ###
##################

### Read in Data ### ----
load("scen6sampled_icc_s1_t1t2_ela01_revision.RData")
icc.comb.scen6sampled <- dplyr::filter(icc.comb.scen6sampled, icc.comb.scen6sampled$covariates == "No covariates")[,c(1,4:5)]
icc.comb.scen6sampled_21 <- icc.comb.scen6sampled[1:865859,]
icc.comb.scen6sampled_2k <- icc.comb.scen6sampled[865860:1731718,]
colnames(icc.comb.scen6sampled_21)[3] <- "ICC(2,1)"
colnames(icc.comb.scen6sampled_2k)[3] <- "ICC(2,k)"
data <- merge(icc.comb.scen6sampled_21, icc.comb.scen6sampled_2k, by = "DNAm.probe")
data <- data[,c(1,3,5)]
Exclude_Probes_Comb <- readRDS("Exclude_Probes_Comb.Rds")
keep <- !(data$DNAm.probe %in% Exclude_Probes_Comb$DNAm.probe)
data <- data[keep,]

### Perform Statistical Testing ### ----
# Paired T-test #
t.test(data$`ICC(2,1)`, data$`ICC(2,k)`, paired = TRUE)
# Proportional bias test #
data$diff <- data$`ICC(2,1)` - data$`ICC(2,k)`
data$mean <- (data$`ICC(2,1)` + data$`ICC(2,k)`)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$`ICC(2,1)`, data$`ICC(2,k)`)
# Chi-square test #
length(data$`ICC(2,1)`[data$`ICC(2,1)` <= 0.01])
length(data$`ICC(2,k)`[data$`ICC(2,k)` <= 0.01])
data_chisq <- data.frame(ICC_21 = c(63725, length(data$`ICC(2,1)`) - 63725),
                         ICC_2k = c(61287, length(data$`ICC(2,k)`) - 61287))
chisq.test(data_chisq)

### Plot Bland-Altman Comparison ### ----
heatscatter(data$mean,data$diff,
            colpal="bl2gr2rd",
            main="Scenario 6",
            xlab = "Mean Reliability",
            ylab = "ICC(2,1) - ICC(2,k)",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

##################
### SCENARIO 7 ###
##################

### Read in Data ### ----
load("scen7sampled_icc_s1_t1t3_ela01_revision.RData")
icc.comb.scen7sampled <- dplyr::filter(icc.comb.scen7sampled, icc.comb.scen7sampled$covariates == "No covariates")[,c(1,4:5)]
icc.comb.scen7sampled_21 <- icc.comb.scen7sampled[1:865859,]
icc.comb.scen7sampled_2k <- icc.comb.scen7sampled[865860:1731718,]
colnames(icc.comb.scen7sampled_21)[3] <- "ICC(2,1)"
colnames(icc.comb.scen7sampled_2k)[3] <- "ICC(2,k)"
data <- merge(icc.comb.scen7sampled_21, icc.comb.scen7sampled_2k, by = "DNAm.probe")
data <- data[,c(1,3,5)]
Exclude_Probes_Comb <- readRDS("Exclude_Probes_Comb.Rds")
keep <- !(data$DNAm.probe %in% Exclude_Probes_Comb$DNAm.probe)
data <- data[keep,]

### Perform Statistical Testing ### ----
# Paired T-test #
t.test(data$`ICC(2,1)`, data$`ICC(2,k)`, paired = TRUE)
# Proportional bias test #
data$diff <- data$`ICC(2,1)` - data$`ICC(2,k)`
data$mean <- (data$`ICC(2,1)` + data$`ICC(2,k)`)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$`ICC(2,1)`, data$`ICC(2,k)`)
# Chi-square test #
length(data$`ICC(2,1)`[data$`ICC(2,1)` <= 0.01])
length(data$`ICC(2,k)`[data$`ICC(2,k)` <= 0.01])
data_chisq <- data.frame(ICC_21 = c(63725, length(data$`ICC(2,1)`) - 63725),
                         ICC_2k = c(61287, length(data$`ICC(2,k)`) - 61287))
chisq.test(data_chisq)

### Plot Bland-Altman Comparison ### ----
heatscatter(data$mean,data$diff,
            colpal="bl2gr2rd",
            main="Scenario 7",
            xlab = "Mean Reliability",
            ylab = "ICC(2,1) - ICC(2,k)",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

##################
### SCENARIO 8 ###
##################

### Read in Data ### ----
load("scen8sampled_icc_s1_t1t4_ela01_revision.RData")
icc.comb.scen8sampled <- dplyr::filter(icc.comb.scen8sampled, icc.comb.scen8sampled$covariates == "No covariates")[,c(1,4:5)]
icc.comb.scen8sampled_21 <- icc.comb.scen8sampled[1:865859,]
icc.comb.scen8sampled_2k <- icc.comb.scen8sampled[865860:1731718,]
colnames(icc.comb.scen8sampled_21)[3] <- "ICC(2,1)"
colnames(icc.comb.scen8sampled_2k)[3] <- "ICC(2,k)"
data <- merge(icc.comb.scen8sampled_21, icc.comb.scen8sampled_2k, by = "DNAm.probe")
data <- data[,c(1,3,5)]
Exclude_Probes_Comb <- readRDS("Exclude_Probes_Comb.Rds")
keep <- !(data$DNAm.probe %in% Exclude_Probes_Comb$DNAm.probe)
data <- data[keep,]

### Perform Statistical Testing ### ----
# Paired T-test #
t.test(data$`ICC(2,1)`, data$`ICC(2,k)`, paired = TRUE)
# Proportional bias test #
data$diff <- data$`ICC(2,1)` - data$`ICC(2,k)`
data$mean <- (data$`ICC(2,1)` + data$`ICC(2,k)`)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$`ICC(2,1)`, data$`ICC(2,k)`)
# Chi-square test #
length(data$`ICC(2,1)`[data$`ICC(2,1)` <= 0.01])
length(data$`ICC(2,k)`[data$`ICC(2,k)` <= 0.01])
data_chisq <- data.frame(ICC_21 = c(63725, length(data$`ICC(2,1)`) - 63725),
                         ICC_2k = c(61287, length(data$`ICC(2,k)`) - 61287))
chisq.test(data_chisq)

### Plot Bland-Altman Comparison ### ----
heatscatter(data$mean,data$diff,
            colpal="bl2gr2rd",
            main="Scenario 8",
            xlab = "Mean Reliability",
            ylab = "ICC(2,1) - ICC(2,k)",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

##################
### SCENARIO 9 ###
##################

### Read in Data ### ----
load("scen9sampled_icc_s1_t3t4_ela01_revision.RData")
icc.comb.scen9sampled <- dplyr::filter(icc.comb.scen9sampled, icc.comb.scen9sampled$covariates == "No covariates")[,c(1,4:5)]
icc.comb.scen9sampled_21 <- icc.comb.scen9sampled[1:865859,]
icc.comb.scen9sampled_2k <- icc.comb.scen9sampled[865860:1731718,]
colnames(icc.comb.scen9sampled_21)[3] <- "ICC(2,1)"
colnames(icc.comb.scen9sampled_2k)[3] <- "ICC(2,k)"
data <- merge(icc.comb.scen9sampled_21, icc.comb.scen9sampled_2k, by = "DNAm.probe")
data <- data[,c(1,3,5)]
Exclude_Probes_Comb <- readRDS("Exclude_Probes_Comb.Rds")
keep <- !(data$DNAm.probe %in% Exclude_Probes_Comb$DNAm.probe)
data <- data[keep,]

### Perform Statistical Testing ### ----
# Paired T-test #
t.test(data$`ICC(2,1)`, data$`ICC(2,k)`, paired = TRUE)
# Proportional bias test #
data$diff <- data$`ICC(2,1)` - data$`ICC(2,k)`
data$mean <- (data$`ICC(2,1)` + data$`ICC(2,k)`)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$`ICC(2,1)`, data$`ICC(2,k)`)
# Chi-square test #
length(data$`ICC(2,1)`[data$`ICC(2,1)` <= 0.01])
length(data$`ICC(2,k)`[data$`ICC(2,k)` <= 0.01])
data_chisq <- data.frame(ICC_21 = c(63725, length(data$`ICC(2,1)`) - 63725),
                         ICC_2k = c(61287, length(data$`ICC(2,k)`) - 61287))
chisq.test(data_chisq)

### Plot Bland-Altman Comparison ### ----
heatscatter(data$mean,data$diff,
            colpal="bl2gr2rd",
            main="Scenario 9",
            xlab = "Mean Reliability",
            ylab = "ICC(2,1) - ICC(2,k)",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

###################
### SCENARIO 10 ###
###################

### Read in Data ### ----
load("scen10_icc_s0s1_t1_ela01_revision.RData")
icc.comb.scen10 <- dplyr::filter(icc.comb.scen10, icc.comb.scen10$covariates == "No covariates")[,c(1,4:5)]
icc.comb.scen10_21 <- icc.comb.scen10[1:865859,]
icc.comb.scen10_2k <- icc.comb.scen10[865860:1731718,]
colnames(icc.comb.scen10_21)[3] <- "ICC(2,1)"
colnames(icc.comb.scen10_2k)[3] <- "ICC(2,k)"
data <- merge(icc.comb.scen10_21, icc.comb.scen10_2k, by = "DNAm.probe")
data <- data[,c(1,3,5)]
Exclude_Probes_Comb <- readRDS("Exclude_Probes_Comb.Rds")
keep <- !(data$DNAm.probe %in% Exclude_Probes_Comb$DNAm.probe)
data <- data[keep,]

### Perform Statistical Testing ### ----
# Paired T-test #
t.test(data$`ICC(2,1)`, data$`ICC(2,k)`, paired = TRUE)
# Proportional bias test #
data$diff <- data$`ICC(2,1)` - data$`ICC(2,k)`
data$mean <- (data$`ICC(2,1)` + data$`ICC(2,k)`)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$`ICC(2,1)`, data$`ICC(2,k)`)
# Chi-square test #
length(data$`ICC(2,1)`[data$`ICC(2,1)` <= 0.01])
length(data$`ICC(2,k)`[data$`ICC(2,k)` <= 0.01])
data_chisq <- data.frame(ICC_21 = c(63725, length(data$`ICC(2,1)`) - 63725),
                         ICC_2k = c(61287, length(data$`ICC(2,k)`) - 61287))
chisq.test(data_chisq)

### Plot Bland-Altman Comparison ### ----
heatscatter(data$mean,data$diff,
            colpal="bl2gr2rd",
            main="Scenario 10",
            xlab = "Mean Reliability",
            ylab = "ICC(2,1) - ICC(2,k)",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

############################################################################
################### Immune Cell Proportion Calculations ####################
############################################################################

##################
### Scenario 1 ###
##################

### Read in Data ### ----
load("scen1_icc_s0_t1t2_ela01_revision.RData")
icc.comb.scen1 <- dplyr::filter(icc.comb.scen1, icc.comb.scen1$type == "ICC(2,1)")[,c(1,2,5)]
icc.comb.scen1_NoAdj <- icc.comb.scen1[1:865859,]
icc.comb.scen1_Adj <- icc.comb.scen1[865860:1731718,]
colnames(icc.comb.scen1_NoAdj)[3] <- "ICC_NoAdj"
colnames(icc.comb.scen1_Adj)[3] <- "ICC_Adj"
data <- merge(icc.comb.scen1_NoAdj, icc.comb.scen1_Adj, by = "DNAm.probe")
data <- data[,c(1,3,5)]
Exclude_Probes_Comb <- readRDS("Exclude_Probes_Comb.Rds")
keep <- !(data$DNAm.probe %in% Exclude_Probes_Comb$DNAm.probe)
data <- data[keep,]

### Perform Statistical Testing ### ----
# Paired T-test #
t.test(data$ICC_NoAdj, data$ICC_Adj, paired = TRUE)
# Proportional bias test #
data$diff <- data$ICC_NoAdj - data$ICC_Adj
data$mean <- (data$ICC_NoAdj + data$ICC_Adj)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$ICC_NoAdj, data$ICC_Adj)
# Chi-square test #
length(data$ICC_NoAdj[data$ICC_NoAdj <= 0.01])
length(data$ICC_Adj[data$ICC_Adj <= 0.01])
data_chisq <- data.frame(ICC_NoAdj = c(63725, length(data$ICC_NoAdj) - 63725),
                           ICC_Adj = c(59371, length(data$ICC_Adj) - 59371))
chisq.test(data_chisq)

### Plot Bland-Altman Comparison ### ----
heatscatter(data$ICC_NoAdj,data$ICC_Adj,
            colpal="bl2gr2rd",
            main="",
            xlab = "No Immune Cell Adjustment",
            ylab = "Immune Cell Adjustment",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

##################
### Scenario 2 ###
##################

### Read in Data ### ----
load("scen2_icc_s0_t1t3_ela01_revision.RData")
icc.comb.scen2 <- dplyr::filter(icc.comb.scen2, icc.comb.scen2$type == "ICC(2,1)")[,c(1,2,5)]
icc.comb.scen2_NoAdj <- icc.comb.scen2[1:865859,]
icc.comb.scen2_Adj <- icc.comb.scen2[865860:1731718,]
colnames(icc.comb.scen2_NoAdj)[3] <- "ICC_NoAdj"
colnames(icc.comb.scen2_Adj)[3] <- "ICC_Adj"
data <- merge(icc.comb.scen2_NoAdj, icc.comb.scen2_Adj, by = "DNAm.probe")
data <- data[,c(1,3,5)]
Exclude_Probes_Comb <- readRDS("Exclude_Probes_Comb.Rds")
keep <- !(data$DNAm.probe %in% Exclude_Probes_Comb$DNAm.probe)
data <- data[keep,]

### Perform Statistical Testing ### ----
# Paired T-test #
t.test(data$ICC_NoAdj, data$ICC_Adj, paired = TRUE)
# Proportional bias test #
data$diff <- data$ICC_NoAdj - data$ICC_Adj
data$mean <- (data$ICC_NoAdj + data$ICC_Adj)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$ICC_NoAdj, data$ICC_Adj)
# Chi-square test #
length(data$ICC_NoAdj[data$ICC_NoAdj <= 0.01])
length(data$ICC_Adj[data$ICC_Adj <= 0.01])
data_chisq <- data.frame(ICC_NoAdj = c(63725, length(data$ICC_NoAdj) - 63725),
                         ICC_Adj = c(59371, length(data$ICC_Adj) - 59371))
chisq.test(data_chisq)

### Plot Bland-Altman Comparison ### ----
heatscatter(data$ICC_NoAdj,data$ICC_Adj,
            colpal="bl2gr2rd",
            main="",
            xlab = "No Immune Cell Adjustment",
            ylab = "Immune Cell Adjustment",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

##################
### Scenario 3 ###
##################

### Read in Data ### ----
load("scen3_icc_s0_t1t4_ela01_revision.RData")
icc.comb.scen3 <- dplyr::filter(icc.comb.scen3, icc.comb.scen3$type == "ICC(2,1)")[,c(1,2,5)]
icc.comb.scen3_NoAdj <- icc.comb.scen3[1:865859,]
icc.comb.scen3_Adj <- icc.comb.scen3[865860:1731718,]
colnames(icc.comb.scen3_NoAdj)[3] <- "ICC_NoAdj"
colnames(icc.comb.scen3_Adj)[3] <- "ICC_Adj"
data <- merge(icc.comb.scen3_NoAdj, icc.comb.scen3_Adj, by = "DNAm.probe")
data <- data[,c(1,3,5)]
Exclude_Probes_Comb <- readRDS("Exclude_Probes_Comb.Rds")
keep <- !(data$DNAm.probe %in% Exclude_Probes_Comb$DNAm.probe)
data <- data[keep,]

### Perform Statistical Testing ### ----
# Paired T-test #
t.test(data$ICC_NoAdj, data$ICC_Adj, paired = TRUE)
# Proportional bias test #
data$diff <- data$ICC_NoAdj - data$ICC_Adj
data$mean <- (data$ICC_NoAdj + data$ICC_Adj)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$ICC_NoAdj, data$ICC_Adj)
# Chi-square test #
length(data$ICC_NoAdj[data$ICC_NoAdj <= 0.01])
length(data$ICC_Adj[data$ICC_Adj <= 0.01])
data_chisq <- data.frame(ICC_NoAdj = c(63725, length(data$ICC_NoAdj) - 63725),
                         ICC_Adj = c(59371, length(data$ICC_Adj) - 59371))
chisq.test(data_chisq)

### Plot Bland-Altman Comparison ### ----
heatscatter(data$ICC_NoAdj,data$ICC_Adj,
            colpal="bl2gr2rd",
            main="",
            xlab = "No Immune Cell Adjustment",
            ylab = "Immune Cell Adjustment",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

##################
### Scenario 4 ###
##################

### Read in Data ### ----
load("scen4_icc_s0_t3t4_ela01_revision.RData")
icc.comb.scen4 <- dplyr::filter(icc.comb.scen4, icc.comb.scen4$type == "ICC(2,1)")[,c(1,2,5)]
icc.comb.scen4_NoAdj <- icc.comb.scen4[1:865859,]
icc.comb.scen4_Adj <- icc.comb.scen4[865860:1731718,]
colnames(icc.comb.scen4_NoAdj)[3] <- "ICC_NoAdj"
colnames(icc.comb.scen4_Adj)[3] <- "ICC_Adj"
data <- merge(icc.comb.scen4_NoAdj, icc.comb.scen4_Adj, by = "DNAm.probe")
data <- data[,c(1,3,5)]
Exclude_Probes_Comb <- readRDS("Exclude_Probes_Comb.Rds")
keep <- !(data$DNAm.probe %in% Exclude_Probes_Comb$DNAm.probe)
data <- data[keep,]

### Perform Statistical Testing ### ----
# Paired T-test #
t.test(data$ICC_NoAdj, data$ICC_Adj, paired = TRUE)
# Proportional bias test #
data$diff <- data$ICC_NoAdj - data$ICC_Adj
data$mean <- (data$ICC_NoAdj + data$ICC_Adj)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$ICC_NoAdj, data$ICC_Adj)
# Chi-square test #
length(data$ICC_NoAdj[data$ICC_NoAdj <= 0.01])
length(data$ICC_Adj[data$ICC_Adj <= 0.01])
data_chisq <- data.frame(ICC_NoAdj = c(63725, length(data$ICC_NoAdj) - 63725),
                         ICC_Adj = c(59371, length(data$ICC_Adj) - 59371))
chisq.test(data_chisq)

### Plot Bland-Altman Comparison ### ----
heatscatter(data$ICC_NoAdj,data$ICC_Adj,
            colpal="bl2gr2rd",
            main="",
            xlab = "No Immune Cell Adjustment",
            ylab = "Immune Cell Adjustment",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

##################
### Scenario 5 ###
##################

### Read in Data ### ----
load("scen5_icc_s0_t1t2t3t4_ela01_revision.RData")
icc.comb.scen5 <- dplyr::filter(icc.comb.scen5, icc.comb.scen5$type == "ICC(2,1)")[,c(1,2,5)]
icc.comb.scen5_NoAdj <- icc.comb.scen5[1:865859,]
icc.comb.scen5_Adj <- icc.comb.scen5[865860:1731718,]
colnames(icc.comb.scen5_NoAdj)[3] <- "ICC_NoAdj"
colnames(icc.comb.scen5_Adj)[3] <- "ICC_Adj"
data <- merge(icc.comb.scen5_NoAdj, icc.comb.scen5_Adj, by = "DNAm.probe")
data <- data[,c(1,3,5)]
Exclude_Probes_Comb <- readRDS("Exclude_Probes_Comb.Rds")
keep <- !(data$DNAm.probe %in% Exclude_Probes_Comb$DNAm.probe)
data <- data[keep,]

### Perform Statistical Testing ### ----
# Paired T-test #
t.test(data$ICC_NoAdj, data$ICC_Adj, paired = TRUE)
# Proportional bias test #
data$diff <- data$ICC_NoAdj - data$ICC_Adj
data$mean <- (data$ICC_NoAdj + data$ICC_Adj)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$ICC_NoAdj, data$ICC_Adj)
# Chi-square test #
length(data$ICC_NoAdj[data$ICC_NoAdj <= 0.01])
length(data$ICC_Adj[data$ICC_Adj <= 0.01])
data_chisq <- data.frame(ICC_NoAdj = c(63725, length(data$ICC_NoAdj) - 63725),
                         ICC_Adj = c(59371, length(data$ICC_Adj) - 59371))
chisq.test(data_chisq)

### Plot Bland-Altman Comparison ### ----
heatscatter(data$ICC_NoAdj,data$ICC_Adj,
            colpal="bl2gr2rd",
            main="",
            xlab = "No Immune Cell Adjustment",
            ylab = "Immune Cell Adjustment",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

##################
### Scenario 6 ###
##################

### Read in Data ### ----
load("scen6sampled_icc_s1_t1t2_ela01_revision.RData")
icc.comb.scen6sampled <- dplyr::filter(icc.comb.scen6sampled, icc.comb.scen6sampled$type == "ICC(2,1)")[,c(1,2,5)]
icc.comb.scen6sampled_NoAdj <- icc.comb.scen6sampled[1:865859,]
icc.comb.scen6sampled_Adj <- icc.comb.scen6sampled[865860:1731718,]
colnames(icc.comb.scen6sampled_NoAdj)[3] <- "ICC_NoAdj"
colnames(icc.comb.scen6sampled_Adj)[3] <- "ICC_Adj"
data <- merge(icc.comb.scen6sampled_NoAdj, icc.comb.scen6sampled_Adj, by = "DNAm.probe")
data <- data[,c(1,3,5)]
Exclude_Probes_Comb <- readRDS("Exclude_Probes_Comb.Rds")
keep <- !(data$DNAm.probe %in% Exclude_Probes_Comb$DNAm.probe)
data <- data[keep,]

### Perform Statistical Testing ### ----
# Paired T-test #
t.test(data$ICC_NoAdj, data$ICC_Adj, paired = TRUE)
# Proportional bias test #
data$diff <- data$ICC_NoAdj - data$ICC_Adj
data$mean <- (data$ICC_NoAdj + data$ICC_Adj)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$ICC_NoAdj, data$ICC_Adj)
# Chi-square test #
length(data$ICC_NoAdj[data$ICC_NoAdj <= 0.01])
length(data$ICC_Adj[data$ICC_Adj <= 0.01])
data_chisq <- data.frame(ICC_NoAdj = c(63725, length(data$ICC_NoAdj) - 63725),
                         ICC_Adj = c(59371, length(data$ICC_Adj) - 59371))
chisq.test(data_chisq)

### Plot Bland-Altman Comparison ### ----
heatscatter(data$ICC_NoAdj,data$ICC_Adj,
            colpal="bl2gr2rd",
            main="",
            xlab = "No Immune Cell Adjustment",
            ylab = "Immune Cell Adjustment",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

##################
### Scenario 7 ###
##################

### Read in Data ### ----
load("scen7sampled_icc_s1_t1t3_ela01_revision.RData")
icc.comb.scen7sampled <- dplyr::filter(icc.comb.scen7sampled, icc.comb.scen7sampled$type == "ICC(2,1)")[,c(1,2,5)]
icc.comb.scen7sampled_NoAdj <- icc.comb.scen7sampled[1:865859,]
icc.comb.scen7sampled_Adj <- icc.comb.scen7sampled[865860:1731718,]
colnames(icc.comb.scen7sampled_NoAdj)[3] <- "ICC_NoAdj"
colnames(icc.comb.scen7sampled_Adj)[3] <- "ICC_Adj"
data <- merge(icc.comb.scen7sampled_NoAdj, icc.comb.scen7sampled_Adj, by = "DNAm.probe")
data <- data[,c(1,3,5)]
Exclude_Probes_Comb <- readRDS("Exclude_Probes_Comb.Rds")
keep <- !(data$DNAm.probe %in% Exclude_Probes_Comb$DNAm.probe)
data <- data[keep,]

### Perform Statistical Testing ### ----
# Paired T-test #
t.test(data$ICC_NoAdj, data$ICC_Adj, paired = TRUE)
# Proportional bias test #
data$diff <- data$ICC_NoAdj - data$ICC_Adj
data$mean <- (data$ICC_NoAdj + data$ICC_Adj)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$ICC_NoAdj, data$ICC_Adj)
# Chi-square test #
length(data$ICC_NoAdj[data$ICC_NoAdj <= 0.01])
length(data$ICC_Adj[data$ICC_Adj <= 0.01])
data_chisq <- data.frame(ICC_NoAdj = c(63725, length(data$ICC_NoAdj) - 63725),
                         ICC_Adj = c(59371, length(data$ICC_Adj) - 59371))
chisq.test(data_chisq)

### Plot Bland-Altman Comparison ### ----
heatscatter(data$ICC_NoAdj,data$ICC_Adj,
            colpal="bl2gr2rd",
            main="",
            xlab = "No Immune Cell Adjustment",
            ylab = "Immune Cell Adjustment",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

##################
### Scenario 8 ###
##################

### Read in Data ### ----
load("scen8sampled_icc_s1_t1t4_ela01_revision.RData")
icc.comb.scen8sampled <- dplyr::filter(icc.comb.scen8sampled, icc.comb.scen8sampled$type == "ICC(2,1)")[,c(1,2,5)]
icc.comb.scen8sampled_NoAdj <- icc.comb.scen8sampled[1:865859,]
icc.comb.scen8sampled_Adj <- icc.comb.scen8sampled[865860:1731718,]
colnames(icc.comb.scen8sampled_NoAdj)[3] <- "ICC_NoAdj"
colnames(icc.comb.scen8sampled_Adj)[3] <- "ICC_Adj"
data <- merge(icc.comb.scen8sampled_NoAdj, icc.comb.scen8sampled_Adj, by = "DNAm.probe")
data <- data[,c(1,3,5)]
Exclude_Probes_Comb <- readRDS("Exclude_Probes_Comb.Rds")
keep <- !(data$DNAm.probe %in% Exclude_Probes_Comb$DNAm.probe)
data <- data[keep,]

### Perform Statistical Testing ### ----
# Paired T-test #
t.test(data$ICC_NoAdj, data$ICC_Adj, paired = TRUE)
# Proportional bias test #
data$diff <- data$ICC_NoAdj - data$ICC_Adj
data$mean <- (data$ICC_NoAdj + data$ICC_Adj)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$ICC_NoAdj, data$ICC_Adj)
# Chi-square test #
length(data$ICC_NoAdj[data$ICC_NoAdj <= 0.01])
length(data$ICC_Adj[data$ICC_Adj <= 0.01])
data_chisq <- data.frame(ICC_NoAdj = c(63725, length(data$ICC_NoAdj) - 63725),
                         ICC_Adj = c(59371, length(data$ICC_Adj) - 59371))
chisq.test(data_chisq)

### Plot Bland-Altman Comparison ### ----
heatscatter(data$ICC_NoAdj,data$ICC_Adj,
            colpal="bl2gr2rd",
            main="",
            xlab = "No Immune Cell Adjustment",
            ylab = "Immune Cell Adjustment",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

##################
### Scenario 9 ###
##################

### Read in Data ### ----
load("scen9sampled_icc_s1_t3t4_ela01_revision.RData")
icc.comb.scen9sampled <- dplyr::filter(icc.comb.scen9sampled, icc.comb.scen9sampled$type == "ICC(2,1)")[,c(1,2,5)]
icc.comb.scen9sampled_NoAdj <- icc.comb.scen9sampled[1:865859,]
icc.comb.scen9sampled_Adj <- icc.comb.scen9sampled[865860:1731718,]
colnames(icc.comb.scen9sampled_NoAdj)[3] <- "ICC_NoAdj"
colnames(icc.comb.scen9sampled_Adj)[3] <- "ICC_Adj"
data <- merge(icc.comb.scen9sampled_NoAdj, icc.comb.scen9sampled_Adj, by = "DNAm.probe")
data <- data[,c(1,3,5)]
Exclude_Probes_Comb <- readRDS("Exclude_Probes_Comb.Rds")
keep <- !(data$DNAm.probe %in% Exclude_Probes_Comb$DNAm.probe)
data <- data[keep,]

### Perform Statistical Testing ### ----
# Paired T-test #
t.test(data$ICC_NoAdj, data$ICC_Adj, paired = TRUE)
# Proportional bias test #
data$diff <- data$ICC_NoAdj - data$ICC_Adj
data$mean <- (data$ICC_NoAdj + data$ICC_Adj)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$ICC_NoAdj, data$ICC_Adj)
# Chi-square test #
length(data$ICC_NoAdj[data$ICC_NoAdj <= 0.01])
length(data$ICC_Adj[data$ICC_Adj <= 0.01])
data_chisq <- data.frame(ICC_NoAdj = c(63725, length(data$ICC_NoAdj) - 63725),
                         ICC_Adj = c(59371, length(data$ICC_Adj) - 59371))
chisq.test(data_chisq)

### Plot Bland-Altman Comparison ### ----
heatscatter(data$ICC_NoAdj,data$ICC_Adj,
            colpal="bl2gr2rd",
            main="",
            xlab = "No Immune Cell Adjustment",
            ylab = "Immune Cell Adjustment",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

###################
### Scenario 10 ###
###################

### Read in Data ### ----
load("scen10_icc_s0s1_t1_ela01_revision.RData")
icc.comb.scen10 <- dplyr::filter(icc.comb.scen10, icc.comb.scen10$type == "ICC(2,1)")[,c(1,2,5)]
icc.comb.scen10_NoAdj <- icc.comb.scen10[1:865859,]
icc.comb.scen10_Adj <- icc.comb.scen10[865860:1731718,]
colnames(icc.comb.scen10_NoAdj)[3] <- "ICC_NoAdj"
colnames(icc.comb.scen10_Adj)[3] <- "ICC_Adj"
data <- merge(icc.comb.scen10_NoAdj, icc.comb.scen10_Adj, by = "DNAm.probe")
data <- data[,c(1,3,5)]
Exclude_Probes_Comb <- readRDS("Exclude_Probes_Comb.Rds")
keep <- !(data$DNAm.probe %in% Exclude_Probes_Comb$DNAm.probe)
data <- data[keep,]

### Perform Statistical Testing ### ----
# Paired T-test #
t.test(data$ICC_NoAdj, data$ICC_Adj, paired = TRUE)
# Proportional bias test #
data$diff <- data$ICC_NoAdj - data$ICC_Adj
data$mean <- (data$ICC_NoAdj + data$ICC_Adj)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$ICC_NoAdj, data$ICC_Adj)
# Chi-square test #
length(data$ICC_NoAdj[data$ICC_NoAdj <= 0.01])
length(data$ICC_Adj[data$ICC_Adj <= 0.01])
data_chisq <- data.frame(ICC_NoAdj = c(63725, length(data$ICC_NoAdj) - 63725),
                         ICC_Adj = c(59371, length(data$ICC_Adj) - 59371))
chisq.test(data_chisq)

### Plot Bland-Altman Comparison ### ----
heatscatter(data$ICC_NoAdj,data$ICC_Adj,
            colpal="bl2gr2rd",
            main="",
            xlab = "No Immune Cell Adjustment",
            ylab = "Immune Cell Adjustment",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

############################################################################
####################### Sample Size Calculations ###########################
############################################################################

### Read in Data ### ----
load("scen6_icc_s1_t1t2_ela01_revision.RData")
icc.comb.scen6 <- dplyr::filter(icc.comb.scen6, icc.comb.scen6$covariates == "Adjusted for DNAm_Monocytes, Methylation_Batch")
icc.comb.scen6 <- dplyr::filter(icc.comb.scen6, icc.comb.scen6$type == "ICC(2,1)")[,c(1,5)]
colnames(icc.comb.scen6)[2] <- "ICC_Scen6"
load("scen7_icc_s1_t1t3_ela01_revision.RData")
icc.comb.scen7 <- dplyr::filter(icc.comb.scen7, icc.comb.scen7$covariates == "Adjusted for DNAm_Monocytes, Methylation_Batch")
icc.comb.scen7 <- dplyr::filter(icc.comb.scen7, icc.comb.scen7$type == "ICC(2,1)")[,c(1,5)]
colnames(icc.comb.scen7)[2] <- "ICC_Scen7"
load("scen8_icc_s1_t1t4_ela01_revision.RData")
icc.comb.scen8 <- dplyr::filter(icc.comb.scen8, icc.comb.scen8$covariates == "Adjusted for DNAm_Monocytes, Methylation_Batch")
icc.comb.scen8 <- dplyr::filter(icc.comb.scen8, icc.comb.scen8$type == "ICC(2,1)")[,c(1,5)]
colnames(icc.comb.scen8)[2] <- "ICC_Scen8"
load("scen9_icc_s1_t3t4_ela01_revision.RData")
icc.comb.scen9 <- dplyr::filter(icc.comb.scen9, icc.comb.scen9$covariates == "Adjusted for DNAm_Monocytes, Methylation_Batch")
icc.comb.scen9 <- dplyr::filter(icc.comb.scen9, icc.comb.scen9$type == "ICC(2,1)")[,c(1,5)]
colnames(icc.comb.scen9)[2] <- "ICC_Scen9"
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

### Merge all Data ### ----
data <- merge(icc.comb.scen6, icc.comb.scen7, by = "DNAm.probe")
data <- merge(data, icc.comb.scen8, by = "DNAm.probe")
data <- merge(data, icc.comb.scen9, by = "DNAm.probe")
data <- merge(data, icc.comb.scen6sampled, by = "DNAm.probe")
data <- merge(data, icc.comb.scen7sampled, by = "DNAm.probe")
data <- merge(data, icc.comb.scen8sampled, by = "DNAm.probe")
data <- merge(data, icc.comb.scen9sampled, by = "DNAm.probe")
Exclude_Probes_Comb <- readRDS("Exclude_Probes_Comb.Rds")
keep <- !(data$DNAm.probe %in% Exclude_Probes_Comb$DNAm.probe)
data <- data[keep,]

# Delete old dataframes #
icc.comb.scen6 <- NULL
icc.comb.scen7 <- NULL
icc.comb.scen8 <- NULL
icc.comb.scen9 <- NULL
icc.comb.scen6sampled <- NULL
icc.comb.scen7sampled <- NULL
icc.comb.scen8sampled <- NULL
icc.comb.scen9sampled <- NULL

### Perform Statistical Testing ### ----
## Scenario 6 ##
# Paired T-test #
t.test(data$ICC_Scen6, data$ICC_Scen6samp, paired = TRUE)
# Proportional bias test #
data$diff <- data$ICC_Scen6 - data$ICC_Scen6samp
data$mean <- (data$ICC_Scen6 + data$ICC_Scen6samp)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$ICC_Scen6, data$ICC_Scen6samp)
# Chi-square test #
length(data$ICC_Scen6[data$ICC_Scen6 <= 0.01])
length(data$ICC_Scen6samp[data$ICC_Scen6samp <= 0.01])
data_chisq <- data.frame(ICC_Scen6 = c(length(data$ICC_Scen6[data$ICC_Scen6 <= 0.01]), 
                                    length(data$ICC_Scen6) - length(data$ICC_Scen6[data$ICC_Scen6 <= 0.01])),
                         ICC_Scen6samp = c(length(data$ICC_Scen6samp[data$ICC_Scen6samp <= 0.01]), 
                                    length(data$ICC_Scen6samp) - length(data$ICC_Scen6samp[data$ICC_Scen6samp <= 0.01])))
chisq.test(data_chisq)

# Plot Bland-Altman Comparison # 
heatscatter(data$ICC_Scen6samp,data$ICC_Scen6,
            colpal="bl2gr2rd",
            main="StressT1-2",
            xlab = "N = 14",
            ylab = "N = 31",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

## Scenario 7 ##
# Paired T-test #
t.test(data$ICC_Scen7, data$ICC_Scen7samp, paired = TRUE)
# Proportional bias test #
data$diff <- data$ICC_Scen7 - data$ICC_Scen7samp
data$mean <- (data$ICC_Scen7 + data$ICC_Scen7samp)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$ICC_Scen7, data$ICC_Scen7samp)
# Chi-square test #
data_chisq <- data.frame(ICC_Scen7 = c(length(data$ICC_Scen7[data$ICC_Scen7 <= 0.01]), 
                                       length(data$ICC_Scen7) - length(data$ICC_Scen7[data$ICC_Scen7 <= 0.01])),
                         ICC_Scen7samp = c(length(data$ICC_Scen7samp[data$ICC_Scen7samp <= 0.01]), 
                                           length(data$ICC_Scen7samp) - length(data$ICC_Scen7samp[data$ICC_Scen7samp <= 0.01])))
chisq.test(data_chisq)

# Plot Bland-Altman Comparison # 
heatscatter(data$ICC_Scen7samp,data$ICC_Scen7,
            colpal="bl2gr2rd",
            main="StressT1-3",
            xlab = "N = 14",
            ylab = "N = 30",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

## Scenario 8 ##
# Paired T-test #
t.test(data$ICC_Scen8, data$ICC_Scen8samp, paired = TRUE)
# Proportional bias test #
data$diff <- data$ICC_Scen8 - data$ICC_Scen8samp
data$mean <- (data$ICC_Scen8 + data$ICC_Scen8samp)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$ICC_Scen8, data$ICC_Scen8samp)
# Chi-square test #
data_chisq <- data.frame(ICC_Scen8 = c(length(data$ICC_Scen8[data$ICC_Scen8 <= 0.01]), 
                                       length(data$ICC_Scen8) - length(data$ICC_Scen8[data$ICC_Scen8 <= 0.01])),
                         ICC_Scen8samp = c(length(data$ICC_Scen8samp[data$ICC_Scen8samp <= 0.01]), 
                                           length(data$ICC_Scen8samp) - length(data$ICC_Scen8samp[data$ICC_Scen8samp <= 0.01])))
chisq.test(data_chisq)

# Plot Bland-Altman Comparison # 
heatscatter(data$ICC_Scen8samp,data$ICC_Scen8,
            colpal="bl2gr2rd",
            main="StressT1-4",
            xlab = "N = 14",
            ylab = "N = 31",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

## Scenario 9 ##
# Paired T-test #
t.test(data$ICC_Scen9, data$ICC_Scen9samp, paired = TRUE)
# Proportional bias test #
data$diff <- data$ICC_Scen9 - data$ICC_Scen9samp
data$mean <- (data$ICC_Scen9 + data$ICC_Scen9samp)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$ICC_Scen9, data$ICC_Scen9samp)
# Chi-square test #
data_chisq <- data.frame(ICC_Scen9 = c(length(data$ICC_Scen9[data$ICC_Scen9 <= 0.01]), 
                                       length(data$ICC_Scen9) - length(data$ICC_Scen9[data$ICC_Scen9 <= 0.01])),
                         ICC_Scen9samp = c(length(data$ICC_Scen9samp[data$ICC_Scen9samp <= 0.01]), 
                                           length(data$ICC_Scen9samp) - length(data$ICC_Scen9samp[data$ICC_Scen9samp <= 0.01])))
chisq.test(data_chisq)

# Plot Bland-Altman Comparison # 
heatscatter(data$ICC_Scen9samp,data$ICC_Scen9,
            colpal="bl2gr2rd",
            main="StressT3-4",
            xlab = "N = 14",
            ylab = "N = 29",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)


############################################################################
##################### Repeated Measures Calculations #######################
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
load("scen12sampled_icc_s1_t1t2t3t4_ela01_revision.Rdata")
icc.comb.scen12sampled <- dplyr::filter(icc.comb.scen12sampled, icc.comb.scen12sampled$covariates == "Adjusted for DNAm_Monocytes, Methylation_Batch")
icc.comb.scen12sampled <- dplyr::filter(icc.comb.scen12sampled, icc.comb.scen12sampled$type == "ICC(2,k)")[,c(1,5)]
colnames(icc.comb.scen12sampled)[2] <- "ICC_Scen12samp"

### Merge all Data ### ----
data <- merge(icc.comb.scen1, icc.comb.scen2, by = "DNAm.probe")
data <- merge(data, icc.comb.scen3, by = "DNAm.probe")
data <- merge(data, icc.comb.scen4, by = "DNAm.probe")
data <- merge(data, icc.comb.scen5, by = "DNAm.probe")
data <- merge(data, icc.comb.scen6sampled, by = "DNAm.probe")
data <- merge(data, icc.comb.scen7sampled, by = "DNAm.probe")
data <- merge(data, icc.comb.scen8sampled, by = "DNAm.probe")
data <- merge(data, icc.comb.scen9sampled, by = "DNAm.probe")
data <- merge(data, icc.comb.scen12sampled, by = "DNAm.probe")
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
icc.comb.scen12sampled <- NULL

### Perform Statistical Testing ### ----
## Scenario 5 and 1 ##
# Paired T-test #
t.test(data$ICC_Scen5, data$ICC_Scen1, paired = TRUE)
# Proportional bias test #
data$diff <- data$ICC_Scen5 - data$ICC_Scen1
data$mean <- (data$ICC_Scen5 + data$ICC_Scen1)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$ICC_Scen5, data$ICC_Scen1)
# Chi-square test #
length(data$ICC_Scen5[data$ICC_Scen5 <= 0.01])
length(data$ICC_Scen1[data$ICC_Scen1 <= 0.01])
data_chisq <- data.frame(ICC_Scen5 = c(length(data$ICC_Scen5[data$ICC_Scen5 <= 0.01]), 
                                       length(data$ICC_Scen5) - length(data$ICC_Scen5[data$ICC_Scen5 <= 0.01])),
                         ICC_Scen1 = c(length(data$ICC_Scen1[data$ICC_Scen1 <= 0.01]), 
                                           length(data$ICC_Scen1) - length(data$ICC_Scen1[data$ICC_Scen1 <= 0.01])))
chisq.test(data_chisq)

# Plot Bland-Altman Comparison # 
heatscatter(data$ICC_Scen5,data$ICC_Scen1,
            colpal="bl2gr2rd",
            main="NoStressT1-2-3-4 vs NoStressT1-2",
            xlab = "4 Repeated Measures",
            ylab = "2 Repeated Measures",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

## Scenario 5 and 2 ##
# Paired T-test #
t.test(data$ICC_Scen5, data$ICC_Scen2, paired = TRUE)
# Proportional bias test #
data$diff <- data$ICC_Scen5 - data$ICC_Scen2
data$mean <- (data$ICC_Scen5 + data$ICC_Scen2)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$ICC_Scen5, data$ICC_Scen2)
# Chi-square test #
data_chisq <- data.frame(ICC_Scen5 = c(length(data$ICC_Scen5[data$ICC_Scen5 <= 0.01]), 
                                       length(data$ICC_Scen5) - length(data$ICC_Scen5[data$ICC_Scen5 <= 0.01])),
                         ICC_Scen2 = c(length(data$ICC_Scen2[data$ICC_Scen2 <= 0.01]), 
                                           length(data$ICC_Scen2) - length(data$ICC_Scen2[data$ICC_Scen2 <= 0.01])))
chisq.test(data_chisq)

# Plot Bland-Altman Comparison # 
heatscatter(data$ICC_Scen5,data$ICC_Scen2,
            colpal="bl2gr2rd",
            main="NoStressT1-2-3-4 vs NoStressT1-3",
            xlab = "4 Repeated Measures",
            ylab = "2 Repeated Measures",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

## Scenario 5 and 3 ##
# Paired T-test #
t.test(data$ICC_Scen5, data$ICC_Scen3, paired = TRUE)
# Proportional bias test #
data$diff <- data$ICC_Scen5 - data$ICC_Scen3
data$mean <- (data$ICC_Scen5 + data$ICC_Scen3)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$ICC_Scen5, data$ICC_Scen3)
# Chi-square test #
data_chisq <- data.frame(ICC_Scen5 = c(length(data$ICC_Scen5[data$ICC_Scen5 <= 0.01]), 
                                       length(data$ICC_Scen5) - length(data$ICC_Scen5[data$ICC_Scen5 <= 0.01])),
                         ICC_Scen3 = c(length(data$ICC_Scen3[data$ICC_Scen3 <= 0.01]), 
                                           length(data$ICC_Scen3) - length(data$ICC_Scen3[data$ICC_Scen3 <= 0.01])))
chisq.test(data_chisq)

# Plot Bland-Altman Comparison # 
heatscatter(data$ICC_Scen5,data$ICC_Scen3,
            colpal="bl2gr2rd",
            main="NoStressT1-2-3-4 vs NoStressT1-4",
            xlab = "4 Repeated Measures",
            ylab = "2 Repeated Measures",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)

## Scenario 5 and 4 ##
# Paired T-test #
t.test(data$ICC_Scen5, data$ICC_Scen4, paired = TRUE)
# Proportional bias test #
data$diff <- data$ICC_Scen5 - data$ICC_Scen4
data$mean <- (data$ICC_Scen5 + data$ICC_Scen4)/2
summary(lm(diff ~ 1 + mean, data))
# KS test #
ks.test(data$ICC_Scen5, data$ICC_Scen4)
# Chi-square test #
data_chisq <- data.frame(ICC_Scen5 = c(length(data$ICC_Scen5[data$ICC_Scen5 <= 0.01]), 
                                       length(data$ICC_Scen5) - length(data$ICC_Scen5[data$ICC_Scen5 <= 0.01])),
                         ICC_Scen4 = c(length(data$ICC_Scen4[data$ICC_Scen4 <= 0.01]), 
                                           length(data$ICC_Scen4) - length(data$ICC_Scen4[data$ICC_Scen4 <= 0.01])))
chisq.test(data_chisq)

# Plot Bland-Altman Comparison # 
heatscatter(data$ICC_Scen5,data$ICC_Scen4,
            colpal="bl2gr2rd",
            main="NoStressT1-2-3-4 vs NoStressT3-4",
            xlab = "4 Repeated Measures",
            ylab = "2 Repeated Measures",
            cor=FALSE, 
            add.contour = FALSE, 
            alpha = 50)
