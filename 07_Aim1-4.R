### Testing Aims 1-4 ### ----

### Load Packages ### ----
library(dplyr)
library(psych)
library(ggplot2)
library(ggpubr)
library(dgof)
library(LSD)

############################################################################
################################### Aim 1 ##################################
############################################################################

##################
### Scenario 1 ###
##################

### Read in Data ### ----
load("scen1_icc_s0_t1t2.RData")
scen1.icc.comb <- dplyr::filter(scen1.icc.comb, scen1.icc.comb$covariates == "No covariates")[,c(1,4:5)]
scen1.icc.comb_21 <- scen1.icc.comb[1:865859,]
scen1.icc.comb_2k <- scen1.icc.comb[865860:1731718,]
colnames(scen1.icc.comb_21)[3] <- "ICC(2,1)"
colnames(scen1.icc.comb_2k)[3] <- "ICC(2,k)"
data <- merge(scen1.icc.comb_21, scen1.icc.comb_2k, by = "DNAm.probe")
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
load("scen2_icc_s0_t1t3.RData")
scen2.icc.comb <- dplyr::filter(scen3.icc.comb, scen3.icc.comb$covariates == "No covariates")[,c(1,4:5)]
scen2.icc.comb_21 <- scen2.icc.comb[1:865859,]
scen2.icc.comb_2k <- scen2.icc.comb[865860:1731718,]
colnames(scen2.icc.comb_21)[3] <- "ICC(2,1)"
colnames(scen2.icc.comb_2k)[3] <- "ICC(2,k)"
data <- merge(scen2.icc.comb_21, scen2.icc.comb_2k, by = "DNAm.probe")
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
load("scen3_icc_s0_t1t4.RData")
scen3.icc.comb <- dplyr::filter(scen11.icc.comb, scen11.icc.comb$covariates == "No covariates")[,c(1,4:5)]
scen3.icc.comb_21 <- scen3.icc.comb[1:865859,]
scen3.icc.comb_2k <- scen3.icc.comb[865860:1731718,]
colnames(scen3.icc.comb_21)[3] <- "ICC(2,1)"
colnames(scen3.icc.comb_2k)[3] <- "ICC(2,k)"
data <- merge(scen3.icc.comb_21, scen3.icc.comb_2k, by = "DNAm.probe")
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
load("scen4_icc_s0_t3t4.RData")
scen4.icc.comb <- dplyr::filter(scen4.icc.comb, scen4.icc.comb$covariates == "No covariates")[,c(1,4:5)]
scen4.icc.comb_21 <- scen4.icc.comb[1:865859,]
scen4.icc.comb_2k <- scen4.icc.comb[865860:1731718,]
colnames(scen4.icc.comb_21)[3] <- "ICC(2,1)"
colnames(scen4.icc.comb_2k)[3] <- "ICC(2,k)"
data <- merge(scen4.icc.comb_21, scen4.icc.comb_2k, by = "DNAm.probe")
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
load("scen5_icc_s0_t1t2t3t4.RData")
scen5.icc.comb <- dplyr::filter(scen5.icc.comb, scen5.icc.comb$covariates == "No covariates")[,c(1,4:5)]
scen5.icc.comb_21 <- scen5.icc.comb[1:865859,]
scen5.icc.comb_2k <- scen5.icc.comb[865860:1731718,]
colnames(scen5.icc.comb_21)[3] <- "ICC(2,1)"
colnames(scen5.icc.comb_2k)[3] <- "ICC(2,k)"
data <- merge(scen5.icc.comb_21, scen5.icc.comb_2k, by = "DNAm.probe")
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
load("scen6_icc_s1_t1t2.RData")
scen6.icc.comb <- dplyr::filter(scen6.icc.comb, scen6.icc.comb$covariates == "No covariates")[,c(1,4:5)]
scen6.icc.comb_21 <- scen6.icc.comb[1:865859,]
scen6.icc.comb_2k <- scen6.icc.comb[865860:1731718,]
colnames(scen6.icc.comb_21)[3] <- "ICC(2,1)"
colnames(scen6.icc.comb_2k)[3] <- "ICC(2,k)"
data <- merge(scen6.icc.comb_21, scen6.icc.comb_2k, by = "DNAm.probe")
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
load("scen7_icc_s1_t1t3.RData")
scen7.icc.comb <- dplyr::filter(scen7.icc.comb, scen7.icc.comb$covariates == "No covariates")[,c(1,4:5)]
scen7.icc.comb_21 <- scen7.icc.comb[1:865859,]
scen7.icc.comb_2k <- scen7.icc.comb[865860:1731718,]
colnames(scen7.icc.comb_21)[3] <- "ICC(2,1)"
colnames(scen7.icc.comb_2k)[3] <- "ICC(2,k)"
data <- merge(scen7.icc.comb_21, scen7.icc.comb_2k, by = "DNAm.probe")
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
load("scen8_icc_s1_t1t4.RData")
scen8.icc.comb <- dplyr::filter(scen8.icc.comb, scen8.icc.comb$covariates == "No covariates")[,c(1,4:5)]
scen8.icc.comb_21 <- scen8.icc.comb[1:865859,]
scen8.icc.comb_2k <- scen8.icc.comb[865860:1731718,]
colnames(scen8.icc.comb_21)[3] <- "ICC(2,1)"
colnames(scen8.icc.comb_2k)[3] <- "ICC(2,k)"
data <- merge(scen8.icc.comb_21, scen8.icc.comb_2k, by = "DNAm.probe")
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
load("scen9_icc_s1_t3t4.RData")
scen9.icc.comb <- dplyr::filter(scen9.icc.comb, scen9.icc.comb$covariates == "No covariates")[,c(1,4:5)]
scen9.icc.comb_21 <- scen9.icc.comb[1:865859,]
scen9.icc.comb_2k <- scen9.icc.comb[865860:1731718,]
colnames(scen9.icc.comb_21)[3] <- "ICC(2,1)"
colnames(scen9.icc.comb_2k)[3] <- "ICC(2,k)"
data <- merge(scen9.icc.comb_21, scen9.icc.comb_2k, by = "DNAm.probe")
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
load("scen10_icc_s1s2_t1.RData")
scen10.icc.comb <- dplyr::filter(scen10.icc.comb, scen10.icc.comb$covariates == "No covariates")[,c(1,4:5)]
scen10.icc.comb_21 <- scen10.icc.comb[1:865859,]
scen10.icc.comb_2k <- scen10.icc.comb[865860:1731718,]
colnames(scen10.icc.comb_21)[3] <- "ICC(2,1)"
colnames(scen10.icc.comb_2k)[3] <- "ICC(2,k)"
data <- merge(scen10.icc.comb_21, scen10.icc.comb_2k, by = "DNAm.probe")
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
################################### Aim 2 ##################################
############################################################################

##################
### Scenario 1 ###
##################

### Read in Data ### ----
load("scen1_icc_s0_t1t2.RData")
scen1.icc.comb <- dplyr::filter(scen1.icc.comb, scen1.icc.comb$type == "ICC(2,1)")[,c(1,2,5)]
scen1.icc.comb_NoAdj <- scen1.icc.comb[1:865859,]
scen1.icc.comb_Adj <- scen1.icc.comb[865860:1731718,]
colnames(scen1.icc.comb_NoAdj)[3] <- "ICC_NoAdj"
colnames(scen1.icc.comb_Adj)[3] <- "ICC_Adj"
data <- merge(scen1.icc.comb_NoAdj, scen1.icc.comb_Adj, by = "DNAm.probe")
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
load("scen2_icc_s0_t1t3.RData")
scen2.icc.comb <- dplyr::filter(scen3.icc.comb, scen3.icc.comb$type == "ICC(2,1)")[,c(1,2,5)]
scen2.icc.comb_NoAdj <- scen2.icc.comb[1:865859,]
scen2.icc.comb_Adj <- scen2.icc.comb[865860:1731718,]
colnames(scen2.icc.comb_NoAdj)[3] <- "ICC_NoAdj"
colnames(scen2.icc.comb_Adj)[3] <- "ICC_Adj"
data <- merge(scen2.icc.comb_NoAdj, scen2.icc.comb_Adj, by = "DNAm.probe")
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
load("scen3_icc_s0_t1t4.RData")
scen3.icc.comb <- dplyr::filter(scen11.icc.comb, scen11.icc.comb$type == "ICC(2,1)")[,c(1,2,5)]
scen3.icc.comb_NoAdj <- scen3.icc.comb[1:865859,]
scen3.icc.comb_Adj <- scen3.icc.comb[865860:1731718,]
colnames(scen3.icc.comb_NoAdj)[3] <- "ICC_NoAdj"
colnames(scen3.icc.comb_Adj)[3] <- "ICC_Adj"
data <- merge(scen3.icc.comb_NoAdj, scen3.icc.comb_Adj, by = "DNAm.probe")
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
load("scen4_icc_s0_t3t4.RData")
scen4.icc.comb <- dplyr::filter(scen4.icc.comb, scen4.icc.comb$type == "ICC(2,1)")[,c(1,2,5)]
scen4.icc.comb_NoAdj <- scen4.icc.comb[1:865859,]
scen4.icc.comb_Adj <- scen4.icc.comb[865860:1731718,]
colnames(scen4.icc.comb_NoAdj)[3] <- "ICC_NoAdj"
colnames(scen4.icc.comb_Adj)[3] <- "ICC_Adj"
data <- merge(scen4.icc.comb_NoAdj, scen4.icc.comb_Adj, by = "DNAm.probe")
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
load("scen5_icc_s0_t1t2t3t4.RData")
scen5.icc.comb <- dplyr::filter(scen5.icc.comb, scen5.icc.comb$type == "ICC(2,1)")[,c(1,2,5)]
scen5.icc.comb_NoAdj <- scen5.icc.comb[1:865859,]
scen5.icc.comb_Adj <- scen5.icc.comb[865860:1731718,]
colnames(scen5.icc.comb_NoAdj)[3] <- "ICC_NoAdj"
colnames(scen5.icc.comb_Adj)[3] <- "ICC_Adj"
data <- merge(scen5.icc.comb_NoAdj, scen5.icc.comb_Adj, by = "DNAm.probe")
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
load("scen6sampled_icc_s1_t1t2.RData")
scen6.icc.comb <- dplyr::filter(scen6sampled.icc.comb, scen6sampled.icc.comb$type == "ICC(2,1)")[,c(1,2,5)]
scen6.icc.comb_NoAdj <- scen6.icc.comb[1:865859,]
scen6.icc.comb_Adj <- scen6.icc.comb[865860:1731718,]
colnames(scen6.icc.comb_NoAdj)[3] <- "ICC_NoAdj"
colnames(scen6.icc.comb_Adj)[3] <- "ICC_Adj"
data <- merge(scen6.icc.comb_NoAdj, scen6.icc.comb_Adj, by = "DNAm.probe")
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
load("scen7sampled_icc_s1_t1t3.RData")
scen7.icc.comb <- dplyr::filter(scen7sampled.icc.comb, scen7sampled.icc.comb$type == "ICC(2,1)")[,c(1,2,5)]
scen7.icc.comb_NoAdj <- scen7.icc.comb[1:865859,]
scen7.icc.comb_Adj <- scen7.icc.comb[865860:1731718,]
colnames(scen7.icc.comb_NoAdj)[3] <- "ICC_NoAdj"
colnames(scen7.icc.comb_Adj)[3] <- "ICC_Adj"
data <- merge(scen7.icc.comb_NoAdj, scen7.icc.comb_Adj, by = "DNAm.probe")
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
load("scen8sampled_icc_s1_t1t4.RData")
scen8.icc.comb <- dplyr::filter(scen8sampled.icc.comb, scen8sampled.icc.comb$type == "ICC(2,1)")[,c(1,2,5)]
scen8.icc.comb_NoAdj <- scen8.icc.comb[1:865859,]
scen8.icc.comb_Adj <- scen8.icc.comb[865860:1731718,]
colnames(scen8.icc.comb_NoAdj)[3] <- "ICC_NoAdj"
colnames(scen8.icc.comb_Adj)[3] <- "ICC_Adj"
data <- merge(scen8.icc.comb_NoAdj, scen8.icc.comb_Adj, by = "DNAm.probe")
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
load("scen9sampled_icc_s1_t3t4.RData")
scen9.icc.comb <- dplyr::filter(scen9sampled.icc.comb, scen9sampled.icc.comb$type == "ICC(2,1)")[,c(1,2,5)]
scen9.icc.comb_NoAdj <- scen9.icc.comb[1:865859,]
scen9.icc.comb_Adj <- scen9.icc.comb[865860:1731718,]
colnames(scen9.icc.comb_NoAdj)[3] <- "ICC_NoAdj"
colnames(scen9.icc.comb_Adj)[3] <- "ICC_Adj"
data <- merge(scen9.icc.comb_NoAdj, scen9.icc.comb_Adj, by = "DNAm.probe")
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
load("scen10_icc_s1s2_t1.RData")
scen10.icc.comb <- dplyr::filter(scen10.icc.comb, scen10.icc.comb$type == "ICC(2,1)")[,c(1,2,5)]
scen10.icc.comb_NoAdj <- scen10.icc.comb[1:865859,]
scen10.icc.comb_Adj <- scen10.icc.comb[865860:1731718,]
colnames(scen10.icc.comb_NoAdj)[3] <- "ICC_NoAdj"
colnames(scen10.icc.comb_Adj)[3] <- "ICC_Adj"
data <- merge(scen10.icc.comb_NoAdj, scen10.icc.comb_Adj, by = "DNAm.probe")
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
################################### Aim 3 ##################################
############################################################################
### Read in Data ### ----
load("scen6_icc_s1_t1t2.RData")
scen6.icc.comb <- dplyr::filter(scen6.icc.comb, scen6.icc.comb$covariates == "Adjusted for DNAm_Monocytes")
scen6.icc.comb <- dplyr::filter(scen6.icc.comb, scen6.icc.comb$type == "ICC(2,1)")[,c(1,5)]
colnames(scen6.icc.comb)[2] <- "ICC_Scen6"
load("scen7_icc_s1_t1t3.RData")
scen7.icc.comb <- dplyr::filter(scen7.icc.comb, scen7.icc.comb$covariates == "Adjusted for DNAm_Monocytes")
scen7.icc.comb <- dplyr::filter(scen7.icc.comb, scen7.icc.comb$type == "ICC(2,1)")[,c(1,5)]
colnames(scen7.icc.comb)[2] <- "ICC_Scen7"
load("scen8_icc_s1_t1t4.RData")
scen8.icc.comb <- dplyr::filter(scen8.icc.comb, scen8.icc.comb$covariates == "Adjusted for DNAm_Monocytes")
scen8.icc.comb <- dplyr::filter(scen8.icc.comb, scen8.icc.comb$type == "ICC(2,1)")[,c(1,5)]
colnames(scen8.icc.comb)[2] <- "ICC_Scen8"
load("scen9_icc_s1_t3t4.RData")
scen9.icc.comb <- dplyr::filter(scen9.icc.comb, scen9.icc.comb$covariates == "Adjusted for DNAm_Monocytes")
scen9.icc.comb <- dplyr::filter(scen9.icc.comb, scen9.icc.comb$type == "ICC(2,1)")[,c(1,5)]
colnames(scen9.icc.comb)[2] <- "ICC_Scen9"
load("scen6sampled_icc_s1_t1t2.RData")
scen6sampled.icc.comb <- dplyr::filter(scen6sampled.icc.comb, scen6sampled.icc.comb$covariates == "Adjusted for DNAm_Monocytes")
scen6sampled.icc.comb <- dplyr::filter(scen6sampled.icc.comb, scen6sampled.icc.comb$type == "ICC(2,k)")[,c(1,5)]
colnames(scen6sampled.icc.comb)[2] <- "ICC_Scen6samp"
load("scen7sampled_icc_s1_t1t3.RData")
scen7sampled.icc.comb <- dplyr::filter(scen7sampled.icc.comb, scen7sampled.icc.comb$covariates == "Adjusted for DNAm_Monocytes")
scen7sampled.icc.comb <- dplyr::filter(scen7sampled.icc.comb, scen7sampled.icc.comb$type == "ICC(2,k)")[,c(1,5)]
colnames(scen7sampled.icc.comb)[2] <- "ICC_Scen7samp"
load("scen8sampled_icc_s1_t1t4.RData")
scen8sampled.icc.comb <- dplyr::filter(scen8sampled.icc.comb, scen8sampled.icc.comb$covariates == "Adjusted for DNAm_Monocytes")
scen8sampled.icc.comb <- dplyr::filter(scen8sampled.icc.comb, scen8sampled.icc.comb$type == "ICC(2,k)")[,c(1,5)]
colnames(scen8sampled.icc.comb)[2] <- "ICC_Scen8samp"
load("scen9sampled_icc_s1_t3t4.RData")
scen9sampled.icc.comb <- dplyr::filter(scen9sampled.icc.comb, scen9sampled.icc.comb$covariates == "Adjusted for DNAm_Monocytes")
scen9sampled.icc.comb <- dplyr::filter(scen9sampled.icc.comb, scen9sampled.icc.comb$type == "ICC(2,k)")[,c(1,5)]
colnames(scen9sampled.icc.comb)[2] <- "ICC_Scen9samp"

### Merge all Data ### ----
data <- merge(scen6.icc.comb, scen7.icc.comb, by = "DNAm.probe")
data <- merge(data, scen8.icc.comb, by = "DNAm.probe")
data <- merge(data, scen9.icc.comb, by = "DNAm.probe")
data <- merge(data, scen6sampled.icc.comb, by = "DNAm.probe")
data <- merge(data, scen7sampled.icc.comb, by = "DNAm.probe")
data <- merge(data, scen8sampled.icc.comb, by = "DNAm.probe")
data <- merge(data, scen9sampled.icc.comb, by = "DNAm.probe")
Exclude_Probes_Comb <- readRDS("Exclude_Probes_Comb.Rds")
keep <- !(data$DNAm.probe %in% Exclude_Probes_Comb$DNAm.probe)
data <- data[keep,]

# Delete old dataframes #
scen6.icc.comb <- NULL
scen7.icc.comb <- NULL
scen8.icc.comb <- NULL
scen9.icc.comb <- NULL
scen6sampled.icc.comb <- NULL
scen7sampled.icc.comb <- NULL
scen8sampled.icc.comb <- NULL
scen9sampled.icc.comb <- NULL

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

### Create Histograms for Comparison ### ----
describe(data)

hist6 <- ggplot(data, aes(x=ICC_Scen6)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.43),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,80000) +
  labs(x = NULL, y = NULL, title = "Scenario 6") +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        plot.title = element_text(hjust = 0.5))  

hist7 <- ggplot(data, aes(x=ICC_Scen7)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.45),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,80000) +
  labs(x = NULL, y = NULL, title = "Scenario 7") +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        plot.title = element_text(hjust = 0.5))  

hist8 <- ggplot(data, aes(x=ICC_Scen8)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.47),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,80000) +
  labs(x = NULL, y = NULL, title = "Scenario 8") +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        plot.title = element_text(hjust = 0.5))  

hist9 <- ggplot(data, aes(x=ICC_Scen9)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.43),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,80000) +
  labs(x = NULL, y = NULL, title = "Scenario 9") +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        plot.title = element_text(hjust = 0.5))  

hist6samp <- ggplot(data, aes(x=ICC_Scen6samp)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.64),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,80000) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        plot.title = element_text(hjust = 0.5))  

hist7samp <- ggplot(data, aes(x=ICC_Scen7samp)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.64),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,80000) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        plot.title = element_text(hjust = 0.5))  

hist8samp <- ggplot(data, aes(x=ICC_Scen8samp)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.65),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,80000) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        plot.title = element_text(hjust = 0.5))  

hist9samp <- ggplot(data, aes(x=ICC_Scen9samp)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.62),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,80000) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), #remove y axis ticks
        plot.title = element_text(hjust = 0.5))  

ggarrange(hist6, hist7, hist8, hist9,
          hist6samp, hist7samp, hist8samp, hist9samp,
          ncol = 4, nrow = 2)

############################################################################
################################### Aim 4 ##################################
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

### Merge all Data ### ----
data <- merge(scen1.icc.comb, scen2.icc.comb, by = "DNAm.probe")
data <- merge(data, scen3.icc.comb, by = "DNAm.probe")
data <- merge(data, scen4.icc.comb, by = "DNAm.probe")
data <- merge(data, scen5.icc.comb, by = "DNAm.probe")
Exclude_Probes_Comb <- readRDS("Exclude_Probes_Comb.Rds")
keep <- !(data$DNAm.probe %in% Exclude_Probes_Comb$DNAm.probe)
data <- data[keep,]

# Delete old dataframes #
scen1.icc.comb <- NULL
scen2.icc.comb <- NULL
scen3.icc.comb <- NULL
scen4.icc.comb <- NULL
scen5.icc.comb <- NULL
scen11.icc.comb <- NULL

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
