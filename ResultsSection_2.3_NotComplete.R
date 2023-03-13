### Results - Section 2.3 ### ----

### Load Packages ### ----
library(dplyr)
library(psych)
library(ggplot2)
library(ggpubr)
library(dgof)
library(fgsea)
library(tidyverse)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)

############################################################################
######################## Epigenetic Aging Clocks ###########################
############################################################################
### Set Seed ### ----
set.seed(100)

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
# Exclude probes with ICC < 0.01 #
data[,2:11][data[,2:11] < 0.01] <- NA

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

######################
### Horvath Probes ###
######################
### Probe Set Enrichment Analysis ### ----
## Scenario 1 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("Horvath_Probes.csv")
keep <- mQTL_probes$DNAm.probe %in% data$DNAm.probe
mQTL_probes <- mQTL_probes[keep,]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S1 <- data[,c(1,2)]
ICC_ranks_S1 <- ICC_ranks_S1 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S1 <- as.numeric(ICC_ranks_S1$ICC_Scen1)
names(ICC_ranks_S1) <- data$DNAm.probe
ICC_ranks_S1 <- sort(ICC_ranks_S1, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S1 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S1,
                           nperm = 10000)

## Scenario 2 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("Horvath_Probes.csv")
keep <- mQTL_probes$DNAm.probe %in% data$DNAm.probe
mQTL_probes <- mQTL_probes[keep,]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S2 <- data[,c(1,3)]
ICC_ranks_S2 <- ICC_ranks_S2 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S2 <- as.numeric(ICC_ranks_S2$ICC_Scen2)
names(ICC_ranks_S2) <- data$DNAm.probe
ICC_ranks_S2 <- sort(ICC_ranks_S2, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S2 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S2,
                           nperm = 10000)

## Scenario 3 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("Horvath_Probes.csv")
keep <- mQTL_probes$DNAm.probe %in% data$DNAm.probe
mQTL_probes <- mQTL_probes[keep,]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S3 <- data[,c(1,4)]
ICC_ranks_S3 <- ICC_ranks_S3 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S3 <- as.numeric(ICC_ranks_S3$ICC_Scen3)
names(ICC_ranks_S3) <- data$DNAm.probe
ICC_ranks_S3 <- sort(ICC_ranks_S3, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S3 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S3,
                           nperm = 10000)

## Scenario 4 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("Horvath_Probes.csv")
keep <- mQTL_probes$DNAm.probe %in% data$DNAm.probe
mQTL_probes <- mQTL_probes[keep,]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S4 <- data[,c(1,5)]
ICC_ranks_S4 <- ICC_ranks_S4 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S4 <- as.numeric(ICC_ranks_S4$ICC_Scen4)
names(ICC_ranks_S4) <- data$DNAm.probe
ICC_ranks_S4 <- sort(ICC_ranks_S4, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S4 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S4,
                           nperm = 10000)

## Scenario 5 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("Horvath_Probes.csv")
keep <- mQTL_probes$DNAm.probe %in% data$DNAm.probe
mQTL_probes <- mQTL_probes[keep,]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S5 <- data[,c(1,6)]
ICC_ranks_S5 <- ICC_ranks_S5 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S5 <- as.numeric(ICC_ranks_S5$ICC_Scen5)
names(ICC_ranks_S5) <- data$DNAm.probe
ICC_ranks_S5 <- sort(ICC_ranks_S5, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S5 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S5,
                           nperm = 10000)

## Scenario 6 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("Horvath_Probes.csv")
keep <- mQTL_probes$DNAm.probe %in% data$DNAm.probe
mQTL_probes <- mQTL_probes[keep,]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S6 <- data[,c(1,7)]
ICC_ranks_S6 <- ICC_ranks_S6 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S6 <- as.numeric(ICC_ranks_S6$ICC_Scen6)
names(ICC_ranks_S6) <- data$DNAm.probe
ICC_ranks_S6 <- sort(ICC_ranks_S6, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S6 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S6,
                           nperm = 10000)

## Scenario 7 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("Horvath_Probes.csv")
keep <- mQTL_probes$DNAm.probe %in% data$DNAm.probe
mQTL_probes <- mQTL_probes[keep,]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S7 <- data[,c(1,8)]
ICC_ranks_S7 <- ICC_ranks_S7 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S7 <- as.numeric(ICC_ranks_S7$ICC_Scen7)
names(ICC_ranks_S7) <- data$DNAm.probe
ICC_ranks_S7 <- sort(ICC_ranks_S7, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S7 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S7,
                           nperm = 10000)

## Scenario 8 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("Horvath_Probes.csv")
keep <- mQTL_probes$DNAm.probe %in% data$DNAm.probe
mQTL_probes <- mQTL_probes[keep,]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S8 <- data[,c(1,9)]
ICC_ranks_S8 <- ICC_ranks_S8 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S8 <- as.numeric(ICC_ranks_S8$ICC_Scen8)
names(ICC_ranks_S8) <- data$DNAm.probe
ICC_ranks_S8 <- sort(ICC_ranks_S8, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S8 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S8,
                           nperm = 10000)

## Scenario 9 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("Horvath_Probes.csv")
keep <- mQTL_probes$DNAm.probe %in% data$DNAm.probe
mQTL_probes <- mQTL_probes[keep,]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S9 <- data[,c(1,10)]
ICC_ranks_S9 <- ICC_ranks_S9 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S9 <- as.numeric(ICC_ranks_S9$ICC_Scen9)
names(ICC_ranks_S9) <- data$DNAm.probe
ICC_ranks_S9 <- sort(ICC_ranks_S9, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S9 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S9,
                           nperm = 10000)

## Scenario 10 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("Horvath_Probes.csv")
keep <- mQTL_probes$DNAm.probe %in% data$DNAm.probe
mQTL_probes <- mQTL_probes[keep,]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S10 <- data[,c(1,11)]
ICC_ranks_S10 <- ICC_ranks_S10 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S10 <- as.numeric(ICC_ranks_S10$ICC_Scen10)
names(ICC_ranks_S10) <- data$DNAm.probe
ICC_ranks_S10 <- sort(ICC_ranks_S10, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S10 <- fgseaSimple(pathways = mQTL_probes, 
                            stats = ICC_ranks_S10,
                            nperm = 10000)


#####################
### Hannum Probes ###
#####################
### Probe Set Enrichment Analysis ### ----
## Scenario 1 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("Hannum_Probes.csv")
keep <- mQTL_probes$DNAm.probe %in% data$DNAm.probe
mQTL_probes <- mQTL_probes[keep,]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S1 <- data[,c(1,2)]
ICC_ranks_S1 <- ICC_ranks_S1 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S1 <- as.numeric(ICC_ranks_S1$ICC_Scen1)
names(ICC_ranks_S1) <- data$DNAm.probe
ICC_ranks_S1 <- sort(ICC_ranks_S1, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S1 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S1,
                           nperm = 10000)

## Scenario 2 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("Hannum_Probes.csv")
keep <- mQTL_probes$DNAm.probe %in% data$DNAm.probe
mQTL_probes <- mQTL_probes[keep,]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S2 <- data[,c(1,3)]
ICC_ranks_S2 <- ICC_ranks_S2 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S2 <- as.numeric(ICC_ranks_S2$ICC_Scen2)
names(ICC_ranks_S2) <- data$DNAm.probe
ICC_ranks_S2 <- sort(ICC_ranks_S2, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S2 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S2,
                           nperm = 10000)

## Scenario 3 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("Hannum_Probes.csv")
keep <- mQTL_probes$DNAm.probe %in% data$DNAm.probe
mQTL_probes <- mQTL_probes[keep,]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S3 <- data[,c(1,4)]
ICC_ranks_S3 <- ICC_ranks_S3 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S3 <- as.numeric(ICC_ranks_S3$ICC_Scen3)
names(ICC_ranks_S3) <- data$DNAm.probe
ICC_ranks_S3 <- sort(ICC_ranks_S3, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S3 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S3,
                           nperm = 10000)

## Scenario 4 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("Hannum_Probes.csv")
keep <- mQTL_probes$DNAm.probe %in% data$DNAm.probe
mQTL_probes <- mQTL_probes[keep,]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S4 <- data[,c(1,5)]
ICC_ranks_S4 <- ICC_ranks_S4 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S4 <- as.numeric(ICC_ranks_S4$ICC_Scen4)
names(ICC_ranks_S4) <- data$DNAm.probe
ICC_ranks_S4 <- sort(ICC_ranks_S4, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S4 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S4,
                           nperm = 10000)

## Scenario 5 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("Hannum_Probes.csv")
keep <- mQTL_probes$DNAm.probe %in% data$DNAm.probe
mQTL_probes <- mQTL_probes[keep,]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S5 <- data[,c(1,6)]
ICC_ranks_S5 <- ICC_ranks_S5 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S5 <- as.numeric(ICC_ranks_S5$ICC_Scen5)
names(ICC_ranks_S5) <- data$DNAm.probe
ICC_ranks_S5 <- sort(ICC_ranks_S5, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S5 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S5,
                           nperm = 10000)

## Scenario 6 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("Hannum_Probes.csv")
keep <- mQTL_probes$DNAm.probe %in% data$DNAm.probe
mQTL_probes <- mQTL_probes[keep,]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S6 <- data[,c(1,7)]
ICC_ranks_S6 <- ICC_ranks_S6 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S6 <- as.numeric(ICC_ranks_S6$ICC_Scen6)
names(ICC_ranks_S6) <- data$DNAm.probe
ICC_ranks_S6 <- sort(ICC_ranks_S6, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S6 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S6,
                           nperm = 10000)

## Scenario 7 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("Hannum_Probes.csv")
keep <- mQTL_probes$DNAm.probe %in% data$DNAm.probe
mQTL_probes <- mQTL_probes[keep,]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S7 <- data[,c(1,8)]
ICC_ranks_S7 <- ICC_ranks_S7 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S7 <- as.numeric(ICC_ranks_S7$ICC_Scen7)
names(ICC_ranks_S7) <- data$DNAm.probe
ICC_ranks_S7 <- sort(ICC_ranks_S7, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S7 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S7,
                           nperm = 10000)

## Scenario 8 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("Hannum_Probes.csv")
keep <- mQTL_probes$DNAm.probe %in% data$DNAm.probe
mQTL_probes <- mQTL_probes[keep,]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S8 <- data[,c(1,9)]
ICC_ranks_S8 <- ICC_ranks_S8 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S8 <- as.numeric(ICC_ranks_S8$ICC_Scen8)
names(ICC_ranks_S8) <- data$DNAm.probe
ICC_ranks_S8 <- sort(ICC_ranks_S8, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S8 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S8,
                           nperm = 10000)

## Scenario 9 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("Hannum_Probes.csv")
keep <- mQTL_probes$DNAm.probe %in% data$DNAm.probe
mQTL_probes <- mQTL_probes[keep,]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S9 <- data[,c(1,10)]
ICC_ranks_S9 <- ICC_ranks_S9 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S9 <- as.numeric(ICC_ranks_S9$ICC_Scen9)
names(ICC_ranks_S9) <- data$DNAm.probe
ICC_ranks_S9 <- sort(ICC_ranks_S9, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S9 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S9,
                           nperm = 10000)

## Scenario 10 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("Hannum_Probes.csv")
keep <- mQTL_probes$DNAm.probe %in% data$DNAm.probe
mQTL_probes <- mQTL_probes[keep,]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S10 <- data[,c(1,11)]
ICC_ranks_S10 <- ICC_ranks_S10 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S10 <- as.numeric(ICC_ranks_S10$ICC_Scen10)
names(ICC_ranks_S10) <- data$DNAm.probe
ICC_ranks_S10 <- sort(ICC_ranks_S10, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S10 <- fgseaSimple(pathways = mQTL_probes, 
                            stats = ICC_ranks_S10,
                            nperm = 10000)

#######################
### PhenoAge Probes ###
#######################
### Probe Set Enrichment Analysis ### ----
## Scenario 1 ##
# Read in significant mQTL probes #
PhenoAge_probes <- read.csv("PhenoAge_Probes.csv")
keep <- PhenoAge_probes$DNAm.probe %in% data$DNAm.probe
PhenoAge_probes <- PhenoAge_probes[keep,]
PhenoAge_probes <- list(mQTL = PhenoAge_probes)
ICC_ranks_S1 <- data[,c(1,2)]
ICC_ranks_S1 <- ICC_ranks_S1 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S1 <- as.numeric(ICC_ranks_S1$ICC_Scen1)
names(ICC_ranks_S1) <- data$DNAm.probe
ICC_ranks_S1 <- sort(ICC_ranks_S1, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S1 <- fgseaSimple(pathways = PhenoAge_probes, 
                           stats = ICC_ranks_S1,
                           nperm = 10000)

## Scenario 2 ##
# Read in significant mQTL probes #
PhenoAge_probes <- read.csv("PhenoAge_Probes.csv")
keep <- PhenoAge_probes$DNAm.probe %in% PhenoAge_probes$DNAm.probe
PhenoAge_probes <- PhenoAge_probes[keep,]
PhenoAge_probes <- list(mQTL = PhenoAge_probes)
ICC_ranks_S2 <- data[,c(1,3)]
ICC_ranks_S2 <- ICC_ranks_S2 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S2 <- as.numeric(ICC_ranks_S2$ICC_Scen2)
names(ICC_ranks_S2) <- data$DNAm.probe
ICC_ranks_S2 <- sort(ICC_ranks_S2, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S2 <- fgseaSimple(pathways = PhenoAge_probes, 
                           stats = ICC_ranks_S2,
                           nperm = 10000)

## Scenario 3 ##
# Read in significant mQTL probes #
PhenoAge_probes <- read.csv("PhenoAge_Probes.csv")
keep <- PhenoAge_probes$DNAm.probe %in% PhenoAge_probes$DNAm.probe
PhenoAge_probes <- PhenoAge_probes[keep,]
PhenoAge_probes <- list(mQTL = PhenoAge_probes)
ICC_ranks_S3 <- data[,c(1,4)]
ICC_ranks_S3 <- ICC_ranks_S3 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S3 <- as.numeric(ICC_ranks_S3$ICC_Scen3)
names(ICC_ranks_S3) <- data$DNAm.probe
ICC_ranks_S3 <- sort(ICC_ranks_S3, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S3 <- fgseaSimple(pathways = PhenoAge_probes, 
                           stats = ICC_ranks_S3,
                           nperm = 10000)

## Scenario 4 ##
# Read in significant mQTL probes #
PhenoAge_probes <- read.csv("PhenoAge_Probes.csv")
keep <- PhenoAge_probes$DNAm.probe %in% PhenoAge_probes$DNAm.probe
PhenoAge_probes <- PhenoAge_probes[keep,]
PhenoAge_probes <- list(mQTL = PhenoAge_probes)
ICC_ranks_S4 <- data[,c(1,5)]
ICC_ranks_S4 <- ICC_ranks_S4 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S4 <- as.numeric(ICC_ranks_S4$ICC_Scen4)
names(ICC_ranks_S4) <- data$DNAm.probe
ICC_ranks_S4 <- sort(ICC_ranks_S4, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S4 <- fgseaSimple(pathways = PhenoAge_probes, 
                           stats = ICC_ranks_S4,
                           nperm = 10000)

## Scenario 5 ##
# Read in significant mQTL probes #
PhenoAge_probes <- read.csv("PhenoAge_Probes.csv")
keep <- PhenoAge_probes$DNAm.probe %in% PhenoAge_probes$DNAm.probe
PhenoAge_probes <- PhenoAge_probes[keep,]
PhenoAge_probes <- list(mQTL = PhenoAge_probes)
ICC_ranks_S5 <- data[,c(1,6)]
ICC_ranks_S5 <- ICC_ranks_S5 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S5 <- as.numeric(ICC_ranks_S5$ICC_Scen5)
names(ICC_ranks_S5) <- data$DNAm.probe
ICC_ranks_S5 <- sort(ICC_ranks_S5, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S5 <- fgseaSimple(pathways = PhenoAge_probes, 
                           stats = ICC_ranks_S5,
                           nperm = 10000)

## Scenario 6 ##
# Read in significant mQTL probes #
PhenoAge_probes <- read.csv("PhenoAge_Probes.csv")
keep <- PhenoAge_probes$DNAm.probe %in% PhenoAge_probes$DNAm.probe
PhenoAge_probes <- PhenoAge_probes[keep,]
PhenoAge_probes <- list(mQTL = PhenoAge_probes)
ICC_ranks_S6 <- data[,c(1,7)]
ICC_ranks_S6 <- ICC_ranks_S6 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S6 <- as.numeric(ICC_ranks_S6$ICC_Scen6)
names(ICC_ranks_S6) <- data$DNAm.probe
ICC_ranks_S6 <- sort(ICC_ranks_S6, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S6 <- fgseaSimple(pathways = PhenoAge_probes, 
                           stats = ICC_ranks_S6,
                           nperm = 10000)

## Scenario 7 ##
# Read in significant mQTL probes #
PhenoAge_probes <- read.csv("PhenoAge_Probes.csv")
keep <- PhenoAge_probes$DNAm.probe %in% PhenoAge_probes$DNAm.probe
PhenoAge_probes <- PhenoAge_probes[keep,]
PhenoAge_probes <- list(mQTL = PhenoAge_probes)
ICC_ranks_S7 <- data[,c(1,8)]
ICC_ranks_S7 <- ICC_ranks_S7 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S7 <- as.numeric(ICC_ranks_S7$ICC_Scen7)
names(ICC_ranks_S7) <- data$DNAm.probe
ICC_ranks_S7 <- sort(ICC_ranks_S7, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S7 <- fgseaSimple(pathways = PhenoAge_probes, 
                           stats = ICC_ranks_S7,
                           nperm = 10000)

## Scenario 8 ##
# Read in significant mQTL probes #
PhenoAge_probes <- read.csv("PhenoAge_Probes.csv")
keep <- PhenoAge_probes$DNAm.probe %in% PhenoAge_probes$DNAm.probe
PhenoAge_probes <- PhenoAge_probes[keep,]
PhenoAge_probes <- list(mQTL = PhenoAge_probes)
ICC_ranks_S8 <- data[,c(1,9)]
ICC_ranks_S8 <- ICC_ranks_S8 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S8 <- as.numeric(ICC_ranks_S8$ICC_Scen8)
names(ICC_ranks_S8) <- data$DNAm.probe
ICC_ranks_S8 <- sort(ICC_ranks_S8, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S8 <- fgseaSimple(pathways = PhenoAge_probes, 
                           stats = ICC_ranks_S8,
                           nperm = 10000)

## Scenario 9 ##
# Read in significant mQTL probes #
PhenoAge_probes <- read.csv("PhenoAge_Probes.csv")
keep <- PhenoAge_probes$DNAm.probe %in% PhenoAge_probes$DNAm.probe
PhenoAge_probes <- PhenoAge_probes[keep,]
PhenoAge_probes <- list(mQTL = PhenoAge_probes)
ICC_ranks_S9 <- data[,c(1,10)]
ICC_ranks_S9 <- ICC_ranks_S9 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S9 <- as.numeric(ICC_ranks_S9$ICC_Scen9)
names(ICC_ranks_S9) <- data$DNAm.probe
ICC_ranks_S9 <- sort(ICC_ranks_S9, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S9 <- fgseaSimple(pathways = PhenoAge_probes, 
                           stats = ICC_ranks_S9,
                           nperm = 10000)

## Scenario 10 ##
# Read in significant mQTL probes #
PhenoAge_probes <- read.csv("PhenoAge_Probes.csv")
keep <- PhenoAge_probes$DNAm.probe %in% PhenoAge_probes$DNAm.probe
PhenoAge_probes <- PhenoAge_probes[keep,]
PhenoAge_probes <- list(mQTL = PhenoAge_probes)
ICC_ranks_S10 <- data[,c(1,11)]
ICC_ranks_S10 <- ICC_ranks_S10 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S10 <- as.numeric(ICC_ranks_S10$ICC_Scen10)
names(ICC_ranks_S10) <- data$DNAm.probe
ICC_ranks_S10 <- sort(ICC_ranks_S10, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S10 <- fgseaSimple(pathways = PhenoAge_probes, 
                            stats = ICC_ranks_S10,
                            nperm = 10000)

#######################
### PC Clock Probes ###
#######################
### Probe Set Enrichment Analysis ### ----
## Scenario 1 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("PC_Probes.csv")
mQTL_probes <- as.data.frame(sample(mQTL_probes$DNAm.probe, 35000))
colnames(mQTL_probes) <- "DNAm.probe"
keep <- mQTL_probes$DNAm.probe %in% data$DNAm.probe
mQTL_probes <- mQTL_probes[keep,]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S1 <- data[,c(1,2)]
ICC_ranks_S1 <- ICC_ranks_S1 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S1 <- as.numeric(ICC_ranks_S1$ICC_Scen1)
names(ICC_ranks_S1) <- data$DNAm.probe
ICC_ranks_S1 <- sort(ICC_ranks_S1, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S1 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S1,
                           nperm = 10000)

## Scenario 2 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("PC_Probes.csv")
mQTL_probes <- as.data.frame(sample(mQTL_probes$DNAm.probe, 35000))
colnames(mQTL_probes) <- "DNAm.probe"
keep <- mQTL_probes$DNAm.probe %in% data$DNAm.probe
mQTL_probes <- mQTL_probes[keep,]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S2 <- data[,c(1,3)]
ICC_ranks_S2 <- ICC_ranks_S2 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S2 <- as.numeric(ICC_ranks_S2$ICC_Scen2)
names(ICC_ranks_S2) <- data$DNAm.probe
ICC_ranks_S2 <- sort(ICC_ranks_S2, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S2 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S2,
                           nperm = 10000)

## Scenario 3 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("PC_Probes.csv")
mQTL_probes <- as.data.frame(sample(mQTL_probes$DNAm.probe, 35000))
colnames(mQTL_probes) <- "DNAm.probe"
keep <- mQTL_probes$DNAm.probe %in% data$DNAm.probe
mQTL_probes <- mQTL_probes[keep,]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S3 <- data[,c(1,4)]
ICC_ranks_S3 <- ICC_ranks_S3 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S3 <- as.numeric(ICC_ranks_S3$ICC_Scen3)
names(ICC_ranks_S3) <- data$DNAm.probe
ICC_ranks_S3 <- sort(ICC_ranks_S3, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S3 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S3,
                           nperm = 10000)

## Scenario 4 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("PC_Probes.csv")
mQTL_probes <- as.data.frame(sample(mQTL_probes$DNAm.probe, 35000))
colnames(mQTL_probes) <- "DNAm.probe"
keep <- mQTL_probes$DNAm.probe %in% data$DNAm.probe
mQTL_probes <- mQTL_probes[keep,]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S4 <- data[,c(1,5)]
ICC_ranks_S4 <- ICC_ranks_S4 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S4 <- as.numeric(ICC_ranks_S4$ICC_Scen4)
names(ICC_ranks_S4) <- data$DNAm.probe
ICC_ranks_S4 <- sort(ICC_ranks_S4, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S4 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S4,
                           nperm = 10000)

## Scenario 5 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("PC_Probes.csv")
mQTL_probes <- as.data.frame(sample(mQTL_probes$DNAm.probe, 35000))
colnames(mQTL_probes) <- "DNAm.probe"
keep <- mQTL_probes$DNAm.probe %in% data$DNAm.probe
mQTL_probes <- mQTL_probes[keep,]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S5 <- data[,c(1,6)]
ICC_ranks_S5 <- ICC_ranks_S5 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S5 <- as.numeric(ICC_ranks_S5$ICC_Scen5)
names(ICC_ranks_S5) <- data$DNAm.probe
ICC_ranks_S5 <- sort(ICC_ranks_S5, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S5 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S5,
                           nperm = 10000)

## Scenario 6 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("PC_Probes.csv")
mQTL_probes <- as.data.frame(sample(mQTL_probes$DNAm.probe, 35000))
colnames(mQTL_probes) <- "DNAm.probe"
keep <- mQTL_probes$DNAm.probe %in% data$DNAm.probe
mQTL_probes <- mQTL_probes[keep,]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S6 <- data[,c(1,7)]
ICC_ranks_S6 <- ICC_ranks_S6 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S6 <- as.numeric(ICC_ranks_S6$ICC_Scen6)
names(ICC_ranks_S6) <- data$DNAm.probe
ICC_ranks_S6 <- sort(ICC_ranks_S6, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S6 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S6,
                           nperm = 10000)

## Scenario 7 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("PC_Probes.csv")
mQTL_probes <- as.data.frame(sample(mQTL_probes$DNAm.probe, 35000))
colnames(mQTL_probes) <- "DNAm.probe"
keep <- mQTL_probes$DNAm.probe %in% data$DNAm.probe
mQTL_probes <- mQTL_probes[keep,]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S7 <- data[,c(1,8)]
ICC_ranks_S7 <- ICC_ranks_S7 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S7 <- as.numeric(ICC_ranks_S7$ICC_Scen7)
names(ICC_ranks_S7) <- data$DNAm.probe
ICC_ranks_S7 <- sort(ICC_ranks_S7, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S7 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S7,
                           nperm = 10000)

## Scenario 8 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("PC_Probes.csv")
mQTL_probes <- as.data.frame(sample(mQTL_probes$DNAm.probe, 35000))
colnames(mQTL_probes) <- "DNAm.probe"
keep <- mQTL_probes$DNAm.probe %in% data$DNAm.probe
mQTL_probes <- mQTL_probes[keep,]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S8 <- data[,c(1,9)]
ICC_ranks_S8 <- ICC_ranks_S8 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S8 <- as.numeric(ICC_ranks_S8$ICC_Scen8)
names(ICC_ranks_S8) <- data$DNAm.probe
ICC_ranks_S8 <- sort(ICC_ranks_S8, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S8 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S8,
                           nperm = 10000)

## Scenario 9 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("PC_Probes.csv")
mQTL_probes <- as.data.frame(sample(mQTL_probes$DNAm.probe, 35000))
colnames(mQTL_probes) <- "DNAm.probe"
keep <- mQTL_probes$DNAm.probe %in% data$DNAm.probe
mQTL_probes <- mQTL_probes[keep,]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S9 <- data[,c(1,10)]
ICC_ranks_S9 <- ICC_ranks_S9 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S9 <- as.numeric(ICC_ranks_S9$ICC_Scen9)
names(ICC_ranks_S9) <- data$DNAm.probe
ICC_ranks_S9 <- sort(ICC_ranks_S9, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S9 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S9,
                           nperm = 10000)

## Scenario 10 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("PC_Probes.csv")
mQTL_probes <- as.data.frame(sample(mQTL_probes$DNAm.probe, 35000))
colnames(mQTL_probes) <- "DNAm.probe"
keep <- mQTL_probes$DNAm.probe %in% data$DNAm.probe
mQTL_probes <- mQTL_probes[keep,]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S10 <- data[,c(1,11)]
ICC_ranks_S10 <- ICC_ranks_S10 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S10 <- as.numeric(ICC_ranks_S10$ICC_Scen10)
names(ICC_ranks_S10) <- data$DNAm.probe
ICC_ranks_S10 <- sort(ICC_ranks_S10, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S10 <- fgseaSimple(pathways = mQTL_probes, 
                            stats = ICC_ranks_S10,
                            nperm = 10000)

###################
### PACE Probes ###
###################
### Probe Set Enrichment Analysis ### ----
## Scenario 1 ##
# Read in significant mQTL probes #
load("sysdata.rda")
keep <- mPOA_Models$model_probes$DunedinPACE %in% data$DNAm.probe
mQTL_probes <- mPOA_Models$model_probes$DunedinPACE[keep]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S1 <- data[,c(1,2)]
ICC_ranks_S1 <- ICC_ranks_S1 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S1 <- as.numeric(ICC_ranks_S1$ICC_Scen1)
names(ICC_ranks_S1) <- data$DNAm.probe
ICC_ranks_S1 <- sort(ICC_ranks_S1, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S1 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S1,
                           nperm = 10000)

## Scenario 2 ##
# Read in significant mQTL probes #
load("sysdata.rda")
keep <- mPOA_Models$model_probes$DunedinPACE %in% data$DNAm.probe
mQTL_probes <- mPOA_Models$model_probes$DunedinPACE[keep]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S2 <- data[,c(1,3)]
ICC_ranks_S2 <- ICC_ranks_S2 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S2 <- as.numeric(ICC_ranks_S2$ICC_Scen2)
names(ICC_ranks_S2) <- data$DNAm.probe
ICC_ranks_S2 <- sort(ICC_ranks_S2, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S2 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S2,
                           nperm = 10000)

## Scenario 3 ##
# Read in significant mQTL probes #
load("sysdata.rda")
keep <- mPOA_Models$model_probes$DunedinPACE %in% data$DNAm.probe
mQTL_probes <- mPOA_Models$model_probes$DunedinPACE[keep]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S3 <- data[,c(1,4)]
ICC_ranks_S3 <- ICC_ranks_S3 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S3 <- as.numeric(ICC_ranks_S3$ICC_Scen3)
names(ICC_ranks_S3) <- data$DNAm.probe
ICC_ranks_S3 <- sort(ICC_ranks_S3, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S3 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S3,
                           nperm = 10000)

## Scenario 4 ##
# Read in significant mQTL probes #
load("sysdata.rda")
keep <- mPOA_Models$model_probes$DunedinPACE %in% data$DNAm.probe
mQTL_probes <- mPOA_Models$model_probes$DunedinPACE[keep]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S4 <- data[,c(1,5)]
ICC_ranks_S4 <- ICC_ranks_S4 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S4 <- as.numeric(ICC_ranks_S4$ICC_Scen4)
names(ICC_ranks_S4) <- data$DNAm.probe
ICC_ranks_S4 <- sort(ICC_ranks_S4, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S4 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S4,
                           nperm = 10000)

## Scenario 5 ##
# Read in significant mQTL probes #
load("sysdata.rda")
keep <- mPOA_Models$model_probes$DunedinPACE %in% data$DNAm.probe
mQTL_probes <- mPOA_Models$model_probes$DunedinPACE[keep]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S5 <- data[,c(1,6)]
ICC_ranks_S5 <- ICC_ranks_S5 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S5 <- as.numeric(ICC_ranks_S5$ICC_Scen5)
names(ICC_ranks_S5) <- data$DNAm.probe
ICC_ranks_S5 <- sort(ICC_ranks_S5, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S5 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S5,
                           nperm = 10000)

## Scenario 6 ##
# Read in significant mQTL probes #
load("sysdata.rda")
keep <- mPOA_Models$model_probes$DunedinPACE %in% data$DNAm.probe
mQTL_probes <- mPOA_Models$model_probes$DunedinPACE[keep]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S6 <- data[,c(1,7)]
ICC_ranks_S6 <- ICC_ranks_S6 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S6 <- as.numeric(ICC_ranks_S6$ICC_Scen6)
names(ICC_ranks_S6) <- data$DNAm.probe
ICC_ranks_S6 <- sort(ICC_ranks_S6, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S6 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S6,
                           nperm = 10000)

## Scenario 7 ##
# Read in significant mQTL probes #
load("sysdata.rda")
keep <- mPOA_Models$model_probes$DunedinPACE %in% data$DNAm.probe
mQTL_probes <- mPOA_Models$model_probes$DunedinPACE[keep]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S7 <- data[,c(1,8)]
ICC_ranks_S7 <- ICC_ranks_S7 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S7 <- as.numeric(ICC_ranks_S7$ICC_Scen7)
names(ICC_ranks_S7) <- data$DNAm.probe
ICC_ranks_S7 <- sort(ICC_ranks_S7, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S7 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S7,
                           nperm = 10000)

## Scenario 8 ##
# Read in significant mQTL probes #
load("sysdata.rda")
keep <- mPOA_Models$model_probes$DunedinPACE %in% data$DNAm.probe
mQTL_probes <- mPOA_Models$model_probes$DunedinPACE[keep]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S8 <- data[,c(1,9)]
ICC_ranks_S8 <- ICC_ranks_S8 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S8 <- as.numeric(ICC_ranks_S8$ICC_Scen8)
names(ICC_ranks_S8) <- data$DNAm.probe
ICC_ranks_S8 <- sort(ICC_ranks_S8, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S8 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S8,
                           nperm = 10000)

## Scenario 9 ##
# Read in significant mQTL probes #
load("sysdata.rda")
keep <- mPOA_Models$model_probes$DunedinPACE %in% data$DNAm.probe
mQTL_probes <- mPOA_Models$model_probes$DunedinPACE[keep]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S9 <- data[,c(1,10)]
ICC_ranks_S9 <- ICC_ranks_S9 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S9 <- as.numeric(ICC_ranks_S9$ICC_Scen9)
names(ICC_ranks_S9) <- data$DNAm.probe
ICC_ranks_S9 <- sort(ICC_ranks_S9, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S9 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S9,
                           nperm = 10000)

## Scenario 10 ##
# Read in significant mQTL probes #
load("sysdata.rda")
keep <- mPOA_Models$model_probes$DunedinPACE %in% data$DNAm.probe
mQTL_probes <- mPOA_Models$model_probes$DunedinPACE[keep]
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S10 <- data[,c(1,11)]
ICC_ranks_S10 <- ICC_ranks_S10 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S10 <- as.numeric(ICC_ranks_S10$ICC_Scen10)
names(ICC_ranks_S10) <- data$DNAm.probe
ICC_ranks_S10 <- sort(ICC_ranks_S10, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S10 <- fgseaSimple(pathways = mQTL_probes, 
                            stats = ICC_ranks_S10,
                            nperm = 10000)

#####################
### Old Salas CpG ###
#####################
### Probe Set Enrichment Analysis ### ----
library(FlowSorted.Blood.EPIC)
library(FlowSorted.BloodExtended.EPIC)
## Scenario 1 ##
# Read in significant mQTL probes #
data(IDOLOptimizedCpGs)
mQTL_probes <- list(mQTL = IDOLOptimizedCpGs)
ICC_ranks_S1 <- data[,c(1,2)]
ICC_ranks_S1 <- ICC_ranks_S1 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S1 <- as.numeric(ICC_ranks_S1$ICC_Scen1)
names(ICC_ranks_S1) <- data$DNAm.probe
ICC_ranks_S1 <- sort(ICC_ranks_S1, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S1 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S1,
                           nperm = 10000)

## Scenario 2 ##
# Read in significant mQTL probes #
data(IDOLOptimizedCpGs)
mQTL_probes <- list(mQTL = IDOLOptimizedCpGs)
ICC_ranks_S2 <- data[,c(1,3)]
ICC_ranks_S2 <- ICC_ranks_S2 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S2 <- as.numeric(ICC_ranks_S2$ICC_Scen2)
names(ICC_ranks_S2) <- data$DNAm.probe
ICC_ranks_S2 <- sort(ICC_ranks_S2, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S2 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S2,
                           nperm = 10000)

## Scenario 3 ##
# Read in significant mQTL probes #
data(IDOLOptimizedCpGs)
mQTL_probes <- list(mQTL = IDOLOptimizedCpGs)
ICC_ranks_S3 <- data[,c(1,4)]
ICC_ranks_S3 <- ICC_ranks_S3 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S3 <- as.numeric(ICC_ranks_S3$ICC_Scen3)
names(ICC_ranks_S3) <- data$DNAm.probe
ICC_ranks_S3 <- sort(ICC_ranks_S3, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S3 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S3,
                           nperm = 10000)

## Scenario 4 ##
# Read in significant mQTL probes #
data(IDOLOptimizedCpGs)
mQTL_probes <- list(mQTL = IDOLOptimizedCpGs)
ICC_ranks_S4 <- data[,c(1,5)]
ICC_ranks_S4 <- ICC_ranks_S4 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S4 <- as.numeric(ICC_ranks_S4$ICC_Scen4)
names(ICC_ranks_S4) <- data$DNAm.probe
ICC_ranks_S4 <- sort(ICC_ranks_S4, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S4 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S4,
                           nperm = 10000)

## Scenario 5 ##
# Read in significant mQTL probes #
data(IDOLOptimizedCpGs)
mQTL_probes <- list(mQTL = IDOLOptimizedCpGs)
ICC_ranks_S5 <- data[,c(1,6)]
ICC_ranks_S5 <- ICC_ranks_S5 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S5 <- as.numeric(ICC_ranks_S5$ICC_Scen5)
names(ICC_ranks_S5) <- data$DNAm.probe
ICC_ranks_S5 <- sort(ICC_ranks_S5, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S5 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S5,
                           nperm = 10000)

## Scenario 6 ##
# Read in significant mQTL probes #
data(IDOLOptimizedCpGs)
mQTL_probes <- list(mQTL = IDOLOptimizedCpGs)
ICC_ranks_S6 <- data[,c(1,7)]
ICC_ranks_S6 <- ICC_ranks_S6 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S6 <- as.numeric(ICC_ranks_S6$ICC_Scen6)
names(ICC_ranks_S6) <- data$DNAm.probe
ICC_ranks_S6 <- sort(ICC_ranks_S6, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S6 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S6,
                           nperm = 10000)

## Scenario 7 ##
# Read in significant mQTL probes #
data(IDOLOptimizedCpGs)
mQTL_probes <- list(mQTL = IDOLOptimizedCpGs)
ICC_ranks_S7 <- data[,c(1,8)]
ICC_ranks_S7 <- ICC_ranks_S7 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S7 <- as.numeric(ICC_ranks_S7$ICC_Scen7)
names(ICC_ranks_S7) <- data$DNAm.probe
ICC_ranks_S7 <- sort(ICC_ranks_S7, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S7 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S7,
                           nperm = 10000)

## Scenario 8 ##
# Read in significant mQTL probes #
data(IDOLOptimizedCpGs)
mQTL_probes <- list(mQTL = IDOLOptimizedCpGs)
ICC_ranks_S8 <- data[,c(1,9)]
ICC_ranks_S8 <- ICC_ranks_S8 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S8 <- as.numeric(ICC_ranks_S8$ICC_Scen8)
names(ICC_ranks_S8) <- data$DNAm.probe
ICC_ranks_S8 <- sort(ICC_ranks_S8, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S8 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S8,
                           nperm = 10000)

## Scenario 9 ##
# Read in significant mQTL probes #
data(IDOLOptimizedCpGs)
mQTL_probes <- list(mQTL = IDOLOptimizedCpGs)
ICC_ranks_S9 <- data[,c(1,10)]
ICC_ranks_S9 <- ICC_ranks_S9 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S9 <- as.numeric(ICC_ranks_S9$ICC_Scen9)
names(ICC_ranks_S9) <- data$DNAm.probe
ICC_ranks_S9 <- sort(ICC_ranks_S9, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S9 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S9,
                           nperm = 10000)

## Scenario 10 ##
# Read in significant mQTL probes #
data(IDOLOptimizedCpGs)
mQTL_probes <- list(mQTL = IDOLOptimizedCpGs)
ICC_ranks_S10 <- data[,c(1,11)]
ICC_ranks_S10 <- ICC_ranks_S10 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S10 <- as.numeric(ICC_ranks_S10$ICC_Scen10)
names(ICC_ranks_S10) <- data$DNAm.probe
ICC_ranks_S10 <- sort(ICC_ranks_S10, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S10 <- fgseaSimple(pathways = mQTL_probes, 
                            stats = ICC_ranks_S10,
                            nperm = 10000)

#####################
### IDOL-Extended ###
#####################
### Probe Set Enrichment Analysis ### ----
## Scenario 1 ##
# Read in significant mQTL probes #
data(IDOLOptimizedCpGsBloodExtended)
mQTL_probes <- list(mQTL = IDOLOptimizedCpGsBloodExtended)
ICC_ranks_S1 <- data[,c(1,2)]
ICC_ranks_S1 <- ICC_ranks_S1 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S1 <- as.numeric(ICC_ranks_S1$ICC_Scen1)
names(ICC_ranks_S1) <- data$DNAm.probe
ICC_ranks_S1 <- sort(ICC_ranks_S1, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S1 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S1,
                           nperm = 10000)

## Scenario 2 ##
# Read in significant mQTL probes #
data(IDOLOptimizedCpGsBloodExtended)
mQTL_probes <- list(mQTL = IDOLOptimizedCpGsBloodExtended)
ICC_ranks_S2 <- data[,c(1,3)]
ICC_ranks_S2 <- ICC_ranks_S2 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S2 <- as.numeric(ICC_ranks_S2$ICC_Scen2)
names(ICC_ranks_S2) <- data$DNAm.probe
ICC_ranks_S2 <- sort(ICC_ranks_S2, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S2 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S2,
                           nperm = 10000)

## Scenario 3 ##
# Read in significant mQTL probes #
data(IDOLOptimizedCpGsBloodExtended)
mQTL_probes <- list(mQTL = IDOLOptimizedCpGsBloodExtended)
ICC_ranks_S3 <- data[,c(1,4)]
ICC_ranks_S3 <- ICC_ranks_S3 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S3 <- as.numeric(ICC_ranks_S3$ICC_Scen3)
names(ICC_ranks_S3) <- data$DNAm.probe
ICC_ranks_S3 <- sort(ICC_ranks_S3, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S3 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S3,
                           nperm = 10000)

## Scenario 4 ##
# Read in significant mQTL probes #
data(IDOLOptimizedCpGsBloodExtended)
mQTL_probes <- list(mQTL = IDOLOptimizedCpGsBloodExtended)
ICC_ranks_S4 <- data[,c(1,5)]
ICC_ranks_S4 <- ICC_ranks_S4 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S4 <- as.numeric(ICC_ranks_S4$ICC_Scen4)
names(ICC_ranks_S4) <- data$DNAm.probe
ICC_ranks_S4 <- sort(ICC_ranks_S4, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S4 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S4,
                           nperm = 10000)

## Scenario 5 ##
# Read in significant mQTL probes #
data(IDOLOptimizedCpGsBloodExtended)
mQTL_probes <- list(mQTL = IDOLOptimizedCpGsBloodExtended)
ICC_ranks_S5 <- data[,c(1,6)]
ICC_ranks_S5 <- ICC_ranks_S5 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S5 <- as.numeric(ICC_ranks_S5$ICC_Scen5)
names(ICC_ranks_S5) <- data$DNAm.probe
ICC_ranks_S5 <- sort(ICC_ranks_S5, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S5 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S5,
                           nperm = 10000)

## Scenario 6 ##
# Read in significant mQTL probes #
data(IDOLOptimizedCpGsBloodExtended)
mQTL_probes <- list(mQTL = IDOLOptimizedCpGsBloodExtended)
ICC_ranks_S6 <- data[,c(1,7)]
ICC_ranks_S6 <- ICC_ranks_S6 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S6 <- as.numeric(ICC_ranks_S6$ICC_Scen6)
names(ICC_ranks_S6) <- data$DNAm.probe
ICC_ranks_S6 <- sort(ICC_ranks_S6, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S6 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S6,
                           nperm = 10000)

## Scenario 7 ##
# Read in significant mQTL probes #
data(IDOLOptimizedCpGsBloodExtended)
mQTL_probes <- list(mQTL = IDOLOptimizedCpGsBloodExtended)
ICC_ranks_S7 <- data[,c(1,8)]
ICC_ranks_S7 <- ICC_ranks_S7 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S7 <- as.numeric(ICC_ranks_S7$ICC_Scen7)
names(ICC_ranks_S7) <- data$DNAm.probe
ICC_ranks_S7 <- sort(ICC_ranks_S7, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S7 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S7,
                           nperm = 10000)

## Scenario 8 ##
# Read in significant mQTL probes #
data(IDOLOptimizedCpGsBloodExtended)
mQTL_probes <- list(mQTL = IDOLOptimizedCpGsBloodExtended)
ICC_ranks_S8 <- data[,c(1,9)]
ICC_ranks_S8 <- ICC_ranks_S8 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S8 <- as.numeric(ICC_ranks_S8$ICC_Scen8)
names(ICC_ranks_S8) <- data$DNAm.probe
ICC_ranks_S8 <- sort(ICC_ranks_S8, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S8 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S8,
                           nperm = 10000)

## Scenario 9 ##
# Read in significant mQTL probes #
data(IDOLOptimizedCpGsBloodExtended)
mQTL_probes <- list(mQTL = IDOLOptimizedCpGsBloodExtended)
ICC_ranks_S9 <- data[,c(1,10)]
ICC_ranks_S9 <- ICC_ranks_S9 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S9 <- as.numeric(ICC_ranks_S9$ICC_Scen9)
names(ICC_ranks_S9) <- data$DNAm.probe
ICC_ranks_S9 <- sort(ICC_ranks_S9, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S9 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S9,
                           nperm = 10000)

## Scenario 10 ##
# Read in significant mQTL probes #
data(IDOLOptimizedCpGsBloodExtended)
mQTL_probes <- list(mQTL = IDOLOptimizedCpGsBloodExtended)
ICC_ranks_S10 <- data[,c(1,11)]
ICC_ranks_S10 <- ICC_ranks_S10 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S10 <- as.numeric(ICC_ranks_S10$ICC_Scen10)
names(ICC_ranks_S10) <- data$DNAm.probe
ICC_ranks_S10 <- sort(ICC_ranks_S10, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S10 <- fgseaSimple(pathways = mQTL_probes, 
                            stats = ICC_ranks_S10,
                            nperm = 10000)



##################################################################################################
################################## Biologically Relevant Probes ##################################
##################################################################################################
### Set Seed ### ----
set.seed(100)

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
# Exclude probes with ICC < 0.01 #
data[,2:11][data[,2:11] < 0.01] <- NA

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


#####################
#### mQTL Probes ####
#####################
### Probe Set Enrichment Analysis ### ----
## Scenario 1 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("mQTL_Probes.csv")
mQTL_probes <- sample(mQTL_probes$DNAm.probe, 35000)
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S1 <- data[,c(1,2)]
ICC_ranks_S1 <- ICC_ranks_S1 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S1 <- as.numeric(ICC_ranks_S1$ICC_Scen1)
names(ICC_ranks_S1) <- data$DNAm.probe
ICC_ranks_S1 <- sort(ICC_ranks_S1, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S1 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S1,
                           nperm = 10000)

## Scenario 2 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("mQTL_Probes.csv")
mQTL_probes <- sample(mQTL_probes$DNAm.probe, 35000)
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S2 <- data[,c(1,3)]
ICC_ranks_S2 <- ICC_ranks_S2 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S2 <- as.numeric(ICC_ranks_S2$ICC_Scen2)
names(ICC_ranks_S2) <- data$DNAm.probe
ICC_ranks_S2 <- sort(ICC_ranks_S2, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S2 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S2,
                           nperm = 10000)

## Scenario 3 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("mQTL_Probes.csv")
mQTL_probes <- sample(mQTL_probes$DNAm.probe, 35000)
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S3 <- data[,c(1,4)]
ICC_ranks_S3 <- ICC_ranks_S3 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S3 <- as.numeric(ICC_ranks_S3$ICC_Scen3)
names(ICC_ranks_S3) <- data$DNAm.probe
ICC_ranks_S3 <- sort(ICC_ranks_S3, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S3 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S3,
                           nperm = 10000)

## Scenario 4 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("mQTL_Probes.csv")
mQTL_probes <- sample(mQTL_probes$DNAm.probe, 35000)
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S4 <- data[,c(1,5)]
ICC_ranks_S4 <- ICC_ranks_S4 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S4 <- as.numeric(ICC_ranks_S4$ICC_Scen4)
names(ICC_ranks_S4) <- data$DNAm.probe
ICC_ranks_S4 <- sort(ICC_ranks_S4, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S4 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S4,
                           nperm = 10000)

## Scenario 5 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("mQTL_Probes.csv")
mQTL_probes <- sample(mQTL_probes$DNAm.probe, 35000)
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S5 <- data[,c(1,6)]
ICC_ranks_S5 <- ICC_ranks_S5 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S5 <- as.numeric(ICC_ranks_S5$ICC_Scen5)
names(ICC_ranks_S5) <- data$DNAm.probe
ICC_ranks_S5 <- sort(ICC_ranks_S5, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S5 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S5,
                           nperm = 10000)

## Scenario 6 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("mQTL_Probes.csv")
mQTL_probes <- sample(mQTL_probes$DNAm.probe, 35000)
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S6 <- data[,c(1,7)]
ICC_ranks_S6 <- ICC_ranks_S6 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S6 <- as.numeric(ICC_ranks_S6$ICC_Scen6)
names(ICC_ranks_S6) <- data$DNAm.probe
ICC_ranks_S6 <- sort(ICC_ranks_S6, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S6 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S6,
                           nperm = 10000)

## Scenario 7 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("mQTL_Probes.csv")
mQTL_probes <- sample(mQTL_probes$DNAm.probe, 35000)
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S7 <- data[,c(1,8)]
ICC_ranks_S7 <- ICC_ranks_S7 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S7 <- as.numeric(ICC_ranks_S7$ICC_Scen7)
names(ICC_ranks_S7) <- data$DNAm.probe
ICC_ranks_S7 <- sort(ICC_ranks_S7, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S7 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S7,
                           nperm = 10000)

## Scenario 8 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("mQTL_Probes.csv")
mQTL_probes <- sample(mQTL_probes$DNAm.probe, 35000)
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S8 <- data[,c(1,9)]
ICC_ranks_S8 <- ICC_ranks_S8 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S8 <- as.numeric(ICC_ranks_S8$ICC_Scen8)
names(ICC_ranks_S8) <- data$DNAm.probe
ICC_ranks_S8 <- sort(ICC_ranks_S8, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S8 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S8,
                           nperm = 10000)

## Scenario 9 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("mQTL_Probes.csv")
mQTL_probes <- sample(mQTL_probes$DNAm.probe, 35000)
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S9 <- data[,c(1,10)]
ICC_ranks_S9 <- ICC_ranks_S9 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S9 <- as.numeric(ICC_ranks_S9$ICC_Scen9)
names(ICC_ranks_S9) <- data$DNAm.probe
ICC_ranks_S9 <- sort(ICC_ranks_S9, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S9 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S9,
                           nperm = 10000)

## Scenario 10 ##
# Read in significant mQTL probes #
mQTL_probes <- read.csv("mQTL_Probes.csv")
mQTL_probes <- sample(mQTL_probes$DNAm.probe, 35000)
mQTL_probes <- list(mQTL = mQTL_probes)
ICC_ranks_S10 <- data[,c(1,11)]
ICC_ranks_S10 <- ICC_ranks_S10 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S10 <- as.numeric(ICC_ranks_S10$ICC_Scen10)
names(ICC_ranks_S10) <- data$DNAm.probe
ICC_ranks_S10 <- sort(ICC_ranks_S10, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S10 <- fgseaSimple(pathways = mQTL_probes, 
                            stats = ICC_ranks_S10,
                            nperm = 10000)

###########
### TSS ###
###########
### Probe Set Enrichment Analysis ### ----
## Scenario 1 ##
# Read in significant mQTL probes #
keep <- grepl("TSS", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
TSS_Probes <- annEPIC@listData$Name[keep]
TSS_Probes <- sample(TSS_Probes, 35000)
mQTL_probes <- list(mQTL = TSS_Probes)
ICC_ranks_S1 <- data[,c(1,2)]
ICC_ranks_S1 <- ICC_ranks_S1 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S1 <- as.numeric(ICC_ranks_S1$ICC_Scen1)
names(ICC_ranks_S1) <- data$DNAm.probe
ICC_ranks_S1 <- sort(ICC_ranks_S1, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S1 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S1,
                           nperm = 10000)

## Scenario 2 ##
# Read in significant mQTL probes #
keep <- grepl("TSS", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
TSS_Probes <- annEPIC@listData$Name[keep]
TSS_Probes <- sample(TSS_Probes, 35000)
mQTL_probes <- list(mQTL = TSS_Probes)
ICC_ranks_S2 <- data[,c(1,3)]
ICC_ranks_S2 <- ICC_ranks_S2 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S2 <- as.numeric(ICC_ranks_S2$ICC_Scen2)
names(ICC_ranks_S2) <- data$DNAm.probe
ICC_ranks_S2 <- sort(ICC_ranks_S2, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S2 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S2,
                           nperm = 10000)

## Scenario 3 ##
# Read in significant mQTL probes #
keep <- grepl("TSS", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
TSS_Probes <- annEPIC@listData$Name[keep]
TSS_Probes <- sample(TSS_Probes, 35000)
mQTL_probes <- list(mQTL = TSS_Probes)
ICC_ranks_S3 <- data[,c(1,4)]
ICC_ranks_S3 <- ICC_ranks_S3 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S3 <- as.numeric(ICC_ranks_S3$ICC_Scen3)
names(ICC_ranks_S3) <- data$DNAm.probe
ICC_ranks_S3 <- sort(ICC_ranks_S3, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S3 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S3,
                           nperm = 10000)

## Scenario 4 ##
# Read in significant mQTL probes #
keep <- grepl("TSS", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
TSS_Probes <- annEPIC@listData$Name[keep]
TSS_Probes <- sample(TSS_Probes, 35000)
mQTL_probes <- list(mQTL = TSS_Probes)
ICC_ranks_S4 <- data[,c(1,5)]
ICC_ranks_S4 <- ICC_ranks_S4 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S4 <- as.numeric(ICC_ranks_S4$ICC_Scen4)
names(ICC_ranks_S4) <- data$DNAm.probe
ICC_ranks_S4 <- sort(ICC_ranks_S4, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S4 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S4,
                           nperm = 10000)

## Scenario 5 ##
# Read in significant mQTL probes #
keep <- grepl("TSS", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
TSS_Probes <- annEPIC@listData$Name[keep]
TSS_Probes <- sample(TSS_Probes, 35000)
mQTL_probes <- list(mQTL = TSS_Probes)
ICC_ranks_S5 <- data[,c(1,6)]
ICC_ranks_S5 <- ICC_ranks_S5 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S5 <- as.numeric(ICC_ranks_S5$ICC_Scen5)
names(ICC_ranks_S5) <- data$DNAm.probe
ICC_ranks_S5 <- sort(ICC_ranks_S5, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S5 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S5,
                           nperm = 10000)

## Scenario 6 ##
# Read in significant mQTL probes #
keep <- grepl("TSS", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
TSS_Probes <- annEPIC@listData$Name[keep]
TSS_Probes <- sample(TSS_Probes, 35000)
mQTL_probes <- list(mQTL = TSS_Probes)
ICC_ranks_S6 <- data[,c(1,7)]
ICC_ranks_S6 <- ICC_ranks_S6 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S6 <- as.numeric(ICC_ranks_S6$ICC_Scen6)
names(ICC_ranks_S6) <- data$DNAm.probe
ICC_ranks_S6 <- sort(ICC_ranks_S6, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S6 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S6,
                           nperm = 10000)

## Scenario 7 ##
# Read in significant mQTL probes #
keep <- grepl("TSS", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
TSS_Probes <- annEPIC@listData$Name[keep]
TSS_Probes <- sample(TSS_Probes, 35000)
mQTL_probes <- list(mQTL = TSS_Probes)
ICC_ranks_S7 <- data[,c(1,8)]
ICC_ranks_S7 <- ICC_ranks_S7 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S7 <- as.numeric(ICC_ranks_S7$ICC_Scen7)
names(ICC_ranks_S7) <- data$DNAm.probe
ICC_ranks_S7 <- sort(ICC_ranks_S7, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S7 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S7,
                           nperm = 10000)

## Scenario 8 ##
# Read in significant mQTL probes #
keep <- grepl("TSS", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
TSS_Probes <- annEPIC@listData$Name[keep]
TSS_Probes <- sample(TSS_Probes, 35000)
mQTL_probes <- list(mQTL = TSS_Probes)
ICC_ranks_S8 <- data[,c(1,9)]
ICC_ranks_S8 <- ICC_ranks_S8 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S8 <- as.numeric(ICC_ranks_S8$ICC_Scen8)
names(ICC_ranks_S8) <- data$DNAm.probe
ICC_ranks_S8 <- sort(ICC_ranks_S8, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S8 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S8,
                           nperm = 10000)

## Scenario 9 ##
# Read in significant mQTL probes #
keep <- grepl("TSS", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
TSS_Probes <- annEPIC@listData$Name[keep]
TSS_Probes <- sample(TSS_Probes, 35000)
mQTL_probes <- list(mQTL = TSS_Probes)
ICC_ranks_S9 <- data[,c(1,10)]
ICC_ranks_S9 <- ICC_ranks_S9 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S9 <- as.numeric(ICC_ranks_S9$ICC_Scen9)
names(ICC_ranks_S9) <- data$DNAm.probe
ICC_ranks_S9 <- sort(ICC_ranks_S9, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S9 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S9,
                           nperm = 10000)

## Scenario 10 ##
# Read in significant mQTL probes #
keep <- grepl("TSS", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
TSS_Probes <- annEPIC@listData$Name[keep]
TSS_Probes <- sample(TSS_Probes, 35000)
mQTL_probes <- list(mQTL = TSS_Probes)
ICC_ranks_S10 <- data[,c(1,11)]
ICC_ranks_S10 <- ICC_ranks_S10 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S10 <- as.numeric(ICC_ranks_S10$ICC_Scen10)
names(ICC_ranks_S10) <- data$DNAm.probe
ICC_ranks_S10 <- sort(ICC_ranks_S10, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S10 <- fgseaSimple(pathways = mQTL_probes, 
                            stats = ICC_ranks_S10,
                            nperm = 10000)

#############
### 5'UTR ###
#############
### Probe Set Enrichment Analysis ### ----
## Scenario 1 ##
# Read in significant mQTL probes #
keep <- grepl("5'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Five_UTR_Probes <- annEPIC@listData$Name[keep]
Five_UTR_Probes <- sample(Five_UTR_Probes, 35000)
mQTL_probes <- list(mQTL = Five_UTR_Probes)
ICC_ranks_S1 <- data[,c(1,2)]
ICC_ranks_S1 <- ICC_ranks_S1 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S1 <- as.numeric(ICC_ranks_S1$ICC_Scen1)
names(ICC_ranks_S1) <- data$DNAm.probe
ICC_ranks_S1 <- sort(ICC_ranks_S1, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S1 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S1,
                           nperm = 10000)

## Scenario 2 ##
# Read in significant mQTL probes #
keep <- grepl("5'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Five_UTR_Probes <- annEPIC@listData$Name[keep]
Five_UTR_Probes <- sample(Five_UTR_Probes, 35000)
mQTL_probes <- list(mQTL = Five_UTR_Probes)
ICC_ranks_S2 <- data[,c(1,3)]
ICC_ranks_S2 <- ICC_ranks_S2 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S2 <- as.numeric(ICC_ranks_S2$ICC_Scen2)
names(ICC_ranks_S2) <- data$DNAm.probe
ICC_ranks_S2 <- sort(ICC_ranks_S2, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S2 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S2,
                           nperm = 10000)

## Scenario 3 ##
# Read in significant mQTL probes #
keep <- grepl("5'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Five_UTR_Probes <- annEPIC@listData$Name[keep]
Five_UTR_Probes <- sample(Five_UTR_Probes, 35000)
mQTL_probes <- list(mQTL = Five_UTR_Probes)
ICC_ranks_S3 <- data[,c(1,4)]
ICC_ranks_S3 <- ICC_ranks_S3 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S3 <- as.numeric(ICC_ranks_S3$ICC_Scen3)
names(ICC_ranks_S3) <- data$DNAm.probe
ICC_ranks_S3 <- sort(ICC_ranks_S3, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S3 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S3,
                           nperm = 10000)

## Scenario 4 ##
# Read in significant mQTL probes #
keep <- grepl("5'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Five_UTR_Probes <- annEPIC@listData$Name[keep]
Five_UTR_Probes <- sample(Five_UTR_Probes, 35000)
mQTL_probes <- list(mQTL = Five_UTR_Probes)
ICC_ranks_S4 <- data[,c(1,5)]
ICC_ranks_S4 <- ICC_ranks_S4 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S4 <- as.numeric(ICC_ranks_S4$ICC_Scen4)
names(ICC_ranks_S4) <- data$DNAm.probe
ICC_ranks_S4 <- sort(ICC_ranks_S4, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S4 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S4,
                           nperm = 10000)

## Scenario 5 ##
# Read in significant mQTL probes #
keep <- grepl("5'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Five_UTR_Probes <- annEPIC@listData$Name[keep]
Five_UTR_Probes <- sample(Five_UTR_Probes, 35000)
mQTL_probes <- list(mQTL = Five_UTR_Probes)
ICC_ranks_S5 <- data[,c(1,6)]
ICC_ranks_S5 <- ICC_ranks_S5 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S5 <- as.numeric(ICC_ranks_S5$ICC_Scen5)
names(ICC_ranks_S5) <- data$DNAm.probe
ICC_ranks_S5 <- sort(ICC_ranks_S5, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S5 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S5,
                           nperm = 10000)

## Scenario 6 ##
# Read in significant mQTL probes #
keep <- grepl("5'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Five_UTR_Probes <- annEPIC@listData$Name[keep]
Five_UTR_Probes <- sample(Five_UTR_Probes, 35000)
mQTL_probes <- list(mQTL = Five_UTR_Probes)
ICC_ranks_S6 <- data[,c(1,7)]
ICC_ranks_S6 <- ICC_ranks_S6 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S6 <- as.numeric(ICC_ranks_S6$ICC_Scen6)
names(ICC_ranks_S6) <- data$DNAm.probe
ICC_ranks_S6 <- sort(ICC_ranks_S6, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S6 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S6,
                           nperm = 10000)

## Scenario 7 ##
# Read in significant mQTL probes #
keep <- grepl("5'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Five_UTR_Probes <- annEPIC@listData$Name[keep]
Five_UTR_Probes <- sample(Five_UTR_Probes, 35000)
mQTL_probes <- list(mQTL = Five_UTR_Probes)
ICC_ranks_S7 <- data[,c(1,8)]
ICC_ranks_S7 <- ICC_ranks_S7 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S7 <- as.numeric(ICC_ranks_S7$ICC_Scen7)
names(ICC_ranks_S7) <- data$DNAm.probe
ICC_ranks_S7 <- sort(ICC_ranks_S7, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S7 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S7,
                           nperm = 10000)

## Scenario 8 ##
# Read in significant mQTL probes #
keep <- grepl("5'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Five_UTR_Probes <- annEPIC@listData$Name[keep]
Five_UTR_Probes <- sample(Five_UTR_Probes, 35000)
mQTL_probes <- list(mQTL = Five_UTR_Probes)
ICC_ranks_S8 <- data[,c(1,9)]
ICC_ranks_S8 <- ICC_ranks_S8 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S8 <- as.numeric(ICC_ranks_S8$ICC_Scen8)
names(ICC_ranks_S8) <- data$DNAm.probe
ICC_ranks_S8 <- sort(ICC_ranks_S8, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S8 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S8,
                           nperm = 10000)

## Scenario 9 ##
# Read in significant mQTL probes #
keep <- grepl("5'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Five_UTR_Probes <- annEPIC@listData$Name[keep]
Five_UTR_Probes <- sample(Five_UTR_Probes, 35000)
mQTL_probes <- list(mQTL = Five_UTR_Probes)
ICC_ranks_S9 <- data[,c(1,10)]
ICC_ranks_S9 <- ICC_ranks_S9 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S9 <- as.numeric(ICC_ranks_S9$ICC_Scen9)
names(ICC_ranks_S9) <- data$DNAm.probe
ICC_ranks_S9 <- sort(ICC_ranks_S9, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S9 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S9,
                           nperm = 10000)

## Scenario 10 ##
# Read in significant mQTL probes #
keep <- grepl("5'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Five_UTR_Probes <- annEPIC@listData$Name[keep]
Five_UTR_Probes <- sample(Five_UTR_Probes, 35000)
mQTL_probes <- list(mQTL = Five_UTR_Probes)
ICC_ranks_S10 <- data[,c(1,11)]
ICC_ranks_S10 <- ICC_ranks_S10 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S10 <- as.numeric(ICC_ranks_S10$ICC_Scen10)
names(ICC_ranks_S10) <- data$DNAm.probe
ICC_ranks_S10 <- sort(ICC_ranks_S10, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S10 <- fgseaSimple(pathways = mQTL_probes, 
                            stats = ICC_ranks_S10,
                            nperm = 10000)

#############
### 3'UTR ###
#############
### Probe Set Enrichment Analysis ### ----
## Scenario 1 ##
# Read in significant mQTL probes #
keep <- grepl("3'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Three_UTR_Probes <- annEPIC@listData$Name[keep]
Three_UTR_Probes <- sample(Three_UTR_Probes, 35000)
mQTL_probes <- list(mQTL = Three_UTR_Probes)
ICC_ranks_S1 <- data[,c(1,2)]
ICC_ranks_S1 <- ICC_ranks_S1 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S1 <- as.numeric(ICC_ranks_S1$ICC_Scen1)
names(ICC_ranks_S1) <- data$DNAm.probe
ICC_ranks_S1 <- sort(ICC_ranks_S1, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S1 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S1,
                           nperm = 10000)

## Scenario 2 ##
# Read in significant mQTL probes #
keep <- grepl("3'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Three_UTR_Probes <- annEPIC@listData$Name[keep]
Three_UTR_Probes <- sample(Three_UTR_Probes, 35000)
mQTL_probes <- list(mQTL = Three_UTR_Probes)
ICC_ranks_S2 <- data[,c(1,3)]
ICC_ranks_S2 <- ICC_ranks_S2 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S2 <- as.numeric(ICC_ranks_S2$ICC_Scen2)
names(ICC_ranks_S2) <- data$DNAm.probe
ICC_ranks_S2 <- sort(ICC_ranks_S2, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S2 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S2,
                           nperm = 10000)

## Scenario 3 ##
# Read in significant mQTL probes #
keep <- grepl("3'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Three_UTR_Probes <- annEPIC@listData$Name[keep]
Three_UTR_Probes <- sample(Three_UTR_Probes, 35000)
mQTL_probes <- list(mQTL = Three_UTR_Probes)
ICC_ranks_S3 <- data[,c(1,4)]
ICC_ranks_S3 <- ICC_ranks_S3 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S3 <- as.numeric(ICC_ranks_S3$ICC_Scen3)
names(ICC_ranks_S3) <- data$DNAm.probe
ICC_ranks_S3 <- sort(ICC_ranks_S3, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S3 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S3,
                           nperm = 10000)

## Scenario 4 ##
# Read in significant mQTL probes #
keep <- grepl("3'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Three_UTR_Probes <- annEPIC@listData$Name[keep]
Three_UTR_Probes <- sample(Three_UTR_Probes, 35000)
mQTL_probes <- list(mQTL = Three_UTR_Probes)
ICC_ranks_S4 <- data[,c(1,5)]
ICC_ranks_S4 <- ICC_ranks_S4 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S4 <- as.numeric(ICC_ranks_S4$ICC_Scen4)
names(ICC_ranks_S4) <- data$DNAm.probe
ICC_ranks_S4 <- sort(ICC_ranks_S4, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S4 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S4,
                           nperm = 10000)

## Scenario 5 ##
# Read in significant mQTL probes #
keep <- grepl("3'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Three_UTR_Probes <- annEPIC@listData$Name[keep]
Three_UTR_Probes <- sample(Three_UTR_Probes, 35000)
mQTL_probes <- list(mQTL = Three_UTR_Probes)
ICC_ranks_S5 <- data[,c(1,6)]
ICC_ranks_S5 <- ICC_ranks_S5 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S5 <- as.numeric(ICC_ranks_S5$ICC_Scen5)
names(ICC_ranks_S5) <- data$DNAm.probe
ICC_ranks_S5 <- sort(ICC_ranks_S5, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S5 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S5,
                           nperm = 10000)

## Scenario 6 ##
# Read in significant mQTL probes #
keep <- grepl("3'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Three_UTR_Probes <- annEPIC@listData$Name[keep]
Three_UTR_Probes <- sample(Three_UTR_Probes, 35000)
mQTL_probes <- list(mQTL = Three_UTR_Probes)
ICC_ranks_S6 <- data[,c(1,7)]
ICC_ranks_S6 <- ICC_ranks_S6 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S6 <- as.numeric(ICC_ranks_S6$ICC_Scen6)
names(ICC_ranks_S6) <- data$DNAm.probe
ICC_ranks_S6 <- sort(ICC_ranks_S6, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S6 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S6,
                           nperm = 10000)

## Scenario 7 ##
# Read in significant mQTL probes #
keep <- grepl("3'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Three_UTR_Probes <- annEPIC@listData$Name[keep]
Three_UTR_Probes <- sample(Three_UTR_Probes, 35000)
mQTL_probes <- list(mQTL = Three_UTR_Probes)
ICC_ranks_S7 <- data[,c(1,8)]
ICC_ranks_S7 <- ICC_ranks_S7 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S7 <- as.numeric(ICC_ranks_S7$ICC_Scen7)
names(ICC_ranks_S7) <- data$DNAm.probe
ICC_ranks_S7 <- sort(ICC_ranks_S7, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S7 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S7,
                           nperm = 10000)

## Scenario 8 ##
# Read in significant mQTL probes #
keep <- grepl("3'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Three_UTR_Probes <- annEPIC@listData$Name[keep]
Three_UTR_Probes <- sample(Three_UTR_Probes, 35000)
mQTL_probes <- list(mQTL = Three_UTR_Probes)
ICC_ranks_S8 <- data[,c(1,9)]
ICC_ranks_S8 <- ICC_ranks_S8 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S8 <- as.numeric(ICC_ranks_S8$ICC_Scen8)
names(ICC_ranks_S8) <- data$DNAm.probe
ICC_ranks_S8 <- sort(ICC_ranks_S8, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S8 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S8,
                           nperm = 10000)

## Scenario 9 ##
# Read in significant mQTL probes #
keep <- grepl("3'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Three_UTR_Probes <- annEPIC@listData$Name[keep]
Three_UTR_Probes <- sample(Three_UTR_Probes, 35000)
mQTL_probes <- list(mQTL = Three_UTR_Probes)
ICC_ranks_S9 <- data[,c(1,10)]
ICC_ranks_S9 <- ICC_ranks_S9 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S9 <- as.numeric(ICC_ranks_S9$ICC_Scen9)
names(ICC_ranks_S9) <- data$DNAm.probe
ICC_ranks_S9 <- sort(ICC_ranks_S9, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S9 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S9,
                           nperm = 10000)

## Scenario 10 ##
# Read in significant mQTL probes #
keep <- grepl("3'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Three_UTR_Probes <- annEPIC@listData$Name[keep]
Three_UTR_Probes <- sample(Three_UTR_Probes, 35000)
mQTL_probes <- list(mQTL = Three_UTR_Probes)
ICC_ranks_S10 <- data[,c(1,11)]
ICC_ranks_S10 <- ICC_ranks_S10 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S10 <- as.numeric(ICC_ranks_S10$ICC_Scen10)
names(ICC_ranks_S10) <- data$DNAm.probe
ICC_ranks_S10 <- sort(ICC_ranks_S10, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S10 <- fgseaSimple(pathways = mQTL_probes, 
                            stats = ICC_ranks_S10,
                            nperm = 10000)

#####################
### Coding Region ###
#####################
### Probe Set Enrichment Analysis ### ----
## Scenario 1 ##
# Read in significant mQTL probes #
keep <- grepl("Exon", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep2 <- grepl("Body", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep <- keep | keep2
CodRegProbes <- annEPIC@listData$Name[keep]
CodRegProbes <- sample(CodRegProbes, 35000)
mQTL_probes <- list(mQTL = CodRegProbes)
ICC_ranks_S1 <- data[,c(1,2)]
ICC_ranks_S1 <- ICC_ranks_S1 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S1 <- as.numeric(ICC_ranks_S1$ICC_Scen1)
names(ICC_ranks_S1) <- data$DNAm.probe
ICC_ranks_S1 <- sort(ICC_ranks_S1, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S1 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S1,
                           nperm = 10000)

## Scenario 2 ##
# Read in significant mQTL probes #
keep <- grepl("Exon", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep2 <- grepl("Body", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep <- keep | keep2
CodRegProbes <- annEPIC@listData$Name[keep]
CodRegProbes <- sample(CodRegProbes, 35000)
mQTL_probes <- list(mQTL = CodRegProbes)
ICC_ranks_S2 <- data[,c(1,3)]
ICC_ranks_S2 <- ICC_ranks_S2 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S2 <- as.numeric(ICC_ranks_S2$ICC_Scen2)
names(ICC_ranks_S2) <- data$DNAm.probe
ICC_ranks_S2 <- sort(ICC_ranks_S2, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S2 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S2,
                           nperm = 10000)

## Scenario 3 ##
# Read in significant mQTL probes #
keep <- grepl("Exon", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep2 <- grepl("Body", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep <- keep | keep2
CodRegProbes <- annEPIC@listData$Name[keep]
CodRegProbes <- sample(CodRegProbes, 35000)
mQTL_probes <- list(mQTL = CodRegProbes)
ICC_ranks_S3 <- data[,c(1,4)]
ICC_ranks_S3 <- ICC_ranks_S3 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S3 <- as.numeric(ICC_ranks_S3$ICC_Scen3)
names(ICC_ranks_S3) <- data$DNAm.probe
ICC_ranks_S3 <- sort(ICC_ranks_S3, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S3 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S3,
                           nperm = 10000)

## Scenario 4 ##
# Read in significant mQTL probes #
keep <- grepl("Exon", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep2 <- grepl("Body", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep <- keep | keep2
CodRegProbes <- annEPIC@listData$Name[keep]
CodRegProbes <- sample(CodRegProbes, 35000)
mQTL_probes <- list(mQTL = CodRegProbes)
ICC_ranks_S4 <- data[,c(1,5)]
ICC_ranks_S4 <- ICC_ranks_S4 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S4 <- as.numeric(ICC_ranks_S4$ICC_Scen4)
names(ICC_ranks_S4) <- data$DNAm.probe
ICC_ranks_S4 <- sort(ICC_ranks_S4, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S4 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S4,
                           nperm = 10000)

## Scenario 5 ##
# Read in significant mQTL probes #
keep <- grepl("Exon", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep2 <- grepl("Body", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep <- keep | keep2
CodRegProbes <- annEPIC@listData$Name[keep]
CodRegProbes <- sample(CodRegProbes, 35000)
mQTL_probes <- list(mQTL = CodRegProbes)
ICC_ranks_S5 <- data[,c(1,6)]
ICC_ranks_S5 <- ICC_ranks_S5 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S5 <- as.numeric(ICC_ranks_S5$ICC_Scen5)
names(ICC_ranks_S5) <- data$DNAm.probe
ICC_ranks_S5 <- sort(ICC_ranks_S5, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S5 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S5,
                           nperm = 10000)

## Scenario 6 ##
# Read in significant mQTL probes #
keep <- grepl("Exon", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep2 <- grepl("Body", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep <- keep | keep2
CodRegProbes <- annEPIC@listData$Name[keep]
CodRegProbes <- sample(CodRegProbes, 35000)
mQTL_probes <- list(mQTL = CodRegProbes)
ICC_ranks_S6 <- data[,c(1,7)]
ICC_ranks_S6 <- ICC_ranks_S6 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S6 <- as.numeric(ICC_ranks_S6$ICC_Scen6)
names(ICC_ranks_S6) <- data$DNAm.probe
ICC_ranks_S6 <- sort(ICC_ranks_S6, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S6 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S6,
                           nperm = 10000)

## Scenario 7 ##
# Read in significant mQTL probes #
keep <- grepl("Exon", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep2 <- grepl("Body", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep <- keep | keep2
CodRegProbes <- annEPIC@listData$Name[keep]
CodRegProbes <- sample(CodRegProbes, 35000)
mQTL_probes <- list(mQTL = CodRegProbes)
ICC_ranks_S7 <- data[,c(1,8)]
ICC_ranks_S7 <- ICC_ranks_S7 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S7 <- as.numeric(ICC_ranks_S7$ICC_Scen7)
names(ICC_ranks_S7) <- data$DNAm.probe
ICC_ranks_S7 <- sort(ICC_ranks_S7, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S7 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S7,
                           nperm = 10000)

## Scenario 8 ##
# Read in significant mQTL probes #
keep <- grepl("Exon", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep2 <- grepl("Body", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep <- keep | keep2
CodRegProbes <- annEPIC@listData$Name[keep]
CodRegProbes <- sample(CodRegProbes, 35000)
mQTL_probes <- list(mQTL = CodRegProbes)
ICC_ranks_S8 <- data[,c(1,9)]
ICC_ranks_S8 <- ICC_ranks_S8 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S8 <- as.numeric(ICC_ranks_S8$ICC_Scen8)
names(ICC_ranks_S8) <- data$DNAm.probe
ICC_ranks_S8 <- sort(ICC_ranks_S8, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S8 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S8,
                           nperm = 10000)

## Scenario 9 ##
# Read in significant mQTL probes #
keep <- grepl("Exon", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep2 <- grepl("Body", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep <- keep | keep2
CodRegProbes <- annEPIC@listData$Name[keep]
CodRegProbes <- sample(CodRegProbes, 35000)
mQTL_probes <- list(mQTL = CodRegProbes)
ICC_ranks_S9 <- data[,c(1,10)]
ICC_ranks_S9 <- ICC_ranks_S9 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S9 <- as.numeric(ICC_ranks_S9$ICC_Scen9)
names(ICC_ranks_S9) <- data$DNAm.probe
ICC_ranks_S9 <- sort(ICC_ranks_S9, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S9 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S9,
                           nperm = 10000)

## Scenario 10 ##
# Read in significant mQTL probes #
keep <- grepl("Exon", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep2 <- grepl("Body", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep <- keep | keep2
CodRegProbes <- annEPIC@listData$Name[keep]
CodRegProbes <- sample(CodRegProbes, 35000)
mQTL_probes <- list(mQTL = CodRegProbes)
ICC_ranks_S10 <- data[,c(1,11)]
ICC_ranks_S10 <- ICC_ranks_S10 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S10 <- as.numeric(ICC_ranks_S10$ICC_Scen10)
names(ICC_ranks_S10) <- data$DNAm.probe
ICC_ranks_S10 <- sort(ICC_ranks_S10, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S10 <- fgseaSimple(pathways = mQTL_probes, 
                            stats = ICC_ranks_S10,
                            nperm = 10000)

#################
### Undefined ###
#################
### Probe Set Enrichment Analysis ### ----
## Scenario 1 ##
# Read in significant mQTL probes #
NoRegProbes <- annEPIC@listData$Name[annEPIC@listData$UCSC_RefGene_Group == ""]
NoRegProbes <- sample(NoRegProbes, 35000)
mQTL_probes <- list(mQTL = NoRegProbes)
ICC_ranks_S1 <- data[,c(1,2)]
ICC_ranks_S1 <- ICC_ranks_S1 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S1 <- as.numeric(ICC_ranks_S1$ICC_Scen1)
names(ICC_ranks_S1) <- data$DNAm.probe
ICC_ranks_S1 <- sort(ICC_ranks_S1, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S1 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S1,
                           nperm = 10000)

## Scenario 2 ##
# Read in significant mQTL probes #
NoRegProbes <- annEPIC@listData$Name[annEPIC@listData$UCSC_RefGene_Group == ""]
NoRegProbes <- sample(NoRegProbes, 35000)
mQTL_probes <- list(mQTL = NoRegProbes)
ICC_ranks_S2 <- data[,c(1,3)]
ICC_ranks_S2 <- ICC_ranks_S2 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S2 <- as.numeric(ICC_ranks_S2$ICC_Scen2)
names(ICC_ranks_S2) <- data$DNAm.probe
ICC_ranks_S2 <- sort(ICC_ranks_S2, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S2 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S2,
                           nperm = 10000)

## Scenario 3 ##
# Read in significant mQTL probes #
NoRegProbes <- annEPIC@listData$Name[annEPIC@listData$UCSC_RefGene_Group == ""]
NoRegProbes <- sample(NoRegProbes, 35000)
mQTL_probes <- list(mQTL = NoRegProbes)
ICC_ranks_S3 <- data[,c(1,4)]
ICC_ranks_S3 <- ICC_ranks_S3 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S3 <- as.numeric(ICC_ranks_S3$ICC_Scen3)
names(ICC_ranks_S3) <- data$DNAm.probe
ICC_ranks_S3 <- sort(ICC_ranks_S3, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S3 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S3,
                           nperm = 10000)

## Scenario 4 ##
# Read in significant mQTL probes #
NoRegProbes <- annEPIC@listData$Name[annEPIC@listData$UCSC_RefGene_Group == ""]
NoRegProbes <- sample(NoRegProbes, 35000)
mQTL_probes <- list(mQTL = NoRegProbes)
ICC_ranks_S4 <- data[,c(1,5)]
ICC_ranks_S4 <- ICC_ranks_S4 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S4 <- as.numeric(ICC_ranks_S4$ICC_Scen4)
names(ICC_ranks_S4) <- data$DNAm.probe
ICC_ranks_S4 <- sort(ICC_ranks_S4, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S4 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S4,
                           nperm = 10000)

## Scenario 5 ##
# Read in significant mQTL probes #
NoRegProbes <- annEPIC@listData$Name[annEPIC@listData$UCSC_RefGene_Group == ""]
NoRegProbes <- sample(NoRegProbes, 35000)
mQTL_probes <- list(mQTL = NoRegProbes)
ICC_ranks_S5 <- data[,c(1,6)]
ICC_ranks_S5 <- ICC_ranks_S5 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S5 <- as.numeric(ICC_ranks_S5$ICC_Scen5)
names(ICC_ranks_S5) <- data$DNAm.probe
ICC_ranks_S5 <- sort(ICC_ranks_S5, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S5 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S5,
                           nperm = 10000)

## Scenario 6 ##
# Read in significant mQTL probes #
NoRegProbes <- annEPIC@listData$Name[annEPIC@listData$UCSC_RefGene_Group == ""]
NoRegProbes <- sample(NoRegProbes, 35000)
mQTL_probes <- list(mQTL = NoRegProbes)
ICC_ranks_S6 <- data[,c(1,7)]
ICC_ranks_S6 <- ICC_ranks_S6 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S6 <- as.numeric(ICC_ranks_S6$ICC_Scen6)
names(ICC_ranks_S6) <- data$DNAm.probe
ICC_ranks_S6 <- sort(ICC_ranks_S6, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S6 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S6,
                           nperm = 10000)

## Scenario 7 ##
# Read in significant mQTL probes #
NoRegProbes <- annEPIC@listData$Name[annEPIC@listData$UCSC_RefGene_Group == ""]
NoRegProbes <- sample(NoRegProbes, 35000)
mQTL_probes <- list(mQTL = NoRegProbes)
ICC_ranks_S7 <- data[,c(1,8)]
ICC_ranks_S7 <- ICC_ranks_S7 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S7 <- as.numeric(ICC_ranks_S7$ICC_Scen7)
names(ICC_ranks_S7) <- data$DNAm.probe
ICC_ranks_S7 <- sort(ICC_ranks_S7, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S7 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S7,
                           nperm = 10000)

## Scenario 8 ##
# Read in significant mQTL probes #
NoRegProbes <- annEPIC@listData$Name[annEPIC@listData$UCSC_RefGene_Group == ""]
NoRegProbes <- sample(NoRegProbes, 35000)
mQTL_probes <- list(mQTL = NoRegProbes)
ICC_ranks_S8 <- data[,c(1,9)]
ICC_ranks_S8 <- ICC_ranks_S8 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S8 <- as.numeric(ICC_ranks_S8$ICC_Scen8)
names(ICC_ranks_S8) <- data$DNAm.probe
ICC_ranks_S8 <- sort(ICC_ranks_S8, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S8 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S8,
                           nperm = 10000)

## Scenario 9 ##
# Read in significant mQTL probes #
NoRegProbes <- annEPIC@listData$Name[annEPIC@listData$UCSC_RefGene_Group == ""]
NoRegProbes <- sample(NoRegProbes, 35000)
mQTL_probes <- list(mQTL = NoRegProbes)
ICC_ranks_S9 <- data[,c(1,10)]
ICC_ranks_S9 <- ICC_ranks_S9 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S9 <- as.numeric(ICC_ranks_S9$ICC_Scen9)
names(ICC_ranks_S9) <- data$DNAm.probe
ICC_ranks_S9 <- sort(ICC_ranks_S9, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S9 <- fgseaSimple(pathways = mQTL_probes, 
                           stats = ICC_ranks_S9,
                           nperm = 10000)

## Scenario 10 ##
# Read in significant mQTL probes #
NoRegProbes <- annEPIC@listData$Name[annEPIC@listData$UCSC_RefGene_Group == ""]
NoRegProbes <- sample(NoRegProbes, 35000)
mQTL_probes <- list(mQTL = NoRegProbes)
ICC_ranks_S10 <- data[,c(1,11)]
ICC_ranks_S10 <- ICC_ranks_S10 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S10 <- as.numeric(ICC_ranks_S10$ICC_Scen10)
names(ICC_ranks_S10) <- data$DNAm.probe
ICC_ranks_S10 <- sort(ICC_ranks_S10, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S10 <- fgseaSimple(pathways = mQTL_probes, 
                            stats = ICC_ranks_S10,
                            nperm = 10000)

### Creating Histogram and Boxplots Combined ### ----
keep <- grepl("TSS", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
TSS_Probes <- annEPIC@listData$Name[keep]
data_TSS <- data[data$DNAm.probe %in% TSS_Probes,]
data_TSS$Category <- "TSS"

keep <- grepl("5'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Five_UTR_Probes <- annEPIC@listData$Name[keep]
data_FiveUTR <- data[data$DNAm.probe %in% Five_UTR_Probes,]
data_FiveUTR$Category <- "3'UTR"

keep <- grepl("3'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Three_UTR_Probes <- annEPIC@listData$Name[keep]
data_ThreeUTR <- data[data$DNAm.probe %in% Three_UTR_Probes,]
data_ThreeUTR$Category <- "5'UTR"

keep <- grepl("Exon", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep2 <- grepl("Body", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep <- keep | keep2
CodRegProbes <- annEPIC@listData$Name[keep]
data_Cod <- data[data$DNAm.probe %in% CodRegProbes,]
data_Cod$Category <- "Coding Region"

NoRegProbes <- annEPIC@listData$Name[annEPIC@listData$UCSC_RefGene_Group == ""]
data_NoCod <- data[data$DNAm.probe %in% NoRegProbes,]
data_NoCod$Category <- "Undefined"

data_regions <- rbind(data_TSS, data_FiveUTR, data_ThreeUTR,
                      data_Cod, data_NoCod)

# Generate Plots #
Hist1 <- ggplot(data, aes(x=ICC_Scen1)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.42),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,120000) +
  labs(title = "Scenario 1", x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())

Box1 <- ggplot(data = data_regions, aes(x=Category, y = ICC_Scen1)) +
  geom_boxplot() + 
  coord_flip() +
  theme_classic() +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        #plot.margin = c(0.25,1,1,1),
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())

Hist2 <- ggplot(data, aes(x=ICC_Scen2)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.43),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,120000) +
  labs(title = "Scenario 2", x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())

Box2 <- ggplot(data = data_regions, aes(x=Category, y = ICC_Scen2)) +
  geom_boxplot() + 
  coord_flip() +
  theme_classic() +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        #plot.margin = c(0.25,1,1,1),
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())

Hist3 <- ggplot(data, aes(x=ICC_Scen3)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.26),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,120000) +
  labs(title = "Scenario 3", x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())

Box3 <- ggplot(data = data_regions, aes(x=Category, y = ICC_Scen3)) +
  geom_boxplot() + 
  coord_flip() +
  theme_classic() +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        #plot.margin = c(0.25,1,1,1),
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())

Hist4 <- ggplot(data, aes(x=ICC_Scen4)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.47),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,120000) +
  labs(title = "Scenario 4", x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())

Box4 <- ggplot(data = data_regions, aes(x=Category, y = ICC_Scen4)) +
  geom_boxplot() + 
  coord_flip() +
  theme_classic() +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        #plot.margin = c(0.25,1,1,1),
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())

Hist5 <- ggplot(data, aes(x=ICC_Scen5)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.42),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,120000) +
  labs(title = "Scenario 5", x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())

Box5 <- ggplot(data = data_regions, aes(x=Category, y = ICC_Scen5)) +
  geom_boxplot() + 
  coord_flip() +
  theme_classic() +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        #plot.margin = c(0.25,1,1,1),
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())

Hist6 <- ggplot(data, aes(x=ICC_Scen6)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.43),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,120000) +
  labs(title = "Scenario 6", x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())

Box6 <- ggplot(data = data_regions, aes(x=Category, y = ICC_Scen6)) +
  geom_boxplot() + 
  coord_flip() +
  theme_classic() +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        #plot.margin = c(0.25,1,1,1),
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())

Hist7 <- ggplot(data, aes(x=ICC_Scen7)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.45),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,120000) +
  labs(title = "Scenario 7", x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())

Box7 <- ggplot(data = data_regions, aes(x=Category, y = ICC_Scen7)) +
  geom_boxplot() + 
  coord_flip() +
  theme_classic() +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        #plot.margin = c(0.25,1,1,1),
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())

Hist8 <- ggplot(data, aes(x=ICC_Scen8)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.47),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,120000) +
  labs(title = "Scenario 8", x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())

Box8 <- ggplot(data = data_regions, aes(x=Category, y = ICC_Scen8)) +
  geom_boxplot() + 
  coord_flip() +
  theme_classic() +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        #plot.margin = c(0.25,1,1,1),
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())

Hist9 <- ggplot(data, aes(x=ICC_Scen9)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.43),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,120000) +
  labs(title = "Scenario 9", x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())

Box9 <- ggplot(data = data_regions, aes(x=Category, y = ICC_Scen9)) +
  geom_boxplot() + 
  coord_flip() +
  theme_classic() +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        #plot.margin = c(0.25,1,1,1),
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())

Hist10 <- ggplot(data, aes(x=ICC_Scen10)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=0.25),
             color="black", linetype="dashed", size=1) +
  theme_bw() +
  ylim(0,120000) +
  labs(title = "Scenario 10", x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())

Box10 <- ggplot(data = data_regions, aes(x=Category, y = ICC_Scen10)) +
  geom_boxplot() + 
  coord_flip() +
  theme_classic() +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        #plot.margin = c(0.25,1,1,1),
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())

### Combine All Plots ### ----
ggarrange(Hist1, Hist2, Hist3, Hist4, Hist5, 
          Box1, Box2, Box3, Box4, Box5,
          Hist6, Hist7, Hist8, Hist9, Hist10,
          Box6, Box7, Box8, Box9, Box10,
          ncol = 5, nrow = 4, align = "hv", heights = c(1, 0.25, 1, 0.25))
