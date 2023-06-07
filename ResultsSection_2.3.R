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
library(FlowSorted.Blood.EPIC)
library(FlowSorted.BloodExtended.EPIC)



#############################################################################
######################## Epigenetic Algorithm Probes ########################
#############################################################################
### Set Seed ### ----
set.seed(100)

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
# Exclude probes with ICC < 0.01 #
data[,2:11][data[,2:11] < 0.01] <- NA

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

######################
### Horvath Probes ###
######################
### Probe Set Enrichment Analysis ### ----
## Scenario 1 ##
# Read in probes #
probe_set <- read.csv("Horvath_Probes.csv")
keep <- probe_set$DNAm.probe %in% data$DNAm.probe
probe_set <- probe_set[keep,]
probe_set <- list(probes = probe_set)
ICC_ranks_S1 <- data[,c(1,2)]
ICC_ranks_S1 <- ICC_ranks_S1 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S1 <- as.numeric(ICC_ranks_S1$ICC_Scen1)
names(ICC_ranks_S1) <- data$DNAm.probe
ICC_ranks_S1 <- sort(ICC_ranks_S1, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S1 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S1,
                           nperm = 10000)

## Scenario 2 ##
# Read in probes #
probe_set <- read.csv("Horvath_Probes.csv")
keep <- probe_set$DNAm.probe %in% data$DNAm.probe
probe_set <- probe_set[keep,]
probe_set <- list(probes = probe_set)
ICC_ranks_S2 <- data[,c(1,3)]
ICC_ranks_S2 <- ICC_ranks_S2 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S2 <- as.numeric(ICC_ranks_S2$ICC_Scen2)
names(ICC_ranks_S2) <- data$DNAm.probe
ICC_ranks_S2 <- sort(ICC_ranks_S2, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S2 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S2,
                           nperm = 10000)

## Scenario 3 ##
# Read in probes #
probe_set <- read.csv("Horvath_Probes.csv")
keep <- probe_set$DNAm.probe %in% data$DNAm.probe
probe_set <- probe_set[keep,]
probe_set <- list(probes = probe_set)
ICC_ranks_S3 <- data[,c(1,4)]
ICC_ranks_S3 <- ICC_ranks_S3 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S3 <- as.numeric(ICC_ranks_S3$ICC_Scen3)
names(ICC_ranks_S3) <- data$DNAm.probe
ICC_ranks_S3 <- sort(ICC_ranks_S3, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S3 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S3,
                           nperm = 10000)

## Scenario 4 ##
# Read in probes #
probe_set <- read.csv("Horvath_Probes.csv")
keep <- probe_set$DNAm.probe %in% data$DNAm.probe
probe_set <- probe_set[keep,]
probe_set <- list(probes = probe_set)
ICC_ranks_S4 <- data[,c(1,5)]
ICC_ranks_S4 <- ICC_ranks_S4 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S4 <- as.numeric(ICC_ranks_S4$ICC_Scen4)
names(ICC_ranks_S4) <- data$DNAm.probe
ICC_ranks_S4 <- sort(ICC_ranks_S4, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S4 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S4,
                           nperm = 10000)

## Scenario 5 ##
# Read in probes #
probe_set <- read.csv("Horvath_Probes.csv")
keep <- probe_set$DNAm.probe %in% data$DNAm.probe
probe_set <- probe_set[keep,]
probe_set <- list(probes = probe_set)
ICC_ranks_S5 <- data[,c(1,6)]
ICC_ranks_S5 <- ICC_ranks_S5 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S5 <- as.numeric(ICC_ranks_S5$ICC_Scen5)
names(ICC_ranks_S5) <- data$DNAm.probe
ICC_ranks_S5 <- sort(ICC_ranks_S5, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S5 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S5,
                           nperm = 10000)

## Scenario 6 ##
# Read in probes #
probe_set <- read.csv("Horvath_Probes.csv")
keep <- probe_set$DNAm.probe %in% data$DNAm.probe
probe_set <- probe_set[keep,]
probe_set <- list(probes = probe_set)
ICC_ranks_S6 <- data[,c(1,7)]
ICC_ranks_S6 <- ICC_ranks_S6 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S6 <- as.numeric(ICC_ranks_S6$ICC_Scen6)
names(ICC_ranks_S6) <- data$DNAm.probe
ICC_ranks_S6 <- sort(ICC_ranks_S6, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S6 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S6,
                           nperm = 10000)

## Scenario 7 ##
# Read in probes #
probe_set <- read.csv("Horvath_Probes.csv")
keep <- probe_set$DNAm.probe %in% data$DNAm.probe
probe_set <- probe_set[keep,]
probe_set <- list(probes = probe_set)
ICC_ranks_S7 <- data[,c(1,8)]
ICC_ranks_S7 <- ICC_ranks_S7 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S7 <- as.numeric(ICC_ranks_S7$ICC_Scen7)
names(ICC_ranks_S7) <- data$DNAm.probe
ICC_ranks_S7 <- sort(ICC_ranks_S7, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S7 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S7,
                           nperm = 10000)

## Scenario 8 ##
# Read in probes #
probe_set <- read.csv("Horvath_Probes.csv")
keep <- probe_set$DNAm.probe %in% data$DNAm.probe
probe_set <- probe_set[keep,]
probe_set <- list(probes = probe_set)
ICC_ranks_S8 <- data[,c(1,9)]
ICC_ranks_S8 <- ICC_ranks_S8 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S8 <- as.numeric(ICC_ranks_S8$ICC_Scen8)
names(ICC_ranks_S8) <- data$DNAm.probe
ICC_ranks_S8 <- sort(ICC_ranks_S8, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S8 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S8,
                           nperm = 10000)

## Scenario 9 ##
# Read in probes #
probe_set <- read.csv("Horvath_Probes.csv")
keep <- probe_set$DNAm.probe %in% data$DNAm.probe
probe_set <- probe_set[keep,]
probe_set <- list(probes = probe_set)
ICC_ranks_S9 <- data[,c(1,10)]
ICC_ranks_S9 <- ICC_ranks_S9 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S9 <- as.numeric(ICC_ranks_S9$ICC_Scen9)
names(ICC_ranks_S9) <- data$DNAm.probe
ICC_ranks_S9 <- sort(ICC_ranks_S9, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S9 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S9,
                           nperm = 10000)

## Scenario 10 ##
# Read in probes #
probe_set <- read.csv("Horvath_Probes.csv")
keep <- probe_set$DNAm.probe %in% data$DNAm.probe
probe_set <- probe_set[keep,]
probe_set <- list(probes = probe_set)
ICC_ranks_S10 <- data[,c(1,11)]
ICC_ranks_S10 <- ICC_ranks_S10 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S10 <- as.numeric(ICC_ranks_S10$ICC_Scen10)
names(ICC_ranks_S10) <- data$DNAm.probe
ICC_ranks_S10 <- sort(ICC_ranks_S10, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S10 <- fgseaSimple(pathways = probe_set, 
                            stats = ICC_ranks_S10,
                            nperm = 10000)


#####################
### Hannum Probes ###
#####################
### Probe Set Enrichment Analysis ### ----
## Scenario 1 ##
# Read in probes #
probe_set <- read.csv("Hannum_Probes.csv")
keep <- probe_set$DNAm.probe %in% data$DNAm.probe
probe_set <- probe_set[keep,]
probe_set <- list(probes = probe_set)
ICC_ranks_S1 <- data[,c(1,2)]
ICC_ranks_S1 <- ICC_ranks_S1 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S1 <- as.numeric(ICC_ranks_S1$ICC_Scen1)
names(ICC_ranks_S1) <- data$DNAm.probe
ICC_ranks_S1 <- sort(ICC_ranks_S1, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S1 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S1,
                           nperm = 10000)

## Scenario 2 ##
# Read in probes #
probe_set <- read.csv("Hannum_Probes.csv")
keep <- probe_set$DNAm.probe %in% data$DNAm.probe
probe_set <- probe_set[keep,]
probe_set <- list(probes = probe_set)
ICC_ranks_S2 <- data[,c(1,3)]
ICC_ranks_S2 <- ICC_ranks_S2 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S2 <- as.numeric(ICC_ranks_S2$ICC_Scen2)
names(ICC_ranks_S2) <- data$DNAm.probe
ICC_ranks_S2 <- sort(ICC_ranks_S2, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S2 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S2,
                           nperm = 10000)

## Scenario 3 ##
# Read in probes #
probe_set <- read.csv("Hannum_Probes.csv")
keep <- probe_set$DNAm.probe %in% data$DNAm.probe
probe_set <- probe_set[keep,]
probe_set <- list(probes = probe_set)
ICC_ranks_S3 <- data[,c(1,4)]
ICC_ranks_S3 <- ICC_ranks_S3 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S3 <- as.numeric(ICC_ranks_S3$ICC_Scen3)
names(ICC_ranks_S3) <- data$DNAm.probe
ICC_ranks_S3 <- sort(ICC_ranks_S3, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S3 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S3,
                           nperm = 10000)

## Scenario 4 ##
# Read in probes #
probe_set <- read.csv("Hannum_Probes.csv")
keep <- probe_set$DNAm.probe %in% data$DNAm.probe
probe_set <- probe_set[keep,]
probe_set <- list(probes = probe_set)
ICC_ranks_S4 <- data[,c(1,5)]
ICC_ranks_S4 <- ICC_ranks_S4 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S4 <- as.numeric(ICC_ranks_S4$ICC_Scen4)
names(ICC_ranks_S4) <- data$DNAm.probe
ICC_ranks_S4 <- sort(ICC_ranks_S4, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S4 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S4,
                           nperm = 10000)

## Scenario 5 ##
# Read in probes #
probe_set <- read.csv("Hannum_Probes.csv")
keep <- probe_set$DNAm.probe %in% data$DNAm.probe
probe_set <- probe_set[keep,]
probe_set <- list(probes = probe_set)
ICC_ranks_S5 <- data[,c(1,6)]
ICC_ranks_S5 <- ICC_ranks_S5 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S5 <- as.numeric(ICC_ranks_S5$ICC_Scen5)
names(ICC_ranks_S5) <- data$DNAm.probe
ICC_ranks_S5 <- sort(ICC_ranks_S5, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S5 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S5,
                           nperm = 10000)

## Scenario 6 ##
# Read in probes #
probe_set <- read.csv("Hannum_Probes.csv")
keep <- probe_set$DNAm.probe %in% data$DNAm.probe
probe_set <- probe_set[keep,]
probe_set <- list(probes = probe_set)
ICC_ranks_S6 <- data[,c(1,7)]
ICC_ranks_S6 <- ICC_ranks_S6 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S6 <- as.numeric(ICC_ranks_S6$ICC_Scen6)
names(ICC_ranks_S6) <- data$DNAm.probe
ICC_ranks_S6 <- sort(ICC_ranks_S6, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S6 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S6,
                           nperm = 10000)

## Scenario 7 ##
# Read in probes #
probe_set <- read.csv("Hannum_Probes.csv")
keep <- probe_set$DNAm.probe %in% data$DNAm.probe
probe_set <- probe_set[keep,]
probe_set <- list(probes = probe_set)
ICC_ranks_S7 <- data[,c(1,8)]
ICC_ranks_S7 <- ICC_ranks_S7 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S7 <- as.numeric(ICC_ranks_S7$ICC_Scen7)
names(ICC_ranks_S7) <- data$DNAm.probe
ICC_ranks_S7 <- sort(ICC_ranks_S7, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S7 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S7,
                           nperm = 10000)

## Scenario 8 ##
# Read in probes #
probe_set <- read.csv("Hannum_Probes.csv")
keep <- probe_set$DNAm.probe %in% data$DNAm.probe
probe_set <- probe_set[keep,]
probe_set <- list(probes = probe_set)
ICC_ranks_S8 <- data[,c(1,9)]
ICC_ranks_S8 <- ICC_ranks_S8 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S8 <- as.numeric(ICC_ranks_S8$ICC_Scen8)
names(ICC_ranks_S8) <- data$DNAm.probe
ICC_ranks_S8 <- sort(ICC_ranks_S8, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S8 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S8,
                           nperm = 10000)

## Scenario 9 ##
# Read in probes #
probe_set <- read.csv("Hannum_Probes.csv")
keep <- probe_set$DNAm.probe %in% data$DNAm.probe
probe_set <- probe_set[keep,]
probe_set <- list(probes = probe_set)
ICC_ranks_S9 <- data[,c(1,10)]
ICC_ranks_S9 <- ICC_ranks_S9 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S9 <- as.numeric(ICC_ranks_S9$ICC_Scen9)
names(ICC_ranks_S9) <- data$DNAm.probe
ICC_ranks_S9 <- sort(ICC_ranks_S9, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S9 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S9,
                           nperm = 10000)

## Scenario 10 ##
# Read in probes #
probe_set <- read.csv("Hannum_Probes.csv")
keep <- probe_set$DNAm.probe %in% data$DNAm.probe
probe_set <- probe_set[keep,]
probe_set <- list(probes = probe_set)
ICC_ranks_S10 <- data[,c(1,11)]
ICC_ranks_S10 <- ICC_ranks_S10 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S10 <- as.numeric(ICC_ranks_S10$ICC_Scen10)
names(ICC_ranks_S10) <- data$DNAm.probe
ICC_ranks_S10 <- sort(ICC_ranks_S10, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S10 <- fgseaSimple(pathways = probe_set, 
                            stats = ICC_ranks_S10,
                            nperm = 10000)

#######################
### PhenoAge Probes ###
#######################
### Probe Set Enrichment Analysis ### ----
## Scenario 1 ##
# Read in probes #
PhenoAge_probes <- read.csv("PhenoAge_Probes.csv")
keep <- PhenoAge_probes$DNAm.probe %in% data$DNAm.probe
PhenoAge_probes <- PhenoAge_probes[keep,]
PhenoAge_probes <- list(probes = PhenoAge_probes)
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
# Read in probes #
PhenoAge_probes <- read.csv("PhenoAge_Probes.csv")
keep <- PhenoAge_probes$DNAm.probe %in% PhenoAge_probes$DNAm.probe
PhenoAge_probes <- PhenoAge_probes[keep,]
PhenoAge_probes <- list(probes = PhenoAge_probes)
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
# Read in probes #
PhenoAge_probes <- read.csv("PhenoAge_Probes.csv")
keep <- PhenoAge_probes$DNAm.probe %in% PhenoAge_probes$DNAm.probe
PhenoAge_probes <- PhenoAge_probes[keep,]
PhenoAge_probes <- list(probes = PhenoAge_probes)
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
# Read in probes #
PhenoAge_probes <- read.csv("PhenoAge_Probes.csv")
keep <- PhenoAge_probes$DNAm.probe %in% PhenoAge_probes$DNAm.probe
PhenoAge_probes <- PhenoAge_probes[keep,]
PhenoAge_probes <- list(probes = PhenoAge_probes)
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
# Read in probes #
PhenoAge_probes <- read.csv("PhenoAge_Probes.csv")
keep <- PhenoAge_probes$DNAm.probe %in% PhenoAge_probes$DNAm.probe
PhenoAge_probes <- PhenoAge_probes[keep,]
PhenoAge_probes <- list(probes = PhenoAge_probes)
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
# Read in probes #
PhenoAge_probes <- read.csv("PhenoAge_Probes.csv")
keep <- PhenoAge_probes$DNAm.probe %in% PhenoAge_probes$DNAm.probe
PhenoAge_probes <- PhenoAge_probes[keep,]
PhenoAge_probes <- list(probes = PhenoAge_probes)
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
# Read in probes #
PhenoAge_probes <- read.csv("PhenoAge_Probes.csv")
keep <- PhenoAge_probes$DNAm.probe %in% PhenoAge_probes$DNAm.probe
PhenoAge_probes <- PhenoAge_probes[keep,]
PhenoAge_probes <- list(probes = PhenoAge_probes)
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
# Read in probes #
PhenoAge_probes <- read.csv("PhenoAge_Probes.csv")
keep <- PhenoAge_probes$DNAm.probe %in% PhenoAge_probes$DNAm.probe
PhenoAge_probes <- PhenoAge_probes[keep,]
PhenoAge_probes <- list(probes = PhenoAge_probes)
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
# Read in probes #
PhenoAge_probes <- read.csv("PhenoAge_Probes.csv")
keep <- PhenoAge_probes$DNAm.probe %in% PhenoAge_probes$DNAm.probe
PhenoAge_probes <- PhenoAge_probes[keep,]
PhenoAge_probes <- list(probes = PhenoAge_probes)
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
# Read in probes #
PhenoAge_probes <- read.csv("PhenoAge_Probes.csv")
keep <- PhenoAge_probes$DNAm.probe %in% PhenoAge_probes$DNAm.probe
PhenoAge_probes <- PhenoAge_probes[keep,]
PhenoAge_probes <- list(probes = PhenoAge_probes)
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
# Read in probes #
probe_set <- read.csv("PC_Probes.csv")
probe_set <- as.data.frame(sample(probe_set$DNAm.probe, 35000))
colnames(probe_set) <- "DNAm.probe"
keep <- probe_set$DNAm.probe %in% data$DNAm.probe
probe_set <- probe_set[keep,]
probe_set <- list(probes = probe_set)
ICC_ranks_S1 <- data[,c(1,2)]
ICC_ranks_S1 <- ICC_ranks_S1 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S1 <- as.numeric(ICC_ranks_S1$ICC_Scen1)
names(ICC_ranks_S1) <- data$DNAm.probe
ICC_ranks_S1 <- sort(ICC_ranks_S1, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S1 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S1,
                           nperm = 10000)

## Scenario 2 ##
# Read in probes #
probe_set <- read.csv("PC_Probes.csv")
probe_set <- as.data.frame(sample(probe_set$DNAm.probe, 35000))
colnames(probe_set) <- "DNAm.probe"
keep <- probe_set$DNAm.probe %in% data$DNAm.probe
probe_set <- probe_set[keep,]
probe_set <- list(probes = probe_set)
ICC_ranks_S2 <- data[,c(1,3)]
ICC_ranks_S2 <- ICC_ranks_S2 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S2 <- as.numeric(ICC_ranks_S2$ICC_Scen2)
names(ICC_ranks_S2) <- data$DNAm.probe
ICC_ranks_S2 <- sort(ICC_ranks_S2, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S2 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S2,
                           nperm = 10000)

## Scenario 3 ##
# Read in probes #
probe_set <- read.csv("PC_Probes.csv")
probe_set <- as.data.frame(sample(probe_set$DNAm.probe, 35000))
colnames(probe_set) <- "DNAm.probe"
keep <- probe_set$DNAm.probe %in% data$DNAm.probe
probe_set <- probe_set[keep,]
probe_set <- list(probes = probe_set)
ICC_ranks_S3 <- data[,c(1,4)]
ICC_ranks_S3 <- ICC_ranks_S3 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S3 <- as.numeric(ICC_ranks_S3$ICC_Scen3)
names(ICC_ranks_S3) <- data$DNAm.probe
ICC_ranks_S3 <- sort(ICC_ranks_S3, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S3 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S3,
                           nperm = 10000)

## Scenario 4 ##
# Read in probes #
probe_set <- read.csv("PC_Probes.csv")
probe_set <- as.data.frame(sample(probe_set$DNAm.probe, 35000))
colnames(probe_set) <- "DNAm.probe"
keep <- probe_set$DNAm.probe %in% data$DNAm.probe
probe_set <- probe_set[keep,]
probe_set <- list(probes = probe_set)
ICC_ranks_S4 <- data[,c(1,5)]
ICC_ranks_S4 <- ICC_ranks_S4 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S4 <- as.numeric(ICC_ranks_S4$ICC_Scen4)
names(ICC_ranks_S4) <- data$DNAm.probe
ICC_ranks_S4 <- sort(ICC_ranks_S4, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S4 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S4,
                           nperm = 10000)

## Scenario 5 ##
# Read in probes #
probe_set <- read.csv("PC_Probes.csv")
probe_set <- as.data.frame(sample(probe_set$DNAm.probe, 35000))
colnames(probe_set) <- "DNAm.probe"
keep <- probe_set$DNAm.probe %in% data$DNAm.probe
probe_set <- probe_set[keep,]
probe_set <- list(probes = probe_set)
ICC_ranks_S5 <- data[,c(1,6)]
ICC_ranks_S5 <- ICC_ranks_S5 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S5 <- as.numeric(ICC_ranks_S5$ICC_Scen5)
names(ICC_ranks_S5) <- data$DNAm.probe
ICC_ranks_S5 <- sort(ICC_ranks_S5, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S5 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S5,
                           nperm = 10000)

## Scenario 6 ##
# Read in probes #
probe_set <- read.csv("PC_Probes.csv")
probe_set <- as.data.frame(sample(probe_set$DNAm.probe, 35000))
colnames(probe_set) <- "DNAm.probe"
keep <- probe_set$DNAm.probe %in% data$DNAm.probe
probe_set <- probe_set[keep,]
probe_set <- list(probes = probe_set)
ICC_ranks_S6 <- data[,c(1,7)]
ICC_ranks_S6 <- ICC_ranks_S6 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S6 <- as.numeric(ICC_ranks_S6$ICC_Scen6)
names(ICC_ranks_S6) <- data$DNAm.probe
ICC_ranks_S6 <- sort(ICC_ranks_S6, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S6 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S6,
                           nperm = 10000)

## Scenario 7 ##
# Read in probes #
probe_set <- read.csv("PC_Probes.csv")
probe_set <- as.data.frame(sample(probe_set$DNAm.probe, 35000))
colnames(probe_set) <- "DNAm.probe"
keep <- probe_set$DNAm.probe %in% data$DNAm.probe
probe_set <- probe_set[keep,]
probe_set <- list(probes = probe_set)
ICC_ranks_S7 <- data[,c(1,8)]
ICC_ranks_S7 <- ICC_ranks_S7 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S7 <- as.numeric(ICC_ranks_S7$ICC_Scen7)
names(ICC_ranks_S7) <- data$DNAm.probe
ICC_ranks_S7 <- sort(ICC_ranks_S7, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S7 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S7,
                           nperm = 10000)

## Scenario 8 ##
# Read in probes #
probe_set <- read.csv("PC_Probes.csv")
probe_set <- as.data.frame(sample(probe_set$DNAm.probe, 35000))
colnames(probe_set) <- "DNAm.probe"
keep <- probe_set$DNAm.probe %in% data$DNAm.probe
probe_set <- probe_set[keep,]
probe_set <- list(probes = probe_set)
ICC_ranks_S8 <- data[,c(1,9)]
ICC_ranks_S8 <- ICC_ranks_S8 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S8 <- as.numeric(ICC_ranks_S8$ICC_Scen8)
names(ICC_ranks_S8) <- data$DNAm.probe
ICC_ranks_S8 <- sort(ICC_ranks_S8, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S8 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S8,
                           nperm = 10000)

## Scenario 9 ##
# Read in probes #
probe_set <- read.csv("PC_Probes.csv")
probe_set <- as.data.frame(sample(probe_set$DNAm.probe, 35000))
colnames(probe_set) <- "DNAm.probe"
keep <- probe_set$DNAm.probe %in% data$DNAm.probe
probe_set <- probe_set[keep,]
probe_set <- list(probes = probe_set)
ICC_ranks_S9 <- data[,c(1,10)]
ICC_ranks_S9 <- ICC_ranks_S9 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S9 <- as.numeric(ICC_ranks_S9$ICC_Scen9)
names(ICC_ranks_S9) <- data$DNAm.probe
ICC_ranks_S9 <- sort(ICC_ranks_S9, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S9 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S9,
                           nperm = 10000)

## Scenario 10 ##
# Read in probes #
probe_set <- read.csv("PC_Probes.csv")
probe_set <- as.data.frame(sample(probe_set$DNAm.probe, 35000))
colnames(probe_set) <- "DNAm.probe"
keep <- probe_set$DNAm.probe %in% data$DNAm.probe
probe_set <- probe_set[keep,]
probe_set <- list(probes = probe_set)
ICC_ranks_S10 <- data[,c(1,11)]
ICC_ranks_S10 <- ICC_ranks_S10 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S10 <- as.numeric(ICC_ranks_S10$ICC_Scen10)
names(ICC_ranks_S10) <- data$DNAm.probe
ICC_ranks_S10 <- sort(ICC_ranks_S10, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S10 <- fgseaSimple(pathways = probe_set, 
                            stats = ICC_ranks_S10,
                            nperm = 10000)

###################
### PACE Probes ###
###################
### Probe Set Enrichment Analysis ### ----
## Scenario 1 ##
# Read in probes #
load("sysdata.rda")
keep <- mPOA_Models$model_probes$DunedinPACE %in% data$DNAm.probe
probe_set <- mPOA_Models$model_probes$DunedinPACE[keep]
probe_set <- list(probes = probe_set)
ICC_ranks_S1 <- data[,c(1,2)]
ICC_ranks_S1 <- ICC_ranks_S1 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S1 <- as.numeric(ICC_ranks_S1$ICC_Scen1)
names(ICC_ranks_S1) <- data$DNAm.probe
ICC_ranks_S1 <- sort(ICC_ranks_S1, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S1 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S1,
                           nperm = 10000)

## Scenario 2 ##
# Read in probes #
load("sysdata.rda")
keep <- mPOA_Models$model_probes$DunedinPACE %in% data$DNAm.probe
probe_set <- mPOA_Models$model_probes$DunedinPACE[keep]
probe_set <- list(probes = probe_set)
ICC_ranks_S2 <- data[,c(1,3)]
ICC_ranks_S2 <- ICC_ranks_S2 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S2 <- as.numeric(ICC_ranks_S2$ICC_Scen2)
names(ICC_ranks_S2) <- data$DNAm.probe
ICC_ranks_S2 <- sort(ICC_ranks_S2, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S2 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S2,
                           nperm = 10000)

## Scenario 3 ##
# Read in probes #
load("sysdata.rda")
keep <- mPOA_Models$model_probes$DunedinPACE %in% data$DNAm.probe
probe_set <- mPOA_Models$model_probes$DunedinPACE[keep]
probe_set <- list(probes = probe_set)
ICC_ranks_S3 <- data[,c(1,4)]
ICC_ranks_S3 <- ICC_ranks_S3 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S3 <- as.numeric(ICC_ranks_S3$ICC_Scen3)
names(ICC_ranks_S3) <- data$DNAm.probe
ICC_ranks_S3 <- sort(ICC_ranks_S3, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S3 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S3,
                           nperm = 10000)

## Scenario 4 ##
# Read in probes #
load("sysdata.rda")
keep <- mPOA_Models$model_probes$DunedinPACE %in% data$DNAm.probe
probe_set <- mPOA_Models$model_probes$DunedinPACE[keep]
probe_set <- list(probes = probe_set)
ICC_ranks_S4 <- data[,c(1,5)]
ICC_ranks_S4 <- ICC_ranks_S4 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S4 <- as.numeric(ICC_ranks_S4$ICC_Scen4)
names(ICC_ranks_S4) <- data$DNAm.probe
ICC_ranks_S4 <- sort(ICC_ranks_S4, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S4 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S4,
                           nperm = 10000)

## Scenario 5 ##
# Read in probes #
load("sysdata.rda")
keep <- mPOA_Models$model_probes$DunedinPACE %in% data$DNAm.probe
probe_set <- mPOA_Models$model_probes$DunedinPACE[keep]
probe_set <- list(probes = probe_set)
ICC_ranks_S5 <- data[,c(1,6)]
ICC_ranks_S5 <- ICC_ranks_S5 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S5 <- as.numeric(ICC_ranks_S5$ICC_Scen5)
names(ICC_ranks_S5) <- data$DNAm.probe
ICC_ranks_S5 <- sort(ICC_ranks_S5, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S5 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S5,
                           nperm = 10000)

## Scenario 6 ##
# Read in probes #
load("sysdata.rda")
keep <- mPOA_Models$model_probes$DunedinPACE %in% data$DNAm.probe
probe_set <- mPOA_Models$model_probes$DunedinPACE[keep]
probe_set <- list(probes = probe_set)
ICC_ranks_S6 <- data[,c(1,7)]
ICC_ranks_S6 <- ICC_ranks_S6 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S6 <- as.numeric(ICC_ranks_S6$ICC_Scen6)
names(ICC_ranks_S6) <- data$DNAm.probe
ICC_ranks_S6 <- sort(ICC_ranks_S6, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S6 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S6,
                           nperm = 10000)

## Scenario 7 ##
# Read in probes #
load("sysdata.rda")
keep <- mPOA_Models$model_probes$DunedinPACE %in% data$DNAm.probe
probe_set <- mPOA_Models$model_probes$DunedinPACE[keep]
probe_set <- list(probes = probe_set)
ICC_ranks_S7 <- data[,c(1,8)]
ICC_ranks_S7 <- ICC_ranks_S7 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S7 <- as.numeric(ICC_ranks_S7$ICC_Scen7)
names(ICC_ranks_S7) <- data$DNAm.probe
ICC_ranks_S7 <- sort(ICC_ranks_S7, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S7 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S7,
                           nperm = 10000)

## Scenario 8 ##
# Read in probes #
load("sysdata.rda")
keep <- mPOA_Models$model_probes$DunedinPACE %in% data$DNAm.probe
probe_set <- mPOA_Models$model_probes$DunedinPACE[keep]
probe_set <- list(probes = probe_set)
ICC_ranks_S8 <- data[,c(1,9)]
ICC_ranks_S8 <- ICC_ranks_S8 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S8 <- as.numeric(ICC_ranks_S8$ICC_Scen8)
names(ICC_ranks_S8) <- data$DNAm.probe
ICC_ranks_S8 <- sort(ICC_ranks_S8, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S8 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S8,
                           nperm = 10000)

## Scenario 9 ##
# Read in probes #
load("sysdata.rda")
keep <- mPOA_Models$model_probes$DunedinPACE %in% data$DNAm.probe
probe_set <- mPOA_Models$model_probes$DunedinPACE[keep]
probe_set <- list(probes = probe_set)
ICC_ranks_S9 <- data[,c(1,10)]
ICC_ranks_S9 <- ICC_ranks_S9 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S9 <- as.numeric(ICC_ranks_S9$ICC_Scen9)
names(ICC_ranks_S9) <- data$DNAm.probe
ICC_ranks_S9 <- sort(ICC_ranks_S9, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S9 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S9,
                           nperm = 10000)

## Scenario 10 ##
# Read in probes #
load("sysdata.rda")
keep <- mPOA_Models$model_probes$DunedinPACE %in% data$DNAm.probe
probe_set <- mPOA_Models$model_probes$DunedinPACE[keep]
probe_set <- list(probes = probe_set)
ICC_ranks_S10 <- data[,c(1,11)]
ICC_ranks_S10 <- ICC_ranks_S10 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S10 <- as.numeric(ICC_ranks_S10$ICC_Scen10)
names(ICC_ranks_S10) <- data$DNAm.probe
ICC_ranks_S10 <- sort(ICC_ranks_S10, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S10 <- fgseaSimple(pathways = probe_set, 
                            stats = ICC_ranks_S10,
                            nperm = 10000)

#####################
### Old Salas CpG ###
#####################
### Probe Set Enrichment Analysis ### ----
library(FlowSorted.Blood.EPIC)
library(FlowSorted.BloodExtended.EPIC)
## Scenario 1 ##
# Read in probes #
data(IDOLOptimizedCpGs)
probe_set <- list(probes = IDOLOptimizedCpGs)
ICC_ranks_S1 <- data[,c(1,2)]
ICC_ranks_S1 <- ICC_ranks_S1 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S1 <- as.numeric(ICC_ranks_S1$ICC_Scen1)
names(ICC_ranks_S1) <- data$DNAm.probe
ICC_ranks_S1 <- sort(ICC_ranks_S1, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S1 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S1,
                           nperm = 10000)

## Scenario 2 ##
# Read in probes #
data(IDOLOptimizedCpGs)
probe_set <- list(probes = IDOLOptimizedCpGs)
ICC_ranks_S2 <- data[,c(1,3)]
ICC_ranks_S2 <- ICC_ranks_S2 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S2 <- as.numeric(ICC_ranks_S2$ICC_Scen2)
names(ICC_ranks_S2) <- data$DNAm.probe
ICC_ranks_S2 <- sort(ICC_ranks_S2, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S2 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S2,
                           nperm = 10000)

## Scenario 3 ##
# Read in probes #
data(IDOLOptimizedCpGs)
probe_set <- list(probes = IDOLOptimizedCpGs)
ICC_ranks_S3 <- data[,c(1,4)]
ICC_ranks_S3 <- ICC_ranks_S3 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S3 <- as.numeric(ICC_ranks_S3$ICC_Scen3)
names(ICC_ranks_S3) <- data$DNAm.probe
ICC_ranks_S3 <- sort(ICC_ranks_S3, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S3 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S3,
                           nperm = 10000)

## Scenario 4 ##
# Read in probes #
data(IDOLOptimizedCpGs)
probe_set <- list(probes = IDOLOptimizedCpGs)
ICC_ranks_S4 <- data[,c(1,5)]
ICC_ranks_S4 <- ICC_ranks_S4 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S4 <- as.numeric(ICC_ranks_S4$ICC_Scen4)
names(ICC_ranks_S4) <- data$DNAm.probe
ICC_ranks_S4 <- sort(ICC_ranks_S4, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S4 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S4,
                           nperm = 10000)

## Scenario 5 ##
# Read in probes #
data(IDOLOptimizedCpGs)
probe_set <- list(probes = IDOLOptimizedCpGs)
ICC_ranks_S5 <- data[,c(1,6)]
ICC_ranks_S5 <- ICC_ranks_S5 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S5 <- as.numeric(ICC_ranks_S5$ICC_Scen5)
names(ICC_ranks_S5) <- data$DNAm.probe
ICC_ranks_S5 <- sort(ICC_ranks_S5, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S5 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S5,
                           nperm = 10000)

## Scenario 6 ##
# Read in probes #
data(IDOLOptimizedCpGs)
probe_set <- list(probes = IDOLOptimizedCpGs)
ICC_ranks_S6 <- data[,c(1,7)]
ICC_ranks_S6 <- ICC_ranks_S6 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S6 <- as.numeric(ICC_ranks_S6$ICC_Scen6)
names(ICC_ranks_S6) <- data$DNAm.probe
ICC_ranks_S6 <- sort(ICC_ranks_S6, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S6 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S6,
                           nperm = 10000)

## Scenario 7 ##
# Read in probes #
data(IDOLOptimizedCpGs)
probe_set <- list(probes = IDOLOptimizedCpGs)
ICC_ranks_S7 <- data[,c(1,8)]
ICC_ranks_S7 <- ICC_ranks_S7 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S7 <- as.numeric(ICC_ranks_S7$ICC_Scen7)
names(ICC_ranks_S7) <- data$DNAm.probe
ICC_ranks_S7 <- sort(ICC_ranks_S7, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S7 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S7,
                           nperm = 10000)

## Scenario 8 ##
# Read in probes #
data(IDOLOptimizedCpGs)
probe_set <- list(probes = IDOLOptimizedCpGs)
ICC_ranks_S8 <- data[,c(1,9)]
ICC_ranks_S8 <- ICC_ranks_S8 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S8 <- as.numeric(ICC_ranks_S8$ICC_Scen8)
names(ICC_ranks_S8) <- data$DNAm.probe
ICC_ranks_S8 <- sort(ICC_ranks_S8, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S8 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S8,
                           nperm = 10000)

## Scenario 9 ##
# Read in probes #
data(IDOLOptimizedCpGs)
probe_set <- list(probes = IDOLOptimizedCpGs)
ICC_ranks_S9 <- data[,c(1,10)]
ICC_ranks_S9 <- ICC_ranks_S9 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S9 <- as.numeric(ICC_ranks_S9$ICC_Scen9)
names(ICC_ranks_S9) <- data$DNAm.probe
ICC_ranks_S9 <- sort(ICC_ranks_S9, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S9 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S9,
                           nperm = 10000)

## Scenario 10 ##
# Read in probes #
data(IDOLOptimizedCpGs)
probe_set <- list(probes = IDOLOptimizedCpGs)
ICC_ranks_S10 <- data[,c(1,11)]
ICC_ranks_S10 <- ICC_ranks_S10 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S10 <- as.numeric(ICC_ranks_S10$ICC_Scen10)
names(ICC_ranks_S10) <- data$DNAm.probe
ICC_ranks_S10 <- sort(ICC_ranks_S10, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S10 <- fgseaSimple(pathways = probe_set, 
                            stats = ICC_ranks_S10,
                            nperm = 10000)

#####################
### IDOL-Extended ###
#####################
### Probe Set Enrichment Analysis ### ----
## Scenario 1 ##
# Read in probes #
data(IDOLOptimizedCpGsBloodExtended)
probe_set <- list(probes = IDOLOptimizedCpGsBloodExtended)
ICC_ranks_S1 <- data[,c(1,2)]
ICC_ranks_S1 <- ICC_ranks_S1 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S1 <- as.numeric(ICC_ranks_S1$ICC_Scen1)
names(ICC_ranks_S1) <- data$DNAm.probe
ICC_ranks_S1 <- sort(ICC_ranks_S1, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S1 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S1,
                           nperm = 10000)

## Scenario 2 ##
# Read in probes #
data(IDOLOptimizedCpGsBloodExtended)
probe_set <- list(probes = IDOLOptimizedCpGsBloodExtended)
ICC_ranks_S2 <- data[,c(1,3)]
ICC_ranks_S2 <- ICC_ranks_S2 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S2 <- as.numeric(ICC_ranks_S2$ICC_Scen2)
names(ICC_ranks_S2) <- data$DNAm.probe
ICC_ranks_S2 <- sort(ICC_ranks_S2, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S2 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S2,
                           nperm = 10000)

## Scenario 3 ##
# Read in probes #
data(IDOLOptimizedCpGsBloodExtended)
probe_set <- list(probes = IDOLOptimizedCpGsBloodExtended)
ICC_ranks_S3 <- data[,c(1,4)]
ICC_ranks_S3 <- ICC_ranks_S3 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S3 <- as.numeric(ICC_ranks_S3$ICC_Scen3)
names(ICC_ranks_S3) <- data$DNAm.probe
ICC_ranks_S3 <- sort(ICC_ranks_S3, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S3 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S3,
                           nperm = 10000)

## Scenario 4 ##
# Read in probes #
data(IDOLOptimizedCpGsBloodExtended)
probe_set <- list(probes = IDOLOptimizedCpGsBloodExtended)
ICC_ranks_S4 <- data[,c(1,5)]
ICC_ranks_S4 <- ICC_ranks_S4 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S4 <- as.numeric(ICC_ranks_S4$ICC_Scen4)
names(ICC_ranks_S4) <- data$DNAm.probe
ICC_ranks_S4 <- sort(ICC_ranks_S4, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S4 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S4,
                           nperm = 10000)

## Scenario 5 ##
# Read in probes #
data(IDOLOptimizedCpGsBloodExtended)
probe_set <- list(probes = IDOLOptimizedCpGsBloodExtended)
ICC_ranks_S5 <- data[,c(1,6)]
ICC_ranks_S5 <- ICC_ranks_S5 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S5 <- as.numeric(ICC_ranks_S5$ICC_Scen5)
names(ICC_ranks_S5) <- data$DNAm.probe
ICC_ranks_S5 <- sort(ICC_ranks_S5, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S5 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S5,
                           nperm = 10000)

## Scenario 6 ##
# Read in probes #
data(IDOLOptimizedCpGsBloodExtended)
probe_set <- list(probes = IDOLOptimizedCpGsBloodExtended)
ICC_ranks_S6 <- data[,c(1,7)]
ICC_ranks_S6 <- ICC_ranks_S6 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S6 <- as.numeric(ICC_ranks_S6$ICC_Scen6)
names(ICC_ranks_S6) <- data$DNAm.probe
ICC_ranks_S6 <- sort(ICC_ranks_S6, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S6 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S6,
                           nperm = 10000)

## Scenario 7 ##
# Read in probes #
data(IDOLOptimizedCpGsBloodExtended)
probe_set <- list(probes = IDOLOptimizedCpGsBloodExtended)
ICC_ranks_S7 <- data[,c(1,8)]
ICC_ranks_S7 <- ICC_ranks_S7 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S7 <- as.numeric(ICC_ranks_S7$ICC_Scen7)
names(ICC_ranks_S7) <- data$DNAm.probe
ICC_ranks_S7 <- sort(ICC_ranks_S7, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S7 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S7,
                           nperm = 10000)

## Scenario 8 ##
# Read in probes #
data(IDOLOptimizedCpGsBloodExtended)
probe_set <- list(probes = IDOLOptimizedCpGsBloodExtended)
ICC_ranks_S8 <- data[,c(1,9)]
ICC_ranks_S8 <- ICC_ranks_S8 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S8 <- as.numeric(ICC_ranks_S8$ICC_Scen8)
names(ICC_ranks_S8) <- data$DNAm.probe
ICC_ranks_S8 <- sort(ICC_ranks_S8, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S8 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S8,
                           nperm = 10000)

## Scenario 9 ##
# Read in probes #
data(IDOLOptimizedCpGsBloodExtended)
probe_set <- list(probes = IDOLOptimizedCpGsBloodExtended)
ICC_ranks_S9 <- data[,c(1,10)]
ICC_ranks_S9 <- ICC_ranks_S9 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S9 <- as.numeric(ICC_ranks_S9$ICC_Scen9)
names(ICC_ranks_S9) <- data$DNAm.probe
ICC_ranks_S9 <- sort(ICC_ranks_S9, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S9 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S9,
                           nperm = 10000)

## Scenario 10 ##
# Read in probes #
data(IDOLOptimizedCpGsBloodExtended)
probe_set <- list(probes = IDOLOptimizedCpGsBloodExtended)
ICC_ranks_S10 <- data[,c(1,11)]
ICC_ranks_S10 <- ICC_ranks_S10 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S10 <- as.numeric(ICC_ranks_S10$ICC_Scen10)
names(ICC_ranks_S10) <- data$DNAm.probe
ICC_ranks_S10 <- sort(ICC_ranks_S10, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S10 <- fgseaSimple(pathways = probe_set, 
                            stats = ICC_ranks_S10,
                            nperm = 10000)

#########################################################
#################### Generate Figure ####################
#########################################################

### Create Boxplots for Each Scenario ### ----
Horvath_probes <- read.csv("Horvath_Probes.csv")
data_All <- data[data$DNAm.probe %in% Horvath_probes$DNAm.probe,]
data_All$Name <- as.factor("Horvath")

Hannum_probes <- read.csv("Hannum_Probes.csv")
data_Hannum <- data[data$DNAm.probe %in% Hannum_probes$DNAm.probe,]
data_Hannum$Name <- as.factor("Hannum")
data_All <- rbind(data_All, data_Hannum)

PhenoAge_probes <- read.csv("PhenoAge_Probes.csv")
data_PhenoAge <- data[data$DNAm.probe %in% PhenoAge_probes$DNAm.probe,]
data_PhenoAge$Name <- as.factor("PhenoAge")
data_All <- rbind(data_All, data_PhenoAge)

PC_Clocks_Probes <- read.csv("PC_Probes.csv")
data_PC_Clocks <- data[data$DNAm.probe %in% PC_Clocks_Probes$DNAm.probe,]
data_PC_Clocks$Name <- as.factor("PC Clocks")
data_All <- rbind(data_All, data_PC_Clocks)

load("sysdata.rda")
data_PACE <- data[data$DNAm.probe %in% mPOA_Models$model_probes$DunedinPACE,]
data_PACE$Name <- as.factor("PACE")
data_All <- rbind(data_All, data_PACE)

data(IDOLOptimizedCpGs)
data_IDOL6 <- data[data$DNAm.probe %in% IDOLOptimizedCpGs,]
data_IDOL6$Name <- as.factor("IDOL-6")
data_All <- rbind(data_All, data_IDOL6)

data(IDOLOptimizedCpGsBloodExtended)
data_IDOLExt <- data[data$DNAm.probe %in% IDOLOptimizedCpGsBloodExtended,]
data_IDOLExt$Name <- as.factor("IDOL-Extended")
data_All <- rbind(data_All, data_IDOLExt)

# Exclude probes with ICC < 0.01 #
data_All[,2:11][data_All[,2:11] < 0.01] <- NA

Box1 <- ggplot(data_All, aes(x=Name, y=ICC_Scen1, fill=Name)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0.47,linetype=2) +
  labs(title="NoStressT1-2", x =NULL, y = NULL) + 
  scale_fill_brewer(palette="Set2") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
Box2 <- ggplot(data_All, aes(x=Name, y=ICC_Scen2, fill=Name)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0.45,linetype=2) +
  labs(title="NoStressT1-3", x =NULL, y = NULL) + 
  scale_fill_brewer(palette="Set2") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
Box3 <- ggplot(data_All, aes(x=Name, y=ICC_Scen3, fill=Name)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0.31,linetype=2) +
  labs(title="NoStressT1-4", x =NULL, y = NULL) + 
  scale_fill_brewer(palette="Set2") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
Box4 <- ggplot(data_All, aes(x=Name, y=ICC_Scen4, fill=Name)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0.50,linetype=2) +
  labs(title="NoStressT3-4", x =NULL, y = NULL) + 
  scale_fill_brewer(palette="Set2") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
Box5 <- ggplot(data_All, aes(x=Name, y=ICC_Scen5, fill=Name)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0.44,linetype=2) +
  labs(title="NoStressT1-2-3-4", x =NULL, y = NULL) + 
  scale_fill_brewer(palette="Set2") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
Box6 <- ggplot(data_All, aes(x=Name, y=ICC_Scen6samp, fill=Name)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0.64,linetype=2) +
  labs(title="StressT1-2", x =NULL, y = NULL) + 
  scale_fill_brewer(palette="Set2") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
Box7 <- ggplot(data_All, aes(x=Name, y=ICC_Scen7samp, fill=Name)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0.64,linetype=2) +
  labs(title="StressT1-3", x =NULL, y = NULL) + 
  scale_fill_brewer(palette="Set2") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
Box8 <- ggplot(data_All, aes(x=Name, y=ICC_Scen8samp, fill=Name)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0.65,linetype=2) +
  labs(title="StressT1-4", x =NULL, y = NULL) + 
  scale_fill_brewer(palette="Set2") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
Box9 <- ggplot(data_All, aes(x=Name, y=ICC_Scen9samp, fill=Name)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0.62,linetype=2) +
  labs(title="StressT3-4", x =NULL, y = NULL) + 
  scale_fill_brewer(palette="Set2") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
Box10 <- ggplot(data_All, aes(x=Name, y=ICC_Scen10, fill=Name)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0.28,linetype=2) +
  labs(title="CrossSessionT1", x =NULL, y = NULL) + 
  scale_fill_brewer(palette="Set2") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

### Combine All Plots ### ----
ggarrange(Box1, Box2, Box3, Box4, Box5,
          Box6, Box7, Box8, Box9, Box10,
          ncol = 5, nrow = 2, legend = "right", common.legend = TRUE)



##################################################################################################
################################## Biologically Relevant Probes ##################################
##################################################################################################
### Set Seed ### ----
set.seed(100)

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
# Exclude probes with ICC < 0.01 #
data[,2:11][data[,2:11] < 0.01] <- NA
# Define annEPIC #
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

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


#####################
#### mQTL Probes ####
#####################
### Probe Set Enrichment Analysis ### ----
## Scenario 1 ##
# Read in probes #
probe_set <- read.csv("probe_set.csv")
probe_set <- sample(probe_set$DNAm.probe, 35000)
probe_set <- list(probes = probe_set)
ICC_ranks_S1 <- data[,c(1,2)]
ICC_ranks_S1 <- ICC_ranks_S1 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S1 <- as.numeric(ICC_ranks_S1$ICC_Scen1)
names(ICC_ranks_S1) <- data$DNAm.probe
ICC_ranks_S1 <- sort(ICC_ranks_S1, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S1 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S1,
                           nperm = 10000)

## Scenario 2 ##
# Read in probes #
probe_set <- read.csv("probe_set.csv")
probe_set <- sample(probe_set$DNAm.probe, 35000)
probe_set <- list(probes = probe_set)
ICC_ranks_S2 <- data[,c(1,3)]
ICC_ranks_S2 <- ICC_ranks_S2 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S2 <- as.numeric(ICC_ranks_S2$ICC_Scen2)
names(ICC_ranks_S2) <- data$DNAm.probe
ICC_ranks_S2 <- sort(ICC_ranks_S2, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S2 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S2,
                           nperm = 10000)

## Scenario 3 ##
# Read in probes #
probe_set <- read.csv("probe_set.csv")
probe_set <- sample(probe_set$DNAm.probe, 35000)
probe_set <- list(probes = probe_set)
ICC_ranks_S3 <- data[,c(1,4)]
ICC_ranks_S3 <- ICC_ranks_S3 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S3 <- as.numeric(ICC_ranks_S3$ICC_Scen3)
names(ICC_ranks_S3) <- data$DNAm.probe
ICC_ranks_S3 <- sort(ICC_ranks_S3, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S3 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S3,
                           nperm = 10000)

## Scenario 4 ##
# Read in probes #
probe_set <- read.csv("probe_set.csv")
probe_set <- sample(probe_set$DNAm.probe, 35000)
probe_set <- list(probes = probe_set)
ICC_ranks_S4 <- data[,c(1,5)]
ICC_ranks_S4 <- ICC_ranks_S4 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S4 <- as.numeric(ICC_ranks_S4$ICC_Scen4)
names(ICC_ranks_S4) <- data$DNAm.probe
ICC_ranks_S4 <- sort(ICC_ranks_S4, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S4 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S4,
                           nperm = 10000)

## Scenario 5 ##
# Read in probes #
probe_set <- read.csv("probe_set.csv")
probe_set <- sample(probe_set$DNAm.probe, 35000)
probe_set <- list(probes = probe_set)
ICC_ranks_S5 <- data[,c(1,6)]
ICC_ranks_S5 <- ICC_ranks_S5 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S5 <- as.numeric(ICC_ranks_S5$ICC_Scen5)
names(ICC_ranks_S5) <- data$DNAm.probe
ICC_ranks_S5 <- sort(ICC_ranks_S5, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S5 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S5,
                           nperm = 10000)

## Scenario 6 ##
# Read in probes #
probe_set <- read.csv("probe_set.csv")
probe_set <- sample(probe_set$DNAm.probe, 35000)
probe_set <- list(probes = probe_set)
ICC_ranks_S6 <- data[,c(1,7)]
ICC_ranks_S6 <- ICC_ranks_S6 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S6 <- as.numeric(ICC_ranks_S6$ICC_Scen6)
names(ICC_ranks_S6) <- data$DNAm.probe
ICC_ranks_S6 <- sort(ICC_ranks_S6, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S6 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S6,
                           nperm = 10000)

## Scenario 7 ##
# Read in probes #
probe_set <- read.csv("probe_set.csv")
probe_set <- sample(probe_set$DNAm.probe, 35000)
probe_set <- list(probes = probe_set)
ICC_ranks_S7 <- data[,c(1,8)]
ICC_ranks_S7 <- ICC_ranks_S7 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S7 <- as.numeric(ICC_ranks_S7$ICC_Scen7)
names(ICC_ranks_S7) <- data$DNAm.probe
ICC_ranks_S7 <- sort(ICC_ranks_S7, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S7 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S7,
                           nperm = 10000)

## Scenario 8 ##
# Read in probes #
probe_set <- read.csv("probe_set.csv")
probe_set <- sample(probe_set$DNAm.probe, 35000)
probe_set <- list(probes = probe_set)
ICC_ranks_S8 <- data[,c(1,9)]
ICC_ranks_S8 <- ICC_ranks_S8 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S8 <- as.numeric(ICC_ranks_S8$ICC_Scen8)
names(ICC_ranks_S8) <- data$DNAm.probe
ICC_ranks_S8 <- sort(ICC_ranks_S8, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S8 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S8,
                           nperm = 10000)

## Scenario 9 ##
# Read in probes #
probe_set <- read.csv("probe_set.csv")
probe_set <- sample(probe_set$DNAm.probe, 35000)
probe_set <- list(probes = probe_set)
ICC_ranks_S9 <- data[,c(1,10)]
ICC_ranks_S9 <- ICC_ranks_S9 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S9 <- as.numeric(ICC_ranks_S9$ICC_Scen9)
names(ICC_ranks_S9) <- data$DNAm.probe
ICC_ranks_S9 <- sort(ICC_ranks_S9, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S9 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S9,
                           nperm = 10000)

## Scenario 10 ##
# Read in probes #
probe_set <- read.csv("probe_set.csv")
probe_set <- sample(probe_set$DNAm.probe, 35000)
probe_set <- list(probes = probe_set)
ICC_ranks_S10 <- data[,c(1,11)]
ICC_ranks_S10 <- ICC_ranks_S10 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S10 <- as.numeric(ICC_ranks_S10$ICC_Scen10)
names(ICC_ranks_S10) <- data$DNAm.probe
ICC_ranks_S10 <- sort(ICC_ranks_S10, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S10 <- fgseaSimple(pathways = probe_set, 
                            stats = ICC_ranks_S10,
                            nperm = 10000)

###########
### TSS ###
###########
### Probe Set Enrichment Analysis ### ----
## Scenario 1 ##
# Read in probes #
keep <- grepl("TSS", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
TSS_Probes <- annEPIC@listData$Name[keep]
TSS_Probes <- sample(TSS_Probes, 35000)
probe_set <- list(probes = TSS_Probes)
ICC_ranks_S1 <- data[,c(1,2)]
ICC_ranks_S1 <- ICC_ranks_S1 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S1 <- as.numeric(ICC_ranks_S1$ICC_Scen1)
names(ICC_ranks_S1) <- data$DNAm.probe
ICC_ranks_S1 <- sort(ICC_ranks_S1, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S1 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S1,
                           nperm = 10000)

## Scenario 2 ##
# Read in probes #
keep <- grepl("TSS", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
TSS_Probes <- annEPIC@listData$Name[keep]
TSS_Probes <- sample(TSS_Probes, 35000)
probe_set <- list(probes = TSS_Probes)
ICC_ranks_S2 <- data[,c(1,3)]
ICC_ranks_S2 <- ICC_ranks_S2 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S2 <- as.numeric(ICC_ranks_S2$ICC_Scen2)
names(ICC_ranks_S2) <- data$DNAm.probe
ICC_ranks_S2 <- sort(ICC_ranks_S2, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S2 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S2,
                           nperm = 10000)

## Scenario 3 ##
# Read in probes #
keep <- grepl("TSS", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
TSS_Probes <- annEPIC@listData$Name[keep]
TSS_Probes <- sample(TSS_Probes, 35000)
probe_set <- list(probes = TSS_Probes)
ICC_ranks_S3 <- data[,c(1,4)]
ICC_ranks_S3 <- ICC_ranks_S3 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S3 <- as.numeric(ICC_ranks_S3$ICC_Scen3)
names(ICC_ranks_S3) <- data$DNAm.probe
ICC_ranks_S3 <- sort(ICC_ranks_S3, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S3 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S3,
                           nperm = 10000)

## Scenario 4 ##
# Read in probes #
keep <- grepl("TSS", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
TSS_Probes <- annEPIC@listData$Name[keep]
TSS_Probes <- sample(TSS_Probes, 35000)
probe_set <- list(probes = TSS_Probes)
ICC_ranks_S4 <- data[,c(1,5)]
ICC_ranks_S4 <- ICC_ranks_S4 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S4 <- as.numeric(ICC_ranks_S4$ICC_Scen4)
names(ICC_ranks_S4) <- data$DNAm.probe
ICC_ranks_S4 <- sort(ICC_ranks_S4, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S4 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S4,
                           nperm = 10000)

## Scenario 5 ##
# Read in probes #
keep <- grepl("TSS", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
TSS_Probes <- annEPIC@listData$Name[keep]
TSS_Probes <- sample(TSS_Probes, 35000)
probe_set <- list(probes = TSS_Probes)
ICC_ranks_S5 <- data[,c(1,6)]
ICC_ranks_S5 <- ICC_ranks_S5 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S5 <- as.numeric(ICC_ranks_S5$ICC_Scen5)
names(ICC_ranks_S5) <- data$DNAm.probe
ICC_ranks_S5 <- sort(ICC_ranks_S5, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S5 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S5,
                           nperm = 10000)

## Scenario 6 ##
# Read in probes #
keep <- grepl("TSS", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
TSS_Probes <- annEPIC@listData$Name[keep]
TSS_Probes <- sample(TSS_Probes, 35000)
probe_set <- list(probes = TSS_Probes)
ICC_ranks_S6 <- data[,c(1,7)]
ICC_ranks_S6 <- ICC_ranks_S6 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S6 <- as.numeric(ICC_ranks_S6$ICC_Scen6)
names(ICC_ranks_S6) <- data$DNAm.probe
ICC_ranks_S6 <- sort(ICC_ranks_S6, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S6 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S6,
                           nperm = 10000)

## Scenario 7 ##
# Read in probes #
keep <- grepl("TSS", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
TSS_Probes <- annEPIC@listData$Name[keep]
TSS_Probes <- sample(TSS_Probes, 35000)
probe_set <- list(probes = TSS_Probes)
ICC_ranks_S7 <- data[,c(1,8)]
ICC_ranks_S7 <- ICC_ranks_S7 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S7 <- as.numeric(ICC_ranks_S7$ICC_Scen7)
names(ICC_ranks_S7) <- data$DNAm.probe
ICC_ranks_S7 <- sort(ICC_ranks_S7, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S7 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S7,
                           nperm = 10000)

## Scenario 8 ##
# Read in probes #
keep <- grepl("TSS", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
TSS_Probes <- annEPIC@listData$Name[keep]
TSS_Probes <- sample(TSS_Probes, 35000)
probe_set <- list(probes = TSS_Probes)
ICC_ranks_S8 <- data[,c(1,9)]
ICC_ranks_S8 <- ICC_ranks_S8 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S8 <- as.numeric(ICC_ranks_S8$ICC_Scen8)
names(ICC_ranks_S8) <- data$DNAm.probe
ICC_ranks_S8 <- sort(ICC_ranks_S8, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S8 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S8,
                           nperm = 10000)

## Scenario 9 ##
# Read in probes #
keep <- grepl("TSS", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
TSS_Probes <- annEPIC@listData$Name[keep]
TSS_Probes <- sample(TSS_Probes, 35000)
probe_set <- list(probes = TSS_Probes)
ICC_ranks_S9 <- data[,c(1,10)]
ICC_ranks_S9 <- ICC_ranks_S9 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S9 <- as.numeric(ICC_ranks_S9$ICC_Scen9)
names(ICC_ranks_S9) <- data$DNAm.probe
ICC_ranks_S9 <- sort(ICC_ranks_S9, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S9 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S9,
                           nperm = 10000)

## Scenario 10 ##
# Read in probes #
keep <- grepl("TSS", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
TSS_Probes <- annEPIC@listData$Name[keep]
TSS_Probes <- sample(TSS_Probes, 35000)
probe_set <- list(probes = TSS_Probes)
ICC_ranks_S10 <- data[,c(1,11)]
ICC_ranks_S10 <- ICC_ranks_S10 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S10 <- as.numeric(ICC_ranks_S10$ICC_Scen10)
names(ICC_ranks_S10) <- data$DNAm.probe
ICC_ranks_S10 <- sort(ICC_ranks_S10, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S10 <- fgseaSimple(pathways = probe_set, 
                            stats = ICC_ranks_S10,
                            nperm = 10000)

#############
### 5'UTR ###
#############
### Probe Set Enrichment Analysis ### ----
## Scenario 1 ##
# Read in probes #
keep <- grepl("5'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Five_UTR_Probes <- annEPIC@listData$Name[keep]
Five_UTR_Probes <- sample(Five_UTR_Probes, 35000)
probe_set <- list(probes = Five_UTR_Probes)
ICC_ranks_S1 <- data[,c(1,2)]
ICC_ranks_S1 <- ICC_ranks_S1 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S1 <- as.numeric(ICC_ranks_S1$ICC_Scen1)
names(ICC_ranks_S1) <- data$DNAm.probe
ICC_ranks_S1 <- sort(ICC_ranks_S1, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S1 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S1,
                           nperm = 10000)

## Scenario 2 ##
# Read in probes #
keep <- grepl("5'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Five_UTR_Probes <- annEPIC@listData$Name[keep]
Five_UTR_Probes <- sample(Five_UTR_Probes, 35000)
probe_set <- list(probes = Five_UTR_Probes)
ICC_ranks_S2 <- data[,c(1,3)]
ICC_ranks_S2 <- ICC_ranks_S2 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S2 <- as.numeric(ICC_ranks_S2$ICC_Scen2)
names(ICC_ranks_S2) <- data$DNAm.probe
ICC_ranks_S2 <- sort(ICC_ranks_S2, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S2 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S2,
                           nperm = 10000)

## Scenario 3 ##
# Read in probes #
keep <- grepl("5'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Five_UTR_Probes <- annEPIC@listData$Name[keep]
Five_UTR_Probes <- sample(Five_UTR_Probes, 35000)
probe_set <- list(probes = Five_UTR_Probes)
ICC_ranks_S3 <- data[,c(1,4)]
ICC_ranks_S3 <- ICC_ranks_S3 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S3 <- as.numeric(ICC_ranks_S3$ICC_Scen3)
names(ICC_ranks_S3) <- data$DNAm.probe
ICC_ranks_S3 <- sort(ICC_ranks_S3, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S3 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S3,
                           nperm = 10000)

## Scenario 4 ##
# Read in probes #
keep <- grepl("5'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Five_UTR_Probes <- annEPIC@listData$Name[keep]
Five_UTR_Probes <- sample(Five_UTR_Probes, 35000)
probe_set <- list(probes = Five_UTR_Probes)
ICC_ranks_S4 <- data[,c(1,5)]
ICC_ranks_S4 <- ICC_ranks_S4 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S4 <- as.numeric(ICC_ranks_S4$ICC_Scen4)
names(ICC_ranks_S4) <- data$DNAm.probe
ICC_ranks_S4 <- sort(ICC_ranks_S4, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S4 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S4,
                           nperm = 10000)

## Scenario 5 ##
# Read in probes #
keep <- grepl("5'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Five_UTR_Probes <- annEPIC@listData$Name[keep]
Five_UTR_Probes <- sample(Five_UTR_Probes, 35000)
probe_set <- list(probes = Five_UTR_Probes)
ICC_ranks_S5 <- data[,c(1,6)]
ICC_ranks_S5 <- ICC_ranks_S5 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S5 <- as.numeric(ICC_ranks_S5$ICC_Scen5)
names(ICC_ranks_S5) <- data$DNAm.probe
ICC_ranks_S5 <- sort(ICC_ranks_S5, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S5 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S5,
                           nperm = 10000)

## Scenario 6 ##
# Read in probes #
keep <- grepl("5'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Five_UTR_Probes <- annEPIC@listData$Name[keep]
Five_UTR_Probes <- sample(Five_UTR_Probes, 35000)
probe_set <- list(probes = Five_UTR_Probes)
ICC_ranks_S6 <- data[,c(1,7)]
ICC_ranks_S6 <- ICC_ranks_S6 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S6 <- as.numeric(ICC_ranks_S6$ICC_Scen6)
names(ICC_ranks_S6) <- data$DNAm.probe
ICC_ranks_S6 <- sort(ICC_ranks_S6, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S6 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S6,
                           nperm = 10000)

## Scenario 7 ##
# Read in probes #
keep <- grepl("5'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Five_UTR_Probes <- annEPIC@listData$Name[keep]
Five_UTR_Probes <- sample(Five_UTR_Probes, 35000)
probe_set <- list(probes = Five_UTR_Probes)
ICC_ranks_S7 <- data[,c(1,8)]
ICC_ranks_S7 <- ICC_ranks_S7 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S7 <- as.numeric(ICC_ranks_S7$ICC_Scen7)
names(ICC_ranks_S7) <- data$DNAm.probe
ICC_ranks_S7 <- sort(ICC_ranks_S7, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S7 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S7,
                           nperm = 10000)

## Scenario 8 ##
# Read in probes #
keep <- grepl("5'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Five_UTR_Probes <- annEPIC@listData$Name[keep]
Five_UTR_Probes <- sample(Five_UTR_Probes, 35000)
probe_set <- list(probes = Five_UTR_Probes)
ICC_ranks_S8 <- data[,c(1,9)]
ICC_ranks_S8 <- ICC_ranks_S8 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S8 <- as.numeric(ICC_ranks_S8$ICC_Scen8)
names(ICC_ranks_S8) <- data$DNAm.probe
ICC_ranks_S8 <- sort(ICC_ranks_S8, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S8 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S8,
                           nperm = 10000)

## Scenario 9 ##
# Read in probes #
keep <- grepl("5'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Five_UTR_Probes <- annEPIC@listData$Name[keep]
Five_UTR_Probes <- sample(Five_UTR_Probes, 35000)
probe_set <- list(probes = Five_UTR_Probes)
ICC_ranks_S9 <- data[,c(1,10)]
ICC_ranks_S9 <- ICC_ranks_S9 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S9 <- as.numeric(ICC_ranks_S9$ICC_Scen9)
names(ICC_ranks_S9) <- data$DNAm.probe
ICC_ranks_S9 <- sort(ICC_ranks_S9, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S9 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S9,
                           nperm = 10000)

## Scenario 10 ##
# Read in probes #
keep <- grepl("5'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Five_UTR_Probes <- annEPIC@listData$Name[keep]
Five_UTR_Probes <- sample(Five_UTR_Probes, 35000)
probe_set <- list(probes = Five_UTR_Probes)
ICC_ranks_S10 <- data[,c(1,11)]
ICC_ranks_S10 <- ICC_ranks_S10 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S10 <- as.numeric(ICC_ranks_S10$ICC_Scen10)
names(ICC_ranks_S10) <- data$DNAm.probe
ICC_ranks_S10 <- sort(ICC_ranks_S10, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S10 <- fgseaSimple(pathways = probe_set, 
                            stats = ICC_ranks_S10,
                            nperm = 10000)

#############
### 3'UTR ###
#############
### Probe Set Enrichment Analysis ### ----
## Scenario 1 ##
# Read in probes #
keep <- grepl("3'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Three_UTR_Probes <- annEPIC@listData$Name[keep]
Three_UTR_Probes <- sample(Three_UTR_Probes, 35000)
probe_set <- list(probes = Three_UTR_Probes)
ICC_ranks_S1 <- data[,c(1,2)]
ICC_ranks_S1 <- ICC_ranks_S1 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S1 <- as.numeric(ICC_ranks_S1$ICC_Scen1)
names(ICC_ranks_S1) <- data$DNAm.probe
ICC_ranks_S1 <- sort(ICC_ranks_S1, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S1 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S1,
                           nperm = 10000)

## Scenario 2 ##
# Read in probes #
keep <- grepl("3'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Three_UTR_Probes <- annEPIC@listData$Name[keep]
Three_UTR_Probes <- sample(Three_UTR_Probes, 35000)
probe_set <- list(probes = Three_UTR_Probes)
ICC_ranks_S2 <- data[,c(1,3)]
ICC_ranks_S2 <- ICC_ranks_S2 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S2 <- as.numeric(ICC_ranks_S2$ICC_Scen2)
names(ICC_ranks_S2) <- data$DNAm.probe
ICC_ranks_S2 <- sort(ICC_ranks_S2, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S2 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S2,
                           nperm = 10000)

## Scenario 3 ##
# Read in probes #
keep <- grepl("3'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Three_UTR_Probes <- annEPIC@listData$Name[keep]
Three_UTR_Probes <- sample(Three_UTR_Probes, 35000)
probe_set <- list(probes = Three_UTR_Probes)
ICC_ranks_S3 <- data[,c(1,4)]
ICC_ranks_S3 <- ICC_ranks_S3 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S3 <- as.numeric(ICC_ranks_S3$ICC_Scen3)
names(ICC_ranks_S3) <- data$DNAm.probe
ICC_ranks_S3 <- sort(ICC_ranks_S3, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S3 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S3,
                           nperm = 10000)

## Scenario 4 ##
# Read in probes #
keep <- grepl("3'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Three_UTR_Probes <- annEPIC@listData$Name[keep]
Three_UTR_Probes <- sample(Three_UTR_Probes, 35000)
probe_set <- list(probes = Three_UTR_Probes)
ICC_ranks_S4 <- data[,c(1,5)]
ICC_ranks_S4 <- ICC_ranks_S4 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S4 <- as.numeric(ICC_ranks_S4$ICC_Scen4)
names(ICC_ranks_S4) <- data$DNAm.probe
ICC_ranks_S4 <- sort(ICC_ranks_S4, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S4 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S4,
                           nperm = 10000)

## Scenario 5 ##
# Read in probes #
keep <- grepl("3'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Three_UTR_Probes <- annEPIC@listData$Name[keep]
Three_UTR_Probes <- sample(Three_UTR_Probes, 35000)
probe_set <- list(probes = Three_UTR_Probes)
ICC_ranks_S5 <- data[,c(1,6)]
ICC_ranks_S5 <- ICC_ranks_S5 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S5 <- as.numeric(ICC_ranks_S5$ICC_Scen5)
names(ICC_ranks_S5) <- data$DNAm.probe
ICC_ranks_S5 <- sort(ICC_ranks_S5, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S5 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S5,
                           nperm = 10000)

## Scenario 6 ##
# Read in probes #
keep <- grepl("3'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Three_UTR_Probes <- annEPIC@listData$Name[keep]
Three_UTR_Probes <- sample(Three_UTR_Probes, 35000)
probe_set <- list(probes = Three_UTR_Probes)
ICC_ranks_S6 <- data[,c(1,7)]
ICC_ranks_S6 <- ICC_ranks_S6 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S6 <- as.numeric(ICC_ranks_S6$ICC_Scen6)
names(ICC_ranks_S6) <- data$DNAm.probe
ICC_ranks_S6 <- sort(ICC_ranks_S6, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S6 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S6,
                           nperm = 10000)

## Scenario 7 ##
# Read in probes #
keep <- grepl("3'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Three_UTR_Probes <- annEPIC@listData$Name[keep]
Three_UTR_Probes <- sample(Three_UTR_Probes, 35000)
probe_set <- list(probes = Three_UTR_Probes)
ICC_ranks_S7 <- data[,c(1,8)]
ICC_ranks_S7 <- ICC_ranks_S7 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S7 <- as.numeric(ICC_ranks_S7$ICC_Scen7)
names(ICC_ranks_S7) <- data$DNAm.probe
ICC_ranks_S7 <- sort(ICC_ranks_S7, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S7 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S7,
                           nperm = 10000)

## Scenario 8 ##
# Read in probes #
keep <- grepl("3'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Three_UTR_Probes <- annEPIC@listData$Name[keep]
Three_UTR_Probes <- sample(Three_UTR_Probes, 35000)
probe_set <- list(probes = Three_UTR_Probes)
ICC_ranks_S8 <- data[,c(1,9)]
ICC_ranks_S8 <- ICC_ranks_S8 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S8 <- as.numeric(ICC_ranks_S8$ICC_Scen8)
names(ICC_ranks_S8) <- data$DNAm.probe
ICC_ranks_S8 <- sort(ICC_ranks_S8, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S8 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S8,
                           nperm = 10000)

## Scenario 9 ##
# Read in probes #
keep <- grepl("3'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Three_UTR_Probes <- annEPIC@listData$Name[keep]
Three_UTR_Probes <- sample(Three_UTR_Probes, 35000)
probe_set <- list(probes = Three_UTR_Probes)
ICC_ranks_S9 <- data[,c(1,10)]
ICC_ranks_S9 <- ICC_ranks_S9 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S9 <- as.numeric(ICC_ranks_S9$ICC_Scen9)
names(ICC_ranks_S9) <- data$DNAm.probe
ICC_ranks_S9 <- sort(ICC_ranks_S9, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S9 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S9,
                           nperm = 10000)

## Scenario 10 ##
# Read in probes #
keep <- grepl("3'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Three_UTR_Probes <- annEPIC@listData$Name[keep]
Three_UTR_Probes <- sample(Three_UTR_Probes, 35000)
probe_set <- list(probes = Three_UTR_Probes)
ICC_ranks_S10 <- data[,c(1,11)]
ICC_ranks_S10 <- ICC_ranks_S10 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S10 <- as.numeric(ICC_ranks_S10$ICC_Scen10)
names(ICC_ranks_S10) <- data$DNAm.probe
ICC_ranks_S10 <- sort(ICC_ranks_S10, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S10 <- fgseaSimple(pathways = probe_set, 
                            stats = ICC_ranks_S10,
                            nperm = 10000)

#####################
### Coding Region ###
#####################
### Probe Set Enrichment Analysis ### ----
## Scenario 1 ##
# Read in probes #
keep <- grepl("Exon", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep2 <- grepl("Body", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep <- keep | keep2
CodRegProbes <- annEPIC@listData$Name[keep]
CodRegProbes <- sample(CodRegProbes, 35000)
probe_set <- list(probes = CodRegProbes)
ICC_ranks_S1 <- data[,c(1,2)]
ICC_ranks_S1 <- ICC_ranks_S1 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S1 <- as.numeric(ICC_ranks_S1$ICC_Scen1)
names(ICC_ranks_S1) <- data$DNAm.probe
ICC_ranks_S1 <- sort(ICC_ranks_S1, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S1 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S1,
                           nperm = 10000)

## Scenario 2 ##
# Read in probes #
keep <- grepl("Exon", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep2 <- grepl("Body", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep <- keep | keep2
CodRegProbes <- annEPIC@listData$Name[keep]
CodRegProbes <- sample(CodRegProbes, 35000)
probe_set <- list(probes = CodRegProbes)
ICC_ranks_S2 <- data[,c(1,3)]
ICC_ranks_S2 <- ICC_ranks_S2 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S2 <- as.numeric(ICC_ranks_S2$ICC_Scen2)
names(ICC_ranks_S2) <- data$DNAm.probe
ICC_ranks_S2 <- sort(ICC_ranks_S2, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S2 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S2,
                           nperm = 10000)

## Scenario 3 ##
# Read in probes #
keep <- grepl("Exon", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep2 <- grepl("Body", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep <- keep | keep2
CodRegProbes <- annEPIC@listData$Name[keep]
CodRegProbes <- sample(CodRegProbes, 35000)
probe_set <- list(probes = CodRegProbes)
ICC_ranks_S3 <- data[,c(1,4)]
ICC_ranks_S3 <- ICC_ranks_S3 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S3 <- as.numeric(ICC_ranks_S3$ICC_Scen3)
names(ICC_ranks_S3) <- data$DNAm.probe
ICC_ranks_S3 <- sort(ICC_ranks_S3, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S3 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S3,
                           nperm = 10000)

## Scenario 4 ##
# Read in probes #
keep <- grepl("Exon", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep2 <- grepl("Body", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep <- keep | keep2
CodRegProbes <- annEPIC@listData$Name[keep]
CodRegProbes <- sample(CodRegProbes, 35000)
probe_set <- list(probes = CodRegProbes)
ICC_ranks_S4 <- data[,c(1,5)]
ICC_ranks_S4 <- ICC_ranks_S4 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S4 <- as.numeric(ICC_ranks_S4$ICC_Scen4)
names(ICC_ranks_S4) <- data$DNAm.probe
ICC_ranks_S4 <- sort(ICC_ranks_S4, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S4 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S4,
                           nperm = 10000)

## Scenario 5 ##
# Read in probes #
keep <- grepl("Exon", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep2 <- grepl("Body", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep <- keep | keep2
CodRegProbes <- annEPIC@listData$Name[keep]
CodRegProbes <- sample(CodRegProbes, 35000)
probe_set <- list(probes = CodRegProbes)
ICC_ranks_S5 <- data[,c(1,6)]
ICC_ranks_S5 <- ICC_ranks_S5 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S5 <- as.numeric(ICC_ranks_S5$ICC_Scen5)
names(ICC_ranks_S5) <- data$DNAm.probe
ICC_ranks_S5 <- sort(ICC_ranks_S5, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S5 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S5,
                           nperm = 10000)

## Scenario 6 ##
# Read in probes #
keep <- grepl("Exon", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep2 <- grepl("Body", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep <- keep | keep2
CodRegProbes <- annEPIC@listData$Name[keep]
CodRegProbes <- sample(CodRegProbes, 35000)
probe_set <- list(probes = CodRegProbes)
ICC_ranks_S6 <- data[,c(1,7)]
ICC_ranks_S6 <- ICC_ranks_S6 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S6 <- as.numeric(ICC_ranks_S6$ICC_Scen6)
names(ICC_ranks_S6) <- data$DNAm.probe
ICC_ranks_S6 <- sort(ICC_ranks_S6, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S6 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S6,
                           nperm = 10000)

## Scenario 7 ##
# Read in probes #
keep <- grepl("Exon", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep2 <- grepl("Body", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep <- keep | keep2
CodRegProbes <- annEPIC@listData$Name[keep]
CodRegProbes <- sample(CodRegProbes, 35000)
probe_set <- list(probes = CodRegProbes)
ICC_ranks_S7 <- data[,c(1,8)]
ICC_ranks_S7 <- ICC_ranks_S7 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S7 <- as.numeric(ICC_ranks_S7$ICC_Scen7)
names(ICC_ranks_S7) <- data$DNAm.probe
ICC_ranks_S7 <- sort(ICC_ranks_S7, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S7 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S7,
                           nperm = 10000)

## Scenario 8 ##
# Read in probes #
keep <- grepl("Exon", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep2 <- grepl("Body", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep <- keep | keep2
CodRegProbes <- annEPIC@listData$Name[keep]
CodRegProbes <- sample(CodRegProbes, 35000)
probe_set <- list(probes = CodRegProbes)
ICC_ranks_S8 <- data[,c(1,9)]
ICC_ranks_S8 <- ICC_ranks_S8 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S8 <- as.numeric(ICC_ranks_S8$ICC_Scen8)
names(ICC_ranks_S8) <- data$DNAm.probe
ICC_ranks_S8 <- sort(ICC_ranks_S8, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S8 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S8,
                           nperm = 10000)

## Scenario 9 ##
# Read in probes #
keep <- grepl("Exon", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep2 <- grepl("Body", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep <- keep | keep2
CodRegProbes <- annEPIC@listData$Name[keep]
CodRegProbes <- sample(CodRegProbes, 35000)
probe_set <- list(probes = CodRegProbes)
ICC_ranks_S9 <- data[,c(1,10)]
ICC_ranks_S9 <- ICC_ranks_S9 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S9 <- as.numeric(ICC_ranks_S9$ICC_Scen9)
names(ICC_ranks_S9) <- data$DNAm.probe
ICC_ranks_S9 <- sort(ICC_ranks_S9, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S9 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S9,
                           nperm = 10000)

## Scenario 10 ##
# Read in probes #
keep <- grepl("Exon", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep2 <- grepl("Body", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep <- keep | keep2
CodRegProbes <- annEPIC@listData$Name[keep]
CodRegProbes <- sample(CodRegProbes, 35000)
probe_set <- list(probes = CodRegProbes)
ICC_ranks_S10 <- data[,c(1,11)]
ICC_ranks_S10 <- ICC_ranks_S10 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S10 <- as.numeric(ICC_ranks_S10$ICC_Scen10)
names(ICC_ranks_S10) <- data$DNAm.probe
ICC_ranks_S10 <- sort(ICC_ranks_S10, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S10 <- fgseaSimple(pathways = probe_set, 
                            stats = ICC_ranks_S10,
                            nperm = 10000)

#################
### Undefined ###
#################
### Probe Set Enrichment Analysis ### ----
## Scenario 1 ##
# Read in probes #
NoRegProbes <- annEPIC@listData$Name[annEPIC@listData$UCSC_RefGene_Group == ""]
NoRegProbes <- sample(NoRegProbes, 35000)
probe_set <- list(probes = NoRegProbes)
ICC_ranks_S1 <- data[,c(1,2)]
ICC_ranks_S1 <- ICC_ranks_S1 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S1 <- as.numeric(ICC_ranks_S1$ICC_Scen1)
names(ICC_ranks_S1) <- data$DNAm.probe
ICC_ranks_S1 <- sort(ICC_ranks_S1, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S1 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S1,
                           nperm = 10000)

## Scenario 2 ##
# Read in probes #
NoRegProbes <- annEPIC@listData$Name[annEPIC@listData$UCSC_RefGene_Group == ""]
NoRegProbes <- sample(NoRegProbes, 35000)
probe_set <- list(probes = NoRegProbes)
ICC_ranks_S2 <- data[,c(1,3)]
ICC_ranks_S2 <- ICC_ranks_S2 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S2 <- as.numeric(ICC_ranks_S2$ICC_Scen2)
names(ICC_ranks_S2) <- data$DNAm.probe
ICC_ranks_S2 <- sort(ICC_ranks_S2, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S2 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S2,
                           nperm = 10000)

## Scenario 3 ##
# Read in probes #
NoRegProbes <- annEPIC@listData$Name[annEPIC@listData$UCSC_RefGene_Group == ""]
NoRegProbes <- sample(NoRegProbes, 35000)
probe_set <- list(probes = NoRegProbes)
ICC_ranks_S3 <- data[,c(1,4)]
ICC_ranks_S3 <- ICC_ranks_S3 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S3 <- as.numeric(ICC_ranks_S3$ICC_Scen3)
names(ICC_ranks_S3) <- data$DNAm.probe
ICC_ranks_S3 <- sort(ICC_ranks_S3, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S3 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S3,
                           nperm = 10000)

## Scenario 4 ##
# Read in probes #
NoRegProbes <- annEPIC@listData$Name[annEPIC@listData$UCSC_RefGene_Group == ""]
NoRegProbes <- sample(NoRegProbes, 35000)
probe_set <- list(probes = NoRegProbes)
ICC_ranks_S4 <- data[,c(1,5)]
ICC_ranks_S4 <- ICC_ranks_S4 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S4 <- as.numeric(ICC_ranks_S4$ICC_Scen4)
names(ICC_ranks_S4) <- data$DNAm.probe
ICC_ranks_S4 <- sort(ICC_ranks_S4, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S4 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S4,
                           nperm = 10000)

## Scenario 5 ##
# Read in probes #
NoRegProbes <- annEPIC@listData$Name[annEPIC@listData$UCSC_RefGene_Group == ""]
NoRegProbes <- sample(NoRegProbes, 35000)
probe_set <- list(probes = NoRegProbes)
ICC_ranks_S5 <- data[,c(1,6)]
ICC_ranks_S5 <- ICC_ranks_S5 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S5 <- as.numeric(ICC_ranks_S5$ICC_Scen5)
names(ICC_ranks_S5) <- data$DNAm.probe
ICC_ranks_S5 <- sort(ICC_ranks_S5, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S5 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S5,
                           nperm = 10000)

## Scenario 6 ##
# Read in probes #
NoRegProbes <- annEPIC@listData$Name[annEPIC@listData$UCSC_RefGene_Group == ""]
NoRegProbes <- sample(NoRegProbes, 35000)
probe_set <- list(probes = NoRegProbes)
ICC_ranks_S6 <- data[,c(1,7)]
ICC_ranks_S6 <- ICC_ranks_S6 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S6 <- as.numeric(ICC_ranks_S6$ICC_Scen6)
names(ICC_ranks_S6) <- data$DNAm.probe
ICC_ranks_S6 <- sort(ICC_ranks_S6, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S6 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S6,
                           nperm = 10000)

## Scenario 7 ##
# Read in probes #
NoRegProbes <- annEPIC@listData$Name[annEPIC@listData$UCSC_RefGene_Group == ""]
NoRegProbes <- sample(NoRegProbes, 35000)
probe_set <- list(probes = NoRegProbes)
ICC_ranks_S7 <- data[,c(1,8)]
ICC_ranks_S7 <- ICC_ranks_S7 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S7 <- as.numeric(ICC_ranks_S7$ICC_Scen7)
names(ICC_ranks_S7) <- data$DNAm.probe
ICC_ranks_S7 <- sort(ICC_ranks_S7, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S7 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S7,
                           nperm = 10000)

## Scenario 8 ##
# Read in probes #
NoRegProbes <- annEPIC@listData$Name[annEPIC@listData$UCSC_RefGene_Group == ""]
NoRegProbes <- sample(NoRegProbes, 35000)
probe_set <- list(probes = NoRegProbes)
ICC_ranks_S8 <- data[,c(1,9)]
ICC_ranks_S8 <- ICC_ranks_S8 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S8 <- as.numeric(ICC_ranks_S8$ICC_Scen8)
names(ICC_ranks_S8) <- data$DNAm.probe
ICC_ranks_S8 <- sort(ICC_ranks_S8, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S8 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S8,
                           nperm = 10000)

## Scenario 9 ##
# Read in probes #
NoRegProbes <- annEPIC@listData$Name[annEPIC@listData$UCSC_RefGene_Group == ""]
NoRegProbes <- sample(NoRegProbes, 35000)
probe_set <- list(probes = NoRegProbes)
ICC_ranks_S9 <- data[,c(1,10)]
ICC_ranks_S9 <- ICC_ranks_S9 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S9 <- as.numeric(ICC_ranks_S9$ICC_Scen9)
names(ICC_ranks_S9) <- data$DNAm.probe
ICC_ranks_S9 <- sort(ICC_ranks_S9, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S9 <- fgseaSimple(pathways = probe_set, 
                           stats = ICC_ranks_S9,
                           nperm = 10000)

## Scenario 10 ##
# Read in probes #
NoRegProbes <- annEPIC@listData$Name[annEPIC@listData$UCSC_RefGene_Group == ""]
NoRegProbes <- sample(NoRegProbes, 35000)
probe_set <- list(probes = NoRegProbes)
ICC_ranks_S10 <- data[,c(1,11)]
ICC_ranks_S10 <- ICC_ranks_S10 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S10 <- as.numeric(ICC_ranks_S10$ICC_Scen10)
names(ICC_ranks_S10) <- data$DNAm.probe
ICC_ranks_S10 <- sort(ICC_ranks_S10, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S10 <- fgseaSimple(pathways = probe_set, 
                            stats = ICC_ranks_S10,
                            nperm = 10000)

#############################
### Circadian Clock Gene  ###
#############################
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ProbeData <- annEPIC@listData
ProbeData_Genes <- annEPIC@listData$UCSC_RefGene_Name
ProbeData_Genes <-gsub(";.*","",ProbeData_Genes)
Circ_Genes <- c("CLOCK", "NPAS2", "BMAL1", "BMAL2", "PER1",
                "PER2", "PER3", "CRY1", "CRY2")
keep <- ProbeData_Genes %in% Circ_Genes
ProbeDataNames <- ProbeData$Name[keep]

## Scenario 1 ##
# Read in probes #
probe_set <- list(probes = ProbeDataNames)
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
# Read in probes #
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
# Read in probes #
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
# Read in probes #
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
# Read in probes #
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
# Read in probes #
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
# Read in probes #
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
# Read in probes #
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
# Read in probes #
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
# Read in probes #
ICC_ranks_S10 <- data[,c(1,11)]
ICC_ranks_S10 <- ICC_ranks_S10 %>% remove_rownames %>% column_to_rownames(var="DNAm.probe")
ICC_ranks_S10 <- as.numeric(ICC_ranks_S10$ICC_Scen10)
names(ICC_ranks_S10) <- data$DNAm.probe
ICC_ranks_S10 <- sort(ICC_ranks_S10, decreasing=TRUE)
# Running probe set enrichment analysis #
fgseaRes_S10 <- fgseaSimple(pathways = mQTL_probes, 
                            stats = ICC_ranks_S10,
                            nperm = 10000)


#############################
###### Creating Figure ######
#############################
### Creating Histogram and Boxplots Combined ### ----
### Load Annotation Data ### ---
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

### Create Boxplots for Each Scenario ### ----
probe_set <- read.csv("probe_set.csv")
data_All <- data[data$DNAm.probe %in% probe_set$DNAm.probe,]
data_All$Name <- as.factor("mQTL")

keep <- grepl("TSS", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
TSS_Probes <- annEPIC@listData$Name[keep]
data_TSS <- data[data$DNAm.probe %in% TSS_Probes,]
data_TSS$Name <- as.factor("TSS")
data_All <- rbind(data_All, data_TSS)

keep <- grepl("5'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Five_UTR_Probes <- annEPIC@listData$Name[keep]
data_FiveUTR <- data[data$DNAm.probe %in% Five_UTR_Probes,]
data_FiveUTR$Name <- as.factor("5'UTR")
data_All <- rbind(data_All, data_FiveUTR)

keep <- grepl("3'UTR", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
Three_UTR_Probes <- annEPIC@listData$Name[keep]
data_ThreeUTR <- data[data$DNAm.probe %in% Three_UTR_Probes,]
data_ThreeUTR$Name <- as.factor("3'UTR")
data_All <- rbind(data_All, data_ThreeUTR)

keep <- grepl("Exon", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep2 <- grepl("Body", annEPIC@listData$UCSC_RefGene_Group, fixed=TRUE)
keep <- keep | keep2
CodRegProbes <- annEPIC@listData$Name[keep]
data_Cod <- data[data$DNAm.probe %in% CodRegProbes,]
data_Cod$Name <- as.factor("Coding Region")
data_All <- rbind(data_All, data_Cod)

NoRegProbes <- annEPIC@listData$Name[annEPIC@listData$UCSC_RefGene_Group == ""]
data_NOR <- data[data$DNAm.probe %in% NoRegProbes,]
data_NOR$Name <- as.factor("Undefined")
data_All <- rbind(data_All, data_NOR)

ProbeData <- annEPIC@listData
ProbeData_Genes <- annEPIC@listData$UCSC_RefGene_Name
ProbeData_Genes <-gsub(";.*","",ProbeData_Genes)
Circ_Genes <- c("CLOCK", "NPAS2", "BMAL1", "BMAL2", "PER1",
                "PER2", "PER3", "CRY1", "CRY2")
keep <- ProbeData_Genes %in% Circ_Genes
CLOCK_Probes <- ProbeData$Name[keep]
data_CLOCK <- data[data$DNAm.probe %in% CLOCK_Probes,]
data_CLOCK$Name <- as.factor("Circadian")
data_All <- rbind(data_All, data_CLOCK)

# Exclude probes with ICC < 0.01 #
data_All[,2:11][data_All[,2:11] < 0.01] <- NA

Box1 <- ggplot(data_All, aes(x=Name, y=ICC_Scen1, fill=Name)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0.47,linetype=2) +
  labs(title="NoStressT1-2", x =NULL, y = NULL) + 
  scale_fill_brewer(palette="Set2") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
Box2 <- ggplot(data_All, aes(x=Name, y=ICC_Scen2, fill=Name)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0.45,linetype=2) +
  labs(title="NoStressT1-3", x =NULL, y = NULL) + 
  scale_fill_brewer(palette="Set2") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
Box3 <- ggplot(data_All, aes(x=Name, y=ICC_Scen3, fill=Name)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0.31,linetype=2) +
  labs(title="NoStressT1-4", x =NULL, y = NULL) + 
  scale_fill_brewer(palette="Set2") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
Box4 <- ggplot(data_All, aes(x=Name, y=ICC_Scen4, fill=Name)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0.50,linetype=2) +
  labs(title="NoStressT3-4", x =NULL, y = NULL) + 
  scale_fill_brewer(palette="Set2") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
Box5 <- ggplot(data_All, aes(x=Name, y=ICC_Scen5, fill=Name)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0.44,linetype=2) +
  labs(title="NoStressT1-2-3-4", x =NULL, y = NULL) + 
  scale_fill_brewer(palette="Set2") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
Box6 <- ggplot(data_All, aes(x=Name, y=ICC_Scen6samp, fill=Name)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0.64,linetype=2) +
  labs(title="StressT1-2", x =NULL, y = NULL) + 
  scale_fill_brewer(palette="Set2") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
Box7 <- ggplot(data_All, aes(x=Name, y=ICC_Scen7samp, fill=Name)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0.64,linetype=2) +
  labs(title="StressT1-3", x =NULL, y = NULL) + 
  scale_fill_brewer(palette="Set2") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
Box8 <- ggplot(data_All, aes(x=Name, y=ICC_Scen8samp, fill=Name)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0.65,linetype=2) +
  labs(title="StressT1-4", x =NULL, y = NULL) + 
  scale_fill_brewer(palette="Set2") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
Box9 <- ggplot(data_All, aes(x=Name, y=ICC_Scen9samp, fill=Name)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0.62,linetype=2) +
  labs(title="StressT3-4", x =NULL, y = NULL) + 
  scale_fill_brewer(palette="Set2") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
Box10 <- ggplot(data_All, aes(x=Name, y=ICC_Scen10, fill=Name)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0.28,linetype=2) +
  labs(title="CrossSessionT1", x =NULL, y = NULL) + 
  scale_fill_brewer(palette="Set2") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

### Combine All Plots ### ----
ggarrange(Box1, Box2, Box3, Box4, Box5,
          Box6, Box7, Box8, Box9, Box10,
          ncol = 5, nrow = 2, legend = "right", common.legend = TRUE)


