### Results - Section 2.4 ### ----

### Load Packages ### ----
library(dplyr)
library(psych)
library(ggplot2)
library(ggpubr)
library(dgof)
library(fgsea)
library(tidyverse)
library('org.Hs.eg.db')
library(pathview)
library(clusterProfiler)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(tibble)
library(lme4)
library(lmerTest)
library(ENmix)

##################################################################################################
################################## Highly Stable Probe Analysis ##################################
##################################################################################################
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

### Filter to Only Keep Probes Highly Reliable (ICC > 0.9), in Scenarios 1-5 and 10 ### ----
HighProbes <- dplyr::filter(data, data$ICC_Scen1 >= 0.9 &
                              data$ICC_Scen2 >= 0.9 &
                              data$ICC_Scen3 >= 0.9 &
                              data$ICC_Scen4 >= 0.9 &
                              data$ICC_Scen5 >= 0.9 &
                              data$ICC_Scen10 >= 0.9)

### Determine Genes Associated with Probes ### ----
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
keep <- grepl("TSS", annEPIC$UCSC_RefGene_Group, fixed=TRUE)
annEPIC <- annEPIC[keep,]

keep <- annEPIC$Name %in% HighProbes$DNAm.probe
annEPIC <- annEPIC[keep,]
Stress_Genes <- annEPIC$UCSC_RefGene_Name
Stress_Genes <- gsub(";.*","",Stress_Genes)

### Keep Only Probes at TSS ### ----
keep <- HighProbes$DNAm.probe %in% annEPIC@listData$Name
HighProbes <- HighProbes[keep,]

### Scatter Plots for Specific Probes ### ----
# Read in data #
bVals <- readRDS("../Methylation/bVals_NoobNorm_SampFilter.RDS")
Phenotypes <- read.csv("Phenotypes.csv")[,-1]
# Keep only probes that are unreliable in all three stress sessions #
Probes_to_Test <- bvals[,1] %in% HighProbes$DNAm.probe
bVals_Test <- as.data.frame(bVals[Probes_to_Test,])
bVals_Test <- data.frame(t(bVals_Test))
bVals_Test <- tibble::rownames_to_column(bVals_Test, "Sample_Name")
# Convert beta to M values #
bVals_Test[,-1] <- log2(bVals_Test[,-1]/(1 - bVals_Test[,-1]))
# Add CpG value to phenotype file #
Phenotypes <- merge(Phenotypes, bVals_Test, by="Sample_Name")
# Keep only stress session values #
Phenotypes <- dplyr::filter(Phenotypes, Phenotypes$Session == 1)

### Multilevel Modeling ### ----
# Define empty dataframe to store results #
MLM_Results <- data.frame(DNAm.probe = colnames(Phenotypes)[30:8195],
                          Intercept = NA,
                          Time2_Beta = NA,
                          Time3_Beta = NA,
                          Time4_Beta = NA,
                          Pval_Int = NA,
                          Pval_T2 = NA,
                          Pval_T3 = NA,
                          Pval_T4 = NA,
                          Pval_Int_BFadj = NA,
                          Pval_T2_BFadj = NA,
                          Pval_T3_BFadj = NA,
                          Pval_T4_BFadj = NA,
                          Pval_Int_FDRadj = NA,
                          Pval_T2_FDRadj = NA,
                          Pval_T3_FDRadj = NA,
                          Pval_T4_FDRadj = NA)
# Run multilevel models #
for (j in 30:ncol(Phenotypes)){
  i = j - 29
  mean_CpG <- mean(Phenotypes[,j])
  sd_CpG <- sd(Phenotypes[,j])
  Cutoff_high <- mean_CpG + (3 * sd_CpG)
  Cutoff_low <- mean_CpG - (3 * sd_CpG)
  for (c in 1:nrow(Phenotypes)){
    if(Phenotypes[c,j] > Cutoff_high | Phenotypes[c,j] < Cutoff_low){Phenotypes[c,j] <- NA}
    #if(Phenotypes[c,j] < Cutoff_low){Phenotypes[c,j] <- NA}
  }
  MLM <- lmerTest::lmer(Phenotypes[,j] ~ 1 + as.factor(Time) + DNAm_Monocytes + as.factor(Methylation_Batch) + (1|Individual), Phenotypes)
  MLM_Results[i,2:5] <- summary(MLM)$coefficients[1:4,1]
  MLM_Results[i,6:9] <- summary(MLM)$coefficients[1:4,5]
  MLM_Results[i,10:13] <- MLM_Results[i,6:9] * nrow(HighProbes)
}
# Add FDR adjusted p-values #
MLM_Results$Pval_Int_FDRadj <- p.adjust(MLM_Results$Pval_Int, method = "fdr")
MLM_Results$Pval_T2_FDRadj <- p.adjust(MLM_Results$Pval_T2, method = "fdr")
MLM_Results$Pval_T3_FDRadj <- p.adjust(MLM_Results$Pval_T3, method = "fdr")
MLM_Results$Pval_T4_FDRadj <- p.adjust(MLM_Results$Pval_T4, method = "fdr")

### Create Manhattan Plot ### ----
# Create data frame for plotting #
MLM_Results_Sorted <- MLM_Results[order(MLM_Results$DNAm.probe),]
Probe_Additional_Info <- data.frame(DNAm.probe = annEPIC@listData$Name,
                                    Chr = annEPIC@listData$chr,
                                    Pos = annEPIC@listData$pos,
                                    Gene = annEPIC@listData$UCSC_RefGene_Name)
# Combine additional info to MLM results #
MLM_Results_Sorted <- merge(MLM_Results_Sorted, Probe_Additional_Info, by = "DNAm.probe")
# Exclude sex chromosome probes #
MLM_Results_Sorted <- dplyr::filter(MLM_Results_Sorted, MLM_Results_Sorted$Chr != "chrX")
MLM_Results_Sorted <- dplyr::filter(MLM_Results_Sorted, MLM_Results_Sorted$Chr != "chrY")
# Save MLM_Results_Sorted #
saveRDS(MLM_Results_Sorted, "MLM_Results_Sorted_Aim15.Rds")
MLM_Results_Sorted <- readRDS("MLM_Results_Sorted_Aim15.Rds")
# Generate Manhattan plots #
# BF sig from sex chromosome excluded #
# Time 2 #
mhtplot(probe=MLM_Results_Sorted$DNAm.probe,
        chr=MLM_Results_Sorted$Chr,
        pos=MLM_Results_Sorted$Pos,
        p=MLM_Results_Sorted$Pval_T2,
        sigthre=0.05,
        sigthre2=(0.05/1293),
        outf="T2vT1_ManhattanPlot",
        #markprobe = c("cg19913448"),
        outfmt="jpg")

# Time 3 #
mhtplot(probe=MLM_Results_Sorted$DNAm.probe,
        chr=MLM_Results_Sorted$Chr,
        pos=MLM_Results_Sorted$Pos,
        p=MLM_Results_Sorted$Pval_T3,
        sigthre=0.05,
        sigthre2=(0.05/1293),
        outf="T3vT1_ManhattanPlot",
        #markprobe = c("cg19913448"),
        outfmt="jpg")

# Time 4 #
mhtplot(probe=MLM_Results_Sorted$DNAm.probe,
        chr=MLM_Results_Sorted$Chr,
        pos=MLM_Results_Sorted$Pos,
        p=MLM_Results_Sorted$Pval_T4,
        sigthre=0.05,
        sigthre2=(0.05/1293),
        outf="T4vT1_ManhattanPlot",
        #markprobe = c("cg19913448"),
        outfmt="jpg")

### Looking at GSR Expression Profiles ### ----
Phenotypes <- read.csv("Phenotypes.csv")[,-1]
GSR_ExpData <- readRDS("GSR_ExpData.Rds")
GSR_ExpData <- t(GSR_ExpData)[-1,]
GSR_ExpData <- data.frame(GSR_Exp = GSR_ExpData)
GSR_ExpData <- tibble::rownames_to_column(GSR_ExpData, "Sample_Name")
GSR_ExpData$GSR_Exp <- as.numeric(GSR_ExpData$GSR_Exp)
Phenotypes <- merge(Phenotypes, GSR_ExpData, by = "Sample_Name")
Phenotypes <- dplyr::filter(Phenotypes, Phenotypes$Session == 1)

GSR_MLM <- lmerTest::lmer(GSR_Exp ~ 1 + as.factor(Time)*ELA + DNAm_Monocytes + (1|Individual), Phenotypes)
summary(GSR_MLM)
