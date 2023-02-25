### Aim 15 ### ----

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

############################################################################
################################## Aim 15 ##################################
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

### Filter to Only Keep Probes Highly Reliable (ICC > 0.9), in Scenarios 1-5 and 10 ### ----
HighProbes <- dplyr::filter(data, data$ICC_Scen1 >= 0.9 &
                              data$ICC_Scen2 >= 0.9 &
                              data$ICC_Scen3 >= 0.9 &
                              data$ICC_Scen4 >= 0.9 &
                              data$ICC_Scen5 >= 0.9 &
                              data$ICC_Scen10 >= 0.9)
Stress_S6_HighProbes <- dplyr::filter(HighProbes, HighProbes$ICC_Scen6 <= 1)
Stress_S7_HighProbes <- dplyr::filter(HighProbes, HighProbes$ICC_Scen7 <= 1)
Stress_S8_HighProbes <- dplyr::filter(HighProbes, HighProbes$ICC_Scen8 <= 1)

### Determine Genes Associated with Probes ### ----
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
keep <- grepl("TSS", annEPIC$UCSC_RefGene_Group, fixed=TRUE)
annEPIC <- annEPIC[keep,]

keep <- annEPIC$Name %in% Stress_S6_HighProbes$DNAm.probe
annEPIC_S6 <- annEPIC[keep,]
Stress_Genes_S6 <- annEPIC_S6$UCSC_RefGene_Name
Stress_Genes_S6 <- gsub(";.*","",Stress_Genes_S6)

keep <- annEPIC$Name %in% Stress_S7_HighProbes$DNAm.probe
annEPIC_S7 <- annEPIC[keep,]
Stress_Genes_S7 <- annEPIC_S7$UCSC_RefGene_Name
Stress_Genes_S7 <- gsub(";.*","",Stress_Genes_S7)

keep <- annEPIC$Name %in% Stress_S8_HighProbes$DNAm.probe
annEPIC_S8 <- annEPIC[keep,]
Stress_Genes_S8 <- annEPIC_S8$UCSC_RefGene_Name
Stress_Genes_S8 <- gsub(";.*","",Stress_Genes_S8)

### Keep Only Probes at TSS ### ----
keep <- Stress_S6_HighProbes$DNAm.probe %in% annEPIC@listData$Name
Stress_S6_HighProbes <- Stress_S6_HighProbes[keep,]

keep <- Stress_S7_HighProbes$DNAm.probe %in% annEPIC@listData$Name
Stress_S7_HighProbes <- Stress_S7_HighProbes[keep,]

keep <- Stress_S8_HighProbes$DNAm.probe %in% annEPIC@listData$Name
Stress_S8_HighProbes <- Stress_S8_HighProbes[keep,]

### Read in Data ### ----
Stress_Genes_S6 <- mapIds(org.Hs.eg.db, Stress_Genes_S6, 'ENTREZID', 'SYMBOL')
Stress_Genes_S7 <- mapIds(org.Hs.eg.db, Stress_Genes_S7, 'ENTREZID', 'SYMBOL')
Stress_Genes_S8 <- mapIds(org.Hs.eg.db, Stress_Genes_S8, 'ENTREZID', 'SYMBOL')

### Gene over-representation analysis ### ----
go_enrich_S6 <- enrichGO(gene = Stress_Genes_S6,
                            OrgDb = org.Hs.eg.db, 
                            keyType = 'ENTREZID',
                            readable = T,
                            ont = "BP",
                            pvalueCutoff = 0.05, 
                            qvalueCutoff = 0.10)

go_enrich_S7 <- enrichGO(gene = Stress_Genes_S7,
                         OrgDb = org.Hs.eg.db, 
                         keyType = 'ENTREZID',
                         readable = T,
                         ont = "BP",
                         pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.10)

go_enrich_S8 <- enrichGO(gene = Stress_Genes_S8,
                         OrgDb = org.Hs.eg.db, 
                         keyType = 'ENTREZID',
                         readable = T,
                         ont = "BP",
                         pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.10)

# Create barplot for top categories #
barplot(go_enrich_S8, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Biological Pathways - Stress",
        font.size = 8)

### Scatter Plots for Specific Probes ### ----
# Read in data #
bVals <- readRDS("../Methylation/bVals_NoobNorm_SampFilter.RDS")
Phenotypes <- read.csv("Phenotypes.csv")[,-1]
# Keep only probes that are unreliable in all three stress sessions #
Probes_to_Test <- unique(c(Stress_S6_HighProbes$DNAm.probe, 
                           Stress_S7_HighProbes$DNAm.probe,
                           Stress_S8_HighProbes$DNAm.probe))
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
MLM_Results <- data.frame(DNAm.probe = colnames(Phenotypes)[30:2218],
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
  MLM_Results[i,10:13] <- MLM_Results[i,6:9] * nrow(Stress_S6_HighProbes)
}
# Add FDR adjusted p-values #
MLM_Results$Pval_Int_FDRadj <- p.adjust(MLM_Results$Pval_Int, method = "fdr")
MLM_Results$Pval_T2_FDRadj <- p.adjust(MLM_Results$Pval_T2, method = "fdr")
MLM_Results$Pval_T3_FDRadj <- p.adjust(MLM_Results$Pval_T3, method = "fdr")
MLM_Results$Pval_T4_FDRadj <- p.adjust(MLM_Results$Pval_T4, method = "fdr")

### Create Manhattan Plot ### ----
# Create data frame for plotting #
MLM_Results_Sorted <- MLM_Results[order(MLM_Results$DNAm.probe),]
Probe_Additional_Info <- data.frame(DNAm.probe = annEPIC_S6@listData$Name,
                                    Chr = annEPIC_S6@listData$chr,
                                    Pos = annEPIC_S6@listData$pos,
                                    Gene = annEPIC_S6@listData$UCSC_RefGene_Name)
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
mhtplot(probe=MLM_Results_Sorted$DNAm.probe,
        chr=MLM_Results_Sorted$Chr,
        pos=MLM_Results_Sorted$Pos,
        p=MLM_Results_Sorted$Pval_T4,
        sigthre=0.05,
        sigthre2=(0.05/1413),
        outf="T4vT1_ManhattanPlot",
        #markprobe = c("cg19913448"),
        outfmt="jpg")

### Gene Set Overrepresentation Analysis ### ----
# Determine which genes to keep #
GSO_Genes <- MLM_Results_Sorted$Gene[MLM_Results_Sorted$Pval_T4 < (0.05/1413)]
GSO_Genes <- gsub(";.*","",GSO_Genes)
# Change gene names to ENTREZID #
GSO_Genes <- mapIds(org.Hs.eg.db, GSO_Genes, 'ENTREZID', 'SYMBOL')
# Gene over-representation analysis #
MLM_Results_GO <- enrichGO(gene = GSO_Genes,
                         OrgDb = org.Hs.eg.db, 
                         keyType = 'ENTREZID',
                         readable = T,
                         ont = "BP",
                         pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.10)

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
