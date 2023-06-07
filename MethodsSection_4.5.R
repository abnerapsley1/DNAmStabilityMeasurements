### Methods - Section 4.5 ### ----
### Adjust ICC for factors using MLM ###----

# Load packages # ----
library(tidyverse)
library(lme4)
library(psych) # ICC() function
library(future) # Parallel processing
plan(multicore)
options(future.globals.maxSize = 1610612736) # Need 9 cores
set.seed(666666)

# Read in the probe data and phenotype data # ----
load("Data/bVals_SampFilt.RData")

bVals.t <- t(bVals)
probe.name <- colnames(bVals.t)
probe.n <- length(probe.name)
targets$Methylation_Batch <- factor(targets$Methylation_Batch)
targets$Sex <- factor(targets$Sex)
targets$MinorityStatus <- factor(targets$MinorityStatus)

# Create new R functions # ----
# Update the function for ICC calculation to be able to control for covariates when calculating ICCs 
# (adapted from the code for the ICC() function)
ICC.n <- function(data, id, item, value, covariates = NULL, alpha = .05, lmer = TRUE) 
{
  cl <- match.call()
  nj <- length(unique(data[, item]))
  n.obs <- length(unique(data[, id]))
  if (lmer) {
    if (is.null(covariates)){
      formula <- as.formula(paste(value, "~ 1 + (1|", id, ") + (1 |", item, ")"))
    }
    else {
      formula <- as.formula(paste(value, "~ 1 + (1|", id, ") + (1 |", item, ")", 
                                  "+", paste0(covariates, collapse = "+")))
    }
    suppressMessages(mod.lmer <- lme4::lmer(formula, 
                                            data = data, na.action = na.omit))
    vc <- lme4::VarCorr(mod.lmer)
    MS_id <- vc[[1]][1,1]
    MS_items <- vc[[2]][1,1]
    MSE <- error <- MS_resid <- (attributes(vc)$sc)^2
    MS.df <- data.frame(variance = c(MS_id, MS_items, MS_resid, 
                                     NA))
    rownames(MS.df) <- c("ID", "Items", "Residual", "Total")
    MS.df["Total", ] <- sum(MS.df[1:3, 1], na.rm = TRUE)
    MS.df["Percent"] <- MS.df/MS.df["Total", 1]
    lmer.MS <- MS.df
    MSB <- nj * MS_id + error
    MSJ <- n.obs * MS_items + error
    MSW <- error + MS_items
    stats <- matrix(NA, ncol = 3, nrow = 5)
    dfB <- n.obs - 1
    dfJ <- nj - 1
    dfE <- (n.obs - 1) * (nj - 1)
    FB <- MSB/MSE
    FJ <- MSJ/MSE
    stats[1, ] <- c(dfB, dfJ, dfE)
    stats[2, ] <- c(MSB * (n.obs - 1), MSJ * (nj - 1), MSE * 
                      (n.obs - 1) * (nj - 1))
    stats[3, ] <- c(MSB, MSJ, MSE)
    stats[4, ] <- c(MSB/MSE, MSJ/MSE, NA)
    stats[5, ] <- c(-expm1(pf(FB, dfB, dfE, log.p = TRUE)), 
                    -expm1(pf(FJ, dfJ, dfE, log.p = TRUE)), NA)
    s.aov <- mod.lmer
  }
  else {
    data[, id] <- as.factor(data[, id])
    data[, item] <- as.factor(data[, item])
    if (is.null(covariates)){
      formula <- as.formula(paste(value, "~", id, "+", item))
      aov.x <- aov(formula, data = data)
      s.aov <- summary(aov.x)
      stats <- matrix(unlist(s.aov), ncol = 3, byrow = TRUE)
    }
    else {
      formula <- as.formula(paste(value, "~", paste0(covariates, collapse = "+"), "+", id, "+", item))
      aov.x <- aov(formula, data = data)
      s.aov <- summary(aov.x)
      stats <- matrix(unlist(s.aov), ncol = 3 + length(covariates), byrow = TRUE)[,-1]
    }
    MSB <- stats[3, 1]
    MSW <- (stats[2, 2] + stats[2, 3])/(stats[1, 2] + stats[1,3])
    MSJ <- stats[3, 2]
    MSE <- stats[3, 3]
    MS.df <- NULL
  }
  colnames(stats) <- c("subjects", "Judges", "Residual")
  rownames(stats) <- c("df", "SumSq", "MS", "F", "p")
  ICC1 <- (MSB - MSW)/(MSB + (nj - 1) * MSW)
  ICC2 <- (MSB - MSE)/(MSB + (nj - 1) * MSE + nj * (MSJ - MSE)/n.obs)
  ICC3 <- (MSB - MSE)/(MSB + (nj - 1) * MSE)
  ICC12 <- (MSB - MSW)/(MSB)
  ICC22 <- (MSB - MSE)/(MSB + (MSJ - MSE)/n.obs)
  ICC32 <- (MSB - MSE)/MSB
  F11 <- MSB/MSW
  df11n <- n.obs - 1
  df11d <- n.obs * (nj - 1)
  p11 <- -expm1(pf(F11, df11n, df11d, log.p = TRUE))
  F21 <- MSB/MSE
  df21n <- n.obs - 1
  df21d <- (n.obs - 1) * (nj - 1)
  p21 <- -expm1(pf(F21, df21n, df21d, log.p = TRUE))
  F31 <- F21
  results <- data.frame(matrix(NA, ncol = 8, nrow = 6))
  colnames(results) <- c("type", "ICC", "F", "df1", "df2", 
                         "p", "lower bound", "upper bound")
  rownames(results) <- c("Single_raters_absolute", "Single_random_raters", 
                         "Single_fixed_raters", "Average_raters_absolute", "Average_random_raters", 
                         "Average_fixed_raters")
  results[, 1] <- c("ICC1", "ICC2", "ICC3", "ICC1k", "ICC2k", 
                    "ICC3k")
  results[, 2] <- c(ICC1, ICC2, ICC3, ICC12, ICC22, ICC32)
  results[, 3] <- c(F11, F21, F21, F11, F21, F21)
  results[, 4] <- c(df11n, df21n, df21n, df11n, df21n, df21n)
  results[, 5] <- c(df11d, df21d, df21d, df11d, df21d, df21d)
  results[, 6] <- c(p11, p21, p21, p11, p21, p21)
  F1L <- F11/qf(1 - alpha/2, df11n, df11d)
  F1U <- F11 * qf(1 - alpha/2, df11d, df11n)
  L1 <- (F1L - 1)/(F1L + (nj - 1))
  U1 <- (F1U - 1)/(F1U + nj - 1)
  F3L <- F31/qf(1 - alpha/2, df21n, df21d)
  F3U <- F31 * qf(1 - alpha/2, df21d, df21n)
  results[1, 7] <- L1
  results[1, 8] <- U1
  results[3, 7] <- (F3L - 1)/(F3L + nj - 1)
  results[3, 8] <- (F3U - 1)/(F3U + nj - 1)
  results[4, 7] <- 1 - 1/F1L
  results[4, 8] <- 1 - 1/F1U
  results[6, 7] <- 1 - 1/F3L
  results[6, 8] <- 1 - 1/F3U
  Fj <- MSJ/MSE
  vn <- (nj - 1) * (n.obs - 1) * ((nj * ICC2 * Fj + n.obs * 
                                     (1 + (nj - 1) * ICC2) - nj * ICC2))^2
  vd <- (n.obs - 1) * nj^2 * ICC2^2 * Fj^2 + (n.obs * (1 + 
                                                         (nj - 1) * ICC2) - nj * ICC2)^2
  v <- vn/vd
  F3U <- qf(1 - alpha/2, n.obs - 1, v)
  F3L <- qf(1 - alpha/2, v, n.obs - 1)
  L3 <- n.obs * (MSB - F3U * MSE)/(F3U * (nj * MSJ + (nj * 
                                                        n.obs - nj - n.obs) * MSE) + n.obs * MSB)
  results[2, 7] <- L3
  U3 <- n.obs * (F3L * MSB - MSE)/(nj * MSJ + (nj * n.obs - 
                                                 nj - n.obs) * MSE + n.obs * F3L * MSB)
  results[2, 8] <- U3
  L3k <- L3 * nj/(1 + L3 * (nj - 1))
  U3k <- U3 * nj/(1 + U3 * (nj - 1))
  results[5, 7] <- L3k
  results[5, 8] <- U3k
  result <- list(results = results, summary = s.aov, stats = stats, 
                 MSW = MSW, lme = MS.df, Call = cl, n.obs = n.obs, n.judge = nj)
  class(result) <- c("psych", "ICC")
  return(result)
}

# Create a function outputting the wanted types of ICC for all the probes
ICC.n.m <- function(probe.data, phenotype.data, subject, time, probe, covariates = NULL, 
                    alpha = .05, lmer = TRUE, 
                    icc.type = c("ICC1", "ICC2", "ICC3", 
                                 "ICC1k", "ICC2k", "ICC3k"), sample.individual = F) {
  res <- data.frame(NULL)
  for(i in 1:length(probe)) {
    
    if (sample.individual == T) {
      samp.ind <- sample(unique(phenotype.data$Individual), size = 14, replace = F)
      phenotype.data.n <- phenotype.data %>%
        filter(Individual %in% samp.ind)
    }
    else {phenotype.data.n <- phenotype.data}
    
    probe.data.n <- probe.data[match(phenotype.data.n$Sample_Name, rownames(probe.data)), i]
    
    data.comb <- cbind(probe.data.n, phenotype.data.n)
    colnames(data.comb)[1] <- probe[i]
    
    covariates2 <- c(NULL)
    for (k in 1:length(covariates)) {
      if (class(data.comb[, covariates[k]]) == "numeric") {
        covariates2 <- c(covariates2, covariates[k])
      } else if (class(data.comb[, covariates[k]]) == "factor") {
        if (length(unique(data.comb[, covariates[k]])) == length(levels(data.comb[, covariates[k]]))) {
          covariates2 <- c(covariates2, covariates[k])
        }
        else {covariates2 <- covariates2}
      }
    }
    res.temp <- ICC.n(data = data.comb, id = subject, item = time, value = probe[i], 
                      covariates = covariates2, alpha = alpha, lmer = lmer)$results
    for (j in 1:length(icc.type)){
      res.ind <- which(res.temp$type == icc.type)
      res.temp.2 <- data.frame(DNAm.probe = probe[i],
                               covariates = ifelse(is.null(covariates), "No covariates", 
                                                   paste0(covariates, collapse = "_")),
                               type.v = rownames(res.temp)[res.ind],
                               res.temp[res.ind, ])
      rownames(res.temp.2) <- NULL
      res <- rbind(res, res.temp.2)

    }
    if (i %% 10000 == 0){
      print(i)
    }
  }
  return(res)
}

# Bland-Altman plot function
### code was adapted from https://cegh.net/article/S2213-3984(21)00139-1/fulltext
BA.plot <- function(m1, m2, note) {
  means <- (m1 + m2) / 2
  diffs <- m1 - m2; mdiff <- mean(diffs); sddiff <- sd(diffs)
  # Positioning of labels on the lines
  xmax=max(means)
  xmin=min(means)
  x1=(xmax-xmin)*0.10
  x2=xmax-x1
  # Adding space for legend
  ylimh <- mdiff + (6* sddiff)
  yliml <- mdiff - (6.5 *sddiff)
  # Plot data
  plot(diffs ~ means, xlab = "Average of two methods", ylab = "Difference of two methods", ylim = c(yliml, ylimh))
  abline(h = mdiff, col='blue', lty = 2)
  text(x2, mdiff, "Mean")
  abline(h = 0)
  # Standard deviations lines & legend
  abline(h = mdiff + 1.96 * sddiff, col='blue', lty = 2)
  text(x2, mdiff + 1.96 * sddiff, "Upper Limit")
  abline(h = mdiff - 1.96 * sddiff, col='blue', lty = 2)
  text(x2, mdiff - 1.96 * sddiff, "Lower Limit")
  paramX =round(mdiff,2); paramY =round(sddiff,2)
  UL=mdiff + 1.96 * sddiff; UL=round(UL,2)
  LL=mdiff - 1.96 * sddiff; LL=round(LL,2)
  expr <- vector("expression", 4)
  expr[[1]] <- bquote(Mean[diff]==.(paramX))
  expr[[2]] <- bquote(SD[diff]==.(paramY))
  expr[[3]] <- bquote(LL==.(LL))
  expr[[4]] <- bquote(UL==.(UL))
  legend("topright", bty = "n",legend = expr)
  title(main = paste("Bland-Altman plot \n", note, sep = ""))}

# ICC calculation # ----
revision <- T
scen <- opt$scenario
if (scen %in% c(1:5, 11)) {
  sess <- 0
  samp <- F
  ela <- c(0, 1)
  tors <- "Time"
  covar <- c("DNAm_Monocytes", "Methylation_Batch")
  comp.list <- case_when(scen == 1 ~ list(c(1, 2)),
                         scen == 2 ~ list(c(2, 3)),
                         scen == 3 ~ list(c(1, 3)),
                         scen == 4 ~ list(c(3, 4)),
                         scen == 5 ~ list(c(1, 2, 3, 4)),
                         scen == 11 ~ list(c(1, 4)))
} else if (scen %in% c(6:9, 12)) {
  sess <- 1
  samp <- T
  ela <- c(0, 1)
  tors <- "Time"
  covar <- c("DNAm_Monocytes", "Methylation_Batch")
  comp.list <- case_when(scen == 6 ~ list(c(1, 2)),
                         scen == 7 ~ list(c(1, 3)),
                         scen == 8 ~ list(c(1, 4)),
                         scen == 9 ~ list(c(3, 4)),
                         scen == 12 ~ list(c(1, 2, 3, 4))) 
  
} else if (scen %in% paste0(rep(c(6:9), each = 2), c("a", "b"))) {
  sess <- 1
  samp <- F
  tors <- "Time"
  covar <- c("DNAm_Monocytes", "Methylation_Batch", "Age", "Sex", "MinorityStatus")
  ela <- ifelse(str_detect(scen, pattern = "a") == T, 0, 1)
  comp.list <- case_when(scen %in% c("6a", "6b") ~ list(c(1, 2)),
                         scen %in% c("7a", "7b") ~ list(c(1, 3)),
                         scen %in% c("8a", "8b") ~ list(c(1, 4)),
                         scen %in% c("9a", "9b") ~ list(c(3, 4)))
} else if (scen == 10) {
  samp <- F
  sess <- c(0, 1)
  comp.list <- list(c(1))
  ela <- c(0, 1)
  tors <- "Session"
  covar <- c("DNAm_Monocytes", "Methylation_Batch")
}


ind <- (targets %>%
          filter(Session %in% sess) %>%
          filter(Time %in% comp.list[[1]]) %>%
          filter(ELA %in% ela) %>%
          group_by(Individual) %>%
          dplyr::count() %>%
          filter(n == ifelse(length(sess) == 1, length(comp.list[[1]]), length(sess))))$Individual

pheno.nt <- targets %>%
  filter(Individual %in% ind) %>%
  filter(Session %in% sess) %>%
  filter(Time %in% comp.list[[1]])

probe.nt <- bVals.t[match(pheno.nt$Sample_Name, rownames(bVals.t)), ]

## Print the time
start.time <- Sys.time()
print(paste("ICC calculation start at", start.time))

## ICC(2,1) and ICC(2,k) for all the probes without adjustment
icc21.1 %<-% {ICC.n.m(probe.data = probe.nt, phenotype.data = pheno.nt, subject = "Individual", time = tors, probe = probe.name[1:round(probe.n/2,0)], 
                      icc.type = "ICC2", sample.individual = samp)} %seed% {seed=TRUE}

icc21.2 %<-% {ICC.n.m(probe.data = probe.nt, phenotype.data = pheno.nt, subject = "Individual", time = tors, probe = probe.name[(round(probe.n/2,0)+1):probe.n], 
                      icc.type = "ICC2", sample.individual = samp)} %seed% {seed=TRUE}


icc2k.1 %<-% {ICC.n.m(probe.data = probe.nt, phenotype.data = pheno.nt, subject = "Individual", time = tors, probe = probe.name[1:round(probe.n/2,0)], 
                      icc.type = "ICC2k", sample.individual = samp)} %seed% {seed=TRUE}

icc2k.2 %<-% {ICC.n.m(probe.data = probe.nt, phenotype.data = pheno.nt, subject = "Individual", time = tors, probe = probe.name[(round(probe.n/2,0)+1):probe.n], 
                      icc.type = "ICC2k", sample.individual = samp)} %seed% {seed=TRUE}


## ICC(2,1) and ICC(2,k) for all the probes adjusting for DNAm monocyte percentage and other covariates
icc21.covar.1 %<-% {ICC.n.m(probe.data = probe.nt, phenotype.data = pheno.nt, subject = "Individual", time = tors, probe = probe.name[1:round(probe.n/2,0)],
                           covariates = covar, icc.type = "ICC2", sample.individual = samp)} %seed% {seed=TRUE}

icc21.covar.2 %<-% {ICC.n.m(probe.data = probe.nt, phenotype.data = pheno.nt, subject = "Individual", time = tors, probe = probe.name[(round(probe.n/2,0)+1):probe.n],
                           covariates = covar, icc.type = "ICC2", sample.individual = samp)} %seed% {seed=TRUE}


icc2k.covar.1 %<-% {ICC.n.m(probe.data = probe.nt, phenotype.data = pheno.nt, subject = "Individual", time = tors, probe = probe.name[1:round(probe.n/2,0)],
                           covariates = covar, icc.type = "ICC2k", sample.individual = samp)} %seed% {seed=TRUE}

icc2k.covar.2 %<-% {ICC.n.m(probe.data = probe.nt, phenotype.data = pheno.nt, subject = "Individual", time = tors, probe = probe.name[(round(probe.n/2,0)+1):probe.n],
                           covariates = covar, icc.type = "ICC2k", sample.individual = samp)} %seed% {seed=TRUE}


## Combine all the results together
icc21 %<-% rbind(icc21.1, icc21.2)
icc2k %<-% rbind(icc2k.1, icc2k.2)
icc21.covar %<-% rbind(icc21.covar.1, icc21.covar.2)
icc2k.covar %<-% rbind(icc2k.covar.1, icc2k.covar.2)

icc.comb <- rbind(icc21,
                  icc2k,
                  icc21.covar,
                  icc2k.covar)

rm(icc21.1)
rm(icc21.2)
rm(icc2k.1)
rm(icc2k.2)
rm(icc21.covar.1)
rm(icc21.covar.2)
rm(icc2k.covar.1)
rm(icc2k.covar.2)

## Check progress
print(paste("ICC calculation and data combination done!", Sys.time()))
print(Sys.time()-start.time)

icc.comb$type_cov <- paste(icc.comb$type, icc.comb$covariates, sep = "_")

icc.comb$type <- factor(icc.comb$type, levels = c("ICC2", "ICC2k"),
                        labels = c("ICC(2,1)", "ICC(2,k)"))

icc.comb$covariates <- factor(icc.comb$covariates,
                              levels = c("No covariates", paste0(covar, collapse = "_")),
                              labels = c("No covariates", paste0("Adjusted for ", paste0(covar, collapse = ", "))))
icc.comb$session <- rep(paste0("s", paste0(sess, collapse = "s")), length(probe.name)*4)
icc.comb$time <- rep(paste0("t", paste0(comp.list[[1]], collapse = "t")), length(probe.name)*4)

print(paste("Variable adding to the combined dataset done!", Sys.time()))
print(Sys.time()-start.time)

## Save ICC data
assign(paste0("icc.comb.scen", scen, ifelse(samp == T, "sampled", "")), icc.comb)
save(list = paste0("icc.comb.scen", scen, ifelse(samp == T, "sampled", "")), 
     file = paste("Results/", "scen", scen, ifelse(samp == T, "sampled", ""), 
                  "_icc", "_s", paste0(sess, collapse = "s"), 
                  "_t", paste0(comp.list[[1]], collapse = "t"), 
                  "_ela", paste0(ela, collapse = ""), 
                  ifelse(revision == T, "_revision", ""), ".RData", sep = ""))
print(paste("Data saving done!", Sys.time()))
print(Sys.time()-start.time)

## Generate and export distribution plots
future({
  p <- ggplot(data = icc.comb, aes(x = ICC)) +
    geom_histogram(fill = ggplot2::alpha("gray", 0.3), color = "black", bins = 30) +
    facet_wrap(facets = covariates~type, nrow = 2, ncol = 2) +
    theme_classic()
  ggsave(file = paste("Results/", "scen", scen, ifelse(samp == T, "sampled", ""), 
                      "_icc_dist", "_s", paste0(sess, collapse = "s"), 
                      "_t", paste0(comp.list[[1]], collapse = "t"), 
                      "_ela", paste0(ela, collapse = ""), 
                      ifelse(revision == T, "_revision", ""), ".pdf", sep = ""),
         plot = p, height = 6, width = 10)
  
  print(paste("Distribution plot done!", Sys.time()))
  print(Sys.time()-start.time)
  
})

## Generate BA plots
pdf(file = paste("Results/", "scen", scen, ifelse(samp == T, "sampled", ""),
                 "_icc_baplot", "_s", paste0(sess, collapse = "s"), 
                 "_t", paste0(comp.list[[1]], collapse = "t"), 
                 "_ela", paste0(ela, collapse = ""), 
                 ifelse(revision == T, "_revision", ""), ".pdf", sep = ""),
    height = 6.5, width = 13.5)
par(mfrow=c(1,2))
BA.plot(icc21$ICC,icc21.covar$ICC, note = "ICC(2,1) vs. ICC(2,1) adjusted for covariates")
BA.plot(icc2k$ICC,icc2k.covar$ICC, note = "ICC(2,k) vs. ICC(2,k) adjusted for covariates")
dev.off()

print(paste("BA plot done!", Sys.time()))
print(Sys.time()-start.time)

## Test if there is a significant proportional bias in the Bland-Altman plots
icc21.baplot %<-% {data.frame(icc21.raw = icc21$ICC,
                              icc21.covariates = icc21.covar$ICC) %>%
    mutate(icc21.mean = (icc21.raw + icc21.covariates)/2,
           icc21.diff = icc21.raw - icc21.covariates)}

icc2k.baplot %<-% {data.frame(icc2k.raw = icc2k$ICC,
                              icc2k.covariates = icc2k.covar$ICC) %>%
    mutate(icc2k.mean = (icc2k.raw + icc2k.covariates)/2,
           icc2k.diff = icc2k.raw - icc2k.covariates)}

ba.biastest.res <- data.frame(ICC.type = c("icc21", "icc2k"),
                              Session = rep(paste0("s", paste0(sess, collapse = "s")), 2),
                              Time.points = rep(paste0("t", paste0(comp.list[[1]], collapse = "t")), 2),
                              rbind(coef(summary(lm(icc21.diff ~ icc21.mean, data = icc21.baplot)))[2,],
                                    coef(summary(lm(icc2k.diff ~ icc2k.mean, data = icc2k.baplot)))[2,]))

colnames(ba.biastest.res) <- c("ICC.type", "Session", "Time.points", "B1.estimate", "Std.error", "t.value", "p.value")

print(paste("BA plot proportional bias test done!", Sys.time()))
print(Sys.time()-start.time)

write_csv(ba.biastest.res, file = paste("Results/", "scen", scen, ifelse(samp == T, "sampled", ""),
                                        "_ba_biastest_res", "_s", paste0(sess, collapse = "s"), 
                                        "_t", paste0(comp.list[[1]], collapse = "t"), 
                                        "_ela", paste0(ela, collapse = ""), 
                                        ifelse(revision == T, "_revision", ""), ".csv",
                                        sep = ""))

print(paste("BA plot proportional bias test results saved!", Sys.time()))
print(Sys.time()-start.time)
print("All done!!!")

