library(TwoSampleMR)
library(MRPRESSO)

exp_dat <- read_exposure_data(
  filename = "FT4_FINAL",
  sep = ",",
  chr_col = "Chr",
  pos_col = "Pos",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P",
  samplesize_col = "TotalSampleSize"
)

#taking out the insignificant p values
exp_dat <- exp_dat[exp_dat$pval.exposure < 0.0000005, ]

snp_list <- exp_dat$SNP

#PTSD dataset as outcome
out_dat <- read_outcome_data(
  snps = exp_dat$SNP,
  filename = "PTSD",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "OR",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P",
  ncase_col = "Nca",
  ncontrol_col = "Nco",
  eaf_col = "FRQ_A_23212"
)

#log(OR) to obtain beta
out_dat$beta.outcome=log(out_dat$beta.outcome)

#time to add SE to the dataframe
#         beta            p-value     tail (2-tailed)
#se <- abs(-0.013896/qnorm(1.48375e-01/2))

se<-c(abs(out_dat$beta.outcome/qnorm(out_dat$pval.outcome/2)))

#add the new column to the df
out_dat$se.outcome <-c(se)


dat2 <- harmonise_data(exp_dat, out_dat, action = 2)

#fake columns removing the infinity in se.exposure SNPs - they don't work
#dat2 <- dat2[-c(95, 401, 582), ]

#renaming id columns to better suit what we are looking at
dat2$id.exposure <- replace(dat2$id.exposure, is.character(dat2$id.exposure), "FT4")
dat2$id.outcome <- replace(dat2$id.outcome, is.character(dat2$id.outcome), "PTSD")

#remove the eaf.exposure column
dat2 = subset(dat2, select = -c(eaf.exposure))

#remove se.outcome 0 values (they cause issues)
dat2$se.outcome[dat2$se.outcome==0] <- NA

dat2 <- dat2 %>%
  na.omit()

#need this package for mr_moe to work
library(dplyr)

load("rf.rdata")
res <- mr_wrapper(dat2)
library("randomForest")
res_moe <- mr_moe(res, rf)
res_moe



# Now you can view the estimates, and see that they have 
# been sorted in order from most likely to least likely to 
# be accurate, based on MOE prediction
res_moe[[1]]$estimates


#only steiger filtered
res_moe_filtered = subset(res_moe$FT4.PTSD$estimates, steiger_filtered== "TRUE")
res_moe_filtered = subset(res_moe_filtered, outlier_filtered == "TRUE")

#estimate MOE values
estimates <- res_moe$FT4.PTSD$estimates

#heterogeneity
heterogeneity <- res_moe$FT4.PTSD$heterogeneity

#directional pleiotropy
directional_pleiotropy <- res_moe$FT4.PTSD$directional_pleiotropy

#info
info <- res_moe$FT4.PTSD$info

#snps retained
snps_retained <- res_moe$FT4.PTSD$snps_retained

###Now we know which methods to use that are best suited to our data

# Perform MR
res <- mr(dat2, method_list = c("mr_weighted_mode", "mr_ivw", "mr_penalised_weighted_median"))

#generate ORs 
ORs <- generate_odds_ratios(res_moe_filtered)

#heterogeneity test - shows large q value (heterogeneity)
#and low p value - indicating substantial difference, none violating of MR assumptions
res_heterogeneity <- mr_heterogeneity(dat2)


plot2 <-mr_scatter_plot(res, dat2)
plot2

#performs the analysis multiple times for each exposure-outcome combination - 
#each time using a different single SNP to perform the analysis.
#forest plots comparing MR values of each SNP
#dat100 <- dat2[1:20,]
res_single <- mr_singlesnp(dat2, all_method = c("mr_weighted_mode", "mr_ivw", "mr_penalised_weighted_median"))

#compare the MR estimates using the different MR methods against the single SNP tests
p2 <- mr_forest_plot(res_single)
p2[[1]]


#calculating horizontal pleiotropy - intercept term in MR Egger regression can be a useful
#indication of whether directional horizontal pleiotropy is driving the results of an MR analysis.
res_pleiotropy <- mr_pleiotropy_test(dat2)


#Asymmetry in a funnel plot is useful for gauging the reliability of a particular MR analysis.
#Funnel plots can be produced using the single SNP results as follows
#funnel plot
p4 <- mr_funnel_plot(res_single)
p4[[1]]

#leave one out analysis, MR but each SNP left out in turn - identifies if a single SNP is driving an association
res_loo <- mr_leaveoneout(dat2, method = mr_ivw)
#can change it with "method" argument - default is IVW
#change method to egger regression

#leave one out plot
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

###MRPRESSO results to confirm
res_presso<- mr_presso(BetaOutcome = "beta.outcome",
                     BetaExposure = "beta.exposure",
                     SdOutcome = "se.outcome",
                     SdExposure = "se.exposure",
                     data = dat2,
                     OUTLIERtest = TRUE,
                     DISTORTIONtest = TRUE
)

#get ORs of PRESSO
b = res_presso$`Main MR results`$`Causal Estimate`
se = res_presso$`Main MR results`$Sd
df <- data.frame(b, se)
presso_ORs <- generate_odds_ratios(df)

###NOW THE OTHER WAY ROUND
#remember you need to log(OR), to get beta
exp_dat <- read_exposure_data(
  filename = "PTSD_diff_sig",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "OR",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P",
  chr_col = "CHR",
  pos_col = "BP",
  ncase_col = "Nca",
  ncontrol_col = "Nco",
  eaf_col = "FRQ_A_23212"
  #clump = TRUE
)

#taking out the insignificant p values - anything <0.05
exp_dat <- exp_dat[exp_dat$pval.exposure < 0.00005, ]

#log(OR) to obtain beta
exp_dat$beta.exposure=log(exp_dat$beta.exposure)

#PLINK clumping method - SNPs in LD within a particular window will be pruned, SNPs with lowest p-value is retained
exp_dat1 <- clump_data(
  exp_dat,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EUR"
)

out_dat <- read.csv("FT4_outcome", sep = ",")

out_dat = subset(out_dat, select = -c(X))



exp_dat$id.exposure <- "PTSD"

out_dat$id.outcome <- "FT4"

out_dat$outcome <- "outcome"
out_dat$mr_keep.outcome <- "TRUE"
out_dat$pval_origin.outcome <- "reported"
out_dat$data_source.outcome <- "textfile"

exp_dat <- subset(exp_dat, select = -c(pos.exposure))
out_dat <- subset(out_dat, select = -c(pos.outcome))


dat2 <- harmonise_data(exp_dat, out_dat, action = 2)

#renaming id columns to better suit what we are looking at
dat2$id.exposure <- replace(dat2$id.exposure, is.character(dat2$id.exposure), "PTSD")
dat2$id.outcome <- replace(dat2$id.outcome, is.character(dat2$id.outcome), "FT4")

#remove se.outcome 0 values (they cause issues)
dat2$se.outcome[dat2$se.outcome==0] <- NA

dat2 <- dat2 %>%
  na.omit()

#need this package for mr_moe to work
library(dplyr)

load("rf.rdata")
res <- mr_wrapper(dat2)
library("randomForest")
res_moe <- mr_moe(res, rf)
res_moe

# Now you can view the estimates, and see that they have 
# been sorted in order from most likely to least likely to 
# be accurate, based on MOE prediction
res_moe[[1]]$estimates


#only steiger filtered
res_moe_filtered = subset(res_moe$PTSD.FT4$estimates, steiger_filtered== "TRUE")
res_moe_filtered = subset(res_moe_filtered, outlier_filtered == "TRUE")

#estimate MOE values
estimates <- res_moe$PTSD.FT4$estimates

#heterogeneity
heterogeneity <- res_moe$PTSD.FT4$heterogeneity

#directional pleiotropy
directional_pleiotropy <- res_moe$PTSD.FT4$directional_pleiotropy

#info
info <- res_moe$PTSD.FT4$info

#snps retained
snps_retained <- res_moe$PTSD.FT4$snps_retained

###Now we know which methods to use that are best suited to our data

# Perform MR
res <- mr(dat2, method_list = c("mr_ivw", "mr_weighted_median", "mr_simple_mode"))

#generate ORs 
ORs <- generate_odds_ratios(res_moe_filtered)

#heterogeneity test - shows large q value (heterogeneity)
#and low p value - indicating substantial difference, none violating of MR assumptions
res_heterogeneity <- mr_heterogeneity(dat2)


plot2 <-mr_scatter_plot(res, dat2)
plot2

#performs the analysis multiple times for each exposure-outcome combination - 
#each time using a different single SNP to perform the analysis.
#forest plots comparing MR values of each SNP
#dat100 <- dat2[1:20,]
res_single <- mr_singlesnp(dat2, all_method = c("mr_ivw", "mr_weighted_median", "mr_simple_mode"))

#compare the MR estimates using the different MR methods against the single SNP tests
p2 <- mr_forest_plot(res_single)
p2[[1]]


#calculating horizontal pleiotropy - intercept term in MR Egger regression can be a useful
#indication of whether directional horizontal pleiotropy is driving the results of an MR analysis.
res_pleiotropy <- mr_pleiotropy_test(dat2)


#Asymmetry in a funnel plot is useful for gauging the reliability of a particular MR analysis.
#Funnel plots can be produced using the single SNP results as follows
#funnel plot
p4 <- mr_funnel_plot(res_single)
p4[[1]]

#leave one out analysis, MR but each SNP left out in turn - identifies if a single SNP is driving an association
res_loo <- mr_leaveoneout(dat2, method = mr_ivw)
#can change it with "method" argument - default is IVW
#change method to egger regression

#leave one out plot
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

###MRPRESSO results to confirm
res_presso<- mr_presso(BetaOutcome = "beta.outcome",
                       BetaExposure = "beta.exposure",
                       SdOutcome = "se.outcome",
                       SdExposure = "se.exposure",
                       data = dat2,
                       OUTLIERtest = TRUE,
                       DISTORTIONtest = TRUE
)

#get ORs of PRESSO
b = res_presso$`Main MR results`$`Causal Estimate`
se = res_presso$`Main MR results`$Sd
df <- data.frame(b, se)
presso_ORs <- generate_odds_ratios(df)

###NO MATCHING SNPS FOR THE REVERSE
