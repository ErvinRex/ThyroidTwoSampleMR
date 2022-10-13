#Get the twosampleMR package
#install.packages("TwoSampleMR-master.zip", repos = NULL, type="source")

###upside down

library(TwoSampleMR)
library(MRPRESSO)



#ThyCa dataset as exposure
exp_dat <- read_exposure_data(
  filename = "ThyCa_SIG",
  sep = ",",
  snp_col = "rsID",
  beta_col = "OR.meta",
  effect_allele_col = "EA",
  other_allele_col = "OA",
  pval_col = "P.meta",
  chr_col = "CHR",
  pos_col = "PosB37",
  eaf_col = "EAF"
)

#taking out the insignificant p values
exp_dat <- exp_dat[exp_dat$pval.exposure < 0.00000005, ]

#log(OR) to obtain beta
exp_dat$beta.exposure=log(exp_dat$beta.exposure)

#time to add SE to the dataframe
#         beta            p-value     tail (2-tailed)
#se <- abs(-0.013896/qnorm(1.48375e-01/2))

se<-c(abs(exp_dat$beta.exposure/qnorm(exp_dat$pval.exposure/2)))

#add the new column to the df
exp_dat$se.exposure <-c(se)


#assuming samplesize of ThyCa is 290,551
#dat2_upside_down$samplesize.exposure <-replace(dat2_upside_down$samplesize.exposure, is.na(dat2_upside_down$samplesize.exposure), 290551)
#this seems to work instead in our flipped version
exp_dat$samplesize.exposure <- 290551

#PLINK clumping method - SNPs in LD within a particular window will be pruned, SNPs with lowest p-value is retained
exp_dat <- clump_data(
  exp_dat,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EUR"
)

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

dat2_upside_down <- harmonise_data(exp_dat, out_dat, action = 2)



#dat2$samplesize.outcome <-replace(dat2$samplesize.outcome, dat2$samplesize.outcome == 338077, 60000)


#dropping eaf.exposure and samplesize.outcome
#dat2 = subset(dat2, select = -c(samplesize.outcome,eaf.exposure) )

#dropping snps with pval.exposure of 0
#dat2 <- dat2[dat2$pval.exposure > 1.052760e-218, ]

#renaming id columns to better suit what we are looking at
dat2_upside_down$id.exposure <- replace(dat2_upside_down$id.exposure, is.character(dat2_upside_down$id.exposure), "Thyroid Cancer")
dat2_upside_down$id.outcome <- replace(dat2_upside_down$id.outcome, is.character(dat2_upside_down$id.outcome), "PTSD")

#remove the eaf.exposure column
#dat2_upside_down = subset(dat2_upside_down, select = -c(eaf.outcome))

#remove se.outcome 0 values (they cause issues)
dat2_upside_down$se.outcome[dat2_upside_down$se.outcome==0] <- NA

dat2_upside_down <- dat2_upside_down %>%
  na.omit()


load("rf.rdata")
res <- mr_wrapper(dat2_upside_down)
library("randomForest")
library("dplyr")
res_moe <- mr_moe(res, rf)
res_moe

#only steiger filtered
res_moe_filtered = subset(res_moe$`Thyroid Cancer.PTSD`$estimates, steiger_filtered== "TRUE")
res_moe_filtered = subset(res_moe_filtered, outlier_filtered == "TRUE")

#estimate MOE values
estimates <- res_moe$`Thyroid Cancer.PTSD`$estimates

#heterogeneity
heterogeneity <- res_moe$`Thyroid Cancer.PTSD`$heterogeneity

#directional pleiotropy
directional_pleiotropy <- res_moe$`Thyroid Cancer.PTSD`$directional_pleiotropy

#info
info <- res_moe$`Thyroid Cancer.PTSD`$info


# Perform MR
res_upside_down <- mr(dat2_upside_down, method_list = c("mr_simple_median", "mr_weighted_median", "mr_penalised_weighted_median"))

#generate ORs 
ORs <- generate_odds_ratios(res_moe_filtered)

#heterogeneity test - shows large q value (heterogeneity)
#and low p value - indicating substantial difference, none violating of MR assumptions
res_heterogeneity <- mr_heterogeneity(dat2_upside_down)

plot2 <-mr_scatter_plot(res_upside_down, dat2_upside_down)
plot2


#calculating horizontal pleiotropy - intercept term in MR Egger regression can be a useful
#indication of whether directional horizontal pleiotropy is driving the results of an MR analysis.
res_pleiotropy <- mr_pleiotropy_test(dat2_upside_down)


#performs the analysis multiple times for each exposure-outcome combination - 
#each time using a different single SNP to perform the analysis.
res_single <- mr_singlesnp(dat2_upside_down, all_method = c("mr_simple_median", "mr_weighted_median", "mr_penalised_weighted_median"))

#forest plots comparing MR values of each SNP
#compare the MR estimates using the different MR methods against the single SNP tests
p2 <- mr_forest_plot(res_single)
p2[[1]]

#Asymmetry in a funnel plot is useful for gauging the reliability of a particular MR analysis.
#Funnel plots can be produced using the single SNP results as follows
#funnel plot
p4 <- mr_funnel_plot(res_single)
p4[[1]]

#leave one out analysis, MR but each SNP left out in turn - identifies if a single SNP is driving an association
#can change it with "method" argument - default is IVW
#change method to egger regression
res_loo <- mr_leaveoneout(dat2_upside_down, method = mr_simple_median)

#leave one out plot
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]


###MR-PRESSO to confirm results

presso_res <- run_mr_presso(dat2_upside_down, NbDistribution = 1000, SignifThreshold = 0.05)

#get ORs of PRESSO
b = presso_res[[1]]$`Main MR results`$`Causal Estimate`
se = presso_res[[1]]$`Main MR results$Sd
df <- data.frame(b, se)
presso_ORs <- generate_odds_ratios(df)
