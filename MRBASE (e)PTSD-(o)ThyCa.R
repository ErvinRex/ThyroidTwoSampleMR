#Get the twosampleMR package
install.packages("TwoSampleMR-master.zip", repos = NULL, type="source")

#get into the workding directory of PTSD-ThyCa
setwd("C:/Users/Ervin/Desktop/Work/QMUL MSc Bioinformatics/Semester B/BIO702P - BIOINFORMATICS RESEARCH PROJECT/MR/PTSD-ThyCa")

library(TwoSampleMR)
#remember you need to log(OR), to get beta
exp_dat <- read_exposure_data(
  filename = "PTSD_SIG",
  sep = ",",
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

#matching the snps to samplesize
exp_dat <- merge(exp_dat, samplesize.exposure, by="SNP")

#log(OR) to obtain beta
exp_dat$beta.exposure=log(exp_dat$beta.exposure)

#PLINK clumping method - SNPs in LD within a particular window will be pruned, SNPs with lowest p-value is retained
exp_dat <- clump_data(
  exp_dat,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EUR"
)

#getting outcome data into viable MR input
#remember to log(OR.meta) to get beta
#problem is matching SNPs here, try get some other chromosomes (+ in out it is written chr1, in exp it is written 1)

out_dat <- read_outcome_data(
  snps = exp_dat$SNP,
  filename = "ThyCa",
  sep = ",",
  snp_col = "rsID",
  beta_col = "OR.meta",
  effect_allele_col = "EA",
  other_allele_col = "OA",
  eaf_col = "EAF",
  pval_col = "P.meta"
)

#taking out the insignificant p values - anything <0.05
#not needed for outcome
#out_dat <- out_dat[out_dat$pval.outcome < 0.05, ]

#log(OR) to obtain beta
out_dat$beta.outcome=log(out_dat$beta.outcome)

#time to add SE to the dataframe
#         beta            p-value     tail (2-tailed)
#se <- abs(-0.013896/qnorm(1.48375e-01/2))

se<-c(abs(out_dat$beta.outcome/qnorm(out_dat$pval.outcome/2)))

#add the new column to the df
out_dat$se.outcome <-c(se)

#set column as numeric to do results
#out_dat$se.outcome <- as.numeric("\\.", "", out_dat$se.outcome)

#At first, the no. of variables was 0 - why is that happening then?
#THEY HAVE TO HAVE OVERLAPPING SNPS, OTHERWISE IT DOESN'T GET SELECTED OUT
#action 3 = Correct strand for non-palindromic SNPs, and drop all palindromic SNPs from the analysis (more conservative). 
#If a single value is passed then this action is applied to all outcomes. But multiple values can be supplied as a vector, each element relating to a different outcome.


dat2 <- harmonise_data(exp_dat, out_dat, action = 2)

#fake columns removing the infinity in se.exposure SNPs - they don't work
#dat2 <- dat2[-c(95, 401, 582), ]

#renaming id columns to better suit what we are looking at
dat2$id.exposure <- replace(dat2$id.exposure, is.character(dat2$id.exposure), "PTSD")
dat2$id.outcome <- replace(dat2$id.outcome, is.character(dat2$id.outcome), "Thyroid Cancer")

# Get effects of instruments on outcome
#outcome_dat <- extract_outcome_data(snps=exp_dat$SNP, outcomes="ieu-b-4914")

#assuming samplesize of ThyCa is 290,551
dat2$samplesize.outcome <-replace(dat2$samplesize.outcome, is.na(dat2$samplesize.outcome), 290551)

#remove the eaf.exposure column
#dat2 = subset(dat2, select = -c(eaf.exposure))

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
res_moe_filtered = subset(res_moe$`PTSD.Thyroid Cancer`$estimates, steiger_filtered== "TRUE")
res_moe_filtered = subset(res_moe_filtered, outlier_filtered == "TRUE")

#merge with the remaining results
res_moe_filtered <- merge(res_moe$`PTSD.Thyroid Cancer`, res_moe_filtered, by="MOE") 

#estimate MOE values
estimates <- res_moe$`PTSD.Thyroid Cancer`$estimates

#heterogeneity
heterogeneity <- res_moe$`PTSD.Thyroid Cancer`$heterogeneity

directional_pleiotropy <- res_moe$`PTSD.Thyroid Cancer`$directional_pleiotropy

#info
info <- res_moe$`PTSD.Thyroid Cancer`$info

#snps retained
snps_retained <- res_moe$`PTSD.Thyroid Cancer`$snps_retained

###Now we know which methods to use that are best suited to our data

# Perform MR
res <- mr(dat2, method_list = c("mr_weighted_median", "mr_ivw", "mr_weighted_mode"))

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
res_single <- mr_singlesnp(dat2, all_method = c("mr_weighted_median", "mr_ivw", "mr_weighted_mode"))

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


###COMBINE SNPS RETAINED WITH EACH DATAFRAME

directionality_test(dat2)
