
# Clean up
rm(list=ls())


# Function to simulate survival data for a GxE experiment and asses power in statistical tests based on LMM an CMM. 1 treatment and 1 control group, where effect size in lnHR. Takes assumed G sd in intercepts and treatments, and corr_IT. Specify n-strains and n-per group. k is number of datasets with exactly same simualted genetic effects. default of 1. Output is a list of length k + 1, where the first k are the datasets and the last is the G-effects sampled
  
# Returns a set reps length statsitcal results on estimates and significance 
# The output is a dataframe containing reps rows of replicated outcomes with the parameters estimated on each replicate
# MD is the estimated mean difference in lifespan for a typical strain from the LMM
# MD_p is the p value for MD
# MD_sigma_G is the estimated among-strain variance in MD
# MD_sigma_G_p is a p-value for MD_sigma_G based on the improvment in model fit test (chi-sq test on change in deviance)
# MD_R is the correlation between strain specific MD in two seperate experiments
# lnHR is the estimated log hazard ratio for a typical strain from the CMM
# lnHR_p is the p value for lnHR
# lnHR_sigma_G is the estimated among-strain variance in lnHR
# lnHR_sigma_G_p is a p-value for lnHR_sigma_G based on the improvment in model fit test (chi-sq test on change in deviance)
# lnHR_R is the correlation between strain specific lnHR in two seperate experiments

# Dependent packages are - make sure they are installed
# (simsurv)
# (lmerTest)
# (coxme)
# (mvtnorm)

power_sim<-function(n_strains, n, lnHR, sigma_IntG, sigma_TrtG, cor_IntTrt, a, b, reps){
	
	# Packages
	require(simsurv)
	require(lmerTest)
	require(coxme)
	require(mvtnorm)
	
	# Object to hole results
	results<-data.frame(replicate=seq(1, reps, 1), MD=NA, MD_p=NA, MD_sigma_G=NA, MD_sigma_G_p=NA, MD_R=NA, lnHR=NA, lnHR_p=NA, lnHR_sigma_G=NA, lnHR_sigma_G_p=NA, lnHR_R=NA)
	
	# model matrix
	ints<-1
	trts<-c(1, 0)
	covs<-data.frame(id=seq(1, n_strains*n*length(trts)), strain=sort(rep(seq(1, n_strains, 1), length(trts)*n)), ints=ints, trts=trts)
	
	# How strong should the main effect of DR be
	main_effect<-lnHR
	
	# Put in to a VCV
	vcv<-array(cor_IntTrt*sigma_IntG*sigma_TrtG, c(2,2))
	diag(vcv)<-c(sigma_IntG, sigma_TrtG)^2
	
	# Loop to repeat for replicates - progress bar as well
	pb<-txtProgressBar(min=0, max=reps, initial=0)
	for(r in 1:reps){
	
		# Replicated using while to make sure all 6 models fit
		model_fit<-F
		while(model_fit == F){
			
			# Simulate strain specific effects for rmvnorm
			strain_specific_effects<-rmvnorm(n_strains, mean=c(0, main_effect), sigma=vcv)
			
			# Add in to the predictor data frame
			pars<-data.frame(ints = rep(strain_specific_effects[,1], each=n*length(trts)), trts = rep(strain_specific_effects[,2], each=n*length(trts)))
			
			# List for output
			out<-list()
			
			# Loop length 2 simulating 2 datasets and adding to a list - we simulate 2 to estimate replicability
			for(i in 1:2){
				# Simulate - parcel up and return
				sim<-simsurv(dist="gompertz", lambdas = a, gammas = b, x=covs, betas=pars)
				dat<-data.frame(covs, sim[,-1])
				dat$diet<-"AL"
				dat$diet[which(dat$trts == 1)]<-"DR"
				out[[i]]<-dat
			}
			
			# Analyse the data via LMM and CMM - note I am supressing warnigs
			oldw<-getOption("warn")
			options(warn=-1)
			lmm<-try(suppressMessages(lmer(eventtime ~ trts + (1 + trts|strain), data=out[[1]])), silent=T)
			cmm<-try(coxme(Surv(eventtime, status) ~ trts + (1 + trts|strain), data=out[[1]]), silent=T)
			
			# Fit null models
			lmm0<-try(suppressMessages(lmer(eventtime ~ trts + (1|strain), data=out[[1]])), silent=T)
			cmm0<-try(coxme(Surv(eventtime, status) ~ trts + (1|strain), data=out[[1]]), silent=T)
			
			# Refit from the second data set to get replicability
			lmm2<-try(suppressMessages(lmer(eventtime ~ trts + (1 + trts|strain), data=out[[2]])), silent=T)
			cmm2<-try(coxme(Surv(eventtime, status) ~ trts + (1 + trts|strain), data=out[[2]]), silent=T)
			options(warn=oldw)
				
			# Assess models for fit and save the results if we have some
			if(class(lmm) != "try-error" & class(cmm) != "try-error" & class(lmm0) != "try-error" & class(cmm0) != "try-error" & class(lmm2) != "try-error" & class(cmm2) != "try-error"){
				
				# Set model fit to T because we are happy
				model_fit<-T
				
				# Save the output from the lmm approach
				results$MD[r]<-lmm@beta[2]
				results$MD_p[r]<-summary(lmm)$coeff[2,5]
				results$MD_sigma_G[r]<-sqrt(summary(lmm)$varcor$strain[2,2])
				results$MD_sigma_G_p[r]<-anova(lmm, lmm0, refit=F)[2,8]
				results$MD_R[r]<-cor(ranef(lmm)$strain[,2], ranef(lmm2)$strain[,2])
				
				# Save the output from the cmm
				results$lnHR[r]<-cmm$coef[1]
				results$lnHR_p[r]<-pnorm(-abs(cmm$coef[1] / sqrt(vcov(cmm))), 0, 1) * 2
				results$lnHR_sigma_G[r]<-sqrt(VarCorr(cmm)$strain[2,2])
				results$lnHR_sigma_G_p[r]<-anova(cmm, cmm0)[2,4]
				results$lnHR_R[r]<-cor(ranef(cmm)$strain[,2], ranef(cmm2)$strain[,2])
				
			}
			
			# Clean up the models
			rm(lmm)
			rm(cmm)
			rm(lmm0)
			rm(cmm0)
			rm(lmm2)
			rm(cmm2)
			
		}
		
		# Update the progress bar
		setTxtProgressBar(pb, r)
	}
	
	# Close the progress bar and return the results
	close(pb)
	return(results)
}

# Here I demonstrate application of the function, based on the parameterisation described in eText 3

# As an example lets get 10 replicates per group with 30 strains (i.e., n = 10 per group)
set.seed(123) # Setting seed for replicability
results<-power_sim(n_strains=30, n=10, lnHR=-0.5, sigma_IntG=0, sigma_TrtG=0.5, cor_IntTrt=0, a=exp(-11.57), b=exp(-4.9), reps=20)

# Check the results
dim(results)
# We have 20 rows - 1 per replicate
head(results)
# All the different columns as per above in the function description

# Median MD is gain 62.8 days IQR goes from 47.8 days to 76 days
quantile(results$MD, c(0.25, 0.5, 0.75))

# Median SD of strains is 70.5 days, IQR is 56.8 to 81.6 days
quantile(results$MD_sigma_G, c(0.25, 0.5, 0.75))

# power to detect among starin variance at alpha = 0.05 is 0.8, this is acceptable, but note this is just 20 reps
mean(results$MD_sigma_G_p < 0.05)

# WHat about CVG for MD?
# Median is 1.41 so a bit of an overestimate
quantile(results$MD_sigma_G / abs(results$MD), c(0.25, 0.5, 0.75))

# For the lnHR based on the CMM we get median SD 0.48, with IQR 0.43 to 0.6 - the true value is 0.5 so not bad
quantile(results$lnHR_sigma_G, c(0.25, 0.5, 0.75))

# Power to detect variance in response is 0.95 
mean(results$lnHR_sigma_G_p < 0.05)

# CVG is an over estimate
quantile(results$lnHR_sigma_G / abs(results$lnHR), c(0.25, 0.5, 0.75))

# OK to be more systemtic lets do 100 reps at n = 10, 15, 20, 25 for 30 and 50 strains
# Starting with 100 reps is a good way to map the space, but confirmation should really do 10,000 reps
# Here I do it long hand to make easier to see what is going on, but below is a loop that automates for potentially more ns

# 30 with all different ns
results30_10<-power_sim(n_strains=30, n=10, lnHR=-0.5, sigma_IntG=0, sigma_TrtG=0.5, cor_IntTrt=0, a=exp(-11.57), b=exp(-4.9), reps=100)
results30_10$n<-10
results30_10$n_strains<-30
results30_15<-power_sim(n_strains=30, n=15, lnHR=-0.5, sigma_IntG=0, sigma_TrtG=0.5, cor_IntTrt=0, a=exp(-11.57), b=exp(-4.9), reps=100)
results30_15$n<-15
results30_15$n_strains<-30
results30_20<-power_sim(n_strains=30, n=20, lnHR=-0.5, sigma_IntG=0, sigma_TrtG=0.5, cor_IntTrt=0, a=exp(-11.57), b=exp(-4.9), reps=100)
results30_20$n<-20
results30_20$n_strains<-30
results30_25<-power_sim(n_strains=30, n=25, lnHR=-0.5, sigma_IntG=0, sigma_TrtG=0.5, cor_IntTrt=0, a=exp(-11.57), b=exp(-4.9), reps=100)
results30_25$n<-25
results30_25$n_strains<-30


# 50 with all different ns
results50_10<-power_sim(n_strains=50, n=10, lnHR=-0.5, sigma_IntG=0, sigma_TrtG=0.5, cor_IntTrt=0, a=exp(-11.57), b=exp(-4.9), reps=100)
results50_10$n<-10
results50_10$n_strains<-50
results50_15<-power_sim(n_strains=50, n=15, lnHR=-0.5, sigma_IntG=0, sigma_TrtG=0.5, cor_IntTrt=0, a=exp(-11.57), b=exp(-4.9), reps=100)
results50_15$n<-15
results50_15$n_strains<-50
results50_20<-power_sim(n_strains=50, n=20, lnHR=-0.5, sigma_IntG=0, sigma_TrtG=0.5, cor_IntTrt=0, a=exp(-11.57), b=exp(-4.9), reps=100)
results50_20$n<-20
results50_20$n_strains<-50
results50_25<-power_sim(n_strains=50, n=25, lnHR=-0.5, sigma_IntG=0, sigma_TrtG=0.5, cor_IntTrt=0, a=exp(-11.57), b=exp(-4.9), reps=100)
results50_25$n<-25
results50_25$n_strains<-50


# Stick them all together
results<-rbind(results30_10, results30_15, results30_20, results30_25, results50_10, results50_15, results50_20, results50_25)
head(results)

# Calculate CVG
results$MD_CVG<-results$MD_sigma_G / abs(results$MD)
results$lnHR_CVG<-results$lnHR_sigma_G / abs(results$lnHR)

# Calculate power to detect
library(plyr)
results$code<-paste0(results$n_strains, "_", results$n)
power<-ddply(results, .(code), summarise, n=n[1], n_strains=n_strains[1], power_MD=mean(MD_sigma_G_p < 0.05), power_lnHR = mean(lnHR_sigma_G_p < 0.05))

# Plot using ggplot
library(ggplot2)

pdf("FS3.pdf", height=4.5, width=5)
p1<-ggplot(power, aes(x=n, y=power_MD, col=n_strains)) + geom_point() + geom_hline(yintercept=0.8) + ylab("Power (LMM)") + labs(title="A.") + ylim(0, 1)
print(p1)
p1<-ggplot(power, aes(x=n, y=power_lnHR, col=n_strains)) + geom_point() + geom_hline(yintercept=0.8) + ylab("Power (CMM)") + labs(title="B.") + ylim(0, 1)
print(p1)
p1<-ggplot(results, aes(x=n, group=paste0(n, n_strains), y=MD_CVG, col=n_strains)) + geom_boxplot() + geom_hline(yintercept=1) + ylab("Esitmated MD CVG") + labs(title="C.")
print(p1)
p1<-ggplot(results, aes(x=n, group=paste0(n, n_strains), y=lnHR_CVG, col=n_strains)) + geom_boxplot() + geom_hline(yintercept=1) + ylab("Esitmated lnHR CVG") + labs(title="D.")
print(p1)
p1<-ggplot(results, aes(x=n, group=paste0(n, n_strains), y=MD_R, col=n_strains)) + geom_boxplot() + geom_hline(yintercept=1) + ylab("Replicability MD") + labs(title="E.")
print(p1)
p1<-ggplot(results, aes(x=n, group=paste0(n, n_strains), y=lnHR_R, col=n_strains)) + geom_boxplot() + geom_hline(yintercept=1) + ylab("Replicability lnHR") + labs(title="F.")
print(p1)
dev.off()

# Done with a loop

# # # Object with the ns we want to assess
# ns<-c(10, 15, 20, 25)

# # Do the first set
# results<-power_sim(n_strains=30, n=ns[1], lnHR=-0.5, sigma_IntG=0, sigma_TrtG=0.5, cor_IntTrt=0, a=exp(-11.57), b=exp(-4.9), reps=100)
# results$n<-ns[1]
# # Now do the others
# for(i in 2:length(ns)){
	# # Get an ith set for that n
	# results_i<-power_sim(n_strains=30, n=ns[i], lnHR=-0.5, sigma_IntG=0, sigma_TrtG=0.5, cor_IntTrt=0, a=exp(-11.57), b=exp(-4.9), reps=100)
	# results_i$n<-ns[i]
	# # Add it in to the results
	# results<-rbind(results, results_i)
# }

# # Lets repeat with n_strains 50 and see what happens
# results_50<-power_sim(n_strains=50, n=ns[1], lnHR=-0.5, sigma_IntG=0, sigma_TrtG=0.5, cor_IntTrt=0, a=exp(-11.57), b=exp(-4.9), reps=100)
# results_50$n<-ns[1]
# # Now do the others
# for(i in 2:length(ns)){
	# # Get an ith set for that n
	# results_i<-power_sim(n_strains=50, n=ns[i], lnHR=-0.5, sigma_IntG=0, sigma_TrtG=0.5, cor_IntTrt=0, a=exp(-11.57), b=exp(-4.9), reps=100)
	# results_i$n<-ns[i]
	# # Add it in to the results
	# results_50<-rbind(results_50, results_i)
# }

# results<-rbind(results, results_50)
