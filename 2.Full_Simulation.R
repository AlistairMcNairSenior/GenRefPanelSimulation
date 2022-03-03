
# Clean up
rm(list=ls())

# Load packages
library(simsurv)
library(lme4)
library(coxme)
library(mvtnorm)

# Get the incoming info from pbs
args<-commandArgs()
job<-as.numeric(args[6])

# Function to simulate survival data for a GxE experiment. 1 treatment and 1 control group, where effect size in lnHR. Takes assemend G sd in intercepts and treatments, and corr_IT. Specify n-strains and n-per group. k is number of datasets with exactly same simualted genetic effects. default of 1. Output is a list of length k + 1, where the first k are the datasets and the last is the G-effects sampled
  
# Data are simulated assuming Gompertz. What to take as the lambda (scale) and gamma (shape)? Simons et al. perform a quasi meta-analysis of DR in rodents giving the ln(a) = -11.57 and ln(b) = -4.90, where a * exp(b * t) is their equation - these are defaults for a and b

my_sim<-function(n_strains, n, lnHR, sigma_IntG, sigma_TrtG, cor_IntTrt, k=1, a=exp(-11.57), b=exp(-4.90)){
	
	# model matrix
	ints<-1
	trts<-c(1, 0)
	covs<-data.frame(id=seq(1, n_strains*n*length(trts)), strain=sort(rep(seq(1, n_strains, 1), length(trts)*n)), ints=ints, trts=trts)
	
	# How strong should the main effect of DR be
	main_effect<-lnHR
	
	# Put in to a VCV
	vcv<-array(cor_IntTrt*sigma_IntG*sigma_TrtG, c(2,2))
	diag(vcv)<-c(sigma_IntG, sigma_TrtG)^2
	
	# Simulate strain specific effects for rmvnorm
	strain_specific_effects<-rmvnorm(n_strains, mean=c(0, main_effect), sigma=vcv)
	
	# Add in to the predictor data frame
	pars<-data.frame(ints = rep(strain_specific_effects[,1], each=n*length(trts)), trts = rep(strain_specific_effects[,2], each=n*length(trts)))
	
	# List for output
	out<-list()
	
	# Loop length k simulating and adding to list
	for(i in 1:k){
		# Simulate - parcel up and return
		sim<-simsurv(dist="gompertz", lambdas = a, gammas = b, x=covs, betas=pars)
		dat<-data.frame(covs, sim[,-1])
		dat$diet<-"AL"
		dat$diet[which(dat$trts == 1)]<-"DR"
		out[[i]]<-dat
	}
	
	# Add in the G effects
	out[[k+1]]<-as.data.frame(strain_specific_effects)
	names(out[[k+1]])<-c("Int_G", "Trt_G")
	
	# Hand it all back
	return(out)

}


# Lets do a series of simulations analysing bias as a function of n and variance using the waterfall plot method and coxmm
ns<-seq(5, 40, 5)
sigma_G<-c(0, 0.2, 0.5, 1)
reps<-1

# List to aggregate the full results from all G
full_results<-list()

# run function for each part across the cores
for(g in 1:length(sigma_G)){
	
	# To hold this set of results
	total_results<-list()
	
	# Loop for each n value
	for(j in 1:length(ns)){
		
		# objects to holds the CVs and keep a running count of significant ps
		raw_CV_np_md<-NULL
		raw_CV_np_hr<-NULL
		raw_CV_pp_md<-NULL
		raw_CV_pp_hr<-NULL
		raw_corrs_np_md<-NULL
		raw_corrs_np_hr<-NULL
		raw_corrs_pp_md<-NULL
		raw_corrs_pp_hr<-NULL
		p<-0
		
		# To do the replicates
		for(k in 1:reps){
			
			# Simulate the two datasets assuming the same GxE Matrix
			sim<-my_sim(n_strain=40, n=ns[j], lnHR=-0.5, sigma_IntG=0, sigma_Tr=sigma_G[g], cor_IntTrt=0, k=2)
			
			# Analyse the first for plotting
			dat<-sim[[1]]
			
			# Now get the strain by strain effect sizes - calculated as t-tests
			effects<-data.frame(strains = seq(1, 40, 1), np_md=0, np_hr=0)
			for(i in 1:nrow(effects)){
				data_i<-dat[which(dat$strain == i),]
				test<-try(lm(eventtime ~ trts, data=data_i), silent=T)
				effects$np_md[i]<-test$coef[2]
				test<-try(coxph(Surv(eventtime, status) ~ trts, data=data_i), silent=T)
				effects$np_hr[i]<-test$coef[1]				
			}
			
			# now we use lmm and coxme
			lmm<-try(lmer(eventtime ~ trts + (1 + trts|strain), data=dat), silent=T)
			coxmm<-try(coxme(Surv(eventtime, status) ~ trts + (1 + trts|strain), data=dat), silent=T)
			
			# Save the CVs
			raw_CV_np_md<-c(raw_CV_np_md, (sd(effects$np_md) / abs(mean(effects$np_md))))
			raw_CV_np_hr<-c(raw_CV_np_hr, (sd(effects$np_hr) / abs(mean(effects$np_hr))))
			raw_CV_pp_md<-c(raw_CV_pp_md, (sqrt(VarCorr(lmm)[1]$strain[2,2]) / abs(lmm@beta[2])))
			raw_CV_pp_hr<-c(raw_CV_pp_hr, (sqrt(coxmm$vcoef$strain[2,2]) / abs(coxmm$coeff)))
								
			# Now the LRT
			coxmm0<-coxme(Surv(eventtime, status) ~ trts + (1|strain), data=dat)
			p<-p + (anova(coxmm0, coxmm)$P[2] < 0.05)
			
			# Now analyse the second and get the correlation between genotypes
			effects$np_md2<-NA
			effects$np_hr2<-NA
			dat<-sim[[2]]
			for(i in 1:nrow(effects)){
				data_i<-dat[which(dat$strain == i),]
				test<-try(lm(eventtime ~ trts, data=data_i), silent=T)
				effects$np_md2[i]<-test$coef[2]
				test<-try(coxph(Surv(eventtime, status) ~ trts, data=data_i), silent=T)
				effects$np_hr2[i]<-test$coef[1]		
			}
						
			# now we use lmm and coxme
			lmm2<-try(lmer(eventtime ~ trts + (1 + trts|strain), data=dat), silent=T)
			coxmm2<-try(coxme(Surv(eventtime, status) ~ trts + (1 + trts|strain), data=dat), silent=T)

			# Correlation between experiments
			raw_corrs_np_md<-c(raw_corrs_np_md, cor(effects$np_md, effects$np_md2))
			raw_corrs_np_hr<-c(raw_corrs_np_hr, cor(effects$np_hr, effects$np_hr2))
			raw_corrs_pp_md<-c(raw_corrs_pp_md, cor(ranef(lmm)$strain[,2], ranef(lmm2)$strain[,2]))
			raw_corrs_pp_hr<-c(raw_corrs_pp_hr, cor(ranef(coxmm)$strain[,2], ranef(coxmm2)$strain[,2]))
			
		}
		
		total_results[[j]]<-cbind(raw_CV_np_md, raw_CV_np_hr, raw_CV_pp_md, raw_CV_pp_hr, p/reps, raw_corrs_np_md, raw_corrs_np_hr, raw_corrs_pp_md, raw_corrs_pp_hr, ns[j], sigma_G[g])
		
	}
	
	full_results[[g]]<-total_results

}



# Parcel Up and save
output<-as.data.frame(t(array(unlist(full_results), c(11,length(ns)*length(sigma_G)))))
names(output)<-c("np_md_CV", "np_hr_CV", "pp_md_CV", "pp_hr_CV", "power", "np_md_r", "np_hr_r", "pp_md_r", "pp_hr_r", "n", "sigma_G")
write.table(output, file=paste0(job, "_Sim_raw.csv"), sep=",", row.names=F, col.names=names(output))

