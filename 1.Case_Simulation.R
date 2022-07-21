
# Clean up
rm(list=ls())

# Load packages
library(simsurv)
library(plyr)
library(lmerTest)
library(coxme)
library(survminer)
library(mvtnorm)
library(doSNOW)


# Set the wd
setwd("/Users/alistairsenior/OneDrive - The University of Sydney (Staff)/Simulation Genetic Mapping/Simulation_Submitted")

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

# Start by simulating data in one strain no genetic effect to show survival there using paramterised a and b, n per group 400. Noting lnHR = -1 from Nakagawa et al.
set.seed(1234)
dat<-my_sim(n_strain=1, n=400, lnHR=-0.5, sigma_IntG=0, sigma_Tr=0, cor_IntTrt=0)[[1]]
fit<-survfit(Surv(eventtime/7, status) ~ trts, data=dat)

# Plot a boxplot + violin
pdf("F1.A.pdf", height=8.2/3, width=8.2/3)
p1<-ggplot(dat, aes(x=eventtime/7, col=diet, fill=diet)) +
	geom_density(alpha=0.2) + theme_bw() + xlab("Age at Death (Weeks)") + labs(title="A") + ylab("Density") +
	theme(legend.position=c(0.23, 0.7)) + xlim(0, 200)
print(p1)
dev.off()

# plot a KM
pdf("F1.B.pdf", height=8.2/3, width=8.2/3)
p2<-ggsurvplot(fit, legend="none", censor=F, title="B", xlab = "Time (Weeks)", ylab="Survivorship")
print(p2)
dev.off()

# Check median survival time and effect
meds<-ddply(dat, .(trts), summarise, median(eventtime), mean(eventtime))
meds
meds[2,2]/meds[1,2]

# Find the modes
dat1<-dat[which(dat$diet == "DR"),]
dat1$eventtime<-round(dat1$eventtime)
mode<-ddply(dat1, .(round(eventtime)), summarise, length(eventtime))
mode[order(mode$..1),]
dat1<-dat[which(dat$diet == "AL"),]
dat1$eventtime<-round(dat1$eventtime)
mode<-ddply(dat1, .(round(eventtime)), summarise, length(eventtime))
mode[order(mode$..1),]

# Add in a histrogram of the mortality data in mice.
Rikke_data<-read.csv("Rikke.csv")
pdf("F1.C.pdf", height=8.2/3, width=8.2/3)
p3<-ggplot(Rikke_data, aes(x=Lifespan/7)) + geom_histogram() + xlab("Age at Death (Weeks)") + labs(title="C") + theme_bw() + ylab("Count")
print(p3)
dev.off()

# Lets simulate some strains, but assume no genetic variance
set.seed(567)
dat<-my_sim(n_strain=40, n=5, lnHR=-0.5, sigma_IntG=0, sigma_Tr=0, cor_IntTrt=0)[[1]]
dat$trts_cat<-as.factor(dat$trts)

# Now get the strain by strain effect sizes - calculated as lm (t-tests, basically)
effects<-data.frame(strains = seq(1, 40, 1), ES=0, SE=0, p=0)
for(i in 1:nrow(effects)){
	data_i<-dat[which(dat$strain == i),]
	test<-lm(eventtime ~ diet, dat=data_i)
	effects$ES[i]<-test$coef[2]
	effects$SE[i]<-summary(test)$coeff[2,2]
	effects$p[i]<-summary(test)$coeff[2,4]
}
effects$signif<-""
effects$signif[which(effects$p <= 0.05)]<-"*"
effects<-effects[order(effects$ES, decreasing = T),]
range(effects$ES)

# Make waterfall plots
pdf("F1.D.pdf", height=8.2/3, width=6)
p3<-ggplot(effects, aes(x=seq(1, 40, 1), y=ES)) + geom_point() + theme_bw() +
	theme(axis.text.x=element_text(angle=90, size=10), axis.ticks.x=element_blank(), panel.grid.minor=element_blank()) + 
	labs(title="D", x = "Strain ID", y = "Diff. Lifespan [DR - AL] - LM") + 
	geom_hline(yintercept=0) + ylim(-300, 300) + 
	annotate("text", x=seq(1, 40, 1), y=280, label=effects$signif, size=10) + 
	scale_x_continuous(breaks=seq(1, 40, 1), labels=effects$strains)
print(p3)
dev.off()

# Let's repeat to see what happens
dat<-my_sim(n_strain=40, n=5, lnHR=-0.5, sigma_IntG=0, sigma_Tr=0, cor_IntTrt=0)[[1]]

# Now get the strain by strain effect sizes - calculated as lm and coxph
effects<-data.frame(strains = seq(1, 40, 1), ES=0, SE=0, p=0, ES_hr=0)
for(i in 1:nrow(effects)){
	data_i<-dat[which(dat$strain == i),]
	test<-lm(eventtime ~ diet, dat=data_i)
	effects$ES[i]<-test$coef[2]
	effects$SE[i]<-summary(test)$coeff[2,2]
	effects$p[i]<-summary(test)$coeff[2,4]
	test<-coxph(Surv(eventtime, status) ~ trts, data=data_i)
	effects$ES_hr[i]<-test$coef
}
effects$signif<-""
effects$signif[which(effects$p <= 0.05)]<-"*"
effects<-effects[order(effects$ES, decreasing = T),]

# Analyses
sd(effects$ES) 
mean(effects$ES) 
sd(effects$ES) / mean(effects$ES)
sd(effects$ES_hr)
mean(effects$ES_hr)
sd(effects$ES_hr) / abs(mean(effects$ES_hr))

# What if we use shrinkage estimate of the variance in effects
lmm<-lmer(eventtime ~ diet + (1+diet|strain), data=dat)
summary(lmm)
# Random effects:
 # Groups   Name        Variance Std.Dev. Corr 
 # strain   (Intercept)  1182.6   34.39        
          # dietDR        988.5   31.44   -1.00
 # Residual             22465.5  149.88        
# Number of obs: 400, groups:  strain, 40

# Fixed effects:
            # Estimate Std. Error     df t value Pr(>|t|)    
# (Intercept)   820.53      11.91  39.30  68.883  < 2e-16 ***
# dietDR         86.39      15.79  93.14   5.471 3.76e-07 ***
31.44 / 86.39


# What if we use coxme - the right model AND shrinkage
coxmm<-coxme(Surv(eventtime, status) ~ trts + (1+trts|strain), data=dat)
summary(coxmm)
# Fixed coefficients
           # coef exp(coef)  se(coef)    z       p
# trts -0.5882912 0.5552754 0.1032758 -5.7 1.2e-08

# Random effects
 # Group  Variable  Std Dev       Variance      Corr         
 # strain Intercept  8.188634e-03  6.705372e-05 -8.394434e-01
        # trts       6.865365e-03  4.713324e-05              

8.188634e-03 / 0.5882912

# Extract the estimates from lmm and coxme and add to the effects for plotting
lmm_ES<-ranef(lmm)$strain
effects$lmm<-lmm_ES[match(effects$strains, row.names(lmm_ES)), 2] + summary(lmm)$coef[2,1]
cox_ES<-ranef(coxmm)$strain
effects$cox<-cox_ES[match(effects$strains, row.names(cox_ES)), 2] + coxmm$coeff

# Make waterfall plots
pdf("F1.E.pdf", height=8.2/3, width=6)
p4<-ggplot(effects, aes(x=seq(1, 40, 1), y=ES)) + geom_point() + theme_bw() +
	theme(axis.text.x=element_text(angle=90, size=10), axis.ticks.x=element_blank(), panel.grid.minor=element_blank()) + 
	labs(title="E", x = "Strain ID", y = "MD Lifespan [DR - AL]") + 
	geom_hline(yintercept=0) +
	annotate("text", x=seq(1, 40, 1), y=280, label=effects$signif, size=10) + 
	scale_x_continuous(breaks=seq(1, 40, 1), labels=effects$strains)	 +
	geom_point(aes(x=seq(1, 40, 1), y=lmm), inherit.aes=F, col="cornflowerblue") +
	scale_y_continuous(limits=c(-300, 300))
print(p4)
dev.off()


# Make waterfall plots of lnHR
effects<-effects[order(effects$ES_hr),]
pdf("F1.F.pdf", height=8.2/3, width=6)
p5<-ggplot(effects, aes(x=seq(1, 40, 1), y=ES_hr)) + geom_point() + theme_bw() +
	theme(axis.text.x=element_text(angle=90, size=10), axis.ticks.x=element_blank(), panel.grid.minor=element_blank()) + 
	labs(title="F", x = "Strain ID", y = "lnHR [DR / AL]") + 
	geom_hline(yintercept=0) +
	scale_x_continuous(breaks=seq(1, 40, 1), labels=effects$strains) + 
	geom_point(aes(x=seq(1, 40, 1), y=cox), inherit.aes=F, col="cornflowerblue") + 
	scale_y_continuous(limits=c(-3, 3))
		
print(p5)
dev.off()
