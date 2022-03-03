
# Clean up
rm(list=ls())

# Load the packages
library(coxme)
library(plyr)
library(lmerTest)
library(ggplot2)
library(MortalityLaws)

# Rikke et al. data
Rikke_data<-read.csv("Rikke.csv")
head(Rikke_data)
Rikke_data$Strain<-as.factor(Rikke_data$Strain)
Rikke_data$trt<-as.numeric(as.factor(Rikke_data$Diet)) - 1
Rikke_data$Status<-1

# No pooling MD
effects<-data.frame(strain=sort(unique(Rikke_data$Strain)), lm=NA, cm=NA)
for(i in 1:nrow(effects)){
	data<-Rikke_data[which(Rikke_data$Strain == effects$strain[i]),]
	effects$lm[i]<-lm(Lifespan ~ trt, data=data)$coef[2]
	effects$cm[i]<-coxph(Surv(Lifespan, Status) ~ trt, data=data)$coef[1]
}

# Fit the mixed-models
lmm<-lmer(Lifespan ~ trt + (1 + trt|Strain), data=Rikke_data)
summary(lmm)
lmm0<-lmer(Lifespan ~ trt + (1|Strain), data=Rikke_data)
summary(lmm0)
anova(lmm0, lmm)

cmm<-coxme(Surv(Lifespan, Status) ~ trt + (1 + trt|Strain), data=Rikke_data)
summary(cmm)
cmm0<-coxme(Surv(Lifespan, Status) ~ trt + (1|Strain), data=Rikke_data)
anova(cmm0, cmm)


# Get the strain specific estimates
effects$lmm<-ranef(lmm)$Strain$trt + lmm@beta[2]
effects$cmm<-ranef(cmm)$Strain[,2] + cmm$coef[1]

# Make the waterfall plots
pdf("F3.AB.pdf", height=4, width=7)

effects<-effects[order(effects$lm, decreasing=T),]
p1<-ggplot(effects, aes(x=seq(1, nrow(effects), 1), y=lm)) + geom_point() + theme_bw() +
	theme(axis.text.x=element_text(angle=90, size=10), axis.ticks.x=element_blank(), panel.grid.minor=element_blank()) + 
	labs(title="A", x = "Strain ID", y = "MD Lifespan [DR - AL]") + 
	geom_hline(yintercept=0) +
	scale_x_continuous(breaks=seq(1, nrow(effects), 1), labels=effects$strain)	 +
	geom_point(aes(x=seq(1, nrow(effects), 1), y=lmm), inherit.aes=F, col="cornflowerblue") +
	scale_y_continuous(limits=c(-600, 600))
print(p1)

effects<-effects[order(effects$cm),]
p1<-ggplot(effects, aes(x=seq(1, nrow(effects), 1), y=cm)) + geom_point() + theme_bw() +
	theme(axis.text.x=element_text(angle=90, size=10), axis.ticks.x=element_blank(), panel.grid.minor=element_blank()) + 
	labs(title="B", x = "Strain ID", y = "lnHR [DR / AL]") + 
	geom_hline(yintercept=0) +
	scale_x_continuous(breaks=seq(1, nrow(effects), 1), labels=effects$strain)	 +
	geom_point(aes(x=seq(1, nrow(effects), 1), y=cmm), inherit.aes=F, col="cornflowerblue") +
	scale_y_continuous(limits=c(-5, 5))
print(p1)

dev.off()

# Here the CVG from the no pooling is
sd(effects$ES) 
mean(effects3$ES)
202  / abs(-5.14)

# From the lmm
183 / abs(-5.17)

# Here we can get no pooling HR
sd(effects3$ES2) 
mean(effects3$ES2)
1.55  / abs(-0.26)

# From the cmm
1.36 / abs(-0.31)

# KM analysis from Rikke et al thinking about proportionality
pdf("FigS2A.pdf")

# Whoah - massive time-dependent effect!
plot(survfit(Surv(Lifespan, Status) ~ trt, data= Rikke_data), col=c(1,2), lwd=2, xlab="Time (Days)", ylab="Survivorship")
legend("topright", c("AL", "DR"), col=c(1,2), lty=1, lwd=2)
abline(v=820)
mtext("A.", at = -50, line=0.5, cex=1.5)

dev.off()

# Can we get a prediction for the mean age at death in a strain

# From Survfit on the raw data
fit<-survfit(Surv(Lifespan, Status) ~ trt, data=Rikke_data)
BL2<-cbind(fit$time, fit$n.risk, fit$n.event, fit$surv, fit$cumhaz)[c(1:298),]

# Try and recalculate
my_att<-ddply(Rikke_data[which(Rikke_data$Diet == "AL"),], .(Lifespan), summarise, n.event=length(Lifespan))
plot(my_att$Lifespan, BL2[,1])
abline(a=0, b=1)
plot(my_att$n.event, BL2[,3])
abline(a=0, b=1)

my_att$n.risk<-(dim(Rikke_data[which(Rikke_data$Diet == "AL"),])[1] - cumsum(my_att$n.event)) + my_att$n.event
plot(my_att$n.risk, BL2[,2])
abline(a=0, b=1)

my_att$haz<-my_att$n.event/my_att$n.risk
my_att$cum_haz<-cumsum(my_att$haz)
plot(my_att$cum_haz, BL2[,5])
abline(a=0, b=1)

my_att$surv<-1 - (cumsum(my_att$n.event) / my_att$n.risk[1])
plot(my_att$surv, BL2[,4])
abline(a=0, b=1)

# OK - easy - know how that is all done

# One option for the BL based on the AL raw data
BL_lt<-LifeTable(x=BL2[,1], lx=BL2[,4])$lt
#mean(Rikke_data[which(Rikke_data$Diet == "AL"),]$Lifespan)
#BL_lt$ex[1] + BL_lt$x[1]
# That looks like a good BL, which generates mean survival close to the whole AL group

# Get a list of the lifetables for each strain and trt
lifespans<-data.frame(Strain=rownames(ranef(cmm)$Strain), AL=NA, DR=NA)
AL_cum_lx<-rep(0, nrow(BL_lt))
DR_cum_lx<-rep(0, nrow(BL_lt))
DR_l95_lx<-rep(0, nrow(BL_lt))
DR_u95_lx<-rep(0, nrow(BL_lt))

for(i in 1:nrow(effects)){
	
	# Pull out the effects for i
	effects_i<-ranef(cmm)$Strain[i,]
	
	# Calculate the hazards for i and get the lifetable
	AL_strain_i<-BL_lt$qx * exp(effects_i[1])
	AL_lt_i<-LifeTable(x=BL_lt$x, qx=AL_strain_i)$lt
	
	# Calculate the hazards for i
	DR_strain_i<-BL_lt$qx * exp(effects_i[1] + cmm$coeff + effects_i[2])
	DR_lt_i<-LifeTable(x=BL_lt$x, qx=DR_strain_i)$lt
	
	# Add in the mean life expectancies
	lifespans$AL[i]<-AL_lt_i$ex[1] + AL_lt_i$x[1]
	lifespans$DR[i]<-DR_lt_i$ex[1] + DR_lt_i$x[1]
	
	# get the weighted average ridix
	AL_wi<-length(which(Rikke_data$trt == 0 & Rikke_data$Strain == row.names(ranef(cmm)$Strain)[i])) / length(which(Rikke_data$trt == 0))
	DR_wi<-length(which(Rikke_data$trt == 1 & Rikke_data$Strain == row.names(ranef(cmm)$Strain)[i])) / length(which(Rikke_data$trt == 1))
	AL_cum_lx<-AL_cum_lx + AL_lt_i$lx * AL_wi
	DR_cum_lx<-DR_cum_lx + DR_lt_i$lx * DR_wi

	# repeat for the low 95		
	DR_strain_i<-BL_lt$qx * exp(effects_i[1] + cmm$coeff - (1.96 * sqrt(vcov(cmm))) + effects_i[2])
	DR_lt_i<-LifeTable(x=BL_lt$x, qx=DR_strain_i)$lt
		
	# get the weighted average ridix
	DR_l95_lx<-DR_l95_lx + DR_lt_i$lx * DR_wi
	
	# repeat for the upper 95		
	DR_strain_i<-BL_lt$qx * exp(effects_i[1] + cmm$coeff + (1.96 * sqrt(vcov(cmm))) + effects_i[2])
	DR_lt_i<-LifeTable(x=BL_lt$x, qx=DR_strain_i)$lt
		
	# get the weighted average ridix
	DR_u95_lx<-DR_u95_lx + DR_lt_i$lx * DR_wi


}

# Plot it
pdf("FigS2B.pdf")

plot(BL_lt$x, AL_cum_lx/AL_cum_lx[1], xlim=c(0, max(Rikke_data$Lifespan)), type="l", lwd=2, xlab="Time (Days)", ylab="Survivorship")
lines(BL_lt$x, DR_cum_lx/DR_cum_lx[1], col=2, lwd=2)
lines(BL_lt$x, DR_l95_lx/DR_l95_lx[1], col=2, lty=2, lwd=2)
lines(BL_lt$x, DR_u95_lx/DR_u95_lx[1], col=2, lty=2, lwd=2)
abline(v=820)
mtext("B.", at = -50, line=0.5, cex=1.5)

dev.off()
