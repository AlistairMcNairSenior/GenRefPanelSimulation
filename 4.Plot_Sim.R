
# Clean up
rm(list=ls())

# Load libraries
library(ggplot2)
library(plyr)

# Read in the data
files<-dir("GenSim")[-1001]
data<-read.csv(paste0("GenSim/", files[1]))
for(i in 2:length(files)){
	data<-rbind(data, read.csv(paste0("GenSim/", files[i])))	
}
dim(data)

# Convert CVG to bias
data$np_md_CV<-data$np_md_CV - (data$sigma_G / 0.5)
data$pp_md_CV<-data$pp_md_CV - (data$sigma_G / 0.5)
data$np_hr_CV<-data$np_hr_CV - (data$sigma_G / 0.5)
data$pp_hr_CV<-data$pp_hr_CV - (data$sigma_G / 0.5)

# Avergaes for each parameter set\
data$param.id<-paste0(data$n, "_", data$sigma_G)
head(data)

data_means<-ddply(data, .(param.id), summarise, n=n[1], sigma_G=sigma_G[1]
	, m.np_md_CV=median(np_md_CV), m.np_hr_CV=median(np_hr_CV), m.pp_md_CV=median(pp_md_CV), m.pp_hr_CV=median(pp_hr_CV)
	, lower_1=quantile(np_md_CV, 0.25), lower_2=quantile(pp_md_CV, 0.25)
	, lower_3=quantile(np_hr_CV, 0.25), lower_4=quantile(pp_hr_CV, 0.25)
	, upper_1=quantile(np_md_CV, 0.75), upper_2=quantile(pp_md_CV, 0.75)
	, upper_3=quantile(np_hr_CV, 0.75), upper_4=quantile(pp_hr_CV, 0.75)
	, power=mean(power), np_md_r=mean(np_md_r), np_hr_r=mean(np_hr_r), pp_md_r=mean(pp_md_r, na.rm=T), pp_hr_r=mean(pp_hr_r))


head(data_means)

pdf("Sim.pdf", height=4, width=4)

sigma_G<-unique(data$sigma_G)

for(g in 1:length(sigma_G)){
	
	# Pull out the right data and scale r
	G1<-data_means[which(data_means$sigma_G == sigma_G[g]),]
	scale_cor<-G1[,c(17:20)] * 3.5
	names(scale_cor)<-c("r1", "r2", "r3", "r4")
	G1<-cbind(G1, scale_cor)
	G1$width<-0
	
	p1<-ggplot(G1, aes(x=n-0.3, y= m.np_md_CV)) +
		geom_point() + geom_errorbar(aes(ymin=lower_1, ymax=upper_1, width=0)) + theme_bw() + geom_hline(yintercept=0) + 
		labs(x = "n per Strain per Diet", y = expression(Bias~~CV[G]~~MD), title=chartr("1234", "ABCD", g)) + 
		geom_point(aes(x=n+0.3, y=m.pp_md_CV), col="cornflowerblue") + geom_errorbar(aes(x=n+0.3, ymin=lower_2, ymax=upper_2, width=0), col="cornflowerblue") +
		geom_point(aes(x=n-0.3, y= r1), shape=4, size=4) +
		geom_point(aes(x=n+0.3, y= r2), shape=4, size=4, col="cornflowerblue")  +
		scale_y_continuous(lim = c(-1, 3.5),  sec.axis=sec_axis(trans= ~./3.5, name = "Replicability", breaks=seq(-1, 1, 0.2))) +
		annotate("text", x = 35, y = -0.9, label=substitute(CV[G]==x, list(x=sigma_G[g]/0.5)), size=5, alpha=0.4)
	print(p1)

}


for(g in 1:length(sigma_G)){
	
	G1<-data_means[which(data_means$sigma_G == sigma_G[g]),]
	scale_cor<-G1[,c(17:20)] * 3.5
	names(scale_cor)<-c("r1", "r2", "r3", "r4")
	G1<-cbind(G1, scale_cor)
	
	p1<-ggplot(G1, aes(x=n-0.3, y= m.np_hr_CV)) +
		geom_point() + geom_errorbar(aes(ymin=lower_3, ymax=upper_3, width=0)) + theme_bw() + geom_hline(yintercept=sigma_G[1]/0.5) + 
		labs(x = "n per Strain per Diet", y = expression(Bias~~CV[G]~~lnHR), title=chartr("1234", "EFGH", g)) + 
		geom_point(aes(x=n+0.3, y=m.pp_hr_CV), col="cornflowerblue") + geom_errorbar(aes(x=n+0.3, ymin=lower_4, ymax=upper_4, width=0), col="cornflowerblue") +
		geom_point(aes(x=n-0.3, y= r3), shape=4, size=4) +
		geom_point(aes(x=n+0.3, y= r4), shape=4, size=4, col="cornflowerblue")  +
		scale_y_continuous(lim = c(-1, 3.5),  sec.axis=sec_axis(trans= ~./3.5, name = "Replicability", breaks=seq(-1, 1, 0.2))) +
		annotate("text", x = 35, y = -0.9, label=substitute(CV[G]==x, list(x=sigma_G[g]/0.5)), size=5, alpha=0.4)
	print(p1)

}

dev.off()


# Check out some power
G1<-data_means[which(data_means$sigma_G == sigma_G[3]),]
G1


