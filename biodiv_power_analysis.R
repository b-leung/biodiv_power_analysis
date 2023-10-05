#biodiversity power analysis: coded by Brian Leung - Oct 4, 2023

d=readRDS("./summary_orig.rds") #summary info from original data - taking systems with < -0.015 decline, and >10 populations (using Bayesian Hierarchical Mixture Model from Leung et al. 2020. Nature 588: 267-271
ns=7 #number of sample years

opt=1 # 0 = increase sampling to increase power to detect declines under status quo conditions (i.e., no policy); 1 = power analysis to detect improvements under policy intervention; 2) power analysis to detect that trends are better than some fixed "reference" thresholds under policy intervention

thin=30 #thin the MCMC for speed
num_samp=1000 # number of points to numerically integrate realized distribution, given each MCMC "reality"

x=rnorm(num_samp) #standard normal - points for numerical integration (will be scaled by standard deviation).
i=0
nnp=c(128)
prob_thresh=c(.7) #Set certainty threshold (essentially fraction of posterior distribution)

severity_thresh=c(-0.015) #set a threshold of severity for comparison (less than this point under status quo; greater than this point under policy intervention). 

#file names to save associated with the 3 options 
nm_out=c("power_status_quo_decline.rds", "power_detect_improvement.rds", "power_ref_threshold.rds")	





if(opt==0)
{
	#power to detect decline
	newdat=as.data.frame(matrix(NA,nrow=nrow(d)*length(nnp)*length(severity_thresh)*length(prob_thresh),ncol=10))
	names(newdat)=c("system","j","npops","m_theta","sd_theta", "nyr","nnewpop","severity","prob_thresh","power")
}else{
	#power to detect change under policy
	#the question might be, how often would you expect to detect a change relationship? 

	policy_theta=0 #the true growth rate, after policy (0 is true stability)	
	
	newdat=as.data.frame(matrix(NA,nrow=nrow(d)*length(nnp)*length(prob_thresh),ncol=10))
	names(newdat)=c("system","j","npops","m_theta","sd_theta","policy_theta","nyr","nnewpop","prob_thresh","power")

}
for(j in 1:nrow(d))
{
	nm=d$system[j]
	nm1=paste("MCMC_runs/MCMC_",gsub(" ","_",nm),".rds",sep="") #read in MCMC 
	f<-readRDS(nm1)
	s=d$pop_fluc_sd[j]
	if(opt<2)
	{	#include real uncertainty
		theta1=mean(f$theta1) #mean and variation of theta's - i.e., the posterior distribution.
		s1=sd(f$theta1)  
	}else{ #compare new populations to a fixed threshold of decline.
		theta1=severity_thresh
		s1=0
	}
	
	mcs=seq(1,length(f$theta1),by=thin) #choose theta, tau from MCMC runs
	len_mcs=length(mcs)

	for(np in nnp) #iterate through each sample size: the probability of detecting a difference from some threshold, given the system as is.
	{
		print(np)
		print(Sys.time())
		s0=((f$tau1[mcs]^2+(s^2)/ns)/np)^.5 #this gives the variance generating new data (i.e., sample means)
		if(opt==0)
		{
			x1=(s0 %*% t(x)) +f$theta1[mcs] #the generation of sample means of populations.
	#		#generates a matrix nrow MCMC, ncol x
			s2=1/(1/s1^2+1/s0^2) #the joint variation - variance - vector has length of s0
			t2=s2*(theta1/s1^2+x1/s0^2) #the mean, given sample mean from x1. From analytic Bayesian - t2 and s2 define the posterior distribution for theta.
			for(thresh in severity_thresh)
			{
				#return back a matrix of Z values. From this, can do all the other statistics on it
				Z=(thresh-t2)/s2^.5 
				for(pt in prob_thresh)
				{
					i=i+1
					#what fraction of realizations (posterior distr) is beyond the threshold probability
					zt=qnorm(pt) 
					power=length(which(Z>zt))/length(Z) #the fraction of realizations, with posterior probability mass beyond threshold of certainty  
					newdat[i,1]=nm
					newdat[i,-1]=c(j,d$pop[j],theta1,s1,ns,np,thresh,pt,power)
			
				}
			}

		}else{ 
		#power to detect improvement post-policy  
		#To calculate change - instead of f$theta1, use alternative hypothesis A to generates a distribution of xbars

			x1=s0 %*% t(x)+policy_theta #this generates xbar, around the alternative hypothesis mean, calculated for each tau; this dimension is s0 (rows) by x (columns)
			theta2=x1-theta1 #so how much does theta change? 
			s2=s0^2+s1^2 #sum variance law
			Z=theta2/s2^.5 #this will allow posterior credible interval to be calculated
			#for each x1, a posterior distribution is generated (defined by theta2 and s2) 
			for(pt in prob_thresh)
			{
				i=i+1
				#what fraction of realizations (posterior distr) is beyond the certainty probability (pt) that an improvement has occurred 
				zt=qnorm(pt) #convert to z threshold - we want to see if policy has positive effect 
				power=length(which(Z>zt))/length(Z) 
				#this give the positions of the matrix that are greater than a value.
				#sum the frequency of those values occurring (pdf1). Since summing across all MCMC, divide by the number of MCMCs
				newdat[i,1]=nm
				newdat[i,-1]=c(j,d$pop[j],theta1,s1,policy_theta,ns,np,pt,power)

			}

		}
	}	
		
}
saveRDS(newdat,nm_out[opt+1])


