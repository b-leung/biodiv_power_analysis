There are several files included.

Data files 

info_2022.rds 
growth_rate_2022.rds 

these are required to run the Bayesian Hierarchical Mixture Model from Leung et al. 2020. Nature 588: 267-271 (see https://doi.org/10.5281/zenodo.3901586 to obtain the code. If used please also cite original paper). 

In script BHM.R, change the code to the file names above (i.e., from "info.rds" and "growth_rates_example.rds", in the BHM.R code). 

These will generate MCMC runs for each country-taxon system, which will need to be saved;
i.e., in function do_fit, include the following lines at the end of the code, to save them into folder MCMC_runs:

		mat=as.data.frame(cbind(theta,tau,frac))
		names(mat)=c(paste("theta",1:K,sep=""),paste("tau",1:K,sep=""),paste("frac",1:K,sep=""))
		saveRDS(mat,paste("./MCMC_runs/MCMC_",gsub(" ","_",u[rep]),".rds",sep="")) 



Alternatively, we directly provide the MCMC runs in the folder "MCMC_runs", for 62 systems with >10 populations and rapid declines (posterior theta < -0.015). 

The R script "biodiversity_power_analysis.R" will draw from these files to run the biodiversity power analysis.


For questions, please contact: Brian Leung (brian.leung2@mcgill.ca)
