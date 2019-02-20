#################################
### simulation_model_function ###
#################################

library(simulateGP)
library(TwoSampleMR)

###########################
## Effect model function ##
###########################

## Gives the snp_exp_effect and plei_effect based on the number of snps, snp exposure variance, pleitropic variance and 
## whether the pleiotropy has no effect, effects the snp exposure independent of size, or effects the snp exposure dependent on size
## i.e. larger pleiotropy effects correspond to smaller snp exposure effects
## This is based on the Omnigenic model, and thinking that smaller exposure effect sizes are more likely to be pleiotropic and together add
## up to have an effect on the outcome.
## Noise has been added to the "dependent" model to make it more realistic whilst still having small effect sizes with large pleitropic effects

effect_model <- function(nsnp, snp_exp_var, plei_var, model=c("none", "independent", "dependent")[1])
{ 
  snp_exp_effect <- abs(choose_effects(nsnp, snp_exp_var)) 
  
  if(model == "dependent")
  {
    plei_effect <- abs(choose_effects(nsnp, plei_var))
    snp_exp_effect <- snp_exp_effect[order(snp_exp_effect, decreasing=FALSE)]
    plei_effect <- plei_effect[order(plei_effect, decreasing=TRUE)]
    noise_effect1 <- choose_effects(nsnp, snp_exp_var)
    noise_effect1 <- noise_effect1/10
    noise_effect2 <- choose_effects(nsnp, plei_var)
    noise_effect2 <- noise_effect2/10
    snp_exp_effect <- rowSums(cbind(snp_exp_effect,noise_effect1))
    plei_effect <- rowSums(cbind(plei_effect,noise_effect2))
  } else if(model == "independent") {
    plei_effect <- abs(choose_effects(nsnp, plei_var))
  } else if (model == "none") {
    plei_effect <- choose_effects(nsnp, 0)
  } else {
    stop("Must specify none, independent or dependent")
  }
  return(list(plei_effect = plei_effect, snp_exp_effect = snp_exp_effect))
} 




###############################
## Simulation model function ##
###############################

## Gives the MR results and heterogeneity test for both one sample and two sample MR based on a theoretical dataset that we simualte
## Simulated dataset is based on number of people(nsample), number of snps found in hypothetical GWAS that relate to exposure,
## snp exposure variance, pleiotropic variance, counfounder exposure variance, exposure outcome variance, counfounder outcome variance
## and the pleitropic model (same as in effect model)

## The "runsim" function uses the "effect_model" function within it so only need to run the "runsim" function to get the results

runsim <- function (nsample, nsnp, snp_exp_var, plei_var, conf_exp_var, exp_out_var, conf_out_var, 
                    pleiotropy_model=c("none", "independent", "dependent"))
{
  g<- make_geno(nsample, nsnp, 0.5)
  g1 <- subset(g[1:(nsample/2),])
  g2 <- subset(g[((nsample/2)+1):nsample,])
  conf1 <- rnorm((nsample/2))
  conf2 <- rnorm((nsample/2))
  pl_effs <- effect_model(nsnp, snp_exp_var, plei_var, model=pleiotropy_model)
  x1 <- make_phen(c(pl_effs$snp_exp_effect, conf_exp_var), cbind(g1, conf1))
  x2 <- make_phen(c(pl_effs$snp_exp_effect, conf_exp_var), cbind(g2, conf2))
  y1 <- make_phen(c(pl_effs$plei_effect, exp_out_var, conf_out_var), cbind(g1, x1, conf1))
  dat1 <- get_effs(x1,y1,g1)
  s1 <- subset(dat1, pval.exposure < 5e-8)
  res1 <- mr(s1, method_list = "mr_ivw") ## one sample MR 
  y2 <- make_phen(c(pl_effs$plei_effect, exp_out_var, conf_out_var), cbind(g2, x2, conf2))
  dat2 <- get_effs(x2,y2,g2)
  s2 <- subset(dat2, SNP %in% s1$SNP)
  res2 <- mr(s2, method_list = "mr_ivw")
  het1 <- mr_heterogeneity(s1, method_list="mr_ivw")
  het1 <- het1[6:8]
  het2 <- mr_heterogeneity(s2, method_list="mr_ivw")
  het2 <- het2[6:8]
  mr_results1 <- cbind(res1, het1)
  mr_results2 <- cbind(res2, het2)
  mr_results1$type <-"one sample"
  mr_results2$type <-"two sample"
  mr_results <- rbind(mr_results1, mr_results2)
  return(mr_results)
}   


########################
## Examples of runsim ##
########################

a <- runsim(40000, 100, 0.2, 0.05, 0.1, 0.2, 0.2, "none")
b <- runsim(40000, 100, 0.2, 0.05, 0.1, 0.2, 0.2, "independent")
c <- runsim(40000, 100, 0.2, 0.05, 0.1, 0.2, 0.2, "dependent")



###############
### runsim2 ###
###############

## The runsim2 function does exactly the same as the runsim function, however the output is only the heterogenetiy.
## This is so we can check that the pleiotropy is working in the way we expect it to

runsim2 <- function (nsample, nsnp, snp_exp_var, plei_var, conf_exp_var, exp_out_var, conf_out_var, 
                     pleiotropy_model=c("none", "independent", "dependent"))
{
  g<- make_geno(nsample, nsnp, 0.5)
  g1 <- subset(g[1:(nsample/2),])
  g2 <- subset(g[((nsample/2)+1):nsample,])
  conf1 <- rnorm((nsample/2))
  conf2 <- rnorm((nsample/2))
  pl_effs <- effect_model(nsnp, snp_exp_var, plei_var, model=pleiotropy_model)
  x1 <- make_phen(c(pl_effs$snp_exp_effect, conf_exp_var), cbind(g1, conf1))
  x2 <- make_phen(c(pl_effs$snp_exp_effect, conf_exp_var), cbind(g2, conf2))
  y1 <- make_phen(c(pl_effs$plei_effect, exp_out_var, conf_out_var), cbind(g1, x1, conf1))
  dat1 <- get_effs(x1,y1,g1)
  s1 <- subset(dat1, pval.exposure < 5e-8)
  res1 <- mr(s1, method_list = "mr_ivw") ## one sample MR 
  y2 <- make_phen(c(pl_effs$plei_effect, exp_out_var, conf_out_var), cbind(g2, x2, conf2))
  dat2 <- get_effs(x2,y2,g2)
  s2 <- subset(dat2, SNP %in% s1$SNP)
  res2 <- mr(s2, method_list = "mr_ivw") ## Two sample MR
  
  het1 <- mr_heterogeneity(s1, method_list="mr_ivw")
  het2 <- mr_heterogeneity(s2, method_list="mr_ivw")
  
  res1$type <- "one sample"
  res2$type <- "two sample"
  
  het1$type <- "one sample"
  het2$type <- "two sample"
  
  mr_results <- rbind(res1=res1, res2=res2)
  
  heterogeneity <- rbind(het1=het1, het2=het2)
  
  return(list( heterogeneity=heterogeneity))
  
}   


################################
## Checking the heterogeneity ##
################################

## The following is to check for heterogeneity in the dependent model to see how often we get significant heterogeneity 
## This was done because it seemed that there was not always evidence of heterogeneity when using the "dependent" model
## it seems that we get heterogenity about 70% of the time with the "dependent" model
## Is this enough??


n.times <- 100
result <- numeric(n.times) 
Q <- numeric(n.times)
for (i in 1:n.times) {
  result[i] <- runsim2(40000, 100, 0.2, 0.05, 0.1, 0.2, 0.2, "dependent")
  Q [i] <- rbind(result[[i]][8])
}

Q <- data.frame(matrix(unlist(Q), nrow=n.times, byrow=T))
colnames(Q)<-c("one_sample", "two_sample")
table(Q[,1] < 0.05) # One sample
table(Q[,2] < 0.05) # Two sample
table(Q <0.05) # one and two sample




#########################
## Multiple parameters ##
#########################

## This allows us to use the "runsim" function with multiple parameters for each variable
## need to be aware of which parameters work together....

params <- expand.grid(
  nsample = c(40000),
  nsnp = c(100),
  snp_exp_var = c(0.2, 0.21),
  plei_var = c(0.05, 0.051),
  conf_exp_var = c(0.1, 0.11),
  exp_out_var = c(0.2, 0.21),
  conf_out_var = c(0.2, 0.21),
  pleiotropy_model = c("none", "independent", "dependent"),
  stringsAsFactors = FALSE
)



## Output 
out <- list()
for(i in 1:nrow(params))
{
  out[[i]] <- runsim(params$nsample[1], params$nsnp[i], params$snp_exp_var[i], params$plei_var[i], params$conf_exp_var[i], 
                     params$exp_out_var[i], params$conf_out_var[1], params$pleiotropy_model[i])
}


