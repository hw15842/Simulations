#################################
### simulation_model_function ###
#################################

library(simulateGP)
library(TwoSampleMR)
library(plyr)

###########################
## Effect model function ##
###########################

## Gives the snp_exp_effect and plei_effect based on the number of snps, snp exposure variance, pleitropic variance and 
## whether the pleiotropy has 1) no effect, 2) effects the snp exposure independent of size, or 3) effects the snp exposure dependent on size
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

## Gives the MR results and heterogeneity test for both one sample and two sample MR based on a theoretical dataset that we simulate
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
  snp_exp_var = c(0.2, 0.5),
  plei_var = c(0.05, 0.1),
  conf_exp_var = c(0.1, 0.3),
  exp_out_var = c(0.2, 0.5),
  conf_out_var = c(0.2, 0.5),
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

out<- ldply(out, data.frame)


#################################################
### Comparing individual Q-stats and F-stats ####
#################################################


## This allows us to look at the individual Q and F stats for each SNP in the MR 
## We want to be able to see if the individual Qstat (amount of heterogeneity that SNP contributes to the whole heterogenetiy)
## is related to the Fstat (fstat shows weak instrument bias, generally if above 10 then an ok instrument, strong enough to use)

runsim_singlesnp_mr <- function (nsample, nsnp, snp_exp_var, plei_var, conf_exp_var, exp_out_var, conf_out_var, 
                                 pleiotropy_model=c("none", "independent", "dependent"))
{
  g<- make_geno(nsample, nsnp, 0.5) # Makes dataset for number of people and number of SNPs specified
  
  g1 <- subset(g[1:(nsample/2),]) # split dataset in half (so can do two sample MR if want to)
  
  g2 <- subset(g[((nsample/2)+1):nsample,])
  
  conf1 <- rnorm((nsample/2)) # create some random confounding
  conf2 <- rnorm((nsample/2))
  
  pl_effs <- effect_model(nsnp, snp_exp_var, plei_var, model=pleiotropy_model) # get the pleiotropic effects from the effect function
  
  x1 <- make_phen(c(pl_effs$snp_exp_effect, conf_exp_var), cbind(g1, conf1)) # Make data for exposure1 using the first half of the dataset(g1) and the conf1
  x2 <- make_phen(c(pl_effs$snp_exp_effect, conf_exp_var), cbind(g2, conf2)) # for two sample MR
  
  y1 <- make_phen(c(pl_effs$plei_effect, exp_out_var, conf_out_var), cbind(g1, x1, conf1)) # make outcome data 
  dat1 <- get_effs(x1,y1,g1) 
  s1 <- subset(dat1, pval.exposure < 5e-8)
  
  res1 <- mr_singlesnp(s1) ## one sample MR - this is doing single snp MR on s1, s1 is the significant values from dat1, dat1 is x on y in g so not measuring snps on x....
  res1$var <- (res1$se)^2
  res1$wj <- 1/(res1$var)
  Bivw_res1 <- sum(res1$wj*res1$b)/sum(res1$wj)
  res1$Qj <- res1$wj*((res1$b - Bivw_res1)^2)
  fstat <- cbind(SNP=s1$SNP,fstat=((s1$beta.exposure/s1$se.exposure)^2))
  res1 <- merge(res1, fstat, id="SNP", all=T)
  
  y2 <- make_phen(c(pl_effs$plei_effect, exp_out_var, conf_out_var), cbind(g2, x2, conf2))
  dat2 <- get_effs(x2,y2,g2)
  s2 <- subset(dat2, SNP %in% s1$SNP)
  res2 <- mr_singlesnp(s2)
  res2$var <- (res2$se)^2
  res2$wj <- 1/(res2$var)
  Bivw_res2 <- sum(res2$wj*res2$b)/sum(res2$wj)
  res2$Qj <- res2$wj*((res2$b - Bivw_res2)^2)
  fstat <- cbind(SNP=s2$SNP,fstat=((s2$beta.exposure/s2$se.exposure)^2))
  res2 <- merge(res2, fstat, id="SNP", all=T)
  
  mr_results <- rbind(res1, res2)
  
  return(mr_results)
  
}   

a <- runsim_singlesnp_mr(40000, 100, 0.2, 0.05, 0.1, 0.2, 0.2, "none")
b <- runsim_singlesnp_mr(40000, 100, 0.2, 0.05, 0.1, 0.2, 0.2, "independent")
c <- runsim_singlesnp_mr(40000, 100, 0.2, 0.05, 0.1, 0.2, 0.2, "dependent")

## Higher pleiotropy ##
e <- runsim_singlesnp_mr(40000, 100, 0.2, 0.2, 0.1, 0.2, 0.2, "none")
f <- runsim_singlesnp_mr(40000, 100, 0.2, 0.2, 0.1, 0.2, 0.2, "independent")
g <- runsim_singlesnp_mr(40000, 100, 0.2, 0.2, 0.1, 0.2, 0.2, "dependent")


plot(a$Qj, a$fstat)
plot(b$Qj, b$fstat)
plot(c$Qj, c$fstat)

plot(e$Qj, e$fstat)
plot(f$Qj, f$fstat)
plot(g$Qj, g$fstat)


plot(a$Qj, a$fstat, xlim=c(0, 200), ylim=c(0, 400))
plot(b$Qj, b$fstat,xlim=c(0, 200), ylim=c(0, 400))
plot(c$Qj, c$fstat,xlim=c(0, 200), ylim=c(0, 400))

plot(e$Qj, e$fstat,xlim=c(0, 200), ylim=c(0, 400))
plot(f$Qj, f$fstat,xlim=c(0, 200), ylim=c(0, 400))
plot(g$Qj, g$fstat,xlim=c(0, 200), ylim=c(0, 400))



### Just the two sample MR ###

runsim_singlesnp_mr_twosample <- function (nsample, nsnp, snp_exp_var, plei_var, conf_exp_var, exp_out_var, conf_out_var, 
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
  
  res1 <- mr_singlesnp(s1) ## one sample MR 
  res1$var <- (res1$se)^2
  res1$wj <- 1/(res1$var)
  Bivw_res1 <- sum(res1$wj*res1$b)/sum(res1$wj)
  res1$Qj <- res1$wj*((res1$b - Bivw_res1)^2)
  fstat <- cbind(SNP=s1$SNP,fstat=((s1$beta.exposure/s1$se.exposure)^2))
  res1 <- merge(res1, fstat, id="SNP", all=T)
  
  y2 <- make_phen(c(pl_effs$plei_effect, exp_out_var, conf_out_var), cbind(g2, x2, conf2))
  dat2 <- get_effs(x2,y2,g2)
  s2 <- subset(dat2, SNP %in% s1$SNP)
  res2 <- mr_singlesnp(s2)
  res2$var <- (res2$se)^2
  res2$wj <- 1/(res2$var)
  Bivw_res2 <- sum(res2$wj*res2$b)/sum(res2$wj)
  res2$Qj <- res2$wj*((res2$b - Bivw_res2)^2)
  fstat <- cbind(SNP=s2$SNP,fstat=((s2$beta.exposure/s2$se.exposure)^2))
  res2 <- merge(res2, fstat, id="SNP", all=T)
  
  #mr_results <- rbind(res1, res2)
  
  return(res2)
  
}   


a_2sampleMR <- runsim_singlesnp_mr_twosample(40000, 100, 0.2, 0.05, 0.1, 0.2, 0.2, "none")
b_2sampleMR <- runsim_singlesnp_mr_twosample(40000, 100, 0.2, 0.05, 0.1, 0.2, 0.2, "independent")
c_2sampleMR <- runsim_singlesnp_mr_twosample(40000, 100, 0.2, 0.05, 0.1, 0.2, 0.2, "dependent")

## Higher pleiotropy ##
e_2sampleMR <- runsim_singlesnp_mr_twosample(40000, 100, 0.2, 0.2, 0.1, 0.2, 0.2, "none")
f_2sampleMR <- runsim_singlesnp_mr_twosample(40000, 100, 0.2, 0.2, 0.1, 0.2, 0.2, "independent")
g_2sampleMR <- runsim_singlesnp_mr_twosample(40000, 100, 0.2, 0.2, 0.1, 0.2, 0.2, "dependent")



plot(a_2sampleMR$Qj, a_2sampleMR$fstat, xlim=c(0, 200), ylim=c(0, 400))
plot(b_2sampleMR$Qj, b_2sampleMR$fstat,xlim=c(0, 200), ylim=c(0, 400))
plot(c_2sampleMR$Qj, c_2sampleMR$fstat,xlim=c(0, 200), ylim=c(0, 400))

plot(e_2sampleMR$Qj, e_2sampleMR$fstat,xlim=c(0, 200), ylim=c(0, 400))
plot(f_2sampleMR$Qj, f_2sampleMR$fstat,xlim=c(0, 200), ylim=c(0, 400))
plot(g_2sampleMR$Qj, g_2sampleMR$fstat,xlim=c(0, 200), ylim=c(0, 400))


plot(a_2sampleMR$Qj, a_2sampleMR$fstat)
plot(b_2sampleMR$Qj, b_2sampleMR$fstat)
plot(c_2sampleMR$Qj, c_2sampleMR$fstat)

plot(e_2sampleMR$Qj, e_2sampleMR$fstat)
plot(f_2sampleMR$Qj, f_2sampleMR$fstat)
plot(g_2sampleMR$Qj, g_2sampleMR$fstat)




