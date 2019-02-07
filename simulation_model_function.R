#################################
### simulation_model_function ###
#################################

library(simulateGP)
library(TwoSampleMR)



effect_model <- function(nsnp, snp_exp_var, plei_var, model=c("none", "independent", "dependent")[1])
{ 
  snp_exp_effect <- abs(choose_effects(nsnp, snp_exp_var)) 
  
  if(model == "dependent")
  {
    plei_effect <- snp_exp_effect[order(snp_exp_effect, decreasing=TRUE)]
    snp_exp_effect <- snp_exp_effect[order(snp_exp_effect, decreasing=FALSE)]
    plei_effect <- plei_effect * plei_var / snp_exp_var
  } else if(model == "independent") {
    plei_effect <- abs(choose_effects(nsnp, plei_var))
  } else if (model == "none") {
    plei_effect <- choose_effects(nsnp, 0)
  } else {
    stop("Must specify none, independent or dependent")
  }
  return(list(plei_effect = plei_effect, snp_exp_effect = snp_exp_effect))
} 








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
  num_sig_snps <- nrow(s2)
  res2 <- mr(s2, method_list = "mr_ivw")
  
  het1 <- mr_heterogeneity(s1, method_list="mr_ivw")
  het2 <- mr_heterogeneity(s2, method_list="mr_ivw")
  
  res1$type <- "one sample"
  res2$type <- "two sample"
  
  het1$type <- "one sample"
  het2$type <- "two sample"
  
  mr_results <- rbind(res1=res1, res2=res2)
  heterogeneity <- rbind(het1=het1, het2=het2)
  
  return(list(mr_results=mr_results, heterogeneity=heterogeneity, num_sig_snps=num_sig_snps))
  
  
}   




a <- runsim(400000, 1000, 0.2, 0.05, 0.1, 0.2, 0.2, "none")
b <- runsim(400000, 1000, 0.2, 0.05, 0.1, 0.2, 0.2, "independent")
c <- runsim(400000, 1000, 0.2, 0.05, 0.1, 0.2, 0.2, "dependent")




