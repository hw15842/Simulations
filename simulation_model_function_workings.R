

library(simulateGP)
library(TwoSampleMR)
nsample <-
nsnp <- 
MAF <-
snp_exp_var <-
plei_var <- 
conf <- rnorm()
conf_exp_var <- 
exp_out_var <-
conf_out_var <-

   

  effect_model <- function(nsnp, snp_exp_var, plei_var, noise_var, model=c("none", "independent", "dependent")[1])
  { 
    snp_exp_effect <- abs(choose_effects(nsnp, snp_exp_var)) 
    
    if(model == "dependent")
    {
      plei_effect <- abs(choose_effects(nsnp, plei_var))
      snp_exp_effect <- snp_exp_effect[order(snp_exp_effect, decreasing=FALSE)]
      plei_effect <- plei_effect[order(plei_effect, decreasing=TRUE)]
      noise_effect <- choose_effects(nsnp, noise_var)
      plei_effect <- rowSums(cbind(plei_effect,noise_effect))
      #plei_effect <- plei_effect * plei_var / snp_exp_var ## Could add noise (multiple by number between0 and 1), make smaller then add noise (x0.8 +0.8)
    } else if(model == "independent") {
      plei_effect <- abs(choose_effects(nsnp, plei_var))
    } else if (model == "none") {
      plei_effect <- choose_effects(nsnp, 0)
    } else {
      stop("Must specify none, independent or dependent")
    }
    return(list(plei_effect = plei_effect, snp_exp_effect = snp_exp_effect))
  } 
   

## so remake a variable with same var as snp_exp, make it around zero, make it smaller, add it on

z <- effect_model(500, 0.2, 0.1, model="none") 
y <- effect_model(500, 0.2, 0.1, model="independent")
x <- effect_model(500, 0.2, 0.1, 0.0005, model="dependent")

noise_effect <- choose_effects(500, 0.0005)  
snp_exp_effect <- abs(choose_effects(500, 0.2)) 
snp_exp_effect <- snp_exp_effect[order(snp_exp_effect, decreasing=FALSE)]
snp_exp_effect <- rowSums(cbind(snp_exp_effect,noise_effect))


  
  ## MAF not that important will always be 0.5, need conf1 and conf2
  ## get out res$IVW and heterogeneity
  
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
    het2 <- mr_heterogeneity(s2, method_list="mr_ivw")
    
    res1$type <- "one sample"
    res2$type <- "two sample"
    
    het1$type <- "one sample"
    het2$type <- "two sample"
    
    mr_results <- rbind(res1=res1, res2=res2)
    
    heterogeneity <- rbind(het1=het1, het2=het2)
    
    #return(list(rbind(res1=res1, res2=res2), (rbind(het1=het1, het2=het2)), num_sig_snps=num_sig_snps))
    
    #return(rbind(res1, res2))
    
    return(list(mr_results=mr_results, heterogeneity=heterogeneity))
    
}   

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
  res2 <- mr(s2, method_list = "mr_ivw")
  
  het1 <- mr_heterogeneity(s1, method_list="mr_ivw")
  het2 <- mr_heterogeneity(s2, method_list="mr_ivw")
  
  res1$type <- "one sample"
  res2$type <- "two sample"
  
  het1$type <- "one sample"
  het2$type <- "two sample"
  
  mr_results <- rbind(res1=res1, res2=res2)
  
  heterogeneity <- rbind(het1=het1, het2=het2)
  
  #return(list(rbind(res1=res1, res2=res2), (rbind(het1=het1, het2=het2)), num_sig_snps=num_sig_snps))
  
  #return(rbind(res1, res2))
  
  return(list( heterogeneity=heterogeneity))
  
}   


mr_method_list()

a <- runsim(40000, 100, 0.2, 0.05, 0.1, 0.2, 0.2, "none")
b <- runsim(40000, 100, 0.2, 0.05, 0.1, 0.2, 0.2, "independent")
c <- runsim2(40000, 100, 0.2, 0.05, 0.1, 0.2, 0.2, "dependent")

c1 <- runsim(40000, 100, 0.2, 0.05, 0.1, 0.2, 0.2, "dependent")
c2 <- runsim(40000, 100, 0.2, 0.05, 0.1, 0.2, 0.2, "dependent")
c3 <- runsim(40000, 100, 0.2, 0.05, 0.1, 0.2, 0.2, "dependent")
c4 <- runsim(40000, 100, 0.2, 0.05, 0.1, 0.2, 0.2, "dependent")

c1$heterogeneity
c2$heterogeneity
c3$heterogeneity
c4$heterogeneity


f<- rbind(c1$heterogeneity$Q_pval, c2$heterogeneity$Q_pval, c3$heterogeneity$Q_pval, c4$heterogeneity$Q_pval)
colnames(f)<-c("one_sample", "two_sample")

table(f[,2] < 0.05)
table(f <0.05)


new_function <- function(num_runs)
{
  for (i in 1:num_runs)
  {
    sim_fun[i] <- runsim2(40000, 100, 0.2, 0.05, 0.1, 0.2, 0.2, "dependent")
  }
  f <- rbind(sim_fun[[i]]$Q_pval)
}

new_function(5)

c <- replicate (10, runsim(40000, 100, 0.2, 0.05, 0.1, 0.2, 0.2, "dependent"))


  
c <- runsim(40000, 100, 0.2, 0.05, 0.1, 0.2, 0.2, "dependent")
het <- c$heterogeneity


for(i in 1:10)
{
  sim_fun[i] <- runsim(40000, 100, 0.2, 0.05, 0.1, 0.2, 0.2, "dependent")
}
## Doesnt give the heterogeneity... 

for(i in 1:10)
{
  sim_fun[i] <- runsim2(40000, 100, 0.2, 0.05, 0.1, 0.2, 0.2, "dependent")
}
## Does give heterogeneity with the runsim2 function




for (i in 10)
{
  c[i] <- runsim(40000, 100, 0.2, 0.05, 0.1, 0.2, 0.2, "dependent")
  het[i] <- c$heterogeneity
}




n.times <- 10
result <- numeric(n.times) ## assuming your function returns numeric
for (i in 1:n.times) {
  result[i] <- runsim2(40000, 100, 0.2, 0.05, 0.1, 0.2, 0.2, "dependent")
}





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
  out[[i]] <- runsim(params$nsample[i], params$nsnp[i], params$snp_exp_var[i], params$plei_var[i], params$conf_exp_var[i], 
                     params$exp_out_var[i], params$conf_out_var[i], params$pleiotropy_model[i])
}


out[[1]]$heterogeneity$Q_pval # Works to get the pvals

out[1:4] ## Works to get everything in runs 1 to 4 but cant figure out how to just get heterogeneity for example

### Needs to have a decent number of SNPs that are significant to work

## At the moment showing that there is no heterogeneity for b or c.... why?? Not enough pepole or SNPs?
## Is it to do with the two sample bit, we are creating res1 and res2 from dat1 and dat2 not s1 and s2 which have been subset
## Not using s1 and s2 at the moment.... why?? changed res1 to be calculated using s1 - so only the significant snps



#### Check run_sim is doing what we want #####

## below works before we start subsetting 

# dependent

o <- effect_model(1000, snp_exp_var=0.5, plei_var=0.05, model="dependent")

g <- make_geno(10000, 1000, 0.5)

x <- make_phen(o[[2]], g)

y <- make_phen(c(o[[1]], 0.3), cbind(g,x))



e <- get_effs(x,y,g)

mr_heterogeneity(e)



m <- mr(e)

mr_scatter_plot(m, e)
mr_forest_plot(mr_singlesnp(m, all_method="mr_ivw"))



# independent

o <- effect_model(1000, snp_exp_var=0.5, plei_var=0.5, model="independent")

g <- make_geno(10000, 1000, 0.5)

x <- make_phen(o[[2]], g)

y <- make_phen(c(o[[1]], 0.3), cbind(g,x))



e <- get_effs(x,y,g)

mr_heterogeneity(e)



m <- mr(e)

mr_scatter_plot(m, e)



# no plei

o <- effect_model(1000, snp_exp_var=0.5, plei_var=0.5, model="none")

g <- make_geno(10000, 1000, 0.5)

x <- make_phen(o[[2]], g)

y <- make_phen(c(o[[1]], 0.3), cbind(g,x))



e <- get_effs(x,y,g)

mr_heterogeneity(e)



m <- mr(e)

mr_scatter_plot(m, e)





## Dependent model is giving no heterogeneity, independent is giving heterogenetiy


snp_exp_effect <- abs(choose_effects(1000, 0.5)) 
plei_effect <- snp_exp_effect[order(snp_exp_effect, decreasing=TRUE)]
snp_exp_effect <- snp_exp_effect[order(snp_exp_effect, decreasing=FALSE)]
plei_effect <- plei_effect * 0.05 / 0.5

plot(snp_exp_effect,plei_effect)


## Tried to add noise gaussian
install.packages("RMThreshold")
library(RMThreshold)
mat <- cbind(plei_effect,snp_exp_effect)
mat <- add.Gaussian.noise(mat, mean = 0, stddev = 1, symm = F)

plot(mat)












## Make a grid with the different inputs

params <- expand.grid(
  snp_exp_var = c(0, 0.01, 0.2),
  bxy = c(0, 0.01, 0.1),
  pleiotropy_model = c("none", "independent", "dependent"),
  one_or_two_sample = c(1,2),
  stringsAsFactors = FALSE
)


## Output 
out <- list()
for(i in 1:nrow(params))
{
  out[[i]] <- runsim(params$snp_exp_var[i], params$bxy[i], params$pleiotropy_model[i], 
                     params$one_or_two_sample[i])
}


#do stuff to `out`



for(i in 1:nrow(params))
{
  Sys.sleep(2)
  print(params[i,])
}


a <- runsim(10000, 100, 0.2, 0.05, 0.1, 0.2, 0.2, "none")
