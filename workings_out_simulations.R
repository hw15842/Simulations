effect_model <- function(nsnp, snp_exp_var, plei_var, model=c("independent", "dependent")[1])
{
  snp_exp_effect <- choose_effects(nsnp, snp_exp_var)
  
  if(model == "dependent")
  {
    plei_effect <- snp_exp_effect[order(snp_exp_effect, decreasing=TRUE)]
    plei_effect <- plei_effect * plei_var / snp_exp_var
  } else if(model == "independent") {
    plei_effect <- choose_effects(nsnp, plei_var)
  } else {
    stop("Must specify independent or dependent")
  }
  return(list(plei_effect = plei_effect, snp_exp_effect = snp_exp_effect))
}


g<- make_geno(1000, 10, 0.5)

g1 <- subset(g[1:(1000/2),])

g2 <- subset(g[((1000/2)+1):1000,])

effs <- choose_effects(10, 0.2)

conf <- rnorm(500)

x1 <- make_phen(c(effs, 0.2), cbind(g1, conf))

x2 <- make_phen(c(effs, 0.2), cbind(g2, conf))

## No pleiotropy

y1 <- make_phen(c(0.2, 0.2), cbind(x1, conf))

dat <- get_effs(x1,y1,g1)

s1 <- subset(dat, pval.exposure < 5e-8)

y2 <- make_phen(c(0.2, 0.2), cbind(x2, conf))

dat2 <- get_effs(x2,y2,g2)

s2 <- subset(dat2, SNP %in% s1$SNP)

mr(s2)


## Independent pleiotropy

pl_effs_independent <- effect_model(10, 0.2, 0.05, model="independent")






pl_effs <- choose_effects(10, 0.05)
pl_effs_independent <- effect_model(10, 0.2, 0.05, model="independent")
pl_effs_dependent <- effect_model(10, 0.2, 0.05, model="dependent")




## I want to match the plei_effects to the snp_exp_effects 
## This is the bit i cant get to work...
y1 = make_phen(c(pl_effs_dependent$plei_effect, 0.2, 0.2), cbind(pl_effs_dependent$snp_exp_effect, x1, conf)) ## This gives me the error message

pl_effs_independent
pl_effs_dependent


## Workings
y3 <- make_phen(c(pl_effs, 0.4, 0.4), cbind(g1[,1:5], x1, conf)) 

c(pl_effs_independent$plei_effect, 0.2, 0.2)
cbind(pl_effs_independent$snp_exp_effect, x1, conf)

## Deconstructing the effect_model funtction
snp_exp_effect <- choose_effects(10, 0.2)

plei_effect <- snp_exp_effect[order(snp_exp_effect, decreasing=TRUE)]
plei_effect <- plei_effect * 0.05 / 0.2

plei_effect_independent <- choose_effects(10, 0.05)



### the make_pheno function takes the effect sizes then the matrix of variables then creates new varible y
### x <- make_phen(effs, g) doesnt have the same length for effs and g but still works
### x1 <- make_phen(c(effs, 0.2), cbind(g1, conf)) makes the effect sizes from combinine effs and adding 0.2 on the end, makes variables from bindin g1 and confounder together
## But how does it know which effect goes with which bit?? 
