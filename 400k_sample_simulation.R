###########################
### 400K sample dataset ###
###########################

### Firstly Simulate a data set with 400K samples, 1K SNPs, split in half for discovery and replication ##
## Have confounding but no pleiotropy in first instance ##

library(simulateGP)

## Simuilate data set, 400k samples 1k SNPs, MAF 0.5
g<- make_geno(400000, 1000, 0.5)

## Split the dataset in half. Is there any reason I cant just take top half and bottom half? It is all randomly simulated anyway? 

g1 <- subset(g[1:200000,])
dim(g1)


g2 <- subset(g[200001:400000,])
dim(g2)

## Choose effects, 
effs <- choose_effects(1000, 0.05)

## Make confounder term

u1 <- rnorm(200000)


## Find SNPs in g1 (discovery sample)

## the 1000 SNPs explain 5% of the variance in X and the confounder (u) expalins 40% of the variance in X

x1 <- make_phen(c(effs, 0.4), cbind(g1, u1))

y1 <- make_phen(c(0.4, 0.4), cbind(x1, u1)) ### X explains 40% of the variance in Y and U explains 40% of the variance in Y


## Find the significant SNPs at <5e-8

dat1 <- get_effs(x1,y1,g1)

table(dat1$pval.exposure < 5e-8)  ## How many SNPs have a pvalue <5e-8 (TRUE is the number of sig SNPs)

# Extract the significant SNPs

s1 <- subset(dat1, pval.exposure < 5e-8)
dim(s1) ### should be the same number of rows as significnat SNPs found above


## Do in sample 2

x2 <- make_phen(c(effs, 0.4), cbind(g2, u1))

y2 <- make_phen(c(0.4, 0.4), cbind(x2, u1))

dat2 <- get_effs(x2,y2,g2)


## Extract the SNPs found in Sample 1 from Sample 2

s2 <- subset(dat2, SNP %in% s1$SNP)
dim(s2) ## Should have the same number of SNPs as s1

## Run MR on these snps

library(TwoSampleMR)
mr(s2)

# make scatter plot of simulated data 
res1 <- mr(s2, method_list="mr_ivw") 
mr_scatter_plot(res1, s2) 

# Forest plot
mr_forest_plot(mr_singlesnp(s2, all_method="mr_ivw"))


## Check for heterogeneity
mr_heterogeneity(s2, method_list="mr_ivw") ## There shouldnt be any heterogeneity as we have not introduced pleiotropy yet


#Sensitivity analysis 
sens_1 <- mr(s2, method_list=c("mr_ivw", "mr_weighted_median", "mr_egger_regression", "mr_weighted_mode"))
sens_1
mr_scatter_plot(sens_1, s2)



########################################
### Adding Pleiotropy to the dataset ###
########################################

pl_effs <- choose_effects(500, 0.05) ## half of the SNPs also have pleiotropc effects 

## Start with the data in g1 (disvorey sample)
# Create X again, NOT adding the pleiotropic effects (same as x1 but for clarity creating a new variable)

x3 <- make_phen(c(effs, 0.4), cbind(g1, u1))


## Add the pleiotropic effects to Y1 (first sample set)

## the pleitropic effects (500 SNPs explain 5% of the variance in Y) affect the first 500 SNPs in g1, they do not effect X as have not added pl_eff term into x3
## 40% of variance in Y3 explined by x1, 40% of variance in Y3 explined by u1

y3 <- make_phen(c(pl_effs, 0.4, 0.4), cbind(g1[,1:500], x3, u1)) 

dat3 <- get_effs(x3,y3,g1)

table(dat3$pval.exposure < 5e-8)

s3 <- subset(dat3, pval.exposure < 5e-8)
dim(s3)


## Do in sample 2

x4 <- make_phen(c(effs, 0.4), cbind(g2, u1))

y4 <- make_phen(c(pl_effs, 0.4, 0.4), cbind(g2[,1:500], x4, u1)) 

dat4 <- get_effs(x4,y4,g2)


## Extract the SNPs found in Sample 1 from Sample 2

s4 <- subset(dat4, SNP %in% s3$SNP)
dim(s4) ## Should have the same number of SNPs as s3

## Run MR on these snps

mr(s4)


# make scatter plot of simulated data 
res2 <- mr(s4, method_list="mr_ivw") 
mr_scatter_plot(res2, s4) 

# Forest plot
mr_forest_plot(mr_singlesnp(s4, all_method="mr_ivw"))


## Check for heterogeneity
mr_heterogeneity(s4, method_list="mr_ivw") ## There should be significant heterogeneity as we have introduced pleiotropy


##############################################################################
### Simulate pleiotropy so the smaller effect sizes have larger pleiotropy ###
##############################################################################
















