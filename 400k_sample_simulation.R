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

s1_sig <- subset(dat1, pval.exposure < 5e-8)
dim(s1_sig) ### should be the same number of rows as significnat SNPs found above


## Do in sample 2

x2 <- make_phen(c(effs, 0.4), cbind(g2, u1))

y2 <- make_phen(c(0.4, 0.4), cbind(x2, u1))

dat2 <- get_effs(x2,y2,g2)


## Extract the SNPs found in Sample 1 from Sample 2

s2_snps <- subset(dat2, SNP %in% s1_sig$SNP)
dim(s2_snps) ## Should have the same number of SNPs as s1_sig

## Run MR on these snps

library(TwoSampleMR)
mr(s2_snps)

# make scatter plot of simulated data 
res <- mr(s2_snps, method_list="mr_ivw") 
mr_scatter_plot(res, s2_snps) 

mr_heterogeneity(s2_snps, method_list="mr_ivw") ## There shouldnt be any heterogeneity as we have not introduced pleiotropy yet







