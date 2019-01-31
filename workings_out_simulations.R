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

## Is this actually doing what I think it is?
## snp_exp_effect is not in any sort of order
## plei_effect is the snp_exp_effect in reverse order and then multiplied by the plei_var and divided by the snp_exp_var (for dependent)
## I then want to make the large pleiotropy associate with the small snp_exp_effect sizes
## However this is still just ordering them? Which basically means negative numbers are small
## But that isnt what we want? We want zero to be small and numbers further away from zero to be considered large?? 
## dont we want it that the numbers further away from zero wnat the larger weighting?
## Could I make it so that they get a weighting depending on how far from zero they are? 
## So the pleiotropy effect sizes get ordered according to how far from zero they are
## Then the snp_exp_effect gets weighted according to how far from zero they are
## And then can reverse match those up? 


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


install.packages("rowr")
library(rowr)

## Need to make the small effect sizes of the snps on x match with the large pleiotropy sizes 
## I want to match the plei_effects to the snp_exp_effects 


y1 = make_phen(c(pl_effs_dependent$plei_effect, 0.2, 0.2), cbind.fill(t(pl_effs_dependent$snp_exp_effect), x1, conf)) 

## transposed the matrix and used cbind to fill, get the two right sizes for vector_of_effect_sizes and matrix_of_variables
## However when run the make_pheno command get a table full of NaN


# make_phen(vector_of_effect_sizes, matrix_of_variables)
# number of rows in the vector_of_effect_sizes needs to equal the number of columns in the matrix_of_variables

vector_of_effect_sizes<- c(pl_effs_dependent$plei_effect, 0.2, 0.2)
matrix_of_variables <- cbind.fill(t(pl_effs_dependent$snp_exp_effect), x1, conf)
# These two are the right sizes and match, so why am I getting NaN? 

pl_effs_independent
pl_effs_dependent


## Workings
y3 <- make_phen(c(pl_effs, 0.4, 0.4), cbind(g1[,1:10], x1, conf)) 


y_test <- make_phen(vector_of_effect_sizes, matrix_of_variables)

c(0.3, 0.4)
cbind(x1,x2)



## Need to make pl_effs_dependent$snp_exp_effect into columns rather than rows then can add it to the matrix


## non dependent pleiotropy example
pl_effs <- choose_effects(10, 0.05)

y3 <- make_phen(c(pl_effs, 0.4, 0.4), cbind(g1[,1:10], x1, conf)) 

c(pl_effs, 0.4, 0.4)
head(cbind(g1[,1:10], x1, conf))



## Deconstructing the effect_model funtction
snp_exp_effect <- choose_effects(10, 0.2)

plei_effect <- snp_exp_effect[order(snp_exp_effect, decreasing=TRUE)]
plei_effect <- plei_effect * 0.05 / 0.2

plei_effect_independent <- choose_effects(10, 0.05)


## make_pheno function figuring it out

make_phen <- function(effs, indep, vy=1, vx=rep(1, length(effs)), my=0)
{
  if(is.null(dim(indep))) indep <- cbind(indep)
  stopifnot(ncol(indep) == length(effs))
  stopifnot(length(vx) == length(effs))
  cors <- effs * vx / sqrt(vx) / sqrt(vy)
  sc <- sum(cors^2)
  if(sc >= 1)
  {
    print(sc)
    stop("effects explain more than 100% of variance")
  }
  cors <- c(cors, sqrt(1-sum(cors^2)))
  indep <- t(t(scale(cbind(indep, rnorm(nrow(indep))))) * cors * c(vx, 1))
  y <- drop(scale(rowSums(indep)) * sqrt(vy)) + my
  return(y)
}

## deconstructiong make_pheno function
effs <- vector_of_effect_sizes 
vx <- rep(1, length(effs)) 
vy=1
indep <- matrix_of_variables 

cors <- effs * vx / sqrt(vx) / sqrt(vy)

indep <- t(t(scale(cbind(indep, rnorm(nrow(indep))))) * cors * c(vx, 1))

#scale, with default settings, will calculate the mean and standard deviation of the entire vector, 
#then "scale" each element by those values by subtracting the mean and dividing by the sd
#scale is calculating the mean and sd of the cbind dataframe (variables and a random norm dist of the same number of variables)
# then it is scaling each element by those values
## Why does it scale it? 

a <- cbind(indep, rnorm(nrow(indep))) ## Works fine
b <- scale(cbind(indep, rnorm(nrow(indep)))) ## throws up the NaN, something about the scale bit isnt working
## Is it because all of them in that column are the same value? so the mean is the same value so getting zero and standar error is zero
## 

c <- t(scale(cbind(indep, rnorm(nrow(indep)))))
d <- t(scale(cbind(indep, rnorm(nrow(indep))))) * cors * c(vx, 1)





### the make_pheno function takes the effect sizes then the matrix of variables then creates new varible y
### x <- make_phen(effs, g) doesnt have the same length for effs and g but still works
### x1 <- make_phen(c(effs, 0.2), cbind(g1, conf)) makes the effect sizes from combinine effs and adding 0.2 on the end, makes variables from bindin g1 and confounder together
## But how does it know which effect goes with which bit?? 
