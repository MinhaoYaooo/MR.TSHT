# MR.TSHT
`MR.TSHT` is an R package to conduct two-sample Mendelian randomization using two-stage hard thresholding. The inputs of `MR.TSHT` are summary statistics which can be obtained in Genome-Wide Association Studies (GWAS). The outputs of `MR.TSHT` are the estimation of causal effect and its corresponding confidence interval, and a set of valid IVs <img src="https://render.githubusercontent.com/render/math?math=\hat{V}">.

# Installation
You can install the development version of `MR.TSHT` from Github via the `devtools` package.
```
devtools::install_github("MinhaoYaooo/MR.TSHT")
```

# Example

We first set the following parameters:

```
n1 <- 10000  # number of sample 1
n2 <- 10000 # number of sample 2
beta <- 1 # magnitude of causal effect
tau <- 0.5 # IV strength
gamma <- tau*c(1,1,1,1,1,-1,-1,-1,-1,-1)  # SNP-exposure effect
pi <- tau*c(0,0,0,0,0,0,1,1,1,1) # invalidity of the candidate IVs
```

Since `MR.TSHT` is used for two-sample Mendelian randomization, we create two independent samples whose genotypes follow the same distribution:

```
pSNP <- runif(10, 0.05, 0.95)  # generate allele frequency

Z1 <- sapply(pSNP, rbinom, n = n1, size = 2) # generate raw genotype of sample 1
Z1 <- apply(Z1, 2, normalize) # normalize it
    
Z2 <- sapply(pSNP, rbinom, n = n2, size = 2) # generate raw genotype of sample 2
Z2 <- apply(Z2, 2, normalize) # normalize it
  
D1 <- Z1 %*% gamma + rnorm(n1, 0, 1) # generate the exposure of sample 1
D2 <- Z2 %*% gamma + rnorm(n2, 0, 1) # generate the exposure of sample 2
Y2 <- D2*beta + Z2 %*% pi + rnorm(n2, 0, 1) # generate the outcome of sample 2
```

After generating two samples and their corresponding phenotypes, we calculate the summary statistics as the inputs:

```
library(bigstatsr)

GWAS1 <- big_univLinReg(X = as_FBM(Z1), y.train = D1) # calculate summary statistics of D~Z
GWAS2 <- big_univLinReg(X = as_FBM(Z2), y.train = Y2) # calculate summary statistics of Y~Z
```

Then with the summary statistics, we can estimate the causal effect with `MR.TSHT`:

```
library(igraph)

ITT_Y = as.numeric(GWAS2$estim); ITT_D = as.numeric(GWAS1$estim) 
SE_Y <- as.numeric(GWAS2$std.err); SE_D <- as.numeric(GWAS1$std.err);
    
TSHT.sum <- TSHT.sumstats(ITT_D, ITT_Y, SE_D, SE_Y, n1, n2, max_clique = TRUE)
```

If we wish to use Searching-and-Sampling method for constructing the robust CI, then we can use the `SearchingSampling.sumstats` function:

```
library(intervals)

SS.sum <- SearchingSampling.sumstats(ITT_D, ITT_Y, SE_D, SE_Y, n1, n2)
```
