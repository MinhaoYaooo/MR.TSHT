#' @title Normalize function
#'
#' @description Standardization of a vector to be 0 mean and unit variance
#' @param x a numeric vector of input data
#' @return a numeric vector of standardized output data
#' @examples
#' x <- rnorm(10, 2, 5);
#' y <- normalize(x);
#' @export
#'
normalize <- function(x){
  (x-mean(x)) / sd(x)
}


#' @title Variance of TSHT estimator
#'
#' @description Calculate the variance of TSHT estimator after getting V_Hat
#' @param gamma a numeric vector of marginal coefficients from IVs to treatment
#' @param Gamma a numeric vector of marginal coefficients from IVs to outcome
#' @param V_gamma a numeric covariance matrix of gamma
#' @param V_Gamma a numeric covariance matrix of Gamma
#' @return a numeric scalar value indicating the variance of the estimator
#' @export

var_TSHT <- function(gamma, Gamma, V_gamma, V_Gamma){

  Var_Gg <- t(gamma) %*% V_Gamma %*% gamma + t(Gamma) %*% V_gamma %*% Gamma
  Var_gg <- 4 * t(gamma) %*% V_gamma %*% gamma
  Cov_Gg <- 2 * t(Gamma) %*% V_gamma %*% gamma

  E_Gg <- sum(gamma*Gamma)
  E_gg <- sum(gamma^2) + sum(diag(V_gamma))

  sig2 <- (E_Gg/E_gg)^2 * (Var_Gg/E_Gg^2 + Var_gg / E_gg^2 - 2*Cov_Gg/(E_Gg*E_gg) )

  return(sig2)
}


#' @title TSHT using summary statistics
#'
#' @description Conduct two-stage hard thresholding using summary statistics
#' @param ITT_Y a numeric vector of GWAS summary statistics of the outcome
#' @param ITT_D a numeric vector of GWAS summary statistics of the treatment
#' @param SE_Y a numeric vector of standard errors of ITT_Y
#' @param SE_D a numeric vector of standard errors of ITT_D
#' @param n1 the sample size of GWAS summary statistics  of the outcome
#' @param n2 the sample size of GWAS summary statistics  of the treatment
#' @param tuning a numeric scalar value tuning parameter for TSHT, with default 1
#' @param max_clique an option to replace the majority and plurality voting procedures with finding maximal clique in the IV voting matrix, with default FALSE
#' @param alpha a numeric scalar value between 0 and 1 indicating the significance level for the confidence interval, with default 0.05
#' @return
#'     \item{\code{VHat}}{a numeric vector denoting the set of valid and relevant IVs}
#'     \item{\code{betaHat}}{a numeric scalar denoting the estimate of treatment effect.}
#'     \item{\code{seHat}}{a numeric scalar denoting the estimated standard error of betaHat.}
#'     \item{\code{ci}}{a two dimensional numeric vector denoting the 1-alpha confidence intervals for betaHat with lower and upper endpoints.}
#' @export

TSHT.sumstats <- function(ITT_Y, ITT_D, SE_Y, SE_D, n1, n2, tuning = 1, max_clique=FALSE, alpha = 0.05){

  pz = length(ITT_Y)

  V_gamma <- matrix(nrow = pz, ncol = pz)
  V_Gamma <- matrix(nrow = pz, ncol = pz)

  for (i in 1:pz) {
    for (j in i:pz) {
      if(i==j){
        V_gamma[i,i] <- SE_D[i]^2
        V_Gamma[i,i] <- SE_Y[i]^2
      }else{
        V_gamma[i,j] <- V_gamma[j,i] <- SE_D[i]*SE_D[j] / n1
        V_Gamma[i,j] <- V_Gamma[j,i] <- SE_Y[i]*SE_Y[j] / n2
      }
    }
  } # Construct V_gamma_hat and V_Gamma_hat based on the formula

  SHat <- (1:pz)[(abs(ITT_D) >= (sqrt(log(n1)*tuning*diag(V_gamma))))]

  SHat.bool = rep(FALSE,pz); SHat.bool[SHat] = TRUE

  nCand = length(SHat)
  VHats.bool = matrix(FALSE,nCand,nCand); colnames(VHats.bool) = rownames(VHats.bool) = SHat

  for(j in SHat) {
    beta.j <- ITT_Y[j] / ITT_D[j]
    pi.j = (ITT_Y * ITT_D[j]  - ITT_D * ITT_Y[j]) / ITT_D[j]
    sigmasq.j = diag(V_Gamma)+(ITT_D/ITT_D[j])^2*V_Gamma[j,j]-2*ITT_D/ITT_D[j]*V_Gamma[,j] +
      beta.j^2*(diag(V_gamma)+(ITT_D/ITT_D[j])^2*V_gamma[j,j]-2*ITT_D/ITT_D[j]*V_gamma[,j])
    PHat.bool.j = abs(pi.j) <= sqrt(sigmasq.j*tuning*log(n2))
    VHat.bool.j = PHat.bool.j * SHat.bool
    VHats.bool[as.character(SHat),as.character(j)] = VHat.bool.j[SHat]
  }

  VHats.boot.sym<-VHats.bool
  for(i in 1:dim(VHats.boot.sym)[1]){
    for(j in 1:dim(VHats.boot.sym)[2]){
      VHats.boot.sym[i,j]<-min(VHats.bool[i,j],VHats.bool[j,i])
    }
  }

  if (max_clique) {
    voting.graph <- as.undirected(graph_from_adjacency_matrix(VHats.boot.sym))
    max_block <- largest_cliques(voting.graph)
    VHat <- as_ids(sort(max_block[[1]]))
    VHat <- as.numeric(VHat) # randomly pick the first one if multiple maximal cliques exist
  } else {
    VM= apply(VHats.boot.sym,1,sum)
    VM.m = rownames(VHats.boot.sym)[VM > (0.5 * length(SHat))] # Majority winners
    VM.p = rownames(VHats.boot.sym)[max(VM) == VM] #Plurality winners
    V.set<-NULL
    for(index in union(VM.m,VM.p)){
      V.set<-union(V.set,names(which(VHats.boot.sym[index,]==1)))
    }
    VHat<-NULL
    for(index in V.set){
      VHat<-union(VHat,names(which(VHats.boot.sym[,index]==1)))
    }
    VHat=as.numeric(VHat)
  }

  # VM= apply(VHats.boot.sym,1,sum)
  # VM.m = rownames(VHats.boot.sym)[VM > (0.5 * length(SHat))] # Majority winners
  # VM.p = rownames(VHats.boot.sym)[max(VM) == VM] #Plurality winners
  #
  # VHat <- as.numeric(union(VM.m, VM.p))

  beta_L <- sum(ITT_D[VHat] * ITT_Y[VHat]) / sum(ITT_D[VHat]^2)

  # se_L <- sqrt(mean(diag(V_Gamma)[VHat] + beta_L^2*diag(V_gamma)[VHat]) / sum(ITT_D[VHat]^2))

  se_L <- sqrt(var_TSHT(ITT_D[VHat], ITT_Y[VHat], V_gamma[VHat,VHat], V_Gamma[VHat,VHat]))

  ci_L = c(beta_L-qnorm(1-alpha/2)*se_L, beta_L+qnorm(1-alpha/2)*se_L)

  return(list(VHat=VHat, beta=beta_L, se=se_L, ci=ci_L))
}


#' @title Handle the possible union of CIs
#'
#' @description Take union of CIs accoding to the specific grid sizes
#' @param CI.matrix a matrix, where each row represents a CI
#' @param true.val a numeric scalar indecating the true value of treatment effect
#' @param grid.size a numeric scalar indicating the size of the grid
#' @return
#'     \item{\code{CI.coverage}}{a numeric scalar denoting whether the CI contains the true value}
#'     \item{\code{CI.length}}{a numeric scalar denoting the length of the CI.}
#'     \item{\code{grid.seq}}{a numeric vector indicating the sequence of grids in CI.matrix}
#' @export

analysis.CI<-function(CI.matrix,true.val,grid.size){
  #### number of rows in the CI.matrix
  d<-dim(CI.matrix)[1]
  CI.coverage<-0
  CI.len<-0
  grid.seq<-NULL
  for (l in 1: d){
    CI.len<-CI.len+CI.matrix[l,2]-CI.matrix[l,1]
    if((CI.matrix[l,2]>true.val)*(CI.matrix[l,1]<true.val)==1){
      CI.coverage<-1
    }
    grid.seq<-c(grid.seq,seq(CI.matrix[l,1],CI.matrix[l,2],by=grid.size))
  }
  return(list(CI.coverage=CI.coverage,CI.len=CI.len,grid.seq=grid.seq))
}


#' @title Searhing method
#'
#' @description Get a possible confidence interval using searching method
#' @param ITT_Y a numeric vector of GWAS summary statistics of the outcome
#' @param ITT_D a numeric vector of GWAS summary statistics of the treatment
#' @param SE_Y a numeric vector of standard errors of ITT_Y
#' @param SE_D a numeric vector of standard errors of ITT_D
#' @param VHat a numeric vector denoting the set of valid and relevant IVs, usually obtained in the results of TSHT.sum
#' @param n1 the sample size of GWAS summary statistics  of the outcome
#' @param n2 the sample size of GWAS summary statistics  of the treatment
#' @param beta.grid a numeric vector of the grid for sampling
#' @return
#'     \item{\code{CI.search}}{a two dimensional numeric vector denoting the 1-alpha confidence intervals for betaHat with lower and upper endpoints constructed by Searching method.}
#'     \item{\code{rule}}{a boolean scalar denoting whether the identification condition is satisfied or not}
#'     \item{\code{valid.grid}}{a numeric vector denoting the candidates satisfying certain thresholding condition}
#' @export

Searching.sumstats <- function(ITT_Y, ITT_D, V_Gamma, V_gamma, VHat, n1, n2, beta.grid = NULL){
  threshold.size<-(length(VHat)/2)
  if(is.null(beta.grid)){
    beta.grid<-seq(-5,5,by=max(n1,500)^{-1})
  }
  n.beta<-length(beta.grid)
  valid.grid<-rep(NA,n.beta)

  Tn<-sqrt(2.005*log(n.beta))
  for(j in 1:n.beta){
    b<-beta.grid[j]
    se.b<-sqrt(diag(V_Gamma[VHat,VHat])+b^2*diag(V_gamma[VHat,VHat]))
    valid.grid[j]<-sum(abs(ITT_Y[VHat]-b*ITT_D[VHat])<Tn*se.b)
  }
  if(length(beta.grid[which(valid.grid>threshold.size)])==0){
    rule=FALSE
    warning("Rule Fails. SS will give misleading CIs, SEs, and p-values.")
    sel.index<-which(valid.grid==max(valid.grid))
  }else{
    rule=TRUE
    sel.index<-which(valid.grid>threshold.size)
  }
  CI<-matrix(NA,nrow=length(sel.index),ncol=2)
  CI[,1]<-beta.grid[sel.index]
  upper.index<-sel.index+1
  upper.index[length(upper.index)]<-min(upper.index[length(upper.index)],n.beta)
  CI[,2]<-beta.grid[upper.index]
  if(dim(as.matrix(CI))[1]==1){
    CI.search<-as.matrix(CI)
  }else{
    uni<- Intervals(CI)
    ###### construct the confidence interval by taking a union
    CI.search<-as.matrix(interval_union(uni))
  }
  return(list(CI.search=CI.search,rule=rule,valid.grid=valid.grid))
}


#' @title Searhing-and-Sampling method
#'
#' @description Get a possible confidence interval using searching method
#' @param ITT_Y a numeric vector of GWAS summary statistics of the outcome
#' @param ITT_D a numeric vector of GWAS summary statistics of the treatment
#' @param SE_Y a numeric vector of standard errors of ITT_Y
#' @param SE_D a numeric vector of standard errors of ITT_D
#' @param VHat a numeric vector denoting the set of valid and relevant IVs, usually obtained in the results of TSHT.sum
#' @param n1 the sample size of GWAS summary statistics  of the outcome
#' @param n2 the sample size of GWAS summary statistics  of the treatment
#' @param beta.grid a numeric vector of the grid for sampling
#' @param rho a numeric scalar denoting thresholding level used in the sampling property
#' @param S a positive integer indicating the number of bootstrap resampling for computing the confidence interval, with default 1000
#' @return
#'     \item{\code{CI.union}}{a two dimensional numeric vector denoting the 1-alpha confidence intervals for betaHat with lower and upper endpoints constructed by Searching-Sampling method.}
#'     \item{\code{rule}}{a boolean scalar denoting whether the identification condition is satisfied or not.}
#'     \item{\code{CI}}{a numeric matrix denoting the confidence intervals for betaHat constructed by valid candidates for betaHat.}
#' @export

Searching.Sampling.sumstats <- function(ITT_Y, ITT_D, V_Gamma, V_gamma, VHat, n1, n2, beta.grid = NULL, rho=NULL, S=1000){
  pz <- length(ITT_D)
  threshold.size<-(length(VHat)/2)
  if(is.null(rho)){
    rho<-(log(n1)/S)^{1/(2*length(VHat))}/6
  }
  if(is.null(beta.grid)){
    beta.grid<-seq(-5,5,by=max(n1,500)^{-1})
  }
  n.beta<-length(beta.grid)
  valid.grid.sample<-matrix(NA,nrow=S,ncol=n.beta)
  #SE.norm<-sqrt(SigmaSqD / n1)
  Tn<-sqrt(2.005*log(n.beta))

  for (s in 1:S) {
    ITT_Y.sample<-ITT_Y-mvrnorm(1,rep(0,pz),V_Gamma)
    ITT_D.sample<-ITT_D-mvrnorm(1,rep(0,pz),V_gamma)
    for(j in 1:n.beta){
      b<-beta.grid[j]
      se.b<-sqrt(diag(V_Gamma[VHat,VHat])+b^2*diag(V_gamma[VHat,VHat]))
      valid.grid.sample[s,j]<-sum(abs(ITT_Y.sample[VHat]-b*ITT_D.sample[VHat])<rho*Tn*se.b)
    }
  }

  CI<-matrix(NA,nrow=S,ncol=2)
  for(s in 1:S){
    if(length(which(valid.grid.sample[s,]>threshold.size))>0){
      CI[s,1]<-min(beta.grid[which(valid.grid.sample[s,]>threshold.size)])
      CI[s,2]<-max(beta.grid[which(valid.grid.sample[s,]>threshold.size)])
    }
  }
  CI<-na.omit(CI)

  while((dim(as.matrix(CI))[1]<min(0.05*S,50)) && (rho<0.5)){
    #print(as.matrix(CI)[1])
    #print(rho)
    rho<-1.25*rho
    for(s in 1:S){
      ITT_Y.sample<-ITT_Y-mvrnorm(1,rep(0,pz),V_Gamma)
      ITT_D.sample<-ITT_D-mvrnorm(1,rep(0,pz),V_gamma)
      for(j in 1:n.beta){
        b<-beta.grid[j]
        se.b<-sqrt(diag(V_Gamma[VHat,VHat])+b^2*diag(V_gamma[VHat,VHat]))
        valid.grid.sample[s,j]<-sum(abs(ITT_Y.sample[VHat]-b*ITT_D.sample[VHat])<rho*Tn*se.b)
      }
    }
    CI<-matrix(NA,nrow=S,ncol=2)
    for(s in 1:S){
      if(length(which(valid.grid.sample[s,]>threshold.size))>0){
        CI[s,1]<-min(beta.grid[which(valid.grid.sample[s,]>threshold.size)])
        CI[s,2]<-max(beta.grid[which(valid.grid.sample[s,]>threshold.size)])
      }
    }
    CI<-na.omit(CI)
  }
  rule<-TRUE
  if(dim(as.matrix(CI))[1]==0){
    rule<-FALSE
    CI.union<-t(as.matrix(c(min(beta.grid),max(beta.grid))))
  }else if(dim(as.matrix(CI))[1]==1){
    CI.union<-as.matrix(CI)
  }else{
    uni<- Intervals(CI)
    ###### construct the confidence interval by taking a union
    CI.union<-as.matrix(interval_union(uni))
  }
  return(list(CI.union=CI.union,rule=rule,CI=CI))
}

