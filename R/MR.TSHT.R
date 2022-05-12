#' @title TSHT using summary statistics
#'
#' @description Conduct two-stage hard thresholding using summary statistics
#' @param ITT_D a numeric vector of GWAS summary statistics of the treatment
#' @param ITT_Y a numeric vector of GWAS summary statistics of the outcome
#' @param SE_D a numeric vector of standard errors of ITT_D
#' @param SE_Y a numeric vector of standard errors of ITT_Y
#' @param n1 the sample size of GWAS summary statistics  of the treatment
#' @param n2 the sample size of GWAS summary statistics  of the outcome
#' @param tuning a numeric scalar value tuning parameter for TSHT, with default 1
#' @param max_clique an option to replace the majority and plurality voting procedures with finding maximal clique in the IV voting matrix, with default FALSE
#' @param alpha a numeric scalar value between 0 and 1 indicating the significance level for the confidence interval, with default 0.05
#' @return
#'     \item{\code{VHat}}{a numeric vector denoting the set of valid and relevant IVs}
#'     \item{\code{betaHat}}{a numeric scalar denoting the estimate of treatment effect.}
#'     \item{\code{seHat}}{a numeric scalar denoting the estimated standard error of betaHat.}
#'     \item{\code{ci}}{a two dimensional numeric vector denoting the 1-alpha confidence intervals for betaHat with lower and upper endpoints.}
#' @export
TSHT.sumstats <- function(ITT_D, ITT_Y, SE_D, SE_Y, n1, n2, tuning = 2.01, max_clique=TRUE, alpha=0.05){

  pz = length(ITT_Y)

  V_gamma <- matrix(nrow = pz, ncol = pz)
  V_Gamma <- matrix(nrow = pz, ncol = pz)

  for (i in 1:pz) {
    for (j in i:pz) {
      if(i==j){
        V_gamma[i,i] <- SE_D[i]^2
        V_Gamma[i,i] <- SE_Y[i]^2
      }else{
        V_gamma[i,j] <- V_gamma[j,i] <- ITT_D[i]*ITT_D[j] / n1
        V_Gamma[i,j] <- V_Gamma[j,i] <- ITT_Y[i]*ITT_Y[j] / n2
      }
    }
  } # Construct V_gamma_hat and V_Gamma_hat based on the formula

  SHat <- (1:pz)[(abs(ITT_D) >= (sqrt(log(pz)*tuning*diag(V_gamma))))]

  SHat.bool = rep(FALSE,pz); SHat.bool[SHat] = TRUE

  nCand = length(SHat)
  VHats.bool = matrix(FALSE,nCand,nCand); colnames(VHats.bool) = rownames(VHats.bool) = SHat

  for(j in SHat) {
    beta.j <- ITT_Y[j] / ITT_D[j]
    pi.j = (ITT_Y * ITT_D[j]  - ITT_D * ITT_Y[j]) / ITT_D[j]
    sigmasq.j = diag(V_Gamma)+(ITT_D/ITT_D[j])^2*V_Gamma[j,j]-2*ITT_D/ITT_D[j]*V_Gamma[,j] +
      beta.j^2*(diag(V_gamma)+(ITT_D/ITT_D[j])^2*V_gamma[j,j]-2*ITT_D/ITT_D[j]*V_gamma[,j])
    PHat.bool.j = abs(pi.j) <= sqrt(sigmasq.j*tuning^2*log(pz))
    VHat.bool.j = PHat.bool.j * SHat.bool
    VHats.bool[as.character(SHat),as.character(j)] = VHat.bool.j[SHat]
  }

  VHats.boot.sym<-VHats.bool
  for(i in 1:dim(VHats.boot.sym)[1]){
    for(j in 1:dim(VHats.boot.sym)[2]){
      VHats.boot.sym[i,j]<-min(VHats.bool[i,j],VHats.bool[j,i])
    }
  }

  diag(VHats.boot.sym) <- 1
  # Voting method
  VM= apply(VHats.boot.sym,1,sum)
  VM.m = rownames(VHats.boot.sym)[VM > (0.5 * length(SHat))] # Majority winners
  VM.p = rownames(VHats.boot.sym)[max(VM) == VM]

  if (max_clique) {
    voting.graph <- igraph::as.undirected(igraph::graph_from_adjacency_matrix(VHats.boot.sym))
    max.clique <- igraph::largest.cliques(voting.graph)
    VHat <- unique(igraph::as_ids(Reduce(c,max.clique))) # take the union if multiple max cliques exist
    VHat <- sort(as.numeric(VHat))
  } else{
    V.set<-NULL
    for(index in VM.p){
      V.set<-union(V.set,names(which(VHats.boot.sym[index,]==1)))
    }
    VHat<-NULL
    for(index in V.set){
      VHat<-union(VHat,names(which(VHats.boot.sym[,index]==1)))
    }
    VHat=sort(as.numeric(VHat))
  }

  if (max_clique) {
    max.clique.mat <- matrix(0,nrow = length(max.clique),ncol = length(max.clique[[1]]))
    CI.temp <- matrix(0,nrow = length(max.clique), ncol = 2)
    beta.temp <- matrix(0,nrow = length(max.clique), ncol = 1)
    betavar.temp <- matrix(0,nrow = length(max.clique), ncol = 1)
    for (i in 1:length(max.clique)) {
      temp <- SHat[sort(as.numeric(max.clique[[i]]))]
      max.clique.mat[i,] <- temp
      betaHat = sum(ITT_D[temp] * ITT_Y[temp]) / sum(ITT_D[temp]^2)
      betaVarHat = var_TSHT(ITT_D[temp], ITT_Y[temp], V_gamma[temp,temp], V_Gamma[temp,temp])
      ci = c(betaHat - qnorm(1-alpha/2) * sqrt(betaVarHat),betaHat + qnorm(1-alpha/2) * sqrt(betaVarHat))
      CI.temp[i,] <- ci
      beta.temp[i,] <- betaHat
      betavar.temp[i,] <- betaVarHat
    }
    uni<- intervals::Intervals(CI.temp)
    ###### construct the confidence interval by taking a union
    CI.union<-as.matrix(intervals::interval_union(uni))

  } else{
    betaHat <- sum(ITT_D[VHat] * ITT_Y[VHat]) / sum(ITT_D[VHat]^2)

    betaVarHat <- var_TSHT(ITT_D[VHat], ITT_Y[VHat], V_gamma[VHat,VHat], V_Gamma[VHat,VHat])

    ci = c(betaHat-qnorm(1-alpha/2)*sqrt(betaVarHat), betaHat+qnorm(1-alpha/2)*sqrt(betaVarHat))
  }



  if (max_clique) {
    returnList <- list(betaHat=betaHat,beta.sdHat = sqrt(betaVarHat),ci=CI.union,SHat=SHat,VHat=VHat,max.clique=max.clique.mat,voting.mat=VHats.boot.sym,
                       beta.clique = beta.temp,beta.sd.clique = sqrt(betavar.temp), CI.clique = CI.temp)
  } else {
    returnList <- list(betaHat=betaHat,beta.sdHat = sqrt(betaVarHat),ci=ci,SHat=SHat,VHat=VHat,voting.mat=VHats.boot.sym)
  }
  return(returnList)
}



#' @title Searhing-and-Sampling method
#'
#' @description Get a possible confidence interval using searching method
#' @param ITT_D a numeric vector of GWAS summary statistics of the treatment
#' @param ITT_Y a numeric vector of GWAS summary statistics of the outcome
#' @param SE_D a numeric vector of standard errors of ITT_D
#' @param SE_Y a numeric vector of standard errors of ITT_Y
#' @param n1 the sample size of GWAS summary statistics  of the treatment
#' @param n2 the sample size of GWAS summary statistics  of the outcome
#' @param CI.init initial interval for beta. If \code{NULL}, it will be generated automatically. (default=\code{NULL})
#' @param a grid size for constructing beta grids (default=0.6)
#' @param Sampling if \code{TRUE}, use the proposed sampling method; else use the proposed searching method. (default=\code{TRUE})
#' @param rho a numeric scalar denoting thresholding level used in the sampling property
#' @param M sampling times. (default=1000)
#' @param prop proportion of intervals kept when sampling. (default=0.1)
#' @return
#'     \item{\code{CI}}{a numeric matrix denoting the confidence intervals for betaHat constructed by valid candidates for betaHat.}
#'     \item{\code{rule}}{a boolean scalar denoting whether the identification condition is satisfied or not.}
#'     \item{\code{VHat}}{valid instruments}
#'     \item{\code{CI.init}}{initial confidence interval for searching and sampling}
#'     \item{\code{TSHT.out}}{standard MR-TSHT output}
#' @export
SearchingSampling.sumstats <- function(ITT_D, ITT_Y, SE_D, SE_Y, n1, n2,
                                       CI.init = NULL,
                                       a=0.6,
                                       Sampling=TRUE,
                                       rho=NULL, M=1000, prop=0.1){

  pz = length(ITT_Y);

  ## Computing (co)variance matrices
  V_gamma <- matrix(nrow = pz, ncol = pz)
  V_Gamma <- matrix(nrow = pz, ncol = pz)

  for (i in 1:pz) {
    for (j in i:pz) {
      if(i==j){
        V_gamma[i,i] <- SE_D[i]^2
        V_Gamma[i,i] <- SE_Y[i]^2
      }else{
        V_gamma[i,j] <- V_gamma[j,i] <- ITT_D[i]*ITT_D[j] / n1
        V_Gamma[i,j] <- V_Gamma[j,i] <- ITT_Y[i]*ITT_Y[j] / n2
      }
    }
  }

  TSHT.out <- TSHT.sumstats(ITT_D, ITT_Y, SE_D, SE_Y, n1, n2)
  V0.hat = sort(TSHT.out$VHat)
  ## Construct range [L, U]
  if(is.vector(CI.init)){
    CI.init.union = t(as.matrix(sort(CI.init)))
  }else{
    ## current method to select initial [L, U]
    var.beta = diag(V_Gamma)/ITT_D^2 + diag(V_gamma)*ITT_Y^2/ITT_D^4
    var.beta = var.beta[V0.hat]
    CI.init = matrix(NA, nrow=length(V0.hat), ncol=2)
    CI.init[,1] = (ITT_Y/ITT_D)[V0.hat] - sqrt(log(n1)*var.beta)
    CI.init[,2] = (ITT_Y/ITT_D)[V0.hat] + sqrt(log(n1)*var.beta)
    uni = Intervals(CI.init)
    CI.init.union = as.matrix(interval_union(uni))
  }

  # Construct beta.grid
  beta.grid = grid.CI(CI.init.union, grid.size=n1^{-a})

  if(Sampling){
    ## Sampling Method
    CI.sampling = Searching.CI.sampling(ITT_D, ITT_Y, V_gamma, V_Gamma, InitiSet=V0.hat,
                                        beta.grid = beta.grid, rho=rho, M=M, prop=prop, filtering=TRUE)
    CI=CI.sampling$CI
    rule=CI.sampling$rule
  }else{
    ## Searching Method
    CI.searching = Searching.CI(ITT_D, ITT_Y, V_gamma, V_Gamma, InitiSet = V0.hat,
                                beta.grid = beta.grid)
    CI=CI.searching$CI
    rule=CI.searching$rule
  }
  returnList <- list(CI=CI, rule=rule, VHat=V0.hat, CI.init=CI.init.union, TSHT.out = TSHT.out)

  return(returnList)
}


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


########## Helpful Functions ##########

var_TSHT <- function(gamma, Gamma, V_gamma, V_Gamma){

  Var_Gg <- t(gamma) %*% V_Gamma %*% gamma + t(Gamma) %*% V_gamma %*% Gamma
  Var_gg <- 4 * t(gamma) %*% V_gamma %*% gamma
  Cov_Gg <- 2 * t(Gamma) %*% V_gamma %*% gamma

  E_Gg <- sum(gamma*Gamma)
  E_gg <- sum(gamma^2) + sum(diag(V_gamma))

  sig2 <- (E_Gg/E_gg)^2 * (Var_Gg/E_Gg^2 + Var_gg / E_gg^2 - 2*Cov_Gg/(E_Gg*E_gg) )

  return(sig2)
}


grid.CI <- function(CI.matrix, grid.size){
  d = dim(CI.matrix)[1]
  grid.seq = NULL
  for(l in 1:d) grid.seq = c(grid.seq, seq(CI.matrix[l, 1], CI.matrix[l, 2], by=grid.size))
  return(grid.seq)
}


Searching.CI <- function(ITT_D, ITT_Y, V_gamma, V_Gamma, InitiSet, beta.grid){
  threshold.size = length(InitiSet)/2
  n.beta = length(beta.grid)
  pz = dim(V_Gamma)[1]

  ## new rho method
  Tn = qnorm(1-0.05/(2*pz))

  ## valid grid
  valid.grid = rep(NA, n.beta)
  for(j in 1:n.beta){
    b = beta.grid[j]
    temp = sqrt(diag(V_Gamma + b^2*V_gamma))
    valid.grid[j] = sum(abs(ITT_Y[InitiSet] - b*ITT_D[InitiSet]) < Tn*temp[InitiSet])
  }

  ## select beta
  if(length(beta.grid[which(valid.grid > threshold.size)])==0){
    rule=FALSE
    warning("Rule Fails. SS will give misleading CIs, SEs, and p-values.")
    sel.index = which(valid.grid==max(valid.grid))
  }else{
    rule=TRUE
    sel.index = which(valid.grid>threshold.size)
  }
  CI = t(as.matrix(c(min(beta.grid[sel.index]), max(beta.grid[sel.index]))))

  return(list(CI=CI, rule=rule))
}


Searching.CI.sampling <- function(ITT_D, ITT_Y, V_gamma, V_Gamma, InitiSet,
                                  beta.grid, rho=NULL, M=1000, prop=0.1, filtering=TRUE){
  threshold.size = length(InitiSet)/2
  n.beta = length(beta.grid)
  pz = dim(V_Gamma)[1]

  ## new rho method
  Tn = qnorm(1-0.05/(2*pz))

  ## Covariance Matrix
  Cov1 <- cbind(V_Gamma, matrix(0, pz, pz))
  Cov2 <- cbind(matrix(0, pz, pz), V_gamma)
  Cov.total <- rbind(Cov1,Cov2)

  valid.grid.sample = matrix(NA, nrow=M, ncol=n.beta)
  Gen.mat = MASS::mvrnorm(M, rep(0, 2*pz), Cov.total)
  if(is.null(rho)) rho = (log(n1)/M)^(1/(2*length(InitiSet)))/6 # initial rho if not specified

  if(filtering){
    temp = abs(t(t(Gen.mat[,(pz+1):(2*pz)]) / sqrt(diag(V_gamma))))
    temp = apply(temp, MARGIN=1, FUN=max)
    temp = (temp <= qnorm(1-0.05/(2*pz)))
    Gen.mat = Gen.mat[temp,]
    M = sum(temp)
  }

  while(rho < 0.5){
    for(m in 1:M){
      ITT_Y.sample = ITT_Y - Gen.mat[m, 1:pz]
      ITT_D.sample = ITT_D - Gen.mat[m, (pz+1):(2*pz)]
      for(j in 1:n.beta){
        b = beta.grid[j]
        temp = sqrt(diag(V_Gamma + b^2*V_gamma))
        valid.grid.sample[m, j] = sum(abs(ITT_Y.sample[InitiSet] - b*ITT_D.sample[InitiSet])<
                                        rho*Tn*temp[InitiSet])
      }
    }
    CI = matrix(NA, nrow=M, ncol=2)
    for(m in 1:M){
      if(length(which(valid.grid.sample[m, ] > threshold.size))>0){
        CI[m, 1] = min(beta.grid[which(valid.grid.sample[m, ] > threshold.size)])
        CI[m, 2] = max(beta.grid[which(valid.grid.sample[m, ] > threshold.size)])
      }
    }
    CI = CI[!rowSums(is.na(CI)), , drop=FALSE] # CI = na.omit(CI)

    ## check CI dims, stop iterations if CI dim is big enough
    if(dim(as.matrix(CI))[1] >= prop*M) break
    rho = 1.25 * rho # increase rho with iterations
  }

  rule = TRUE
  if(dim(as.matrix(CI))[1] < prop*M){
    warning("Sampling Criterion not met, trasfer to Searching Method.")
    CI.searching = Searching.CI(ITT_Y, ITT_D, V_Gamma, V_gamma, InitiSet, beta.grid)
    rule = CI.searching$rule
    CI = CI.searching$CI
  }else{
    uni = Intervals(CI)
    CI = as.matrix(interval_union(uni))
    CI = t(as.matrix(c(min(CI[,1]), max(CI[,2]))))
  }

  return(list(CI=CI, rule=rule))
}



















