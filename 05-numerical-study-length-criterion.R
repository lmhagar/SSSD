## confirmation analysis for Table 1 in the paper
require(foreach)
require(doParallel)
require(doSNOW)
require(MASS)
require(mvtnorm)

## helper function that outputs HDI length
outputLen <- function(alpha, sampleSize, alpha1, beta1, alpha2, beta2, mu1_0 = 2, tau1_0 = 0.1, 
         kappa1_0 = 2, lambda1_0 = 0.1, mu2_0 = 2, tau2_0 = 0.1, 
         kappa2_0 = 2, lambda2_0 = 0.1, symmetric = FALSE){

  ## alpha: coverage for HDI
  ## sampleSize: sample size n that determines how many observations are randomly generated from each group
  ## alpha1, beta1, alpha2, and beta2 are the relevant design values
  ## alpha_j has a Gamma(muj_0, tauj_0) prior, where tauj_0 is a rate
  ## beta_j has a Gamma(kappaj_0, lambdaj_0) prior, where lambdaj_0 is a rate
  ## symmetric is a binary variable to denote whether the credible inteval should be symmetric (FALSE for HDI)
  
  y1 <- rgamma(sampleSize, alpha1, beta1)
  y2 <- rgamma(sampleSize, alpha2, beta2)
  
  theta.diff <- GammaPostSIR(1000000, 100000, y1, y2, alpha1, beta1, alpha2, beta2, mu1_0, tau1_0,
                             kappa1_0, lambda1_0, mu2_0, tau2_0, kappa2_0, lambda2_0,
                             tau = 4.29)
  theta.diff[which(!is.finite(theta.diff))] <- 10^7
  theta.diff[which(is.na(theta.diff))] <- 10^7
  
  if (symmetric == TRUE)
  {
    avg.length <- unname(quantile(theta.diff, 1 - alpha/2) - quantile(theta.diff, alpha/2))
  }
  else
  {
    optInfo <- optimize(intervalLength, c(0, alpha),
                        data = theta.diff, alpha = alpha)
    avg.length <- unname(optInfo$objective)
  }
  
  return(avg.length)
}

## helper function to optimize to find HDI length
intervalLength <- function(data, base, alpha)
{
  quantile(data, base + (1 - alpha)) - quantile(data, base)
}

## helper function to approximate the posterior of theta = log(theta_1) - log(theta_2)
## similar to the function in files 02 and 03, but the output is different
GammaPostSIR <- function(k, l, y1, y2, alpha1_0, beta1_0, alpha2_0, beta2_0, mu1, tau1, kappa1, lambda1,
                         mu2, tau2, kappa2, lambda2, metric = "tail", tau = 4.29){

  ## k is the number of points generated from sampling distribution
  ## l is the number of resampled points returned by algorithm
  ## y1 and y2 are the food expenditure observations in each group
  ## alphaj_0 and betaj_0 are the design values
  ## alpha_j has a Gamma(mu_j, tau_j) prior, where tau_j is a rate
  ## beta_j has a Gamma(kappa_j, lambda_j) prior, where lambda_j is a rate
  ## metric can be "tail" (tail probability), "mean" (mean), "cv" (coefficient
  ## of variation), or a number between 0 and 1 that denotes a quantile

  ## helper function to return value that is proportional to the posterior likelihood;
  ## used to create importance sampling weights
  propLog <- function(alpha, beta, dlogprod, dsum, mu, tau, kappa, lambda, n){
    return(ifelse(alpha > 0 & beta > 0,
                  (alpha*n + kappa - 1)*log(beta) - n*lgamma(alpha) + (alpha - 1)*dlogprod
                  - beta*(dsum + lambda) - tau*alpha + (mu - 1)*log(alpha), -Inf))
  }

  n1 <- length(y1)
  n2 <- length(y2)
  sum1 <- sum(y1)
  sum2 <- sum(y2)
  logprod1 <- sum(log(y1))
  logprod2 <- sum(log(y2))

  ## use method of moments as starting point to solve for gamma MLEs
  if (n1 > 1){
    momalpha <- mean(y1)^2/var(y1)
    mombeta <- var(y1)/mean(y1)

    gmll <- function(theta,datta)
    {
      a <- theta[1]; b <- theta[2]
      n <- length(datta); sumd <- sum(datta); sumlogd <- sum(log(datta))
      gmll <- n*a*log(b) + n*lgamma(a) + sumd/b - (a-1)*sumlogd
      gmll
    }

    ## solve MLEs using nlm and convert the scale to a rate (for beta)
    gammasearch <- tryCatch(nlm(gmll,c(momalpha,mombeta),hessian=T,datta=y1)$estimate,
                            error = function(e) {c(-1,-1)})

    if (gammasearch[1] > 0 & gammasearch[2] > 0){
      alpha1_0 <- gammasearch[1]
      beta1_0 <- 1/gammasearch[2]
    }
  }

  ## specify covariance matrix for first group and generate sample from proposal
  mat1 <- solve(matrix(c(trigamma(alpha1_0), -1/beta1_0, -1/beta1_0, alpha1_0/beta1_0^2), nrow = 2))
  samp1 <- mvrnorm(n = k, mu = c(alpha1_0, beta1_0), Sigma = 1.5*mat1/n1)

  ## form importance sampling weights on log-scale
  w1 <- propLog(samp1[,1], samp1[,2], sum(log(y1)), sum(y1), mu1, tau1, kappa1, lambda1, length(y1)) -
    dmvnorm(samp1, mean = c(alpha1_0, beta1_0), sigma = 1.5*mat1/n1, log = TRUE)
  w1 <- exp(w1 - max(w1))
  w1 <- w1/sum(w1)
  ## resample to create approximate sample from posterior
  inds1 <- sample(seq(1,k,by = 1), size = l, replace = TRUE, prob = w1)
  samp1_post <- samp1[inds1,]
  alpha1 <- samp1_post[,1]
  beta1 <- samp1_post[,2]

  if (n2 > 1){
    ## use method of moments as starting point to solve for gamma MLEs
    momalpha <- mean(y2)^2/var(y2)
    mombeta <- var(y2)/mean(y2)

    ## solve MLEs using nlm and convert the scale to a rate (for beta)
    gammasearch <- tryCatch(nlm(gmll,c(momalpha,mombeta),hessian=T,datta=y2)$estimate,
                            error = function(e) {c(-1,-1)})

    if (gammasearch[1] > 0 & gammasearch[2] > 0){
      alpha2_0 <- gammasearch[1]
      beta2_0 <- 1/gammasearch[2]
    }
  }

  ## specify covariance matrix for first group and generate sample from proposal
  mat2 <- solve(matrix(c(trigamma(alpha2_0), -1/beta2_0, -1/beta2_0, alpha2_0/beta2_0^2), nrow = 2))
  samp2 <- mvrnorm(n = k, mu = c(alpha2_0, beta2_0), Sigma = 1.5*mat2/n2)

  ## form importance sampling weights on log-scale
  w2 <- propLog(samp2[,1], samp2[,2], sum(log(y2)), sum(y2), mu2, tau2, kappa2, lambda2, length(y2)) -
    dmvnorm(samp2, mean = c(alpha2_0, beta2_0), sigma = 1.5*mat2/n2, log = TRUE)
  w2 <- exp(w2 - max(w2))
  w2 <- w2/sum(w2)
  ## resample to create approximate sample from posterior
  inds2 <- sample(seq(1,k,by = 1), size = l, replace = TRUE, prob = w2)
  samp2_post <- samp2[inds2,]
  alpha2 <- samp2_post[,1]
  beta2 <- samp2_post[,2]

  ## output vectors of interest
  if (metric == "tail"){
    theta1 <- 1 - pgamma(tau, alpha1, beta1)
    theta2 <- 1 - pgamma(tau, alpha2, beta2)
  }
  else if (metric == "mean"){
    theta1 <- alpha1/beta1
    theta2 <- alpha2/beta2
  }
  else if (metric == "cv"){
    theta1 <- 1/sqrt(alpha1)
    theta2 <- 1/sqrt(alpha2)
  }
  else {
    theta1 <- qgamma(0.01*metric, alpha1, beta1)
    theta2 <- qgamma(0.01*metric, alpha2, beta2)
  }
  ## return posterior draws needed to compute figure
  return(log(theta1) - log(theta2))
}

## parameters for conviction thresholds, design values, etc
## as defined in paper
convictions <- c(0.5, 0.9, 0.8)
powers <- c(0.6, 0.7, 0.8)
equivalences <- c(0.25, 0.3, 0.15)
alphas <- 1 - convictions

informs <- read.csv("informs_aguas.csv")

## get hyperparameters for uninformative and informative settings
mu1s <- c(2, informs[1,1])
tau1s <- c(0.1, informs[1,2])
kappa1s <- c(2, informs[2,1])
lambda1s <- c(0.1, informs[2,2])
mu2s <- c(2, informs[3,1])
tau2s <- c(0.1, informs[3,2])
kappa2s <- c(2, informs[4,1])
lambda2s <- c(0.1, informs[4,2])

alp1 <- 2.11
bet1 <- 0.69
alp2 <- 2.43
bet2 <- 0.79

## parallelize simulations
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

registerDoSNOW(cl)
pb <- txtProgressBar(max = 10000, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## implement simulation for each combination of conviction threshold, power, delta_* (i)
## and each prior specification (k) 100 times for each of the 100 estimated SSSDs from file 04
## output results from each percentile to a .csv file
for (i in 1:3){
  for (k in 1:2){
    params <- NULL
    params$mu <- rep(read.csv(paste0("ba_results_rest_", k, i, ".csv"))$mu_final, each = 100)
    params$sigma <- rep(read.csv(paste0("ba_results_rest_", k, i, ".csv"))$sigma_final, each = 100)
    Sys.time()
    finalMatrix <- foreach(j=1:10000, .combine=rbind, .packages = c("MASS", "mvtnorm"),
                           .options.snow=opts) %dopar% {
                             ## if sample size is negative, deem length criterion unsatisfied 
                             ## (does not occur for log(theta_1) - log(theta_2))
                             if (ceiling(qnorm(0.05, params$mu[j], params$sigma[j])) < 1){
                               tempMatrix = 1000
                             } else {tempMatrix = outputLen(alphas[i],max(1,ceiling(qnorm(0.05, params$mu[j], params$sigma[j]))),
                                                            alp1, bet1,alp2, bet2,
                                                            mu1s[k], tau1s[k], kappa1s[k], lambda1s[k], mu2s[k], tau2s[k],
                                                            kappa2s[k], lambda2s[k]) #calling a function
                             }
                             
                             tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
                           }
    new_samples <- as.numeric(finalMatrix)
    write.csv(t(new_samples),paste0("ba_gamma_p05_matrix_adj", k, i, ".csv"), row.names = FALSE)
    Sys.time()
    
    Sys.time()
    finalMatrix <- foreach(j=1:10000, .combine=rbind, .packages = c("MASS", "mvtnorm"),
                           .options.snow=opts) %dopar% {
                             tempMatrix = outputLen(alphas[i],max(1,ceiling(qnorm(0.1, params$mu[j], params$sigma[j]))),
                                                    alp1, bet1,alp2, bet2,
                                                    mu1s[k], tau1s[k], kappa1s[k], lambda1s[k], mu2s[k], tau2s[k],
                                                    kappa2s[k], lambda2s[k]) #calling a function
                             
                             tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
                           }
    new_samples <- as.numeric(finalMatrix)
    write.csv(t(new_samples),paste0("ba_gamma_p10_matrix_adj", k, i, ".csv"), row.names = FALSE)
    Sys.time()
    
    Sys.time()
    finalMatrix <- foreach(j=1:10000, .combine=rbind, .packages = c("MASS", "mvtnorm"),
                           .options.snow=opts) %dopar% {
                             tempMatrix = outputLen(alphas[i],max(1,ceiling(qnorm(0.25, params$mu[j], params$sigma[j]))),
                                                    alp1, bet1,alp2, bet2,
                                                    mu1s[k], tau1s[k], kappa1s[k], lambda1s[k], mu2s[k], tau2s[k],
                                                    kappa2s[k], lambda2s[k]) #calling a function
                             
                             tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
                           }
    new_samples <- as.numeric(finalMatrix)
    write.csv(t(new_samples),paste0("ba_gamma_p25_matrix_adj", k, i, ".csv"), row.names = FALSE)
    Sys.time()
    
    Sys.time()
    finalMatrix <- foreach(j=1:10000, .combine=rbind, .packages = c("MASS", "mvtnorm"),
                           .options.snow=opts) %dopar% {
                             tempMatrix = outputLen(alphas[i],max(1,ceiling(qnorm(0.50, params$mu[j], params$sigma[j]))),
                                                    alp1, bet1,alp2, bet2,
                                                    mu1s[k], tau1s[k], kappa1s[k], lambda1s[k], mu2s[k], tau2s[k],
                                                    kappa2s[k], lambda2s[k]) #calling a function
                             
                             tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
                           }
    new_samples <- as.numeric(finalMatrix)
    write.csv(t(new_samples),paste0("ba_gamma_p50_matrix_adj", k, i, ".csv"), row.names = FALSE)
    Sys.time()
    
    Sys.time()
    finalMatrix <- foreach(j=1:10000, .combine=rbind, .packages = c("MASS", "mvtnorm"),
                           .options.snow=opts) %dopar% {
                             tempMatrix = outputLen(alphas[i],max(1,ceiling(qnorm(0.75, params$mu[j], params$sigma[j]))),
                                                    alp1, bet1,alp2, bet2,
                                                    mu1s[k], tau1s[k], kappa1s[k], lambda1s[k], mu2s[k], tau2s[k],
                                                    kappa2s[k], lambda2s[k]) #calling a function
                             
                             tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
                           }
    new_samples <- as.numeric(finalMatrix)
    write.csv(t(new_samples),paste0("ba_gamma_p75_matrix_adj", k, i, ".csv"), row.names = FALSE)
    Sys.time()
    
    Sys.time()
    finalMatrix <- foreach(j=1:10000, .combine=rbind, .packages = c("MASS", "mvtnorm"),
                           .options.snow=opts) %dopar% {
                             tempMatrix = outputLen(alphas[i],max(1,ceiling(qnorm(0.9, params$mu[j], params$sigma[j]))),
                                                    alp1, bet1,alp2, bet2,
                                                    mu1s[k], tau1s[k], kappa1s[k], lambda1s[k], mu2s[k], tau2s[k], 
                                                    kappa2s[k], lambda2s[k]) #calling a function
                             
                             tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
                           }
    new_samples <- as.numeric(finalMatrix)
    write.csv(t(new_samples),paste0("ba_gamma_p90_matrix_adj", k, i, ".csv"), row.names = FALSE)
    Sys.time()
    
    Sys.time()
    finalMatrix <- foreach(j=1:10000, .combine=rbind, .packages = c("MASS", "mvtnorm"),
                           .options.snow=opts) %dopar% {
                             tempMatrix = outputLen(alphas[i],max(1,ceiling(qnorm(0.95, params$mu[j], params$sigma[j]))),
                                                    alp1, bet1,alp2, bet2,
                                                    mu1s[k], tau1s[k], kappa1s[k], lambda1s[k], mu2s[k], tau2s[k], 
                                                    kappa2s[k], lambda2s[k]) #calling a function
                             
                             tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
                           }
    new_samples <- as.numeric(finalMatrix)
    write.csv(t(new_samples),paste0("ba_gamma_p95_matrix_adj", k, i, ".csv"), row.names = FALSE)
    Sys.time()
  }
}

## compute values in Table 1 
k = 2
i = 1
len <- rep(read.csv(paste0("ba_results_rest_", k, i, ".csv"))$criterion_l, each = 100)
paste0(mean(ifelse(read.csv(paste0("ba_gamma_p05_matrix_adj",k,i,".csv"))[1,] <= len,1,0)), " & ",
       mean(ifelse(read.csv(paste0("ba_gamma_p10_matrix_adj",k,i,".csv"))[1,] <= len,1,0))," & ",
       mean(ifelse(read.csv(paste0("ba_gamma_p25_matrix_adj",k,i,".csv"))[1,] <= len,1,0))," & ",
       mean(ifelse(read.csv(paste0("ba_gamma_p50_matrix_adj",k,i,".csv"))[1,] <= len,1,0))," & ",
       mean(ifelse(read.csv(paste0("ba_gamma_p75_matrix_adj",k,i,".csv"))[1,] <= len,1,0))," & ",
       mean(ifelse(read.csv(paste0("ba_gamma_p90_matrix_adj",k,i,".csv"))[1,] <= len,1,0))," & ",
       mean(ifelse(read.csv(paste0("ba_gamma_p95_matrix_adj",k,i,".csv"))[1,] <= len,1,0)))