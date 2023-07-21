require(foreach)
require(doParallel)
require(doSNOW)
require(numDeriv)
require(MASS)
require(mvtnorm)
require(ggplot2)
require(cowplot)
require(ggpubr)

## file to approximate the power curves using simulated data (red curves in
## Figure 2 of the paper)
## settings for sampling-resampling algorithm; k_app is number of sampled
## points and l_app is the number of resampled points
k_app <- 750000; l_app <- 90000

## settings for simulation (values to compute power)
rrr_lower <- c(50,100,100)
rrr_upper <- c(250,600,1750)
rrr_inc <- c(10, 20, 50)

## settings for three gamma scenarios defined in paper
convictions <- c(0.5, 0.9, 0.8)
powers <- c(0.6, 0.7, 0.8)
equivalences <- c(0.25, 0.3, 0.15)

gamma_alpha.1 <- 2.11
gamma_beta.1 <- 0.69
gamma_alpha.2 <- 2.43
gamma_beta.2 <- 0.79
metric <- "tail"
threshold <- 4.29

## read in informative hyperparameters
informs <- read.csv("informs_aguas.csv")

## function to do sampling-resampling for the posterior of theta1/theta2; this is nearly the same as
## the function from file 05 but the output is for the posterior of theta, NOT log(theta)
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
  theta <- theta1/theta2
  theta <- ifelse(is.na(theta), Inf, theta)
  return(theta)
}

## set up parallelization
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

registerDoSNOW(cl)
pb <- txtProgressBar(max = 10000, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## output power curve for informative priors
for (k in 1:3){
  rrr <- seq(rrr_lower[k], rrr_upper[k], rrr_inc[k])
  if (k == 1){
    ## add more granularity for smallest sample size setting
    rrr <- c(seq(10,45,5), rrr)
  }
  results_rrr <- rep(0, length(rrr))
  
  ## obtain parameters for particular setting
  eq <- equivalences[k]
  alpha <- 1 - convictions[k]
  pwer <- convictions[k]
  
  for (i in 1:length(rrr)){
    print(rrr[i])
    
    pwr_rep <- foreach(j=1:10000, .combine='rbind', .packages = c("MASS", "mvtnorm"),
                         .options.snow=opts) %dopar% {
                           
                           ## simulate data; approximate posterior of theta1/theta2; compute PoA
                           y_star1 <- rgamma(rrr[i], gamma_alpha.1, gamma_beta.1)
                           y_star2 <- rgamma(rrr[i], gamma_alpha.2, gamma_beta.2)
                           
                           theta.diff <- GammaPostSIR(k_app, l_app, y_star1, y_star2, gamma_alpha.1, gamma_beta.1,
                                                      gamma_alpha.2, gamma_beta.2, informs[1,1], informs[1,2],
                                                      informs[2,1], informs[2,2], informs[3,1], informs[3,2],
                                                      informs[4,1], informs[4,2],
                                                      metric, threshold)
                           
                           mean(ifelse(theta.diff > (1 + eq)^(-1), theta.diff < (1 + eq),0))
                           
                         }
    
    pwr_pwr3 = as.numeric(pwr_rep)
    results_rrr[i] <- mean(ifelse(pwr_rep>=pwer,1,0))
    write.csv(data.frame(rrr = rrr, results_rrr = results_rrr), paste0("confirm_results_2", k, ".csv"), row.names = FALSE)
  }
}

## output power curve for uninformative priors
for (k in 1:3){
  rrr <- seq(rrr_lower[k], rrr_upper[k], rrr_inc[k])
  if (k == 1){
    ## add more granularity for smallest sample size setting
    rrr <- c(seq(10,45,5), rrr)
  }
  results_rrr <- rep(0, length(rrr))
  
  ## obtain parameters for particular setting
  eq <- equivalences[k]
  alpha <- 1 - convictions[k]
  pwer <- convictions[k]
  
  for (i in 1:length(rrr)){
    print(rrr[i])
    
    pwr_rep <- foreach(j=1:10000, .combine='rbind', .packages = c("MASS", "mvtnorm"),
                         .options.snow=opts) %dopar% {
                           
                           ## simulate data; approximate posterior of theta1/theta2; compute PoA
                           y_star1 <- rgamma(rrr[i], gamma_alpha.1, gamma_beta.1)
                           y_star2 <- rgamma(rrr[i], gamma_alpha.2, gamma_beta.2)
                           
                           theta.diff <- GammaPostSIR(k_app, l_app, y_star1, y_star2, gamma_alpha.1, gamma_beta.1,
                                                      gamma_alpha.2, gamma_beta.2,
                                                      2, 0.1, 2, 0.1,
                                                      2, 0.1, 2, 0.1,
                                                      metric, threshold)
                           
                           mean(ifelse(theta.diff > (1 + eq)^(-1), theta.diff < (1 + eq),0))
                           
                         }
    
    pwr_rep <- as.numeric(pwr_rep)
    results_rrr[i] <- mean(ifelse(pwr_rep>=pwer,1,0))
    write.csv(data.frame(rrr = rrr, results_rrr = results_rrr), paste0("confirm_results_1", k, ".csv"), row.names = FALSE)
  }
}

## Make Figure 2 in the paper
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

powers <- c(0.6, 0.7, 0.8)
settings <- c("a", "b", "c")
jj <- 0
for (k in c(1,2,3)){
  jj <- jj + 1
  for (ii in 1:2){
    ## read in data from estimated SSSDs in file 04 and the confirmation
    ## simulations conducted earlier in this file
    res_1 <- read.csv(paste0("ba_samps_pwr_", ii, k, ".csv"))
    res_2 <- read.csv(paste0("ba_results_rest_", ii, k, ".csv"))
    rrr <- read.csv(paste0("confirm_results_",ii, k, ".csv"))$rrr
    results_rrr <- read.csv(paste0("confirm_results_", ii, k, ".csv"))$results_rrr
    power <- powers[jj]
    y_mat <- NULL
    z_mat <- NULL
    for (i in 1:100){
      mu_l <- res_2$mu_l[i]
      sigma_l <- res_2$sigma_l[i]
      mu_v3 <- res_2$mu_v3[i]
      mu_final <- res_2$mu_final[i]
      sigma_final <- res_2$sigma_final[i]
      rough_pwr <- res_2$rough_pwr[i]
      
      funecdf <- ecdf(as.numeric(res_1[i,]))
      
      ## compute scaling factor for inverse of power curve
      scl <- as.numeric(mu_v3/quantile(res_1[i,], rough_pwr))
      xxx <- seq(0, max(res_1[i,]), length.out = 600)
      
      ## adjust the sample sizes by the scaling factor
      y_mat <- rbind(y_mat, data.frame(n = matrix(scl*xxx, ncol = 1), curve = rep(i + 100, length(xxx))))
      z_mat <- rbind(z_mat, data.frame(power = matrix(funecdf(xxx), ncol = 1)))
      
    }
    
    ## these curves will be used to compute the blue portions of the power curves
    ## that are in the probable domain of the prior adjusted SSSD
    for (i in 1:100){
      mu_l <- res_2$mu_l[i]
      sigma_l <- res_2$sigma_l[i]
      mu_v3 <- res_2$mu_v3[i]
      mu_final <- res_2$mu_final[i]
      sigma_final <- res_2$sigma_final[i]
      rough_pwr <- res_2$rough_pwr[i]
      
      funecdf <- ecdf(as.numeric(res_1[i,]))
      
      scl <- as.numeric(mu_v3/quantile(res_1[i,], rough_pwr))
      xxx <- seq(qnorm(0.001, mu_final, sigma_final), qnorm(0.999, mu_final, sigma_final), length.out = 200)
      
      y_mat <- rbind(y_mat, data.frame(n = matrix(scl*xxx, ncol = 1), curve = rep(i + 200, length(xxx))))
      z_mat <- rbind(z_mat, data.frame(power = matrix(funecdf(xxx), ncol = 1)))
      
    }
    
    assign(paste0("power", ii, k), power)
    
    assign(paste0("data_full",ii, k), data.frame(n = y_mat$n, power = z_mat$power, curve = y_mat$curve))
    assign(paste0("data_full",ii, k), subset(get(paste0("data_full",ii, k)), get(paste0("data_full",ii, k))$n > 0))
    assign(paste0("data_full",ii, k), subset(get(paste0("data_full",ii, k)), get(paste0("data_full",ii, k))$power > 0))
    
    assign(paste0("data_sim",ii, k), data.frame(n = rrr, power = results_rrr))
    assign(paste0("title", ii, k), paste0("Setting ", ii, settings[jj]))
    
    assign(paste0("plot",ii, k), ggplot(get(paste0("data_full",ii, k)), aes(x=n)) + theme_bw() +
             geom_line(aes(y = power, color=as.factor(curve), alpha = 0.25, ), size = 1) +
             labs(title=get(paste0("title", ii, k))) +
             labs(x= bquote(italic(n)), y= bquote('Power')) +
             theme(plot.title = element_text(size=20,face="bold",
                                             margin = margin(t = 0, 0, 5, 0))) +
             theme(axis.text=element_text(size=16),
                   axis.title=element_text(size=18)) +
             theme(legend.position="none") +
             scale_color_manual(name = " ", 
                                values = c(rep("grey", 100), rep(cbbPalette[3], 100), "firebrick")) +
             theme(legend.text=element_text(size=18)) +
             ylim(0,1) + 
             xlim(0, max(get(paste0("data_sim",ii, k))$n)) + 
             theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0))) +
             theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0))) +
             geom_line(data = get(paste0("data_sim",ii, k)), aes(x = n, y = power, color="firebrick"), size = 2) +
             geom_hline(yintercept=get(paste0("power", ii, k)), lty = 2)
    )
  }
}

## compile subfigures in grid form
figp.row1 <- plot_grid(plot11 + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")), 
                       plot21 + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")),
                       rel_widths = c(1, 1))
figp.row2 <- plot_grid(plot12 + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")), 
                       plot22 + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")),
                       rel_widths = c(1, 1))
figp.row3 <- plot_grid(plot13 + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")), 
                       plot23 + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")),
                       rel_widths = c(1, 1))
figp <- plot_grid(figp.row1, figp.row2, figp.row3, nrow = 3)

# output as .pdf file for the article
pdf(file = "FigTwoBA.pdf",   # The directory you want to save the file in
    width = 10.5, # The width of the plot in inches (12.41)
    height = 8.65) # The height of the plot in inches (10.7)

figp

dev.off()