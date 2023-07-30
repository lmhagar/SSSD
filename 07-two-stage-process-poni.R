require(numDeriv)
require(MASS)
require(mvtnorm)
require(qrng)
require(nleqslv)

## this function implements the two-stage procedure to estimate the SSSD for equivalence tests facilitated
## via the probability of noninferiority with respect to group 2
sssdGammaPONI <- function(conviction, power, equivalence, gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2,  
                      length = NULL, alpha = NULL, mu0.1 = 2, tau0.1 = 0.1,
                      kappa0.1 = 2, lambda0.1 = 0.1, mu0.2 = 2, tau0.2 = 0.1, kappa0.2 = 2, lambda0.2 = 0.1,
                      threshold = NULL, rep_pwr = 1024, rep_var = 256){
  
  ## conviction: the conviction threshold in [0.5, 1), denoted gamma in the paper
  ## power: the target power, denoted capital Gamma in the paper
  ## equivalence: the extension to the equivalence margin for ratio-based comparisons; delta_* in the paper
  ## gamma_alpha.j: the design value for the shape parameter of the gamma distribution for group j = 1, 2
  ## gamma_beta.j: the design value for the rate parameter of the gamma distribution for group j = 1, 2
  ## length: used to specify a target HDI length that is not calibrated with target power; NULL by default
  ## alpha: used to specify a target HDI coverage that is not calibrated with target power; NULL by default
  ## alpha_j has a Gamma(mu0.j, tau0.j) prior, where tau0.j is a rate for group j = 1, 2
  ## beta_j has a Gamma(kappa0.j, lambda0.j) prior, where lambda0.j is a rate for group j = 1, 2
  ## threshold: the threshold for the tail probability in each group, denoted tau in the paper
  ## rep_pwr: length of Sobol' sequence to approximate power curve, denoted n_sob in the paper
  ## rep_var: length of Sobol' sequence to compute final power estimate, denoted n_var in the paper
  
  ## find good upper bound for the sample size to generate the asymptotic power curve
  ROPE_len <- (1 + equivalence) - (1 + equivalence)^{-1}
  
  alpha1 <- gamma_alpha.1; alpha2 <- gamma_alpha.2; beta1 <- gamma_beta.1; beta2 <- gamma_beta.2
  ## define theta metrics in terms of design values
  theta1 <- 1 - pgamma(threshold, alpha1, beta1)
  theta2 <- 1 - pgamma(threshold, alpha2, beta2)
  
  ## compute partial derivatives of theta with respect to alpha and beta for each group
  int_1 <- integrate(function(alpha, beta, x){dgamma(x, alpha, beta)*log(beta*x)},
                     lower = 0, upper = threshold, alpha = alpha1, beta = beta1)$value
  
  d_alpha1 <- digamma(alpha1)*pgamma(threshold,alpha1, beta1) - int_1
  
  d_beta1 <- (alpha1/beta1)*(pgamma(threshold, alpha1+1, beta1) - 
                               pgamma(threshold, alpha1, beta1))
  
  int_2 <- integrate(function(alpha, beta, x){dgamma(x, alpha, beta)*log(beta*x)},
                     lower = 0, upper = threshold, alpha = alpha2, beta = beta2)$value
  
  d_alpha2 <- digamma(alpha2)*pgamma(threshold, alpha2, beta2) - int_2
  
  d_beta2 <- (alpha2/beta2)*(pgamma(threshold, alpha2+1, beta2) -
                               pgamma(threshold, alpha2,beta2))
  
  ## compute Fisher information for each gamma model
  Fish1 <- matrix(c(trigamma(alpha1), -1/beta1, - 1/beta1, alpha1/beta1^2), nrow = 2)
  
  Fish2 <- matrix(c(trigamma(alpha2), -1/beta2, - 1/beta2, alpha2/beta2^2), nrow = 2)
  
  ## apply the delta method to get the limiting variance for each theta_j metric
  avar1 <- t(c(d_alpha1, d_beta1))%*%solve(Fish1)%*%c(d_alpha1, d_beta1)
  
  avar2 <- t(c(d_alpha2, d_beta2))%*%solve(Fish2)%*%c(d_alpha2, d_beta2)
  
  ## apply the delta method to get the limiting variance for theta = logtheta1 - logtheta2
  Fish_ratio_mu <- avar1/theta1^2 + avar2/theta2^2
  
  ## return initial upper bound for the root finding algorithm
  mu_start <- (4*qnorm(1 - (1 - conviction)/2)*sqrt(Fish_ratio_mu)/ROPE_len)^2
  
  ## generate a Sobol' sequence to find asymptotic power curve
  sob_pwr <- sobol(rep_pwr, d = 4, randomize = "digital.shift")
  
  ## function to input into simplified uniroot function
  ## u is point from Sobol' sequence, targetP is the target power
  ## for the root-finding algorithm to obtain, and n_val is the 
  ## sample size presently explored by the root-finding algorithm
  targetPower <- function(u, targetP, params, delta_star, n_val){
    
    ## return negative power if sample size is not positive
    if (n_val <= 0){return(-1.5)}
    
    gamma_alpha.1 <- params[1]
    gamma_beta.1 <- params[2]
    gamma_alpha.2 <- params[3]
    gamma_beta.2 <- params[4]
    threshold <- params[5]
    
    ## generate approximately normal MLEs for group 1 using delta method
    rho1 <- 1/sqrt(gamma_alpha.1*trigamma(gamma_alpha.1))
    mat1 <- matrix(c(1/gamma_alpha.1, 1/gamma_alpha.1, 1/gamma_alpha.1, 
                     trigamma(gamma_alpha.1))/(trigamma(gamma_alpha.1)*gamma_alpha.1 - 1), nrow = 2)
    a1 <- qnorm(u[1], log(gamma_alpha.1), sqrt(mat1[1,1]/n_val))
    b1 <- qnorm(u[2], log(gamma_beta.1) + rho1*(a1 - log(gamma_alpha.1))*(sqrt(mat1[2,2])/sqrt(mat1[1,1])), sqrt(1- rho1^2)*sqrt(mat1[2,2]/n_val))
    
    ## generate approximately normal MLEs for group 2 using delta method
    rho2 <- 1/sqrt(gamma_alpha.2*trigamma(gamma_alpha.2))
    mat2 <- matrix(c(1/gamma_alpha.2, 1/gamma_alpha.2, 1/gamma_alpha.2, 
                     trigamma(gamma_alpha.2))/(trigamma(gamma_alpha.2)*gamma_alpha.2 - 1), nrow = 2)
    a2 <- qnorm(u[3], log(gamma_alpha.2), sqrt(mat2[1,1]/n_val))
    b2 <- qnorm(u[4], log(gamma_beta.2) + rho2*(a2 - log(gamma_alpha.2))*(sqrt(mat2[2,2])/sqrt(mat2[1,1])), sqrt(1- rho2^2)*sqrt(mat2[2,2]/n_val))
    
    ## exponentiate MLEs
    a1 <- exp(a1); b1 <- exp(b1); a2 <- exp(a2); b2 <- exp(b2)
    
    ## ensure no MLEs underflow to 0
    if(max(a1 <= 0, b1 <= 0, a2 <= 0, b2<= 0)){return(-1.5)}
    
    ## define theta metrics in terms of design values
    theta1 <- 1 - pgamma(threshold, a1, b1)
    theta2 <- 1 - pgamma(threshold, a2, b2)
    
    ## compute partial derivatives of theta with respect to alpha and beta for each group
    int_1 <- tryCatch(integrate(function(alpha, beta, x){dgamma(x, alpha, beta)*log(beta*x)},
                                lower = 0, upper = threshold, alpha = a1, beta = b1)$value,
                      error = function(e) {-10001})
    
    ## say power criterion is not satisfied if the integral does not converge
    if (int_1 < -10000){return(-1.5)}
    
    d_alpha1 <- digamma(a1)*pgamma(threshold,a1, b1) - int_1
    
    d_beta1 <- (a1/b1)*(pgamma(threshold, a1+1, b1) - 
                          pgamma(threshold, a1, b1))
    
    int_2 <- tryCatch(integrate(function(alpha, beta, x){dgamma(x, alpha, beta)*log(beta*x)},
                                lower = 0, upper = threshold, alpha = a2, beta = b2)$value,
                      error = function(e) {-10001})
    
    if (int_2 < -10000){return(-1.5)}
    
    d_alpha2 <- digamma(a2)*pgamma(threshold, a2, b2) - int_2
    
    d_beta2 <- (a2/b2)*(pgamma(threshold, a2+1, b2) -
                          pgamma(threshold, a2,b2))
    
    ## compute Fisher information for each gamma model
    Fish1 <- matrix(c(trigamma(a1), -1/b1, - 1/b1, a1/b1^2), nrow = 2)
    
    Fish2 <- matrix(c(trigamma(a2), -1/b2, - 1/b2, a2/b2^2), nrow = 2)
    
    ## apply the delta method to get the limiting variance for each theta_j metric
    avar1 <- t(c(d_alpha1, d_beta1))%*%solve(Fish1)%*%c(d_alpha1, d_beta1)
    
    avar2 <- t(c(d_alpha2, d_beta2))%*%solve(Fish2)%*%c(d_alpha2, d_beta2)
    
    ## apply the delta method to get the limiting variance for theta = logtheta1 - logtheta2
    Fish_ratio_mu <- tryCatch(avar1/theta1^2 + avar2/theta2^2, error = function(e) {-10001}) 
    
    ## return negative power if division causes NA/Inf values
    if (is.na(Fish_ratio_mu)){return(-1.5)}
    if (Fish_ratio_mu < -10000){return(-1.5)}
    
    ## return power based on normal approximation induced by Bernstein-von Mises
    realP <- 1 - pnorm(log(1/(1 + delta_star)), log(theta1/theta2), sqrt(Fish_ratio_mu)/sqrt(n_val))
    
    ## return estimated power less target power (for root-finding algorithm)
    return(realP - targetP)
  }
  
  ## uniroot function without the error checking
  uu <- function(f, lower, upper, tol = 1e-4, maxiter =1000L, ...) {
    f.lower <- f(lower, ...)
    f.upper <- f(upper, ...)
    val <- .External2(stats:::C_zeroin2, function(arg) f(arg, ...),
                      lower, upper, f.lower, f.upper, tol, as.integer(maxiter))
    return(val[1])
  }
  
  ## find some better endpoints for the root-finding algorithm (closer to the 
  ## sample size that should be returned by the root-finding algorithm).
  ## we use a target power of 0 here (that is subtracted from the power at
  ## that sample size) to obtain the actual power
  endpoints1_vec <- NULL
  endpoints2_vec <- NULL
  endpoints3_vec <- NULL
  for (i in 1:nrow(sob_pwr)){
    endpoints1_vec[i] <- targetPower(targetP = 0, n_val = ceiling(0.5*mu_start),
                                     params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                                     delta_star = equivalence, u = sob_pwr[i,])
    endpoints2_vec[i] <- targetPower(targetP = 0, n_val = ceiling(mu_start),
                                     params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                                     delta_star = equivalence, u = sob_pwr[i,])
    endpoints3_vec[i] <- targetPower(targetP = 0, n_val = ceiling(2*mu_start),
                                     params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                                     delta_star = equivalence, u = sob_pwr[i,])
  }
  
  ## we categorize appropriate endpoints based on the smallest sample size where power > conviction
  endpoints_cat <- ifelse(endpoints1_vec >= conviction, 1,
                          ifelse(endpoints2_vec >= conviction, 2,
                                 ifelse(endpoints3_vec >= conviction, 3, 4)))
  
  ## find larger upper bound if the initial one is not sufficient
  ## (i.e., if some points still do not satisfy power criterion at 2*mu_start)
  last_group <- which(endpoints_cat == 4)
  if (length(last_group) == 0){
    upper_c <- 2
  } else{
    upper_c <- 2
    while(length(last_group) > 0){
      if (upper_c > 32){
        last_group <- NULL
      }
      upper_c <- 2*upper_c
      endpoints4_vec <- NULL
      for (i in 1:length(last_group)){
        endpoints4_vec[i] <- targetPower(targetP = 0, n_val = ceiling(upper_c*mu_start),
                                         params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                                         delta_star = equivalence, u = sob_pwr[i,])
      }
      keep_vec <- ifelse(endpoints4_vec >= conviction, FALSE, TRUE)
      ## only keep points that still do not satisfy power criterion after increasing
      ## the upper bound for the sample size
      last_group <- last_group[keep_vec]
    }
  }
  
  ## implement the root-finding algorithm for each point in the Sobol' sequence
  ## with different endpoints depending on the last section of code
  samps_pwr <- NULL
  for (i in 1:nrow(sob_pwr)){
    if (endpoints_cat[i] == 1){
      samps_pwr[i] <- uu(targetPower, lower =1, upper = ceiling(0.5*mu_start), targetP = conviction, 
                         params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                         delta_star = equivalence, u = sob_pwr[i,])
    }
    else if (endpoints_cat[i] == 2){
      samps_pwr[i] <- uu(targetPower, lower =ceiling(0.5*mu_start), upper = ceiling(mu_start), targetP = conviction, 
                         params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                         delta_star = equivalence, u = sob_pwr[i,])
    }
    else if (endpoints_cat[i] == 3){
      samps_pwr[i] <- uu(targetPower, lower =ceiling(mu_start), upper = ceiling(2*mu_start), targetP = conviction, 
                         params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                         delta_star = equivalence, u = sob_pwr[i,])
    }
    else{
      samps_pwr[i] <- uu(targetPower, lower =ceiling(2*mu_start), upper = ceiling(upper_c*mu_start), targetP = conviction, 
                         params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                         delta_star = equivalence, u = sob_pwr[i,])
    }
  }
  
  ## approximate the power curve with the empirical CDF of these sample sizes
  funecdf <- ecdf(samps_pwr)
  
  ## align coverage target with conviction threshold if not alternatively specified
  criterion_alpha <- 1 - conviction
  
  ## allow for specification of non-power calibrated target HDI coverage
  if (is.null(alpha)) {alpha <- criterion_alpha}
  
  ## find the Gamma-quantile of the approximate power curve to obtain mu_l
  ecdf_root <- function(quant, pwr){return(funecdf(quant) - pwr)}
  mu_l <- uu(ecdf_root, lower = quantile(samps_pwr, 0.5*power), upper = quantile(samps_pwr, power + 0.5*(1-power)), pwr = power)
  
  ## solve for l; rearrange (3.4) from main text
  criterion_l <- 2*qnorm(1 - criterion_alpha/2)*sqrt(Fish_ratio_mu)/sqrt(mu_l)
  
  ## allow for specification of non-power calibrated target HDI length
  if (is.null(length)) {length <- criterion_l}
  
  ## define helper function to compute partial derivatives of the square root
  ## of the inverse Fisher information (used to compute A as referenced in text)
  Fisher_ratio <- function(alpha1, beta1, alpha2, beta2, kappa){
    theta1 <- 1 - pgamma(kappa, alpha1, beta1)
    theta2 <- 1 - pgamma(kappa, alpha2, beta2)
    
    int_1 <- integrate(function(alpha, beta, x){dgamma(x, alpha, beta)*log(beta*x)},
                       lower = 0, upper = kappa, alpha = alpha1, beta = beta1)$value
    
    d_alpha1 <- digamma(alpha1)*pgamma(kappa,alpha1, beta1) - int_1
    
    d_beta1 <- (alpha1/beta1)*(pgamma(kappa, alpha1+1, beta1) - 
                                 pgamma(kappa, alpha1, beta1))
    
    int_2 <- integrate(function(alpha, beta, x){dgamma(x, alpha, beta)*log(beta*x)},
                       lower = 0, upper = kappa, alpha = alpha2, beta = beta2)$value
    
    d_alpha2 <- digamma(alpha2)*pgamma(kappa, alpha2, beta2) - int_2
    
    d_beta2 <- (alpha2/beta2)*(pgamma(kappa, alpha2+1, beta2) -
                                 pgamma(kappa, alpha2,beta2))
    
    Fish1 <- matrix(c(trigamma(alpha1), -1/beta1, - 1/beta1, alpha1/beta1^2), nrow = 2)
    
    Fish2 <- matrix(c(trigamma(alpha2), -1/beta2, - 1/beta2, alpha2/beta2^2), nrow = 2)
    
    avar1 <- t(c(d_alpha1, d_beta1))%*%solve(Fish1)%*%c(d_alpha1, d_beta1)
    
    avar2 <- t(c(d_alpha2, d_beta2))%*%solve(Fish2)%*%c(d_alpha2, d_beta2)
    
    ## asymptotic variance of theta = logtheta1 - logtheta2
    Fish_ratio <- avar1/theta1^2 + avar2/theta2^2
    return(sqrt(Fish_ratio))
  }
  
  ## define wrapper functions for each design value
  ## grad function requires a function with a single input for which the partial derivative
  ## will be taken with respect to that parameter
  Fisher_alpha1 <- function(alpha1){
    Fisher_ratio(alpha1, beta1 = beta1, alpha2 = alpha2, beta2 = beta2, kappa = threshold)
  }
  
  Fisher_beta1 <- function(beta1){
    Fisher_ratio(alpha1 = alpha1, beta1, alpha2 = alpha2, beta2 = beta2, kappa = threshold)
  }
  
  Fisher_alpha2 <- function(alpha2){
    Fisher_ratio(alpha1 = alpha1, beta1 = beta1, alpha2, beta2 = beta2, kappa = threshold)
  }
  
  Fisher_beta2 <- function(beta2){
    Fisher_ratio(alpha1 = alpha1, beta1 = beta1, alpha2 = alpha2, beta2, kappa = threshold)
  }
  
  ## compute the partial derivatives numerically
  dd_alpha1 <- grad(Fisher_alpha1, x = alpha1)
  dd_beta1 <- grad(Fisher_beta1, x = beta1)
  dd_alpha2 <- grad(Fisher_alpha2, x = alpha2)
  dd_beta2 <- grad(Fisher_beta2, x = beta2)
  
  ## apply the delta method again
  avar_dd1 <- t(c(dd_alpha1, dd_beta1))%*%solve(Fish1)%*%c(dd_alpha1, dd_beta1)
  
  avar_dd2 <- t(c(dd_alpha2, dd_beta2))%*%solve(Fish2)%*%c(dd_alpha2, dd_beta2)
  
  Fish_ratio_sd <- avar_dd1 + avar_dd2
  
  ## find standard deviation using the formula in the article
  sigma_l <- 4*qnorm(1-criterion_alpha/2)*sqrt(Fish_ratio_sd)/criterion_l
  
  ## function to solve for the posterior mode of each group
  ## x is the point (logalpha, logbeta) at which this function is evaluated at
  ## y = c(length(data), sum(log(data)), sum(data)); it is more efficient to
  ## compute this outside of the function that is to be optimized over
  ## other parameters are hyperparameters for the prior
  fn_grad <- function(x, y, mu, tau, kappa, lambda) {
    
    res1 <- exp(x[1])*y[1]*x[2] - y[1]*digamma(exp(x[1]))*exp(x[1]) + exp(x[1])*(y[2] - tau) + mu
    
    res2 <- y[1]*exp(x[1]) - (y[3] + lambda)*exp(x[2]) + kappa
    
    return(c(res1, res2))
  }
  
  ## calculate the covariance matrix for the normal approximation to the 
  ## posterior of the gamma model according to equation (4.2) on page 84
  ## of Gelman et al. (2013)
  calc_covar <- function(u, v, y_star, hyper){
    a <- exp(u); b <- exp(v); n <- length(y_star)
    tau <- hyper[1]; lambda <- hyper[2]
    d11 <- a*(n*digamma(a) + tau - n*v - sum(log(y_star))) + n*a^2*trigamma(a)
    d12 <- -1*a*n
    d22 <- b*(sum(y_star) + lambda)
    mat <- rbind(c(d11, d12), c(d12, d22))
    return(solve(mat))
  }
  
  ## helper function to obtain the mean and standard deviation of the normal approximation
  ## to the posterior of theta = log(theta_1) - log(theta_2)
  ## all input are defined as in targetPower() helper function; 
  ## hyper is a matrix of all hyperparameters (for the two groups)
  getNormal <- function(u, targetP, params, delta_star, n_val, hyper){
    
    if (n_val <= 0){return(-1.5)}
    
    gamma_alpha.1 <- params[1]
    gamma_beta.1 <- params[2]
    gamma_alpha.2 <- params[3]
    gamma_beta.2 <- params[4]
    threshold <- params[5]
    
    ## get approximately normal MLEs for logalpha and logbeta
    rho1 <- 1/sqrt(gamma_alpha.1*trigamma(gamma_alpha.1))
    mat1 <- matrix(c(1/gamma_alpha.1, 1/gamma_alpha.1, 1/gamma_alpha.1, 
                     trigamma(gamma_alpha.1))/(trigamma(gamma_alpha.1)*gamma_alpha.1 - 1), nrow = 2)
    a1 <- qnorm(u[1], log(gamma_alpha.1), sqrt(mat1[1,1]/n_val))
    b1 <- qnorm(u[2], log(gamma_beta.1) + rho1*(a1 - log(gamma_alpha.1))*(sqrt(mat1[2,2])/sqrt(mat1[1,1])), sqrt(1- rho1^2)*sqrt(mat1[2,2]/n_val))
    
    rho2 <- 1/sqrt(gamma_alpha.2*trigamma(gamma_alpha.2))
    mat2 <- matrix(c(1/gamma_alpha.2, 1/gamma_alpha.2, 1/gamma_alpha.2, 
                     trigamma(gamma_alpha.2))/(trigamma(gamma_alpha.2)*gamma_alpha.2 - 1), nrow = 2)
    a2 <- qnorm(u[3], log(gamma_alpha.2), sqrt(mat2[1,1]/n_val))
    b2 <- qnorm(u[4], log(gamma_beta.2) + rho2*(a2 - log(gamma_alpha.2))*(sqrt(mat2[2,2])/sqrt(mat2[1,1])), sqrt(1- rho2^2)*sqrt(mat2[2,2]/n_val))
    
    ## exponentiate to get MLEs
    a1 <- exp(a1); b1 <- exp(b1); a2 <- exp(a2); b2 <- exp(b2)
    
    if(max(a1 <= 0, b1 <= 0, a2 <= 0, b2<= 0)){return(-1.5)}
    
    gamma_alpha.1 <- a1; gamma_beta.1 <- b1; gamma_alpha.2 <- a2; gamma_beta.2 <- b2
    
    ## generate a non-random sample from each group
    LHS <- seq(0.5, ceiling(as.numeric(n_val)) - 0.5)/ceiling(as.numeric(n_val))
    y_star1 <- qgamma(LHS, gamma_alpha.1, gamma_beta.1)
    y_star2 <- qgamma(LHS, gamma_alpha.2, gamma_beta.2)
    
    ## summarize information from first group of data (faster computation)
    yy_star1 <- c(length(y_star1), sum(log(y_star1)), sum(y_star1))
    ## find posterior modes for the first group (logalpha and logbeta)
    modes <- nleqslv(log(c(gamma_alpha.1, gamma_beta.1)), fn_grad, y = yy_star1, mu = hyper[1,1], tau = hyper[1,2],
                     kappa = hyper[2,1], lambda = hyper[2,2] )$x
    
    ## create covariance matrix for approximately normal posterior that accounts for priors
    mat1_new <- calc_covar(modes[1], modes[2], y_star1, c(hyper[1,2], hyper[2,2]))
    ## exponentiate modes to return to standard scale
    modes1 <- exp(modes)
    
    ## repeat all steps for the second group
    yy_star2 <- c(length(y_star2), sum(log(y_star2)), sum(y_star2))
    modes <- nleqslv(log(c(gamma_alpha.2, gamma_beta.2)), fn_grad, y = yy_star2, mu = hyper[3,1], tau = hyper[3,2],
                     kappa = hyper[4,1], lambda = hyper[4,2] )$x
    
    mat2_new <- calc_covar(modes[1], modes[2], y_star2, c(hyper[3,2], hyper[4,2]))
    modes2 <- exp(modes)
    a1 <- modes1[1]; b1 <- modes1[2]; a2 <- modes2[1]; b2 <- modes2[2]
    
    ## ensure no modes are 0 due to underflow errors
    if(max(a1 <= 0, b1 <= 0, a2 <= 0, b2<= 0)){return(-1.5)}
    
    ## define theta metrics in terms of design values
    theta1 <- 1 - pgamma(threshold, a1, b1)
    theta2 <- 1 - pgamma(threshold, a2, b2)
    
    ## compute partial derivatives of theta with respect to logalpha and logbeta for each group
    ## this is different from the previous function that computes the derivatives with respect
    ## to alpha and beta
    int_1 <- tryCatch(integrate(function(alpha, beta, x){dgamma(x, alpha, beta)*log(beta*x)},
                                lower = 0, upper = threshold, alpha = a1, beta = b1)$value,
                      error = function(e) {-10001})
    
    if (int_1 < -10000){return(-1.5)}
    
    d_alpha1 <- (digamma(a1)*pgamma(threshold,a1, b1) - int_1)*a1
    
    d_beta1 <- ((a1/b1)*(pgamma(threshold, a1+1, b1) - 
                           pgamma(threshold, a1, b1)))*b1
    
    int_2 <- tryCatch(integrate(function(alpha, beta, x){dgamma(x, alpha, beta)*log(beta*x)},
                                lower = 0, upper = threshold, alpha = a2, beta = b2)$value,
                      error = function(e) {-10001})
    
    if (int_2 < -10000){return(-1.5)}
    
    d_alpha2 <- (digamma(a2)*pgamma(threshold, a2, b2) - int_2)*a2
    
    d_beta2 <- ((a2/b2)*(pgamma(threshold, a2+1, b2) -
                           pgamma(threshold, a2,b2)))*b2
  
    ## apply the delta method to get the limiting variance for each theta_j metric
    avar1 <- t(c(d_alpha1, d_beta1))%*%mat1_new%*%c(d_alpha1, d_beta1)
    
    avar2 <- t(c(d_alpha2, d_beta2))%*%mat2_new%*%c(d_alpha2, d_beta2)
    
    ## apply the delta method to get the limiting variance for theta = logtheta1 - logtheta2
    Fish_ratio_mu <- avar1/theta1^2 + avar2/theta2^2
    return(c(log(theta1/theta2), sqrt(Fish_ratio_mu)))
  }
  
  ## compute initial slope and intercept given in Theorem 1
  b0a <- -1*(qnorm(1 - alpha/2)*sqrt(Fish_ratio_mu)/mu_l^(0.5))
  b1a <- (qnorm(1 - alpha/2)*sqrt(Fish_ratio_mu)/mu_l^(1.5))
  ## print limiting mean and standard deviation
  print(c(mu_l, sigma_l))
  
  ## compute estimate for expected HDI length at mu_l (generate non-random sample)
  LHS <- seq(0.5, ceiling(as.numeric(mu_l)) - 0.5)/ceiling(as.numeric(mu_l))
  y_star1 <- qgamma(LHS, gamma_alpha.1, gamma_beta.1)
  y_star2 <- qgamma(LHS, gamma_alpha.2, gamma_beta.2)
  
  ## concatenate all hyperparameters in matrix form (used in getNormal() function)
  hyper_mat <- rbind(c(mu0.1, tau0.1), c(kappa0.1, lambda0.1), c(mu0.2, tau0.2), c(kappa0.2, lambda0.2))
  
  ## get fast normal approximation to the posterior for non-random sample
  normp1 <- getNormal(u = rep(0.5,4), targetP = conviction, delta_star = equivalence, hyper = hyper_mat,
                      params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold), n_val = ceiling(as.numeric(mu_l)))
  
  ## expected HDI length
  l_star1 <- 2*qnorm(1 - alpha/2)*normp1[2]
  
  ## find a closer sample size; denoted dot mu_l in the paper
  mu_v2 <- -1*(length-l_star1)/b1a + mu_l
  
  ## compute expected HDI length for second sample size
  LHS <- seq(0.5, ceiling(as.numeric(mu_v2)) - 0.5)/ceiling(as.numeric(mu_v2))
  y_star1 <- qgamma(LHS, gamma_alpha.1, gamma_beta.1)
  y_star2 <- qgamma(LHS, gamma_alpha.2, gamma_beta.2)
  
  normp2 <- getNormal(u = rep(0.5,4), targetP = conviction, delta_star = equivalence, hyper = hyper_mat,
                      params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold), n_val = ceiling(as.numeric(mu_v2)))
  ## expected HDI length
  l_star2 <- 2*qnorm(1 - alpha/2)*normp2[2]
  
  ## find a closer sample size (tilde mu_l from the paper is mu_v3)
  slope2 <- (l_star1-l_star2)/(mu_v2 - mu_l)

  ## this sample size is used to estimate the variability in HDI length and estimate power
  mu_v3 <- -1*(length-l_star2)/b1a + mu_v2
  
  ## generate a sample of approximately normal MLEs for logalpha and logbeta in each group
  rho1 <- 1/sqrt(gamma_alpha.1*trigamma(gamma_alpha.1))
  sob <- sobol(rep_var, d = 4, randomize = "digital.shift")
  mat1 <- matrix(c(1/gamma_alpha.1, 1/gamma_alpha.1, 1/gamma_alpha.1, 
                   trigamma(gamma_alpha.1))/(trigamma(gamma_alpha.1)*gamma_alpha.1 - 1), nrow = 2)
  a1 <- qnorm(sob[,1], log(gamma_alpha.1), sqrt(mat1[1,1]/ceiling(mu_v3)))
  b1 <- qnorm(sob[,2], log(gamma_beta.1) + rho1*(a1 - log(gamma_alpha.1))*(sqrt(mat1[2,2])/sqrt(mat1[1,1])), sqrt(1- rho1^2)*sqrt(mat1[2,2]/ceiling(mu_v3)))
  
  rho2 <- 1/sqrt(gamma_alpha.2*trigamma(gamma_alpha.2))
  mat2 <- matrix(c(1/gamma_alpha.2, 1/gamma_alpha.2, 1/gamma_alpha.2, 
                   trigamma(gamma_alpha.2))/(trigamma(gamma_alpha.2)*gamma_alpha.2 - 1), nrow = 2)
  a2 <- qnorm(sob[,3], log(gamma_alpha.2), sqrt(mat2[1,1]/ceiling(mu_v3)))
  b2 <- qnorm(sob[,4], log(gamma_beta.2) + rho2*(a2 - log(gamma_alpha.2))*(sqrt(mat2[2,2])/sqrt(mat2[1,1])), sqrt(1- rho2^2)*sqrt(mat2[2,2]/ceiling(mu_v3)))
  
  ## exponentiate to get MLEs on regular scale
  a1s <- exp(a1); b1s <- exp(b1); a2s <- exp(a2); b2s <- exp(b2)
  
  ## generate a non-random sample for each MLE
  ## compute the expected HDI length and power for each MLE
  var_results <- NULL
  for (i in 1:rep_var){
    var_temp <- getNormal(u = rep(0.5,4), targetP = conviction, delta_star = equivalence, hyper = hyper_mat,
                          params = c(a1s[i], b1s[i], a2s[i], b2s[i], threshold), n_val = ceiling(as.numeric(mu_v3)))
    
    var_results <- rbind(var_results, c(var_temp, 2*qnorm(1 - alpha/2)*var_temp[2], 
                                        1 - pnorm(log(1/(1 + equivalence)), var_temp[1], var_temp[2])))
  }
  
  ## obtain rough estimate for power
  l_var_pwr <- var_results[,4]
  rough_pwr <- mean(ifelse(l_var_pwr>=conviction,1,0)) 

  ## obtain estimate for variance of error terms (probit model)
  l_var_rep <- var_results[,3]
  l_sd <- sd(l_var_rep)
  
  ## calculate intercept of initial line used to find mu_v3
  int2 <- (length-l_star2) - slope2*mu_v2
  
  ## obtain penultimate estimates for SSSD parameters to get sample sizes at
  ## which to evaluate expected HDI length to fit linear model
  mu_pen <- -int2/slope2
  sigma_pen <- 10/(qnorm(pnorm(int2 + slope2*(mu_pen + 10),0,l_sd)))
  ## make adjustments if initial standard deviation estimate is unstable
  sigma_pen <- min(sigma_pen, sqrt(1.5)*sigma_l)
  
  ## find the 45th and 55th percentiles of this distribution
  mu_low <- qnorm(0.45, mu_pen, sigma_pen)
  mu_high <- qnorm(0.55, mu_pen, sigma_pen)
  ## add and subtract 5 from mean if estimated standard deviation is negative
  ## this would only occur if our estimated slope2 was negative (should not occur)
  if (sigma_pen < 0){
    mu_low <- mu_pen - 5
    mu_high <- mu_pen + 5
  }
  
  ## compute expected HDI length at three sample sizes (low, pen = penultimate mean, high)
  normplow <- getNormal(u = rep(0.5,4), targetP = conviction, delta_star = equivalence, hyper = hyper_mat,
                        params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold), n_val = ceiling(as.numeric(mu_low)))
  normppen <- getNormal(u = rep(0.5,4), targetP = conviction, delta_star = equivalence, hyper = hyper_mat,
                        params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold), n_val = ceiling(as.numeric(mu_pen)))
  normphigh <- getNormal(u = rep(0.5,4), targetP = conviction, delta_star = equivalence, hyper = hyper_mat,
                         params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold), n_val = ceiling(as.numeric(mu_high)))
  l_starlow <- 2*qnorm(1 - alpha/2)*normplow[2]
  l_starpen <- 2*qnorm(1 - alpha/2)*normppen[2]
  l_starhigh <- 2*qnorm(1 - alpha/2)*normphigh[2]
  
  ## fit a linear model to get coefficients
  l_y <- as.numeric(length) - c(as.numeric(l_starlow), as.numeric(l_starpen), as.numeric(l_starhigh), as.numeric(l_star2))
  l_x <- c(ceiling(mu_low), ceiling(mu_pen), ceiling(mu_high), ceiling(mu_v2))
  
  l_mod <- summary(lm(l_y ~ l_x))
  
  b1_final <- l_mod$coefficients[2,1]
  b0_final <- l_mod$coefficients[1,1]
  
  ## obtain final estimates for SSSD parameters
  mu_final <- -b0_final/b1_final
  sigma_final <- 10/(qnorm(pnorm(b0_final + b1_final*(mu_final + 10),0,l_sd)))
  print(c(mu_final,sigma_final))
  
  ## compute scaling factor to proportionally adjust the power curve for prior information
  scl <- as.numeric(mu_v3/quantile(samps_pwr, rough_pwr))
  x_plot <- seq(qnorm(0.001, mu_final, sigma_final), qnorm(0.999, mu_final, sigma_final), length.out = 200)
  y_plot <- funecdf(x_plot)
  
  ## obtain final sample size recommendation
  z_plot <- scl*x_plot
  fapp <- approxfun(x = y_plot, y = z_plot, ties = "ordered")
  samp_final <- fapp(power)
  
  ## compute final sample size recommendation for extreme percentiles
  if(is.na(samp_final)){
    if (rough_pwr > power){
      x_plot <- seq(0.01, qnorm(0.001, mu_final, sigma_final), length.out = 200)
      y_plot <- funecdf(x_plot)
      
      z_plot <- scl*x_plot
      fapp <- approxfun(x = y_plot, y = z_plot)
      samp_final <- fapp(power)
    }
    else{
      x_plot <- seq(qnorm(0.999, mu_final, sigma_final), 2*qnorm(0.999, mu_final, sigma_final), length.out = 200)
      y_plot <- funecdf(x_plot)
      
      z_plot <- scl*x_plot
      fapp <- approxfun(x = y_plot, y = z_plot)
      samp_final <- fapp(power)
    }
  }
  
  ## compute percentile of prior-adjusted SSSD
  pctl_final <- pnorm(samp_final, mu_final, sigma_final)
  print(pctl_final)
  print(samp_final)
  
  ## save results
  results_vec <- data.frame(criterion_l = length, mu_l = mu_l, sigma_l = sigma_l, mu_v3 = mu_v3, rough_pwr = rough_pwr,
                            mu_final = mu_final, sigma_final = sigma_final, pctl_final = pctl_final, 
                            samp_final = samp_final)
  
  return(list(samps_pwr, results_vec))
}

## define scenarios as in main text, but only run analysis for settings 1b and 2b (k = 2)
informs <- read.csv("informs_aguas.csv")
convictions <- c(0.5, 0.9, 0.8)
powers <- c(0.6, 0.7, 0.8)
equivalences <- c(0.25, 0.3, 0.15)
tic <- Sys.time()
for (k in 2){
  res_1 <- NULL
  res_2 <- NULL
  for (j in 1:100){
    print(j)
    res_temp <- sssdGammaPONI(conviction = convictions[k], power = powers[k], equivalence = equivalences[k], gamma_alpha.1 = 2.11, gamma_beta.1 = 0.69,
                              gamma_alpha.2 = 2.43, gamma_beta.2 = 0.79, threshold = 4.29, mu0.1 = informs[1,1], tau0.1 = informs[1,2],
                              kappa0.1 = informs[2,1], lambda0.1 = informs[2,2], mu0.2 = informs[3,1], tau0.2 = informs[3,2],
                              kappa0.2 = informs[4,1], lambda0.2 = informs[4,2], rep_pwr = 1024, rep_var = 256)
    
    res_1 <- rbind(res_1, res_temp[[1]])
    res_2 <- rbind(res_2, res_temp[[2]])
    
    write.csv(res_1, paste0("ba_poni_samps_pwr_2", k, ".csv"), row.names = FALSE)
    write.csv(res_2, paste0("ba_poni_results_rest_2", k, ".csv"), row.names = FALSE)
  }
}
toc <- Sys.time()
toc - tic

tic <- Sys.time()
for (k in 2){
  res_1 <- NULL
  res_2 <- NULL
  for (j in 1:100){
    print(j)
    res_temp <- sssdGammaPONI(conviction = convictions[k], power = powers[k], equivalence = equivalences[k], gamma_alpha.1 = 2.11, gamma_beta.1 = 0.69,
                              gamma_alpha.2 = 2.43, gamma_beta.2 = 0.79, threshold = 4.29, mu0.1 = 2, tau0.1 = 0.1,
                              kappa0.1 = 2, lambda0.1 = 0.1, mu0.2 = 2, tau0.2 = 0.1, kappa0.2 = 2, lambda0.2 = 0.1, rep_pwr = 1024, rep_var = 256)
    
    res_1 <- rbind(res_1, res_temp[[1]])
    res_2 <- rbind(res_2, res_temp[[2]])
    
    write.csv(res_1, paste0("ba_poni_samps_pwr_1", k, ".csv"), row.names = FALSE)
    write.csv(res_2, paste0("ba_poni_results_rest_1", k, ".csv"), row.names = FALSE)
  }
}
toc <- Sys.time()
toc - tic

## get the appropriate means for Table E2
for (i in 1:2){
  for (j in 2){
    dat <- read.csv(paste0("ba_poni_results_rest_",i,j,".csv"))
    print(c(i,j, round(mean(dat$mu_final),2), round(mean(dat$sigma_final),2), round(mean(dat$criterion_l),3)))
  }
}