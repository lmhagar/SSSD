## code to determine design values and informative priors
## this procedure involves processing data from the ENIGH 2018 survey

## BEGIN SETUP ##

## Process to obtain data file from ENIGH website
## Step 1: Go to this website: https://www.inegi.org.mx/programas/enigh/nc/2018/#Datos_abiertos
## Step 2: Click the "Download CSV" button to download a .zip file
## Step 3: Open .zip file
## Step 4: Open the "conjunto_de_datos_concentradohogar_enigh_2018_ns" subfolder
## Step 5: Open the "conjunto_de_datos" subfolder
## Step 6: Extract the lone .csv file: "conjunto_de_datos_concentradohogar_enigh_2018_ns.csv" and move it to the working directory

require(mvtnorm)
require(MASS)
## read data set
enigh <- read.csv("conjunto_de_datos_concentradohogar_enigh_2018_ns.csv")

## create a numeric character for region
enigh$region <- as.numeric(substr(as.character(enigh$ubica_geo),1,nchar(as.character(enigh$ubica_geo))-3))

## keep only the households from Aguascalientes (Region = 1):
aguas <- enigh[enigh$region %in% c(1),]

## keep only columns of interest
aguas <- aguas[c("region","factor","est_socio","educa_jefe", "tot_integ", "alimentos", "sexo_jefe")]

## convert column titles to English
names(aguas)[4] <- "education"
names(aguas)[5] <- "total_people"
names(aguas)[6] <- "total_food"
names(aguas)[7] <- "sex"

## remove households with 0 quarterly expenditure on food
aguas_full <- aguas ## save all data for later
aguas <- aguas[aguas$total_food > 0,]

## keep only individuals with estimated socioeconomic class 2
aguas <- aguas[aguas$est_socio ==2 ,]

## create simplified weighting factor 
aguas$factor2 <- round(aguas$factor/75)

## repeat the observations according to new weighting factor
aguas_long <- aguas[1,]
for (i in 1:nrow(aguas)){
  if (i %% 100 == 0){
    print(i)
  }
  for (j in 1:aguas$factor2[i]){
    aguas_long <- rbind(aguas_long, aguas[i,])
  }
}
aguas_long <- aguas_long[-1,]

## calculate food expense per person in thousands of pesos (adjusted for 2 percent inflation)
aguas_long$food <- 0.001*aguas_long$total_food/aguas_long$total_people*(1.02)^2

## split based on sex of main household provider
aguas_F <- aguas_long[aguas_long$sex == 2,]
aguas_M <- aguas_long[aguas_long$sex ==1,]

## remove households with more than 20000 pesos per person per month;
## this is three households with female providers, and 2 with male providers
aguas_F <- subset(aguas_F, aguas_F$food <= 20)
aguas_M <- subset(aguas_M, aguas_M$food <= 20)

## save food expenditure data for both groups
write.csv(aguas_F[c("food")], "aguas2018_food_F.csv", row.names = FALSE)
write.csv(aguas_M[c("food")], "aguas2018_food_M.csv", row.names = FALSE)

## find the median food expenditure per person in upper income household from same region
aguas_soc4 <- aguas_full[aguas_full$est_socio == 4,]

## exclude households with no quarterly food expenditure (none for this example)
aguas_soc4 <- aguas_soc4[aguas_soc4$total_food > 0,]

## calculate food expense per person in thousands of pesos (adjusted for 2 percent inflation)
aguas_soc4$food <- 0.001*aguas_soc4$total_food/aguas_soc4$total_people*(1.02)^2

food4_rep <- NULL
for (i in 1:nrow(aguas_soc4)){
  food4_rep <- c(food4_rep, rep(aguas_soc4$food[i], aguas_soc4$factor[i]))
}

## confirmation that this expense is 4.29 thousand pesos (design value for threshold)
median(food4_rep)

## we now obtain design values for alpha1, beta1, alpha2, and beta2
## first load the .csv files just created
## save the food data in vector form
y1 <- read.csv("aguas2018_food_F.csv")$food 
y2 <- read.csv("aguas2018_food_M.csv")$food

GammaPostSIR <- function(k, l, y1, y2, mu1, tau1, kappa1, lambda1, 
                         mu2, tau2, kappa2, lambda2, metric = "tail", tau = 4.29){
  
  ## k is the number of points generated from sampling distribution
  ## l is the number of resampled points returned by algorithm
  ## y1 and y2 are the food expenditure observations in each group
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
  ## return posterior draws needed to compute figure
  return(list(alpha1, beta1, alpha2, beta2, theta))
}

## approximate the posterior (this may several seconds)
set.seed(2)
design.data <- GammaPostSIR(1000000, 100000, y1, y2, 2, 0.1, 2, 0.1, 
                             2, 0.1, 2, 0.1, metric = "tail", tau = 4.29)

## extract the marginal posteriors from list form
theta <- design.data[[5]]
alpha.1 <- design.data[[1]]
beta.1 <- design.data[[2]]
alpha.2 <- design.data[[3]]
beta.2 <- design.data[[4]]

## find posterior means for design values
mean(alpha.1) ## should be 2.11
mean(beta.1) ## should be 0.69
mean(alpha.2) ## should be 2.43
mean(beta.2) ## should be 0.79

## save posterior draws to .csv files for reference
write.csv(alpha.1, "alpha1s_2018.csv", row.names = FALSE)
write.csv(beta.1, "beta1s_2018.csv", row.names = FALSE)
write.csv(alpha.2, "alpha2s_2018.csv", row.names = FALSE)
write.csv(beta.2, "beta2s_2018.csv", row.names = FALSE)

## now we inflate the variances of these approximately gamma
## marginal posteriors to find informative priors for the numerical study

## define factor by which to inflate the variance
var_inf <- 10

## consider alpha1
## find the posterior mode and variance
a1dens <- density(alpha.1, n = 1024, bw = "SJ")
a1_mode <- a1dens$x[which.max(a1dens$y)]
a1_var <- var(alpha.1)

## solve system of equations to get desired gamma distribution
b1_star <- (-a1_mode - sqrt(a1_mode^2 + 4*var_inf*a1_var))/(-2*var_inf*a1_var); b1_star
a1_star <- 1 + a1_mode*b1_star; a1_star

## should be GAMMA(33.79, 15.66)
informs <- c(a1_star, b1_star)

## consider beta1
## find the posterior mode and variance
a1dens <- density(beta.1, n = 1024, bw = "SJ")
a1_mode <- a1dens$x[which.max(a1dens$y)]
a1_var <- var(beta.1)

## solve system of equations to get desired gamma distribution
b1_star <- (-a1_mode - sqrt(a1_mode^2 + 4*var_inf*a1_var))/(-2*var_inf*a1_var); b1_star
a1_star <- 1 + a1_mode*b1_star; a1_star

## should be GAMMA(26.96, 37.92)
informs <- rbind(informs, c(a1_star, b1_star))

## consider alpha2
## find the posterior mode and variance
a1dens <- density(alpha.2, n = 1024, bw = "SJ")
a1_mode <- a1dens$x[which.max(a1dens$y)]
a1_var <- var(alpha.2)

## solve system of equations to get desired gamma distribution
b1_star <- (-a1_mode - sqrt(a1_mode^2 + 4*var_inf*a1_var))/(-2*var_inf*a1_var); b1_star
a1_star <- 1 + a1_mode*b1_star; a1_star

## should be GAMMA(105.53, 42.97)
informs <- rbind(informs, c(a1_star, b1_star))

## consider beta2
## find the posterior mode and variance
a1dens <- density(beta.2, n = 1024, bw = "SJ")
a1_mode <- a1dens$x[which.max(a1dens$y)]
a1_var <- var(beta.2)

## solve system of equations to get desired gamma distribution
b1_star <- (-a1_mode - sqrt(a1_mode^2 + 4*var_inf*a1_var))/(-2*var_inf*a1_var); b1_star
a1_star <- 1 + a1_mode*b1_star; a1_star

## should be GAMMA(85.43, 106.31)
informs <- rbind(informs, c(a1_star, b1_star))

## output informative gamma distribution hyperparameters to .csv file
write.csv(informs, "informs_aguas.csv", row.names = FALSE)