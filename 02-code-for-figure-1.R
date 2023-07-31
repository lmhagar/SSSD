## code to reproduce Figure 1 in the main text
## this file should be run after "01-food-data-2020.R"

## BEGIN SETUP ##

## load necessary packages
require(ggplot2)
require(cowplot)
require(mvtnorm)
require(MASS)

## load the .csv files created using "01-food-data-2020.R" 
## save the food data in vector form
y1 <- read.csv("aguas_food_F.csv")$food 
y2 <- read.csv("aguas_food_M.csv")$food

## save the data as data frame (for plotting later)
foodF <- read.csv("aguas_food_F.csv")
foodM <- read.csv("aguas_food_M.csv")

GammaPostSIR <- function(k, l, y1, y2, mu1, tau1, kappa1, lambda1, 
         mu2, tau2, kappa2, lambda2, metric = "tail", tau = 4.86){
  
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
  
  ## summarize quantities from data
  n1 <- length(y1)
  n2 <- length(y2)
  sum1 <- sum(y1)
  sum2 <- sum(y2)
  logprod1 <- sum(log(y1))
  logprod2 <- sum(log(y2))
  
  ## use method of moments as starting point to solve for gamma MLEs
  ## MLE is used in the proposal distribution
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
set.seed(1)
figure1.data <- GammaPostSIR(1000000, 100000, y1, y2, 2, 0.1, 2, 0.1, 
                             2, 0.1, 2, 0.1, metric = "tail", tau = 4.82)

## extract marginal posteriors from list form
theta <- figure1.data[[5]]
alpha.1 <- figure1.data[[1]]
beta.1 <- figure1.data[[2]]
alpha.2 <- figure1.data[[3]]
beta.2 <- figure1.data[[4]]

## find good histogram breaks based on the data
breaks1 <- hist(y1, plot=FALSE)$breaks
breaks2 <- hist(y2, plot=FALSE)$breaks
binwidth2 <-0.64
bins2 <- seq(0, 19.2, binwidth2)

## create data frames for plotting
x_fit1 <- seq(floor(qgamma(0.001,mean(alpha.1),mean(beta.1))),
              ceiling(qgamma(0.999,mean(alpha.1),mean(beta.1))), by = 0.01)
fit1 <- data.frame(xfit = x_fit1, yfit = dgamma(x_fit1,mean(alpha.1),mean(beta.1)))

x_fit2 <- seq(floor(qgamma(0.001,mean(alpha.2),mean(beta.2))),
              ceiling(qgamma(0.999,mean(alpha.2),mean(beta.2))), by = 0.01)
fit2 <- data.frame(xfit = x_fit2, yfit = dgamma(x_fit2,mean(alpha.2),mean(beta.2)))

hist_lower2 <- 0
hist_upper2 <- 20

group1_name <- "Female Household Provider"
n1 <- length(y1)

use_pctl <- 1
metric_plot <- paste("Tail Probability", sep="")
metric_value1 <- round(mean(ifelse(foodF$food >= 4.82,1,0)),3)
percentile <- 0.9

## make the plot for the female provider group data (top left corner)
plot1a <-
  ggplot(data=foodF, aes(foodF$food)) + theme_bw() +
  geom_histogram(aes(y = (stat(count) / sum(count))/binwidth2), breaks = bins2,
                 col="cadetblue4", 
                 fill="cadetblue2", 
                 alpha = 0, size = 1) + 
  coord_cartesian(xlim = c(hist_lower2, hist_upper2)) +
  labs(title=paste(group1_name, "\n")) +
  labs(x=bquote(atop(' ', atop(textstyle('Food Expenditure per Person (MXN $1000)'),
                               textstyle('n = '*.(n1)*", "*.(metric_plot)*' ('*hat(theta)[1]*  
                                           ') = '*.(metric_value1))))), y="Density\n") + 
  theme(plot.title = element_text(hjust = 0.5,size=16,face="bold")) + 
  theme(plot.subtitle = element_text(hjust = 0.5,size=14)) +
  geom_area(data = fit1, aes(x=xfit, y=yfit), fill="cadetblue2", col="cadetblue4", alpha=0.45, size = 1) +
  geom_segment(aes(x = 4.82, y = 0, 
                   xend=4.82, yend = (Inf*(use_pctl))), color="grey16", size=1.5) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) + ylim(0, 0.3)

group2_name <- "Male Household Provider"
n2 <- length(y2)
metric_value2 <- round(mean(ifelse(foodM$food >= 4.82,1,0)),3)

## make the plot for the male provider group (bottom left corner)
plot1c <- ggplot(data=foodM, aes(foodM$food)) + theme_bw() +
  geom_histogram(aes(y = (stat(count) / sum(count))/binwidth2), breaks = bins2,
                 col="darkolivegreen", 
                 fill="darkolivegreen2", 
                 alpha = 0, size = 1) + 
  coord_cartesian(xlim = c(hist_lower2, hist_upper2)) +
  labs(title=paste(group2_name,"\n")) +
  labs(x=bquote(atop(' ', atop(textstyle('Food Expenditure per Person (MXN $1000)'),
                               textstyle('n = '*.(n2)*", "*.(metric_plot)*' ('*hat(theta)[2]*  
                                           ') = '*.(metric_value2))))), y="Density\n") + 
  theme(plot.title = element_text(hjust = 0.5,size=16,face="bold")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  geom_area(data = fit2, aes(x=xfit, y=yfit), fill="darkolivegreen2", col="darkolivegreen", alpha=0.45, size = 1) + 
  geom_segment(aes(x = 4.82, y = 0, 
                   xend=4.82, yend = (Inf*(use_pctl))), color="grey16", size=1.5) + ylim(0, 0.3)

## approximate the posterior using nonparametric density estimation
d <- density(theta)
d_frame <- data.frame(x = d$x, y = d$y)

## make plot with no red shading first (additional layers to be added later)
plot1e <- ggplot(data=d_frame, aes(x=d_frame$x)) + theme_bw() +
  geom_polygon(aes(y=d_frame$y), col="gray10", fill="snow1", size=0.75, alpha=0.9) +
  labs(title = bquote(bold("Posterior pdf of "~ theta[1] ~ '\u00F7' ~ theta[2]))) +
  labs(subtitle = " ") +
  labs(x=bquote(bold('\n Posterior pdf of'~ theta[1] ~ '\u00F7' ~ theta[2])), y=" ") + 
  theme(plot.title = element_text(hjust = 0.5,size=18,face="bold", colour = "#FFFFFF")) +
  theme(plot.subtitle = element_text(hjust = 0.5,size=16,face="bold")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold")) + 
  theme(axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")))

## calculate CPMs
delta.1 <- 1/1.1
delta.2 <- 1.1

CPM_4 <- round(mean(ifelse(theta < delta.2, theta > delta.1,0)),4)
CPM_title <- bquote(bold(italic(Pr)*'( '* 1.1^-1 < theta[1]~'\u00F7'~ theta[2] < 1.1*' | '*italic(data)*') =' ~ .(CPM_4)))

x_title <- bquote(bold('\n Posterior pdf of'~ theta[1] ~ '\u00F7' ~ theta[2]))

CPM_subset = subset(d_frame, x > delta.1 & x < delta.2)
x_new = 
  c(max(delta.1,min(CPM_subset$x)), unlist(CPM_subset$x), 
    min(delta.2,max(CPM_subset$x)))
y_new = c(0,unlist(CPM_subset$y),0)
CPM_frame = data.frame(x = x_new,y = y_new)

## make plot for probability of agreement (top right corner)
plot1b <- 
  plot1e + labs(title = NULL) + labs(x = NULL) +
  labs(title=CPM_title) +
  labs(x= x_title) +
  geom_polygon(data = CPM_frame, aes(x=x, y=y), fill="red", col="gray10", size = 0.75, alpha = 0.75) +
  theme(plot.title = element_text(hjust = 0.5,size=18,face="bold", colour = "#000000"))

delta.1 <- 1/1.1
delta.2 <- Inf

CPM_4 <- round(mean(theta > delta.1),4)
CPM_title <- bquote(bold(italic(Pr)*'( '* theta[1] > 1.1^-1*theta[2]*' | '*italic(data)*') =' ~ .(CPM_4)))

x_title <- bquote(bold('\n Posterior pdf of'~ theta[1] ~ '\u00F7' ~ theta[2]))

CPM_subset = subset(d_frame, x > delta.1 & x < delta.2)
x_new = 
  c(max(delta.1,min(CPM_subset$x)), unlist(CPM_subset$x), 
    min(delta.2,max(CPM_subset$x)))
y_new = c(0,unlist(CPM_subset$y),0)
CPM_frame = data.frame(x = x_new,y = y_new)

## make plot for probability of noninferiority (bottom right corner)
plot1d <- 
  plot1e + labs(title = NULL) + labs(x = NULL) +
  labs(title=CPM_title) +
  labs(x= x_title) +
  geom_polygon(data = CPM_frame, aes(x=x, y=y), fill="red", col="gray10", size = 0.75, alpha = 0.75) +
  theme(plot.title = element_text(hjust = 0.5,size=18,face="bold", colour = "#000000"))

## combine into grid
fig1.row1 <- plot_grid(plot1a + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")), 
                       plot1b + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")),
                       rel_widths = c(1.45, 1.15))
fig1.row2 <- plot_grid(plot1c + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")), 
                       plot1d + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")),
                       rel_widths = c(1.45, 1.15))
fig1 <- plot_grid(fig1.row1, fig1.row2, nrow = 2)

## output as .pdf file for the article
pdf(file = "Figure1BA.pdf",   # The directory you want to save the file in
    width = 11.7895, # The width of the plot in inches (12.41)
    height = 10.165) # The height of the plot in inches (10.7)

fig1

dev.off()