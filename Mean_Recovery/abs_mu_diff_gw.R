# Parameters

alpha=2;beta=5;
scale_val = 0.05;df=3;
shape=2;location=2;
mu = 0;sigma=1;
mul = 1;
df=3;
shape_p = 4;

beta_mean = alpha / (alpha + beta)
gamma_mean = alpha / beta
normal_mean = mu
logistic_mean = 0
t_mean = 0
laplace_mean = mu
chisq_mean = df
lnorm_mean = mul
weibull_mean = gamma(1 + 1/shape)
pareto_mean = (shape_p * location) / (shape_p - 1)


#----------------------------

G1 <- function(y) pnorm(y)
G2 <- function(y) dnorm(y)

EMemp <- function(x, grid, sigma=1, start=0, max_step=1000, eps2=1e-10, eps1=1e-5){
  
  n <- length(x)
  pn <- (hist(x, breaks=grid)$counts)/n
  
  mu_new <- start
  mu_vec <- mu_new
  delta  <- 10
  
  for(i in 1:max_step){
    
    mu_old <- mu_new
    
    alpha <- (grid-mu_old)/sigma
    temp  <- (diff(G2(alpha))+eps2)/(diff(G1(alpha))+eps2)
    
    mu_new <- mu_old - sigma*sum(pn*temp)
    mu_vec <- c(mu_new, mu_vec)
    
    delta <- abs(mu_new-mu_old)
    
    if(delta<eps1) break
    
  }
  
  return(list(mu_hat = mu_new, mu_vec=mu_vec))
  
}


EMtheo <- function(p0, grid, sigma=1, start=0, max_step=1000, eps2=1e-10, eps1=1e-5){
  
  
  pn <- p0
  
  mu_new <- start
  mu_vec <- mu_new
  delta  <- 10
  
  for(i in 1:max_step){
    
    mu_old <- mu_new
    
    alpha <- (grid-mu_old)/sigma
    temp  <- (diff(G2(alpha))+eps2)/(diff(G1(alpha))+eps2)
    
    mu_new <- mu_old - sigma*sum(pn*temp)
    mu_vec <- c(mu_new, mu_vec)
    
    delta <- abs(mu_new-mu_old)
    
    if(delta<eps1) break
    
  }
  
  return(list(mu_hat = mu_new, mu_vec=mu_vec))
  
}


theo_mu_k_variation = function(cdf0,true_mean,min_value_grid,max_value_grid,max_step = 3000){
  
  mu_theo = c()
  delta_vec = c()
  k_vec = c()
  
  for (i in 1:max_step) {
    
    grid <- seq(min_value_grid,max_value_grid, length.out=i+2)
    p0 <- diff(cdf0(grid))
    interval_centers <- (head(grid, -1) + tail(grid, -1)) / 2
    mid_point <- sum(interval_centers * p0) 
    weighted_variance <- sum(p0 * (interval_centers - mid_point)^2)
    
    mu_theo[i] = EMtheo(p0=p0, grid, start=mid_point,sigma = sqrt(weighted_variance))$mu_hat
    delta_vec[i] = min(diff(grid))
    k_vec[i] = i + 1
    
    if ( delta_vec[i] <= 1e-04) {
      
      break
      
    }
    
  }
  
  return(list("Theoretical_mu" = mu_theo,"Grid_Width" = delta_vec,"Number_of_bins" = k_vec))
  
}

theo_mu_k_check = function(cdf0,true_mean,min_value_grid,max_value_grid,max_step = 3000){
  
  mu_theo = c()
  delta_vec = c()
  k_vec = c()
  check = TRUE
  
  for (i in 1:max_step) {
    
    grid <- seq(min_value_grid,max_value_grid, length.out=i+2)
    p0 <- diff(cdf0(grid))
    interval_centers <- (head(grid, -1) + tail(grid, -1)) / 2
    mid_point <- sum(interval_centers * p0) 
    weighted_variance <- sum(p0 * (interval_centers - mid_point)^2)
    
    mu_theo[i] = EMtheo(p0=p0, grid, start=mid_point,sigma = sqrt(weighted_variance))$mu_hat
    delta_vec[i] = min(diff(grid))
    k_vec[i] = i + 1
    
    if ( delta_vec[i] <= 1e-04) {
      
      break
      
    }
    
  }
  
  if( any(abs(mu_theo - true_mean) > 2*delta_vec) == TRUE){
    check = F
  }
  
  return(check)
  
}

# grDevices::pdf("Theo_mu-delta.pdf",width = 16, height = 14)

par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(5,2),cex.lab=1.5,cex.axis=1.5,
    font.axis=1,cex.main=1.5)

# Beta

cdf_beta = function(x){pbeta(x,alpha,beta)}

beta_2delta = theo_mu_k_variation(cdf0 = cdf_beta,true_mean = beta_mean,min_value_grid = 0,max_value_grid = 1)
theo_mu_k_check(cdf0 = cdf_beta,true_mean = beta_mean,min_value_grid = 0,max_value_grid = 1)

# Provided data
theoretical_mu <- beta_2delta$Theoretical_mu
grid_width <- beta_2delta$Grid_Width

# Round grid width for display
grid_width_rounded <- round(grid_width, 2)

# Calculate the upper and lower bounds for the shaded area
upper_bounds <- beta_mean + 2 * grid_width
lower_bounds <- beta_mean - 2 * grid_width

# Plot setup
plot(grid_width_rounded, theoretical_mu, type = "o", pch = 16, col = "#E69F00", 
     xlab = expression(delta), ylab = expression(widehat(mu)[0]), 
     ylim = c(min(lower_bounds), max(upper_bounds)),xlim = c(0.018,max(grid_width_rounded)),
     cex = 2, lwd = 3)

x_start <- 0
x_end <- max(grid_width_rounded)  # or set a specific endpoint if needed

# Draw the horizontal line from x = 0 to x = x_end
segments(x0 = x_start, y0 = beta_mean, x1 = x_end, y1 = beta_mean, col = "red", lty = 2)

# Shading the area within ±2 grid widths around the true mean
polygon(c(grid_width_rounded,0,rev(grid_width_rounded)), 
        c(upper_bounds,beta_mean,rev(lower_bounds)), 
        col="#E69F00", border = NA,angle = 45,density = 10)

polygon(c(grid_width_rounded,0,rev(grid_width_rounded)), 
        c(upper_bounds,beta_mean, rev(lower_bounds)), 
        , col=rgb(230, 159, 0, alpha=100, maxColorValue=255),
        border=NA,angle = 45)

# Gamma

cdf_gamma = function(x){pgamma(x,alpha,beta)}
MM <- ceiling(qgamma(0.9999, alpha, beta))+5

gamma_2delta = theo_mu_k_variation(cdf0 = cdf_gamma,true_mean = gamma_mean,min_value_grid = 0,max_value_grid = MM)
theo_mu_k_check(cdf0 = cdf_gamma,true_mean = gamma_mean,min_value_grid = 0,max_value_grid = MM)

# Provided data
theoretical_mu <- gamma_2delta$Theoretical_mu
grid_width <- gamma_2delta$Grid_Width

# Round grid width for display
grid_width_rounded <- round(grid_width, 2)

# Calculate the upper and lower bounds for the shaded area
upper_bounds <- gamma_mean + 2 * grid_width
lower_bounds <- gamma_mean - 2 * grid_width

# Plot setup
plot(grid_width_rounded, theoretical_mu, type = "o", pch = 16, col = "#009E73", 
     xlab = expression(delta), ylab = expression(widehat(mu)[0]), 
     ylim = c(min(lower_bounds), max(upper_bounds)),xlim = c(0.14,max(grid_width_rounded)),
     cex = 2, lwd = 3)

x_start <- 0
x_end <- max(grid_width_rounded)  # or set a specific endpoint if needed

# Draw the horizontal line from x = 0 to x = x_end
segments(x0 = x_start, y0 = gamma_mean, x1 = x_end, y1 = gamma_mean, col = "red", lty = 2)

# Shading the area within ±2 grid widths around the true mean
polygon(c(grid_width_rounded,0,rev(grid_width_rounded)), 
        c(upper_bounds,gamma_mean,rev(lower_bounds)), 
        col="#009E73", border = NA,angle = 45,density = 10)

polygon(c(grid_width_rounded,0,rev(grid_width_rounded)), 
        c(upper_bounds,gamma_mean, rev(lower_bounds)), 
        , col=rgb(0, 158, 115, alpha=100, maxColorValue=255),
        border=NA,angle = 45)

# Normal

cdf_normal = function(x){pnorm(x,mu,sigma)}
MM <- ceiling(qnorm(0.9999, mu, sigma))+5

normal_2delta = theo_mu_k_variation(cdf0 = cdf_normal,true_mean = normal_mean,min_value_grid = -MM-2,max_value_grid = MM)
theo_mu_k_check(cdf0 = cdf_normal,true_mean = normal_mean,min_value_grid = -MM-2,max_value_grid = MM)

# Provided data
theoretical_mu <- normal_2delta$Theoretical_mu
grid_width <- normal_2delta$Grid_Width

# Round grid width for display
grid_width_rounded <- round(grid_width, 2)

# Calculate the upper and lower bounds for the shaded area
upper_bounds <- normal_mean + 2 * grid_width
lower_bounds <- normal_mean - 2 * grid_width

# Plot setup
plot(grid_width_rounded, theoretical_mu, type = "o", pch = 16, col = "#56B4E9", 
     xlab = expression(delta), ylab = expression(widehat(mu)[0]), 
     ylim = c(min(lower_bounds), max(upper_bounds)),xlim = c(0.32,max(grid_width_rounded)),
     cex = 2, lwd = 3)

x_start <- 0
x_end <- max(grid_width_rounded)  # or set a specific endpoint if needed

# Draw the horizontal line from x = 0 to x = x_end
segments(x0 = x_start-0.02, y0 = normal_mean, x1 = x_end, y1 = normal_mean, col = "red", lty = 2)

# Shading the area within ±2 grid widths around the true mean
polygon(c(grid_width_rounded,0,rev(grid_width_rounded)), 
        c(upper_bounds,normal_mean,rev(lower_bounds)), 
        col="#56B4E9", border = NA,angle = 45,density = 10)

polygon(c(grid_width_rounded,0,rev(grid_width_rounded)), 
        c(upper_bounds,normal_mean, rev(lower_bounds)), 
        , col=rgb(86, 180, 233, alpha=100, maxColorValue=255),
        border=NA,angle = 45)

# Logistic

cdf_logistic <- function(x) {plogis(x)}

MM <- ceiling(qlogis(0.9999))+5

logistic_2delta = theo_mu_k_variation(cdf0 = cdf_logistic,true_mean = logistic_mean,min_value_grid = -MM-2,max_value_grid = MM)
theo_mu_k_check(cdf0 = cdf_logistic,true_mean = logistic_mean,min_value_grid = -MM-2,max_value_grid = MM)

# Provided data
theoretical_mu <- logistic_2delta$Theoretical_mu
grid_width <- logistic_2delta$Grid_Width

# Round grid width for display
grid_width_rounded <- round(grid_width, 2)

# Calculate the upper and lower bounds for the shaded area
upper_bounds <- logistic_mean + 2 * grid_width
lower_bounds <- logistic_mean - 2 * grid_width

# Plot setup
plot(grid_width_rounded, theoretical_mu, type = "o", pch = 16, col = "#CC79A7", 
     xlab = expression(delta), ylab = expression(widehat(mu)[0]), 
     ylim = c(min(lower_bounds), max(upper_bounds)),xlim = c(0.57,max(grid_width_rounded)),
     cex = 2, lwd = 3)

x_start <- 0
x_end <- max(grid_width_rounded)  # or set a specific endpoint if needed

# Draw the horizontal line from x = 0 to x = x_end
segments(x0 = x_start-0.02, y0 = logistic_mean, x1 = x_end, y1 = logistic_mean, col = "red", lty = 2)

# Shading the area within ±2 grid widths around the true mean
polygon(c(grid_width_rounded,0,rev(grid_width_rounded)), 
        c(upper_bounds,logistic_mean,rev(lower_bounds)), 
        col="#CC79A7", border = NA,angle = 45,density = 10)

polygon(c(grid_width_rounded,0,rev(grid_width_rounded)), 
        c(upper_bounds,logistic_mean, rev(lower_bounds)), 
        , col = rgb(204, 121, 167, alpha=100, maxColorValue=255),
        border=NA,angle = 45)
# t

cdf_t = function(x){pt(x,df)}
MM <- ceiling(qt(0.9999, df))+5

t_2delta = theo_mu_k_variation(cdf0 = cdf_t,true_mean = t_mean,min_value_grid = -MM-2,max_value_grid = MM)
theo_mu_k_check(cdf0 = cdf_t,true_mean = t_mean,min_value_grid = -MM-2,max_value_grid = MM)

# Provided data
theoretical_mu <- t_2delta$Theoretical_mu
grid_width <- t_2delta$Grid_Width

# Round grid width for display
grid_width_rounded <- round(grid_width, 2)

# Calculate the upper and lower bounds for the shaded area
upper_bounds <- t_mean + 2 * grid_width
lower_bounds <- t_mean - 2 * grid_width

# Plot setup
plot(grid_width_rounded, theoretical_mu, type = "o", pch = 16, col = "#D55E00", 
     xlab = expression(delta), ylab = expression(widehat(mu)[0]), 
     ylim = c(min(lower_bounds), max(upper_bounds)),xlim = c(1,max(grid_width_rounded)),
     cex = 2, lwd = 3)

x_start <- 0
x_end <- max(grid_width_rounded)  # or set a specific endpoint if needed

# Draw the horizontal line from x = 0 to x = x_end
segments(x0 = x_start-0.02, y0 = t_mean, x1 = x_end, y1 = t_mean, col = "red", lty = 2)

# Shading the area within ±2 grid widths around the true mean
polygon(c(grid_width_rounded,0,rev(grid_width_rounded)), 
        c(upper_bounds,t_mean,rev(lower_bounds)), 
        col="#D55E00", border = NA,angle = 45,density = 10)

polygon(c(grid_width_rounded,0,rev(grid_width_rounded)), 
        c(upper_bounds,t_mean, rev(lower_bounds)), 
        , col = rgb(213, 94, 0, alpha=100, maxColorValue=255),
        border=NA,angle = 45)


library(jmuOutlier)
library(EnvStats)

# Laplace

cdf_laplace <- function(x) {plaplace(x,mu,sd = sigma)}
MM <- ceiling(qlaplace(0.9999, mu, sigma))+5

laplace_2delta = theo_mu_k_variation(cdf0 = cdf_laplace,true_mean = laplace_mean,min_value_grid = -MM-2,max_value_grid = MM)
theo_mu_k_check(cdf0 = cdf_laplace,true_mean = laplace_mean,min_value_grid = -MM-2,max_value_grid = MM)


# Provided data
theoretical_mu <- laplace_2delta$Theoretical_mu
grid_width <- laplace_2delta$Grid_Width

# Round grid width for display
grid_width_rounded <- round(grid_width, 2)

# Calculate the upper and lower bounds for the shaded area
upper_bounds <- laplace_mean + 2 * grid_width
lower_bounds <- laplace_mean - 2 * grid_width

# Plot setup
plot(grid_width_rounded, theoretical_mu, type = "o", pch = 16, col = "#FDB462", 
     xlab = expression(delta), ylab = expression(widehat(mu)[0]), 
     ylim = c(min(lower_bounds), max(upper_bounds)),xlim = c(0.42,max(grid_width_rounded)),
     cex = 2, lwd = 3)

x_start <- 0
x_end <- max(grid_width_rounded)  # or set a specific endpoint if needed

# Draw the horizontal line from x = 0 to x = x_end
segments(x0 = x_start-0.02, y0 = laplace_mean, x1 = x_end, y1 = laplace_mean, col = "red", lty = 2)

# Shading the area within ±2 grid widths around the true mean
polygon(c(grid_width_rounded,0,rev(grid_width_rounded)), 
        c(upper_bounds,laplace_mean,rev(lower_bounds)), 
        col="#FDB462", border = NA,angle = 45,density = 10)

polygon(c(grid_width_rounded,0,rev(grid_width_rounded)), 
        c(upper_bounds,laplace_mean, rev(lower_bounds)), 
        , col = rgb(253, 180, 98, alpha=100, maxColorValue=255),
        border=NA,angle = 45)


# Chi-Square 

cdf_chisq <- function(x) {pchisq(x,df)}
MM <- ceiling(qchisq(0.9999, df))+5

chisq_2delta = theo_mu_k_variation(cdf0 = cdf_chisq,true_mean = chisq_mean,min_value_grid = -MM-2,max_value_grid = MM)
theo_mu_k_check(cdf0 = cdf_chisq,true_mean = chisq_mean,min_value_grid = -MM-2,max_value_grid = MM)

# Provided data
theoretical_mu <- chisq_2delta$Theoretical_mu
grid_width <- chisq_2delta$Grid_Width

# Round grid width for display
grid_width_rounded <- round(grid_width, 2)

# Calculate the upper and lower bounds for the shaded area
upper_bounds <- chisq_mean + 2 * grid_width
lower_bounds <- chisq_mean - 2 * grid_width

# Plot setup
plot(grid_width_rounded, theoretical_mu, type = "o", pch = 16, col = "#B3DE69", 
     xlab = expression(delta), ylab = expression(widehat(mu)[0]), 
     ylim = c(min(lower_bounds), max(upper_bounds)),xlim = c(1,max(grid_width_rounded)),
     cex = 2, lwd = 3)

x_start <- 0
x_end <- max(grid_width_rounded)  # or set a specific endpoint if needed

# Draw the horizontal line from x = 0 to x = x_end
segments(x0 = x_start-0.02, y0 = chisq_mean, x1 = x_end, y1 = chisq_mean, col = "red", lty = 2)

# Shading the area within ±2 grid widths around the true mean
polygon(c(grid_width_rounded,0,rev(grid_width_rounded)), 
        c(upper_bounds,chisq_mean,rev(lower_bounds)), 
        col="#B3DE69", border = NA,angle = 45,density = 10)

polygon(c(grid_width_rounded,0,rev(grid_width_rounded)), 
        c(upper_bounds,chisq_mean, rev(lower_bounds)), 
        , col = rgb(179, 222, 105, alpha=100, maxColorValue=255),
        border=NA,angle = 45)

# Log-Normal

cdf_lnorm <- function(x) {plnorm(x,meanlog = mul,sdlog = sigma)}
MM <- ceiling(qlnorm(0.9999, mul, sigma))+5

lnorm_2delta = theo_mu_k_variation(cdf0 = cdf_lnorm,true_mean = lnorm_mean,min_value_grid = -MM-2,max_value_grid = MM)
theo_mu_k_check(cdf0 = cdf_lnorm,true_mean = lnorm_mean,min_value_grid = -MM-2,max_value_grid = MM)

# Provided data
theoretical_mu <- lnorm_2delta$Theoretical_mu
grid_width <- lnorm_2delta$Grid_Width

# Round grid width for display
grid_width_rounded <- round(grid_width, 2)

# Calculate the upper and lower bounds for the shaded area
upper_bounds <- lnorm_mean + 2 * grid_width
lower_bounds <- lnorm_mean - 2 * grid_width

# Plot setup
plot(grid_width_rounded, theoretical_mu, type = "o", pch = 16, col = "#BC80BD", 
     xlab = expression(delta), ylab = expression(widehat(mu)[0]), 
     ylim = c(min(lower_bounds), max(upper_bounds)),xlim = c(4,max(grid_width_rounded)),
     cex = 2, lwd = 3)

x_start <- 0
x_end <- max(grid_width_rounded)  # or set a specific endpoint if needed

# Draw the horizontal line from x = 0 to x = x_end
segments(x0 = x_start-0.02, y0 = lnorm_mean, x1 = x_end, y1 = lnorm_mean, col = "red", lty = 2)

# Shading the area within ±2 grid widths around the true mean
polygon(c(grid_width_rounded,0,rev(grid_width_rounded)), 
        c(upper_bounds,lnorm_mean,rev(lower_bounds)), 
        col="#BC80BD", border = NA,angle = 45,density = 10)

polygon(c(grid_width_rounded,0,rev(grid_width_rounded)), 
        c(upper_bounds,lnorm_mean, rev(lower_bounds)), 
        , col = rgb(188, 128, 189, alpha=100, maxColorValue=255),
        border=NA,angle = 45)

# Weibull 

cdf_weibull <- function(x) {pweibull(x,shape)}
MM <- ceiling(qweibull(0.9999, shape))+5

weibull_2delta = theo_mu_k_variation(cdf0 = cdf_weibull,true_mean = weibull_mean,min_value_grid = -MM-2,max_value_grid = MM)
theo_mu_k_check(cdf0 = cdf_weibull,true_mean = weibull_mean,min_value_grid = -MM-2,max_value_grid = MM)

# Provided data
theoretical_mu <- weibull_2delta$Theoretical_mu
grid_width <- weibull_2delta$Grid_Width

# Round grid width for display
grid_width_rounded <- round(grid_width, 2)

# Calculate the upper and lower bounds for the shaded area
upper_bounds <- weibull_mean + 2 * grid_width
lower_bounds <- weibull_mean - 2 * grid_width

# Plot setup
plot(grid_width_rounded, theoretical_mu, type = "o", pch = 16, col = "cyan4", 
     xlab = expression(delta), ylab = expression(widehat(mu)[0]), 
     ylim = c(min(lower_bounds), max(upper_bounds)),xlim = c(0.35,max(grid_width_rounded)),
     cex = 2, lwd = 3)

x_start <- 0
x_end <- max(grid_width_rounded)  # or set a specific endpoint if needed

# Draw the horizontal line from x = 0 to x = x_end
segments(x0 = x_start-0.02, y0 = weibull_mean, x1 = x_end, y1 = weibull_mean, col = "red", lty = 2)

# Shading the area within ±2 grid widths around the true mean
polygon(c(grid_width_rounded,0,rev(grid_width_rounded)), 
        c(upper_bounds,weibull_mean,rev(lower_bounds)), 
        col="cyan4", border = NA,angle = 45,density = 10)

polygon(c(grid_width_rounded,0,rev(grid_width_rounded)), 
        c(upper_bounds,weibull_mean, rev(lower_bounds)), 
        , col=rgb(0, 139, 139, alpha=100, maxColorValue=255),
        border=NA,angle = 45)


# Pareto

cdf_pareto <- function(x) {ppareto(x,location = location,shape = shape_p)}
MM <- ceiling(qpareto(0.9999, location, shape_p))+5

pareto_2delta = theo_mu_k_variation(cdf0 = cdf_pareto,true_mean = pareto_mean,min_value_grid = -MM-2,max_value_grid = MM)
theo_mu_k_check(cdf0 = cdf_pareto,true_mean = pareto_mean,min_value_grid = -MM-2,max_value_grid = MM)

# Provided data
theoretical_mu <- pareto_2delta$Theoretical_mu
grid_width <- pareto_2delta$Grid_Width

# Round grid width for display
grid_width_rounded <- round(grid_width, 2)

# Calculate the upper and lower bounds for the shaded area
upper_bounds <- pareto_mean + 2 * grid_width
lower_bounds <- pareto_mean - 2 * grid_width

# Plot setup
plot(grid_width_rounded, theoretical_mu, type = "o", pch = 16, col = "#FB8072", 
     xlab = expression(delta), ylab = expression(widehat(mu)[0]), 
     ylim = c(min(lower_bounds), max(upper_bounds)),xlim = c(0.95,max(grid_width_rounded)),
     cex = 2, lwd = 3)

x_start <- 0
x_end <- max(grid_width_rounded)  # or set a specific endpoint if needed

# Draw the horizontal line from x = 0 to x = x_end
segments(x0 = x_start-0.02, y0 = pareto_mean, x1 = x_end, y1 = pareto_mean, col = "red", lty = 2)

# Shading the area within ±2 grid widths around the true mean
polygon(c(grid_width_rounded,0,rev(grid_width_rounded)), 
        c(upper_bounds,pareto_mean,rev(lower_bounds)), 
        col="#FB8072", border = NA,angle = 45,density = 10)

polygon(c(grid_width_rounded,0,rev(grid_width_rounded)), 
        c(upper_bounds,pareto_mean, rev(lower_bounds)), 
        ,col = rgb(251, 128, 114, alpha=100, maxColorValue=255),
        border=NA,angle = 45)

# dev.off()




