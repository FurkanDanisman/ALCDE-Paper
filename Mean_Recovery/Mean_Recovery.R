G1 <- function(y) pnorm(y)
G2 <- function(y) dnorm(y)
G3 <- function(y) dnorm(y)*y 

EMemp_mu <- function(x, grid, sigma, start, max_step=3000, eps2=1e-10, eps1=1e-5){
  
  n <- length(x)
  pn <- (hist(x, breaks=grid,plot = F)$counts)/n
  
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


EMemp_var <- function(x, grid, mu_new, mu_old,sigma, max_step=3000, eps2=1e-10, eps1=1e-5){
  
  n <- length(x)
  pn <- (hist(x, breaks=grid,plot = F)$counts)/n
  
  alpha <- (grid-mu_old)/sigma
  
  temp1 = (diff(G3(alpha))+eps2)/(diff(G1(alpha))+eps2)
  temp2 = (diff(G2(alpha))+eps2)/(diff(G1(alpha))+eps2)
  
  var_new = sigma^2*sum((1-temp1)*pn) + (mu_new - mu_old)^2 + 2*sigma*(mu_new - mu_old)*sum(temp2*pn)
  
  return(var_new)
  
}

EMtheo_mu <- function(p0, grid, sigma, start, max_step=3000, eps2=1e-10, eps1=1e-5){
  
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


EMtheo_var <- function(p0, grid, mu_new, mu_old, sigma){
  
  pn <- p0
  
  alpha <- (grid-mu_old)/sigma
  
  temp1 = (diff(G3(alpha))+eps2)/(diff(G1(alpha))+eps2)
  temp2 = (diff(G2(alpha))+eps2)/(diff(G1(alpha))+eps2)
  
  var_new = sigma^2*sum((1-temp1)*pn) + (mu_new - mu_old)^2 + 2*sigma*(mu_new - mu_old)*sum(temp2*pn)
  
  return(var_new)
  
}


EM_emp <- function(x,grid,mu_0,var_0,max_step=5000,eps2=1e-10, eps1=1e-5){
  
  flag<- 0
  mu_new = mu_0
  var_new = var_0
  
  for (i in 1:max_step){
    
    mu_old = mu_new
    var_old = var_new
    mu_new <- EMemp_mu(x, grid, start=mu_old,sigma = sqrt(var_old))$mu_hat
    var_new <- EMemp_var(x, grid, sigma=sqrt(var_old),mu_new = mu_new,mu_old = mu_old)
    
    if(abs(mu_old-mu_new)<eps1 &
       abs(var_old-var_new)<eps1){
      flag<-1 ;break
    }
    
  }
  
  if(!flag){
    warning("Didn't Converge \n")
  }
  
  updateres<- base::list("mu_estimate" = mu_new,
                         "var_estimate" = var_new)
  
  return(updateres)
}

EM_theo <- function(p0,grid,mu_0,var_0,max_step=5000,eps2=1e-10, eps1=1e-5){
  
  flag<- 0
  mu_new = mu_0
  var_new = var_0
  
  for (i in 1:max_step){
    
    mu_old = mu_new
    var_old = var_new
    mu_new <- EMtheo_mu(p0, grid, start=mu_old,sigma = sqrt(var_old))$mu_hat
    var_new <- EMtheo_var(p0, grid, sigma=sqrt(var_old),mu_new = mu_new,mu_old = mu_old)
    
    if(abs(mu_old-mu_new)<eps1 &
       abs(var_old-var_new)<eps1){
      flag<-1 ;break
    }
    
  }
  
  if(!flag){
    warning("Didn't Converge \n")
  }
  
  updateres<- base::list("mu_estimate" = mu_new,
                         "var_estimate" = var_new)
  
  return(updateres)
}

n <- 1000
grid <- seq(-5,5, length.out=21)
x <- rnorm(n,0,1)
p0 <- diff(pnorm(grid,0,1))
pn <- (hist(x, breaks=grid,plot = F)$counts)/n

# Empirical Results 

EMemp_mu(x, grid, start=2,sigma=200)$mu_hat
EM_emp(x, grid,mu_0 = 2,var_0 = 2.4)

# Theoretical Results

EMtheo_mu(p0=p0, grid, start=2,sigma=2)$mu_hat
EM_theo(p0=p0, grid,mu_0 = 2,var_0 = 2)

n <- 1000
grid <- seq(0,1, length.out=10)
a1 <- 2
a2 <- 5
x <- rbeta(n,a1,a2)
p0 <- diff(pbeta(grid,alpha,beta))
pn <- (hist(x, breaks=grid)$counts)/n

# Empirical Results 

EMemp_mu(x, grid, start=2,sigma=2)$mu_hat
EM_emp(x, grid,mu_0 = 2,var_0 = 2)

beta_sd = sqrt(alpha*beta / ((alpha + beta)^2 * (alpha + beta +1 ) ) )
beta_var = beta_sd^2

# Theoretical Results

EMtheo_mu(p0=p0, grid, start=2,sigma=2)$mu_hat
EM_theo(p0=p0, grid,mu_0 = 2,var_0 = 2)


# Plot-1 

# Function 

k_bound_calculator = function(cdf0,true_mean,min_value_grid,max_value_grid,max_step = 3000){
  
  k = 0 
  
  for (i in 1:max_step) {
    
    grid <- seq(min_value_grid,max_value_grid, length.out=i+2)
    p0 <- diff(cdf0(grid))
    interval_centers <- (head(grid, -1) + tail(grid, -1)) / 2
    mid_point <- sum(interval_centers * p0)  # Weighted average
    weighted_variance <- sum(p0 * (interval_centers - mid_point)^2)
    
    mu_theo = EMtheo(p0=p0, grid,start = mid_point,sigma = sqrt(weighted_variance))$mu_hat
    
    if (abs(mu_theo - true_mean) < 1e-4) {
      
      k = (i+1)
      
      break
      
    }
    
  }
  
  return(k)
  
}

k_bound_calculator_ver2 = function(cdf0,true_mean,true_sigma,min_value_grid,max_value_grid,max_step = 40){
  
  delta_vec = c()
  mu_theo = c()
  delta_over_sd = c()
  k = c()
  
  for (i in 1:max_step) {
    
    delta_over_sd[i] = i/4
    
    delta_vec[i] = (i/4) * true_sigma
    
    if (delta_vec[i] > max_value_grid - min_value_grid) {
      mu_theo[i] = NA
      next
    }
    grid <- seq(min_value_grid,max_value_grid, by = delta_vec[i])
    
    p0 <- diff(cdf0(grid))
    
    interval_centers <- (head(grid, -1) + tail(grid, -1)) / 2
    mid_point <- sum(interval_centers * p0)  # Weighted average
    weighted_variance <- sum(p0 * (interval_centers - mid_point)^2)
    
    mu_theo[i] = EMtheo(p0=p0, grid,start = mid_point,sigma = sqrt(weighted_variance))$mu_hat
    
    k[i] = length(diff(grid))
    
  }
  
  return(list("mu_diff" = abs(true_mean-mu_theo),"Delta" = delta_vec,"delta_over_sd" = delta_over_sd,"k" = k))
  
}

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

# Standard deviations
beta_sd <- sqrt((alpha * beta) / ((alpha + beta)^2 * (alpha + beta + 1)))
gamma_sd <- sqrt(alpha) / beta
normal_sd <- sigma
logistic_sd <- sqrt((pi^2) / 3) * 1
t_sd <- ifelse(df > 2, sqrt(df / (df - 2)), Inf)
laplace_sd <- sqrt(2) * sigma
chisq_sd <- sqrt(2 * df)
lnorm_sd <- sqrt((exp(sigma^2) - 1) * exp(2 * mul + sigma^2))
weibull_sd <- sqrt(gamma(1 + 2 / shape) - (gamma(1 + 1 / shape)^2))
pareto_sd <- ifelse(shape_p > 2, sqrt((shape_p * location^2) / ((shape_p - 1)^2 * (shape_p - 2))), Inf) # Valid only for shape_p > 2

k_bound_vector = rep(0,10)

# Beta

cdf_beta = function(x){pbeta(x,alpha,beta)}

k_bound_vector[1] = k_bound_calculator(cdf0 = cdf_beta,true_mean = beta_mean,min_value_grid = 0,max_value_grid = 1)

grid <- seq(0,1, length.out=10)
beta_gw = round(min(diff(grid)),2)

Beta_mu_diff_list = k_bound_calculator_ver2(cdf0 = cdf_beta,true_mean = beta_mean,true_sigma = beta_sd,min_value_grid = 0,max_value_grid = 1)

round(Beta_mu_diff_list$mu_diff[which(Beta_mu_diff_list$delta_over_sd %in% c(0.5,1,1.5,2,2.5,3,3.5,4,4.5,5))],4)

round(Beta_mu_diff_list$Delta[which(Beta_mu_diff_list$delta_over_sd %in% c(0.5,1,1.5,2,2.5))],4)

round(Beta_mu_diff_list$k[which(Beta_mu_diff_list$delta_over_sd %in% c(0.5,1,1.5,2,2.5))],4)

# Gamma

cdf_gamma = function(x){pgamma(x,alpha,beta)}
MM <- ceiling(qgamma(0.9999, alpha, beta))+5

k_bound_vector[2] = k_bound_calculator(cdf0 = cdf_gamma,true_mean = gamma_mean,min_value_grid = 0,max_value_grid = MM)

grid <- seq(-MM-2,MM, length.out=46)
gamma_gw = round(min(diff(grid)),2)

Gamma_mu_diff_list = k_bound_calculator_ver2(cdf0 = cdf_gamma,true_mean = gamma_mean,true_sigma = gamma_sd,min_value_grid = 0,max_value_grid = MM)

round(Gamma_mu_diff_list$mu_diff[which(Gamma_mu_diff_list$delta_over_sd %in% c(0.5,1,1.5,2,2.5,3,3.5,4,4.5,5))],4)

round(Gamma_mu_diff_list$Delta[which(Gamma_mu_diff_list$delta_over_sd %in% c(0.5,1,1.5,2,2.5))],4)

round(Gamma_mu_diff_list$k[which(Gamma_mu_diff_list$delta_over_sd %in% c(0.5,1,1.5,2,2.5))],4)

# Normal

cdf_normal = function(x){pnorm(x,mu,sigma)}
MM <- ceiling(qnorm(0.9999, mu, sigma))+5

k_bound_vector[3] = k_bound_calculator(cdf0 = cdf_normal,true_mean = normal_mean,min_value_grid = -MM-2,max_value_grid = MM)

grid <- seq(-MM-2,MM, length.out=11)
normal_gw = round(min(diff(grid)),2)

Normal_mu_diff_list = k_bound_calculator_ver2(cdf0 = cdf_normal,true_mean = normal_mean,true_sigma = normal_sd,min_value_grid = -MM-2,max_value_grid = MM)

round(Normal_mu_diff_list$mu_diff[which(Normal_mu_diff_list$delta_over_sd %in% c(0.5,1,1.5,2,2.5,3,3.5,4,4.5,5))],4)

round(Normal_mu_diff_list$Delta[which(Normal_mu_diff_list$delta_over_sd %in% c(0.5,1,1.5,2,2.5))],4)

round(Normal_mu_diff_list$k[which(Normal_mu_diff_list$delta_over_sd %in% c(0.5,1,1.5,2,2.5))],4)

# Logistic

cdf_logistic <- function(x) {plogis(x)}

MM <- ceiling(qlogis(0.9999))+5

k_bound_vector[4] = k_bound_calculator(cdf0 = cdf_logistic,true_mean = logistic_mean,min_value_grid = -MM-2,max_value_grid = MM)

grid <- seq(-MM-2,MM, length.out=17)
logistic_gw = round(min(diff(grid)),2)

Logistic_mu_diff_list = k_bound_calculator_ver2(cdf0 = cdf_logistic,true_mean = logistic_mean,true_sigma = logistic_sd,min_value_grid = -MM-2,max_value_grid = MM)

round(Logistic_mu_diff_list$mu_diff[which(Logistic_mu_diff_list$delta_over_sd %in% c(0.5,1,1.5,2,2.5,3,3.5,4,4.5,5))],4)

round(Logistic_mu_diff_list$Delta[which(Logistic_mu_diff_list$delta_over_sd %in% c(0.5,1,1.5,2,2.5))],4)

round(Logistic_mu_diff_list$k[which(Logistic_mu_diff_list$delta_over_sd %in% c(0.5,1,1.5,2,2.5))],4)

# t - No Convergence

cdf_t = function(x){pt(x,df)}
MM <- ceiling(qt(0.9999, df))+5

k_bound_vector[5] = k_bound_calculator(cdf0 = cdf_t,true_mean = t_mean,min_value_grid = -MM-2,max_value_grid = MM)

grid <- seq(-MM-2,MM, length.out=18)
t_gw = round(min(diff(grid)),2)

t_mu_diff_list = k_bound_calculator_ver2(cdf0 = cdf_t,true_mean = t_mean,true_sigma = t_sd,min_value_grid = -MM-2,max_value_grid = MM)

round(t_mu_diff_list$mu_diff[which(t_mu_diff_list$delta_over_sd %in% c(0.5,1,1.5,2,2.5,3,3.5,4,4.5,5))],4)

round(t_mu_diff_list$Delta[which(t_mu_diff_list$delta_over_sd %in% c(0.5,1,1.5,2,2.5))],4)

round(t_mu_diff_list$k[which(t_mu_diff_list$delta_over_sd %in% c(0.5,1,1.5,2,2.5))],4)

library(jmuOutlier)
library(EnvStats)

# Laplace

cdf_laplace <- function(x) {plaplace(x,mu,sd = sigma)}
MM <- ceiling(qlaplace(0.9999, mu, sigma))+5

k_bound_vector[6] = k_bound_calculator(cdf0 = cdf_laplace,true_mean = laplace_mean,min_value_grid = -MM-2,max_value_grid = MM)

grid <- seq(-MM-2,MM, length.out=14)
laplace_gw = round(min(diff(grid)),2)

Laplace_mu_diff_list = k_bound_calculator_ver2(cdf0 = cdf_laplace,true_mean = laplace_mean,true_sigma = laplace_sd,min_value_grid = -MM-2,max_value_grid = MM)

round(Laplace_mu_diff_list$mu_diff[which(Laplace_mu_diff_list$delta_over_sd %in% c(0.5,1,1.5,2,2.5,3,3.5,4,4.5,5))],4)

round(Laplace_mu_diff_list$Delta[which(Laplace_mu_diff_list$delta_over_sd %in% c(0.5,1,1.5,2,2.5))],4)

round(Laplace_mu_diff_list$k[which(Laplace_mu_diff_list$delta_over_sd %in% c(0.5,1,1.5,2,2.5))],4)

# Chi-Square - No Convergence

cdf_chisq <- function(x) {pchisq(x,df)}
MM <- ceiling(qchisq(0.9999, df))+5

k_bound_vector[7] = k_bound_calculator(cdf0 = cdf_chisq,true_mean = chisq_mean,min_value_grid = 0,max_value_grid = MM)

grid <- seq(0,MM, length.out=38)
chi_sqr_gw = round(min(diff(grid)),2)

Chisqr_mu_diff_list = k_bound_calculator_ver2(cdf0 = cdf_chisq,true_mean = chisq_mean,true_sigma=chisq_sd,min_value_grid = 0,max_value_grid = MM)

round(Chisqr_mu_diff_list$mu_diff[which(Chisqr_mu_diff_list$delta_over_sd %in% c(0.5,1,1.5,2,2.5,3,3.5,4,4.5,5))],4)

round(Chisqr_mu_diff_list$Delta[which(Chisqr_mu_diff_list$delta_over_sd %in% c(0.5,1,1.5,2,2.5))],4)

round(Chisqr_mu_diff_list$k[which(Chisqr_mu_diff_list$delta_over_sd %in% c(0.5,1,1.5,2,2.5))],4)

# Log-Normal - No Convergence

cdf_lnorm <- function(x) {plnorm(x,meanlog = mul,sdlog = sigma)}
MM <- ceiling(qlnorm(0.9999, mul, sigma))+5

k_bound_vector[8] = k_bound_calculator(cdf0 = cdf_lnorm,true_mean = exp(mul + (sigma^2 / 2)), min_value_grid = 0,max_value_grid = MM)
grid <- seq(0,MM, length.out=38)
lnorm_gw = round(min(diff(grid)),2)

Lnorm_mu_diff_list = k_bound_calculator_ver2(cdf0 = cdf_lnorm,true_mean = exp(mul + (sigma^2 / 2)),true_sigma = lnorm_sd, min_value_grid = 0,max_value_grid = MM)

round(Lnorm_mu_diff_list$mu_diff[which(Lnorm_mu_diff_list$delta_over_sd %in% c(0.5,1,1.5,2,2.5,3,3.5,4,4.5,5))],4)

round(Lnorm_mu_diff_list$Delta[which(Lnorm_mu_diff_list$delta_over_sd %in% c(0.5,1,1.5,2,2.5))],4)

round(Lnorm_mu_diff_list$k[which(Lnorm_mu_diff_list$delta_over_sd %in% c(0.5,1,1.5,2,2.5))],4)

# Weibull 

cdf_weibull <- function(x) {pweibull(x,shape)}
MM <- ceiling(qweibull(0.9999, shape))+5

k_bound_vector[9] = k_bound_calculator(cdf0 = cdf_weibull,true_mean = weibull_mean,min_value_grid = 0,max_value_grid = MM)

Weibull_mu_diff_list = k_bound_calculator_ver2(cdf0 = cdf_weibull,true_mean = weibull_mean,true_sigma = weibull_sd,min_value_grid = 0,max_value_grid = MM)

round(Weibull_mu_diff_list$mu_diff[which(Weibull_mu_diff_list$delta_over_sd %in% c(0.5,1,1.5,2,2.5,3,3.5,4,4.5,5))],4)

round(Weibull_mu_diff_list$Delta[which(Weibull_mu_diff_list$delta_over_sd %in% c(0.5,1,1.5,2,2.5))],4)

round(Weibull_mu_diff_list$k[which(Weibull_mu_diff_list$delta_over_sd %in% c(0.5,1,1.5,2,2.5))],4)

# Pareto - No Convergence

cdf_pareto <- function(x) {ppareto(x,location = location,shape = shape_p)}
MM <- ceiling(qpareto(0.9999, location, shape_p))+5

k_bound_vector[10] = k_bound_calculator(cdf0 = cdf_pareto,true_mean = pareto_mean,min_value_grid = location,max_value_grid = MM)
grid <- seq(location,MM, length.out=38)
pareto_gw = round(min(diff(grid)),2)

Pareto_mu_diff_list = k_bound_calculator_ver2(cdf0 = cdf_pareto,true_mean = pareto_mean,true_sigma = pareto_sd,min_value_grid = location,max_value_grid = MM)

round(Pareto_mu_diff_list$mu_diff[which(Pareto_mu_diff_list$delta_over_sd %in% c(0.5,1,1.5,2,2.5,3,3.5,4,4.5,5))],4)

round(Pareto_mu_diff_list$Delta[which(Pareto_mu_diff_list$delta_over_sd %in% c(0.5,1,1.5,2,2.5))],4)

round(Pareto_mu_diff_list$k[which(Pareto_mu_diff_list$delta_over_sd %in% c(0.5,1,1.5,2,2.5))],4)

n = 10000
x_norm = rnorm(n,mu,sigma)
x_logis = rlogis(n)
x_laplace = rlaplace(x,mu,sd = sigma)
x_t = rt(n,df)
x_gamma = rgamma(n,alpha,beta)
x_beta = rbeta(n,alpha,beta)
x_weibull = rweibull(n,shape)

# Plot 

library(ggplot2)

dist_name = c("Beta","Normal","Laplace","Logistic","Gamma","Weibull")
df_k_bound = data.frame(Distribution = dist_name, k = sort(k_bound_vector[k_bound_vector!=0]))

# Ensure the order in the plot matches the data frame order
df_k_bound$Distribution <- factor(df_k_bound$Distribution, levels = df_k_bound$Distribution)

k_bound_plot = ggplot(df_k_bound, aes(x = Distribution, y = k)) +
  geom_point(size = 4) +  # Enlarged and colored points
  geom_text(aes(label = k), vjust = -1, size = 6, show.legend = FALSE) +  # Add y-values as text and remove from legend
  theme_light() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 22),  # Center and bold title
    axis.title.x = element_blank(),  # Increase x-axis label size
    axis.title.y = element_text(size = 20),  # Increase y-axis label size
    axis.text = element_text(size = 20),     # Increase axis text size
    axis.text.x = element_text(size = 20),   # Increase x-axis text size
    axis.text.y = element_text(size = 20),   # Increase y-axis text size
  ) +
  scale_y_continuous(breaks = seq(0, max(k_bound) + 5, by = 5)) +  # Set Y-axis breaks every 5
  guides(color = guide_legend(override.aes = list(label = ""))) +  # Ensure no labels (like "a") appear in the legend
  ggplot2::labs(
    title = "Minimum Number of Bins Required",   # Title of the plot
    y = "Number of Bins"  # Y-axis label
  )

# Display the plot
k_bound_plot

c(beta_gw / beta_sd,normal_gw / normal_sd,laplace_gw / laplace_sd,logistic_gw / logistic_sd,
  gamma_gw / gamma_sd,Weibull_gw / weibull_sd)



# idea (gw / sd range should be 0.5,1,1.5,...,5) Maybe end it around 2.5 as well but start should be 0.5 I think? 


# Plot-1 

# Function 

n_bound_calculator = function(n,rg,p0,grid,max_step = 3000){
  
  for (i in 1:max_step) {
    
    x = rg(n)
    
    interval_centers <- (head(grid, -1) + tail(grid, -1)) / 2
    mid_point <- sum(interval_centers * p0)  # Weighted average
    weighted_variance <- sum(p0 * (interval_centers - mid_point)^2)
    
    
    mu_emp = EM_emp(x, grid,mu_0 = mid_point,var_0 = weighted_variance)$mu_estimate
    
    
    mu_theo = EM_theo(p0=p0, grid,mu_0 = mid_point,var_0 = weighted_variance)$mu_estimate
    
    if (abs(mu_emp - mu_theo) < 1e-4) {
      break
    }
    
    n = n + 100
    
  }
  
  return(n)
  
}

# Normal

B = 20
n_bound_norm = rep(0,B)

norm_rg = function(n) rnorm(n)

for (i in 1:B) {
  
  MM <- ceiling(qnorm(0.9999))+5
  grid <- seq(-MM,MM, length.out=i + 3)
  p0 <- diff(pnorm(grid))
  n_bound_norm[i] = n_bound_calculator(n = 1000,rg = norm_rg,p0 = p0,grid = grid)
  
}


# Beta

B = 20
n_bound_beta = rep(0,B)
beta_rg = function(n) rbeta(n,alpha,beta)


for (i in 1:B) {
  
  grid <- seq(0,1, length.out=i+3)
  p0 <- diff(pbeta(grid,alpha,beta))
  n_bound_beta[i] = n_bound_calculator(n = 1000,rg = beta_rg,p0 = p0,grid = grid)
  
}


# Gamma

B = 20
n_bound_gamma = rep(0,B)

gamma_rg = function(n) rgamma(n,alpha,beta)

for (i in 1:B) {
  
  MM <- ceiling(qgamma(0.9999,alpha,beta))+5
  grid <- seq(0,MM, length.out=i+3)
  p0 <- diff(pgamma(grid,alpha,beta))
  n_bound_gamma[i] = n_bound_calculator(n = 1000,rg = gamma_rg,p0 = p0,grid = grid)
  
}


# Plot 

library(ggplot2)

# Sample data (replace k_bound_vector with your actual data)
k_num = seq(3,22)
df_n_bound = data.frame(n = n_bound, k = k_num)

# Plot
n_bound_plot_norm = ggplot(df_n_bound, aes(x = k_num, y = n)) +
  geom_point(size = 4, color = "blue") +  # Enlarged and colored points
  theme_light() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 22),  # Center and bold title
    axis.title.x = element_blank(),  # Increase x-axis label size
    axis.title.y = element_text(size = 20),  # Increase y-axis label size
    axis.text = element_text(size = 20),     # Increase axis text size
    axis.text.x = element_text(size = 20),   # Increase x-axis text size
    axis.text.y = element_text(size = 20),   # Increase y-axis text size
  ) +
  ggplot2::labs(
    title = "Normal",   # Title of the plot
    y = "Min(n)"  # Y-axis label
  ) + 
  ylim(0,40000)

# Display the plot
n_bound_plot_norm


# Sample data (replace k_bound_vector with your actual data)
k_num = seq(3,22)
df_n_bound = data.frame(n = n_bound_beta, k = k_num)

# Plot
n_bound_plot_beta = ggplot(df_n_bound, aes(x = k_num, y = n)) +
  geom_point(size = 4, color = "blue") +  # Enlarged and colored points
  theme_light() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 22),  # Center and bold title
    axis.title.x = element_blank(),  # Increase x-axis label size
    axis.title.y = element_text(size = 20),  # Increase y-axis label size
    axis.text = element_text(size = 20),     # Increase axis text size
    axis.text.x = element_text(size = 20),   # Increase x-axis text size
    axis.text.y = element_text(size = 20),   # Increase y-axis text size
  ) +
  ggplot2::labs(
    title = "Beta",   # Title of the plot
    y = "Min(n)"  # Y-axis label
  ) + 
  ylim(0,40000)

# Display the plot
n_bound_plot_beta


# Sample data (replace k_bound_vector with your actual data)
k_num = seq(3,22)
df_n_bound = data.frame(n = n_bound_gamma, k = k_num)

# Plot
n_bound_plot_gamma = ggplot(df_n_bound, aes(x = k_num, y = n)) +
  geom_point(size = 4, color = "blue") +  # Enlarged and colored points
  theme_light() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 22),  # Center and bold title
    axis.title.x = element_blank(),  # Increase x-axis label size
    axis.title.y = element_text(size = 20),  # Increase y-axis label size
    axis.text = element_text(size = 20),     # Increase axis text size
    axis.text.x = element_text(size = 20),   # Increase x-axis text size
    axis.text.y = element_text(size = 20),   # Increase y-axis text size
  ) +
  ggplot2::labs(
    title = "Gamma",   # Title of the plot
    y = "Min(n)"  # Y-axis label
  ) + 
  ylim(0,40000)

# Display the plot
n_bound_plot_gamma


# Plots for non-convergence

NoConvergence = function(cdf0,min_value_grid,max_value_grid,mu_0,var_0,max_step = 200){
  
  mu_theo_vec = c()
  k_vec = c()
  
  index = 1
  
  for (i in 1:max_step) {
    
    if (i %% 4 == 0) {
      
      grid <- seq(min_value_grid,max_value_grid, length.out=i+2)
      p0 <- diff(cdf0(grid))
      interval_centers <- (head(grid, -1) + tail(grid, -1)) / 2
      mid_point <- sum(interval_centers * p0)  # Weighted average
      weighted_variance <- sum(p0 * (interval_centers - mid_point)^2)
      
      mu_theo_vec[index] = EM_theo(p0=p0, grid,mu_0 = mid_point,var_0 = weighted_variance)$mu_estimate
      k_vec[index] = i + 1 
      
      index = index + 1
      
    }
    
  }
  
  return(list("mu_theo_vec" = mu_theo_vec, "k_vec" = k_vec))
  
}

MM <- ceiling(qchisq(0.9999, df))+5
chisq_theo_mu_k = NoConvergence(cdf0 = cdf_chisq,min_value_grid = -MM,max_value_grid = MM,mu_0 = 2,var_0 = 2)
chisq_theo_mu_k = as.data.frame(chisq_theo_mu_k)

MM <- ceiling(qlnorm(0.9999, mul, sigma))+5
lnorm_theo_mu_k = NoConvergence(cdf0 = cdf_lnorm, min_value_grid = -MM,max_value_grid = MM,mu_0 = 2,var_0 = 2)
lnorm_mean = exp(mul + (sigma^2 / 2))
lnorm_theo_mu_k = as.data.frame(lnorm_theo_mu_k)

MM <- ceiling(qpareto(0.9999, location, shape_p))+5
pareto_theo_mu_k = NoConvergence(cdf0 = cdf_pareto,min_value_grid = -MM,max_value_grid = MM,mu_0 = 2,var_0 = 2)
pareto_theo_mu_k = as.data.frame(pareto_theo_mu_k)

# Plot
NC_plot_chisq = ggplot(chisq_theo_mu_k, aes(x = k_vec, y = mu_theo_vec)) +
  geom_point(size = 4, color = "blue") +  # Enlarged and colored points
  theme_light() + 
  geom_hline(yintercept = chisq_mean, linetype = "dashed", color = "red") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22),  # Center and bold title 
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20), 
    axis.text = element_text(size = 20),     # Increase axis text size
    axis.text.x = element_text(size = 20),   # Increase x-axis text size
    axis.text.y = element_text(size = 20),   # Increase y-axis text size
  ) +
  ggplot2::labs(
    title = "Chi-Square",   # Title of the plot
    y = "Theoretical Mean",  # Y-axis label
    x = "k"  # Y-axis label
  ) + 
  ylim(2.5,3.25)

NC_plot_lnorm = ggplot(lnorm_theo_mu_k, aes(x = k_vec, y = mu_theo_vec)) +
  geom_point(size = 4, color = "blue") +  # Enlarged and colored points
  theme_light() + 
  geom_hline(yintercept = lnorm_mean, linetype = "dashed", color = "red") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22),  # Center and bold title 
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20), 
    axis.text = element_text(size = 20),     # Increase axis text size
    axis.text.x = element_text(size = 20),   # Increase x-axis text size
    axis.text.y = element_text(size = 20),   # Increase y-axis text size
  ) +
  ggplot2::labs(
    title = "Log-Normal",   # Title of the plot
    y = "Theoretical Mean",  # Y-axis label
    x = "k"  # Y-axis label
  ) + 
  ylim(2,7.5)

NC_plot_pareto = ggplot(pareto_theo_mu_k, aes(x = k_vec, y = mu_theo_vec)) +
  geom_point(size = 4, color = "blue") +  # Enlarged and colored points
  theme_light() + 
  geom_hline(yintercept = pareto_mean, linetype = "dashed", color = "red") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22),  # Center and bold title 
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20), 
    axis.text = element_text(size = 20),     # Increase axis text size
    axis.text.x = element_text(size = 20),   # Increase x-axis text size
    axis.text.y = element_text(size = 20),   # Increase y-axis text size
  ) +
  ggplot2::labs(
    title = "Pareto",   # Title of the plot
    y = "Theoretical Mean",  # Y-axis label
    x = "k"  # Y-axis label
  ) +
  ylim(1,4.5)


library(patchwork)

grDevices::pdf("Minimum-k.pdf",width = 12, height = 10)
k_bound_plot
dev.off()

grDevices::pdf("Non-Convergence-k.pdf",width = 12, height = 10)
NC_plot_chisq / NC_plot_lnorm / NC_plot_pareto
dev.off()

grDevices::pdf("Minimum-n.pdf",width = 12, height = 10)
n_bound_plot_norm / n_bound_plot_beta / n_bound_plot_gamma
dev.off()


# Play with the sigma on estimated and theoretical mu 

Emp_mu_sigma_relationship = function(p0,x,grid,mu_0,sigma,real_mu,range = 10,range_diff = 0.5){
  
  mu_est = rep(0,range*2 / range_diff)
  mu_theo = rep(0,range*2 / range_diff)
  sigma_0 = rep(0,range*2 / range_diff)
  
  for (i in 1:(range*2 / range_diff)) {
    
    sigma_0[i] = sigma - range + range_diff*i
    mu_est[i] = EMemp(x = x, grid = grid, start=mu_0,sigma = sigma_0[i])$mu_hat
    mu_theo[i] = EMtheo(p0, grid, start=mu_0,sigma = sigma_0[i])$mu_hat
    
  }
  
  return(list("Mu_est" = mu_est,"Mu_theo" = mu_theo, "Sigma_0" = sigma_0))
  
}

# K = 3 (N = 500)

# normal example

n <- 500
grid <- seq(-15,15, length.out=4)
x <- rnorm(n,sd = 3)
p0 <- diff(pnorm(grid,sd = 3))
delta <- min(diff(grid))
y <- floor(x/delta)
y = y*delta
mu_0 = (mean(y/delta)+0.5)*delta
sigma_0 = sd(y)
normal_mean = 0
normal_sigma = 3

mu_list = Emp_mu_sigma_relationship(p0,x,grid,mu_0,sigma = normal_sigma,real_mu = normal_mean,range = 40,range_diff = 2)
mu_emp_norm_n1 = mu_list$Mu_est
mu_theo_norm = mu_list$Mu_theo
diff_sigma_values_norm_n1 = mu_list$Sigma_0

data_norm <- data.frame(
  Sigma = diff_sigma_values_norm_n1,
  MuEstn1 = mu_emp_norm_n1,
  MuTheo = mu_theo_norm
)

library(grid)  # Needed for unit() function
library(dplyr)  # Needed for the filter function

# Filter out rows where Sigma is equal to 0 (the point where it's undefined)
data_norm_filtered <- data_norm %>%
  filter(Sigma != 0)

# Create the plot with the filtered data
plot_data_norm_k1 <- ggplot(data_norm_filtered) +
  geom_point(aes(x = Sigma, y = MuEstn1, color = "Empirical Mu (N=500)"),size = 3) +
  geom_point(aes(x = Sigma, y = MuTheo, color = "Theoretical Mu"), size = 3) +
  geom_hline(yintercept = normal_mean, color = "red", linetype = "dashed") +
  geom_vline(xintercept = normal_sigma, color = "navy", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  ggplot2::labs(
    y = "Mu",
    x = "Sigma"
  ) +
  scale_color_manual(
    values = c(
      "Empirical Mu (N=500)" = "darkgreen",
      "Theoretical Mu" = "brown3"
    ),
    breaks = c("Empirical Mu (N=500)", "Theoretical Mu")
  ) +               
  theme_minimal() +
  theme(
    legend.position = "none"
  ) +
  ylim(-5,5)
# Add annotation for the navy vertical line (True Sigma)
# annotate("text", x = normal_sigma, y = 3.5, label = "True Sigma", color = "navy", hjust = -0.1, size = 5)

plot_data_norm_k1


#----------------------------

# beta example

set.seed(10)
n <- 500
alpha <- 2
beta <- 5
grid <- seq(0,1.3, by=0.32)
x <- rbeta(n,alpha,beta)
delta <- min(diff(grid))
y <- floor(x/delta)
y = y*delta
mu_0 = (mean(y/delta)+0.5)*delta
sd = sd(y)
beta_mean = alpha / (alpha + beta)
beta_sd = sqrt(alpha*beta / ((alpha + beta)^2 * (alpha + beta +1 ) ) )

p0 <- diff(pbeta(grid,alpha,beta))

mu_list = Emp_mu_sigma_relationship(p0 = p0,x = x,grid = grid,mu_0 = mu_0,sigma = beta_sd,real_mu = beta_mean,range = 2,range_diff = 0.01)
mu_emp_beta = mu_list$Mu_est
mu_theo_beta = mu_list$Mu_theo
sigma_values_beta = mu_list$Sigma_0

data_beta <- data.frame(
  Sigma = sigma_values_beta,
  MuEstn1 = mu_emp_beta,
  MuTheo = mu_theo_beta
)

library(grid)  # Needed for unit() function
library(dplyr)  # Needed for the filter function

# Filter out rows where Sigma is equal to 0 (the point where it's undefined)
data_beta_filtered <- data_beta %>%
  filter(Sigma > 0)

# Create the plot with the filtered data
plot_data_beta_k1 <- ggplot(data_beta_filtered) +
  geom_line(aes(x = Sigma, y = MuEstn1, color = "Empirical Mu (N=500)"),linewidth = 1.5) +
  geom_line(aes(x = Sigma, y = MuTheo, color = "Theoretical Mu"), linewidth = 1.5) +
  geom_hline(yintercept = beta_mean, color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = beta_sd, color = "navy", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = delta/2, color = "orange", linetype = "dashed", linewidth = 1) +
  ggplot2::labs(
    y = "",
    x = ""
  ) +
  scale_color_manual(
    values = c(
      "Empirical Mu (N=500)" = "darkgreen",
      "Theoretical Mu" = "brown3"
    ),
    breaks = c("Empirical Mu (N=500)", "Theoretical Mu")
  ) +               
  theme_minimal() +
  theme(
    legend.position = "none"
  ) +
  ylim(0.223,0.35) + 
  theme(legend.position = "none",  # Increase title size
        axis.title.x = element_text(size = 22),  # Increase x-axis title size
        axis.title.y = element_text(size = 22),  # Increase y-axis title size
        axis.text.x = element_text(size = 22),  # Increase x-axis text size
        axis.text.y = element_text(size = 22),
        strip.text = element_text(size = 22))

# Add annotation for the navy vertical line (True Sigma)
# annotate("text", x = beta_sd, y = 0.34, label = "True Sigma", color = "navy", hjust = -0.1, size = 5)

plot_data_beta_k1

#----------------------------

# gamma example

set.seed(10)
n <- 500
alpha <- 2
beta <- 5
MM <- ceiling(qgamma(0.9999, alpha, beta))+5
grid <- seq(0,MM, by = 0.57)
x <- rgamma(n,alpha,beta)
delta <- min(diff(grid))
y <- floor(x/delta)
y = y*delta
mu_0 = (mean(y/delta)+0.5)*delta
sd = sd(y)
gamma_mean = alpha / beta
gamma_sd = sqrt(alpha) / beta

p0 <- diff(pgamma(grid,alpha,beta))

mu_list = Emp_mu_sigma_relationship(p0 = p0,x = x,grid = grid,mu_0 = mu_0,sigma = gamma_sd,real_mu = gamma_mean,range = 5,range_diff = 0.01)
mu_emp_gamma = mu_list$Mu_est
mu_theo_gamma = mu_list$Mu_theo
sigma_values_gamma = mu_list$Sigma_0

data_gamma <- data.frame(
  Sigma = sigma_values_gamma,
  MuEstn1 = mu_emp_gamma,
  MuTheo = mu_theo_gamma
)

library(grid)  # Needed for unit() function
library(dplyr)  # Needed for the filter function

# Filter out rows where Sigma is equal to 0 (the point where it's undefined)
data_gamma_filtered <- data_gamma %>%
  filter(Sigma > 0)

# Create the plot with the filtered data
plot_data_gamma_k1 <- ggplot(data_gamma_filtered) +
  geom_line(aes(x = Sigma, y = MuEstn1, color = "Empirical Mu (N=500)"),linewidth = 1.5) +
  geom_line(aes(x = Sigma, y = MuTheo, color = "Theoretical Mu"), linewidth = 1.5) +
  geom_hline(yintercept = gamma_mean, color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = gamma_sd, color = "navy", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = delta/2, color = "orange", linetype = "dashed", linewidth = 1) +
  ggplot2::labs(
    y = "",
    x = ""
  ) +
  scale_color_manual(
    values = c(
      "Empirical Mu (N=500)" = "darkgreen",
      "Theoretical Mu" = "brown3"
    ),
    breaks = c("Empirical Mu (N=500)", "Theoretical Mu")
  ) +               
  theme_minimal() +
  theme(
    legend.position = "none"
  ) + 
  theme(legend.position = "none",  # Increase title size
        axis.title.x = element_text(size = 22),  # Increase x-axis title size
        axis.title.y = element_text(size = 22),  # Increase y-axis title size
        axis.text.x = element_text(size = 22),  # Increase x-axis text size
        axis.text.y = element_text(size = 22),
        strip.text = element_text(size = 22))

# Add annotation for the navy vertical line (True Sigma)
# annotate("text", x = gamma_sd, y = 10.34, label = "True Sigma", color = "navy", hjust = -0.1, size = 5)

plot_data_gamma_k1

# K = 9 (N = 500)

# normal example

n <- 500
grid <- seq(-15,15, length.out=10)
x <- rnorm(n,sd = 3)
p0 <- diff(pnorm(grid,sd = 3))
delta <- min(diff(grid))
y <- floor(x/delta)
y = y*delta
mu_0 = (mean(y/delta)+0.5)*delta
sigma_0 = sd(y)
normal_mean = 0
normal_sigma = 3

mu_list = Emp_mu_sigma_relationship(p0,x,grid,mu_0,sigma = normal_sigma,real_mu = normal_mean,range = 40,range_diff = 2)
mu_emp_norm_n1 = mu_list$Mu_est
mu_theo_norm = mu_list$Mu_theo
diff_sigma_values_norm_n1 = mu_list$Sigma_0

data_norm <- data.frame(
  Sigma = diff_sigma_values_norm_n1,
  MuEstn1 = mu_emp_norm_n1,
  MuTheo = mu_theo_norm
)

library(grid)  # Needed for unit() function
library(dplyr)  # Needed for the filter function

# Filter out rows where Sigma is equal to 0 (the point where it's undefined)
data_norm_filtered <- data_norm %>%
  filter(Sigma != 0)

# Create the plot with the filtered data
plot_data_norm_k2 <- ggplot(data_norm_filtered) +
  geom_point(aes(x = Sigma, y = MuEstn1, color = "Empirical Mu (N=500)"),size = 3) +
  geom_point(aes(x = Sigma, y = MuTheo, color = "Theoretical Mu"), size = 3) +
  geom_hline(yintercept = normal_mean, color = "red", linetype = "dashed") +
  geom_vline(xintercept = normal_sigma, color = "navy", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  ggplot2::labs(
    y = "Mu",
    x = "Sigma"
  ) +
  scale_color_manual(
    values = c(
      "Empirical Mu (N=500)" = "darkgreen",
      "Theoretical Mu" = "brown3"
    ),
    breaks = c("Empirical Mu (N=500)", "Theoretical Mu")
  ) +               
  theme_minimal() +
  theme(
    legend.position = "none"
  ) +
  ylim(-0.3,0.3)
# Add annotation for the navy vertical line (True Sigma)
# annotate("text", x = normal_sigma, y = 0.2, label = "True Sigma", color = "navy", hjust = -0.1, size = 5)

plot_data_norm_k2


#----------------------------

# beta example

set.seed(10)
n <- 500
alpha <- 2
beta <- 5
grid <- seq(0,1.2, by=0.16)
x <- rbeta(n,alpha,beta)
delta <- min(diff(grid))
y <- floor(x/delta)
y = y*delta
mu_0 = (mean(y/delta)+0.5)*delta
sd = sd(y)
beta_mean = alpha / (alpha + beta)
beta_sd = sqrt(alpha*beta / ((alpha + beta)^2 * (alpha + beta +1 ) ) )

p0 <- diff(pbeta(grid,alpha,beta))

mu_list = Emp_mu_sigma_relationship(p0 = p0,x = x,grid = grid,mu_0 = mu_0,sigma = beta_sd,real_mu = beta_mean,range = 2,range_diff = 0.01)
mu_emp_beta = mu_list$Mu_est
mu_theo_beta = mu_list$Mu_theo
sigma_values_beta = mu_list$Sigma_0

data_beta <- data.frame(
  Sigma = sigma_values_beta,
  MuEstn1 = mu_emp_beta,
  MuTheo = mu_theo_beta
)

library(grid)  # Needed for unit() function
library(dplyr)  # Needed for the filter function

# Filter out rows where Sigma is equal to 0 (the point where it's undefined)
data_beta_filtered <- data_beta %>%
  filter(Sigma > 0)

# Create the plot with the filtered data
plot_data_beta_k2 <- ggplot(data_beta_filtered) +
  geom_line(aes(x = Sigma, y = MuEstn1, color = "Empirical Mu (N=500)"),linewidth = 1.5) +
  geom_line(aes(x = Sigma, y = MuTheo, color = "Theoretical Mu"), linewidth = 1.5) +
  geom_hline(yintercept = beta_mean, color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = beta_sd, color = "navy", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = delta/2, color = "orange", linetype = "dashed", linewidth = 1) +
  ggplot2::labs(
    y = "",
    x = ""
  ) +
  scale_color_manual(
    values = c(
      "Empirical Mu (N=500)" = "darkgreen",
      "Theoretical Mu" = "brown3"
    ),
    breaks = c("Empirical Mu (N=500)", "Theoretical Mu")
  ) +               
  theme_minimal() +
  theme(
    legend.position = "none"
  ) +
  ylim(0.223,0.35) + 
  theme(legend.position = "none",  # Increase title size
        axis.title.x = element_text(size = 22),  # Increase x-axis title size
        axis.title.y = element_text(size = 22),  # Increase y-axis title size
        axis.text.x = element_text(size = 22),  # Increase x-axis text size
        axis.text.y = element_text(size = 22),
        strip.text = element_text(size = 22))
# Add annotation for the navy vertical line (True Sigma)
# annotate("text", x = beta_sd, y = 0.5, label = "True Sigma", color = "navy", hjust = -0.1, size = 5)

plot_data_beta_k2

#----------------------------

# gamma example

set.seed(10)

n <- 500
alpha <- 2
beta <- 5
MM <- ceiling(qgamma(0.9999, alpha, beta))+5
grid <- seq(0,MM, by = 0.28)
x <- rgamma(n,alpha,beta)
delta <- min(diff(grid))
y <- floor(x/delta)
y = y*delta
mu_0 = (mean(y/delta)+0.5)*delta
sd = sd(y)
gamma_mean = alpha / beta
gamma_sd = sqrt(alpha) / beta

p0 <- diff(pgamma(grid,alpha,beta))

mu_list = Emp_mu_sigma_relationship(p0 = p0,x = x,grid = grid,mu_0 = mu_0,sigma = gamma_sd,real_mu = gamma_mean,range = 5,range_diff = 0.1)
mu_emp_gamma = mu_list$Mu_est
mu_theo_gamma = mu_list$Mu_theo
sigma_values_gamma = mu_list$Sigma_0

data_gamma <- data.frame(
  Sigma = sigma_values_gamma,
  MuEstn1 = mu_emp_gamma,
  MuTheo = mu_theo_gamma
)

library(grid)  # Needed for unit() function
library(dplyr)  # Needed for the filter function

# Filter out rows where Sigma is equal to 0 (the point where it's undefined)
data_gamma_filtered <- data_gamma %>%
  filter(Sigma > 0)

# Create the plot with the filtered data
plot_data_gamma_k2 <- ggplot(data_gamma_filtered) +
  geom_line(aes(x = Sigma, y = MuEstn1, color = "Empirical Mu (N=500)"),linewidth = 1.5) +
  geom_line(aes(x = Sigma, y = MuTheo, color = "Theoretical Mu"), linewidth = 1.5) +
  geom_hline(yintercept = gamma_mean, color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = gamma_sd, color = "navy", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = delta/2, color = "orange", linetype = "dashed", linewidth = 1) +
  ggplot2::labs(
    y = "",
    x = ""
  ) +
  scale_color_manual(
    values = c(
      "Empirical Mu (N=500)" = "darkgreen",
      "Theoretical Mu" = "brown3"
    ),
    breaks = c("Empirical Mu (N=500)", "Theoretical Mu")
  ) +               
  theme_minimal() +
  theme(
    legend.position = "none"
  ) + 
  theme(legend.position = "none",  # Increase title size
        axis.title.x = element_text(size = 22),  # Increase x-axis title size
        axis.title.y = element_text(size = 22),  # Increase y-axis title size
        axis.text.x = element_text(size = 22),  # Increase x-axis text size
        axis.text.y = element_text(size = 22),
        strip.text = element_text(size = 22))

# Add annotation for the navy vertical line (True Sigma)
# annotate("text", x = gamma_sd, y = 10.34, label = "True Sigma", color = "navy", hjust = -0.1, size = 5)

plot_data_gamma_k2


# K = 15 (N = 500)

# normal example

n <- 500
grid <- seq(-15,15, length.out=16)
x <- rnorm(n,sd = 3)
p0 <- diff(pnorm(grid,sd = 3))
delta <- min(diff(grid))
y <- floor(x/delta)
y = y*delta
mu_0 = (mean(y/delta)+0.5)*delta
sigma_0 = sd(y)
normal_mean = 0
normal_sigma = 3

mu_list = Emp_mu_sigma_relationship(p0,x,grid,mu_0,sigma = normal_sigma,real_mu = normal_mean,range = 40,range_diff = 2)
mu_emp_norm_n1 = mu_list$Mu_est
mu_theo_norm = mu_list$Mu_theo
diff_sigma_values_norm_n1 = mu_list$Sigma_0

data_norm <- data.frame(
  Sigma = diff_sigma_values_norm_n1,
  MuEstn1 = mu_emp_norm_n1,
  MuTheo = mu_theo_norm
)

library(grid)  # Needed for unit() function
library(dplyr)  # Needed for the filter function

# Filter out rows where Sigma is equal to 0 (the point where it's undefined)
data_norm_filtered <- data_norm %>%
  filter(Sigma != 0)

# Create the plot with the filtered data
plot_data_norm_k3 <- ggplot(data_norm_filtered) +
  geom_point(aes(x = Sigma, y = MuEstn1, color = "Empirical Mu (N=500)"),size = 3) +
  geom_point(aes(x = Sigma, y = MuTheo, color = "Theoretical Mu"), size = 3) +
  geom_hline(yintercept = normal_mean, color = "red", linetype = "dashed") +
  geom_vline(xintercept = normal_sigma, color = "navy", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  ggplot2::labs(
    y = "Mu",
    x = "Sigma"
  ) +
  scale_color_manual(
    values = c(
      "Empirical Mu (N=500)" = "darkgreen",
      "Theoretical Mu" = "brown3"
    ),
    breaks = c("Empirical Mu (N=500)", "Theoretical Mu")
  ) +               
  theme_minimal() +
  theme(
    legend.position = "none"
  ) +
  ylim(-0.3,0.3)
# Add annotation for the navy vertical line (True Sigma)
# annotate("text", x = normal_sigma, y = 0.2, label = "True Sigma", color = "navy", hjust = -0.1, size = 5)

plot_data_norm_k3


#----------------------------

# beta example

n <- 500
alpha <- 2
beta <- 5
grid <- seq(0,1, length.out=16)
x <- rbeta(n,alpha,beta)
delta <- min(diff(grid))
y <- floor(x/delta)
y = y*delta
mu_0 = (mean(y/delta)+0.5)*delta
sd = sd(y)
beta_mean = alpha / (alpha + beta)
beta_sd = sqrt(alpha*beta / ((alpha + beta)^2 * (alpha + beta +1 ) ) )

p0 <- diff(pbeta(grid,alpha,beta))

mu_list = Emp_mu_sigma_relationship(p0 = p0,x = x,grid = grid,mu_0 = mu_0,sigma = beta_sd,real_mu = beta_mean,range = 2,range_diff = 0.12)
mu_emp_beta = mu_list$Mu_est
mu_theo_beta = mu_list$Mu_theo
sigma_values_beta = mu_list$Sigma_0

data_beta <- data.frame(
  Sigma = sigma_values_beta,
  MuEstn1 = mu_emp_beta,
  MuTheo = mu_theo_beta
)

library(grid)  # Needed for unit() function
library(dplyr)  # Needed for the filter function

# Filter out rows where Sigma is equal to 0 (the point where it's undefined)
data_beta_filtered <- data_beta %>%
  filter(Sigma != 0)

# Create the plot with the filtered data
plot_data_beta_k3 <- ggplot(data_beta_filtered) +
  geom_point(aes(x = Sigma, y = MuEstn1, color = "Empirical Mu (N=500)"),size = 3) +
  geom_point(aes(x = Sigma, y = MuTheo, color = "Theoretical Mu"), size = 3) +
  geom_hline(yintercept = beta_mean, color = "red", linetype = "dashed") +
  geom_vline(xintercept = beta_sd, color = "navy", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  ggplot2::labs(
    y = "Mu",
    x = "Sigma"
  ) +
  scale_color_manual(
    values = c(
      "Empirical Mu (N=500)" = "darkgreen",
      "Theoretical Mu" = "brown3"
    ),
    breaks = c("Empirical Mu (N=500)", "Theoretical Mu")
  ) +               
  theme_minimal() +
  theme(
    legend.position = "none"
  ) +
  ylim(0.223,0.35)
# Add annotation for the navy vertical line (True Sigma)
# annotate("text", x = beta_sd, y = 0.5, label = "True Sigma", color = "navy", hjust = -0.1, size = 5)

plot_data_beta_k3


#----------------------------

# gamma example

n <- 500
alpha <- 4
beta <- 0.5
MM <- ceiling(qgamma(0.9999, alpha, beta))+5
grid <- seq(0,MM, length.out=16)
x <- rgamma(n,alpha,beta)
delta <- min(diff(grid))
y <- floor(x/delta)
y = y*delta
mu_0 = (mean(y/delta)+0.5)*delta
sd = sd(y)
gamma_mean = alpha / beta
gamma_sd = sqrt(alpha) / beta

p0 <- diff(pgamma(grid,alpha,beta))

mu_list = Emp_mu_sigma_relationship(p0 = p0,x = x,grid = grid,mu_0 = mu_0,sigma = gamma_sd,real_mu = gamma_mean,range = 30,range_diff = 2)
mu_emp_gamma = mu_list$Mu_est
mu_theo_gamma = mu_list$Mu_theo
sigma_values_gamma = mu_list$Sigma_0

data_gamma <- data.frame(
  Sigma = sigma_values_gamma,
  MuEstn1 = mu_emp_gamma,
  MuTheo = mu_theo_gamma
)

library(grid)  # Needed for unit() function
library(dplyr)  # Needed for the filter function

# Filter out rows where Sigma is equal to 0 (the point where it's undefined)
data_gamma_filtered <- data_gamma %>%
  filter(Sigma != 0)

# Create the plot with the filtered data
plot_data_gamma_k3 <- ggplot(data_gamma_filtered) +
  geom_point(aes(x = Sigma, y = MuEstn1, color = "Empirical Mu (N=500)"),size = 3) +
  geom_point(aes(x = Sigma, y = MuTheo, color = "Theoretical Mu"), size = 3) +
  geom_hline(yintercept = gamma_mean, color = "red", linetype = "dashed") +
  geom_vline(xintercept = gamma_sd, color = "navy", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  ggplot2::labs(
    y = "Mu",
    x = "Sigma"
  ) +
  scale_color_manual(
    values = c(
      "Empirical Mu (N=500)" = "darkgreen",
      "Theoretical Mu" = "brown3"
    ),
    breaks = c("Empirical Mu (N=500)", "Theoretical Mu")
  ) +               
  theme_minimal() +
  theme(
    legend.position = "none"
  ) +
  ylim(4,12)
# Add annotation for the navy vertical line (True Sigma)
#annotate("text", x = gamma_sd, y = 10.34, label = "True Sigma", color = "navy", hjust = -0.1, size = 5)

plot_data_gamma_k3


# Arrange the plots in a 3x3 grid
final_plot <- (plot_data_norm_k1 | plot_data_norm_k2 | plot_data_norm_k3) /
  (plot_data_beta_k1 | plot_data_beta_k2 | plot_data_beta_k3) /
  (plot_data_gamma_k1 | plot_data_gamma_k2 | plot_data_gamma_k3)

# Display the final combined plot
final_plot

library(ggplot2)
library(cowplot)

# Create label plots for each row and column to insert into the grid
label_k3 <- ggdraw() + draw_label("  k = 3", fontface = 'bold', size = 15)
label_k7 <- ggdraw() + draw_label("  k = 9", fontface = 'bold', size = 15)
label_k11 <- ggdraw() + draw_label("  k = 15", fontface = 'bold', size = 15)

label_normal <- ggdraw() + draw_label("      Normal", fontface = 'bold', size = 12, angle = 90)
label_beta <- ggdraw() + draw_label("    Beta", fontface = 'bold', size = 12, angle = 90)
label_gamma <- ggdraw() + draw_label("      Gamma", fontface = 'bold', size = 12, angle = 90)

# Combine the plots in a 3x3 grid with row and column labels
final_plot <- plot_grid(
  NULL, label_k3, label_k7, label_k11,               # Top row with column labels
  label_normal, plot_data_norm_k1, plot_data_norm_k2, plot_data_norm_k3, # First row with 'Normal'
  label_beta, plot_data_beta_k1, plot_data_beta_k2, plot_data_beta_k3,   # Second row with 'Beta'
  label_gamma, plot_data_gamma_k1, plot_data_gamma_k2, plot_data_gamma_k3, # Third row with 'Gamma'
  ncol = 4, nrow = 4,                                # Define a 4x4 grid
  rel_widths = c(0.1, 1, 1, 1),                      # Control widths
  rel_heights = c(0.1, 1, 1, 1)                      # Control heights
)

# Display the final combined plot
print(final_plot)


grDevices::pdf("mu-sigma-k.pdf",width = 14, height = 10)
final_plot
dev.off()

grDevices::pdf("Beta-STDelta-1-mu-sigma-k.pdf",width = 12, height = 10)
plot_data_beta_k1
dev.off()

grDevices::pdf("Beta-STDelta-2-mu-sigma-k.pdf",width = 12, height = 10)
plot_data_beta_k2
dev.off()

grDevices::pdf("Gamma-STDelta-1-mu-sigma-k.pdf",width = 12, height = 10)
plot_data_gamma_k1
dev.off()

grDevices::pdf("Gamma-STDelta-2-mu-sigma-k.pdf",width = 12, height = 10)
plot_data_gamma_k2
dev.off()


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

grDevices::pdf("Theo_mu-delta.pdf",width = 16, height = 14)

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


# Chi-Square - No Convergence

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

# Log-Normal - No Convergence

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


# Pareto - No Convergence

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

dev.off()

n <- 1000
grid <- seq(0,60, length.out=121)
x <- rgamma(n, shape=1, scale=5)
mu_0 = 1
gamma_mean = 5
gamma_sd = sqrt(1) * 5 

# Trial-123

EMemp(x,grid,start = mu_0,sigma = -24)$mu_hat
EMemp(x,grid,start = mu_0,sigma = 0.4)$mu_hat
EMemp(x,grid,start = mu_0,sigma = 0.8)$mu_hat
EMemp(x,grid,start = mu_0,sigma = 1)$mu_hat
EMemp(x,grid,start = mu_0,sigma = 2)$mu_hat
EMemp(x,grid,start = mu_0,sigma = 3)$mu_hat
EMemp(x,grid,start = mu_0,sigma = 4)$mu_hat
EMemp(x,grid,start = mu_0,sigma = 5)$mu_hat
EMemp(x,grid,start = mu_0,sigma = 6)$mu_hat
EMemp(x,grid,start = mu_0,sigma = 24)$mu_hat
EMemp(x,grid,start = mu_0,sigma = 50)$mu_hat
EMemp(x,grid,start = mu_0,sigma = 60)$mu_hat
EMemp(x,grid,start = mu_0,sigma = 100)$mu_hat
EMemp(x,grid,start = mu_0,sigma = gamma_sd)$mu_hat


# Gamma in my example 

