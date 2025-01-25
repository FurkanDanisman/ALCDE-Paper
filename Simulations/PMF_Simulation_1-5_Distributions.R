# All distributions together n = 1000, n = 10000

# Each Distribution, Sample size increasing, grid width increasing 

# Parameters 

alpha=2;beta=5;
scale_val = 0.05;df=3;
shape=2;location=2;
mu = 0;sigma=1;
mul = 1;
df=3;
shape_p = 4;

# CDFs 

cdf_beta <- function(x) {
  pbeta(x,alpha,beta)
}

cdf_norm <- function(x) {
  pnorm(x, mean = mu, sd = sigma)
}

cdf_gamma <- function(x) {
  pgamma(x, alpha, beta)
}

cdf_logistic <- function(x) {
  plogis(x)
}

cdf_cauchy <- function(x) {
  pcauchy(x,scale = scale_val)
}

cdf_t <- function(x) {
  pt(x,df)
}

# Standardized Delta Values # 

norm_deltas  = c(0.5, 1.0, 1.5, 2.0, 2.5)
beta_deltas  = c(0.0799, 0.1597, 0.2396, 0.3194, 0.3993)
gamma_deltas = c(0.1414, 0.2828, 0.4243, 0.5657, 0.7071)
logistic_deltas = c(0.9069, 1.8138, 2.7207, 3.6276, 4.5345)
t_deltas = c(0.8660, 1.7321, 2.5981, 3.4641, 4.3301)

# Corresponding k values to Standardized Delta Values # 

norm_ks  = c(12, 6, 4, 3, 2)
beta_ks  = c(22, 11, 7, 5, 4)
gamma_ks = c(56, 28, 18, 14, 11)
logistic_ks = c(35, 17, 11, 8, 7)
t_ks = c(66, 33, 22, 16, 13)

# Simulations Within the Distribution 

colors <- c("#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#D55E00")

# Colors for Normal Distribution 

c11 = darken("#92CAEB",factor = 1.1)
c21 = darken(c11,factor = 1.1)
c31 = darken(c21,factor = 1.1)
c41 = darken(c31,factor = 1.1)
colors1 = c("#92CAEB",c11,c21,c31,c41)

# Colors for Beta Distribution

c11 = darken("#F3CF70",factor = 1.1)
c21 = darken(c11,factor = 1.1)
c31 = darken(c21,factor = 1.1)
c41 = darken(c31,factor = 1.1)
colors2 = c("#F3CF70",c11,c21,c31,c41)

# Colors for Gamma Distribution 

c11 = darken("#66D7A5",factor = 1.1)
c21 = darken(c11,factor = 1.1)
c31 = darken(c21,factor = 1.1)
c41 = darken(c31,factor = 1.1)
colors3 = c("#66D7A5",c11,c21,c31,c41)

# Colors for Logistic Distribution

c11 = darken("#E6A4C6",factor = 1.1)
c21 = darken(c11,factor = 1.1)
c31 = darken(c21,factor = 1.1)
c41 = darken(c31,factor = 1.1)
colors4 = c("#E6A4C6",c11,c21,c31,c41)

# Colors for Student's t Distribution

c11 = darken("#D55E00",factor = 1.1)
c21 = darken(c11,factor = 1.1)
c31 = darken(c21,factor = 1.1)
c41 = darken(c31,factor = 1.1)
colors5 = c("#D55E00",c11,c21,c31,c41)


n = c(10^1,10^2,10^3,10^4,10^5)
B = 100

# Normal Distribution - Sample Size Simulation 

MM = qnorm(0.9999999,mu,sigma) + 3
norm_grid <- seq(-MM,MM, by = norm_deltas[1])

res_norm_pmf_n_11 = rep(0,B)
res_norm_pmf_n_21 = rep(0,B)
res_norm_pmf_n_31 = rep(0,B)
res_norm_pmf_n_41 = rep(0,B)
res_norm_pmf_n_51 = rep(0,B)

res_norm_pmf_n_12 = rep(0,B)
res_norm_pmf_n_22 = rep(0,B)
res_norm_pmf_n_32 = rep(0,B)
res_norm_pmf_n_42 = rep(0,B)
res_norm_pmf_n_52 = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  # Normal Distribution - n1
  
  x <- sort(rnorm(n[1], mean = mu,sd=sigma))
  n_i <- hist(x, breaks = norm_grid, plot = F)$counts
  
  res_norm_pmf_n_11[i] = L1_L2_pmf(n_i,norm_grid,cdf = cdf_norm)$L1_pmf
  res_norm_pmf_n_12[i] = L1_L2_pmf(n_i,norm_grid,cdf = cdf_norm)$L2_pmf
  
  # Normal Distribution - n2
  
  x <- sort(rnorm(n[2], mean = mu,sd=sigma))
  n_i <- hist(x, breaks = norm_grid, plot = F)$counts
  
  res_norm_pmf_n_21[i] = L1_L2_pmf(n_i,norm_grid,cdf = cdf_norm)$L1_pmf
  res_norm_pmf_n_22[i] = L1_L2_pmf(n_i,norm_grid,cdf = cdf_norm)$L2_pmf
  
  # Normal Distribution - n3
  
  x <- sort(rnorm(n[3], mean = mu,sd=sigma))
  n_i <- hist(x, breaks = norm_grid, plot = F)$counts
  
  res_norm_pmf_n_31[i] = L1_L2_pmf(n_i,norm_grid,cdf = cdf_norm)$L1_pmf
  res_norm_pmf_n_32[i] = L1_L2_pmf(n_i,norm_grid,cdf = cdf_norm)$L2_pmf
  
  # Normal Distribution - n4
  
  x <- sort(rnorm(n[4], mean = mu,sd=sigma))
  n_i <- hist(x, breaks = norm_grid, plot = F)$counts
  
  res_norm_pmf_n_41[i] = L1_L2_pmf(n_i,norm_grid,cdf = cdf_norm)$L1_pmf
  res_norm_pmf_n_42[i] = L1_L2_pmf(n_i,norm_grid,cdf = cdf_norm)$L2_pmf
  
  # Normal Distribution - n5
  
  x <- sort(rnorm(n[5], mean = mu,sd=sigma))
  n_i <- hist(x, breaks = norm_grid, plot = F)$counts
  
  res_norm_pmf_n_51[i] = L1_L2_pmf(n_i,norm_grid,cdf = cdf_norm)$L1_pmf
  res_norm_pmf_n_52[i] = L1_L2_pmf(n_i,norm_grid,cdf = cdf_norm)$L2_pmf
  
}

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L1_n/L1-PMF-1-norm-n-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_norm_pmf_n_11,res_norm_pmf_n_21,res_norm_pmf_n_31,res_norm_pmf_n_41,res_norm_pmf_n_51,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors1)
segments(x0 = 0.29, y0 =0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L2_n/L2-PMF-1-norm-n-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_norm_pmf_n_12,res_norm_pmf_n_22,res_norm_pmf_n_32,res_norm_pmf_n_42,res_norm_pmf_n_52,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors1)
segments(x0 = 0.29, y0 =0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

# Beta Distribution - Sample Size Simulation 

beta_grid <- seq(0,1.04, by = beta_deltas[1])

res_beta_pmf_n_11 = rep(0,B)
res_beta_pmf_n_21 = rep(0,B)
res_beta_pmf_n_31 = rep(0,B)
res_beta_pmf_n_41 = rep(0,B)
res_beta_pmf_n_51 = rep(0,B)

res_beta_pmf_n_12 = rep(0,B)
res_beta_pmf_n_22 = rep(0,B)
res_beta_pmf_n_32 = rep(0,B)
res_beta_pmf_n_42 = rep(0,B)
res_beta_pmf_n_52 = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  # Beta Distribution - n1
  
  x   <- sort(rbeta(n[1],alpha,beta))
  n_i <- hist(x, breaks = beta_grid, plot = F)$counts
  
  res_beta_pmf_n_11[i] = L1_L2_pmf(n_i,beta_grid,cdf = cdf_beta)$L1_pmf
  res_beta_pmf_n_12[i] = L1_L2_pmf(n_i,beta_grid,cdf = cdf_beta)$L2_pmf
  
  # Beta Distribution - n2
  
  x   <- sort(rbeta(n[2],alpha,beta))
  n_i <- hist(x, breaks = beta_grid, plot = F)$counts
  
  
  res_beta_pmf_n_21[i] = L1_L2_pmf(n_i,beta_grid,cdf = cdf_beta)$L1_pmf
  res_beta_pmf_n_22[i] = L1_L2_pmf(n_i,beta_grid,cdf = cdf_beta)$L2_pmf
  
  # Beta Distribution - n3
  
  x   <- sort(rbeta(n[3],alpha,beta))
  n_i <- hist(x, breaks = beta_grid, plot = F)$counts
  
  res_beta_pmf_n_31[i] = L1_L2_pmf(n_i,beta_grid,cdf = cdf_beta)$L1_pmf
  res_beta_pmf_n_32[i] = L1_L2_pmf(n_i,beta_grid,cdf = cdf_beta)$L2_pmf
  
  # Beta Distribution - n4
  
  x   <- sort(rbeta(n[4],alpha,beta))
  n_i <- hist(x, breaks = beta_grid, plot = F)$counts
  
  res_beta_pmf_n_41[i] = L1_L2_pmf(n_i,beta_grid,cdf = cdf_beta)$L1_pmf
  res_beta_pmf_n_42[i] = L1_L2_pmf(n_i,beta_grid,cdf = cdf_beta)$L2_pmf
  
  # Beta Distribution - n5
  
  x   <- sort(rbeta(n[5],alpha,beta))
  n_i <- hist(x, breaks = beta_grid, plot = F)$counts
  
  res_beta_pmf_n_51[i] = L1_L2_pmf(n_i,beta_grid,cdf = cdf_beta)$L1_pmf
  res_beta_pmf_n_52[i] = L1_L2_pmf(n_i,beta_grid,cdf = cdf_beta)$L2_pmf
  
}

grDevices::pdf("L1-PMF-1-Beta-n-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pmf_n_11,res_beta_pmf_n_21,res_beta_pmf_n_31,res_beta_pmf_n_41,res_beta_pmf_n_51,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2,ylim=c(0,1))
segments(x0 = 0.29, y0 =0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L2_n/L2-PMF-1-Beta-n-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pmf_n_12,res_beta_pmf_n_22,res_beta_pmf_n_32,res_beta_pmf_n_42,res_beta_pmf_n_52,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2,ylim=c(0,0.6))
segments(x0 = 0.29, y0 =0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

# Gamma Distribution - Sample Size Simulation 

MM = qgamma(0.9999999,alpha,beta) + 3
gamma_grid <- seq(0,MM, by = gamma_deltas[1])

res_gamma_pmf_n_11 = rep(0,B)
res_gamma_pmf_n_21 = rep(0,B)
res_gamma_pmf_n_31 = rep(0,B)
res_gamma_pmf_n_41 = rep(0,B)
res_gamma_pmf_n_51 = rep(0,B)

res_gamma_pmf_n_12 = rep(0,B)
res_gamma_pmf_n_22 = rep(0,B)
res_gamma_pmf_n_32 = rep(0,B)
res_gamma_pmf_n_42 = rep(0,B)
res_gamma_pmf_n_52 = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  # Gamma Distribution - n1
  
  x   <- sort(rgamma(n[1],alpha,beta))
  n_i <- hist(x, breaks = gamma_grid, plot = F)$counts
  
  res_gamma_pmf_n_11[i] = L1_L2_pmf(n_i,gamma_grid,cdf = cdf_gamma)$L1_pmf
  res_gamma_pmf_n_12[i] = L1_L2_pmf(n_i,gamma_grid,cdf = cdf_gamma)$L2_pmf
  
  # Gamma Distribution - n2
  
  x   <- sort(rgamma(n[2],alpha,beta))
  n_i <- hist(x, breaks = gamma_grid, plot = F)$counts
  
  res_gamma_pmf_n_21[i] = L1_L2_pmf(n_i,gamma_grid,cdf = cdf_gamma)$L1_pmf
  res_gamma_pmf_n_22[i] = L1_L2_pmf(n_i,gamma_grid,cdf = cdf_gamma)$L2_pmf
  
  # Gamma Distribution - n3
  
  x   <- sort(rgamma(n[3],alpha,beta))
  n_i <- hist(x, breaks = gamma_grid, plot = F)$counts
  
  res_gamma_pmf_n_31[i] = L1_L2_pmf(n_i,gamma_grid,cdf = cdf_gamma)$L1_pmf
  res_gamma_pmf_n_32[i] = L1_L2_pmf(n_i,gamma_grid,cdf = cdf_gamma)$L2_pmf
  
  # Gamma Distribution - n4
  
  x   <- sort(rgamma(n[4],alpha,beta))
  n_i <- hist(x, breaks = gamma_grid, plot = F)$counts
  
  res_gamma_pmf_n_41[i] = L1_L2_pmf(n_i,gamma_grid,cdf = cdf_gamma)$L1_pmf
  res_gamma_pmf_n_42[i] = L1_L2_pmf(n_i,gamma_grid,cdf = cdf_gamma)$L2_pmf
  
  # Gamma Distribution - n5
  
  x   <- sort(rgamma(n[5],alpha,beta))
  n_i <- hist(x, breaks = gamma_grid, plot = F)$counts
  
  res_gamma_pmf_n_51[i] = L1_L2_pmf(n_i,gamma_grid,cdf = cdf_gamma)$L1_pmf
  res_gamma_pmf_n_52[i] = L1_L2_pmf(n_i,gamma_grid,cdf = cdf_gamma)$L2_pmf
  
}

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L1_n/L1-PMF-1-Gamma-n-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pmf_n_11,res_gamma_pmf_n_21,res_gamma_pmf_n_31,res_gamma_pmf_n_41,res_gamma_pmf_n_51,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3,ylim=c(0,1.2))
segments(x0 = 0.29, y0 =0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L2_n/L2-PMF-1-Gamma-n-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pmf_n_12,res_gamma_pmf_n_22,res_gamma_pmf_n_32,res_gamma_pmf_n_42,res_gamma_pmf_n_52,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3,ylim=c(0,0.8))
segments(x0 = 0.29, y0 =0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

# Logistic Distribution - Sample Size Simulation 

MM = qlogis(0.9999999,mu,sigma) + 10
logistic_grid <- seq(-MM,MM, by = logistic_deltas[1])

res_logis_pmf_n_11 = rep(0,B)
res_logis_pmf_n_21 = rep(0,B)
res_logis_pmf_n_31 = rep(0,B)
res_logis_pmf_n_41 = rep(0,B)
res_logis_pmf_n_51 = rep(0,B)

res_logis_pmf_n_12 = rep(0,B)
res_logis_pmf_n_22 = rep(0,B)
res_logis_pmf_n_32 = rep(0,B)
res_logis_pmf_n_42 = rep(0,B)
res_logis_pmf_n_52 = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  # Logistic Distribution - n1
  
  x = sort(rlogis(n[1]))
  n_i <- hist(x, breaks = logistic_grid, plot = F)$counts
  
  res_logis_pmf_n_11[i] = L1_L2_pmf(n_i,logistic_grid,cdf = cdf_logistic)$L1_pmf
  res_logis_pmf_n_12[i] = L1_L2_pmf(n_i,logistic_grid,cdf = cdf_logistic)$L2_pmf
  
  # Logistic Distribution - n2
  
  x = sort(rlogis(n[2]))
  n_i <- hist(x, breaks = logistic_grid, plot = F)$counts
  
  res_logis_pmf_n_21[i] = L1_L2_pmf(n_i,logistic_grid,cdf = cdf_logistic)$L1_pmf
  res_logis_pmf_n_22[i] = L1_L2_pmf(n_i,logistic_grid,cdf = cdf_logistic)$L2_pmf
  
  # Logistic Distribution - n3
  
  x = sort(rlogis(n[3]))
  n_i <- hist(x, breaks = logistic_grid, plot = F)$counts
  
  res_logis_pmf_n_31[i] = L1_L2_pmf(n_i,logistic_grid,cdf = cdf_logistic)$L1_pmf
  res_logis_pmf_n_32[i] = L1_L2_pmf(n_i,logistic_grid,cdf = cdf_logistic)$L2_pmf
  
  # Logistic Distribution - n4
  
  x = sort(rlogis(n[4]))
  n_i <- hist(x, breaks = logistic_grid, plot = F)$counts
  
  res_logis_pmf_n_41[i] = L1_L2_pmf(n_i,logistic_grid,cdf = cdf_logistic)$L1_pmf
  res_logis_pmf_n_42[i] = L1_L2_pmf(n_i,logistic_grid,cdf = cdf_logistic)$L2_pmf
  
  # Logistic Distribution - n5
  
  x = sort(rlogis(n[5]))
  n_i <- hist(x, breaks = logistic_grid, plot = F)$counts
  
  res_logis_pmf_n_51[i] = L1_L2_pmf(n_i,logistic_grid,cdf = cdf_logistic)$L1_pmf
  res_logis_pmf_n_52[i] = L1_L2_pmf(n_i,logistic_grid,cdf = cdf_logistic)$L2_pmf
  
}

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L1_n/L1-PMF-1-Logistic-n-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pmf_n_11,res_logis_pmf_n_21,res_logis_pmf_n_31,res_logis_pmf_n_41,res_logis_pmf_n_51,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors4,ylim=c(0,1.2))
segments(x0 = 0.29, y0 =0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L2_n/L2-PMF-1-Logistic-n-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pmf_n_12,res_logis_pmf_n_22,res_logis_pmf_n_32,res_logis_pmf_n_42,res_logis_pmf_n_52,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors4,ylim=c(0,0.8))
segments(x0 = 0.29, y0 =0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

# Student's t Distribution - Sample Size Simulation 

MM <- ceiling(qt(0.99999999, df))+100
t_grid <- seq(-MM,MM, by = t_deltas[1])

res_t_pmf_n_11 = rep(0,B)
res_t_pmf_n_21 = rep(0,B)
res_t_pmf_n_31 = rep(0,B)
res_t_pmf_n_41 = rep(0,B)
res_t_pmf_n_51 = rep(0,B)

res_t_pmf_n_12 = rep(0,B)
res_t_pmf_n_22 = rep(0,B)
res_t_pmf_n_32 = rep(0,B)
res_t_pmf_n_42 = rep(0,B)
res_t_pmf_n_52 = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  # t Distribution - n1
  
  x = sort(rt(n[1],df = df))
  n_i <- hist(x, breaks = t_grid, plot = F)$counts
  
  res_t_pmf_n_11[i] = L1_L2_pmf(n_i,t_grid,cdf = cdf_t)$L1_pmf
  res_t_pmf_n_12[i] = L1_L2_pmf(n_i,t_grid,cdf = cdf_t)$L2_pmf
  
  # t Distribution - n2
  
  x = sort(rt(n[2],df = df))
  n_i <- hist(x, breaks = t_grid, plot = F)$counts
  
  res_t_pmf_n_21[i] = L1_L2_pmf(n_i,t_grid,cdf = cdf_t)$L1_pmf
  res_t_pmf_n_22[i] = L1_L2_pmf(n_i,t_grid,cdf = cdf_t)$L2_pmf
  
  # t Distribution - n3
  
  x = sort(rt(n[3],df = df))
  n_i <- hist(x, breaks = t_grid, plot = F)$counts
  
  res_t_pmf_n_31[i] = L1_L2_pmf(n_i,t_grid,cdf = cdf_t)$L1_pmf
  res_t_pmf_n_32[i] = L1_L2_pmf(n_i,t_grid,cdf = cdf_t)$L2_pmf
  
  # t Distribution - n4
  
  x = sort(rt(n[4],df = df))
  n_i <- hist(x, breaks = t_grid, plot = F)$counts
  
  res_t_pmf_n_41[i] = L1_L2_pmf(n_i,t_grid,cdf = cdf_t)$L1_pmf
  res_t_pmf_n_42[i] = L1_L2_pmf(n_i,t_grid,cdf = cdf_t)$L2_pmf
  
  # t Distribution - n5
  
  x = sort(rt(n[5],df = df))
  n_i <- hist(x, breaks = t_grid, plot = F)$counts
  
  res_t_pmf_n_51[i] = L1_L2_pmf(n_i,t_grid,cdf = cdf_t)$L1_pmf
  res_t_pmf_n_52[i] = L1_L2_pmf(n_i,t_grid,cdf = cdf_t)$L2_pmf
  
}

t_L1_pmf_limit = 0
t_L2_pmf_limit = 0

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L1_n/L1-PMF-1-t-n-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pmf_n_11,res_t_pmf_n_21,res_t_pmf_n_31,res_t_pmf_n_41,res_t_pmf_n_51,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors5,ylim=c(0,0.9))
segments(x0 = 0.29, y0 = median(res_t_pmf_n_51), x1 = 10^5, y1 = t_L1_pmf_limit, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L2_n/L2-PMF-1-t-n-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pmf_n_12,res_t_pmf_n_22,res_t_pmf_n_32,res_t_pmf_n_42,res_t_pmf_n_52,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors5,ylim=c(0,0.8))
segments(x0 = 0.29, y0 = median(res_t_pmf_n_52), x1 = 10^5, y1 = t_L2_pmf_limit, col = "red", lty = 2)
dev.off()

# Normal Distribution - Grid-Width Simulation 

n = c(10^3)
B = 100

MM = qnorm(0.9999999,mu,sigma) + 3
norm_grids = list(seq(-MM,MM, by = norm_deltas[5]),
                  seq(-MM,MM, by = norm_deltas[4]),
                  seq(-MM,MM, by = norm_deltas[3]),
                  seq(-MM,MM, by = norm_deltas[2]),
                  seq(-MM,MM, by = norm_deltas[1]))

res_norm_pmf_gw_11 = rep(0,B)
res_norm_pmf_gw_21 = rep(0,B)
res_norm_pmf_gw_31 = rep(0,B)
res_norm_pmf_gw_41 = rep(0,B)
res_norm_pmf_gw_51 = rep(0,B)

res_norm_pmf_gw_12 = rep(0,B)
res_norm_pmf_gw_22 = rep(0,B)
res_norm_pmf_gw_32 = rep(0,B)
res_norm_pmf_gw_42 = rep(0,B)
res_norm_pmf_gw_52 = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  x <- sort(rnorm(n, mean = mu,sd=sigma))
  
  # Normal Distribution - GW1

  n_i <- hist(x, breaks = norm_grids[[1]], plot = F)$counts
  
  res_norm_pmf_gw_11[i] = L1_L2_pmf(n_i,norm_grids[[1]],cdf = cdf_norm)$L1_pmf
  res_norm_pmf_gw_12[i] = L1_L2_pmf(n_i,norm_grids[[1]],cdf = cdf_norm)$L2_pmf
  
  # Normal Distribution - GW2

  n_i <- hist(x, breaks = norm_grids[[2]], plot = F)$counts
  
  res_norm_pmf_gw_21[i] = L1_L2_pmf(n_i,norm_grids[[2]],cdf = cdf_norm)$L1_pmf
  res_norm_pmf_gw_22[i] = L1_L2_pmf(n_i,norm_grids[[2]],cdf = cdf_norm)$L2_pmf
  
  # Normal Distribution - GW3

  n_i <- hist(x, breaks = norm_grids[[3]], plot = F)$counts
  
  res_norm_pmf_gw_31[i] = L1_L2_pmf(n_i,norm_grids[[3]],cdf = cdf_norm)$L1_pmf
  res_norm_pmf_gw_32[i] = L1_L2_pmf(n_i,norm_grids[[3]],cdf = cdf_norm)$L2_pmf
  
  # Normal Distribution - GW4

  n_i <- hist(x, breaks = norm_grids[[4]], plot = F)$counts
  
  res_norm_pmf_gw_41[i] = L1_L2_pmf(n_i,norm_grids[[4]],cdf = cdf_norm)$L1_pmf
  res_norm_pmf_gw_42[i] = L1_L2_pmf(n_i,norm_grids[[4]],cdf = cdf_norm)$L2_pmf
  
  # Normal Distribution - GW5

  n_i <- hist(x, breaks = norm_grids[[5]], plot = F)$counts
  
  res_norm_pmf_gw_51[i] = L1_L2_pmf(n_i,norm_grids[[5]],cdf = cdf_norm)$L1_pmf
  res_norm_pmf_gw_52[i] = L1_L2_pmf(n_i,norm_grids[[5]],cdf = cdf_norm)$L2_pmf
  
}

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L1_gw/L1-PMF-1-norm-gw-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_norm_pmf_gw_11,res_norm_pmf_gw_21,res_norm_pmf_gw_31,res_norm_pmf_gw_41,res_norm_pmf_gw_51,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,0.15))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L2_gw/L2-PMF-1-norm-gw-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_norm_pmf_gw_12,res_norm_pmf_gw_22,res_norm_pmf_gw_32,res_norm_pmf_gw_42,res_norm_pmf_gw_52,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,0.06))
dev.off()

# Beta Distribution - Grid Width Simulation 

beta_grids = list(seq(0,1.2, by = beta_deltas[5]),
                  seq(0,1.3, by = beta_deltas[4]),
                  seq(0,1.2, by = beta_deltas[3]),
                  seq(0,1.2, by = beta_deltas[2]),
                  seq(0,1.1, by = beta_deltas[1]))

res_beta_pmf_gw_11 = rep(0,B)
res_beta_pmf_gw_21 = rep(0,B)
res_beta_pmf_gw_31 = rep(0,B)
res_beta_pmf_gw_41 = rep(0,B)
res_beta_pmf_gw_51 = rep(0,B)

res_beta_pmf_gw_12 = rep(0,B)
res_beta_pmf_gw_22 = rep(0,B)
res_beta_pmf_gw_32 = rep(0,B)
res_beta_pmf_gw_42 = rep(0,B)
res_beta_pmf_gw_52 = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  x <- sort(rbeta(n,alpha,beta))
  
  # Beta Distribution - GW1
  
  n_i <- hist(x, breaks = beta_grids[[1]], plot = F)$counts
  
  res_beta_pmf_gw_11[i] = L1_L2_pmf(n_i,beta_grids[[1]],cdf = cdf_beta)$L1_pmf
  res_beta_pmf_gw_12[i] = L1_L2_pmf(n_i,beta_grids[[1]],cdf = cdf_beta)$L2_pmf
  
  # Beta Distribution - GW2

  n_i <- hist(x, breaks = beta_grids[[2]], plot = F)$counts
  
  res_beta_pmf_gw_21[i] = L1_L2_pmf(n_i,beta_grids[[2]],cdf = cdf_beta)$L1_pmf
  res_beta_pmf_gw_22[i] = L1_L2_pmf(n_i,beta_grids[[2]],cdf = cdf_beta)$L2_pmf
  
  # Beta Distribution - GW3

  n_i <- hist(x, breaks = beta_grids[[3]], plot = F)$counts
  
  res_beta_pmf_gw_31[i] = L1_L2_pmf(n_i,beta_grids[[3]],cdf = cdf_beta)$L1_pmf
  res_beta_pmf_gw_32[i] = L1_L2_pmf(n_i,beta_grids[[3]],cdf = cdf_beta)$L2_pmf
  
  # Beta Distribution - GW4

  n_i <- hist(x, breaks = beta_grids[[4]], plot = F)$counts
  
  res_beta_pmf_gw_41[i] = L1_L2_pmf(n_i,beta_grids[[4]],cdf = cdf_beta)$L1_pmf
  res_beta_pmf_gw_42[i] = L1_L2_pmf(n_i,beta_grids[[4]],cdf = cdf_beta)$L2_pmf
  
  # Beta Distribution - GW5

  n_i <- hist(x, breaks = beta_grids[[5]], plot = F)$counts
  
  res_beta_pmf_gw_51[i] = L1_L2_pmf(n_i,beta_grids[[5]],cdf = cdf_beta)$L1_pmf
  res_beta_pmf_gw_52[i] = L1_L2_pmf(n_i,beta_grids[[5]],cdf = cdf_beta)$L2_pmf
  
}

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L1_gw/L1-PMF-1-Beta-gw-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pmf_gw_11,res_beta_pmf_gw_21,res_beta_pmf_gw_31,res_beta_pmf_gw_41,res_beta_pmf_gw_51,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.1))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L2_gw/L2-PMF-1-Beta-gw-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pmf_gw_12,res_beta_pmf_gw_22,res_beta_pmf_gw_32,res_beta_pmf_gw_42,res_beta_pmf_gw_52,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.06))
dev.off()

# Gamma Distribution - Grid Width Simulation 

MM = qgamma(0.9999999,alpha,beta) + 3
gamma_grids = list(seq(0,MM, by = gamma_deltas[5]),
                   seq(0,MM, by = gamma_deltas[4]),
                   seq(0,MM, by = gamma_deltas[3]),
                   seq(0,MM, by = gamma_deltas[2]),
                   seq(0,MM, by = gamma_deltas[1]))

res_gamma_pmf_gw_11 = rep(0,B)
res_gamma_pmf_gw_21 = rep(0,B)
res_gamma_pmf_gw_31 = rep(0,B)
res_gamma_pmf_gw_41 = rep(0,B)
res_gamma_pmf_gw_51 = rep(0,B)

res_gamma_pmf_gw_12 = rep(0,B)
res_gamma_pmf_gw_22 = rep(0,B)
res_gamma_pmf_gw_32 = rep(0,B)
res_gamma_pmf_gw_42 = rep(0,B)
res_gamma_pmf_gw_52 = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  x  <- sort(rgamma(n,alpha,beta))
  
  # Gamma Distribution - GW1

  n_i <- hist(x, breaks = gamma_grids[[1]], plot = F)$counts
  
  res_gamma_pmf_gw_11[i] = L1_L2_pmf(n_i,gamma_grids[[1]],cdf = cdf_gamma)$L1_pmf
  res_gamma_pmf_gw_12[i] = L1_L2_pmf(n_i,gamma_grids[[1]],cdf = cdf_gamma)$L2_pmf
  
  # Gamma Distribution - GW2

  n_i <- hist(x, breaks = gamma_grids[[2]], plot = F)$counts
  
  res_gamma_pmf_gw_21[i] = L1_L2_pmf(n_i,gamma_grids[[2]],cdf = cdf_gamma)$L1_pmf
  res_gamma_pmf_gw_22[i] = L1_L2_pmf(n_i,gamma_grids[[2]],cdf = cdf_gamma)$L2_pmf
  
  # Gamma Distribution - GW3

  n_i <- hist(x, breaks = gamma_grids[[3]], plot = F)$counts
  
  res_gamma_pmf_gw_31[i] = L1_L2_pmf(n_i,gamma_grids[[3]],cdf = cdf_gamma)$L1_pmf
  res_gamma_pmf_gw_32[i] = L1_L2_pmf(n_i,gamma_grids[[3]],cdf = cdf_gamma)$L2_pmf
  
  # Gamma Distribution - GW4

  n_i <- hist(x, breaks = gamma_grids[[4]], plot = F)$counts
  
  res_gamma_pmf_gw_41[i] = L1_L2_pmf(n_i,gamma_grids[[4]],cdf = cdf_gamma)$L1_pmf
  res_gamma_pmf_gw_42[i] = L1_L2_pmf(n_i,gamma_grids[[4]],cdf = cdf_gamma)$L2_pmf
  
  # Gamma Distribution - GW5

  n_i <- hist(x, breaks = gamma_grids[[5]], plot = F)$counts
  
  res_gamma_pmf_gw_51[i] = L1_L2_pmf(n_i,gamma_grids[[5]],cdf = cdf_gamma)$L1_pmf
  res_gamma_pmf_gw_52[i] = L1_L2_pmf(n_i,gamma_grids[[5]],cdf = cdf_gamma)$L2_pmf
  
}

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L1_gw/L1-PMF-1-Gamma-gw-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pmf_gw_11,res_gamma_pmf_gw_21,res_gamma_pmf_gw_31,res_gamma_pmf_gw_41,res_gamma_pmf_gw_51,
        main = "",
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,0.1))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L2_gw/L2-PMF-1-Gamma-gw-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pmf_gw_12,res_gamma_pmf_gw_22,res_gamma_pmf_gw_32,res_gamma_pmf_gw_42,res_gamma_pmf_gw_52,
        main = "",
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,0.06))
dev.off()

# Logistic Distribution - Grid Width Simulation 

MM = qlogis(0.9999999,mu,sigma) + 3
logistic_grids = list(seq(-MM,MM, by = logistic_deltas[5]),
                      seq(-MM,MM, by = logistic_deltas[4]),
                      seq(-MM,MM, by = logistic_deltas[3]),
                      seq(-MM,MM, by = logistic_deltas[2]),
                      seq(-MM,MM, by = logistic_deltas[1]))

res_logis_pmf_gw_11 = rep(0,B)
res_logis_pmf_gw_21 = rep(0,B)
res_logis_pmf_gw_31 = rep(0,B)
res_logis_pmf_gw_41 = rep(0,B)
res_logis_pmf_gw_51 = rep(0,B)

res_logis_pmf_gw_12 = rep(0,B)
res_logis_pmf_gw_22 = rep(0,B)
res_logis_pmf_gw_32 = rep(0,B)
res_logis_pmf_gw_42 = rep(0,B)
res_logis_pmf_gw_52 = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  # Logistic Distribution - GW1
  
  x = sort(rlogis(n))
  n_i <- hist(x, breaks = logistic_grids[[1]], plot = F)$counts
  
  res_logis_pmf_gw_11[i] = L1_L2_pmf(n_i,logistic_grids[[1]],cdf = cdf_logistic)$L1_pmf
  res_logis_pmf_gw_12[i] = L1_L2_pmf(n_i,logistic_grids[[1]],cdf = cdf_logistic)$L2_pmf
  
  # Logistic Distribution - GW2
  
  x = sort(rlogis(n))
  n_i <- hist(x, breaks = logistic_grids[[2]], plot = F)$counts
  
  res_logis_pmf_gw_21[i] = L1_L2_pmf(n_i,logistic_grids[[2]],cdf = cdf_logistic)$L1_pmf
  res_logis_pmf_gw_22[i] = L1_L2_pmf(n_i,logistic_grids[[2]],cdf = cdf_logistic)$L2_pmf
  
  # Logistic Distribution - GW3
  
  x = sort(rlogis(n))
  n_i <- hist(x, breaks = logistic_grids[[3]], plot = F)$counts
  
  res_logis_pmf_gw_31[i] = L1_L2_pmf(n_i,logistic_grids[[3]],cdf = cdf_logistic)$L1_pmf
  res_logis_pmf_gw_32[i] = L1_L2_pmf(n_i,logistic_grids[[3]],cdf = cdf_logistic)$L2_pmf
  
  # Logistic Distribution - GW4
  
  x = sort(rlogis(n))
  n_i <- hist(x, breaks = logistic_grids[[4]], plot = F)$counts
  
  res_logis_pmf_gw_41[i] = L1_L2_pmf(n_i,logistic_grids[[4]],cdf = cdf_logistic)$L1_pmf
  res_logis_pmf_gw_42[i] = L1_L2_pmf(n_i,logistic_grids[[4]],cdf = cdf_logistic)$L2_pmf
  
  # Logistic Distribution - GW5
  
  x = sort(rlogis(n))
  n_i <- hist(x, breaks = logistic_grids[[5]], plot = F)$counts
  
  res_logis_pmf_gw_51[i] = L1_L2_pmf(n_i,logistic_grids[[5]],cdf = cdf_logistic)$L1_pmf
  res_logis_pmf_gw_52[i] = L1_L2_pmf(n_i,logistic_grids[[5]],cdf = cdf_logistic)$L2_pmf
  
}

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L1_gw/L1-PMF-1-Logistic-gw-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pmf_gw_11,res_logis_pmf_gw_21,res_logis_pmf_gw_31,res_logis_pmf_gw_41,res_logis_pmf_gw_51,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,0.1))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L2_gw/L2-PMF-1-Logistic-gw-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pmf_gw_12,res_logis_pmf_gw_22,res_logis_pmf_gw_32,res_logis_pmf_gw_42,res_logis_pmf_gw_52,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,0.06))
dev.off()

# Student's t Distribution - Grid Width Simulation 

MM = qt(0.9999999,df) + 3
t_grids = list(seq(-MM,MM, by = t_deltas[5]),
               seq(-MM,MM, by = t_deltas[4]),
               seq(-MM,MM, by = t_deltas[3]),
               seq(-MM,MM, by = t_deltas[2]),
               seq(-MM,MM, by = t_deltas[1]))

res_t_pmf_gw_11 = rep(0,B)
res_t_pmf_gw_21 = rep(0,B)
res_t_pmf_gw_31 = rep(0,B)
res_t_pmf_gw_41 = rep(0,B)
res_t_pmf_gw_51 = rep(0,B)

res_t_pmf_gw_12 = rep(0,B)
res_t_pmf_gw_22 = rep(0,B)
res_t_pmf_gw_32 = rep(0,B)
res_t_pmf_gw_42 = rep(0,B)
res_t_pmf_gw_52 = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  # t Distribution - GW1
  
  x = sort(rt(n,df = df))
  n_i <- hist(x, breaks = t_grids[[1]], plot = F)$counts
  
  res_t_pmf_gw_11[i] = L1_L2_pmf(n_i,t_grids[[1]],cdf = cdf_t)$L1_pmf
  res_t_pmf_gw_12[i] = L1_L2_pmf(n_i,t_grids[[1]],cdf = cdf_t)$L2_pmf
  
  # t Distribution - GW2
  
  x = sort(rt(n,df = df))
  n_i <- hist(x, breaks = t_grids[[2]], plot = F)$counts
  
  res_t_pmf_gw_21[i] = L1_L2_pmf(n_i,t_grids[[2]],cdf = cdf_t)$L1_pmf
  res_t_pmf_gw_22[i] = L1_L2_pmf(n_i,t_grids[[2]],cdf = cdf_t)$L2_pmf
  
  # t Distribution - GW3
  
  x = sort(rt(n,df = df))
  n_i <- hist(x, breaks = t_grids[[3]], plot = F)$counts
  
  res_t_pmf_gw_31[i] = L1_L2_pmf(n_i,t_grids[[3]],cdf = cdf_t)$L1_pmf
  res_t_pmf_gw_32[i] = L1_L2_pmf(n_i,t_grids[[3]],cdf = cdf_t)$L2_pmf
  
  # t Distribution - GW4
  
  x = sort(rt(n,df = df))
  n_i <- hist(x, breaks = t_grids[[4]], plot = F)$counts
  
  res_t_pmf_gw_41[i] = L1_L2_pmf(n_i,t_grids[[4]],cdf = cdf_t)$L1_pmf
  res_t_pmf_gw_42[i] = L1_L2_pmf(n_i,t_grids[[4]],cdf = cdf_t)$L2_pmf
  
  # t Distribution - GW5
  
  x = sort(rt(n,df = df))
  n_i <- hist(x, breaks = t_grids[[5]], plot = F)$counts
  
  res_t_pmf_gw_51[i] = L1_L2_pmf(n_i,t_grids[[5]],cdf = cdf_t)$L1_pmf
  res_t_pmf_gw_52[i] = L1_L2_pmf(n_i,t_grids[[5]],cdf = cdf_t)$L1_pmf
  
}

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L1_gw/L1-PMF-1-t-gw-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pmf_gw_11,res_t_pmf_gw_21,res_t_pmf_gw_31,res_t_pmf_gw_41,res_t_pmf_gw_51,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,0.1))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L2_gw/L2-PMF-1-t-gw-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pmf_gw_12,res_t_pmf_gw_22,res_t_pmf_gw_32,res_t_pmf_gw_42,res_t_pmf_gw_52,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,0.15))
dev.off()
