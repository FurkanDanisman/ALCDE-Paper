# All distributions together n = 1000, n = 10000

# Each Distribution, Sample size increasing, grid width increasing 

# Simulation Study 

# CDFs

cdf_laplace <- function(x) {
  plaplace(x,mean = mu,sd = sigma)
}

cdf_chi <- function(x) {
  pchisq(x,df)
}

cdf_lnorm <- function(x) {
  plnorm(x,meanlog = mul,sdlog = sigma)
}

cdf_weibull <- function(x) {
  pweibull(x,shape)
}

cdf_pareto <- function(x) {
  ppareto(x,location = location,shape = shape_p)
}


n <- 1000
B = 100

res_pmf_laplace_L1 = rep(0,B)
res_pmf_chisq_L1 = rep(0,B)
res_pmf_lnorm_L1 = rep(0,B)
res_pmf_wei_L1 = rep(0,B)
res_pmf_pareto_L1 = rep(0,B)

res_pmf_laplace_L2 = rep(0,B)
res_pmf_chisq_L2 = rep(0,B)
res_pmf_lnorm_L2 = rep(0,B)
res_pmf_wei_L2 = rep(0,B)
res_pmf_pareto_L2 = rep(0,B)


library(jmuOutlier)
library(EnvStats)

# Parameters

shape=2;location=1;
mu = 10;sigma=1;
mul = 1;
df=3;
shape_p = 4;

# Standardized Delta Values # 

laplace_deltas  = c(0.7071, 1.4142, 2.1213, 2.8284, 3.5355)
chisqr_deltas  = c(1.2247, 2.4495, 3.6742, 4.8990, 6.1237)
lnorm_deltas = c(2.9374, 5.8747, 8.8121, 11.7495,14.6869)
weibull_deltas = c(0.2316, 0.4633, 0.6949, 0.9265, 1.1581)
pareto_deltas = c(0.4714, 0.9428, 1.4142, 1.8856, 2.3570)

# Corresponding k values to Standardized Delta Values # 

laplace_ks  = c(36, 18, 12, 9, 7)
chisqr_ks  = c(22, 11, 7, 5, 4)
lnorm_ks = c(40, 20, 13, 10, 8)
weibull_ks = c(38, 19, 12, 9,  7)
pareto_ks = c(50, 25, 16, 12, 10)

# Simulations Within the Distributions 

c11 = darken("#FDB462", factor = 1.1)
c21 = darken(c11, factor = 1.1)
c31 = darken(c21, factor = 1.1)
c41 = darken(c31, factor = 1.1)
colors1 = c("#FDB462", c11, c21, c31, c41)

# Colors for  Distribution

c11 = darken("#B3DE69", factor = 1.1)
c21 = darken(c11, factor = 1.1)
c31 = darken(c21, factor = 1.1)
c41 = darken(c31, factor = 1.1)
colors2 = c("#B3DE69", c11, c21, c31, c41)

# Colors for Beta Distribution

c11 = darken("#BC80BD", factor = 1.1)
c21 = darken(c11, factor = 1.1)
c31 = darken(c21, factor = 1.1)
c41 = darken(c31, factor = 1.1)
colors3 = c("#BC80BD", c11, c21, c31, c41)

# Colors for Bernoulli Distribution

c11 = darken("cyan", factor = 1.1)
c21 = darken(c11, factor = 1.1)
c31 = darken(c21, factor = 1.1)
c41 = darken(c31, factor = 1.1)
colors4 = c("cyan", c11, c21, c31, c41)

# Colors for Poisson Distribution

c11 = darken("#FB8072", factor = 1.1)
c21 = darken(c11, factor = 1.1)
c31 = darken(c21, factor = 1.1)
c41 = darken(c31, factor = 1.1)
colors5 = c("#FB8072", c11, c21, c31, c41)

n = c(10^2,10^3,10^4,10^5,10^6)
B = 100

# Laplace Distribution - Sample Size Simulation 

MM = qlaplace(0.9999999,mu,sigma) + 3
laplace_grid <- seq(-MM,MM, by = laplace_deltas[1])

res_lap_pmf_n_11 = rep(0,B)
res_lap_pmf_n_21 = rep(0,B)
res_lap_pmf_n_31 = rep(0,B)
res_lap_pmf_n_41 = rep(0,B)
res_lap_pmf_n_51 = rep(0,B)

res_lap_pmf_n_12 = rep(0,B)
res_lap_pmf_n_22 = rep(0,B)
res_lap_pmf_n_32 = rep(0,B)
res_lap_pmf_n_42 = rep(0,B)
res_lap_pmf_n_52 = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  # Laplace Distribution - n1
  
  x <- sort(rlaplace(n[1], mean = mu,sd=sigma))
  n_i <- hist(x, breaks = laplace_grid, plot = F)$counts
  
  res_lap_pmf_n_11[i] = L1_L2_pmf(n_i,laplace_grid,cdf = cdf_laplace)$L1_pmf
  res_lap_pmf_n_12[i] = L1_L2_pmf(n_i,laplace_grid,cdf = cdf_laplace)$L2_pmf
  
  # Laplace Distribution - n2
  
  x <- sort(rlaplace(n[2], mean = mu,sd=sigma))
  n_i <- hist(x, breaks = laplace_grid, plot = F)$counts
  
  res_lap_pmf_n_21[i] = L1_L2_pmf(n_i,laplace_grid,cdf = cdf_laplace)$L1_pmf
  res_lap_pmf_n_22[i] = L1_L2_pmf(n_i,laplace_grid,cdf = cdf_laplace)$L2_pmf
  
  # Laplace Distribution - n3
  
  x <- sort(rlaplace(n[3], mean = mu,sd=sigma))
  n_i <- hist(x, breaks = laplace_grid, plot = F)$counts
  
  res_lap_pmf_n_31[i] = L1_L2_pmf(n_i,laplace_grid,cdf = cdf_laplace)$L1_pmf
  res_lap_pmf_n_32[i] = L1_L2_pmf(n_i,laplace_grid,cdf = cdf_laplace)$L2_pmf
  
  # Laplace Distribution - n4
  
  x <- sort(rlaplace(n[4], mean = mu,sd=sigma))
  n_i <- hist(x, breaks = laplace_grid, plot = F)$counts
  
  res_lap_pmf_n_41[i] = L1_L2_pmf(n_i,laplace_grid,cdf = cdf_laplace)$L1_pmf
  res_lap_pmf_n_42[i] = L1_L2_pmf(n_i,laplace_grid,cdf = cdf_laplace)$L2_pmf
  
  # Laplace Distribution - n5
  
  x <- sort(rlaplace(n[5], mean = mu,sd=sigma))
  n_i <- hist(x, breaks = laplace_grid, plot = F)$counts
  
  res_lap_pmf_n_51[i] = L1_L2_pmf(n_i,laplace_grid,cdf = cdf_laplace)$L1_pmf
  res_lap_pmf_n_52[i] = L1_L2_pmf(n_i,laplace_grid,cdf = cdf_laplace)$L2_pmf
  
}

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L1_n/L1-PMF-2-Laplace-n-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_lap_pmf_n_11,res_lap_pmf_n_21,res_lap_pmf_n_31,res_lap_pmf_n_41,res_lap_pmf_n_51,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors1,ylim=c(0,0.4))
segments(x0 = 0.29, y0 =0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L2_n/L2-PMF-2-Laplace-n-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_lap_pmf_n_12,res_lap_pmf_n_22,res_lap_pmf_n_32,res_lap_pmf_n_42,res_lap_pmf_n_52,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors1,ylim=c(0,0.15))
segments(x0 = 0.29, y0 =0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

# Chi-Square Distribution - Sample Size Simulation 

MM = qchisq(0.99999999,df) + 20
chisq_grid <- seq(0,MM, by = chisqr_deltas[1])

res_chisq_pmf_n_11 = rep(0,B)
res_chisq_pmf_n_21 = rep(0,B)
res_chisq_pmf_n_31 = rep(0,B)
res_chisq_pmf_n_41 = rep(0,B)
res_chisq_pmf_n_51 = rep(0,B)

res_chisq_pmf_n_12 = rep(0,B)
res_chisq_pmf_n_22 = rep(0,B)
res_chisq_pmf_n_32 = rep(0,B)
res_chisq_pmf_n_42 = rep(0,B)
res_chisq_pmf_n_52 = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  # Chi-square Distribution - n1
  
  x   <- sort(rchisq(n[1],df))
  n_i <- hist(x, breaks = chisq_grid, plot = F)$counts
  
  res_chisq_pmf_n_11[i] = L1_L2_pmf(n_i,chisq_grid,cdf = cdf_chi)$L1_pmf
  res_chisq_pmf_n_12[i] = L1_L2_pmf(n_i,chisq_grid,cdf = cdf_chi)$L2_pmf
  
  # Chi-square Distribution - n2
  
  x   <- sort(rchisq(n[2],df))
  n_i <- hist(x, breaks = chisq_grid, plot = F)$counts
  
  res_chisq_pmf_n_21[i] = L1_L2_pmf(n_i,chisq_grid,cdf = cdf_chi)$L1_pmf
  res_chisq_pmf_n_22[i] = L1_L2_pmf(n_i,chisq_grid,cdf = cdf_chi)$L2_pmf
  
  # Chi-square Distribution - n3
  
  x   <- sort(rchisq(n[3],df))
  n_i <- hist(x, breaks = chisq_grid, plot = F)$counts
  
  res_chisq_pmf_n_31[i] = L1_L2_pmf(n_i,chisq_grid,cdf = cdf_chi)$L1_pmf
  res_chisq_pmf_n_32[i] = L1_L2_pmf(n_i,chisq_grid,cdf = cdf_chi)$L2_pmf
  
  # Chi-square Distribution - n4
  
  x   <- sort(rchisq(n[4],df))
  n_i <- hist(x, breaks = chisq_grid, plot = F)$counts
  
  res_chisq_pmf_n_41[i] = L1_L2_pmf(n_i,chisq_grid,cdf = cdf_chi)$L1_pmf
  res_chisq_pmf_n_42[i] = L1_L2_pmf(n_i,chisq_grid,cdf = cdf_chi)$L2_pmf
  
  # Chi-Square Distribution - n5
  
  x   <- sort(rchisq(n[5],df))
  n_i <- hist(x, breaks = chisq_grid, plot = F)$counts
  
  res_chisq_pmf_n_51[i] = L1_L2_pmf(n_i,chisq_grid,cdf = cdf_chi)$L1_pmf
  res_chisq_pmf_n_52[i] = L1_L2_pmf(n_i,chisq_grid,cdf = cdf_chi)$L2_pmf
  
}

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L1_n/L1-PMF-2-Chisq-n-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_chisq_pmf_n_11,res_chisq_pmf_n_21,res_chisq_pmf_n_31,res_chisq_pmf_n_41,res_chisq_pmf_n_51,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2,ylim=c(0,0.3))
segments(x0 = 0.29, y0 =0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L2_n/L2-PMF-2-Chisq-n-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_chisq_pmf_n_12,res_chisq_pmf_n_22,res_chisq_pmf_n_32,res_chisq_pmf_n_42,res_chisq_pmf_n_52,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2,ylim=c(0,0.15))
segments(x0 = 0.29, y0 =0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

# Log-Normal Distribution - Sample Size Simulation 

MM = qlnorm(0.9999999,mu,sigma) + 3
lnrom_grid <- seq(0,MM, by = lnorm_deltas[1])

res_logn_pmf_n_11 = rep(0,B)
res_logn_pmf_n_21 = rep(0,B)
res_logn_pmf_n_31 = rep(0,B)
res_logn_pmf_n_41 = rep(0,B)
res_logn_pmf_n_51 = rep(0,B)

res_logn_pmf_n_12 = rep(0,B)
res_logn_pmf_n_22 = rep(0,B)
res_logn_pmf_n_32 = rep(0,B)
res_logn_pmf_n_42 = rep(0,B)
res_logn_pmf_n_52 = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  # Log-Normal Distribution - n1
  
  x   <- sort(rlnorm(n[1],mul,sigma))
  n_i <- hist(x, breaks = lnrom_grid, plot = F)$counts
  
  res_logn_pmf_n_11[i] = L1_L2_pmf(n_i,lnrom_grid,cdf = cdf_lnorm)$L1_pmf
  res_logn_pmf_n_12[i] = L1_L2_pmf(n_i,lnrom_grid,cdf = cdf_lnorm)$L2_pmf
  
  # Log-Normal Distribution - n2
  
  x   <- sort(rlnorm(n[2],mul,sigma))
  n_i <- hist(x, breaks = lnrom_grid, plot = F)$counts
  
  res_logn_pmf_n_21[i] = L1_L2_pmf(n_i,lnrom_grid,cdf = cdf_lnorm)$L1_pmf
  res_logn_pmf_n_22[i] = L1_L2_pmf(n_i,lnrom_grid,cdf = cdf_lnorm)$L2_pmf
  
  # Log-Normal Distribution - n3
  
  x   <- sort(rlnorm(n[3],mul,sigma))
  n_i <- hist(x, breaks = lnrom_grid, plot = F)$counts
  
  res_logn_pmf_n_31[i] = L1_L2_pmf(n_i,lnrom_grid,cdf = cdf_lnorm)$L1_pmf
  res_logn_pmf_n_32[i] = L1_L2_pmf(n_i,lnrom_grid,cdf = cdf_lnorm)$L2_pmf
  
  # Log-Normal Distribution - n4
  
  x   <- sort(rlnorm(n[4],mul,sigma))
  n_i <- hist(x, breaks = lnrom_grid, plot = F)$counts
  
  res_logn_pmf_n_41[i] = L1_L2_pmf(n_i,lnrom_grid,cdf = cdf_lnorm)$L1_pmf
  res_logn_pmf_n_42[i] = L1_L2_pmf(n_i,lnrom_grid,cdf = cdf_lnorm)$L2_pmf
  
  # Log-Normal Distribution - n5
  
  x   <- sort(rlnorm(n[5],mul,sigma))
  n_i <- hist(x, breaks = lnrom_grid, plot = F)$counts
  
  res_logn_pmf_n_51[i] = L1_L2_pmf(n_i,lnrom_grid,cdf = cdf_lnorm)$L1_pmf
  res_logn_pmf_n_52[i] = L1_L2_pmf(n_i,lnrom_grid,cdf = cdf_lnorm)$L2_pmf
  
}

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L1_n/L1-PMF-2-Lnorm-n-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logn_pmf_n_11,res_logn_pmf_n_21,res_logn_pmf_n_31,res_logn_pmf_n_41,res_logn_pmf_n_51,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3,ylim=c(0,0.3))
segments(x0 = 0.29, y0 = median(res_logn_pmf_n_51), x1 = 10^5, y1 = mean(res_logn_pmf_n_51), col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L2_n/L2-PMF-2-Lnorm-n-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logn_pmf_n_12,res_logn_pmf_n_22,res_logn_pmf_n_32,res_logn_pmf_n_42,res_logn_pmf_n_52,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3,ylim=c(0,0.15))
segments(x0 = 0.29, y0 = median(res_logn_pmf_n_52), x1 = 10^5, y1 = mean(res_logn_pmf_n_52), col = "red", lty = 2)
dev.off()

# Weibull Distribution - Sample Size Simulation 

MM = qweibull(0.9999999,shape) + 3
weibull_grid <- seq(-MM,MM, by = weibull_deltas[1])

res_wei_pmf_n_11 = rep(0,B)
res_wei_pmf_n_21 = rep(0,B)
res_wei_pmf_n_31 = rep(0,B)
res_wei_pmf_n_41 = rep(0,B)
res_wei_pmf_n_51 = rep(0,B)

res_wei_pmf_n_12 = rep(0,B)
res_wei_pmf_n_22 = rep(0,B)
res_wei_pmf_n_32 = rep(0,B)
res_wei_pmf_n_42 = rep(0,B)
res_wei_pmf_n_52 = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  # Weibull Distribution - n1
  
  x = sort(rweibull(n[1],shape))
  n_i <- hist(x, breaks = weibull_grid, plot = F)$counts
  
  res_wei_pmf_n_11[i] = L1_L2_pmf(n_i,weibull_grid,cdf = cdf_weibull)$L1_pmf
  res_wei_pmf_n_12[i] = L1_L2_pmf(n_i,weibull_grid,cdf = cdf_weibull)$L2_pmf
  
  # Weibull Distribution - n2
  
  x = sort(rweibull(n[2],shape))
  n_i <- hist(x, breaks = weibull_grid, plot = F)$counts
  
  res_wei_pmf_n_21[i] = L1_L2_pmf(n_i,weibull_grid,cdf = cdf_weibull)$L1_pmf
  res_wei_pmf_n_22[i] = L1_L2_pmf(n_i,weibull_grid,cdf = cdf_weibull)$L2_pmf
  
  # Weibull Distribution - n3
  
  x = sort(rweibull(n[3],shape))
  n_i <- hist(x, breaks = weibull_grid, plot = F)$counts
  
  res_wei_pmf_n_31[i] = L1_L2_pmf(n_i,weibull_grid,cdf = cdf_weibull)$L1_pmf
  res_wei_pmf_n_32[i] = L1_L2_pmf(n_i,weibull_grid,cdf = cdf_weibull)$L2_pmf
  
  # Weibull Distribution - n4
  
  x = sort(rweibull(n[4],shape))
  n_i <- hist(x, breaks = weibull_grid, plot = F)$counts
  
  res_wei_pmf_n_41[i] = L1_L2_pmf(n_i,weibull_grid,cdf = cdf_weibull)$L1_pmf
  res_wei_pmf_n_42[i] = L1_L2_pmf(n_i,weibull_grid,cdf = cdf_weibull)$L2_pmf
  
  # Weibull Distribution - n5
  
  x = sort(rweibull(n[5],shape))
  n_i <- hist(x, breaks = weibull_grid, plot = F)$counts
  
  res_wei_pmf_n_51[i] = L1_L2_pmf(n_i,weibull_grid,cdf = cdf_weibull)$L1_pmf
  res_wei_pmf_n_52[i] = L1_L2_pmf(n_i,weibull_grid,cdf = cdf_weibull)$L2_pmf
  
}

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L1_n/L1-PMF-2-Wei-n-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_wei_pmf_n_11,res_wei_pmf_n_21,res_wei_pmf_n_31,res_wei_pmf_n_41,res_wei_pmf_n_51,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors4,ylim=c(0,0.3))
segments(x0 = 0.29, y0 =0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L2_n/L2-PMF-2-Wei-n-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_wei_pmf_n_12,res_wei_pmf_n_22,res_wei_pmf_n_32,res_wei_pmf_n_42,res_wei_pmf_n_52,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors4,ylim=c(0,0.15))
segments(x0 = 0.29, y0 =0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

# Pareto Distribution - Sample Size Simulation 

MM = qpareto(0.99999999999999,location,shape_p) + 3
pareto_grid <- seq(location,MM, by = pareto_deltas[1])

res_par_pmf_n_11 = rep(0,B)
res_par_pmf_n_21 = rep(0,B)
res_par_pmf_n_31 = rep(0,B)
res_par_pmf_n_41 = rep(0,B)
res_par_pmf_n_51 = rep(0,B)

res_par_pmf_n_12 = rep(0,B)
res_par_pmf_n_22 = rep(0,B)
res_par_pmf_n_32 = rep(0,B)
res_par_pmf_n_42 = rep(0,B)
res_par_pmf_n_52 = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  # Pareto Distribution - n1
  
  x = sort(rpareto(n[1],location,shape_p))
  n_i <- hist(x, breaks = pareto_grid, plot = F)$counts
  
  res_par_pmf_n_11[i] = L1_L2_pmf(n_i,pareto_grid,cdf = cdf_pareto)$L1_pmf
  res_par_pmf_n_12[i] = L1_L2_pmf(n_i,pareto_grid,cdf = cdf_pareto)$L2_pmf
  
  # Pareto Distribution - n2
  
  x = sort(rpareto(n[2],location,shape_p))
  n_i <- hist(x, breaks = pareto_grid, plot = F)$counts
  
  res_par_pmf_n_21[i] = L1_L2_pmf(n_i,pareto_grid,cdf = cdf_pareto)$L1_pmf
  res_par_pmf_n_22[i] = L1_L2_pmf(n_i,pareto_grid,cdf = cdf_pareto)$L2_pmf
  
  # Pareto Distribution - n3
  
  x = sort(rpareto(n[3],location,shape_p))
  n_i <- hist(x, breaks = pareto_grid, plot = F)$counts
  
  res_par_pmf_n_31[i] = L1_L2_pmf(n_i,pareto_grid,cdf = cdf_pareto)$L1_pmf
  res_par_pmf_n_32[i] = L1_L2_pmf(n_i,pareto_grid,cdf = cdf_pareto)$L2_pmf
  
  # Pareto Distribution - n4
  
  x = sort(rpareto(n[4],location,shape_p))
  n_i <- hist(x, breaks = pareto_grid, plot = F)$counts
  
  res_par_pmf_n_41[i] = L1_L2_pmf(n_i,pareto_grid,cdf = cdf_pareto)$L1_pmf
  res_par_pmf_n_42[i] = L1_L2_pmf(n_i,pareto_grid,cdf = cdf_pareto)$L2_pmf
  
  # Pareto Distribution - n5
  
  x = sort(rpareto(n[5],location,shape_p))
  n_i <- hist(x, breaks = pareto_grid, plot = F)$counts
  
  res_par_pmf_n_51[i] = L1_L2_pmf(n_i,pareto_grid,cdf = cdf_pareto)$L1_pmf
  res_par_pmf_n_52[i] = L1_L2_pmf(n_i,pareto_grid,cdf = cdf_pareto)$L2_pmf
  
}


grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L1_n/L1-PMF-2-Pareto-n-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_par_pmf_n_11,res_par_pmf_n_21,res_par_pmf_n_31,res_par_pmf_n_41,res_par_pmf_n_51,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors5,ylim=c(0,0.3))
segments(x0 = 0.29, y0 = median(res_par_pmf_n_51), x1 = 10^5, y1 = median(res_par_pmf_n_51), col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L2_n/L2-PMF-2-Pareto-n-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_par_pmf_n_12,res_par_pmf_n_22,res_par_pmf_n_32,res_par_pmf_n_42,res_par_pmf_n_52,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors5,ylim=c(0,0.2))
segments(x0 = 0.29, y0 = median(res_par_pmf_n_52), x1 = 10^5, y1 = mean(res_par_pmf_n_52), col = "red", lty = 2)
dev.off()

n = c(10^3)
B = 500

# Laplace Distribution - Grid-Width Simulation 

MM = qlaplace(0.9999999,mu,sigma) + 3
laplace_grids = list(seq(-MM,MM, by = laplace_deltas[5]),
                     seq(-MM,MM, by = laplace_deltas[4]),
                     seq(-MM,MM, by = laplace_deltas[3]),
                     seq(-MM,MM, by = laplace_deltas[2]),
                     seq(-MM,MM, by = laplace_deltas[1]))

res_lap_pmf_gw_11 = rep(0,B)
res_lap_pmf_gw_21 = rep(0,B)
res_lap_pmf_gw_31 = rep(0,B)
res_lap_pmf_gw_41 = rep(0,B)
res_lap_pmf_gw_51 = rep(0,B)

res_lap_pmf_gw_12 = rep(0,B)
res_lap_pmf_gw_22 = rep(0,B)
res_lap_pmf_gw_32 = rep(0,B)
res_lap_pmf_gw_42 = rep(0,B)
res_lap_pmf_gw_52 = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  x <- sort(rlaplace(n, mean = mu,sd=sigma))
  
  # Laplace Distribution - GW1
  
  n_i <- hist(x, breaks = laplace_grids[[1]], plot = F)$counts
  
  res_lap_pmf_gw_11[i] = L1_L2_pmf(n_i,laplace_grids[[1]],cdf = cdf_laplace)$L1_pmf
  res_lap_pmf_gw_12[i] = L1_L2_pmf(n_i,laplace_grids[[1]],cdf = cdf_laplace)$L2_pmf
  
  # Laplace Distribution - GW2
  
  n_i <- hist(x, breaks = laplace_grids[[2]], plot = F)$counts
  
  res_lap_pmf_gw_21[i] = L1_L2_pmf(n_i,laplace_grids[[2]],cdf = cdf_laplace)$L1_pmf
  res_lap_pmf_gw_22[i] = L1_L2_pmf(n_i,laplace_grids[[2]],cdf = cdf_laplace)$L2_pmf
  
  # Laplace Distribution - GW3
  
  n_i <- hist(x, breaks = laplace_grids[[3]], plot = F)$counts
  
  res_lap_pmf_gw_31[i] = L1_L2_pmf(n_i,laplace_grids[[3]],cdf = cdf_laplace)$L1_pmf
  res_lap_pmf_gw_32[i] = L1_L2_pmf(n_i,laplace_grids[[3]],cdf = cdf_laplace)$L2_pmf
  
  # Laplace Distribution - GW4
  
  n_i <- hist(x, breaks = laplace_grids[[4]], plot = F)$counts
  
  res_lap_pmf_gw_41[i] = L1_L2_pmf(n_i,laplace_grids[[4]],cdf = cdf_laplace)$L1_pmf
  res_lap_pmf_gw_42[i] = L1_L2_pmf(n_i,laplace_grids[[4]],cdf = cdf_laplace)$L2_pmf
  
  # Laplace Distribution - GW5
  
  n_i <- hist(x, breaks = laplace_grids[[5]], plot = F)$counts
  
  res_lap_pmf_gw_51[i] = L1_L2_pmf(n_i,laplace_grids[[5]],cdf = cdf_laplace)$L1_pmf
  res_lap_pmf_gw_52[i] = L1_L2_pmf(n_i,laplace_grids[[5]],cdf = cdf_laplace)$L2_pmf
  
}

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L1_gw/L1-PMF-2-Laplace-gw-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_lap_pmf_gw_11,res_lap_pmf_gw_21,res_lap_pmf_gw_31,res_lap_pmf_gw_41,res_lap_pmf_gw_51,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,0.3))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L2_gw/L2-PMF-2-Laplace-gw-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_lap_pmf_gw_12,res_lap_pmf_gw_22,res_lap_pmf_gw_32,res_lap_pmf_gw_42,res_lap_pmf_gw_52,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,0.25))
dev.off()

# Chi-Square Distribution - Grid Width Simulation 

MM = qchisq(0.9999999,df) + 3
chisqr_grids = list(seq(0,MM, by = chisqr_deltas[5]),
                    seq(0,MM, by = chisqr_deltas[4]),
                    seq(0,MM, by = chisqr_deltas[3]),
                    seq(0,MM, by = chisqr_deltas[2]),
                    seq(0,MM, by = chisqr_deltas[1]))

res_chisq_pmf_gw_11 = rep(0,B)
res_chisq_pmf_gw_21 = rep(0,B)
res_chisq_pmf_gw_31 = rep(0,B)
res_chisq_pmf_gw_41 = rep(0,B)
res_chisq_pmf_gw_51 = rep(0,B)

res_chisq_pmf_gw_12 = rep(0,B)
res_chisq_pmf_gw_22 = rep(0,B)
res_chisq_pmf_gw_32 = rep(0,B)
res_chisq_pmf_gw_42 = rep(0,B)
res_chisq_pmf_gw_52 = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  x <- sort(rchisq(n,df))
  
  # Chi-Square Distribution - GW1
  
  n_i <- hist(x, breaks = chisqr_grids[[1]], plot = F)$counts
  
  res_chisq_pmf_gw_11[i] = L1_L2_pmf(n_i,chisqr_grids[[1]],cdf = cdf_chi)$L1_pmf
  res_chisq_pmf_gw_12[i] = L1_L2_pmf(n_i,chisqr_grids[[1]],cdf = cdf_chi)$L2_pmf
  
  # Chi-Square Distribution - GW2
  
  n_i <- hist(x, breaks = chisqr_grids[[2]], plot = F)$counts
  
  res_chisq_pmf_gw_21[i] = L1_L2_pmf(n_i,chisqr_grids[[2]],cdf = cdf_chi)$L1_pmf
  res_chisq_pmf_gw_22[i] = L1_L2_pmf(n_i,chisqr_grids[[2]],cdf = cdf_chi)$L2_pmf
  
  # Chi-Square Distribution - GW3
  
  n_i <- hist(x, breaks = chisqr_grids[[3]], plot = F)$counts
  
  res_chisq_pmf_gw_31[i] = L1_L2_pmf(n_i,chisqr_grids[[3]],cdf = cdf_chi)$L1_pmf
  res_chisq_pmf_gw_32[i] = L1_L2_pmf(n_i,chisqr_grids[[3]],cdf = cdf_chi)$L2_pmf
  
  # Chi-Square Distribution - GW4
  
  n_i <- hist(x, breaks = chisqr_grids[[4]], plot = F)$counts
  
  res_chisq_pmf_gw_41[i] = L1_L2_pmf(n_i,chisqr_grids[[4]],cdf = cdf_chi)$L1_pmf
  res_chisq_pmf_gw_42[i] = L1_L2_pmf(n_i,chisqr_grids[[4]],cdf = cdf_chi)$L2_pmf
  
  # Chi-Square Distribution - GW5
  
  n_i <- hist(x, breaks = chisqr_grids[[5]], plot = F)$counts
  
  res_chisq_pmf_gw_51[i] = L1_L2_pmf(n_i,chisqr_grids[[5]],cdf = cdf_chi)$L1_pmf
  res_chisq_pmf_gw_52[i] = L1_L2_pmf(n_i,chisqr_grids[[5]],cdf = cdf_chi)$L2_pmf
  
}

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L1_gw/L1-PMF-2-Chisq-gw-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_chisq_pmf_gw_11,res_chisq_pmf_gw_21,res_chisq_pmf_gw_31,res_chisq_pmf_gw_41,res_chisq_pmf_gw_51,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.1))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L2_gw/L2-PMF-2-Chisq-gw-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_chisq_pmf_gw_12,res_chisq_pmf_gw_22,res_chisq_pmf_gw_32,res_chisq_pmf_gw_42,res_chisq_pmf_gw_52,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.06))
dev.off()

# Log-Normal Distribution - Grid Width Simulation 

MM = qlnorm(0.9999999,mu,sigma) + 3
lnorm_grids = list(seq(0,MM, by = lnorm_deltas[5]),
                   seq(0,MM, by = lnorm_deltas[4]),
                   seq(0,MM, by = lnorm_deltas[3]),
                   seq(0,MM, by = lnorm_deltas[2]),
                   seq(0,MM, by = lnorm_deltas[1]))

res_logn_pmf_gw_11 = rep(0,B)
res_logn_pmf_gw_21 = rep(0,B)
res_logn_pmf_gw_31 = rep(0,B)
res_logn_pmf_gw_41 = rep(0,B)
res_logn_pmf_gw_51 = rep(0,B)

res_logn_pmf_gw_12 = rep(0,B)
res_logn_pmf_gw_22 = rep(0,B)
res_logn_pmf_gw_32 = rep(0,B)
res_logn_pmf_gw_42 = rep(0,B)
res_logn_pmf_gw_52 = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  x <- sort(rlnorm(n,mul,sigma))
  
  # Log-Normal Distribution - GW1

  n_i <- hist(x, breaks = lnorm_grids[[1]], plot = F)$counts
  
  res_logn_pmf_gw_11[i] = L1_L2_pmf(n_i,lnorm_grids[[1]],cdf = cdf_lnorm)$L1_pmf
  res_logn_pmf_gw_12[i] = L1_L2_pmf(n_i,lnorm_grids[[1]],cdf = cdf_lnorm)$L2_pmf
  
  # Log-Normal Distribution - GW2
  
  n_i <- hist(x, breaks = lnorm_grids[[2]], plot = F)$counts
  
  res_logn_pmf_gw_21[i] = L1_L2_pmf(n_i,lnorm_grids[[2]],cdf = cdf_lnorm)$L1_pmf
  res_logn_pmf_gw_22[i] = L1_L2_pmf(n_i,lnorm_grids[[2]],cdf = cdf_lnorm)$L2_pmf
  
  # Log-Normal Distribution - GW3
  
  n_i <- hist(x, breaks = lnorm_grids[[3]], plot = F)$counts
  
  res_logn_pmf_gw_31[i] = L1_L2_pmf(n_i,lnorm_grids[[3]],cdf = cdf_lnorm)$L1_pmf
  res_logn_pmf_gw_32[i] = L1_L2_pmf(n_i,lnorm_grids[[3]],cdf = cdf_lnorm)$L2_pmf
  
  # Log-Normal Distribution - GW4
  
  n_i <- hist(x, breaks = lnorm_grids[[4]], plot = F)$counts
  
  res_logn_pmf_gw_41[i] = L1_L2_pmf(n_i,lnorm_grids[[4]],cdf = cdf_lnorm)$L1_pmf
  res_logn_pmf_gw_42[i] = L1_L2_pmf(n_i,lnorm_grids[[4]],cdf = cdf_lnorm)$L2_pmf
  
  # Log-Normal Distribution - GW5
  
  n_i <- hist(x, breaks = lnorm_grids[[5]], plot = F)$counts
  
  res_logn_pmf_gw_51[i] = L1_L2_pmf(n_i,lnorm_grids[[5]],cdf = cdf_lnorm)$L1_pmf
  res_logn_pmf_gw_52[i] = L1_L2_pmf(n_i,lnorm_grids[[5]],cdf = cdf_lnorm)$L2_pmf
  
}

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L1_gw/L1-PMF-2-Lnorm-gw-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logn_pmf_gw_11,res_logn_pmf_gw_21,res_logn_pmf_gw_31,res_logn_pmf_gw_41,res_logn_pmf_gw_51,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,0.2))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L2_gw/L2-PMF-2-Lnorm-gw-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logn_pmf_gw_12,res_logn_pmf_gw_22,res_logn_pmf_gw_32,res_logn_pmf_gw_42,res_logn_pmf_gw_52,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,0.1))
dev.off()

# Weibull Distribution - Grid Width Simulation 

MM = qweibull(0.9999999,shape) + 3
weibull_grids = list(seq(-MM,MM, by = weibull_deltas[5]),
                     seq(-MM,MM, by = weibull_deltas[4]),
                     seq(-MM,MM, by = weibull_deltas[3]),
                     seq(-MM,MM, by = weibull_deltas[2]),
                     seq(-MM,MM, by = weibull_deltas[1]))

res_wei_pmf_gw_11 = rep(0,B)
res_wei_pmf_gw_21 = rep(0,B)
res_wei_pmf_gw_31 = rep(0,B)
res_wei_pmf_gw_41 = rep(0,B)
res_wei_pmf_gw_51 = rep(0,B)

res_wei_pmf_gw_12 = rep(0,B)
res_wei_pmf_gw_22 = rep(0,B)
res_wei_pmf_gw_32 = rep(0,B)
res_wei_pmf_gw_42 = rep(0,B)
res_wei_pmf_gw_52 = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  x = sort(rweibull(n,shape))
  
  # Weibull Distribution - GW1

  n_i <- hist(x, breaks = weibull_grids[[1]], plot = F)$counts
  
  res_wei_pmf_gw_11[i] = L1_L2_pmf(n_i,weibull_grids[[1]],cdf = cdf_weibull)$L1_pmf
  res_wei_pmf_gw_12[i] = L1_L2_pmf(n_i,weibull_grids[[1]],cdf = cdf_weibull)$L2_pmf
  
  # Weibull Distribution - GW2
  
  n_i <- hist(x, breaks = weibull_grids[[2]], plot = F)$counts
  
  res_wei_pmf_gw_21[i] = L1_L2_pmf(n_i,weibull_grids[[2]],cdf = cdf_weibull)$L1_pmf
  res_wei_pmf_gw_22[i] = L1_L2_pmf(n_i,weibull_grids[[2]],cdf = cdf_weibull)$L2_pmf
  
  # Weibull Distribution - GW3
  
  n_i <- hist(x, breaks = weibull_grids[[3]], plot = F)$counts
  
  res_wei_pmf_gw_31[i] = L1_L2_pmf(n_i,weibull_grids[[3]],cdf = cdf_weibull)$L1_pmf
  res_wei_pmf_gw_32[i] = L1_L2_pmf(n_i,weibull_grids[[3]],cdf = cdf_weibull)$L2_pmf
  
  # Weibull Distribution - GW4
  
  n_i <- hist(x, breaks = weibull_grids[[4]], plot = F)$counts
  
  res_wei_pmf_gw_41[i] = L1_L2_pmf(n_i,weibull_grids[[4]],cdf = cdf_weibull)$L1_pmf
  res_wei_pmf_gw_42[i] = L1_L2_pmf(n_i,weibull_grids[[4]],cdf = cdf_weibull)$L2_pmf
  
  # Weibull Distribution - GW5
  
  n_i <- hist(x, breaks = weibull_grids[[5]], plot = F)$counts
  
  res_wei_pmf_gw_51[i] = L1_L2_pmf(n_i,weibull_grids[[5]],cdf = cdf_weibull)$L1_pmf
  res_wei_pmf_gw_52[i] = L1_L2_pmf(n_i,weibull_grids[[5]],cdf = cdf_weibull)$L2_pmf
  
}

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L1_gw/L1-PMF-2-Wei-gw-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_wei_pmf_gw_11,res_wei_pmf_gw_21,res_wei_pmf_gw_31,res_wei_pmf_gw_41,res_wei_pmf_gw_51,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,0.1))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L2_gw/L2-PMF-2-Wei-gw-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_wei_pmf_gw_12,res_wei_pmf_gw_22,res_wei_pmf_gw_32,res_wei_pmf_gw_42,res_wei_pmf_gw_52,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,0.06))
dev.off()

# Pareto Distribution - Grid Width Simulation 

MM = qpareto(0.9999999999,location,shape_p) + 3
pareto_grids = list(seq(location,MM, by = pareto_deltas[1]),
                    seq(location,MM, by = pareto_deltas[2]),
                    seq(location,MM, by = pareto_deltas[3]),
                    seq(location,MM, by = pareto_deltas[4]),
                    seq(location,MM, by = pareto_deltas[5]))

res_par_pmf_gw_11 = rep(0,B)
res_par_pmf_gw_21 = rep(0,B)
res_par_pmf_gw_31 = rep(0,B)
res_par_pmf_gw_41 = rep(0,B)
res_par_pmf_gw_51 = rep(0,B)

res_par_pmf_gw_12 = rep(0,B)
res_par_pmf_gw_22 = rep(0,B)
res_par_pmf_gw_32 = rep(0,B)
res_par_pmf_gw_42 = rep(0,B)
res_par_pmf_gw_52 = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  x = sort(rpareto(n,location,shape_p))
  
  # Pareto Distribution - GW1
  
  n_i <- hist(x, breaks = pareto_grids[[1]], plot = F)$counts
  
  res_par_pmf_gw_11[i] = L1_L2_pmf(n_i,pareto_grids[[1]],cdf = cdf_pareto)$L1_pmf
  res_par_pmf_gw_12[i] = L1_L2_pmf(n_i,pareto_grids[[1]],cdf = cdf_pareto)$L2_pmf
  
  # Pareto Distribution - GW2
  
  n_i <- hist(x, breaks = pareto_grids[[2]], plot = F)$counts
  
  res_par_pmf_gw_21[i] = L1_L2_pmf(n_i,pareto_grids[[2]],cdf = cdf_pareto)$L1_pmf
  res_par_pmf_gw_22[i] = L1_L2_pmf(n_i,pareto_grids[[2]],cdf = cdf_pareto)$L2_pmf
  
  # Pareto Distribution - GW3
  
  n_i <- hist(x, breaks = pareto_grids[[3]], plot = F)$counts
  
  res_par_pmf_gw_31[i] = L1_L2_pmf(n_i,pareto_grids[[3]],cdf = cdf_pareto)$L1_pmf
  res_par_pmf_gw_32[i] = L1_L2_pmf(n_i,pareto_grids[[3]],cdf = cdf_pareto)$L2_pmf
  
  # Pareto Distribution - GW4
  
  n_i <- hist(x, breaks = pareto_grids[[4]], plot = F)$counts
  
  res_par_pmf_gw_41[i] = L1_L2_pmf(n_i,pareto_grids[[4]],cdf = cdf_pareto)$L1_pmf
  res_par_pmf_gw_42[i] = L1_L2_pmf(n_i,pareto_grids[[4]],cdf = cdf_pareto)$L2_pmf

  # Pareto Distribution - GW5
  
  n_i <- hist(x, breaks = pareto_grids[[5]], plot = F)$counts
  
  res_par_pmf_gw_51[i] = L1_L2_pmf(n_i,pareto_grids[[5]],cdf = cdf_pareto)$L1_pmf
  res_par_pmf_gw_52[i] = L1_L2_pmf(n_i,pareto_grids[[5]],cdf = cdf_pareto)$L2_pmf
  
}

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L1_gw/L1-PMF-2-Pareto-gw-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_par_pmf_gw_11,res_par_pmf_gw_21,res_par_pmf_gw_31,res_par_pmf_gw_41,res_par_pmf_gw_51,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,1))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PMF/L2_gw/L2-PMF-2-Pareto-gw-update.pdf",width = 16, height = 8)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_par_pmf_gw_12,res_par_pmf_gw_22,res_par_pmf_gw_32,res_par_pmf_gw_42,res_par_pmf_gw_52,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,1))
dev.off()

