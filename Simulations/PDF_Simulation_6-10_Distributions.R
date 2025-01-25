library(jmuOutlier)
library(EnvStats)

# Simulation Study 

# Range 

range_laplace = c(-Inf,Inf)
range_chisq = c(0,Inf)
range_lnorm = c(0,Inf)
range_wei = c(0,Inf)
range_pareto = c(0,Inf)

# Pdfs 

pdf_laplace <- function(x) {
  dlaplace(x,mu,sd = sigma)
}

pdf_chisq <- function(x) {
  dchisq(x,df)
}

pdf_lnorm <- function(x) {
  dlnorm(x,meanlog = mul,sdlog = sigma)
}

pdf_weibull <- function(x) {
  dweibull(x,shape)
}

pdf_pareto <- function(x) {
  dpareto(x,location = location,shape = shape_p)
}

n <- 1000
B = 100

res_pdf_laplace_L1 = rep(0,B)
res_pdf_chisq_L1 = rep(0,B)
res_pdf_lnorm_L1 = rep(0,B)
res_pdf_wei_L1 = rep(0,B)
res_pdf_pareto_L1 = rep(0,B)

res_pdf_laplace_L2 = rep(0,B)
res_pdf_chisq_L2 = rep(0,B)
res_pdf_lnorm_L2 = rep(0,B)
res_pdf_wei_L2 = rep(0,B)
res_pdf_pareto_L2 = rep(0,B)

# Parameters

alpha=2;beta=5;
scale_val = 0.05;df=3;
shape=2;location=2;
mu = 0;sigma=1;
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

# Convergence Limits 

pdf_L1_laplace_limit = mean(res_pdf_laplace_L1)
pdf_L2_laplace_limit = mean(res_pdf_laplace_L2)
pdf_L1_lnorm_limit = mean(res_pdf_lnorm_L1)
pdf_L2_lnorm_limit = mean(res_pdf_lnorm_L2)
pdf_L1_chisq_limit = mean(res_pdf_chisq_L1) 
pdf_L2_chisq_limit = mean(res_pdf_chisq_L2)
pdf_L1_wei_limit = mean(res_pdf_wei_L1)
pdf_L2_wei_limit = mean(res_pdf_wei_L2)
pdf_L1_pareto_limit = mean(res_pdf_pareto_L1)
pdf_L2_pareto_limit = mean(res_pdf_pareto_L2)

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

# Laplace Distribution - Sample Size Simulation 

n = c(10^1,10^2,10^3,10^4,10^5)
B = 100

MM = qlaplace(0.9999999,mu,sigma) + 3
laplace_grid <- seq(-MM,MM, by = laplace_deltas[1])

# Method-1 # 

res_lap_pdf_n_11_method1 = rep(0,B)
res_lap_pdf_n_21_method1 = rep(0,B)
res_lap_pdf_n_31_method1 = rep(0,B)
res_lap_pdf_n_41_method1 = rep(0,B)
res_lap_pdf_n_51_method1 = rep(0,B)

res_lap_pdf_n_12_method1 = rep(0,B)
res_lap_pdf_n_22_method1 = rep(0,B)
res_lap_pdf_n_32_method1 = rep(0,B)
res_lap_pdf_n_42_method1 = rep(0,B)
res_lap_pdf_n_52_method1 = rep(0,B)

# Method-2 # 

res_lap_pdf_n_11_method2 = rep(0,B)
res_lap_pdf_n_21_method2 = rep(0,B)
res_lap_pdf_n_31_method2 = rep(0,B)
res_lap_pdf_n_41_method2 = rep(0,B)
res_lap_pdf_n_51_method2 = rep(0,B)

res_lap_pdf_n_12_method2 = rep(0,B)
res_lap_pdf_n_22_method2 = rep(0,B)
res_lap_pdf_n_32_method2 = rep(0,B)
res_lap_pdf_n_42_method2 = rep(0,B)
res_lap_pdf_n_52_method2 = rep(0,B)

# Method-3 # 

res_lap_pdf_n_11_method3 = rep(0,B)
res_lap_pdf_n_21_method3 = rep(0,B)
res_lap_pdf_n_31_method3 = rep(0,B)
res_lap_pdf_n_41_method3 = rep(0,B)
res_lap_pdf_n_51_method3 = rep(0,B)

res_lap_pdf_n_12_method3 = rep(0,B)
res_lap_pdf_n_22_method3 = rep(0,B)
res_lap_pdf_n_32_method3 = rep(0,B)
res_lap_pdf_n_42_method3 = rep(0,B)
res_lap_pdf_n_52_method3 = rep(0,B)

# EM # 

res_lap_pdf_n_11_EM = rep(0,B)
res_lap_pdf_n_21_EM = rep(0,B)
res_lap_pdf_n_31_EM = rep(0,B)
res_lap_pdf_n_41_EM = rep(0,B)
res_lap_pdf_n_51_EM = rep(0,B)

res_lap_pdf_n_12_EM = rep(0,B)
res_lap_pdf_n_22_EM = rep(0,B)
res_lap_pdf_n_32_EM = rep(0,B)
res_lap_pdf_n_42_EM = rep(0,B)
res_lap_pdf_n_52_EM = rep(0,B)

set.seed(9)

for (i in 87:B) {
  
  # Laplace Distribution - n1
  
  x <- sort(rlaplace(n[1], mean = mu,sd=sigma))
  n_i <- hist(x, breaks = laplace_grid, plot = F)$counts
  
  res_lap_pdf_n_11_method1[i] = L1_Distance_method1(n_i,laplace_grid,pdf_laplace,range_laplace)
  res_lap_pdf_n_12_method1[i] = L2_Distance_method1(n_i,laplace_grid,pdf_laplace,range_laplace)
  
  res_lap_pdf_n_11_method2[i] = L1_Distance_method2(n_i,laplace_grid,pdf_laplace,range_laplace)
  res_lap_pdf_n_12_method2[i] = L2_Distance_method2(n_i,laplace_grid,pdf_laplace,range_laplace)
  
  res_lap_pdf_n_11_method3[i] = L1_Distance_method3(n_i,laplace_grid,pdf_laplace,range_laplace)
  res_lap_pdf_n_12_method3[i] = L2_Distance_method3(n_i,laplace_grid,pdf_laplace,range_laplace)
  
  res_lap_pdf_n_11_EM[i] = EM_L1_L2(n_i,laplace_grid,pdf_laplace,range_laplace)$L1
  res_lap_pdf_n_12_EM[i] = EM_L1_L2(n_i,laplace_grid,pdf_laplace,range_laplace)$L2
  
  # Laplace Distribution - n2
  
  x <- sort(rlaplace(n[2], mean = mu,sd=sigma))
  n_i <- hist(x, breaks = laplace_grid, plot = F)$counts
  
  res_lap_pdf_n_21_method1[i] = L1_Distance_method1(n_i,laplace_grid,pdf_laplace,range_laplace)
  res_lap_pdf_n_22_method1[i] = L2_Distance_method1(n_i,laplace_grid,pdf_laplace,range_laplace)
  
  res_lap_pdf_n_21_method2[i] = L1_Distance_method2(n_i,laplace_grid,pdf_laplace,range_laplace)
  res_lap_pdf_n_22_method2[i] = L2_Distance_method2(n_i,laplace_grid,pdf_laplace,range_laplace)
  
  res_lap_pdf_n_21_method3[i] = L1_Distance_method3(n_i,laplace_grid,pdf_laplace,range_laplace)
  res_lap_pdf_n_22_method3[i] = L2_Distance_method3(n_i,laplace_grid,pdf_laplace,range_laplace)
  
  res_lap_pdf_n_21_EM[i] = EM_L1_L2(n_i,laplace_grid,pdf_laplace,range_laplace)$L1
  res_lap_pdf_n_22_EM[i] = EM_L1_L2(n_i,laplace_grid,pdf_laplace,range_laplace)$L2
  
  # Laplace Distribution - n3
  
  x <- sort(rlaplace(n[3], mean = mu,sd=sigma))
  n_i <- hist(x, breaks = laplace_grid, plot = F)$counts
  
  res_lap_pdf_n_31_method1[i] = L1_Distance_method1(n_i,laplace_grid,pdf_laplace,range_laplace)
  res_lap_pdf_n_32_method1[i] = L2_Distance_method1(n_i,laplace_grid,pdf_laplace,range_laplace)
  
  res_lap_pdf_n_31_method2[i] = L1_Distance_method2(n_i,laplace_grid,pdf_laplace,range_laplace)
  res_lap_pdf_n_32_method2[i] = L2_Distance_method2(n_i,laplace_grid,pdf_laplace,range_laplace)

  res_lap_pdf_n_31_method3[i] = L1_Distance_method3(n_i,laplace_grid,pdf_laplace,range_laplace)
  res_lap_pdf_n_32_method3[i] = L2_Distance_method3(n_i,laplace_grid,pdf_laplace,range_laplace)
  
  res_lap_pdf_n_31_EM[i] = EM_L1_L2(n_i,laplace_grid,pdf_laplace,range_laplace)$L1
  res_lap_pdf_n_32_EM[i] = EM_L1_L2(n_i,laplace_grid,pdf_laplace,range_laplace)$L2
  
  # Laplace Distribution - n4
  
  x <- sort(rlaplace(n[4], mean = mu,sd=sigma))
  n_i <- hist(x, breaks = laplace_grid, plot = F)$counts
  
  res_lap_pdf_n_41_method1[i] = L1_Distance_method1(n_i,laplace_grid,pdf_laplace,range_laplace)
  res_lap_pdf_n_42_method1[i] = L2_Distance_method1(n_i,laplace_grid,pdf_laplace,range_laplace)
  
  res_lap_pdf_n_41_method2[i] = L1_Distance_method2(n_i,laplace_grid,pdf_laplace,range_laplace)
  res_lap_pdf_n_42_method2[i] = L2_Distance_method2(n_i,laplace_grid,pdf_laplace,range_laplace)
  
  res_lap_pdf_n_41_method3[i] = L1_Distance_method3(n_i,laplace_grid,pdf_laplace,range_laplace)
  res_lap_pdf_n_42_method3[i] = L2_Distance_method3(n_i,laplace_grid,pdf_laplace,range_laplace)
  
  res_lap_pdf_n_41_EM[i] = EM_L1_L2(n_i,laplace_grid,pdf_laplace,range_laplace)$L1
  res_lap_pdf_n_42_EM[i] = EM_L1_L2(n_i,laplace_grid,pdf_laplace,range_laplace)$L2
  
  # Laplace Distribution - n5
  
  x <- sort(rlaplace(n[5], mean = mu,sd=sigma))
  n_i <- hist(x, breaks = laplace_grid, plot = F)$counts
  
  res_lap_pdf_n_51_method1[i] = L1_Distance_method1(n_i,laplace_grid,pdf_laplace,range_laplace)
  res_lap_pdf_n_52_method1[i] = L2_Distance_method1(n_i,laplace_grid,pdf_laplace,range_laplace)
  
  res_lap_pdf_n_51_method2[i] = L1_Distance_method2(n_i,laplace_grid,pdf_laplace,range_laplace)
  res_lap_pdf_n_52_method2[i] = L2_Distance_method2(n_i,laplace_grid,pdf_laplace,range_laplace)
  
  res_lap_pdf_n_51_method3[i] = L1_Distance_method3(n_i,laplace_grid,pdf_laplace,range_laplace)
  res_lap_pdf_n_52_method3[i] = L2_Distance_method3(n_i,laplace_grid,pdf_laplace,range_laplace)
  
  res_lap_pdf_n_51_EM[i] = EM_L1_L2(n_i,laplace_grid,pdf_laplace,range_laplace)$L1
  res_lap_pdf_n_52_EM[i] = EM_L1_L2(n_i,laplace_grid,pdf_laplace,range_laplace)$L2
  
}


# Method-1 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L1-n/L1-PDF-2-Laplace-n-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_lap_pdf_n_11_method1,res_lap_pdf_n_21_method1,res_lap_pdf_n_31_method1,res_lap_pdf_n_41_method1,res_lap_pdf_n_51_method1,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors1,ylim=c(0,1.2))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L2-n/L2-PDF-2-Laplace-n-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_lap_pdf_n_12_method1,res_lap_pdf_n_22_method1,res_lap_pdf_n_32_method1,res_lap_pdf_n_42_method1,res_lap_pdf_n_52_method1,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors1,ylim=c(0,0.7))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

# Method-2 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L1-n/L1-PDF-2-Laplace-n-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_lap_pdf_n_11_method2,res_lap_pdf_n_21_method2,res_lap_pdf_n_31_method2,res_lap_pdf_n_41_method2,res_lap_pdf_n_51_method2,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors1,ylim=c(0,1.2))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L2-n/L2-PDF-2-Laplace-n-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_lap_pdf_n_12_method2,res_lap_pdf_n_22_method2,res_lap_pdf_n_32_method2,res_lap_pdf_n_42_method2,res_lap_pdf_n_52_method2,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors1,ylim=c(0,0.7))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

# Method-3 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L1-n/L1-PDF-2-Laplace-n-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_lap_pdf_n_11_method3,res_lap_pdf_n_21_method3,res_lap_pdf_n_31_method3,res_lap_pdf_n_41_method3,res_lap_pdf_n_51_method3,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors1,ylim=c(0,1.2))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L2-n/L2-PDF-2-Laplace-n-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_lap_pdf_n_12_method3,res_lap_pdf_n_22_method3,res_lap_pdf_n_32_method3,res_lap_pdf_n_42_method3,res_lap_pdf_n_52_method3,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors1,ylim=c(0,0.7))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

# EM # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L1-n/L1-PDF-2-Laplace-n-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_lap_pdf_n_11_EM,res_lap_pdf_n_21_EM,res_lap_pdf_n_31_EM,res_lap_pdf_n_41_EM,res_lap_pdf_n_51_EM,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors1,ylim=c(0,1.2))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L2-n/L2-PDF-2-Laplace-n-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_lap_pdf_n_12_EM,res_lap_pdf_n_22_EM,res_lap_pdf_n_32_EM,res_lap_pdf_n_42_EM,res_lap_pdf_n_52_EM,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors1,ylim=c(0,0.7))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

# Chi-Square Distribution - Sample Size Simulation 

MM = qchisq(0.999999999999,df) + 30
chisq_grid <- seq(0,MM, by = chisqr_deltas[1])

# Method-1 # 

res_chisq_pdf_n_11_method1 = rep(0,B)
res_chisq_pdf_n_21_method1 = rep(0,B)
res_chisq_pdf_n_31_method1 = rep(0,B)
res_chisq_pdf_n_41_method1 = rep(0,B)
res_chisq_pdf_n_51_method1 = rep(0,B)

res_chisq_pdf_n_12_method1 = rep(0,B)
res_chisq_pdf_n_22_method1 = rep(0,B)
res_chisq_pdf_n_32_method1 = rep(0,B)
res_chisq_pdf_n_42_method1 = rep(0,B)
res_chisq_pdf_n_52_method1 = rep(0,B)

# Method-2 # 

res_chisq_pdf_n_11_method2 = rep(0,B)
res_chisq_pdf_n_21_method2 = rep(0,B)
res_chisq_pdf_n_31_method2 = rep(0,B)
res_chisq_pdf_n_41_method2 = rep(0,B)
res_chisq_pdf_n_51_method2 = rep(0,B)

res_chisq_pdf_n_12_method2 = rep(0,B)
res_chisq_pdf_n_22_method2 = rep(0,B)
res_chisq_pdf_n_32_method2 = rep(0,B)
res_chisq_pdf_n_42_method2 = rep(0,B)
res_chisq_pdf_n_52_method2 = rep(0,B)

# Method-3 # 

res_chisq_pdf_n_11_method3 = rep(0,B)
res_chisq_pdf_n_21_method3 = rep(0,B)
res_chisq_pdf_n_31_method3 = rep(0,B)
res_chisq_pdf_n_41_method3 = rep(0,B)
res_chisq_pdf_n_51_method3 = rep(0,B)

res_chisq_pdf_n_12_method3 = rep(0,B)
res_chisq_pdf_n_22_method3 = rep(0,B)
res_chisq_pdf_n_32_method3 = rep(0,B)
res_chisq_pdf_n_42_method3 = rep(0,B)
res_chisq_pdf_n_52_method3 = rep(0,B)

# EM # 

res_chisq_pdf_n_11_EM = rep(0,B)
res_chisq_pdf_n_21_EM = rep(0,B)
res_chisq_pdf_n_31_EM = rep(0,B)
res_chisq_pdf_n_41_EM = rep(0,B)
res_chisq_pdf_n_51_EM = rep(0,B)

res_chisq_pdf_n_12_EM = rep(0,B)
res_chisq_pdf_n_22_EM = rep(0,B)
res_chisq_pdf_n_32_EM = rep(0,B)
res_chisq_pdf_n_42_EM = rep(0,B)
res_chisq_pdf_n_52_EM = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  # Chi-square Distribution - n1
  
  x   <- sort(rchisq(n[1],df))
  n_i <- hist(x, breaks = chisq_grid, plot = F)$counts
  
  res_chisq_pdf_n_11_method1[i] = L1_Distance_method1(n_i,chisq_grid,pdf_chisq,range_chisq)
  res_chisq_pdf_n_12_method1[i] = L2_Distance_method1(n_i,chisq_grid,pdf_chisq,range_chisq)
  
  res_chisq_pdf_n_11_method2[i] = L1_Distance_method2(n_i,chisq_grid,pdf_chisq,range_chisq)
  res_chisq_pdf_n_12_method2[i] = L2_Distance_method2(n_i,chisq_grid,pdf_chisq,range_chisq)
  
  res_chisq_pdf_n_11_method3[i] = L1_Distance_method3(n_i,chisq_grid,pdf_chisq,range_chisq)
  res_chisq_pdf_n_12_method3[i] = L2_Distance_method3(n_i,chisq_grid,pdf_chisq,range_chisq)
  
  res_chisq_pdf_n_11_EM[i] = EM_L1_L2(n_i,chisq_grid,pdf_chisq,range_chisq)$L1
  res_chisq_pdf_n_12_EM[i] = EM_L1_L2(n_i,chisq_grid,pdf_chisq,range_chisq)$L2
  
  # Chi-square Distribution - n2
  
  x   <- sort(rchisq(n[2],df))
  n_i <- hist(x, breaks = chisq_grid, plot = F)$counts
  
  res_chisq_pdf_n_21_method1[i] = L1_Distance_method1(n_i,chisq_grid,pdf_chisq,range_chisq)
  res_chisq_pdf_n_22_method1[i] = L2_Distance_method1(n_i,chisq_grid,pdf_chisq,range_chisq)
  
  res_chisq_pdf_n_21_method2[i] = L1_Distance_method2(n_i,chisq_grid,pdf_chisq,range_chisq)
  res_chisq_pdf_n_22_method2[i] = L2_Distance_method2(n_i,chisq_grid,pdf_chisq,range_chisq)
  
  res_chisq_pdf_n_21_method3[i] = L1_Distance_method3(n_i,chisq_grid,pdf_chisq,range_chisq)
  res_chisq_pdf_n_22_method3[i] = L2_Distance_method3(n_i,chisq_grid,pdf_chisq,range_chisq)
  
  res_chisq_pdf_n_21_EM[i] = EM_L1_L2(n_i,chisq_grid,pdf_chisq,range_chisq)$L1
  res_chisq_pdf_n_22_EM[i] = EM_L1_L2(n_i,chisq_grid,pdf_chisq,range_chisq)$L2
  
  # Chi-square Distribution - n3
  
  x   <- sort(rchisq(n[3],df))
  n_i <- hist(x, breaks = chisq_grid, plot = F)$counts
  
  res_chisq_pdf_n_31_method1[i] = L1_Distance_method1(n_i,chisq_grid,pdf_chisq,range_chisq)
  res_chisq_pdf_n_32_method1[i] = L2_Distance_method1(n_i,chisq_grid,pdf_chisq,range_chisq)
  
  res_chisq_pdf_n_31_method2[i] = L1_Distance_method2(n_i,chisq_grid,pdf_chisq,range_chisq)
  res_chisq_pdf_n_32_method2[i] = L2_Distance_method2(n_i,chisq_grid,pdf_chisq,range_chisq)
  
  res_chisq_pdf_n_31_method3[i] = L1_Distance_method3(n_i,chisq_grid,pdf_chisq,range_chisq)
  res_chisq_pdf_n_32_method3[i] = L2_Distance_method3(n_i,chisq_grid,pdf_chisq,range_chisq)
  
  res_chisq_pdf_n_31_EM[i] = EM_L1_L2(n_i,chisq_grid,pdf_chisq,range_chisq)$L1
  res_chisq_pdf_n_32_EM[i] = EM_L1_L2(n_i,chisq_grid,pdf_chisq,range_chisq)$L2
  
  # Chi-square Distribution - n4
  
  x   <- sort(rchisq(n[4],df))
  n_i <- hist(x, breaks = chisq_grid, plot = F)$counts
  
  res_chisq_pdf_n_41_method1[i] = L1_Distance_method1(n_i,chisq_grid,pdf_chisq,range_chisq)
  res_chisq_pdf_n_42_method1[i] = L2_Distance_method1(n_i,chisq_grid,pdf_chisq,range_chisq)
  
  res_chisq_pdf_n_41_method2[i] = L1_Distance_method2(n_i,chisq_grid,pdf_chisq,range_chisq)
  res_chisq_pdf_n_42_method2[i] = L2_Distance_method2(n_i,chisq_grid,pdf_chisq,range_chisq)
  
  res_chisq_pdf_n_41_method3[i] = L1_Distance_method3(n_i,chisq_grid,pdf_chisq,range_chisq)
  res_chisq_pdf_n_42_method3[i] = L2_Distance_method3(n_i,chisq_grid,pdf_chisq,range_chisq)
  
  res_chisq_pdf_n_41_EM[i] = EM_L1_L2(n_i,chisq_grid,pdf_chisq,range_chisq)$L1
  res_chisq_pdf_n_42_EM[i] = EM_L1_L2(n_i,chisq_grid,pdf_chisq,range_chisq)$L2
  
  # Chi-Square Distribution - n5
  
  x   <- sort(rchisq(n[5],df))
  n_i <- hist(x, breaks = chisq_grid, plot = F)$counts
  
  res_chisq_pdf_n_51_method1[i] = L1_Distance_method1(n_i,chisq_grid,pdf_chisq,range_chisq)
  res_chisq_pdf_n_52_method1[i] = L2_Distance_method1(n_i,chisq_grid,pdf_chisq,range_chisq)
  
  res_chisq_pdf_n_51_method2[i] = L1_Distance_method2(n_i,chisq_grid,pdf_chisq,range_chisq)
  res_chisq_pdf_n_52_method2[i] = L2_Distance_method2(n_i,chisq_grid,pdf_chisq,range_chisq)
  
  res_chisq_pdf_n_51_method3[i] = L1_Distance_method3(n_i,chisq_grid,pdf_chisq,range_chisq)
  res_chisq_pdf_n_52_method3[i] = L2_Distance_method3(n_i,chisq_grid,pdf_chisq,range_chisq)
  
  res_chisq_pdf_n_51_EM[i] = EM_L1_L2(n_i,chisq_grid,pdf_chisq,range_chisq)$L1
  res_chisq_pdf_n_52_EM[i] = EM_L1_L2(n_i,chisq_grid,pdf_chisq,range_chisq)$L2
  
}

# Method-1 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L1-n/L1-PDF-2-Chisq-n-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_chisq_pdf_n_11_method1,res_chisq_pdf_n_21_method1,res_chisq_pdf_n_31_method1,res_chisq_pdf_n_41_method1,res_chisq_pdf_n_51_method1,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2,ylim=c(0,1.0))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L2-n/L2-PDF-2-Chisq-n-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_chisq_pdf_n_12_method1,res_chisq_pdf_n_22_method1,res_chisq_pdf_n_32_method1,res_chisq_pdf_n_42_method1,res_chisq_pdf_n_52_method1,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2,ylim=c(0,0.5))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

# Method-2 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L1-n/L1-PDF-2-Chisq-n-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_chisq_pdf_n_11_method2,res_chisq_pdf_n_21_method2,res_chisq_pdf_n_31_method2,res_chisq_pdf_n_41_method2,res_chisq_pdf_n_51_method2,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2,ylim=c(0,1.0))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L2-n/L2-PDF-2-Chisq-n-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_chisq_pdf_n_12_method2,res_chisq_pdf_n_22_method2,res_chisq_pdf_n_32_method2,res_chisq_pdf_n_42_method2,res_chisq_pdf_n_52_method2,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2,ylim=c(0,0.5))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

# Method-3 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L1-n/L1-PDF-2-Chisq-n-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_chisq_pdf_n_11_method3,res_chisq_pdf_n_21_method3,res_chisq_pdf_n_31_method3,res_chisq_pdf_n_41_method3,res_chisq_pdf_n_51_method3,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2,ylim=c(0,1.0))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L2-n/L2-PDF-2-Chisq-n-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_chisq_pdf_n_12_method3,res_chisq_pdf_n_22_method3,res_chisq_pdf_n_32_method3,res_chisq_pdf_n_42_method3,res_chisq_pdf_n_52_method3,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2,ylim=c(0,0.5))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

# EM # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L1-n/L1-PDF-2-Chisq-n-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_chisq_pdf_n_11_EM,res_chisq_pdf_n_21_EM,res_chisq_pdf_n_31_EM,res_chisq_pdf_n_41_EM,res_chisq_pdf_n_51_EM,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2,ylim=c(0,1.0))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L2-n/L2-PDF-2-Chisq-n-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_chisq_pdf_n_12_EM,res_chisq_pdf_n_22_EM,res_chisq_pdf_n_32_EM,res_chisq_pdf_n_42_EM,res_chisq_pdf_n_52_EM,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2,ylim=c(0,0.5))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

# Log-Normal Distribution - Sample Size Simulation 

MM = qlnorm(0.99999999999999,mu,sigma) + 30
lnrom_grid <- seq(0,MM, by = lnorm_deltas[1])

# Method-1 # 

res_logn_pdf_n_11_method1 = rep(0,B)
res_logn_pdf_n_21_method1 = rep(0,B)
res_logn_pdf_n_31_method1 = rep(0,B)
res_logn_pdf_n_41_method1 = rep(0,B)
res_logn_pdf_n_51_method1 = rep(0,B)

res_logn_pdf_n_12_method1 = rep(0,B)
res_logn_pdf_n_22_method1 = rep(0,B)
res_logn_pdf_n_32_method1 = rep(0,B)
res_logn_pdf_n_42_method1 = rep(0,B)
res_logn_pdf_n_52_method1 = rep(0,B)

# Method-2 # 

res_logn_pdf_n_11_method2 = rep(0,B)
res_logn_pdf_n_21_method2 = rep(0,B)
res_logn_pdf_n_31_method2 = rep(0,B)
res_logn_pdf_n_41_method2 = rep(0,B)
res_logn_pdf_n_51_method2 = rep(0,B)

res_logn_pdf_n_12_method2 = rep(0,B)
res_logn_pdf_n_22_method2 = rep(0,B)
res_logn_pdf_n_32_method2 = rep(0,B)
res_logn_pdf_n_42_method2 = rep(0,B)
res_logn_pdf_n_52_method2 = rep(0,B)

# Method-3 # 

res_logn_pdf_n_11_method3 = rep(0,B)
res_logn_pdf_n_21_method3 = rep(0,B)
res_logn_pdf_n_31_method3 = rep(0,B)
res_logn_pdf_n_41_method3 = rep(0,B)
res_logn_pdf_n_51_method3 = rep(0,B)

res_logn_pdf_n_12_method3 = rep(0,B)
res_logn_pdf_n_22_method3 = rep(0,B)
res_logn_pdf_n_32_method3 = rep(0,B)
res_logn_pdf_n_42_method3 = rep(0,B)
res_logn_pdf_n_52_method3 = rep(0,B)

# EM # 

res_logn_pdf_n_11_EM = rep(0,B)
res_logn_pdf_n_21_EM = rep(0,B)
res_logn_pdf_n_31_EM = rep(0,B)
res_logn_pdf_n_41_EM = rep(0,B)
res_logn_pdf_n_51_EM = rep(0,B)

res_logn_pdf_n_12_EM = rep(0,B)
res_logn_pdf_n_22_EM = rep(0,B)
res_logn_pdf_n_32_EM = rep(0,B)
res_logn_pdf_n_42_EM = rep(0,B)
res_logn_pdf_n_52_EM = rep(0,B)

set.seed(9)

for (i in 99:B) {
  
  # Log-Normal Distribution - n1
  
  x   <- sort(rlnorm(n[1],mul,sigma))
  n_i <- hist(x, breaks = lnrom_grid, plot = F)$counts
  
  res_logn_pdf_n_11_method1[i] = L1_Distance_method1(n_i,lnrom_grid,pdf_lnorm,range_lnorm)
  res_logn_pdf_n_12_method1[i] = L2_Distance_method1(n_i,lnrom_grid,pdf_lnorm,range_lnorm)
  
  res_logn_pdf_n_11_method2[i] = L1_Distance_method2(n_i,lnrom_grid,pdf_lnorm,range_lnorm)
  res_logn_pdf_n_12_method2[i] = L2_Distance_method2(n_i,lnrom_grid,pdf_lnorm,range_lnorm)
  
  res_logn_pdf_n_11_method3[i] = L1_Distance_method3(n_i,lnrom_grid,pdf_lnorm,range_lnorm)
  res_logn_pdf_n_12_method3[i] = L2_Distance_method3(n_i,lnrom_grid,pdf_lnorm,range_lnorm)
  
  res_logn_pdf_n_11_EM[i] = EM_L1_L2(n_i,lnrom_grid,pdf_lnorm,range_lnorm)$L1
  res_logn_pdf_n_12_EM[i] = EM_L1_L2(n_i,lnrom_grid,pdf_lnorm,range_lnorm)$L2
  
  # Log-Normal Distribution - n2
  
  x   <- sort(rlnorm(n[2],mul,sigma))
  n_i <- hist(x, breaks = lnrom_grid, plot = F)$counts
  
  res_logn_pdf_n_21_method1[i] = L1_Distance_method1(n_i,lnrom_grid,pdf_lnorm,range_lnorm)
  res_logn_pdf_n_22_method1[i] = L2_Distance_method1(n_i,lnrom_grid,pdf_lnorm,range_lnorm)
  
  res_logn_pdf_n_21_method2[i] = L1_Distance_method2(n_i,lnrom_grid,pdf_lnorm,range_lnorm)
  res_logn_pdf_n_22_method2[i] = L2_Distance_method2(n_i,lnrom_grid,pdf_lnorm,range_lnorm)
  
  res_logn_pdf_n_21_method3[i] = L1_Distance_method3(n_i,lnrom_grid,pdf_lnorm,range_lnorm)
  res_logn_pdf_n_22_method3[i] = L2_Distance_method3(n_i,lnrom_grid,pdf_lnorm,range_lnorm)
  
  res_logn_pdf_n_21_EM[i] = EM_L1_L2(n_i,lnrom_grid,pdf_lnorm,range_lnorm)$L1
  res_logn_pdf_n_22_EM[i] = EM_L1_L2(n_i,lnrom_grid,pdf_lnorm,range_lnorm)$L2
  
  # Log-Normal Distribution - n3
  
  x   <- sort(rlnorm(n[3],mul,sigma))
  n_i <- hist(x, breaks = lnrom_grid, plot = F)$counts
  
  res_logn_pdf_n_31_method1[i] = L1_Distance_method1(n_i,lnrom_grid,pdf_lnorm,range_lnorm)
  res_logn_pdf_n_32_method1[i] = L2_Distance_method1(n_i,lnrom_grid,pdf_lnorm,range_lnorm)
  
  res_logn_pdf_n_31_method2[i] = L1_Distance_method2(n_i,lnrom_grid,pdf_lnorm,range_lnorm)
  res_logn_pdf_n_32_method2[i] = L2_Distance_method2(n_i,lnrom_grid,pdf_lnorm,range_lnorm)
  
  res_logn_pdf_n_31_method3[i] = L1_Distance_method3(n_i,lnrom_grid,pdf_lnorm,range_lnorm)
  res_logn_pdf_n_32_method3[i] = L2_Distance_method3(n_i,lnrom_grid,pdf_lnorm,range_lnorm)
  
  res_logn_pdf_n_31_EM[i] = EM_L1_L2(n_i,lnrom_grid,pdf_lnorm,range_lnorm)$L1
  res_logn_pdf_n_32_EM[i] = EM_L1_L2(n_i,lnrom_grid,pdf_lnorm,range_lnorm)$L2
  
  # Log-Normal Distribution - n4
  
  x   <- sort(rlnorm(n[4],mul,sigma))
  n_i <- hist(x, breaks = lnrom_grid, plot = F)$counts
  
  res_logn_pdf_n_41_method1[i] = L1_Distance_method1(n_i,lnrom_grid,pdf_lnorm,range_lnorm)
  res_logn_pdf_n_42_method1[i] = L2_Distance_method1(n_i,lnrom_grid,pdf_lnorm,range_lnorm)
  
  res_logn_pdf_n_41_method2[i] = L1_Distance_method2(n_i,lnrom_grid,pdf_lnorm,range_lnorm)
  res_logn_pdf_n_42_method2[i] = L2_Distance_method2(n_i,lnrom_grid,pdf_lnorm,range_lnorm)
  
  res_logn_pdf_n_41_method3[i] = L1_Distance_method3(n_i,lnrom_grid,pdf_lnorm,range_lnorm)
  res_logn_pdf_n_42_method3[i] = L2_Distance_method3(n_i,lnrom_grid,pdf_lnorm,range_lnorm)
  
  res_logn_pdf_n_41_EM[i] = EM_L1_L2(n_i,lnrom_grid,pdf_lnorm,range_lnorm)$L1
  res_logn_pdf_n_42_EM[i] = EM_L1_L2(n_i,lnrom_grid,pdf_lnorm,range_lnorm)$L2
  
  # Log-Normal Distribution - n5
  
  x   <- sort(rlnorm(n[5],mul,sigma))
  n_i <- hist(x, breaks = lnrom_grid, plot = F)$counts
  
  res_logn_pdf_n_51_method1[i] = L1_Distance_method1(n_i,lnrom_grid,pdf_lnorm,range_lnorm)
  res_logn_pdf_n_52_method1[i] = L2_Distance_method1(n_i,lnrom_grid,pdf_lnorm,range_lnorm)
  
  res_logn_pdf_n_51_method2[i] = L1_Distance_method2(n_i,lnrom_grid,pdf_lnorm,range_lnorm)
  res_logn_pdf_n_52_method2[i] = L2_Distance_method2(n_i,lnrom_grid,pdf_lnorm,range_lnorm)
  
  res_logn_pdf_n_51_method3[i] = L1_Distance_method3(n_i,lnrom_grid,pdf_lnorm,range_lnorm)
  res_logn_pdf_n_52_method3[i] = L2_Distance_method3(n_i,lnrom_grid,pdf_lnorm,range_lnorm)
  
  res_logn_pdf_n_51_EM[i] = EM_L1_L2(n_i,lnrom_grid,pdf_lnorm,range_lnorm)$L1
  res_logn_pdf_n_52_EM[i] = EM_L1_L2(n_i,lnrom_grid,pdf_lnorm,range_lnorm)$L2
  
}

# Method-1 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L1-n/L1-PDF-2-Lnorm-n-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logn_pdf_n_11_method1,res_logn_pdf_n_21_method1,res_logn_pdf_n_31_method1,res_logn_pdf_n_41_method1,res_logn_pdf_n_51_method1,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3,ylim=c(0,1.6))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L2-n/L2-PDF-2-Lnorm-n-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logn_pdf_n_12_method1,res_logn_pdf_n_22_method1,res_logn_pdf_n_32_method1,res_logn_pdf_n_42_method1,res_logn_pdf_n_52_method1,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3,ylim=c(0,0.4))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

# Method-2 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L1-n/L1-PDF-2-Lnorm-n-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logn_pdf_n_11_method2,res_logn_pdf_n_21_method2,res_logn_pdf_n_31_method2,res_logn_pdf_n_41_method2,res_logn_pdf_n_51_method2,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3,ylim=c(0,1.6))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L2-n/L2-PDF-2-Lnorm-n-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logn_pdf_n_12_method2,res_logn_pdf_n_22_method2,res_logn_pdf_n_32_method2,res_logn_pdf_n_42_method2,res_logn_pdf_n_52_method2,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3,ylim=c(0,0.4))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

# Method-3 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L1-n/L1-PDF-2-Lnorm-n-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logn_pdf_n_11_method3,res_logn_pdf_n_21_method3,res_logn_pdf_n_31_method3,res_logn_pdf_n_41_method3,res_logn_pdf_n_51_method3,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3,ylim=c(0,1.6))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L2-n/L2-PDF-2-Lnorm-n-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logn_pdf_n_12_method3,res_logn_pdf_n_22_method3,res_logn_pdf_n_32_method3,res_logn_pdf_n_42_method3,res_logn_pdf_n_52_method3,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3,ylim=c(0,0.4))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

# EM # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L1-n/L1-PDF-2-Lnorm-n-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logn_pdf_n_11_EM,res_logn_pdf_n_21_EM,res_logn_pdf_n_31_EM,res_logn_pdf_n_41_EM,res_logn_pdf_n_51_EM,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3,ylim=c(0,1.6))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L2-n/L2-PDF-2-Lnorm-n-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logn_pdf_n_12_EM,res_logn_pdf_n_22_EM,res_logn_pdf_n_32_EM,res_logn_pdf_n_42_EM,res_logn_pdf_n_52_EM,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3,ylim=c(0,0.4))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

# Weibull Distribution - Sample Size Simulation 

MM = qweibull(0.9999999,shape) + 3
weibull_grid <- seq(-MM,MM, by = weibull_deltas[1])

# Method-1 # 

res_wei_pdf_n_11_method1 = rep(0,B)
res_wei_pdf_n_21_method1 = rep(0,B)
res_wei_pdf_n_31_method1 = rep(0,B)
res_wei_pdf_n_41_method1 = rep(0,B)
res_wei_pdf_n_51_method1 = rep(0,B)

res_wei_pdf_n_12_method1 = rep(0,B)
res_wei_pdf_n_22_method1 = rep(0,B)
res_wei_pdf_n_32_method1 = rep(0,B)
res_wei_pdf_n_42_method1 = rep(0,B)
res_wei_pdf_n_52_method1 = rep(0,B)

# Method-2 # 

res_wei_pdf_n_11_method2 = rep(0,B)
res_wei_pdf_n_21_method2 = rep(0,B)
res_wei_pdf_n_31_method2 = rep(0,B)
res_wei_pdf_n_41_method2 = rep(0,B)
res_wei_pdf_n_51_method2 = rep(0,B)

res_wei_pdf_n_12_method2 = rep(0,B)
res_wei_pdf_n_22_method2 = rep(0,B)
res_wei_pdf_n_32_method2 = rep(0,B)
res_wei_pdf_n_42_method2 = rep(0,B)
res_wei_pdf_n_52_method2 = rep(0,B)

# Method-3 # 

res_wei_pdf_n_11_method3 = rep(0,B)
res_wei_pdf_n_21_method3 = rep(0,B)
res_wei_pdf_n_31_method3 = rep(0,B)
res_wei_pdf_n_41_method3 = rep(0,B)
res_wei_pdf_n_51_method3 = rep(0,B)

res_wei_pdf_n_12_method3 = rep(0,B)
res_wei_pdf_n_22_method3 = rep(0,B)
res_wei_pdf_n_32_method3 = rep(0,B)
res_wei_pdf_n_42_method3 = rep(0,B)
res_wei_pdf_n_52_method3 = rep(0,B)

# EM # 

res_wei_pdf_n_11_EM = rep(0,B)
res_wei_pdf_n_21_EM = rep(0,B)
res_wei_pdf_n_31_EM = rep(0,B)
res_wei_pdf_n_41_EM = rep(0,B)
res_wei_pdf_n_51_EM = rep(0,B)

res_wei_pdf_n_12_EM = rep(0,B)
res_wei_pdf_n_22_EM = rep(0,B)
res_wei_pdf_n_32_EM = rep(0,B)
res_wei_pdf_n_42_EM = rep(0,B)
res_wei_pdf_n_52_EM = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  # Weibull Distribution - n1
  
  x = sort(rweibull(n[1],shape))
  n_i <- hist(x, breaks = weibull_grid, plot = F)$counts
  
  res_wei_pdf_n_11_method1[i] = L1_Distance_method1(n_i,weibull_grid,pdf_weibull,range_wei)
  res_wei_pdf_n_12_method1[i] = L2_Distance_method1(n_i,weibull_grid,pdf_weibull,range_wei)
  
  res_wei_pdf_n_11_method2[i] = L1_Distance_method2(n_i,weibull_grid,pdf_weibull,range_wei)
  res_wei_pdf_n_12_method2[i] = L2_Distance_method2(n_i,weibull_grid,pdf_weibull,range_wei)
  
  res_wei_pdf_n_11_method3[i] = L1_Distance_method3(n_i,weibull_grid,pdf_weibull,range_wei)
  res_wei_pdf_n_12_method3[i] = L2_Distance_method3(n_i,weibull_grid,pdf_weibull,range_wei)
  
  res_wei_pdf_n_11_EM[i] = EM_L1_L2(n_i,weibull_grid,pdf_weibull,range_wei)$L1
  res_wei_pdf_n_12_EM[i] = EM_L1_L2(n_i,weibull_grid,pdf_weibull,range_wei)$L2
  
  # Weibull Distribution - n2
  
  x = sort(rweibull(n[2],shape))
  n_i <- hist(x, breaks = weibull_grid, plot = F)$counts
  
  res_wei_pdf_n_21_method1[i] = L1_Distance_method1(n_i,weibull_grid,pdf_weibull,range_wei)
  res_wei_pdf_n_22_method1[i] = L2_Distance_method1(n_i,weibull_grid,pdf_weibull,range_wei)
  
  res_wei_pdf_n_21_method2[i] = L1_Distance_method2(n_i,weibull_grid,pdf_weibull,range_wei)
  res_wei_pdf_n_22_method2[i] = L2_Distance_method2(n_i,weibull_grid,pdf_weibull,range_wei)
  
  res_wei_pdf_n_21_method3[i] = L1_Distance_method3(n_i,weibull_grid,pdf_weibull,range_wei)
  res_wei_pdf_n_22_method3[i] = L2_Distance_method3(n_i,weibull_grid,pdf_weibull,range_wei)
  
  res_wei_pdf_n_21_EM[i] = EM_L1_L2(n_i,weibull_grid,pdf_weibull,range_wei)$L1
  res_wei_pdf_n_22_EM[i] = EM_L1_L2(n_i,weibull_grid,pdf_weibull,range_wei)$L2
  
  # Weibull Distribution - n3
  
  x = sort(rweibull(n[3],shape))
  n_i <- hist(x, breaks = weibull_grid, plot = F)$counts
  
  res_wei_pdf_n_31_method1[i] = L1_Distance_method1(n_i,weibull_grid,pdf_weibull,range_wei)
  res_wei_pdf_n_32_method1[i] = L2_Distance_method1(n_i,weibull_grid,pdf_weibull,range_wei)
  
  res_wei_pdf_n_31_method2[i] = L1_Distance_method2(n_i,weibull_grid,pdf_weibull,range_wei)
  res_wei_pdf_n_32_method2[i] = L2_Distance_method2(n_i,weibull_grid,pdf_weibull,range_wei)
  
  res_wei_pdf_n_31_method3[i] = L1_Distance_method3(n_i,weibull_grid,pdf_weibull,range_wei)
  res_wei_pdf_n_32_method3[i] = L2_Distance_method3(n_i,weibull_grid,pdf_weibull,range_wei)
  
  res_wei_pdf_n_31_EM[i] = EM_L1_L2(n_i,weibull_grid,pdf_weibull,range_wei)$L1
  res_wei_pdf_n_32_EM[i] = EM_L1_L2(n_i,weibull_grid,pdf_weibull,range_wei)$L2
  
  # Weibull Distribution - n4
  
  x = sort(rweibull(n[4],shape))
  n_i <- hist(x, breaks = weibull_grid, plot = F)$counts
  
  res_wei_pdf_n_41_method1[i] = L1_Distance_method1(n_i,weibull_grid,pdf_weibull,range_wei)
  res_wei_pdf_n_42_method1[i] = L2_Distance_method1(n_i,weibull_grid,pdf_weibull,range_wei)
  
  res_wei_pdf_n_41_method2[i] = L1_Distance_method2(n_i,weibull_grid,pdf_weibull,range_wei)
  res_wei_pdf_n_42_method2[i] = L2_Distance_method2(n_i,weibull_grid,pdf_weibull,range_wei)
  
  res_wei_pdf_n_41_method3[i] = L1_Distance_method3(n_i,weibull_grid,pdf_weibull,range_wei)
  res_wei_pdf_n_42_method3[i] = L2_Distance_method3(n_i,weibull_grid,pdf_weibull,range_wei)
  
  res_wei_pdf_n_41_EM[i] = EM_L1_L2(n_i,weibull_grid,pdf_weibull,range_wei)$L1
  res_wei_pdf_n_42_EM[i] = EM_L1_L2(n_i,weibull_grid,pdf_weibull,range_wei)$L2
  
  # Weibull Distribution - n5
  
  x = sort(rweibull(n[5],shape))
  n_i <- hist(x, breaks = weibull_grid, plot = F)$counts
  
  res_wei_pdf_n_51_method1[i] = L1_Distance_method1(n_i,weibull_grid,pdf_weibull,range_wei)
  res_wei_pdf_n_52_method1[i] = L2_Distance_method1(n_i,weibull_grid,pdf_weibull,range_wei)
  
  res_wei_pdf_n_51_method2[i] = L1_Distance_method2(n_i,weibull_grid,pdf_weibull,range_wei)
  res_wei_pdf_n_52_method2[i] = L2_Distance_method2(n_i,weibull_grid,pdf_weibull,range_wei)
  
  res_wei_pdf_n_51_method3[i] = L1_Distance_method3(n_i,weibull_grid,pdf_weibull,range_wei)
  res_wei_pdf_n_52_method3[i] = L2_Distance_method3(n_i,weibull_grid,pdf_weibull,range_wei)
  
  res_wei_pdf_n_51_EM[i] = EM_L1_L2(n_i,weibull_grid,pdf_weibull,range_wei)$L1
  res_wei_pdf_n_52_EM[i] = EM_L1_L2(n_i,weibull_grid,pdf_weibull,range_wei)$L2
  
}


# Method-1 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L1-n/L1-PDF-2-Wei-n-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_wei_pdf_n_11_method1,res_wei_pdf_n_21_method1,res_wei_pdf_n_31_method1,res_wei_pdf_n_41_method1,res_wei_pdf_n_51_method1,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors4,ylim=c(0,1.0))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L2-n/L2-PDF-2-Wei-n-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_wei_pdf_n_12_method1,res_wei_pdf_n_22_method1,res_wei_pdf_n_32_method1,res_wei_pdf_n_42_method1,res_wei_pdf_n_52_method1,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors4,ylim=c(0,0.8))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 =  0, col = "red", lty = 2)
dev.off()

# Method-2 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L1-n/L1-PDF-2-Wei-n-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_wei_pdf_n_11_method2,res_wei_pdf_n_21_method2,res_wei_pdf_n_31_method2,res_wei_pdf_n_41_method2,res_wei_pdf_n_51_method2,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors4,ylim=c(0,1.0))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L2-n/L2-PDF-2-Wei-n-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_wei_pdf_n_12_method2,res_wei_pdf_n_22_method2,res_wei_pdf_n_32_method2,res_wei_pdf_n_42_method2,res_wei_pdf_n_52_method2,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors4,ylim=c(0,0.8))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 =  0, col = "red", lty = 2)
dev.off()

# Method-3 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L1-n/L1-PDF-2-Wei-n-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_wei_pdf_n_11_method3,res_wei_pdf_n_21_method3,res_wei_pdf_n_31_method3,res_wei_pdf_n_41_method3,res_wei_pdf_n_51_method3,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors4,ylim=c(0,1))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L2-n/L2-PDF-2-Wei-n-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_wei_pdf_n_12_method3,res_wei_pdf_n_22_method3,res_wei_pdf_n_32_method3,res_wei_pdf_n_42_method3,res_wei_pdf_n_52_method3,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors4,ylim=c(0,0.8))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 =  0, col = "red", lty = 2)
dev.off()

# EM # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L1-n/L1-PDF-2-Wei-n-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_wei_pdf_n_11_EM,res_wei_pdf_n_21_EM,res_wei_pdf_n_31_EM,res_wei_pdf_n_41_EM,res_wei_pdf_n_51_EM,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors4,ylim=c(0,1))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L2-n/L2-PDF-2-Wei-n-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_wei_pdf_n_12_EM,res_wei_pdf_n_22_EM,res_wei_pdf_n_32_EM,res_wei_pdf_n_42_EM,res_wei_pdf_n_52_EM,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors4,ylim=c(0,0.8))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 =  0, col = "red", lty = 2)
dev.off()

# Pareto Distribution - Sample Size Simulation 

MM = qpareto(0.999999999,location,shape_p) + 3
pareto_grid <- seq(location,MM, by = pareto_deltas[1])

# Method-1 # 

res_pareto_pdf_n_11_method1 = rep(0,B)
res_pareto_pdf_n_21_method1 = rep(0,B)
res_pareto_pdf_n_31_method1 = rep(0,B)
res_pareto_pdf_n_41_method1 = rep(0,B)
res_pareto_pdf_n_51_method1 = rep(0,B)

res_pareto_pdf_n_12_method1 = rep(0,B)
res_pareto_pdf_n_22_method1 = rep(0,B)
res_pareto_pdf_n_32_method1 = rep(0,B)
res_pareto_pdf_n_42_method1 = rep(0,B)
res_pareto_pdf_n_52_method1 = rep(0,B)

# Method-2 # 

res_pareto_pdf_n_11_method2 = rep(0,B)
res_pareto_pdf_n_21_method2 = rep(0,B)
res_pareto_pdf_n_31_method2 = rep(0,B)
res_pareto_pdf_n_41_method2 = rep(0,B)
res_pareto_pdf_n_51_method2 = rep(0,B)

res_pareto_pdf_n_12_method2 = rep(0,B)
res_pareto_pdf_n_22_method2 = rep(0,B)
res_pareto_pdf_n_32_method2 = rep(0,B)
res_pareto_pdf_n_42_method2 = rep(0,B)
res_pareto_pdf_n_52_method2 = rep(0,B)

# Method-3 # 

res_pareto_pdf_n_11_method3 = rep(0,B)
res_pareto_pdf_n_21_method3 = rep(0,B)
res_pareto_pdf_n_31_method3 = rep(0,B)
res_pareto_pdf_n_41_method3 = rep(0,B)
res_pareto_pdf_n_51_method3 = rep(0,B)

res_pareto_pdf_n_12_method3 = rep(0,B)
res_pareto_pdf_n_22_method3 = rep(0,B)
res_pareto_pdf_n_32_method3 = rep(0,B)
res_pareto_pdf_n_42_method3 = rep(0,B)
res_pareto_pdf_n_52_method3 = rep(0,B)

# EM # 

res_pareto_pdf_n_11_EM = rep(0,B)
res_pareto_pdf_n_21_EM = rep(0,B)
res_pareto_pdf_n_31_EM = rep(0,B)
res_pareto_pdf_n_41_EM = rep(0,B)
res_pareto_pdf_n_51_EM = rep(0,B)

res_pareto_pdf_n_12_EM = rep(0,B)
res_pareto_pdf_n_22_EM = rep(0,B)
res_pareto_pdf_n_32_EM = rep(0,B)
res_pareto_pdf_n_42_EM = rep(0,B)
res_pareto_pdf_n_52_EM = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  # Pareto Distribution - n1
  
  x = sort(rpareto(n[1],location,shape_p))
  n_i <- hist(x, breaks = pareto_grid, plot = F)$counts
  
  res_pareto_pdf_n_11_method1[i] = L1_Distance_method1(n_i,pareto_grid,pdf_pareto,range_pareto)
  res_pareto_pdf_n_12_method1[i] = L2_Distance_method1(n_i,pareto_grid,pdf_pareto,range_pareto)
  
  res_pareto_pdf_n_11_method2[i] = L1_Distance_method2(n_i,pareto_grid,pdf_pareto,range_pareto)
  res_pareto_pdf_n_12_method2[i] = L2_Distance_method2(n_i,pareto_grid,pdf_pareto,range_pareto)
  
  res_pareto_pdf_n_11_method3[i] = L1_Distance_method3(n_i,pareto_grid,pdf_pareto,range_pareto)
  res_pareto_pdf_n_12_method3[i] = L2_Distance_method3(n_i,pareto_grid,pdf_pareto,range_pareto)
  
  res_pareto_pdf_n_11_EM[i] = EM_L1_L2(n_i,pareto_grid,pdf_pareto,range_pareto)$L1
  res_pareto_pdf_n_12_EM[i] = EM_L1_L2(n_i,pareto_grid,pdf_pareto,range_pareto)$L2
  
  # Pareto Distribution - n2
  
  x = sort(rpareto(n[2],location,shape_p))
  n_i <- hist(x, breaks = pareto_grid, plot = F)$counts
  
  res_pareto_pdf_n_21_method1[i] = L1_Distance_method1(n_i,pareto_grid,pdf_pareto,range_pareto)
  res_pareto_pdf_n_22_method1[i] = L2_Distance_method1(n_i,pareto_grid,pdf_pareto,range_pareto)
  
  res_pareto_pdf_n_21_method2[i] = L1_Distance_method2(n_i,pareto_grid,pdf_pareto,range_pareto)
  res_pareto_pdf_n_22_method2[i] = L2_Distance_method2(n_i,pareto_grid,pdf_pareto,range_pareto)
  
  res_pareto_pdf_n_21_method3[i] = L1_Distance_method3(n_i,pareto_grid,pdf_pareto,range_pareto)
  res_pareto_pdf_n_22_method3[i] = L2_Distance_method3(n_i,pareto_grid,pdf_pareto,range_pareto)
  
  res_pareto_pdf_n_21_EM[i] = EM_L1_L2(n_i,pareto_grid,pdf_pareto,range_pareto)$L1
  res_pareto_pdf_n_22_EM[i] = EM_L1_L2(n_i,pareto_grid,pdf_pareto,range_pareto)$L2
  
  # Pareto Distribution - n3
  
  x = sort(rpareto(n[3],location,shape_p))
  n_i <- hist(x, breaks = pareto_grid, plot = F)$counts
  
  res_pareto_pdf_n_31_method1[i] = L1_Distance_method1(n_i,pareto_grid,pdf_pareto,range_pareto)
  res_pareto_pdf_n_32_method1[i] = L2_Distance_method1(n_i,pareto_grid,pdf_pareto,range_pareto)
  
  res_pareto_pdf_n_31_method2[i] = L1_Distance_method2(n_i,pareto_grid,pdf_pareto,range_pareto)
  res_pareto_pdf_n_32_method2[i] = L2_Distance_method2(n_i,pareto_grid,pdf_pareto,range_pareto)
  
  res_pareto_pdf_n_31_method3[i] = L1_Distance_method3(n_i,pareto_grid,pdf_pareto,range_pareto)
  res_pareto_pdf_n_32_method3[i] = L2_Distance_method3(n_i,pareto_grid,pdf_pareto,range_pareto)
  
  res_pareto_pdf_n_31_EM[i] = EM_L1_L2(n_i,pareto_grid,pdf_pareto,range_pareto)$L1
  res_pareto_pdf_n_32_EM[i] = EM_L1_L2(n_i,pareto_grid,pdf_pareto,range_pareto)$L2
  
  # Pareto Distribution - n4
  
  x = sort(rpareto(n[4],location,shape_p))
  n_i <- hist(x, breaks = pareto_grid, plot = F)$counts
  
  res_pareto_pdf_n_41_method1[i] = L1_Distance_method1(n_i,pareto_grid,pdf_pareto,range_pareto)
  res_pareto_pdf_n_42_method1[i] = L2_Distance_method1(n_i,pareto_grid,pdf_pareto,range_pareto)
  
  res_pareto_pdf_n_41_method2[i] = L1_Distance_method2(n_i,pareto_grid,pdf_pareto,range_pareto)
  res_pareto_pdf_n_42_method2[i] = L2_Distance_method2(n_i,pareto_grid,pdf_pareto,range_pareto)
  
  res_pareto_pdf_n_41_method3[i] = L1_Distance_method3(n_i,pareto_grid,pdf_pareto,range_pareto)
  res_pareto_pdf_n_42_method3[i] = L2_Distance_method3(n_i,pareto_grid,pdf_pareto,range_pareto)
  
  res_pareto_pdf_n_41_EM[i] = EM_L1_L2(n_i,pareto_grid,pdf_pareto,range_pareto)$L1
  res_pareto_pdf_n_42_EM[i] = EM_L1_L2(n_i,pareto_grid,pdf_pareto,range_pareto)$L2
  
  # Pareto Distribution - n5
  
  x = sort(rpareto(n[5],location,shape_p))
  n_i <- hist(x, breaks = pareto_grid, plot = F)$counts
  
  res_pareto_pdf_n_51_method1[i] = L1_Distance_method1(n_i,pareto_grid,pdf_pareto,range_pareto)
  res_pareto_pdf_n_52_method1[i] = L2_Distance_method1(n_i,pareto_grid,pdf_pareto,range_pareto)
  
  res_pareto_pdf_n_51_method2[i] = L1_Distance_method2(n_i,pareto_grid,pdf_pareto,range_pareto)
  res_pareto_pdf_n_52_method2[i] = L2_Distance_method2(n_i,pareto_grid,pdf_pareto,range_pareto)
  
  res_pareto_pdf_n_51_method3[i] = L1_Distance_method3(n_i,pareto_grid,pdf_pareto,range_pareto)
  res_pareto_pdf_n_52_method3[i] = L2_Distance_method3(n_i,pareto_grid,pdf_pareto,range_pareto)
  
  res_pareto_pdf_n_51_EM[i] = EM_L1_L2(n_i,pareto_grid,pdf_pareto,range_pareto)$L1
  res_pareto_pdf_n_52_EM[i] = EM_L1_L2(n_i,pareto_grid,pdf_pareto,range_pareto)$L2
  
}

# Method-1 #  

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L1-n/L1-PDF-2-Pareto-n-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_pareto_pdf_n_11_method1,res_pareto_pdf_n_21_method1,res_pareto_pdf_n_31_method1,res_pareto_pdf_n_41_method1,res_pareto_pdf_n_51_method1,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors5,ylim=c(0,1.6))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L2-n/L2-PDF-2-Pareto-n-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_pareto_pdf_n_12_method1,res_pareto_pdf_n_22_method1,res_pareto_pdf_n_32_method1,res_pareto_pdf_n_42_method1,res_pareto_pdf_n_52_method1,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors5,ylim=c(0,1.6))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

# Method-2 #  

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L1-n/L1-PDF-2-Pareto-n-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_pareto_pdf_n_11_method2,res_pareto_pdf_n_21_method2,res_pareto_pdf_n_31_method2,res_pareto_pdf_n_41_method2,res_pareto_pdf_n_51_method2,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors5,ylim=c(0,1.6))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L2-n/L2-PDF-2-Pareto-n-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_pareto_pdf_n_12_method2,res_pareto_pdf_n_22_method2,res_pareto_pdf_n_32_method2,res_pareto_pdf_n_42_method2,res_pareto_pdf_n_52_method2,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors5,ylim=c(0,1.6))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

# Method-3 #  

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L1-n/L1-PDF-2-Pareto-n-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_pareto_pdf_n_11_method3,res_pareto_pdf_n_21_method3,res_pareto_pdf_n_31_method3,res_pareto_pdf_n_41_method3,res_pareto_pdf_n_51_method3,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors5,ylim=c(0,1.6))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L2-n/L2-PDF-2-Pareto-n-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_pareto_pdf_n_12_method3,res_pareto_pdf_n_22_method3,res_pareto_pdf_n_32_method3,res_pareto_pdf_n_42_method3,res_pareto_pdf_n_52_method3,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors5,ylim=c(0,1.6))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

# EM #  

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L1-n/L1-PDF-2-Pareto-n-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_pareto_pdf_n_11_EM,res_pareto_pdf_n_21_EM,res_pareto_pdf_n_31_EM,res_pareto_pdf_n_41_EM,res_pareto_pdf_n_51_EM,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors5,ylim=c(0,1.6))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L2-n/L2-PDF-2-Pareto-n-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_pareto_pdf_n_12_EM,res_pareto_pdf_n_22_EM,res_pareto_pdf_n_32_EM,res_pareto_pdf_n_42_EM,res_pareto_pdf_n_52_EM,
        names = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors5,ylim=c(0,1.6))
#segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
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

# Method-1 # 

res_lap_pdf_gw_11_method1 = rep(0,B)
res_lap_pdf_gw_21_method1 = rep(0,B)
res_lap_pdf_gw_31_method1 = rep(0,B)
res_lap_pdf_gw_41_method1 = rep(0,B)
res_lap_pdf_gw_51_method1 = rep(0,B)

res_lap_pdf_gw_12_method1 = rep(0,B)
res_lap_pdf_gw_22_method1 = rep(0,B)
res_lap_pdf_gw_32_method1 = rep(0,B)
res_lap_pdf_gw_42_method1 = rep(0,B)
res_lap_pdf_gw_52_method1 = rep(0,B)

# Method-2 # 

res_lap_pdf_gw_11_method2 = rep(0,B)
res_lap_pdf_gw_21_method2 = rep(0,B)
res_lap_pdf_gw_31_method2 = rep(0,B)
res_lap_pdf_gw_41_method2 = rep(0,B)
res_lap_pdf_gw_51_method2 = rep(0,B)

res_lap_pdf_gw_12_method2 = rep(0,B)
res_lap_pdf_gw_22_method2 = rep(0,B)
res_lap_pdf_gw_32_method2 = rep(0,B)
res_lap_pdf_gw_42_method2 = rep(0,B)
res_lap_pdf_gw_52_method2 = rep(0,B)

# Method-3 # 

res_lap_pdf_gw_11_method3 = rep(0,B)
res_lap_pdf_gw_21_method3 = rep(0,B)
res_lap_pdf_gw_31_method3 = rep(0,B)
res_lap_pdf_gw_41_method3 = rep(0,B)
res_lap_pdf_gw_51_method3 = rep(0,B)

res_lap_pdf_gw_12_method3 = rep(0,B)
res_lap_pdf_gw_22_method3 = rep(0,B)
res_lap_pdf_gw_32_method3 = rep(0,B)
res_lap_pdf_gw_42_method3 = rep(0,B)
res_lap_pdf_gw_52_method3 = rep(0,B)

# EM # 

res_lap_pdf_gw_11_EM = rep(0,B)
res_lap_pdf_gw_21_EM = rep(0,B)
res_lap_pdf_gw_31_EM = rep(0,B)
res_lap_pdf_gw_41_EM = rep(0,B)
res_lap_pdf_gw_51_EM = rep(0,B)

res_lap_pdf_gw_12_EM = rep(0,B)
res_lap_pdf_gw_22_EM = rep(0,B)
res_lap_pdf_gw_32_EM = rep(0,B)
res_lap_pdf_gw_42_EM = rep(0,B)
res_lap_pdf_gw_52_EM = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  x <- sort(rlaplace(n, mean = mu,sd=sigma))
  
  # Laplace Distribution - GW1
  
  n_i <- hist(x, breaks = laplace_grids[[1]], plot = F)$counts
  
  res_lap_pdf_gw_11_method1[i] = L1_Distance_method1(n_i,laplace_grids[[1]],pdf_laplace,range_laplace)
  res_lap_pdf_gw_12_method1[i] = L2_Distance_method1(n_i,laplace_grids[[1]],pdf_laplace,range_laplace)
  
  res_lap_pdf_gw_11_method2[i] = L1_Distance_method2(n_i,laplace_grids[[1]],pdf_laplace,range_laplace)
  res_lap_pdf_gw_12_method2[i] = L2_Distance_method2(n_i,laplace_grids[[1]],pdf_laplace,range_laplace)
  
  res_lap_pdf_gw_11_method3[i] = L1_Distance_method3(n_i,laplace_grids[[1]],pdf_laplace,range_laplace)
  res_lap_pdf_gw_12_method3[i] = L2_Distance_method3(n_i,laplace_grids[[1]],pdf_laplace,range_laplace)
  
  res_lap_pdf_gw_11_EM[i] = EM_L1_L2(n_i,laplace_grids[[1]],pdf_laplace,range_laplace)$L1
  res_lap_pdf_gw_12_EM[i] = EM_L1_L2(n_i,laplace_grids[[1]],pdf_laplace,range_laplace)$L2
  
  # Laplace Distribution - GW2
  
  n_i <- hist(x, breaks = laplace_grids[[2]], plot = F)$counts
  
  res_lap_pdf_gw_21_method1[i] = L1_Distance_method1(n_i,laplace_grids[[2]],pdf_laplace,range_laplace)
  res_lap_pdf_gw_22_method1[i] = L2_Distance_method1(n_i,laplace_grids[[2]],pdf_laplace,range_laplace)
  
  res_lap_pdf_gw_21_method2[i] = L1_Distance_method2(n_i,laplace_grids[[2]],pdf_laplace,range_laplace)
  res_lap_pdf_gw_22_method2[i] = L2_Distance_method2(n_i,laplace_grids[[2]],pdf_laplace,range_laplace)
  
  res_lap_pdf_gw_21_method3[i] = L1_Distance_method3(n_i,laplace_grids[[2]],pdf_laplace,range_laplace)
  res_lap_pdf_gw_22_method3[i] = L2_Distance_method3(n_i,laplace_grids[[2]],pdf_laplace,range_laplace)
  
  res_lap_pdf_gw_21_EM[i] = EM_L1_L2(n_i,laplace_grids[[2]],pdf_laplace,range_laplace)$L1
  res_lap_pdf_gw_22_EM[i] = EM_L1_L2(n_i,laplace_grids[[2]],pdf_laplace,range_laplace)$L2
  
  # Laplace Distribution - GW3
  
  n_i <- hist(x, breaks = laplace_grids[[3]], plot = F)$counts
  
  res_lap_pdf_gw_31_method1[i] = L1_Distance_method1(n_i,laplace_grids[[3]],pdf_laplace,range_laplace)
  res_lap_pdf_gw_32_method1[i] = L2_Distance_method1(n_i,laplace_grids[[3]],pdf_laplace,range_laplace)
  
  res_lap_pdf_gw_31_method2[i] = L1_Distance_method2(n_i,laplace_grids[[3]],pdf_laplace,range_laplace)
  res_lap_pdf_gw_32_method2[i] = L2_Distance_method2(n_i,laplace_grids[[3]],pdf_laplace,range_laplace)
  
  res_lap_pdf_gw_31_method3[i] = L1_Distance_method3(n_i,laplace_grids[[3]],pdf_laplace,range_laplace)
  res_lap_pdf_gw_32_method3[i] = L2_Distance_method3(n_i,laplace_grids[[3]],pdf_laplace,range_laplace)
  
  res_lap_pdf_gw_31_EM[i] = EM_L1_L2(n_i,laplace_grids[[3]],pdf_laplace,range_laplace)$L1
  res_lap_pdf_gw_32_EM[i] = EM_L1_L2(n_i,laplace_grids[[3]],pdf_laplace,range_laplace)$L2
  
  # Laplace Distribution - GW4
  
  n_i <- hist(x, breaks = laplace_grids[[4]], plot = F)$counts
  
  res_lap_pdf_gw_41_method1[i] = L1_Distance_method1(n_i,laplace_grids[[4]],pdf_laplace,range_laplace)
  res_lap_pdf_gw_42_method1[i] = L2_Distance_method1(n_i,laplace_grids[[4]],pdf_laplace,range_laplace)
  
  res_lap_pdf_gw_41_method2[i] = L1_Distance_method2(n_i,laplace_grids[[4]],pdf_laplace,range_laplace)
  res_lap_pdf_gw_42_method2[i] = L2_Distance_method2(n_i,laplace_grids[[4]],pdf_laplace,range_laplace)
  
  res_lap_pdf_gw_41_method3[i] = L1_Distance_method3(n_i,laplace_grids[[4]],pdf_laplace,range_laplace)
  res_lap_pdf_gw_42_method3[i] = L2_Distance_method3(n_i,laplace_grids[[4]],pdf_laplace,range_laplace)
  
  res_lap_pdf_gw_41_EM[i] = EM_L1_L2(n_i,laplace_grids[[4]],pdf_laplace,range_laplace)$L1
  res_lap_pdf_gw_42_EM[i] = EM_L1_L2(n_i,laplace_grids[[4]],pdf_laplace,range_laplace)$L2
  
  # Laplace Distribution - GW5
  
  n_i <- hist(x, breaks = laplace_grids[[5]], plot = F)$counts
  
  res_lap_pdf_gw_51_method1[i] = L1_Distance_method1(n_i,laplace_grids[[5]],pdf_laplace,range_laplace)
  res_lap_pdf_gw_52_method1[i] = L2_Distance_method1(n_i,laplace_grids[[5]],pdf_laplace,range_laplace)
  
  res_lap_pdf_gw_51_method2[i] = L1_Distance_method2(n_i,laplace_grids[[5]],pdf_laplace,range_laplace)
  res_lap_pdf_gw_52_method2[i] = L2_Distance_method2(n_i,laplace_grids[[5]],pdf_laplace,range_laplace)
  
  res_lap_pdf_gw_51_method3[i] = L1_Distance_method3(n_i,laplace_grids[[5]],pdf_laplace,range_laplace)
  res_lap_pdf_gw_52_method3[i] = L2_Distance_method3(n_i,laplace_grids[[5]],pdf_laplace,range_laplace)
  
  res_lap_pdf_gw_51_EM[i] = EM_L1_L2(n_i,laplace_grids[[5]],pdf_laplace,range_laplace)$L1
  res_lap_pdf_gw_52_EM[i] = EM_L1_L2(n_i,laplace_grids[[5]],pdf_laplace,range_laplace)$L2
  
}

# Method-1 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L1-GW/L1-PDF-2-Laplace-gw-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_lap_pdf_gw_11_method1,res_lap_pdf_gw_21_method1,res_lap_pdf_gw_31_method1,res_lap_pdf_gw_41_method1,res_lap_pdf_gw_51_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,1.7))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L2-GW/L2-PDF-2-Laplace-gw-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_lap_pdf_gw_12_method1,res_lap_pdf_gw_22_method1,res_lap_pdf_gw_32_method1,res_lap_pdf_gw_42_method1,res_lap_pdf_gw_52_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,1.2))
dev.off()

# Method-2 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L1-GW/L1-PDF-2-Laplace-gw-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_lap_pdf_gw_11_method1,res_lap_pdf_gw_21_method1,res_lap_pdf_gw_31_method1,res_lap_pdf_gw_41_method1,res_lap_pdf_gw_51_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,1.7))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L2-GW/L2-PDF-2-Laplace-gw-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_lap_pdf_gw_12_method1,res_lap_pdf_gw_22_method1,res_lap_pdf_gw_32_method1,res_lap_pdf_gw_42_method1,res_lap_pdf_gw_52_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,1.2))
dev.off()

# Method-3 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L1-GW/L1-PDF-2-Laplace-gw-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_lap_pdf_gw_11_method3,res_lap_pdf_gw_21_method3,res_lap_pdf_gw_31_method3,res_lap_pdf_gw_41_method3,res_lap_pdf_gw_51_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,1.7))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L2-GW/L2-PDF-2-Laplace-gw-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_lap_pdf_gw_12_method3,res_lap_pdf_gw_22_method3,res_lap_pdf_gw_32_method3,res_lap_pdf_gw_42_method3,res_lap_pdf_gw_52_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,1.2))
dev.off()

# EM # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L1-GW/L1-PDF-2-Laplace-gw-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_lap_pdf_gw_11_EM,res_lap_pdf_gw_21_EM,res_lap_pdf_gw_31_EM,res_lap_pdf_gw_41_EM,res_lap_pdf_gw_51_EM,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,1.7))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L2-GW/L2-PDF-2-Laplace-gw-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_lap_pdf_gw_12_EM,res_lap_pdf_gw_22_EM,res_lap_pdf_gw_32_EM,res_lap_pdf_gw_42_EM,res_lap_pdf_gw_52_EM,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,1.2))
dev.off()

# Chi-Square Distribution - Grid Width Simulation 

MM = qchisq(0.9999999,df) + 3
chisqr_grids = list(seq(0,MM, by = chisqr_deltas[5]),
                    seq(0,MM, by = chisqr_deltas[4]),
                    seq(0,MM, by = chisqr_deltas[3]),
                    seq(0,MM, by = chisqr_deltas[2]),
                    seq(0,MM, by = chisqr_deltas[1]))

# Method-1 # 

res_chisq_pdf_gw_11_method1 = rep(0,B)
res_chisq_pdf_gw_21_method1 = rep(0,B)
res_chisq_pdf_gw_31_method1 = rep(0,B)
res_chisq_pdf_gw_41_method1 = rep(0,B)
res_chisq_pdf_gw_51_method1 = rep(0,B)

res_chisq_pdf_gw_12_method1 = rep(0,B)
res_chisq_pdf_gw_22_method1 = rep(0,B)
res_chisq_pdf_gw_32_method1 = rep(0,B)
res_chisq_pdf_gw_42_method1 = rep(0,B)
res_chisq_pdf_gw_52_method1 = rep(0,B)

# Method-2 # 

res_chisq_pdf_gw_11_method2 = rep(0,B)
res_chisq_pdf_gw_21_method2 = rep(0,B)
res_chisq_pdf_gw_31_method2 = rep(0,B)
res_chisq_pdf_gw_41_method2 = rep(0,B)
res_chisq_pdf_gw_51_method2 = rep(0,B)

res_chisq_pdf_gw_12_method2 = rep(0,B)
res_chisq_pdf_gw_22_method2 = rep(0,B)
res_chisq_pdf_gw_32_method2 = rep(0,B)
res_chisq_pdf_gw_42_method2 = rep(0,B)
res_chisq_pdf_gw_52_method2 = rep(0,B)

# Method-3 # 

res_chisq_pdf_gw_11_method3 = rep(0,B)
res_chisq_pdf_gw_21_method3 = rep(0,B)
res_chisq_pdf_gw_31_method3 = rep(0,B)
res_chisq_pdf_gw_41_method3 = rep(0,B)
res_chisq_pdf_gw_51_method3 = rep(0,B)

res_chisq_pdf_gw_12_method3 = rep(0,B)
res_chisq_pdf_gw_22_method3 = rep(0,B)
res_chisq_pdf_gw_32_method3 = rep(0,B)
res_chisq_pdf_gw_42_method3 = rep(0,B)
res_chisq_pdf_gw_52_method3 = rep(0,B)

# EM # 

res_chisq_pdf_gw_11_EM = rep(0,B)
res_chisq_pdf_gw_21_EM = rep(0,B)
res_chisq_pdf_gw_31_EM = rep(0,B)
res_chisq_pdf_gw_41_EM = rep(0,B)
res_chisq_pdf_gw_51_EM = rep(0,B)

res_chisq_pdf_gw_12_EM = rep(0,B)
res_chisq_pdf_gw_22_EM = rep(0,B)
res_chisq_pdf_gw_32_EM = rep(0,B)
res_chisq_pdf_gw_42_EM = rep(0,B)
res_chisq_pdf_gw_52_EM = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  x <- sort(rchisq(n,df))
  
  # Chi-Square Distribution - GW1
  
  n_i <- hist(x, breaks = chisqr_grids[[1]], plot = F)$counts
  
  res_chisq_pdf_gw_11_method1[i] = L1_Distance_method1(n_i,chisqr_grids[[1]],pdf_chisq,range_chisq)
  res_chisq_pdf_gw_12_method1[i] = L2_Distance_method1(n_i,chisqr_grids[[1]],pdf_chisq,range_chisq)
  
  res_chisq_pdf_gw_11_method2[i] = L1_Distance_method2(n_i,chisqr_grids[[1]],pdf_chisq,range_chisq)
  res_chisq_pdf_gw_12_method2[i] = L2_Distance_method2(n_i,chisqr_grids[[1]],pdf_chisq,range_chisq)
  
  res_chisq_pdf_gw_11_method3[i] = L1_Distance_method3(n_i,chisqr_grids[[1]],pdf_chisq,range_chisq)
  res_chisq_pdf_gw_12_method3[i] = L2_Distance_method3(n_i,chisqr_grids[[1]],pdf_chisq,range_chisq)
  
  res_chisq_pdf_gw_11_EM[i] = EM_L1_L2(n_i,chisqr_grids[[1]],pdf_chisq,range_chisq)$L1
  res_chisq_pdf_gw_12_EM[i] = EM_L1_L2(n_i,chisqr_grids[[1]],pdf_chisq,range_chisq)$L2
  
  # Chi-Square Distribution - GW2
  
  n_i <- hist(x, breaks = chisqr_grids[[2]], plot = F)$counts
  
  res_chisq_pdf_gw_21_method1[i] = L1_Distance_method1(n_i,chisqr_grids[[2]],pdf_chisq,range_chisq)
  res_chisq_pdf_gw_22_method1[i] = L2_Distance_method1(n_i,chisqr_grids[[2]],pdf_chisq,range_chisq)
  
  res_chisq_pdf_gw_21_method2[i] = L1_Distance_method2(n_i,chisqr_grids[[2]],pdf_chisq,range_chisq)
  res_chisq_pdf_gw_22_method2[i] = L2_Distance_method2(n_i,chisqr_grids[[2]],pdf_chisq,range_chisq)
  
  res_chisq_pdf_gw_21_method3[i] = L1_Distance_method3(n_i,chisqr_grids[[2]],pdf_chisq,range_chisq)
  res_chisq_pdf_gw_22_method3[i] = L2_Distance_method3(n_i,chisqr_grids[[2]],pdf_chisq,range_chisq)
  
  res_chisq_pdf_gw_21_EM[i] = EM_L1_L2(n_i,chisqr_grids[[2]],pdf_chisq,range_chisq)$L1
  res_chisq_pdf_gw_22_EM[i] = EM_L1_L2(n_i,chisqr_grids[[2]],pdf_chisq,range_chisq)$L2
  
  # Chi-Square Distribution - GW3
  
  n_i <- hist(x, breaks = chisqr_grids[[3]], plot = F)$counts
  
  res_chisq_pdf_gw_31_method1[i] = L1_Distance_method1(n_i,chisqr_grids[[3]],pdf_chisq,range_chisq)
  res_chisq_pdf_gw_32_method1[i] = L2_Distance_method1(n_i,chisqr_grids[[3]],pdf_chisq,range_chisq)
  
  res_chisq_pdf_gw_31_method2[i] = L1_Distance_method2(n_i,chisqr_grids[[3]],pdf_chisq,range_chisq)
  res_chisq_pdf_gw_32_method2[i] = L2_Distance_method2(n_i,chisqr_grids[[3]],pdf_chisq,range_chisq)
  
  res_chisq_pdf_gw_31_method3[i] = L1_Distance_method3(n_i,chisqr_grids[[3]],pdf_chisq,range_chisq)
  res_chisq_pdf_gw_32_method3[i] = L2_Distance_method3(n_i,chisqr_grids[[3]],pdf_chisq,range_chisq)
  
  res_chisq_pdf_gw_31_EM[i] = EM_L1_L2(n_i,chisqr_grids[[3]],pdf_chisq,range_chisq)$L1
  res_chisq_pdf_gw_32_EM[i] = EM_L1_L2(n_i,chisqr_grids[[3]],pdf_chisq,range_chisq)$L2
  
  # Chi-Square Distribution - GW4
  
  n_i <- hist(x, breaks = chisqr_grids[[4]], plot = F)$counts
  
  res_chisq_pdf_gw_41_method1[i] = L1_Distance_method1(n_i,chisqr_grids[[4]],pdf_chisq,range_chisq)
  res_chisq_pdf_gw_42_method1[i] = L2_Distance_method1(n_i,chisqr_grids[[4]],pdf_chisq,range_chisq)
  
  res_chisq_pdf_gw_41_method2[i] = L1_Distance_method2(n_i,chisqr_grids[[4]],pdf_chisq,range_chisq)
  res_chisq_pdf_gw_42_method2[i] = L2_Distance_method2(n_i,chisqr_grids[[4]],pdf_chisq,range_chisq)
  
  res_chisq_pdf_gw_41_method3[i] = L1_Distance_method3(n_i,chisqr_grids[[4]],pdf_chisq,range_chisq)
  res_chisq_pdf_gw_42_method3[i] = L2_Distance_method3(n_i,chisqr_grids[[4]],pdf_chisq,range_chisq)
  
  res_chisq_pdf_gw_41_EM[i] = EM_L1_L2(n_i,chisqr_grids[[4]],pdf_chisq,range_chisq)$L1
  res_chisq_pdf_gw_42_EM[i] = EM_L1_L2(n_i,chisqr_grids[[4]],pdf_chisq,range_chisq)$L2
  
  # Chi-Square Distribution - GW5
  
  n_i <- hist(x, breaks = chisqr_grids[[5]], plot = F)$counts
  
  res_chisq_pdf_gw_51_method1[i] = L1_Distance_method1(n_i,chisqr_grids[[5]],pdf_chisq,range_chisq)
  res_chisq_pdf_gw_52_method1[i] = L2_Distance_method1(n_i,chisqr_grids[[5]],pdf_chisq,range_chisq)

  res_chisq_pdf_gw_51_method2[i] = L1_Distance_method2(n_i,chisqr_grids[[5]],pdf_chisq,range_chisq)
  res_chisq_pdf_gw_52_method2[i] = L2_Distance_method2(n_i,chisqr_grids[[5]],pdf_chisq,range_chisq)
  
  res_chisq_pdf_gw_51_method3[i] = L1_Distance_method3(n_i,chisqr_grids[[5]],pdf_chisq,range_chisq)
  res_chisq_pdf_gw_52_method3[i] = L2_Distance_method3(n_i,chisqr_grids[[5]],pdf_chisq,range_chisq)
  
  res_chisq_pdf_gw_51_EM[i] = EM_L1_L2(n_i,chisqr_grids[[5]],pdf_chisq,range_chisq)$L1
  res_chisq_pdf_gw_52_EM[i] = EM_L1_L2(n_i,chisqr_grids[[5]],pdf_chisq,range_chisq)$L2
  
}

# Method-1 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L1-GW/L1-PDF-2-Chisq-gw-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_chisq_pdf_gw_11_method1,res_chisq_pdf_gw_21_method1,res_chisq_pdf_gw_31_method1,res_chisq_pdf_gw_41_method1,res_chisq_pdf_gw_51_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.9))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L2-GW/L2-PDF-2-Chisq-gw-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_chisq_pdf_gw_12_method1,res_chisq_pdf_gw_22_method1,res_chisq_pdf_gw_32_method1,res_chisq_pdf_gw_42_method1,res_chisq_pdf_gw_52_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.5))
dev.off()

# Method-2 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L1-GW/L1-PDF-2-Chisq-gw-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_chisq_pdf_gw_11_method2,res_chisq_pdf_gw_21_method2,res_chisq_pdf_gw_31_method2,res_chisq_pdf_gw_41_method2,res_chisq_pdf_gw_51_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.9))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L2-GW/L2-PDF-2-Chisq-gw-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_chisq_pdf_gw_12_method2,res_chisq_pdf_gw_22_method2,res_chisq_pdf_gw_32_method2,res_chisq_pdf_gw_42_method2,res_chisq_pdf_gw_52_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.5))
dev.off()

# Method-3 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L1-GW/L1-PDF-2-Chisq-gw-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_chisq_pdf_gw_11_method3,res_chisq_pdf_gw_21_method3,res_chisq_pdf_gw_31_method3,res_chisq_pdf_gw_41_method3,res_chisq_pdf_gw_51_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.9))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L2-GW/L2-PDF-2-Chisq-gw-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_chisq_pdf_gw_12_method3,res_chisq_pdf_gw_22_method3,res_chisq_pdf_gw_32_method3,res_chisq_pdf_gw_42_method3,res_chisq_pdf_gw_52_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.5))
dev.off()

# EM # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L1-GW/L1-PDF-2-Chisq-gw-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_chisq_pdf_gw_11_EM,res_chisq_pdf_gw_21_EM,res_chisq_pdf_gw_31_EM,res_chisq_pdf_gw_41_EM,res_chisq_pdf_gw_51_EM,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.9))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L2-GW/L2-PDF-2-Chisq-gw-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_chisq_pdf_gw_12_EM,res_chisq_pdf_gw_22_EM,res_chisq_pdf_gw_32_EM,res_chisq_pdf_gw_42_EM,res_chisq_pdf_gw_52_EM,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.5))
dev.off()

# Log-Normal Distribution - Grid Width Simulation 

MM = qlnorm(0.99999999999,mu,sigma) + 30
lnorm_grids = list(seq(0,MM, by = lnorm_deltas[5]),
                   seq(0,MM, by = lnorm_deltas[4]),
                   seq(0,MM, by = lnorm_deltas[3]),
                   seq(0,MM, by = lnorm_deltas[2]),
                   seq(0,MM, by = lnorm_deltas[1]))

# Method-1 # 

res_lnorm_pdf_gw_11_method1 = rep(0,B)
res_lnorm_pdf_gw_21_method1 = rep(0,B)
res_lnorm_pdf_gw_31_method1 = rep(0,B)
res_lnorm_pdf_gw_41_method1 = rep(0,B)
res_lnorm_pdf_gw_51_method1 = rep(0,B)

res_lnorm_pdf_gw_12_method1 = rep(0,B)
res_lnorm_pdf_gw_22_method1 = rep(0,B)
res_lnorm_pdf_gw_32_method1 = rep(0,B)
res_lnorm_pdf_gw_42_method1 = rep(0,B)
res_lnorm_pdf_gw_52_method1 = rep(0,B)

# Method-2 # 

res_lnorm_pdf_gw_11_method2 = rep(0,B)
res_lnorm_pdf_gw_21_method2 = rep(0,B)
res_lnorm_pdf_gw_31_method2 = rep(0,B)
res_lnorm_pdf_gw_41_method2 = rep(0,B)
res_lnorm_pdf_gw_51_method2 = rep(0,B)

res_lnorm_pdf_gw_12_method2 = rep(0,B)
res_lnorm_pdf_gw_22_method2 = rep(0,B)
res_lnorm_pdf_gw_32_method2 = rep(0,B)
res_lnorm_pdf_gw_42_method2 = rep(0,B)
res_lnorm_pdf_gw_52_method2 = rep(0,B)

# Method-3 # 

res_lnorm_pdf_gw_11_method3 = rep(0,B)
res_lnorm_pdf_gw_21_method3 = rep(0,B)
res_lnorm_pdf_gw_31_method3 = rep(0,B)
res_lnorm_pdf_gw_41_method3 = rep(0,B)
res_lnorm_pdf_gw_51_method3 = rep(0,B)

res_lnorm_pdf_gw_12_method3 = rep(0,B)
res_lnorm_pdf_gw_22_method3 = rep(0,B)
res_lnorm_pdf_gw_32_method3 = rep(0,B)
res_lnorm_pdf_gw_42_method3 = rep(0,B)
res_lnorm_pdf_gw_52_method3 = rep(0,B)

# EM # 

res_lnorm_pdf_gw_11_EM = rep(0,B)
res_lnorm_pdf_gw_21_EM = rep(0,B)
res_lnorm_pdf_gw_31_EM = rep(0,B)
res_lnorm_pdf_gw_41_EM = rep(0,B)
res_lnorm_pdf_gw_51_EM = rep(0,B)

res_lnorm_pdf_gw_12_EM = rep(0,B)
res_lnorm_pdf_gw_22_EM = rep(0,B)
res_lnorm_pdf_gw_32_EM = rep(0,B)
res_lnorm_pdf_gw_42_EM = rep(0,B)
res_lnorm_pdf_gw_52_EM = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  x <- sort(rlnorm(n,mul,sigma))
  
  # Log-Normal Distribution - GW1
  
  n_i <- hist(x, breaks = lnorm_grids[[1]], plot = F)$counts
  
  res_lnorm_pdf_gw_11_method1[i] = L1_Distance_method1(n_i,lnorm_grids[[1]],pdf_lnorm,range_lnorm)
  res_lnorm_pdf_gw_12_method1[i] = L2_Distance_method1(n_i,lnorm_grids[[1]],pdf_lnorm,range_lnorm)
  
  res_lnorm_pdf_gw_11_method2[i] = L1_Distance_method2(n_i,lnorm_grids[[1]],pdf_lnorm,range_lnorm)
  res_lnorm_pdf_gw_12_method2[i] = L2_Distance_method2(n_i,lnorm_grids[[1]],pdf_lnorm,range_lnorm)
  
  res_lnorm_pdf_gw_11_method3[i] = L1_Distance_method3(n_i,lnorm_grids[[1]],pdf_lnorm,range_lnorm)
  res_lnorm_pdf_gw_12_method3[i] = L2_Distance_method3(n_i,lnorm_grids[[1]],pdf_lnorm,range_lnorm)
  
  res_lnorm_pdf_gw_11_EM[i] = EM_L1_L2(n_i,lnorm_grids[[1]],pdf_lnorm,range_lnorm)$L1
  res_lnorm_pdf_gw_12_EM[i] = EM_L1_L2(n_i,lnorm_grids[[1]],pdf_lnorm,range_lnorm)$L2
  
  # Log-Normal Distribution - GW2
  
  n_i <- hist(x, breaks = lnorm_grids[[2]], plot = F)$counts
  
  res_lnorm_pdf_gw_21_method1[i] = L1_Distance_method1(n_i,lnorm_grids[[2]],pdf_lnorm,range_lnorm)
  res_lnorm_pdf_gw_22_method1[i] = L2_Distance_method1(n_i,lnorm_grids[[2]],pdf_lnorm,range_lnorm)
  
  res_lnorm_pdf_gw_21_method2[i] = L1_Distance_method2(n_i,lnorm_grids[[2]],pdf_lnorm,range_lnorm)
  res_lnorm_pdf_gw_22_method2[i] = L2_Distance_method2(n_i,lnorm_grids[[2]],pdf_lnorm,range_lnorm)
  
  res_lnorm_pdf_gw_21_method3[i] = L1_Distance_method3(n_i,lnorm_grids[[2]],pdf_lnorm,range_lnorm)
  res_lnorm_pdf_gw_22_method3[i] = L2_Distance_method3(n_i,lnorm_grids[[2]],pdf_lnorm,range_lnorm)
  
  res_lnorm_pdf_gw_21_EM[i] = EM_L1_L2(n_i,lnorm_grids[[2]],pdf_lnorm,range_lnorm)$L1
  res_lnorm_pdf_gw_22_EM[i] = EM_L1_L2(n_i,lnorm_grids[[2]],pdf_lnorm,range_lnorm)$L2
  
  # Log-Normal Distribution - GW3
  
  n_i <- hist(x, breaks = lnorm_grids[[3]], plot = F)$counts
  
  res_lnorm_pdf_gw_31_method1[i] = L1_Distance_method1(n_i,lnorm_grids[[3]],pdf_lnorm,range_lnorm)
  res_lnorm_pdf_gw_32_method1[i] = L2_Distance_method1(n_i,lnorm_grids[[3]],pdf_lnorm,range_lnorm)
  
  res_lnorm_pdf_gw_31_method2[i] = L1_Distance_method2(n_i,lnorm_grids[[3]],pdf_lnorm,range_lnorm)
  res_lnorm_pdf_gw_32_method2[i] = L2_Distance_method2(n_i,lnorm_grids[[3]],pdf_lnorm,range_lnorm)
  
  res_lnorm_pdf_gw_31_method3[i] = L1_Distance_method3(n_i,lnorm_grids[[3]],pdf_lnorm,range_lnorm)
  res_lnorm_pdf_gw_32_method3[i] = L2_Distance_method3(n_i,lnorm_grids[[3]],pdf_lnorm,range_lnorm)
  
  res_lnorm_pdf_gw_31_EM[i] = EM_L1_L2(n_i,lnorm_grids[[3]],pdf_lnorm,range_lnorm)$L1
  res_lnorm_pdf_gw_32_EM[i] = EM_L1_L2(n_i,lnorm_grids[[3]],pdf_lnorm,range_lnorm)$L2
  
  # Log-Normal Distribution - GW4
  
  n_i <- hist(x, breaks = lnorm_grids[[4]], plot = F)$counts
  
  res_lnorm_pdf_gw_41_method1[i] = L1_Distance_method1(n_i,lnorm_grids[[4]],pdf_lnorm,range_lnorm)
  res_lnorm_pdf_gw_42_method1[i] = L2_Distance_method1(n_i,lnorm_grids[[4]],pdf_lnorm,range_lnorm)
  
  res_lnorm_pdf_gw_41_method2[i] = L1_Distance_method2(n_i,lnorm_grids[[4]],pdf_lnorm,range_lnorm)
  res_lnorm_pdf_gw_42_method2[i] = L2_Distance_method2(n_i,lnorm_grids[[4]],pdf_lnorm,range_lnorm)
  
  res_lnorm_pdf_gw_41_method3[i] = L1_Distance_method3(n_i,lnorm_grids[[4]],pdf_lnorm,range_lnorm)
  res_lnorm_pdf_gw_42_method3[i] = L2_Distance_method3(n_i,lnorm_grids[[4]],pdf_lnorm,range_lnorm)
  
  res_lnorm_pdf_gw_41_EM[i] = EM_L1_L2(n_i,lnorm_grids[[4]],pdf_lnorm,range_lnorm)$L1
  res_lnorm_pdf_gw_42_EM[i] = EM_L1_L2(n_i,lnorm_grids[[4]],pdf_lnorm,range_lnorm)$L2
  
  # Log-Normal Distribution - GW5
  
  n_i <- hist(x, breaks = lnorm_grids[[5]], plot = F)$counts
  
  res_lnorm_pdf_gw_51_method1[i] = L1_Distance_method1(n_i,lnorm_grids[[5]],pdf_lnorm,range_lnorm)
  res_lnorm_pdf_gw_52_method1[i] = L2_Distance_method1(n_i,lnorm_grids[[5]],pdf_lnorm,range_lnorm)
  
  res_lnorm_pdf_gw_51_method2[i] = L1_Distance_method2(n_i,lnorm_grids[[5]],pdf_lnorm,range_lnorm)
  res_lnorm_pdf_gw_52_method2[i] = L2_Distance_method2(n_i,lnorm_grids[[5]],pdf_lnorm,range_lnorm)
  
  res_lnorm_pdf_gw_51_method3[i] = L1_Distance_method3(n_i,lnorm_grids[[5]],pdf_lnorm,range_lnorm)
  res_lnorm_pdf_gw_52_method3[i] = L2_Distance_method3(n_i,lnorm_grids[[5]],pdf_lnorm,range_lnorm)
  
  res_lnorm_pdf_gw_51_EM[i] = EM_L1_L2(n_i,lnorm_grids[[5]],pdf_lnorm,range_lnorm)$L1
  res_lnorm_pdf_gw_52_EM[i] = EM_L1_L2(n_i,lnorm_grids[[5]],pdf_lnorm,range_lnorm)$L2
  
}

# Method-1 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L1-GW/L1-PDF-2-Lnorm-gw-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_lnorm_pdf_gw_11_method1,res_lnorm_pdf_gw_21_method1,res_lnorm_pdf_gw_31_method1,res_lnorm_pdf_gw_41_method1,res_lnorm_pdf_gw_51_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,1.5))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L2-GW/L2-PDF-2-Lnorm-gw-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_lnorm_pdf_gw_12_method1,res_lnorm_pdf_gw_22_method1,res_lnorm_pdf_gw_32_method1,res_lnorm_pdf_gw_42_method1,res_lnorm_pdf_gw_52_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,0.8))
dev.off()

# Method-2 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L1-GW/L1-PDF-2-Lnorm-gw-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_lnorm_pdf_gw_11_method2,res_lnorm_pdf_gw_21_method2,res_lnorm_pdf_gw_31_method2,res_lnorm_pdf_gw_41_method2,res_lnorm_pdf_gw_51_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,1.5))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L2-GW/L2-PDF-2-Lnorm-gw-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_lnorm_pdf_gw_12_method2,res_lnorm_pdf_gw_22_method2,res_lnorm_pdf_gw_32_method2,res_lnorm_pdf_gw_42_method2,res_lnorm_pdf_gw_52_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,0.8))
dev.off()

# Method-3 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L1-GW/L1-PDF-2-Lnorm-gw-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_lnorm_pdf_gw_11_method3,res_lnorm_pdf_gw_21_method3,res_lnorm_pdf_gw_31_method3,res_lnorm_pdf_gw_41_method3,res_lnorm_pdf_gw_51_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,1.5))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L2-GW/L2-PDF-2-Lnorm-gw-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_lnorm_pdf_gw_12_method3,res_lnorm_pdf_gw_22_method3,res_lnorm_pdf_gw_32_method3,res_lnorm_pdf_gw_42_method3,res_lnorm_pdf_gw_52_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,0.8))
dev.off()

# EM # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L1-GW/L1-PDF-2-Lnorm-gw-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_lnorm_pdf_gw_11_EM,res_lnorm_pdf_gw_21_EM,res_lnorm_pdf_gw_31_EM,res_lnorm_pdf_gw_41_EM,res_lnorm_pdf_gw_51_EM,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,1.5))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L2-GW/L2-PDF-2-Lnorm-gw-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_lnorm_pdf_gw_12_EM,res_lnorm_pdf_gw_22_EM,res_lnorm_pdf_gw_32_EM,res_lnorm_pdf_gw_42_EM,res_lnorm_pdf_gw_52_EM,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,0.8))
dev.off()

# Weibull Distribution - Grid Width Simulation 

MM = qweibull(0.9999999,shape) + 3
weibull_grids = list(seq(-MM,MM, by = weibull_deltas[5]),
                     seq(-MM,MM, by = weibull_deltas[4]),
                     seq(-MM,MM, by = weibull_deltas[3]),
                     seq(-MM,MM, by = weibull_deltas[2]),
                     seq(-MM,MM, by = weibull_deltas[1]))

# Method-1 # 

res_wei_pdf_gw_11_method1 = rep(0,B)
res_wei_pdf_gw_21_method1 = rep(0,B)
res_wei_pdf_gw_31_method1 = rep(0,B)
res_wei_pdf_gw_41_method1 = rep(0,B)
res_wei_pdf_gw_51_method1 = rep(0,B)

res_wei_pdf_gw_12_method1 = rep(0,B)
res_wei_pdf_gw_22_method1 = rep(0,B)
res_wei_pdf_gw_32_method1 = rep(0,B)
res_wei_pdf_gw_42_method1 = rep(0,B)
res_wei_pdf_gw_52_method1 = rep(0,B)

# Method-2 # 

res_wei_pdf_gw_11_method2 = rep(0,B)
res_wei_pdf_gw_21_method2 = rep(0,B)
res_wei_pdf_gw_31_method2 = rep(0,B)
res_wei_pdf_gw_41_method2 = rep(0,B)
res_wei_pdf_gw_51_method2 = rep(0,B)

res_wei_pdf_gw_12_method2 = rep(0,B)
res_wei_pdf_gw_22_method2 = rep(0,B)
res_wei_pdf_gw_32_method2 = rep(0,B)
res_wei_pdf_gw_42_method2 = rep(0,B)
res_wei_pdf_gw_52_method2 = rep(0,B)

# Method-3 # 

res_wei_pdf_gw_11_method3 = rep(0,B)
res_wei_pdf_gw_21_method3 = rep(0,B)
res_wei_pdf_gw_31_method3 = rep(0,B)
res_wei_pdf_gw_41_method3 = rep(0,B)
res_wei_pdf_gw_51_method3 = rep(0,B)

res_wei_pdf_gw_12_method3 = rep(0,B)
res_wei_pdf_gw_22_method3 = rep(0,B)
res_wei_pdf_gw_32_method3 = rep(0,B)
res_wei_pdf_gw_42_method3 = rep(0,B)
res_wei_pdf_gw_52_method3 = rep(0,B)

# EM # 

res_wei_pdf_gw_11_EM = rep(0,B)
res_wei_pdf_gw_21_EM = rep(0,B)
res_wei_pdf_gw_31_EM = rep(0,B)
res_wei_pdf_gw_41_EM = rep(0,B)
res_wei_pdf_gw_51_EM = rep(0,B)

res_wei_pdf_gw_12_EM = rep(0,B)
res_wei_pdf_gw_22_EM = rep(0,B)
res_wei_pdf_gw_32_EM = rep(0,B)
res_wei_pdf_gw_42_EM = rep(0,B)
res_wei_pdf_gw_52_EM = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  x = sort(rweibull(n,shape))
  
  # Weibull Distribution - GW1
  
  n_i <- hist(x, breaks = weibull_grids[[1]], plot = F)$counts
  
  res_wei_pdf_gw_11_method1[i] = L1_Distance_method1(n_i,weibull_grids[[1]],pdf_weibull,range_wei)
  res_wei_pdf_gw_12_method1[i] = L2_Distance_method1(n_i,weibull_grids[[1]],pdf_weibull,range_wei)
  
  res_wei_pdf_gw_11_method2[i] = L1_Distance_method2(n_i,weibull_grids[[1]],pdf_weibull,range_wei)
  res_wei_pdf_gw_12_method2[i] = L2_Distance_method2(n_i,weibull_grids[[1]],pdf_weibull,range_wei)
  
  res_wei_pdf_gw_11_method3[i] = L1_Distance_method3(n_i,weibull_grids[[1]],pdf_weibull,range_wei)
  res_wei_pdf_gw_12_method3[i] = L2_Distance_method3(n_i,weibull_grids[[1]],pdf_weibull,range_wei)
  
  res_wei_pdf_gw_11_EM[i] = EM_L1_L2(n_i,weibull_grids[[1]],pdf_weibull,range_wei)$L1
  res_wei_pdf_gw_12_EM[i] = EM_L1_L2(n_i,weibull_grids[[1]],pdf_weibull,range_wei)$L2
  
  # Weibull Distribution - GW2
  
  n_i <- hist(x, breaks = weibull_grids[[2]], plot = F)$counts
  
  res_wei_pdf_gw_21_method1[i] = L1_Distance_method1(n_i,weibull_grids[[2]],pdf_weibull,range_wei)
  res_wei_pdf_gw_22_method1[i] = L2_Distance_method1(n_i,weibull_grids[[2]],pdf_weibull,range_wei)
  
  res_wei_pdf_gw_21_method2[i] = L1_Distance_method2(n_i,weibull_grids[[2]],pdf_weibull,range_wei)
  res_wei_pdf_gw_22_method2[i] = L2_Distance_method2(n_i,weibull_grids[[2]],pdf_weibull,range_wei)
  
  res_wei_pdf_gw_21_method3[i] = L1_Distance_method3(n_i,weibull_grids[[2]],pdf_weibull,range_wei)
  res_wei_pdf_gw_22_method3[i] = L2_Distance_method3(n_i,weibull_grids[[2]],pdf_weibull,range_wei)
  
  res_wei_pdf_gw_21_EM[i] = EM_L1_L2(n_i,weibull_grids[[2]],pdf_weibull,range_wei)$L1
  res_wei_pdf_gw_22_EM[i] = EM_L1_L2(n_i,weibull_grids[[2]],pdf_weibull,range_wei)$L2
  
  # Weibull Distribution - GW3
  
  n_i <- hist(x, breaks = weibull_grids[[3]], plot = F)$counts
  
  res_wei_pdf_gw_31_method1[i] = L1_Distance_method1(n_i,weibull_grids[[3]],pdf_weibull,range_wei)
  res_wei_pdf_gw_32_method1[i] = L2_Distance_method1(n_i,weibull_grids[[3]],pdf_weibull,range_wei)
  
  res_wei_pdf_gw_31_method2[i] = L1_Distance_method2(n_i,weibull_grids[[3]],pdf_weibull,range_wei)
  res_wei_pdf_gw_32_method2[i] = L2_Distance_method2(n_i,weibull_grids[[3]],pdf_weibull,range_wei)
  
  res_wei_pdf_gw_31_method3[i] = L1_Distance_method3(n_i,weibull_grids[[3]],pdf_weibull,range_wei)
  res_wei_pdf_gw_32_method3[i] = L2_Distance_method3(n_i,weibull_grids[[3]],pdf_weibull,range_wei)
  
  res_wei_pdf_gw_31_EM[i] = EM_L1_L2(n_i,weibull_grids[[3]],pdf_weibull,range_wei)$L1
  res_wei_pdf_gw_32_EM[i] = EM_L1_L2(n_i,weibull_grids[[3]],pdf_weibull,range_wei)$L2
  
  # Weibull Distribution - GW4
  
  n_i <- hist(x, breaks = weibull_grids[[4]], plot = F)$counts
  
  res_wei_pdf_gw_41_method1[i] = L1_Distance_method1(n_i,weibull_grids[[4]],pdf_weibull,range_wei)
  res_wei_pdf_gw_42_method1[i] = L2_Distance_method1(n_i,weibull_grids[[4]],pdf_weibull,range_wei)
  
  res_wei_pdf_gw_41_method2[i] = L1_Distance_method2(n_i,weibull_grids[[4]],pdf_weibull,range_wei)
  res_wei_pdf_gw_42_method2[i] = L2_Distance_method2(n_i,weibull_grids[[4]],pdf_weibull,range_wei)
  
  res_wei_pdf_gw_41_method3[i] = L1_Distance_method3(n_i,weibull_grids[[4]],pdf_weibull,range_wei)
  res_wei_pdf_gw_42_method3[i] = L2_Distance_method3(n_i,weibull_grids[[4]],pdf_weibull,range_wei)
  
  res_wei_pdf_gw_41_EM[i] = EM_L1_L2(n_i,weibull_grids[[4]],pdf_weibull,range_wei)$L1
  res_wei_pdf_gw_42_EM[i] = EM_L1_L2(n_i,weibull_grids[[4]],pdf_weibull,range_wei)$L2
  
  # Weibull Distribution - GW5
  
  n_i <- hist(x, breaks = weibull_grids[[5]], plot = F)$counts
  
  res_wei_pdf_gw_51_method1[i] = L1_Distance_method1(n_i,weibull_grids[[5]],pdf_weibull,range_wei)
  res_wei_pdf_gw_52_method1[i] = L2_Distance_method1(n_i,weibull_grids[[5]],pdf_weibull,range_wei)
  
  res_wei_pdf_gw_51_method2[i] = L1_Distance_method2(n_i,weibull_grids[[5]],pdf_weibull,range_wei)
  res_wei_pdf_gw_52_method2[i] = L2_Distance_method2(n_i,weibull_grids[[5]],pdf_weibull,range_wei)
  
  res_wei_pdf_gw_51_method3[i] = L1_Distance_method3(n_i,weibull_grids[[5]],pdf_weibull,range_wei)
  res_wei_pdf_gw_52_method3[i] = L2_Distance_method3(n_i,weibull_grids[[5]],pdf_weibull,range_wei)
  
  res_wei_pdf_gw_51_EM[i] = EM_L1_L2(n_i,weibull_grids[[5]],pdf_weibull,range_wei)$L1
  res_wei_pdf_gw_52_EM[i] = EM_L1_L2(n_i,weibull_grids[[5]],pdf_weibull,range_wei)$L2
  
}

# Method-1 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L1-GW/L1-PDF-2-Wei-gw-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_wei_pdf_gw_11_method1,res_wei_pdf_gw_21_method1,res_wei_pdf_gw_31_method1,res_wei_pdf_gw_41_method1,res_wei_pdf_gw_51_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,0.9))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L2-GW/L2-PDF-2-Wei-gw-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_wei_pdf_gw_12_method1,res_wei_pdf_gw_22_method1,res_wei_pdf_gw_32_method1,res_wei_pdf_gw_42_method1,res_wei_pdf_gw_52_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,0.7))
dev.off()

# Method-2 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L1-GW/L1-PDF-2-Wei-gw-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_wei_pdf_gw_11_method2,res_wei_pdf_gw_21_method2,res_wei_pdf_gw_31_method2,res_wei_pdf_gw_41_method2,res_wei_pdf_gw_51_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,0.9))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L2-GW/L2-PDF-2-Wei-gw-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_wei_pdf_gw_12_method2,res_wei_pdf_gw_22_method2,res_wei_pdf_gw_32_method2,res_wei_pdf_gw_42_method2,res_wei_pdf_gw_52_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,0.4))
dev.off()

# Method-3 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L1-GW/L1-PDF-2-Wei-gw-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_wei_pdf_gw_11_method3,res_wei_pdf_gw_21_method3,res_wei_pdf_gw_31_method3,res_wei_pdf_gw_41_method3,res_wei_pdf_gw_51_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,0.9))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L2-GW/L2-PDF-2-Wei-gw-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_wei_pdf_gw_12_method3,res_wei_pdf_gw_22_method3,res_wei_pdf_gw_32_method3,res_wei_pdf_gw_42_method3,res_wei_pdf_gw_52_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,0.7))
dev.off()

# EM # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L1-GW/L1-PDF-2-Wei-gw-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_wei_pdf_gw_11_EM,res_wei_pdf_gw_21_EM,res_wei_pdf_gw_31_EM,res_wei_pdf_gw_41_EM,res_wei_pdf_gw_51_EM,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,0.9))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L2-GW/L2-PDF-2-Wei-gw-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_wei_pdf_gw_12_EM,res_wei_pdf_gw_22_EM,res_wei_pdf_gw_32_EM,res_wei_pdf_gw_42_EM,res_wei_pdf_gw_52_EM,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,0.7))
dev.off()

# Pareto Distribution - Grid Width Simulation 

MM = qpareto(0.9999999,location,shape_p) + 3
pareto_grids = list(seq(location,MM, by = pareto_deltas[5]),
                    seq(location,MM, by = pareto_deltas[4]),
                    seq(location,MM, by = pareto_deltas[3]),
                    seq(location,MM, by = pareto_deltas[2]),
                    seq(location,MM, by = pareto_deltas[1]))

# Method-1 # 

res_pareto_pdf_gw_11_method1 = rep(0,B)
res_pareto_pdf_gw_21_method1 = rep(0,B)
res_pareto_pdf_gw_31_method1 = rep(0,B)
res_pareto_pdf_gw_41_method1 = rep(0,B)
res_pareto_pdf_gw_51_method1 = rep(0,B)

res_pareto_pdf_gw_12_method1 = rep(0,B)
res_pareto_pdf_gw_22_method1 = rep(0,B)
res_pareto_pdf_gw_32_method1 = rep(0,B)
res_pareto_pdf_gw_42_method1 = rep(0,B)
res_pareto_pdf_gw_52_method1 = rep(0,B)

# Method-2 # 

res_pareto_pdf_gw_11_method2 = rep(0,B)
res_pareto_pdf_gw_21_method2 = rep(0,B)
res_pareto_pdf_gw_31_method2 = rep(0,B)
res_pareto_pdf_gw_41_method2 = rep(0,B)
res_pareto_pdf_gw_51_method2 = rep(0,B)

res_pareto_pdf_gw_12_method2 = rep(0,B)
res_pareto_pdf_gw_22_method2 = rep(0,B)
res_pareto_pdf_gw_32_method2 = rep(0,B)
res_pareto_pdf_gw_42_method2 = rep(0,B)
res_pareto_pdf_gw_52_method2 = rep(0,B)

# Method-3 # 

res_pareto_pdf_gw_11_method3 = rep(0,B)
res_pareto_pdf_gw_21_method3 = rep(0,B)
res_pareto_pdf_gw_31_method3 = rep(0,B)
res_pareto_pdf_gw_41_method3 = rep(0,B)
res_pareto_pdf_gw_51_method3 = rep(0,B)

res_pareto_pdf_gw_12_method3 = rep(0,B)
res_pareto_pdf_gw_22_method3 = rep(0,B)
res_pareto_pdf_gw_32_method3 = rep(0,B)
res_pareto_pdf_gw_42_method3 = rep(0,B)
res_pareto_pdf_gw_52_method3 = rep(0,B)

# EM # 

res_pareto_pdf_gw_11_EM = rep(0,B)
res_pareto_pdf_gw_21_EM = rep(0,B)
res_pareto_pdf_gw_31_EM = rep(0,B)
res_pareto_pdf_gw_41_EM = rep(0,B)
res_pareto_pdf_gw_51_EM = rep(0,B)

res_pareto_pdf_gw_12_EM = rep(0,B)
res_pareto_pdf_gw_22_EM = rep(0,B)
res_pareto_pdf_gw_32_EM = rep(0,B)
res_pareto_pdf_gw_42_EM = rep(0,B)
res_pareto_pdf_gw_52_EM = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  x = sort(rpareto(n,location,shape_p))
  
  # Pareto Distribution - GW1
  
  n_i <- hist(x, breaks = pareto_grids[[1]], plot = F)$counts
  
  res_pareto_pdf_gw_11_method1[i] = L1_Distance_method1(n_i,pareto_grids[[1]],pdf_pareto,range_pareto)
  res_pareto_pdf_gw_12_method1[i] = L2_Distance_method1(n_i,pareto_grids[[1]],pdf_pareto,range_pareto)
  
  res_pareto_pdf_gw_11_method2[i] = L1_Distance_method2(n_i,pareto_grids[[1]],pdf_pareto,range_pareto)
  res_pareto_pdf_gw_12_method2[i] = L2_Distance_method2(n_i,pareto_grids[[1]],pdf_pareto,range_pareto)
  
  res_pareto_pdf_gw_11_method3[i] = L1_Distance_method3(n_i,pareto_grids[[1]],pdf_pareto,range_pareto)
  res_pareto_pdf_gw_12_method3[i] = L2_Distance_method3(n_i,pareto_grids[[1]],pdf_pareto,range_pareto)
  
  res_pareto_pdf_gw_11_EM[i] = EM_L1_L2(n_i,pareto_grids[[1]],pdf_pareto,range_pareto)$L1
  res_pareto_pdf_gw_12_EM[i] = EM_L1_L2(n_i,pareto_grids[[1]],pdf_pareto,range_pareto)$L2
  
  # Pareto Distribution - GW2
  
  n_i <- hist(x, breaks = pareto_grids[[2]], plot = F)$counts
  
  res_pareto_pdf_gw_21_method1[i] = L1_Distance_method1(n_i,pareto_grids[[2]],pdf_pareto,range_pareto)
  res_pareto_pdf_gw_22_method1[i] = L2_Distance_method1(n_i,pareto_grids[[2]],pdf_pareto,range_pareto)

  res_pareto_pdf_gw_21_method2[i] = L1_Distance_method2(n_i,pareto_grids[[2]],pdf_pareto,range_pareto)
  res_pareto_pdf_gw_22_method2[i] = L2_Distance_method2(n_i,pareto_grids[[2]],pdf_pareto,range_pareto)
  
  res_pareto_pdf_gw_21_method3[i] = L1_Distance_method3(n_i,pareto_grids[[2]],pdf_pareto,range_pareto)
  res_pareto_pdf_gw_22_method3[i] = L2_Distance_method3(n_i,pareto_grids[[2]],pdf_pareto,range_pareto)
  
  res_pareto_pdf_gw_21_EM[i] = EM_L1_L2(n_i,pareto_grids[[2]],pdf_pareto,range_pareto)$L1
  res_pareto_pdf_gw_22_EM[i] = EM_L1_L2(n_i,pareto_grids[[2]],pdf_pareto,range_pareto)$L2
  
  # Pareto Distribution - GW3
  
  n_i <- hist(x, breaks = pareto_grids[[3]], plot = F)$counts
  
  res_pareto_pdf_gw_31_method1[i] = L1_Distance_method1(n_i,pareto_grids[[3]],pdf_pareto,range_pareto)
  res_pareto_pdf_gw_32_method1[i] = L2_Distance_method1(n_i,pareto_grids[[3]],pdf_pareto,range_pareto)
  
  res_pareto_pdf_gw_31_method2[i] = L1_Distance_method2(n_i,pareto_grids[[3]],pdf_pareto,range_pareto)
  res_pareto_pdf_gw_32_method2[i] = L2_Distance_method2(n_i,pareto_grids[[3]],pdf_pareto,range_pareto)
  
  res_pareto_pdf_gw_31_method3[i] = L1_Distance_method3(n_i,pareto_grids[[3]],pdf_pareto,range_pareto)
  res_pareto_pdf_gw_32_method3[i] = L2_Distance_method3(n_i,pareto_grids[[3]],pdf_pareto,range_pareto)
  
  res_pareto_pdf_gw_31_EM[i] = EM_L1_L2(n_i,pareto_grids[[3]],pdf_pareto,range_pareto)$L1
  res_pareto_pdf_gw_32_EM[i] = EM_L1_L2(n_i,pareto_grids[[3]],pdf_pareto,range_pareto)$L2
  
  # Pareto Distribution - GW4
  
  n_i <- hist(x, breaks = pareto_grids[[4]], plot = F)$counts
  
  res_pareto_pdf_gw_41_method1[i] = L1_Distance_method1(n_i,pareto_grids[[4]],pdf_pareto,range_pareto)
  res_pareto_pdf_gw_42_method1[i] = L2_Distance_method1(n_i,pareto_grids[[4]],pdf_pareto,range_pareto)
  
  res_pareto_pdf_gw_41_method2[i] = L1_Distance_method2(n_i,pareto_grids[[4]],pdf_pareto,range_pareto)
  res_pareto_pdf_gw_42_method2[i] = L2_Distance_method2(n_i,pareto_grids[[4]],pdf_pareto,range_pareto)
  
  res_pareto_pdf_gw_41_method3[i] = L1_Distance_method3(n_i,pareto_grids[[4]],pdf_pareto,range_pareto)
  res_pareto_pdf_gw_42_method3[i] = L2_Distance_method3(n_i,pareto_grids[[4]],pdf_pareto,range_pareto)
  
  res_pareto_pdf_gw_41_EM[i] = EM_L1_L2(n_i,pareto_grids[[4]],pdf_pareto,range_pareto)$L1
  res_pareto_pdf_gw_42_EM[i] = EM_L1_L2(n_i,pareto_grids[[4]],pdf_pareto,range_pareto)$L2
  
  # Pareto Distribution - GW5
  
  n_i <- hist(x, breaks = pareto_grids[[5]], plot = F)$counts
  
  res_pareto_pdf_gw_51_method1[i] = L1_Distance_method1(n_i,pareto_grids[[5]],pdf_pareto,range_pareto)
  res_pareto_pdf_gw_52_method1[i] = L2_Distance_method1(n_i,pareto_grids[[5]],pdf_pareto,range_pareto)

  res_pareto_pdf_gw_51_method2[i] = L1_Distance_method2(n_i,pareto_grids[[5]],pdf_pareto,range_pareto)
  res_pareto_pdf_gw_52_method2[i] = L2_Distance_method2(n_i,pareto_grids[[5]],pdf_pareto,range_pareto)

  res_pareto_pdf_gw_51_method3[i] = L1_Distance_method3(n_i,pareto_grids[[5]],pdf_pareto,range_pareto)
  res_pareto_pdf_gw_52_method3[i] = L2_Distance_method3(n_i,pareto_grids[[5]],pdf_pareto,range_pareto)
  
  res_pareto_pdf_gw_51_EM[i] = EM_L1_L2(n_i,pareto_grids[[5]],pdf_pareto,range_pareto)$L1
  res_pareto_pdf_gw_52_EM[i] = EM_L1_L2(n_i,pareto_grids[[5]],pdf_pareto,range_pareto)$L2
  
}

# Method-1 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L1-GW/L1-PDF-2-Pareto-gw-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_pareto_pdf_gw_11_method1,res_pareto_pdf_gw_21_method1,res_pareto_pdf_gw_31_method1,res_pareto_pdf_gw_41_method1,res_pareto_pdf_gw_51_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,1.7))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L2-GW/L2-PDF-2-Pareto-gw-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_pareto_pdf_gw_12_method1,res_pareto_pdf_gw_22_method1,res_pareto_pdf_gw_32_method1,res_pareto_pdf_gw_42_method1,res_pareto_pdf_gw_52_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,4.2))
dev.off()

# Method-2 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L1-GW/L1-PDF-2-Pareto-gw-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_pareto_pdf_gw_11_method2,res_pareto_pdf_gw_21_method2,res_pareto_pdf_gw_31_method2,res_pareto_pdf_gw_41_method2,res_pareto_pdf_gw_51_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,1.7))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L2-GW/L2-PDF-2-Pareto-gw-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_pareto_pdf_gw_12_method2,res_pareto_pdf_gw_22_method2,res_pareto_pdf_gw_32_method2,res_pareto_pdf_gw_42_method2,res_pareto_pdf_gw_52_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,4.2))
dev.off()

# Method-3 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L1-GW/L1-PDF-2-Pareto-gw-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_pareto_pdf_gw_11_method3,res_pareto_pdf_gw_21_method3,res_pareto_pdf_gw_31_method3,res_pareto_pdf_gw_41_method3,res_pareto_pdf_gw_51_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,1.7))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L2-GW/L2-PDF-2-Pareto-gw-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_pareto_pdf_gw_12_method3,res_pareto_pdf_gw_22_method3,res_pareto_pdf_gw_32_method3,res_pareto_pdf_gw_42_method3,res_pareto_pdf_gw_52_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,4.2))
dev.off()

# EM # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L1-GW/L1-PDF-2-Pareto-gw-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_pareto_pdf_gw_11_EM,res_pareto_pdf_gw_21_EM,res_pareto_pdf_gw_31_EM,res_pareto_pdf_gw_41_EM,res_pareto_pdf_gw_51_EM,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,1.7))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L2-GW/L2-PDF-2-Pareto-gw-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_pareto_pdf_gw_12_EM,res_pareto_pdf_gw_22_EM,res_pareto_pdf_gw_32_EM,res_pareto_pdf_gw_42_EM,res_pareto_pdf_gw_52_EM,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,4.2))
dev.off()

