# Simulation Study 

n <- 1000000
B = 500

res_pdf_L1-norm = rep(0,B)
res_pdf_L1_beta = rep(0,B)
res_pdf_L1_gamma = rep(0,B)
res_pdf_L1_logis = rep(0,B)
res_pdf_L1_t = rep(0,B)

res_pdf_L2-norm = rep(0,B)
res_pdf_L2_beta = rep(0,B)
res_pdf_L2_gamma = rep(0,B)
res_pdf_L2_logis = rep(0,B)
res_pdf_L2_t = rep(0,B)

# Parameters

alpha=2;beta=5;
scale_val = 0.05;df=3;
shape=2;location=2;
mu = 0;sigma=1;
mul = 1;
df=3;
shape_p = 4;

# Ranges

range_norm = c(-Inf,Inf)
range_beta = c(0,1)
range_gamma = c(0,Inf)
range_logistic = c(-Inf,Inf)
range_t = c(-Inf,Inf)

# Pdfs 

pdf_beta <- function(x) {
  dbeta(x,alpha,beta)
}

pdf_norm <- function(x) {
  dnorm(x, mean = mu, sd = sigma)
}

pdf_gamma <- function(x) {
  dgamma(x, alpha, beta)
}

pdf_logistic <- function(x) {
  dlogis(x)
}

pdf_t <- function(x) {
  dt(x,df)
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

# Convergence Limits

pdf_L1-norm_limit = mean(res_pdf_L1-norm)
pdf_L2-norm_limit = mean(res_pdf_L2-norm)
pdf_L1_gamma_limit = mean(res_pdf_L1_gamma)
pdf_L2_gamma_limit = mean(res_pdf_L2_gamma)
pdf_L1_beta_limit = mean(res_pdf_L1_beta)
pdf_L2_beta_limit = mean(res_pdf_L2_beta)
pdf_L1_logis_limit = mean(res_pdf_L1_logis)
pdf_L2_logis_limit = mean(res_pdf_L2_logis)
pdf_L1_t_limit = mean(res_pdf_L1_t)
pdf_L2_t_limit = mean(res_pdf_L2_t)

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

# Normal Distribution - Sample Size Simulation 

n = c(10^1,10^2,10^3,10^4,10^5)
B = 100

MM = qnorm(0.9999999,mu,sigma) + 3
norm_grid <- seq(-MM,MM, by = norm_deltas[1])

# Method-1 Residuals # 

res_norm_pdf_n_11_method1 = rep(0,B)
res_norm_pdf_n_21_method1 = rep(0,B)
res_norm_pdf_n_31_method1 = rep(0,B)
res_norm_pdf_n_41_method1 = rep(0,B)
res_norm_pdf_n_51_method1 = rep(0,B)

res_norm_pdf_n_12_method1 = rep(0,B)
res_norm_pdf_n_22_method1 = rep(0,B)
res_norm_pdf_n_32_method1 = rep(0,B)
res_norm_pdf_n_42_method1 = rep(0,B)
res_norm_pdf_n_52_method1 = rep(0,B)

# Method-2 Residuals # 

res_norm_pdf_n_11_method2 = rep(0,B)
res_norm_pdf_n_21_method2 = rep(0,B)
res_norm_pdf_n_31_method2 = rep(0,B)
res_norm_pdf_n_41_method2 = rep(0,B)
res_norm_pdf_n_51_method2 = rep(0,B)

res_norm_pdf_n_12_method2 = rep(0,B)
res_norm_pdf_n_22_method2 = rep(0,B)
res_norm_pdf_n_32_method2 = rep(0,B)
res_norm_pdf_n_42_method2 = rep(0,B)
res_norm_pdf_n_52_method2 = rep(0,B)

# Method-3 Residuals # 

res_norm_pdf_n_11_method3 = rep(0,B)
res_norm_pdf_n_21_method3 = rep(0,B)
res_norm_pdf_n_31_method3 = rep(0,B)
res_norm_pdf_n_41_method3 = rep(0,B)
res_norm_pdf_n_51_method3 = rep(0,B)

res_norm_pdf_n_12_method3 = rep(0,B)
res_norm_pdf_n_22_method3 = rep(0,B)
res_norm_pdf_n_32_method3 = rep(0,B)
res_norm_pdf_n_42_method3 = rep(0,B)
res_norm_pdf_n_52_method3 = rep(0,B)

# EM # 

res_norm_pdf_n_11_EM = rep(0,B)
res_norm_pdf_n_21_EM = rep(0,B)
res_norm_pdf_n_31_EM = rep(0,B)
res_norm_pdf_n_41_EM = rep(0,B)
res_norm_pdf_n_51_EM = rep(0,B)

res_norm_pdf_n_12_EM = rep(0,B)
res_norm_pdf_n_22_EM = rep(0,B)
res_norm_pdf_n_32_EM = rep(0,B)
res_norm_pdf_n_42_EM = rep(0,B)
res_norm_pdf_n_52_EM = rep(0,B)

set.seed(10)

for (i in 51:B) {
  
  # Normal Distribution - n1
  
  x <- sort(rnorm(n[1], mean = mu,sd=sigma))
  n_i <- hist(x, breaks = norm_grid, plot = F)$counts
  
  res_norm_pdf_n_11_method1[i] = L1_Distance_method1(n_i,norm_grid,pdf_norm,range_norm)
  res_norm_pdf_n_12_method1[i] = L2_Distance_method1(n_i,norm_grid,pdf_norm,range_norm)
  
  res_norm_pdf_n_11_method2[i] = L1_Distance_method2(n_i,norm_grid,pdf_norm,range_norm)
  res_norm_pdf_n_12_method2[i] = L2_Distance_method2(n_i,norm_grid,pdf_norm,range_norm)
  
  res_norm_pdf_n_11_method3[i] = L1_Distance_method3(n_i,norm_grid,pdf_norm,range_norm)
  res_norm_pdf_n_12_method3[i] = L2_Distance_method3(n_i,norm_grid,pdf_norm,range_norm)
  
  res_norm_pdf_n_11_EM[i] = EM_L1_L2(n_i,norm_grid,pdf_norm,range_norm)$L1
  res_norm_pdf_n_12_EM[i] = EM_L1_L2(n_i,norm_grid,pdf_norm,range_norm)$L2
  
  # Normal Distribution - n2
  
  x <- sort(rnorm(n[2], mean = mu,sd=sigma))
  n_i <- hist(x, breaks = norm_grid, plot = F)$counts
  
  res_norm_pdf_n_21_method1[i] = L1_Distance_method1(n_i,norm_grid,pdf_norm,range_norm)
  res_norm_pdf_n_22_method1[i] = L2_Distance_method1(n_i,norm_grid,pdf_norm,range_norm)
  
  res_norm_pdf_n_21_method2[i] = L1_Distance_method2(n_i,norm_grid,pdf_norm,range_norm)
  res_norm_pdf_n_22_method2[i] = L2_Distance_method2(n_i,norm_grid,pdf_norm,range_norm)
  
  res_norm_pdf_n_21_method3[i] = L1_Distance_method3(n_i,norm_grid,pdf_norm,range_norm)
  res_norm_pdf_n_22_method3[i] = L2_Distance_method3(n_i,norm_grid,pdf_norm,range_norm)
  
  res_norm_pdf_n_21_EM[i] = EM_L1_L2(n_i,norm_grid,pdf_norm,range_norm)$L1
  res_norm_pdf_n_22_EM[i] = EM_L1_L2(n_i,norm_grid,pdf_norm,range_norm)$L2
  
  # Normal Distribution - n3
  
  x <- sort(rnorm(n[3], mean = mu,sd=sigma))
  n_i <- hist(x, breaks = norm_grid, plot = F)$counts
  
  res_norm_pdf_n_31_method1[i] = L1_Distance_method1(n_i,norm_grid,pdf_norm,range_norm)
  res_norm_pdf_n_32_method1[i] = L2_Distance_method1(n_i,norm_grid,pdf_norm,range_norm)
  
  res_norm_pdf_n_31_method2[i] = L1_Distance_method2(n_i,norm_grid,pdf_norm,range_norm)
  res_norm_pdf_n_32_method2[i] = L2_Distance_method2(n_i,norm_grid,pdf_norm,range_norm)
  
  res_norm_pdf_n_31_method3[i] = L1_Distance_method3(n_i,norm_grid,pdf_norm,range_norm)
  res_norm_pdf_n_32_method3[i] = L2_Distance_method3(n_i,norm_grid,pdf_norm,range_norm)
  
  res_norm_pdf_n_31_EM[i] = EM_L1_L2(n_i,norm_grid,pdf_norm,range_norm)$L1
  res_norm_pdf_n_32_EM[i] = EM_L1_L2(n_i,norm_grid,pdf_norm,range_norm)$L2
  
  # Normal Distribution - n4
  
  x <- sort(rnorm(n[4], mean = mu,sd=sigma))
  n_i <- hist(x, breaks = norm_grid, plot = F)$counts
  
  res_norm_pdf_n_41_method1[i] = L1_Distance_method1(n_i,norm_grid,pdf_norm,range_norm)
  res_norm_pdf_n_42_method1[i] = L2_Distance_method1(n_i,norm_grid,pdf_norm,range_norm)
  
  res_norm_pdf_n_41_method2[i] = L1_Distance_method2(n_i,norm_grid,pdf_norm,range_norm)
  res_norm_pdf_n_42_method2[i] = L2_Distance_method2(n_i,norm_grid,pdf_norm,range_norm)
  
  res_norm_pdf_n_41_method3[i] = L1_Distance_method3(n_i,norm_grid,pdf_norm,range_norm)
  res_norm_pdf_n_42_method3[i] = L2_Distance_method3(n_i,norm_grid,pdf_norm,range_norm)
  
  res_norm_pdf_n_41_EM[i] = EM_L1_L2(n_i,norm_grid,pdf_norm,range_norm)$L1
  res_norm_pdf_n_42_EM[i] = EM_L1_L2(n_i,norm_grid,pdf_norm,range_norm)$L2
  
  # Normal Distribution - n5
  
  x <- sort(rnorm(n[5], mean = mu,sd=sigma))
  n_i <- hist(x, breaks = norm_grid, plot = F)$counts
  
  res_norm_pdf_n_51_method1[i] = L1_Distance_method1(n_i,norm_grid,pdf_norm,range_norm)
  res_norm_pdf_n_52_method1[i] = L2_Distance_method1(n_i,norm_grid,pdf_norm,range_norm)
  
  res_norm_pdf_n_51_method2[i] = L1_Distance_method2(n_i,norm_grid,pdf_norm,range_norm)
  res_norm_pdf_n_52_method2[i] = L2_Distance_method2(n_i,norm_grid,pdf_norm,range_norm)
  
  res_norm_pdf_n_51_method3[i] = L1_Distance_method3(n_i,norm_grid,pdf_norm,range_norm)
  res_norm_pdf_n_52_method3[i] = L2_Distance_method3(n_i,norm_grid,pdf_norm,range_norm)
  
  res_norm_pdf_n_51_EM[i] = EM_L1_L2(n_i,norm_grid,pdf_norm,range_norm)$L1
  res_norm_pdf_n_52_EM[i] = EM_L1_L2(n_i,norm_grid,pdf_norm,range_norm)$L2
  
  print(i)
  
}

# Method-1 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L1-n/L1-PDF-1-Normal-n-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_norm_pdf_n_11_method1,res_norm_pdf_n_21_method1,res_norm_pdf_n_31_method1,res_norm_pdf_n_41_method1,res_norm_pdf_n_51_method1,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors1,ylim=c(0,0.8))
#segments(x0 = 0.29, y0 = pdf_L1-norm_limit, x1 = 10^5, y1 = pdf_L1-norm_limit, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L2-n/L2-PDF-1-Normal-n-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_norm_pdf_n_12_method1,res_norm_pdf_n_22_method1,res_norm_pdf_n_32_method1,res_norm_pdf_n_42_method1,res_norm_pdf_n_52_method1,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors1,ylim=c(0,0.6))
# segments(x0 = 0.29, y0 = pdf_L2-norm_limit, x1 = 10^5, y1 = pdf_L2-norm_limit, col = "red", lty = 2)
dev.off()

# Method-2 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L1-n/L1-PDF-1-Normal-n-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_norm_pdf_n_11_method2,res_norm_pdf_n_21_method2,res_norm_pdf_n_31_method2,res_norm_pdf_n_41_method2,res_norm_pdf_n_51_method2,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors1,ylim=c(0,0.8))
# segments(x0 = 0.29, y0 = pdf_L1-norm_limit, x1 = 10^5, y1 = pdf_L1-norm_limit, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L2-n/L2-PDF-1-Normal-n-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_norm_pdf_n_12_method2,res_norm_pdf_n_22_method2,res_norm_pdf_n_32_method2,res_norm_pdf_n_42_method2,res_norm_pdf_n_52_method2,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors1,ylim=c(0,0.6))
# segments(x0 = 0.29, y0 = pdf_L2-norm_limit, x1 = 10^5, y1 = pdf_L2-norm_limit, col = "red", lty = 2)
dev.off()

# Method-3 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L1-n/L1-PDF-1-Normal-n-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_norm_pdf_n_11_method3,res_norm_pdf_n_21_method3,res_norm_pdf_n_31_method3,res_norm_pdf_n_41_method3,res_norm_pdf_n_51_method3,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors1,ylim=c(0,0.8))
# segments(x0 = 0.29, y0 = pdf_L1-norm_limit, x1 = 10^5, y1 = pdf_L1-norm_limit, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L2-n/L2-PDF-1-Normal-n-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_norm_pdf_n_12_method3,res_norm_pdf_n_22_method3,res_norm_pdf_n_32_method3,res_norm_pdf_n_42_method3,res_norm_pdf_n_52_method3,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors1,ylim=c(0,0.6))
# segments(x0 = 0.29, y0 = pdf_L2-norm_limit, x1 = 10^5, y1 = pdf_L2-norm_limit, col = "red", lty = 2)
dev.off()

# EM # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L1-n/L1-PDF-1-Normal-n-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_norm_pdf_n_11_EM,res_norm_pdf_n_21_EM,res_norm_pdf_n_31_EM,res_norm_pdf_n_41_EM,res_norm_pdf_n_51_EM,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors1,ylim=c(0,0.8))
# segments(x0 = 0.29, y0 = pdf_L1-norm_limit, x1 = 10^5, y1 = pdf_L1-norm_limit, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L2-n/L2-PDF-1-Normal-n-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_norm_pdf_n_12_EM,res_norm_pdf_n_22_EM,res_norm_pdf_n_32_EM,res_norm_pdf_n_42_EM,res_norm_pdf_n_52_EM,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors1,ylim=c(0,0.6))
# segments(x0 = 0.29, y0 = pdf_L2-norm_limit, x1 = 10^5, y1 = pdf_L2-norm_limit, col = "red", lty = 2)
dev.off()

# Beta Distribution - Sample Size Simulation 

beta_grid <- seq(0,1.04, by = beta_deltas[1])

# Method 1 #

res_beta_pdf_n_11_method1 = rep(0,B)
res_beta_pdf_n_21_method1 = rep(0,B)
res_beta_pdf_n_31_method1 = rep(0,B)
res_beta_pdf_n_41_method1 = rep(0,B)
res_beta_pdf_n_51_method1 = rep(0,B)

res_beta_pdf_n_12_method1 = rep(0,B)
res_beta_pdf_n_22_method1 = rep(0,B)
res_beta_pdf_n_32_method1 = rep(0,B)
res_beta_pdf_n_42_method1 = rep(0,B)
res_beta_pdf_n_52_method1 = rep(0,B)

# Method 2 #

res_beta_pdf_n_11_method2 = rep(0,B)
res_beta_pdf_n_21_method2 = rep(0,B)
res_beta_pdf_n_31_method2 = rep(0,B)
res_beta_pdf_n_41_method2 = rep(0,B)
res_beta_pdf_n_51_method2 = rep(0,B)

res_beta_pdf_n_12_method2 = rep(0,B)
res_beta_pdf_n_22_method2 = rep(0,B)
res_beta_pdf_n_32_method2 = rep(0,B)
res_beta_pdf_n_42_method2 = rep(0,B)
res_beta_pdf_n_52_method2 = rep(0,B)

# Method 3 #

res_beta_pdf_n_11_method3 = rep(0,B)
res_beta_pdf_n_21_method3 = rep(0,B)
res_beta_pdf_n_31_method3 = rep(0,B)
res_beta_pdf_n_41_method3 = rep(0,B)
res_beta_pdf_n_51_method3 = rep(0,B)

res_beta_pdf_n_12_method3 = rep(0,B)
res_beta_pdf_n_22_method3 = rep(0,B)
res_beta_pdf_n_32_method3 = rep(0,B)
res_beta_pdf_n_42_method3 = rep(0,B)
res_beta_pdf_n_52_method3 = rep(0,B)

# EM #

res_beta_pdf_n_11_EM = rep(0,B)
res_beta_pdf_n_21_EM = rep(0,B)
res_beta_pdf_n_31_EM = rep(0,B)
res_beta_pdf_n_41_EM = rep(0,B)
res_beta_pdf_n_51_EM = rep(0,B)

res_beta_pdf_n_12_EM = rep(0,B)
res_beta_pdf_n_22_EM = rep(0,B)
res_beta_pdf_n_32_EM = rep(0,B)
res_beta_pdf_n_42_EM = rep(0,B)
res_beta_pdf_n_52_EM = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  # Beta Distribution - n1
  
  x   <- sort(rbeta(n[1],alpha,beta))
  n_i <- hist(x, breaks = beta_grid, plot = F)$counts
  
  res_beta_pdf_n_11_method1[i] = L1_Distance_method1(n_i,beta_grid,pdf_beta,range_beta)
  res_beta_pdf_n_12_method1[i] = L2_Distance_method1(n_i,beta_grid,pdf_beta,range_beta)
  
  res_beta_pdf_n_11_method2[i] = L1_Distance_method2(n_i,beta_grid,pdf_beta,range_beta)
  res_beta_pdf_n_12_method2[i] = L2_Distance_method2(n_i,beta_grid,pdf_beta,range_beta)
  
  res_beta_pdf_n_11_method3[i] = L1_Distance_method3(n_i,beta_grid,pdf_beta,range_beta)
  res_beta_pdf_n_12_method3[i] = L2_Distance_method3(n_i,beta_grid,pdf_beta,range_beta)
  
  res_beta_pdf_n_11_EM[i] = EM_L1_L2(n_i,beta_grid,pdf_beta,range_beta)$L1
  res_beta_pdf_n_12_EM[i] = EM_L1_L2(n_i,beta_grid,pdf_beta,range_beta)$L2
  
  # Beta Distribution - n2
  
  x   <- sort(rbeta(n[2],alpha,beta))
  n_i <- hist(x, breaks = beta_grid, plot = F)$counts
  
  res_beta_pdf_n_21_method1[i] = L1_Distance_method1(n_i,beta_grid,pdf_beta,range_beta)
  res_beta_pdf_n_22_method1[i] = L2_Distance_method1(n_i,beta_grid,pdf_beta,range_beta)
  
  res_beta_pdf_n_21_method2[i] = L1_Distance_method2(n_i,beta_grid,pdf_beta,range_beta)
  res_beta_pdf_n_22_method2[i] = L2_Distance_method2(n_i,beta_grid,pdf_beta,range_beta)
  
  res_beta_pdf_n_21_method3[i] = L1_Distance_method3(n_i,beta_grid,pdf_beta,range_beta)
  res_beta_pdf_n_22_method3[i] = L2_Distance_method3(n_i,beta_grid,pdf_beta,range_beta)
  
  res_beta_pdf_n_21_EM[i] = EM_L1_L2(n_i,beta_grid,pdf_beta,range_beta)$L1
  res_beta_pdf_n_22_EM[i] = EM_L1_L2(n_i,beta_grid,pdf_beta,range_beta)$L2
  
  # Beta Distribution - n3
  
  x   <- sort(rbeta(n[3],alpha,beta))
  n_i <- hist(x, breaks = beta_grid, plot = F)$counts
  
  res_beta_pdf_n_31_method1[i] = L1_Distance_method1(n_i,beta_grid,pdf_beta,range_beta)
  res_beta_pdf_n_32_method1[i] = L2_Distance_method1(n_i,beta_grid,pdf_beta,range_beta)

  res_beta_pdf_n_31_method2[i] = L1_Distance_method2(n_i,beta_grid,pdf_beta,range_beta)
  res_beta_pdf_n_32_method2[i] = L2_Distance_method2(n_i,beta_grid,pdf_beta,range_beta)
  
  res_beta_pdf_n_31_method3[i] = L1_Distance_method3(n_i,beta_grid,pdf_beta,range_beta)
  res_beta_pdf_n_32_method3[i] = L2_Distance_method3(n_i,beta_grid,pdf_beta,range_beta)
  
  res_beta_pdf_n_31_EM[i] = EM_L1_L2(n_i,beta_grid,pdf_beta,range_beta)$L1
  res_beta_pdf_n_32_EM[i] = EM_L1_L2(n_i,beta_grid,pdf_beta,range_beta)$L2
  
  # Beta Distribution - n4
  
  x   <- sort(rbeta(n[4],alpha,beta))
  n_i <- hist(x, breaks = beta_grid, plot = F)$counts
  
  res_beta_pdf_n_41_method1[i] = L1_Distance_method1(n_i,beta_grid,pdf_beta,range_beta)
  res_beta_pdf_n_42_method1[i] = L2_Distance_method1(n_i,beta_grid,pdf_beta,range_beta)
  
  res_beta_pdf_n_41_method2[i] = L1_Distance_method2(n_i,beta_grid,pdf_beta,range_beta)
  res_beta_pdf_n_42_method2[i] = L2_Distance_method2(n_i,beta_grid,pdf_beta,range_beta)
  
  res_beta_pdf_n_41_method3[i] = L1_Distance_method3(n_i,beta_grid,pdf_beta,range_beta)
  res_beta_pdf_n_42_method3[i] = L2_Distance_method3(n_i,beta_grid,pdf_beta,range_beta)
  
  res_beta_pdf_n_41_EM[i] = EM_L1_L2(n_i,beta_grid,pdf_beta,range_beta)$L1
  res_beta_pdf_n_42_EM[i] = EM_L1_L2(n_i,beta_grid,pdf_beta,range_beta)$L2
  
  # Beta Distribution - n5
  
  x   <- sort(rbeta(n[5],alpha,beta))
  n_i <- hist(x, breaks = beta_grid, plot = F)$counts
  
  res_beta_pdf_n_51_method1[i] = L1_Distance_method1(n_i,beta_grid,pdf_beta,range_beta)
  res_beta_pdf_n_52_method1[i] = L2_Distance_method1(n_i,beta_grid,pdf_beta,range_beta)
  
  res_beta_pdf_n_51_method2[i] = L1_Distance_method2(n_i,beta_grid,pdf_beta,range_beta)
  res_beta_pdf_n_52_method2[i] = L2_Distance_method2(n_i,beta_grid,pdf_beta,range_beta)
  
  res_beta_pdf_n_51_method3[i] = L1_Distance_method3(n_i,beta_grid,pdf_beta,range_beta)
  res_beta_pdf_n_52_method3[i] = L2_Distance_method3(n_i,beta_grid,pdf_beta,range_beta)
  
  res_beta_pdf_n_51_EM[i] = EM_L1_L2(n_i,beta_grid,pdf_beta,range_beta)$L1
  res_beta_pdf_n_52_EM[i] = EM_L1_L2(n_i,beta_grid,pdf_beta,range_beta)$L2
  
  print(i)
  
}

pdf_L1_beta_limit = 0
pdf_L2_beta_limit = 0

# Method-1 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L1-n/L1-PDF-1-Beta-n-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pdf_n_11_method1,res_beta_pdf_n_21_method1,res_beta_pdf_n_31_method1,res_beta_pdf_n_41_method1,res_beta_pdf_n_51_method1,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2,ylim=c(0,0.8))
# segments(x0 = 0.29, y0 = pdf_L1_beta_limit, x1 = 10^5, y1 = pdf_L1_beta_limit, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L2-n/L2-PDF-1-Beta-n-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pdf_n_12_method1,res_beta_pdf_n_22_method1,res_beta_pdf_n_32_method1,res_beta_pdf_n_42_method1,res_beta_pdf_n_52_method1,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2,ylim=c(0,1.2))
# segments(x0 = 0.29, y0 = pdf_L2_beta_limit, x1 = 10^5, y1 = pdf_L2_beta_limit, col = "red", lty = 2)
dev.off()

# Method-2 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L1-n/L1-PDF-1-Beta-n-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pdf_n_11_method2,res_beta_pdf_n_21_method2,res_beta_pdf_n_31_method2,res_beta_pdf_n_41_method2,res_beta_pdf_n_51_method2,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2,ylim=c(0,0.8))
segments(x0 = 0.29, y0 = pdf_L1_beta_limit, x1 = 10^5, y1 = pdf_L1_beta_limit, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L2-n/L2-PDF-1-Beta-n-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pdf_n_12_method2,res_beta_pdf_n_22_method2,res_beta_pdf_n_32_method2,res_beta_pdf_n_42_method2,res_beta_pdf_n_52_method2,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2,ylim=c(0,1.2))
# segments(x0 = 0.29, y0 = pdf_L2_beta_limit, x1 = 10^5, y1 = pdf_L2_beta_limit, col = "red", lty = 2)
dev.off()

# Method-3 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L1-n/L1-PDF-1-Beta-n-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pdf_n_11_method3,res_beta_pdf_n_21_method3,res_beta_pdf_n_31_method3,res_beta_pdf_n_41_method3,res_beta_pdf_n_51_method3,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2,ylim=c(0,0.8))
# segments(x0 = 0.29, y0 = pdf_L1_beta_limit, x1 = 10^5, y1 = pdf_L1_beta_limit, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L2-n/L2-PDF-1-Beta-n-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pdf_n_12_method3,res_beta_pdf_n_22_method3,res_beta_pdf_n_32_method3,res_beta_pdf_n_42_method3,res_beta_pdf_n_52_method3,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2,ylim=c(0,1.2))
# segments(x0 = 0.29, y0 = pdf_L2_beta_limit, x1 = 10^5, y1 = pdf_L2_beta_limit, col = "red", lty = 2)
dev.off()

# EM # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L1-n/L1-PDF-1-Beta-n-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pdf_n_11_EM,res_beta_pdf_n_21_EM,res_beta_pdf_n_31_EM,res_beta_pdf_n_41_EM,res_beta_pdf_n_51_EM,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2,ylim=c(0,0.8))
# segments(x0 = 0.29, y0 = pdf_L1_beta_limit, x1 = 10^5, y1 = pdf_L1_beta_limit, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L2-n/L2-PDF-1-Beta-n-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pdf_n_12_EM,res_beta_pdf_n_22_EM,res_beta_pdf_n_32_EM,res_beta_pdf_n_42_EM,res_beta_pdf_n_52_EM,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors2,ylim=c(0,1.2))
# segments(x0 = 0.29, y0 = pdf_L2_beta_limit, x1 = 10^5, y1 = pdf_L2_beta_limit, col = "red", lty = 2)
dev.off()

# Gamma Distribution - Sample Size Simulation 

MM = qgamma(0.9999999,alpha,beta) + 3
gamma_grid <- seq(0,MM, by = gamma_deltas[1])

# Method-1 #

res_gamma_pdf_n_11_method1 = rep(0,B)
res_gamma_pdf_n_21_method1 = rep(0,B)
res_gamma_pdf_n_31_method1 = rep(0,B)
res_gamma_pdf_n_41_method1 = rep(0,B)
res_gamma_pdf_n_51_method1 = rep(0,B)

res_gamma_pdf_n_12_method1 = rep(0,B)
res_gamma_pdf_n_22_method1 = rep(0,B)
res_gamma_pdf_n_32_method1 = rep(0,B)
res_gamma_pdf_n_42_method1 = rep(0,B)
res_gamma_pdf_n_52_method1 = rep(0,B)

# Method-2 #

res_gamma_pdf_n_11_method2 = rep(0,B)
res_gamma_pdf_n_21_method2 = rep(0,B)
res_gamma_pdf_n_31_method2 = rep(0,B)
res_gamma_pdf_n_41_method2 = rep(0,B)
res_gamma_pdf_n_51_method2 = rep(0,B)

res_gamma_pdf_n_12_method2 = rep(0,B)
res_gamma_pdf_n_22_method2 = rep(0,B)
res_gamma_pdf_n_32_method2 = rep(0,B)
res_gamma_pdf_n_42_method2 = rep(0,B)
res_gamma_pdf_n_52_method2 = rep(0,B)

# Method-3 #

res_gamma_pdf_n_11_method3 = rep(0,B)
res_gamma_pdf_n_21_method3 = rep(0,B)
res_gamma_pdf_n_31_method3 = rep(0,B)
res_gamma_pdf_n_41_method3 = rep(0,B)
res_gamma_pdf_n_51_method3 = rep(0,B)

res_gamma_pdf_n_12_method3 = rep(0,B)
res_gamma_pdf_n_22_method3 = rep(0,B)
res_gamma_pdf_n_32_method3 = rep(0,B)
res_gamma_pdf_n_42_method3 = rep(0,B)
res_gamma_pdf_n_52_method3 = rep(0,B)

# EM #

res_gamma_pdf_n_11_EM = rep(0,B)
res_gamma_pdf_n_21_EM = rep(0,B)
res_gamma_pdf_n_31_EM = rep(0,B)
res_gamma_pdf_n_41_EM = rep(0,B)
res_gamma_pdf_n_51_EM = rep(0,B)

res_gamma_pdf_n_12_EM = rep(0,B)
res_gamma_pdf_n_22_EM = rep(0,B)
res_gamma_pdf_n_32_EM = rep(0,B)
res_gamma_pdf_n_42_EM = rep(0,B)
res_gamma_pdf_n_52_EM = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  # Gamma Distribution - n1
  
  x   <- sort(rgamma(n[1],alpha,beta))
  n_i <- hist(x, breaks = gamma_grid, plot = F)$counts

  res_gamma_pdf_n_11_method1[i] = L1_Distance_method1(n_i,gamma_grid,pdf_gamma,range_gamma)
  res_gamma_pdf_n_12_method1[i] = L2_Distance_method1(n_i,gamma_grid,pdf_gamma,range_gamma)
  
  res_gamma_pdf_n_11_method2[i] = L1_Distance_method2(n_i,gamma_grid,pdf_gamma,range_gamma)
  res_gamma_pdf_n_12_method2[i] = L2_Distance_method2(n_i,gamma_grid,pdf_gamma,range_gamma)
  
  res_gamma_pdf_n_11_method3[i] = L1_Distance_method3(n_i,gamma_grid,pdf_gamma,range_gamma)
  res_gamma_pdf_n_12_method3[i] = L2_Distance_method3(n_i,gamma_grid,pdf_gamma,range_gamma)
  
  res_gamma_pdf_n_11_EM[i] = EM_L1_L2(n_i,gamma_grid,pdf_gamma,range_gamma)$L1
  res_gamma_pdf_n_12_EM[i] = EM_L1_L2(n_i,gamma_grid,pdf_gamma,range_gamma)$L2
  
  # Gamma Distribution - n2
  
  x   <- sort(rgamma(n[2],alpha,beta))
  n_i <- hist(x, breaks = gamma_grid, plot = F)$counts
  
  res_gamma_pdf_n_21_method1[i] = L1_Distance_method1(n_i,gamma_grid,pdf_gamma,range_gamma)
  res_gamma_pdf_n_22_method1[i] = L2_Distance_method1(n_i,gamma_grid,pdf_gamma,range_gamma)
  
  res_gamma_pdf_n_21_method2[i] = L1_Distance_method2(n_i,gamma_grid,pdf_gamma,range_gamma)
  res_gamma_pdf_n_22_method2[i] = L2_Distance_method2(n_i,gamma_grid,pdf_gamma,range_gamma)
  
  res_gamma_pdf_n_21_method3[i] = L1_Distance_method3(n_i,gamma_grid,pdf_gamma,range_gamma)
  res_gamma_pdf_n_22_method3[i] = L2_Distance_method3(n_i,gamma_grid,pdf_gamma,range_gamma)
  
  res_gamma_pdf_n_21_EM[i] = EM_L1_L2(n_i,gamma_grid,pdf_gamma,range_gamma)$L1
  res_gamma_pdf_n_22_EM[i] = EM_L1_L2(n_i,gamma_grid,pdf_gamma,range_gamma)$L2
  
  # Gamma Distribution - n3
  
  x   <- sort(rgamma(n[3],alpha,beta))
  n_i <- hist(x, breaks = gamma_grid, plot = F)$counts
  
  res_gamma_pdf_n_31_method1[i] = L1_Distance_method1(n_i,gamma_grid,pdf_gamma,range_gamma)
  res_gamma_pdf_n_32_method1[i] = L2_Distance_method1(n_i,gamma_grid,pdf_gamma,range_gamma)
  
  res_gamma_pdf_n_31_method2[i] = L1_Distance_method2(n_i,gamma_grid,pdf_gamma,range_gamma)
  res_gamma_pdf_n_32_method2[i] = L2_Distance_method2(n_i,gamma_grid,pdf_gamma,range_gamma)
  
  res_gamma_pdf_n_31_method3[i] = L1_Distance_method3(n_i,gamma_grid,pdf_gamma,range_gamma)
  res_gamma_pdf_n_32_method3[i] = L2_Distance_method3(n_i,gamma_grid,pdf_gamma,range_gamma)
  
  res_gamma_pdf_n_31_EM[i] = EM_L1_L2(n_i,gamma_grid,pdf_gamma,range_gamma)$L1
  res_gamma_pdf_n_32_EM[i] = EM_L1_L2(n_i,gamma_grid,pdf_gamma,range_gamma)$L2
  
  # Gamma Distribution - n4
  
  x   <- sort(rgamma(n[4],alpha,beta))
  n_i <- hist(x, breaks = gamma_grid, plot = F)$counts
  
  res_gamma_pdf_n_41_method1[i] = L1_Distance_method1(n_i,gamma_grid,pdf_gamma,range_gamma)
  res_gamma_pdf_n_42_method1[i] = L2_Distance_method1(n_i,gamma_grid,pdf_gamma,range_gamma)
  
  res_gamma_pdf_n_41_method2[i] = L1_Distance_method2(n_i,gamma_grid,pdf_gamma,range_gamma)
  res_gamma_pdf_n_42_method2[i] = L2_Distance_method2(n_i,gamma_grid,pdf_gamma,range_gamma)
  
  res_gamma_pdf_n_41_method3[i] = L1_Distance_method3(n_i,gamma_grid,pdf_gamma,range_gamma)
  res_gamma_pdf_n_42_method3[i] = L2_Distance_method3(n_i,gamma_grid,pdf_gamma,range_gamma)
  
  res_gamma_pdf_n_41_EM[i] = EM_L1_L2(n_i,gamma_grid,pdf_gamma,range_gamma)$L1
  res_gamma_pdf_n_42_EM[i] = EM_L1_L2(n_i,gamma_grid,pdf_gamma,range_gamma)$L2
  
  # Gamma Distribution - n5
  
  x   <- sort(rgamma(n[5],alpha,beta))
  n_i <- hist(x, breaks = gamma_grid, plot = F)$counts
  
  res_gamma_pdf_n_51_method1[i] = L1_Distance_method1(n_i,gamma_grid,pdf_gamma,range_gamma)
  res_gamma_pdf_n_52_method1[i] = L2_Distance_method1(n_i,gamma_grid,pdf_gamma,range_gamma)
  
  res_gamma_pdf_n_51_method2[i] = L1_Distance_method2(n_i,gamma_grid,pdf_gamma,range_gamma)
  res_gamma_pdf_n_52_method2[i] = L2_Distance_method2(n_i,gamma_grid,pdf_gamma,range_gamma)
  
  res_gamma_pdf_n_51_method3[i] = L1_Distance_method3(n_i,gamma_grid,pdf_gamma,range_gamma)
  res_gamma_pdf_n_52_method3[i] = L2_Distance_method3(n_i,gamma_grid,pdf_gamma,range_gamma)
  
  res_gamma_pdf_n_51_EM[i] = EM_L1_L2(n_i,gamma_grid,pdf_gamma,range_gamma)$L1
  res_gamma_pdf_n_52_EM[i] = EM_L1_L2(n_i,gamma_grid,pdf_gamma,range_gamma)$L2
  
  print(i)
  
}

pdf_L1_gamma_limit = 0
pdf_L2_gamma_limit = 0

# Method-1 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L1-n/L1-PDF-1-Gamma-n-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pdf_n_11_method1,res_gamma_pdf_n_21_method1,res_gamma_pdf_n_31_method1,res_gamma_pdf_n_41_method1,res_gamma_pdf_n_51_method1,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3,ylim=c(0,1.2))
# segments(x0 = 0.29, y0 = pdf_L1_gamma_limit, x1 = 10^5, y1 = pdf_L1_gamma_limit, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L2-n/L2-PDF-1-Gamma-n-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pdf_n_12_method1,res_gamma_pdf_n_22_method1,res_gamma_pdf_n_32_method1,res_gamma_pdf_n_42_method1,res_gamma_pdf_n_52_method1,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3,ylim=c(0,1.3))
# segments(x0 = 0.29, y0 = pdf_L2_gamma_limit, x1 = 10^5, y1 = pdf_L2_gamma_limit, col = "red", lty = 2)
dev.off()

# Method-2 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L1-n/L1-PDF-1-Gamma-n-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pdf_n_11_method2,res_gamma_pdf_n_21_method2,res_gamma_pdf_n_31_method2,res_gamma_pdf_n_41_method2,res_gamma_pdf_n_51_method2,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3,ylim=c(0,1.2))
# segments(x0 = 0.29, y0 = pdf_L1_gamma_limit, x1 = 10^5, y1 = pdf_L1_gamma_limit, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L2-n/L2-PDF-1-Gamma-n-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pdf_n_12_method2,res_gamma_pdf_n_22_method2,res_gamma_pdf_n_32_method2,res_gamma_pdf_n_42_method2,res_gamma_pdf_n_52_method2,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3,ylim=c(0,1.3))
#segments(x0 = 0.29, y0 = pdf_L2_gamma_limit, x1 = 10^5, y1 = pdf_L2_gamma_limit, col = "red", lty = 2)
dev.off()

# Method-3 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L1-n/L1-PDF-1-Gamma-n-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pdf_n_11_method3,res_gamma_pdf_n_21_method3,res_gamma_pdf_n_31_method3,res_gamma_pdf_n_41_method3,res_gamma_pdf_n_51_method3,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3,ylim=c(0,1.2))
#segments(x0 = 0.29, y0 = pdf_L1_gamma_limit, x1 = 10^5, y1 = pdf_L1_gamma_limit, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L2-n/L2-PDF-1-Gamma-n-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pdf_n_12_method3,res_gamma_pdf_n_22_method3,res_gamma_pdf_n_32_method3,res_gamma_pdf_n_42_method3,res_gamma_pdf_n_52_method3,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3,ylim=c(0,1.3))
# segments(x0 = 0.29, y0 = pdf_L2_gamma_limit, x1 = 10^5, y1 = pdf_L2_gamma_limit, col = "red", lty = 2)
dev.off()

# EM # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L1-n/L1-PDF-1-Gamma-n-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pdf_n_11_EM,res_gamma_pdf_n_21_EM,res_gamma_pdf_n_31_EM,res_gamma_pdf_n_41_EM,res_gamma_pdf_n_51_EM,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3,ylim=c(0,1.2))
# segments(x0 = 0.29, y0 = pdf_L1_gamma_limit, x1 = 10^5, y1 = pdf_L1_gamma_limit, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L2-n/L2-PDF-1-Gamma-n-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pdf_n_12_EM,res_gamma_pdf_n_22_EM,res_gamma_pdf_n_32_EM,res_gamma_pdf_n_42_EM,res_gamma_pdf_n_52_EM,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors3,ylim=c(0,1.3))
# segments(x0 = 0.29, y0 = pdf_L2_gamma_limit, x1 = 10^5, y1 = pdf_L2_gamma_limit, col = "red", lty = 2)
dev.off()

# Logistic Distribution - Sample Size Simulation 

MM = qlogis(0.9999999,mu,sigma) + 10
logistic_grid <- seq(-MM,MM, by = logistic_deltas[1])

# Method-1 # 

res_logis_pdf_n_11_method1 = rep(0,B)
res_logis_pdf_n_21_method1 = rep(0,B)
res_logis_pdf_n_31_method1 = rep(0,B)
res_logis_pdf_n_41_method1 = rep(0,B)
res_logis_pdf_n_51_method1 = rep(0,B)

res_logis_pdf_n_12_method1 = rep(0,B)
res_logis_pdf_n_22_method1 = rep(0,B)
res_logis_pdf_n_32_method1 = rep(0,B)
res_logis_pdf_n_42_method1 = rep(0,B)
res_logis_pdf_n_52_method1 = rep(0,B)

# Method-2 # 

res_logis_pdf_n_11_method2 = rep(0,B)
res_logis_pdf_n_21_method2 = rep(0,B)
res_logis_pdf_n_31_method2 = rep(0,B)
res_logis_pdf_n_41_method2 = rep(0,B)
res_logis_pdf_n_51_method2 = rep(0,B)

res_logis_pdf_n_12_method2 = rep(0,B)
res_logis_pdf_n_22_method2 = rep(0,B)
res_logis_pdf_n_32_method2 = rep(0,B)
res_logis_pdf_n_42_method2 = rep(0,B)
res_logis_pdf_n_52_method2 = rep(0,B)

# Method-3 # 

res_logis_pdf_n_11_method3 = rep(0,B)
res_logis_pdf_n_21_method3 = rep(0,B)
res_logis_pdf_n_31_method3 = rep(0,B)
res_logis_pdf_n_41_method3 = rep(0,B)
res_logis_pdf_n_51_method3 = rep(0,B)

res_logis_pdf_n_12_method3 = rep(0,B)
res_logis_pdf_n_22_method3 = rep(0,B)
res_logis_pdf_n_32_method3 = rep(0,B)
res_logis_pdf_n_42_method3 = rep(0,B)
res_logis_pdf_n_52_method3 = rep(0,B)

# EM # 

res_logis_pdf_n_11_EM = rep(0,B)
res_logis_pdf_n_21_EM = rep(0,B)
res_logis_pdf_n_31_EM = rep(0,B)
res_logis_pdf_n_41_EM = rep(0,B)
res_logis_pdf_n_51_EM = rep(0,B)

res_logis_pdf_n_12_EM = rep(0,B)
res_logis_pdf_n_22_EM = rep(0,B)
res_logis_pdf_n_32_EM = rep(0,B)
res_logis_pdf_n_42_EM = rep(0,B)
res_logis_pdf_n_52_EM = rep(0,B)

set.seed(10)
set.seed(11)
set.seed(12)

for (i in 94:B) {
  
  # Logistic Distribution - n1
  
  x = sort(rlogis(n[1]))
  n_i <- hist(x, breaks = logistic_grid, plot = F)$counts
  
  res_logis_pdf_n_11_method1[i] = L1_Distance_method1(n_i,logistic_grid,pdf_logistic,range_logistic)
  res_logis_pdf_n_12_method1[i] = L2_Distance_method1(n_i,logistic_grid,pdf_logistic,range_logistic)
  
  res_logis_pdf_n_11_method2[i] = L1_Distance_method2(n_i,logistic_grid,pdf_logistic,range_logistic)
  res_logis_pdf_n_12_method2[i] = L2_Distance_method2(n_i,logistic_grid,pdf_logistic,range_logistic)
  
  res_logis_pdf_n_11_method3[i] = L1_Distance_method3(n_i,logistic_grid,pdf_logistic,range_logistic)
  res_logis_pdf_n_12_method3[i] = L2_Distance_method3(n_i,logistic_grid,pdf_logistic,range_logistic)
  
  res_logis_pdf_n_11_EM[i] = EM_L1_L2(n_i,logistic_grid,pdf_logistic,range_logistic)$L1
  res_logis_pdf_n_12_EM[i] = EM_L1_L2(n_i,logistic_grid,pdf_logistic,range_logistic)$L2
  
  # Logistic Distribution - n2
  
  x = sort(rlogis(n[2]))
  n_i <- hist(x, breaks = logistic_grid, plot = F)$counts
  
  res_logis_pdf_n_21_method1[i] = L1_Distance_method1(n_i,logistic_grid,pdf_logistic,range_logistic)
  res_logis_pdf_n_22_method1[i] = L2_Distance_method1(n_i,logistic_grid,pdf_logistic,range_logistic)
  
  res_logis_pdf_n_21_method2[i] = L1_Distance_method2(n_i,logistic_grid,pdf_logistic,range_logistic)
  res_logis_pdf_n_22_method2[i] = L2_Distance_method2(n_i,logistic_grid,pdf_logistic,range_logistic)
  
  res_logis_pdf_n_21_method3[i] = L1_Distance_method3(n_i,logistic_grid,pdf_logistic,range_logistic)
  res_logis_pdf_n_22_method3[i] = L2_Distance_method3(n_i,logistic_grid,pdf_logistic,range_logistic)
  
  res_logis_pdf_n_21_EM[i] = EM_L1_L2(n_i,logistic_grid,pdf_logistic,range_logistic)$L1
  res_logis_pdf_n_22_EM[i] = EM_L1_L2(n_i,logistic_grid,pdf_logistic,range_logistic)$L2
  
  # Logistic Distribution - n3
  
  x = sort(rlogis(n[3]))
  n_i <- hist(x, breaks = logistic_grid, plot = F)$counts
  
  res_logis_pdf_n_31_method1[i] = L1_Distance_method1(n_i,logistic_grid,pdf_logistic,range_logistic)
  res_logis_pdf_n_32_method1[i] = L2_Distance_method1(n_i,logistic_grid,pdf_logistic,range_logistic)
  
  res_logis_pdf_n_31_method2[i] = L1_Distance_method2(n_i,logistic_grid,pdf_logistic,range_logistic)
  res_logis_pdf_n_32_method2[i] = L2_Distance_method2(n_i,logistic_grid,pdf_logistic,range_logistic)
  
  res_logis_pdf_n_31_method3[i] = L1_Distance_method3(n_i,logistic_grid,pdf_logistic,range_logistic)
  res_logis_pdf_n_32_method3[i] = L2_Distance_method3(n_i,logistic_grid,pdf_logistic,range_logistic)
  
  res_logis_pdf_n_31_EM[i] = EM_L1_L2(n_i,logistic_grid,pdf_logistic,range_logistic)$L1
  res_logis_pdf_n_32_EM[i] = EM_L1_L2(n_i,logistic_grid,pdf_logistic,range_logistic)$L2
  
  # Logistic Distribution - n4
  
  x = sort(rlogis(n[4]))
  n_i <- hist(x, breaks = logistic_grid, plot = F)$counts
  
  res_logis_pdf_n_41_method1[i] = L1_Distance_method1(n_i,logistic_grid,pdf_logistic,range_logistic)
  res_logis_pdf_n_42_method1[i] = L2_Distance_method1(n_i,logistic_grid,pdf_logistic,range_logistic)
  
  res_logis_pdf_n_41_method2[i] = L1_Distance_method2(n_i,logistic_grid,pdf_logistic,range_logistic)
  res_logis_pdf_n_42_method2[i] = L2_Distance_method2(n_i,logistic_grid,pdf_logistic,range_logistic)
  
  res_logis_pdf_n_41_method3[i] = L1_Distance_method3(n_i,logistic_grid,pdf_logistic,range_logistic)
  res_logis_pdf_n_42_method3[i] = L2_Distance_method3(n_i,logistic_grid,pdf_logistic,range_logistic)
  
  res_logis_pdf_n_41_EM[i] = EM_L1_L2(n_i,logistic_grid,pdf_logistic,range_logistic)$L1
  res_logis_pdf_n_42_EM[i] = EM_L1_L2(n_i,logistic_grid,pdf_logistic,range_logistic)$L2
  
  # Logistic Distribution - n5
  
  x = sort(rlogis(n[5]))
  n_i <- hist(x, breaks = logistic_grid, plot = F)$counts
  
  res_logis_pdf_n_51_method1[i] = L1_Distance_method1(n_i,logistic_grid,pdf_logistic,range_logistic)
  res_logis_pdf_n_52_method1[i] = L2_Distance_method1(n_i,logistic_grid,pdf_logistic,range_logistic)
  
  res_logis_pdf_n_51_method2[i] = L1_Distance_method2(n_i,logistic_grid,pdf_logistic,range_logistic)
  res_logis_pdf_n_52_method2[i] = L2_Distance_method2(n_i,logistic_grid,pdf_logistic,range_logistic)
  
  res_logis_pdf_n_51_method3[i] = L1_Distance_method3(n_i,logistic_grid,pdf_logistic,range_logistic)
  res_logis_pdf_n_52_method3[i] = L2_Distance_method3(n_i,logistic_grid,pdf_logistic,range_logistic)
  
  res_logis_pdf_n_51_EM[i] = EM_L1_L2(n_i,logistic_grid,pdf_logistic,range_logistic)$L1
  res_logis_pdf_n_52_EM[i] = EM_L1_L2(n_i,logistic_grid,pdf_logistic,range_logistic)$L2
  
  print(i)
  
}

# Method-1 #

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L1-n/L1-PDF-1-Logistic-n-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pdf_n_11_method1,res_logis_pdf_n_21_method1,res_logis_pdf_n_31_method1,res_logis_pdf_n_41_method1,res_logis_pdf_n_51_method1,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors4,ylim=c(0,1.2))
# segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L2-n/L2-PDF-1-Logistic-n-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pdf_n_12_method1,res_logis_pdf_n_22_method1,res_logis_pdf_n_32_method1,res_logis_pdf_n_42_method1,res_logis_pdf_n_52_method1,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors4,ylim=c(0,0.6))
# segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

# Method-2 #

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L1-n/L1-PDF-1-Logistic-n-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pdf_n_11_method2,res_logis_pdf_n_21_method2,res_logis_pdf_n_31_method2,res_logis_pdf_n_41_method2,res_logis_pdf_n_51_method2,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors4,ylim=c(0,1.2))
# segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L2-n/L2-PDF-1-Logistic-n-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pdf_n_12_method2,res_logis_pdf_n_22_method2,res_logis_pdf_n_32_method2,res_logis_pdf_n_42_method2,res_logis_pdf_n_52_method2,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors4,ylim=c(0,0.6))
# segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

# Method-3 #

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L1-n/L1-PDF-1-Logistic-n-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pdf_n_11_method3,res_logis_pdf_n_21_method3,res_logis_pdf_n_31_method3,res_logis_pdf_n_41_method3,res_logis_pdf_n_51_method3,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors4,ylim=c(0,1.2))
# segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L2-n/L2-PDF-1-Logistic-n-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pdf_n_12_method3,res_logis_pdf_n_22_method3,res_logis_pdf_n_32_method3,res_logis_pdf_n_42_method3,res_logis_pdf_n_52_method3,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors4,ylim=c(0,0.6))
# segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

# EM #

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L1-n/L1-PDF-1-Logistic-n-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pdf_n_11_EM,res_logis_pdf_n_21_EM,res_logis_pdf_n_31_EM,res_logis_pdf_n_41_EM,res_logis_pdf_n_51_EM,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors4,ylim=c(0,1.2))
# segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L2-n/L2-PDF-1-Logistic-n-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pdf_n_12_EM,res_logis_pdf_n_22_EM,res_logis_pdf_n_32_EM,res_logis_pdf_n_42_EM,res_logis_pdf_n_52_EM,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors4,ylim=c(0,0.6))
# segments(x0 = 0.29, y0 = 0, x1 = 10^5, y1 = 0, col = "red", lty = 2)
dev.off()

# Student's t Distribution - Sample Size Simulation 

MM <- ceiling(qt(0.99999999, df))+50

t_grid <- seq(-MM,MM, by = t_deltas[1])

# Method-1 # 

res_t_pdf_n_11_method1 = rep(0,B)
res_t_pdf_n_21_method1 = rep(0,B)
res_t_pdf_n_31_method1 = rep(0,B)
res_t_pdf_n_41_method1 = rep(0,B)
res_t_pdf_n_51_method1 = rep(0,B)

res_t_pdf_n_12_method1 = rep(0,B)
res_t_pdf_n_22_method1 = rep(0,B)
res_t_pdf_n_32_method1 = rep(0,B)
res_t_pdf_n_42_method1 = rep(0,B)
res_t_pdf_n_52_method1 = rep(0,B)

# Method-2 # 

res_t_pdf_n_11_method2 = rep(0,B)
res_t_pdf_n_21_method2 = rep(0,B)
res_t_pdf_n_31_method2 = rep(0,B)
res_t_pdf_n_41_method2 = rep(0,B)
res_t_pdf_n_51_method2 = rep(0,B)

res_t_pdf_n_12_method2 = rep(0,B)
res_t_pdf_n_22_method2 = rep(0,B)
res_t_pdf_n_32_method2 = rep(0,B)
res_t_pdf_n_42_method2 = rep(0,B)
res_t_pdf_n_52_method2 = rep(0,B)

# Method-3 # 

res_t_pdf_n_11_method3 = rep(0,B)
res_t_pdf_n_21_method3 = rep(0,B)
res_t_pdf_n_31_method3 = rep(0,B)
res_t_pdf_n_41_method3 = rep(0,B)
res_t_pdf_n_51_method3 = rep(0,B)

res_t_pdf_n_12_method3 = rep(0,B)
res_t_pdf_n_22_method3 = rep(0,B)
res_t_pdf_n_32_method3 = rep(0,B)
res_t_pdf_n_42_method3 = rep(0,B)
res_t_pdf_n_52_method3 = rep(0,B)

# EM # 

res_t_pdf_n_11_EM = rep(0,B)
res_t_pdf_n_21_EM = rep(0,B)
res_t_pdf_n_31_EM = rep(0,B)
res_t_pdf_n_41_EM = rep(0,B)
res_t_pdf_n_51_EM = rep(0,B)

res_t_pdf_n_12_EM = rep(0,B)
res_t_pdf_n_22_EM = rep(0,B)
res_t_pdf_n_32_EM = rep(0,B)
res_t_pdf_n_42_EM = rep(0,B)
res_t_pdf_n_52_EM = rep(0,B)

set.seed(11)
set.seed(12)

for (i in 99:B) {
  
  # t Distribution - n1
  
  x = sort(rt(n[1],df = df))
  n_i <- hist(x, breaks = t_grid, plot = F)$counts
  
  res_t_pdf_n_11_method1[i] = L1_Distance_method1(n_i,t_grid,pdf_t,range_t)
  res_t_pdf_n_12_method1[i] = L2_Distance_method1(n_i,t_grid,pdf_t,range_t)
  
  res_t_pdf_n_11_method2[i] = L1_Distance_method2(n_i,t_grid,pdf_t,range_t)
  res_t_pdf_n_12_method2[i] = L2_Distance_method2(n_i,t_grid,pdf_t,range_t)
  
  res_t_pdf_n_11_method3[i] = L1_Distance_method3(n_i,t_grid,pdf_t,range_t)
  res_t_pdf_n_12_method3[i] = L2_Distance_method3(n_i,t_grid,pdf_t,range_t)
  
  res_t_pdf_n_11_EM[i] = EM_L1_L2(n_i,t_grid,pdf_t,range_t)$L1
  res_t_pdf_n_12_EM[i] = EM_L1_L2(n_i,t_grid,pdf_t,range_t)$L2
  
  # t Distribution - n2
  
  x = sort(rt(n[2],df = df))
  n_i <- hist(x, breaks = t_grid, plot = F)$counts
  
  res_t_pdf_n_21_method1[i] = L1_Distance_method1(n_i,t_grid,pdf_t,range_t)
  res_t_pdf_n_22_method1[i] = L2_Distance_method1(n_i,t_grid,pdf_t,range_t)
  
  res_t_pdf_n_21_method2[i] = L1_Distance_method2(n_i,t_grid,pdf_t,range_t)
  res_t_pdf_n_22_method2[i] = L2_Distance_method2(n_i,t_grid,pdf_t,range_t)
  
  res_t_pdf_n_21_method3[i] = L1_Distance_method3(n_i,t_grid,pdf_t,range_t)
  res_t_pdf_n_22_method3[i] = L2_Distance_method3(n_i,t_grid,pdf_t,range_t)
  
  res_t_pdf_n_21_EM[i] = EM_L1_L2(n_i,t_grid,pdf_t,range_t)$L1
  res_t_pdf_n_22_EM[i] = EM_L1_L2(n_i,t_grid,pdf_t,range_t)$L2
  
  # t Distribution - n3
  
  x = sort(rt(n[3],df = df))
  n_i <- hist(x, breaks = t_grid, plot = F)$counts
  
  res_t_pdf_n_31_method1[i] = L1_Distance_method1(n_i,t_grid,pdf_t,range_t)
  res_t_pdf_n_32_method1[i] = L2_Distance_method1(n_i,t_grid,pdf_t,range_t)
  
  res_t_pdf_n_31_method2[i] = L1_Distance_method2(n_i,t_grid,pdf_t,range_t)
  res_t_pdf_n_32_method2[i] = L2_Distance_method2(n_i,t_grid,pdf_t,range_t)
  
  res_t_pdf_n_31_method3[i] = L1_Distance_method3(n_i,t_grid,pdf_t,range_t)
  res_t_pdf_n_32_method3[i] = L2_Distance_method3(n_i,t_grid,pdf_t,range_t)
  
  res_t_pdf_n_31_EM[i] = EM_L1_L2(n_i,t_grid,pdf_t,range_t)$L1
  res_t_pdf_n_32_EM[i] = EM_L1_L2(n_i,t_grid,pdf_t,range_t)$L2
  
  # t Distribution - n4
  
  x = sort(rt(n[4],df = df))
  n_i <- hist(x, breaks = t_grid, plot = F)$counts
  
  res_t_pdf_n_41_method1[i] = L1_Distance_method1(n_i,t_grid,pdf_t,range_t)
  res_t_pdf_n_42_method1[i] = L2_Distance_method1(n_i,t_grid,pdf_t,range_t)
  
  res_t_pdf_n_41_method2[i] = L1_Distance_method2(n_i,t_grid,pdf_t,range_t)
  res_t_pdf_n_42_method2[i] = L2_Distance_method2(n_i,t_grid,pdf_t,range_t)
  
  res_t_pdf_n_41_method3[i] = L1_Distance_method3(n_i,t_grid,pdf_t,range_t)
  res_t_pdf_n_42_method3[i] = L2_Distance_method3(n_i,t_grid,pdf_t,range_t)
  
  res_t_pdf_n_41_EM[i] = EM_L1_L2(n_i,t_grid,pdf_t,range_t)$L1
  res_t_pdf_n_42_EM[i] = EM_L1_L2(n_i,t_grid,pdf_t,range_t)$L2
  
  # t Distribution - n5
  
  x = sort(rt(n[5],df = df))
  n_i <- hist(x, breaks = t_grid, plot = F)$counts
  
  res_t_pdf_n_51_method1[i] = L1_Distance_method1(n_i,t_grid,pdf_t,range_t)
  res_t_pdf_n_52_method1[i] = L2_Distance_method1(n_i,t_grid,pdf_t,range_t)
  
  res_t_pdf_n_51_method2[i] = L1_Distance_method2(n_i,t_grid,pdf_t,range_t)
  res_t_pdf_n_52_method2[i] = L2_Distance_method2(n_i,t_grid,pdf_t,range_t)
  
  res_t_pdf_n_51_method3[i] = L1_Distance_method3(n_i,t_grid,pdf_t,range_t)
  res_t_pdf_n_52_method3[i] = L2_Distance_method3(n_i,t_grid,pdf_t,range_t)
  
  res_t_pdf_n_51_EM[i] = EM_L1_L2(n_i,t_grid,pdf_t,range_t)$L1
  res_t_pdf_n_52_EM[i] = EM_L1_L2(n_i,t_grid,pdf_t,range_t)$L2
  
  print(i)
  
}

pdf_L1_t_limit = 0
pdf_L2_t_limit = 0

# Method-1 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L1-n/L1-PDF-1-t-n-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pdf_n_11_method1,res_t_pdf_n_21_method1,res_t_pdf_n_31_method1,res_t_pdf_n_41_method1,res_t_pdf_n_51_method1,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors5,ylim=c(0,1.2))
# segments(x0 = 0.29, y0 = pdf_L1_t_limit, x1 = 10^5, y1 = pdf_L1_t_limit, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L2-n/L2-PDF-1-t-n-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pdf_n_12_method1,res_t_pdf_n_22_method1,res_t_pdf_n_32_method1,res_t_pdf_n_42_method1,res_t_pdf_n_52_method1,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors5,ylim=c(0,0.6))
# segments(x0 = 0.29, y0 = pdf_L2_t_limit, x1 = 10^5, y1 = pdf_L2_t_limit, col = "red", lty = 2)
dev.off()

# Method-2 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L1-n/L1-PDF-1-t-n-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pdf_n_11_method2,res_t_pdf_n_21_method2,res_t_pdf_n_31_method2,res_t_pdf_n_41_method2,res_t_pdf_n_51_method2,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors5,ylim=c(0,2.2))
# segments(x0 = 0.29, y0 = pdf_L1_t_limit, x1 = 10^5, y1 = pdf_L1_t_limit, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L2-n/L2-PDF-1-t-n-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pdf_n_12_method2,res_t_pdf_n_22_method2,res_t_pdf_n_32_method2,res_t_pdf_n_42_method2,res_t_pdf_n_52_method2,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors5,ylim=c(0,0.6))
# segments(x0 = 0.29, y0 = pdf_L2_t_limit, x1 = 10^5, y1 = pdf_L2_t_limit, col = "red", lty = 2)
dev.off()

# Method-3 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L1-n/L1-PDF-1-t-n-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pdf_n_11_method3,res_t_pdf_n_21_method3,res_t_pdf_n_31_method3,res_t_pdf_n_41_method3,res_t_pdf_n_51_method3,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors5,ylim=c(0,1.2))
# segments(x0 = 0.29, y0 = pdf_L1_t_limit, x1 = 10^5, y1 = pdf_L1_t_limit, col = "red", lty = 2)
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L2-n/L2-PDF-1-t-n-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pdf_n_12_method3,res_t_pdf_n_22_method3,res_t_pdf_n_32_method3,res_t_pdf_n_42_method3,res_t_pdf_n_52_method3,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors5,ylim=c(0,0.6))
# segments(x0 = 0.29, y0 = pdf_L2_t_limit, x1 = 10^5, y1 = pdf_L2_t_limit, col = "red", lty = 2)
dev.off()

# EM # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L1-n/L1-PDF-1-t-n-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pdf_n_11_EM,res_t_pdf_n_21_EM,res_t_pdf_n_31_EM,res_t_pdf_n_41_EM,res_t_pdf_n_51_EM,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors5,ylim=c(0,1.2))
# segments(x0 = 0.29, y0 = pdf_L1_t_limit, x1 = 10^5, y1 = pdf_L1_t_limit, col = "red", lty = 2)
dev.off()


grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L2-n/L2-PDF-1-t-n-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pdf_n_12_EM,res_t_pdf_n_22_EM,res_t_pdf_n_32_EM,res_t_pdf_n_42_EM,res_t_pdf_n_52_EM,
        names = c(expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5)),
        col = colors5,ylim=c(0,0.6))
# segments(x0 = 0.29, y0 = pdf_L2_t_limit, x1 = 10^5, y1 = pdf_L2_t_limit, col = "red", lty = 2)
dev.off()


n = c(10^3)
B = 500


# Normal Distribution - Grid-Width Simulation 

MM = qnorm(0.9999999,mu,sigma) + 3
norm_grids = list(seq(-MM,MM, by = norm_deltas[5]),
                  seq(-MM,MM, by = norm_deltas[4]),
                  seq(-MM,MM, by = norm_deltas[3]),
                  seq(-MM,MM, by = norm_deltas[2]),
                  seq(-MM,MM, by = norm_deltas[1]))

# Method-1 # 

res_norm_pdf_gw_11_method1 = rep(0,B)
res_norm_pdf_gw_21_method1 = rep(0,B)
res_norm_pdf_gw_31_method1 = rep(0,B)
res_norm_pdf_gw_41_method1 = rep(0,B)
res_norm_pdf_gw_51_method1 = rep(0,B)

res_norm_pdf_gw_12_method1 = rep(0,B)
res_norm_pdf_gw_22_method1 = rep(0,B)
res_norm_pdf_gw_32_method1 = rep(0,B)
res_norm_pdf_gw_42_method1 = rep(0,B)
res_norm_pdf_gw_52_method1 = rep(0,B)

# Method-2 # 

res_norm_pdf_gw_11_method2 = rep(0,B)
res_norm_pdf_gw_21_method2 = rep(0,B)
res_norm_pdf_gw_31_method2 = rep(0,B)
res_norm_pdf_gw_41_method2 = rep(0,B)
res_norm_pdf_gw_51_method2 = rep(0,B)

res_norm_pdf_gw_12_method2 = rep(0,B)
res_norm_pdf_gw_22_method2 = rep(0,B)
res_norm_pdf_gw_32_method2 = rep(0,B)
res_norm_pdf_gw_42_method2 = rep(0,B)
res_norm_pdf_gw_52_method2 = rep(0,B)

# Method-3 # 

res_norm_pdf_gw_11_method3 = rep(0,B)
res_norm_pdf_gw_21_method3 = rep(0,B)
res_norm_pdf_gw_31_method3 = rep(0,B)
res_norm_pdf_gw_41_method3 = rep(0,B)
res_norm_pdf_gw_51_method3 = rep(0,B)

res_norm_pdf_gw_12_method3 = rep(0,B)
res_norm_pdf_gw_22_method3 = rep(0,B)
res_norm_pdf_gw_32_method3 = rep(0,B)
res_norm_pdf_gw_42_method3 = rep(0,B)
res_norm_pdf_gw_52_method3 = rep(0,B)

# EM # 

res_norm_pdf_gw_11_EM = rep(0,B)
res_norm_pdf_gw_21_EM = rep(0,B)
res_norm_pdf_gw_31_EM = rep(0,B)
res_norm_pdf_gw_41_EM = rep(0,B)
res_norm_pdf_gw_51_EM = rep(0,B)

res_norm_pdf_gw_12_EM = rep(0,B)
res_norm_pdf_gw_22_EM = rep(0,B)
res_norm_pdf_gw_32_EM = rep(0,B)
res_norm_pdf_gw_42_EM = rep(0,B)
res_norm_pdf_gw_52_EM = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  x <- sort(rnorm(n, mean = mu,sd=sigma))
  
  # Normal Distribution - GW1
  
  n_i <- hist(x, breaks = norm_grids[[1]], plot = F)$counts
  
  res_norm_pdf_gw_11_method1[i] = L1_Distance_method1(n_i,norm_grids[[1]],pdf_norm,range_norm)
  res_norm_pdf_gw_12_method1[i] = L2_Distance_method1(n_i,norm_grids[[1]],pdf_norm,range_norm)
  
  res_norm_pdf_gw_11_method2[i] = L1_Distance_method2(n_i,norm_grids[[1]],pdf_norm,range_norm)
  res_norm_pdf_gw_12_method2[i] = L2_Distance_method2(n_i,norm_grids[[1]],pdf_norm,range_norm)
  
  res_norm_pdf_gw_11_method3[i] = L1_Distance_method3(n_i,norm_grids[[1]],pdf_norm,range_norm)
  res_norm_pdf_gw_12_method3[i] = L2_Distance_method3(n_i,norm_grids[[1]],pdf_norm,range_norm)
  
  res_norm_pdf_gw_11_EM[i] = EM_L1_L2(n_i,norm_grids[[1]],pdf_norm,range_norm)$L1
  res_norm_pdf_gw_12_EM[i] = EM_L1_L2(n_i,norm_grids[[1]],pdf_norm,range_norm)$L2
  
  # Normal Distribution - GW2
  
  n_i <- hist(x, breaks = norm_grids[[2]], plot = F)$counts
  
  res_norm_pdf_gw_21_method1[i] = L1_Distance_method1(n_i,norm_grids[[2]],pdf_norm,range_norm)
  res_norm_pdf_gw_22_method1[i] = L2_Distance_method1(n_i,norm_grids[[2]],pdf_norm,range_norm)
  
  res_norm_pdf_gw_21_method2[i] = L1_Distance_method2(n_i,norm_grids[[2]],pdf_norm,range_norm)
  res_norm_pdf_gw_22_method2[i] = L2_Distance_method2(n_i,norm_grids[[2]],pdf_norm,range_norm)
  
  res_norm_pdf_gw_21_method3[i] = L1_Distance_method3(n_i,norm_grids[[2]],pdf_norm,range_norm)
  res_norm_pdf_gw_22_method3[i] = L2_Distance_method3(n_i,norm_grids[[2]],pdf_norm,range_norm)
  
  res_norm_pdf_gw_21_EM[i] = EM_L1_L2(n_i,norm_grids[[2]],pdf_norm,range_norm)$L1
  res_norm_pdf_gw_22_EM[i] = EM_L1_L2(n_i,norm_grids[[2]],pdf_norm,range_norm)$L2
  
  # Normal Distribution - GW3
  
  n_i <- hist(x, breaks = norm_grids[[3]], plot = F)$counts
  
  res_norm_pdf_gw_31_method1[i] = L1_Distance_method1(n_i,norm_grids[[3]],pdf_norm,range_norm)
  res_norm_pdf_gw_32_method1[i] = L2_Distance_method1(n_i,norm_grids[[3]],pdf_norm,range_norm)
  
  res_norm_pdf_gw_31_method2[i] = L1_Distance_method2(n_i,norm_grids[[3]],pdf_norm,range_norm)
  res_norm_pdf_gw_32_method2[i] = L2_Distance_method2(n_i,norm_grids[[3]],pdf_norm,range_norm)
  
  res_norm_pdf_gw_31_method3[i] = L1_Distance_method3(n_i,norm_grids[[3]],pdf_norm,range_norm)
  res_norm_pdf_gw_32_method3[i] = L2_Distance_method3(n_i,norm_grids[[3]],pdf_norm,range_norm)
  
  res_norm_pdf_gw_31_EM[i] = EM_L1_L2(n_i,norm_grids[[3]],pdf_norm,range_norm)$L1
  res_norm_pdf_gw_32_EM[i] = EM_L1_L2(n_i,norm_grids[[3]],pdf_norm,range_norm)$L2
  
  # Normal Distribution - GW4
  
  n_i <- hist(x, breaks = norm_grids[[4]], plot = F)$counts
  
  res_norm_pdf_gw_41_method1[i] = L1_Distance_method1(n_i,norm_grids[[4]],pdf_norm,range_norm)
  res_norm_pdf_gw_42_method1[i] = L2_Distance_method1(n_i,norm_grids[[4]],pdf_norm,range_norm)
  
  res_norm_pdf_gw_41_method2[i] = L1_Distance_method2(n_i,norm_grids[[4]],pdf_norm,range_norm)
  res_norm_pdf_gw_42_method2[i] = L2_Distance_method2(n_i,norm_grids[[4]],pdf_norm,range_norm)
  
  res_norm_pdf_gw_41_method3[i] = L1_Distance_method3(n_i,norm_grids[[4]],pdf_norm,range_norm)
  res_norm_pdf_gw_42_method3[i] = L2_Distance_method3(n_i,norm_grids[[4]],pdf_norm,range_norm)
  
  res_norm_pdf_gw_41_EM[i] = EM_L1_L2(n_i,norm_grids[[4]],pdf_norm,range_norm)$L1
  res_norm_pdf_gw_42_EM[i] = EM_L1_L2(n_i,norm_grids[[4]],pdf_norm,range_norm)$L2
  
  # Normal Distribution - GW5
  
  n_i <- hist(x, breaks = norm_grids[[5]], plot = F)$counts
  
  res_norm_pdf_gw_51_method1[i] = L1_Distance_method1(n_i,norm_grids[[5]],pdf_norm,range_norm)
  res_norm_pdf_gw_52_method1[i] = L2_Distance_method1(n_i,norm_grids[[5]],pdf_norm,range_norm)
  
  res_norm_pdf_gw_51_method2[i] = L1_Distance_method2(n_i,norm_grids[[5]],pdf_norm,range_norm)
  res_norm_pdf_gw_52_method2[i] = L2_Distance_method2(n_i,norm_grids[[5]],pdf_norm,range_norm)
  
  res_norm_pdf_gw_51_method3[i] = L1_Distance_method3(n_i,norm_grids[[5]],pdf_norm,range_norm)
  res_norm_pdf_gw_52_method3[i] = L2_Distance_method3(n_i,norm_grids[[5]],pdf_norm,range_norm)
  
  res_norm_pdf_gw_51_EM[i] = EM_L1_L2(n_i,norm_grids[[5]],pdf_norm,range_norm)$L1
  res_norm_pdf_gw_52_EM[i] = EM_L1_L2(n_i,norm_grids[[5]],pdf_norm,range_norm)$L2
  
}

# Method-1 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L1-GW/L1-PDF-1-Normal-gw-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_norm_pdf_gw_11_method1,res_norm_pdf_gw_21_method1,res_norm_pdf_gw_31_method1,res_norm_pdf_gw_41_method1,res_norm_pdf_gw_51_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,0.7))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L2-GW/L2-PDF-1-Normal-gw-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_norm_pdf_gw_12_method1,res_norm_pdf_gw_22_method1,res_norm_pdf_gw_32_method1,res_norm_pdf_gw_42_method1,res_norm_pdf_gw_52_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,0.4))
dev.off()

# Method-2 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L1-GW/L1-PDF-1-Normal-gw-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_norm_pdf_gw_11_method2,res_norm_pdf_gw_21_method2,res_norm_pdf_gw_31_method2,res_norm_pdf_gw_41_method2,res_norm_pdf_gw_51_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,0.7))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L2-GW/L2-PDF-1-Normal-gw-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_norm_pdf_gw_12_method2,res_norm_pdf_gw_22_method2,res_norm_pdf_gw_32_method2,res_norm_pdf_gw_42_method2,res_norm_pdf_gw_52_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,0.4))
dev.off()

# Method-3 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L1-GW/L1-PDF-1-Normal-gw-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_norm_pdf_gw_11_method3,res_norm_pdf_gw_21_method3,res_norm_pdf_gw_31_method3,res_norm_pdf_gw_41_method3,res_norm_pdf_gw_51_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,0.7))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L2-GW/L2-PDF-1-Normal-gw-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_norm_pdf_gw_12_method3,res_norm_pdf_gw_22_method3,res_norm_pdf_gw_32_method3,res_norm_pdf_gw_42_method3,res_norm_pdf_gw_52_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,0.4))
dev.off()

# EM # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L1-GW/L1-PDF-1-Normal-gw-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_norm_pdf_gw_11_EM,res_norm_pdf_gw_21_EM,res_norm_pdf_gw_31_EM,res_norm_pdf_gw_41_EM,res_norm_pdf_gw_51_EM,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,0.7))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L2-GW/L2-PDF-1-Normal-gw-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_norm_pdf_gw_12_EM,res_norm_pdf_gw_22_EM,res_norm_pdf_gw_32_EM,res_norm_pdf_gw_42_EM,res_norm_pdf_gw_52_EM,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors1,ylim=c(0,0.4))
dev.off()


# Beta Distribution - Grid Width Simulation 

beta_grids = list(seq(0,1.2, by = beta_deltas[5]),
                  seq(0,1.3, by = beta_deltas[4]),
                  seq(0,1.2, by = beta_deltas[3]),
                  seq(0,1.2, by = beta_deltas[2]),
                  seq(0,1.1, by = beta_deltas[1]))

# Method-1 #

res_beta_pdf_gw_11_method1 = rep(0,B)
res_beta_pdf_gw_21_method1 = rep(0,B)
res_beta_pdf_gw_31_method1 = rep(0,B)
res_beta_pdf_gw_41_method1 = rep(0,B)
res_beta_pdf_gw_51_method1 = rep(0,B)

res_beta_pdf_gw_12_method1 = rep(0,B)
res_beta_pdf_gw_22_method1 = rep(0,B)
res_beta_pdf_gw_32_method1 = rep(0,B)
res_beta_pdf_gw_42_method1 = rep(0,B)
res_beta_pdf_gw_52_method1 = rep(0,B)

# Method-2 #

res_beta_pdf_gw_11_method2 = rep(0,B)
res_beta_pdf_gw_21_method2 = rep(0,B)
res_beta_pdf_gw_31_method2 = rep(0,B)
res_beta_pdf_gw_41_method2 = rep(0,B)
res_beta_pdf_gw_51_method2 = rep(0,B)

res_beta_pdf_gw_12_method2 = rep(0,B)
res_beta_pdf_gw_22_method2 = rep(0,B)
res_beta_pdf_gw_32_method2 = rep(0,B)
res_beta_pdf_gw_42_method2 = rep(0,B)
res_beta_pdf_gw_52_method2 = rep(0,B)

# Method-3 #

res_beta_pdf_gw_11_method3 = rep(0,B)
res_beta_pdf_gw_21_method3 = rep(0,B)
res_beta_pdf_gw_31_method3 = rep(0,B)
res_beta_pdf_gw_41_method3 = rep(0,B)
res_beta_pdf_gw_51_method3 = rep(0,B)

res_beta_pdf_gw_12_method3 = rep(0,B)
res_beta_pdf_gw_22_method3 = rep(0,B)
res_beta_pdf_gw_32_method3 = rep(0,B)
res_beta_pdf_gw_42_method3 = rep(0,B)
res_beta_pdf_gw_52_method3 = rep(0,B)

# EM #

res_beta_pdf_gw_11_EM = rep(0,B)
res_beta_pdf_gw_21_EM = rep(0,B)
res_beta_pdf_gw_31_EM = rep(0,B)
res_beta_pdf_gw_41_EM = rep(0,B)
res_beta_pdf_gw_51_EM = rep(0,B)

res_beta_pdf_gw_12_EM = rep(0,B)
res_beta_pdf_gw_22_EM = rep(0,B)
res_beta_pdf_gw_32_EM = rep(0,B)
res_beta_pdf_gw_42_EM = rep(0,B)
res_beta_pdf_gw_52_EM = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  x <- sort(rbeta(n,alpha,beta))
  
  # Beta Distribution - GW1
  
  n_i <- hist(x, breaks = beta_grids[[1]], plot = F)$counts
  
  res_beta_pdf_gw_11_method1[i] = L1_Distance_method1(n_i,beta_grids[[1]],pdf_beta,range_beta)
  res_beta_pdf_gw_12_method1[i] = L2_Distance_method1(n_i,beta_grids[[1]],pdf_beta,range_beta)
  
  res_beta_pdf_gw_11_method2[i] = L1_Distance_method2(n_i,beta_grids[[1]],pdf_beta,range_beta)
  res_beta_pdf_gw_12_method2[i] = L2_Distance_method2(n_i,beta_grids[[1]],pdf_beta,range_beta)
  
  res_beta_pdf_gw_11_method3[i] = L1_Distance_method3(n_i,beta_grids[[1]],pdf_beta,range_beta)
  res_beta_pdf_gw_12_method3[i] = L2_Distance_method3(n_i,beta_grids[[1]],pdf_beta,range_beta)
  
  res_beta_pdf_gw_11_EM[i] = EM_L1_L2(n_i,beta_grids[[1]],pdf_beta,range_beta)$L1
  res_beta_pdf_gw_12_EM[i] = EM_L1_L2(n_i,beta_grids[[1]],pdf_beta,range_beta)$L2
  
  # Beta Distribution - GW2

  n_i <- hist(x, breaks = beta_grids[[2]], plot = F)$counts
  
  res_beta_pdf_gw_21_method1[i] = L1_Distance_method1(n_i,beta_grids[[2]],pdf_beta,range_beta)
  res_beta_pdf_gw_22_method1[i] = L2_Distance_method1(n_i,beta_grids[[2]],pdf_beta,range_beta)
  
  res_beta_pdf_gw_21_method2[i] = L1_Distance_method2(n_i,beta_grids[[2]],pdf_beta,range_beta)
  res_beta_pdf_gw_22_method2[i] = L2_Distance_method2(n_i,beta_grids[[2]],pdf_beta,range_beta)
  
  res_beta_pdf_gw_21_method3[i] = L1_Distance_method3(n_i,beta_grids[[2]],pdf_beta,range_beta)
  res_beta_pdf_gw_22_method3[i] = L2_Distance_method3(n_i,beta_grids[[2]],pdf_beta,range_beta)
  
  res_beta_pdf_gw_21_EM[i] = EM_L1_L2(n_i,beta_grids[[2]],pdf_beta,range_beta)$L1
  res_beta_pdf_gw_22_EM[i] = EM_L1_L2(n_i,beta_grids[[2]],pdf_beta,range_beta)$L2
  
  # Beta Distribution - GW3

  n_i <- hist(x, breaks = beta_grids[[3]], plot = F)$counts
  
  res_beta_pdf_gw_31_method1[i] = L1_Distance_method1(n_i,beta_grids[[3]],pdf_beta,range_beta)
  res_beta_pdf_gw_32_method1[i] = L2_Distance_method1(n_i,beta_grids[[3]],pdf_beta,range_beta)
  
  res_beta_pdf_gw_31_method2[i] = L1_Distance_method2(n_i,beta_grids[[3]],pdf_beta,range_beta)
  res_beta_pdf_gw_32_method2[i] = L2_Distance_method2(n_i,beta_grids[[3]],pdf_beta,range_beta)
  
  res_beta_pdf_gw_31_method3[i] = L1_Distance_method3(n_i,beta_grids[[3]],pdf_beta,range_beta)
  res_beta_pdf_gw_32_method3[i] = L2_Distance_method3(n_i,beta_grids[[3]],pdf_beta,range_beta)
  
  res_beta_pdf_gw_31_EM[i] = EM_L1_L2(n_i,beta_grids[[3]],pdf_beta,range_beta)$L1
  res_beta_pdf_gw_32_EM[i] = EM_L1_L2(n_i,beta_grids[[3]],pdf_beta,range_beta)$L2
  
  # Beta Distribution - GW4

  n_i <- hist(x, breaks = beta_grids[[4]], plot = F)$counts
  
  res_beta_pdf_gw_41_method1[i] = L1_Distance_method1(n_i,beta_grids[[4]],pdf_beta,range_beta)
  res_beta_pdf_gw_42_method1[i] = L2_Distance_method1(n_i,beta_grids[[4]],pdf_beta,range_beta)
  
  res_beta_pdf_gw_41_method2[i] = L1_Distance_method2(n_i,beta_grids[[4]],pdf_beta,range_beta)
  res_beta_pdf_gw_42_method2[i] = L2_Distance_method2(n_i,beta_grids[[4]],pdf_beta,range_beta)
  
  res_beta_pdf_gw_41_method3[i] = L1_Distance_method3(n_i,beta_grids[[4]],pdf_beta,range_beta)
  res_beta_pdf_gw_42_method3[i] = L2_Distance_method3(n_i,beta_grids[[4]],pdf_beta,range_beta)
  
  res_beta_pdf_gw_41_EM[i] = EM_L1_L2(n_i,beta_grids[[4]],pdf_beta,range_beta)$L1
  res_beta_pdf_gw_42_EM[i] = EM_L1_L2(n_i,beta_grids[[4]],pdf_beta,range_beta)$L2
  
  # Beta Distribution - GW5

  n_i <- hist(x, breaks = beta_grids[[5]], plot = F)$counts
  
  res_beta_pdf_gw_51_method1[i] = L1_Distance_method1(n_i,beta_grids[[5]],pdf_beta,range_beta)
  res_beta_pdf_gw_52_method1[i] = L2_Distance_method1(n_i,beta_grids[[5]],pdf_beta,range_beta)
  
  res_beta_pdf_gw_51_method2[i] = L1_Distance_method2(n_i,beta_grids[[5]],pdf_beta,range_beta)
  res_beta_pdf_gw_52_method2[i] = L2_Distance_method2(n_i,beta_grids[[5]],pdf_beta,range_beta)
  
  res_beta_pdf_gw_51_method3[i] = L1_Distance_method3(n_i,beta_grids[[5]],pdf_beta,range_beta)
  res_beta_pdf_gw_52_method3[i] = L2_Distance_method3(n_i,beta_grids[[5]],pdf_beta,range_beta)
  
  res_beta_pdf_gw_51_EM[i] = EM_L1_L2(n_i,beta_grids[[5]],pdf_beta,range_beta)$L1
  res_beta_pdf_gw_52_EM[i] = EM_L1_L2(n_i,beta_grids[[5]],pdf_beta,range_beta)$L2
  
}

# Method-1 #

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L1-GW/L1-PDF-1-Beta-gw-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pdf_gw_11_method1,res_beta_pdf_gw_21_method1,res_beta_pdf_gw_31_method1,res_beta_pdf_gw_41_method1,res_beta_pdf_gw_51_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.7))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L2-GW/L2-PDF-1-Beta-gw-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pdf_gw_12_method1,res_beta_pdf_gw_22_method1,res_beta_pdf_gw_32_method1,res_beta_pdf_gw_42_method1,res_beta_pdf_gw_52_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.7))
dev.off()

# Method-2 #

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L1-GW/L1-PDF-1-Beta-gw-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pdf_gw_11_method2,res_beta_pdf_gw_21_method2,res_beta_pdf_gw_31_method2,res_beta_pdf_gw_41_method2,res_beta_pdf_gw_51_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.7))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L2-GW/L2-PDF-1-Beta-gw-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pdf_gw_12_method2,res_beta_pdf_gw_22_method2,res_beta_pdf_gw_32_method2,res_beta_pdf_gw_42_method2,res_beta_pdf_gw_52_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.7))
dev.off()


# Method-3 #

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L1-GW/L1-PDF-1-Beta-gw-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pdf_gw_11_method3,res_beta_pdf_gw_21_method3,res_beta_pdf_gw_31_method3,res_beta_pdf_gw_41_method3,res_beta_pdf_gw_51_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.7))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L2-GW/L2-PDF-1-Beta-gw-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pdf_gw_12_method3,res_beta_pdf_gw_22_method3,res_beta_pdf_gw_32_method3,res_beta_pdf_gw_42_method3,res_beta_pdf_gw_52_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.7))
dev.off()

# EM #

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L1-GW/L1-PDF-1-Beta-gw-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pdf_gw_11_EM,res_beta_pdf_gw_21_EM,res_beta_pdf_gw_31_EM,res_beta_pdf_gw_41_EM,res_beta_pdf_gw_51_EM,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.7))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L2-GW/L2-PDF-1-Beta-gw-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_beta_pdf_gw_12_EM,res_beta_pdf_gw_22_EM,res_beta_pdf_gw_32_EM,res_beta_pdf_gw_42_EM,res_beta_pdf_gw_52_EM,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors2,ylim=c(0,0.7))
dev.off()


# Gamma Distribution - Grid Width Simulation 

MM = qgamma(0.9999999,alpha,beta) + 3
gamma_grids = list(seq(0,MM, by = gamma_deltas[5]),
                   seq(0,MM, by = gamma_deltas[4]),
                   seq(0,MM, by = gamma_deltas[3]),
                   seq(0,MM, by = gamma_deltas[2]),
                   seq(0,MM, by = gamma_deltas[1]))

# Method-1 # 

res_gamma_pdf_gw_11_method1 = rep(0,B)
res_gamma_pdf_gw_21_method1 = rep(0,B)
res_gamma_pdf_gw_31_method1 = rep(0,B)
res_gamma_pdf_gw_41_method1 = rep(0,B)
res_gamma_pdf_gw_51_method1 = rep(0,B)

res_gamma_pdf_gw_12_method1 = rep(0,B)
res_gamma_pdf_gw_22_method1 = rep(0,B)
res_gamma_pdf_gw_32_method1 = rep(0,B)
res_gamma_pdf_gw_42_method1 = rep(0,B)
res_gamma_pdf_gw_52_method1 = rep(0,B)

# Method-2 # 

res_gamma_pdf_gw_11_method2 = rep(0,B)
res_gamma_pdf_gw_21_method2 = rep(0,B)
res_gamma_pdf_gw_31_method2 = rep(0,B)
res_gamma_pdf_gw_41_method2 = rep(0,B)
res_gamma_pdf_gw_51_method2 = rep(0,B)

res_gamma_pdf_gw_12_method2 = rep(0,B)
res_gamma_pdf_gw_22_method2 = rep(0,B)
res_gamma_pdf_gw_32_method2 = rep(0,B)
res_gamma_pdf_gw_42_method2 = rep(0,B)
res_gamma_pdf_gw_52_method2 = rep(0,B)

# Method-3 # 

res_gamma_pdf_gw_11_method3 = rep(0,B)
res_gamma_pdf_gw_21_method3 = rep(0,B)
res_gamma_pdf_gw_31_method3 = rep(0,B)
res_gamma_pdf_gw_41_method3 = rep(0,B)
res_gamma_pdf_gw_51_method3 = rep(0,B)

res_gamma_pdf_gw_12_method3 = rep(0,B)
res_gamma_pdf_gw_22_method3 = rep(0,B)
res_gamma_pdf_gw_32_method3 = rep(0,B)
res_gamma_pdf_gw_42_method3 = rep(0,B)
res_gamma_pdf_gw_52_method3 = rep(0,B)

# EM # 

res_gamma_pdf_gw_11_EM = rep(0,B)
res_gamma_pdf_gw_21_EM = rep(0,B)
res_gamma_pdf_gw_31_EM = rep(0,B)
res_gamma_pdf_gw_41_EM = rep(0,B)
res_gamma_pdf_gw_51_EM = rep(0,B)

res_gamma_pdf_gw_12_EM = rep(0,B)
res_gamma_pdf_gw_22_EM = rep(0,B)
res_gamma_pdf_gw_32_EM = rep(0,B)
res_gamma_pdf_gw_42_EM = rep(0,B)
res_gamma_pdf_gw_52_EM = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  x <- sort(rgamma(n,alpha,beta))
  
  # Gamma Distribution - GW1
  
  n_i <- hist(x, breaks = gamma_grids[[1]], plot = F)$counts
  
  res_gamma_pdf_gw_11_method1[i] = L1_Distance_method1(n_i,gamma_grids[[1]],pdf_gamma,range_gamma)
  res_gamma_pdf_gw_12_method1[i] = L2_Distance_method1(n_i,gamma_grids[[1]],pdf_gamma,range_gamma)
  
  res_gamma_pdf_gw_11_method2[i] = L1_Distance_method2(n_i,gamma_grids[[1]],pdf_gamma,range_gamma)
  res_gamma_pdf_gw_12_method2[i] = L2_Distance_method2(n_i,gamma_grids[[1]],pdf_gamma,range_gamma)
  
  res_gamma_pdf_gw_11_method3[i] = L1_Distance_method3(n_i,gamma_grids[[1]],pdf_gamma,range_gamma)
  res_gamma_pdf_gw_12_method3[i] = L2_Distance_method3(n_i,gamma_grids[[1]],pdf_gamma,range_gamma)
  
  res_gamma_pdf_gw_11_EM[i] = EM_L1_L2(n_i,gamma_grids[[1]],pdf_gamma,range_gamma)$L1
  res_gamma_pdf_gw_12_EM[i] = EM_L1_L2(n_i,gamma_grids[[1]],pdf_gamma,range_gamma)$L2
  
  # Gamma Distribution - GW2
  
  n_i <- hist(x, breaks = gamma_grids[[2]], plot = F)$counts
  
  res_gamma_pdf_gw_21_method1[i] = L1_Distance_method1(n_i,gamma_grids[[2]],pdf_gamma,range_gamma)
  res_gamma_pdf_gw_22_method1[i] = L2_Distance_method1(n_i,gamma_grids[[2]],pdf_gamma,range_gamma)
  
  res_gamma_pdf_gw_21_method2[i] = L1_Distance_method2(n_i,gamma_grids[[2]],pdf_gamma,range_gamma)
  res_gamma_pdf_gw_22_method2[i] = L2_Distance_method2(n_i,gamma_grids[[2]],pdf_gamma,range_gamma)
  
  res_gamma_pdf_gw_21_method3[i] = L1_Distance_method3(n_i,gamma_grids[[2]],pdf_gamma,range_gamma)
  res_gamma_pdf_gw_22_method3[i] = L2_Distance_method3(n_i,gamma_grids[[2]],pdf_gamma,range_gamma)
  
  res_gamma_pdf_gw_21_EM[i] = EM_L1_L2(n_i,gamma_grids[[2]],pdf_gamma,range_gamma)$L1
  res_gamma_pdf_gw_22_EM[i] = EM_L1_L2(n_i,gamma_grids[[2]],pdf_gamma,range_gamma)$L2
  
  # Gamma Distribution - GW3
  
  n_i <- hist(x, breaks = gamma_grids[[3]], plot = F)$counts
  
  res_gamma_pdf_gw_31_method1[i] = L1_Distance_method1(n_i,gamma_grids[[3]],pdf_gamma,range_gamma)
  res_gamma_pdf_gw_32_method1[i] = L2_Distance_method1(n_i,gamma_grids[[3]],pdf_gamma,range_gamma)
  
  res_gamma_pdf_gw_31_method2[i] = L1_Distance_method2(n_i,gamma_grids[[3]],pdf_gamma,range_gamma)
  res_gamma_pdf_gw_32_method2[i] = L2_Distance_method2(n_i,gamma_grids[[3]],pdf_gamma,range_gamma)
  
  res_gamma_pdf_gw_31_method3[i] = L1_Distance_method3(n_i,gamma_grids[[3]],pdf_gamma,range_gamma)
  res_gamma_pdf_gw_32_method3[i] = L2_Distance_method3(n_i,gamma_grids[[3]],pdf_gamma,range_gamma)
  
  res_gamma_pdf_gw_31_EM[i] = EM_L1_L2(n_i,gamma_grids[[3]],pdf_gamma,range_gamma)$L1
  res_gamma_pdf_gw_32_EM[i] = EM_L1_L2(n_i,gamma_grids[[3]],pdf_gamma,range_gamma)$L2
  
  # Gamma Distribution - GW4
  
  n_i <- hist(x, breaks = gamma_grids[[4]], plot = F)$counts
  
  res_gamma_pdf_gw_41_method1[i] = L1_Distance_method1(n_i,gamma_grids[[4]],pdf_gamma,range_gamma)
  res_gamma_pdf_gw_42_method1[i] = L2_Distance_method1(n_i,gamma_grids[[4]],pdf_gamma,range_gamma)
  
  res_gamma_pdf_gw_41_method2[i] = L1_Distance_method2(n_i,gamma_grids[[4]],pdf_gamma,range_gamma)
  res_gamma_pdf_gw_42_method2[i] = L2_Distance_method2(n_i,gamma_grids[[4]],pdf_gamma,range_gamma)
  
  res_gamma_pdf_gw_41_method3[i] = L1_Distance_method3(n_i,gamma_grids[[4]],pdf_gamma,range_gamma)
  res_gamma_pdf_gw_42_method3[i] = L2_Distance_method3(n_i,gamma_grids[[4]],pdf_gamma,range_gamma)
  
  res_gamma_pdf_gw_41_EM[i] = EM_L1_L2(n_i,gamma_grids[[4]],pdf_gamma,range_gamma)$L1
  res_gamma_pdf_gw_42_EM[i] = EM_L1_L2(n_i,gamma_grids[[4]],pdf_gamma,range_gamma)$L2
  
  # Gamma Distribution - GW5
  
  n_i <- hist(x, breaks = gamma_grids[[5]], plot = F)$counts
  
  res_gamma_pdf_gw_51_method1[i] = L1_Distance_method1(n_i,gamma_grids[[5]],pdf_gamma,range_gamma)
  res_gamma_pdf_gw_52_method1[i] = L2_Distance_method1(n_i,gamma_grids[[5]],pdf_gamma,range_gamma)
  
  res_gamma_pdf_gw_51_method2[i] = L1_Distance_method2(n_i,gamma_grids[[5]],pdf_gamma,range_gamma)
  res_gamma_pdf_gw_52_method2[i] = L2_Distance_method2(n_i,gamma_grids[[5]],pdf_gamma,range_gamma)
  
  res_gamma_pdf_gw_51_method3[i] = L1_Distance_method3(n_i,gamma_grids[[5]],pdf_gamma,range_gamma)
  res_gamma_pdf_gw_52_method3[i] = L2_Distance_method3(n_i,gamma_grids[[5]],pdf_gamma,range_gamma)
  
  res_gamma_pdf_gw_51_EM[i] = EM_L1_L2(n_i,gamma_grids[[5]],pdf_gamma,range_gamma)$L1
  res_gamma_pdf_gw_52_EM[i] = EM_L1_L2(n_i,gamma_grids[[5]],pdf_gamma,range_gamma)$L2
  
}

# Method-1 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L1-GW/L1-PDF-1-Gamma-gw-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pdf_gw_11_method1,res_gamma_pdf_gw_21_method1,res_gamma_pdf_gw_31_method1,res_gamma_pdf_gw_41_method1,res_gamma_pdf_gw_51_method1,
        names = c("1 (7)", "0.8 (8)", "0.6 (11)", "0.4 (16)", "0.2 (31)"),
        col = colors3,ylim=c(0,0.9))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L2-GW/L2-PDF-1-Gamma-gw-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pdf_gw_12_method1,res_gamma_pdf_gw_22_method1,res_gamma_pdf_gw_32_method1,res_gamma_pdf_gw_42_method1,res_gamma_pdf_gw_52_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,0.8))
dev.off()

# Method-2 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L1-GW/L1-PDF-1-Gamma-gw-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pdf_gw_11_method2,res_gamma_pdf_gw_21_method2,res_gamma_pdf_gw_31_method2,res_gamma_pdf_gw_41_method2,res_gamma_pdf_gw_51_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,0.9))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L2-GW/L2-PDF-1-Gamma-gw-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pdf_gw_12_method2,res_gamma_pdf_gw_22_method2,res_gamma_pdf_gw_32_method2,res_gamma_pdf_gw_42_method2,res_gamma_pdf_gw_52_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,0.8))
dev.off()

# Method-3 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L1-GW/L1-PDF-1-Gamma-gw-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pdf_gw_11_method3,res_gamma_pdf_gw_21_method3,res_gamma_pdf_gw_31_method3,res_gamma_pdf_gw_41_method3,res_gamma_pdf_gw_51_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,0.9))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L2-GW/L2-PDF-1-Gamma-gw-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pdf_gw_12_method3,res_gamma_pdf_gw_22_method3,res_gamma_pdf_gw_32_method3,res_gamma_pdf_gw_42_method3,res_gamma_pdf_gw_52_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,0.8))
dev.off()

# EM # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L1-GW/L1-PDF-1-Gamma-gw-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pdf_gw_11_EM,res_gamma_pdf_gw_21_EM,res_gamma_pdf_gw_31_EM,res_gamma_pdf_gw_41_EM,res_gamma_pdf_gw_51_EM,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,0.9))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L2-GW/L2-PDF-1-Gamma-gw-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_gamma_pdf_gw_12_EM,res_gamma_pdf_gw_22_EM,res_gamma_pdf_gw_32_EM,res_gamma_pdf_gw_42_EM,res_gamma_pdf_gw_52_EM,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors3,ylim=c(0,0.8))
dev.off()

# Logistic Distribution - Grid Width Simulation 

MM = qlogis(0.9999999,mu,sigma) + 3
logistic_grids = list(seq(-MM,MM, by = logistic_deltas[5]),
                      seq(-MM,MM, by = logistic_deltas[4]),
                      seq(-MM,MM, by = logistic_deltas[3]),
                      seq(-MM,MM, by = logistic_deltas[2]),
                      seq(-MM,MM, by = logistic_deltas[1]))


# Method-1 # 

res_logis_pdf_gw_11_method1 = rep(0,B)
res_logis_pdf_gw_21_method1 = rep(0,B)
res_logis_pdf_gw_31_method1 = rep(0,B)
res_logis_pdf_gw_41_method1 = rep(0,B)
res_logis_pdf_gw_51_method1 = rep(0,B)

res_logis_pdf_gw_12_method1 = rep(0,B)
res_logis_pdf_gw_22_method1 = rep(0,B)
res_logis_pdf_gw_32_method1 = rep(0,B)
res_logis_pdf_gw_42_method1 = rep(0,B)
res_logis_pdf_gw_52_method1 = rep(0,B)

# Method-2 # 

res_logis_pdf_gw_11_method2 = rep(0,B)
res_logis_pdf_gw_21_method2 = rep(0,B)
res_logis_pdf_gw_31_method2 = rep(0,B)
res_logis_pdf_gw_41_method2 = rep(0,B)
res_logis_pdf_gw_51_method2 = rep(0,B)

res_logis_pdf_gw_12_method2 = rep(0,B)
res_logis_pdf_gw_22_method2 = rep(0,B)
res_logis_pdf_gw_32_method2 = rep(0,B)
res_logis_pdf_gw_42_method2 = rep(0,B)
res_logis_pdf_gw_52_method2 = rep(0,B)

# Method-3 # 

res_logis_pdf_gw_11_method3 = rep(0,B)
res_logis_pdf_gw_21_method3 = rep(0,B)
res_logis_pdf_gw_31_method3 = rep(0,B)
res_logis_pdf_gw_41_method3 = rep(0,B)
res_logis_pdf_gw_51_method3 = rep(0,B)

res_logis_pdf_gw_12_method3 = rep(0,B)
res_logis_pdf_gw_22_method3 = rep(0,B)
res_logis_pdf_gw_32_method3 = rep(0,B)
res_logis_pdf_gw_42_method3 = rep(0,B)
res_logis_pdf_gw_52_method3 = rep(0,B)

# EM # 

res_logis_pdf_gw_11_EM = rep(0,B)
res_logis_pdf_gw_21_EM = rep(0,B)
res_logis_pdf_gw_31_EM = rep(0,B)
res_logis_pdf_gw_41_EM = rep(0,B)
res_logis_pdf_gw_51_EM = rep(0,B)

res_logis_pdf_gw_12_EM = rep(0,B)
res_logis_pdf_gw_22_EM = rep(0,B)
res_logis_pdf_gw_32_EM = rep(0,B)
res_logis_pdf_gw_42_EM = rep(0,B)
res_logis_pdf_gw_52_EM = rep(0,B)

set.seed(9)

for (i in 1:B) {
  
  x = sort(rlogis(n))
  
  # Logistic Distribution - GW1
  
  n_i <- hist(x, breaks = logistic_grids[[1]], plot = F)$counts
  
  res_logis_pdf_gw_11_method1[i] = L1_Distance_method1(n_i,logistic_grids[[1]],pdf_logistic,range_logistic)
  res_logis_pdf_gw_12_method1[i] = L2_Distance_method1(n_i,logistic_grids[[1]],pdf_logistic,range_logistic)
  
  res_logis_pdf_gw_11_method2[i] = L1_Distance_method2(n_i,logistic_grids[[1]],pdf_logistic,range_logistic)
  res_logis_pdf_gw_12_method2[i] = L2_Distance_method2(n_i,logistic_grids[[1]],pdf_logistic,range_logistic)
  
  res_logis_pdf_gw_11_method3[i] = L1_Distance_method3(n_i,logistic_grids[[1]],pdf_logistic,range_logistic)
  res_logis_pdf_gw_12_method3[i] = L2_Distance_method3(n_i,logistic_grids[[1]],pdf_logistic,range_logistic)
  
  res_logis_pdf_gw_11_EM[i] = EM_L1_L2(n_i,logistic_grids[[1]],pdf_logistic,range_logistic)$L1
  res_logis_pdf_gw_12_EM[i] = EM_L1_L2(n_i,logistic_grids[[1]],pdf_logistic,range_logistic)$L2
  
  # Logistic Distribution - GW2
  
  n_i <- hist(x, breaks = logistic_grids[[2]], plot = F)$counts
  
  res_logis_pdf_gw_21_method1[i] = L1_Distance_method1(n_i,logistic_grids[[2]],pdf_logistic,range_logistic)
  res_logis_pdf_gw_22_method1[i] = L2_Distance_method1(n_i,logistic_grids[[2]],pdf_logistic,range_logistic)
  
  res_logis_pdf_gw_21_method2[i] = L1_Distance_method2(n_i,logistic_grids[[2]],pdf_logistic,range_logistic)
  res_logis_pdf_gw_22_method2[i] = L2_Distance_method2(n_i,logistic_grids[[2]],pdf_logistic,range_logistic)
  
  res_logis_pdf_gw_21_method3[i] = L1_Distance_method3(n_i,logistic_grids[[2]],pdf_logistic,range_logistic)
  res_logis_pdf_gw_22_method3[i] = L2_Distance_method3(n_i,logistic_grids[[2]],pdf_logistic,range_logistic)
  
  res_logis_pdf_gw_21_EM[i] = EM_L1_L2(n_i,logistic_grids[[2]],pdf_logistic,range_logistic)$L1
  res_logis_pdf_gw_22_EM[i] = EM_L1_L2(n_i,logistic_grids[[2]],pdf_logistic,range_logistic)$L2
  
  # Logistic Distribution - GW3
  
  n_i <- hist(x, breaks = logistic_grids[[3]], plot = F)$counts
  
  res_logis_pdf_gw_31_method1[i] = L1_Distance_method1(n_i,logistic_grids[[3]],pdf_logistic,range_logistic)
  res_logis_pdf_gw_32_method1[i] = L2_Distance_method1(n_i,logistic_grids[[3]],pdf_logistic,range_logistic)
  
  res_logis_pdf_gw_31_method2[i] = L1_Distance_method2(n_i,logistic_grids[[3]],pdf_logistic,range_logistic)
  res_logis_pdf_gw_32_method2[i] = L2_Distance_method2(n_i,logistic_grids[[3]],pdf_logistic,range_logistic)
  
  res_logis_pdf_gw_31_method3[i] = L1_Distance_method3(n_i,logistic_grids[[3]],pdf_logistic,range_logistic)
  res_logis_pdf_gw_32_method3[i] = L2_Distance_method3(n_i,logistic_grids[[3]],pdf_logistic,range_logistic)  
  
  res_logis_pdf_gw_31_EM[i] = EM_L1_L2(n_i,logistic_grids[[3]],pdf_logistic,range_logistic)$L1
  res_logis_pdf_gw_32_EM[i] = EM_L1_L2(n_i,logistic_grids[[3]],pdf_logistic,range_logistic)$L2
  
  # Logistic Distribution - GW4
  
  n_i <- hist(x, breaks = logistic_grids[[4]], plot = F)$counts
  
  res_logis_pdf_gw_41_method1[i] = L1_Distance_method1(n_i,logistic_grids[[4]],pdf_logistic,range_logistic)
  res_logis_pdf_gw_42_method1[i] = L2_Distance_method1(n_i,logistic_grids[[4]],pdf_logistic,range_logistic)
  
  res_logis_pdf_gw_41_method2[i] = L1_Distance_method2(n_i,logistic_grids[[4]],pdf_logistic,range_logistic)
  res_logis_pdf_gw_42_method2[i] = L2_Distance_method2(n_i,logistic_grids[[4]],pdf_logistic,range_logistic)
  
  res_logis_pdf_gw_41_method3[i] = L1_Distance_method3(n_i,logistic_grids[[4]],pdf_logistic,range_logistic)
  res_logis_pdf_gw_42_method3[i] = L2_Distance_method3(n_i,logistic_grids[[4]],pdf_logistic,range_logistic)
  
  res_logis_pdf_gw_41_EM[i] = EM_L1_L2(n_i,logistic_grids[[4]],pdf_logistic,range_logistic)$L1
  res_logis_pdf_gw_42_EM[i] = EM_L1_L2(n_i,logistic_grids[[4]],pdf_logistic,range_logistic)$L2
  
  # Logistic Distribution - GW5
  
  n_i <- hist(x, breaks = logistic_grids[[5]], plot = F)$counts
  
  res_logis_pdf_gw_51_method1[i] = L1_Distance_method1(n_i,logistic_grids[[5]],pdf_logistic,range_logistic)
  res_logis_pdf_gw_52_method1[i] = L2_Distance_method1(n_i,logistic_grids[[5]],pdf_logistic,range_logistic)

  res_logis_pdf_gw_51_method2[i] = L1_Distance_method2(n_i,logistic_grids[[5]],pdf_logistic,range_logistic)
  res_logis_pdf_gw_52_method2[i] = L2_Distance_method2(n_i,logistic_grids[[5]],pdf_logistic,range_logistic)
  
  res_logis_pdf_gw_51_method3[i] = L1_Distance_method3(n_i,logistic_grids[[5]],pdf_logistic,range_logistic)
  res_logis_pdf_gw_52_method3[i] = L2_Distance_method3(n_i,logistic_grids[[5]],pdf_logistic,range_logistic)
  
  res_logis_pdf_gw_51_EM[i] = EM_L1_L2(n_i,logistic_grids[[5]],pdf_logistic,range_logistic)$L1
  res_logis_pdf_gw_52_EM[i] = EM_L1_L2(n_i,logistic_grids[[5]],pdf_logistic,range_logistic)$L2
  
}

# Method-1 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L1-GW/L1-PDF-1-Logistic-gw-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pdf_gw_11_method1,res_logis_pdf_gw_21_method1,res_logis_pdf_gw_31_method1,res_logis_pdf_gw_41_method1,res_logis_pdf_gw_51_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,1.2))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L2-GW/L2-PDF-1-Logistic-gw-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pdf_gw_12_method1,res_logis_pdf_gw_22_method1,res_logis_pdf_gw_32_method1,res_logis_pdf_gw_42_method1,res_logis_pdf_gw_52_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,0.4))
dev.off()

# Method-2 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L1-GW/L1-PDF-1-Logistic-gw-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pdf_gw_11_method2,res_logis_pdf_gw_21_method2,res_logis_pdf_gw_31_method2,res_logis_pdf_gw_41_method2,res_logis_pdf_gw_51_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,1.2))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L2-GW/L2-PDF-1-Logistic-gw-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pdf_gw_12_method2,res_logis_pdf_gw_22_method2,res_logis_pdf_gw_32_method2,res_logis_pdf_gw_42_method2,res_logis_pdf_gw_52_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,0.4))
dev.off()

# Method-3 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L1-GW/L1-PDF-1-Logistic-gw-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pdf_gw_11_method3,res_logis_pdf_gw_21_method3,res_logis_pdf_gw_31_method3,res_logis_pdf_gw_41_method3,res_logis_pdf_gw_51_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,1.2))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L2-GW/L2-PDF-1-Logistic-gw-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pdf_gw_12_method3,res_logis_pdf_gw_22_method3,res_logis_pdf_gw_32_method3,res_logis_pdf_gw_42_method3,res_logis_pdf_gw_52_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,0.4))
dev.off()

# EM # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L1-GW/L1-PDF-1-Logistic-gw-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pdf_gw_11_EM,res_logis_pdf_gw_21_EM,res_logis_pdf_gw_31_EM,res_logis_pdf_gw_41_EM,res_logis_pdf_gw_51_EM,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,1.2))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L2-GW/L2-PDF-1-Logistic-gw-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_logis_pdf_gw_12_EM,res_logis_pdf_gw_22_EM,res_logis_pdf_gw_32_EM,res_logis_pdf_gw_42_EM,res_logis_pdf_gw_52_EM,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors4,ylim=c(0,0.4))
dev.off()

# Student's t Distribution - Grid Width Simulation 

MM = qt(0.9999999,df) + 3
t_grids = list(seq(-MM,MM, by = t_deltas[5]),
               seq(-MM,MM, by = t_deltas[4]),
               seq(-MM,MM, by = t_deltas[3]),
               seq(-MM,MM, by = t_deltas[2]),
               seq(-MM,MM, by = t_deltas[1]))

# Method-1 # 

res_t_pdf_gw_11_method1 = rep(0,B)
res_t_pdf_gw_21_method1 = rep(0,B)
res_t_pdf_gw_31_method1 = rep(0,B)
res_t_pdf_gw_41_method1 = rep(0,B)
res_t_pdf_gw_51_method1 = rep(0,B)

res_t_pdf_gw_12_method1 = rep(0,B)
res_t_pdf_gw_22_method1 = rep(0,B)
res_t_pdf_gw_32_method1 = rep(0,B)
res_t_pdf_gw_42_method1 = rep(0,B)
res_t_pdf_gw_52_method1 = rep(0,B)

# Method-2 # 

res_t_pdf_gw_11_method2 = rep(0,B)
res_t_pdf_gw_21_method2 = rep(0,B)
res_t_pdf_gw_31_method2 = rep(0,B)
res_t_pdf_gw_41_method2 = rep(0,B)
res_t_pdf_gw_51_method2 = rep(0,B)

res_t_pdf_gw_12_method2 = rep(0,B)
res_t_pdf_gw_22_method2 = rep(0,B)
res_t_pdf_gw_32_method2 = rep(0,B)
res_t_pdf_gw_42_method2 = rep(0,B)
res_t_pdf_gw_52_method2 = rep(0,B)

# Method-3 # 

res_t_pdf_gw_11_method3 = rep(0,B)
res_t_pdf_gw_21_method3 = rep(0,B)
res_t_pdf_gw_31_method3 = rep(0,B)
res_t_pdf_gw_41_method3 = rep(0,B)
res_t_pdf_gw_51_method3 = rep(0,B)

res_t_pdf_gw_12_method3 = rep(0,B)
res_t_pdf_gw_22_method3 = rep(0,B)
res_t_pdf_gw_32_method3 = rep(0,B)
res_t_pdf_gw_42_method3 = rep(0,B)
res_t_pdf_gw_52_method3 = rep(0,B)

# EM # 

res_t_pdf_gw_11_EM = rep(0,B)
res_t_pdf_gw_21_EM = rep(0,B)
res_t_pdf_gw_31_EM = rep(0,B)
res_t_pdf_gw_41_EM = rep(0,B)
res_t_pdf_gw_51_EM = rep(0,B)

res_t_pdf_gw_12_EM = rep(0,B)
res_t_pdf_gw_22_EM = rep(0,B)
res_t_pdf_gw_32_EM = rep(0,B)
res_t_pdf_gw_42_EM = rep(0,B)
res_t_pdf_gw_52_EM = rep(0,B)

set.seed(10)

for (i in 1:B) {
  
  x = sort(rt(n,df = df))
  
  # t Distribution - GW1

  n_i <- hist(x, breaks = t_grids[[1]], plot = F)$counts
  
  res_t_pdf_gw_11_method1[i] = L1_Distance_method1(n_i,t_grids[[1]],pdf_t,range_t)
  res_t_pdf_gw_12_method1[i] = L2_Distance_method1(n_i,t_grids[[1]],pdf_t,range_t)
  
  res_t_pdf_gw_11_method2[i] = L1_Distance_method2(n_i,t_grids[[1]],pdf_t,range_t)
  res_t_pdf_gw_12_method2[i] = L2_Distance_method2(n_i,t_grids[[1]],pdf_t,range_t)
  
  res_t_pdf_gw_11_method3[i] = L1_Distance_method3(n_i,t_grids[[1]],pdf_t,range_t)
  res_t_pdf_gw_12_method3[i] = L2_Distance_method3(n_i,t_grids[[1]],pdf_t,range_t)
  
  res_t_pdf_gw_11_EM[i] = EM_L1_L2(n_i,t_grids[[1]],pdf_t,range_t)$L1
  res_t_pdf_gw_12_EM[i] = EM_L1_L2(n_i,t_grids[[1]],pdf_t,range_t)$L2
  
  # t Distribution - GW2
  
  n_i <- hist(x, breaks = t_grids[[2]], plot = F)$counts
  
  res_t_pdf_gw_21_method1[i] = L1_Distance_method1(n_i,t_grids[[2]],pdf_t,range_t)
  res_t_pdf_gw_22_method1[i] = L2_Distance_method1(n_i,t_grids[[2]],pdf_t,range_t)
  
  res_t_pdf_gw_21_method2[i] = L1_Distance_method2(n_i,t_grids[[2]],pdf_t,range_t)
  res_t_pdf_gw_22_method2[i] = L2_Distance_method2(n_i,t_grids[[2]],pdf_t,range_t)
  
  res_t_pdf_gw_21_method3[i] = L1_Distance_method3(n_i,t_grids[[2]],pdf_t,range_t)
  res_t_pdf_gw_22_method3[i] = L2_Distance_method3(n_i,t_grids[[2]],pdf_t,range_t)
  
  res_t_pdf_gw_21_EM[i] = EM_L1_L2(n_i,t_grids[[2]],pdf_t,range_t)$L1
  res_t_pdf_gw_22_EM[i] = EM_L1_L2(n_i,t_grids[[2]],pdf_t,range_t)$L2
  
  # t Distribution - GW3
  
  n_i <- hist(x, breaks = t_grids[[3]], plot = F)$counts
  
  res_t_pdf_gw_31_method1[i] = L1_Distance_method1(n_i,t_grids[[3]],pdf_t,range_t)
  res_t_pdf_gw_32_method1[i] = L2_Distance_method1(n_i,t_grids[[3]],pdf_t,range_t)
  
  res_t_pdf_gw_31_method2[i] = L1_Distance_method2(n_i,t_grids[[3]],pdf_t,range_t)
  res_t_pdf_gw_32_method2[i] = L2_Distance_method2(n_i,t_grids[[3]],pdf_t,range_t)
  
  res_t_pdf_gw_31_method3[i] = L1_Distance_method3(n_i,t_grids[[3]],pdf_t,range_t)
  res_t_pdf_gw_32_method3[i] = L2_Distance_method3(n_i,t_grids[[3]],pdf_t,range_t)
  
  res_t_pdf_gw_31_EM[i] = EM_L1_L2(n_i,t_grids[[3]],pdf_t,range_t)$L1
  res_t_pdf_gw_32_EM[i] = EM_L1_L2(n_i,t_grids[[3]],pdf_t,range_t)$L2
  
  # t Distribution - GW4
  
  n_i <- hist(x, breaks = t_grids[[4]], plot = F)$counts
  
  res_t_pdf_gw_41_method1[i] = L1_Distance_method1(n_i,t_grids[[4]],pdf_t,range_t)
  res_t_pdf_gw_42_method1[i] = L2_Distance_method1(n_i,t_grids[[4]],pdf_t,range_t)
  
  res_t_pdf_gw_41_method2[i] = L1_Distance_method2(n_i,t_grids[[4]],pdf_t,range_t)
  res_t_pdf_gw_42_method2[i] = L2_Distance_method2(n_i,t_grids[[4]],pdf_t,range_t)
  
  res_t_pdf_gw_41_method3[i] = L1_Distance_method3(n_i,t_grids[[4]],pdf_t,range_t)
  res_t_pdf_gw_42_method3[i] = L2_Distance_method3(n_i,t_grids[[4]],pdf_t,range_t)
  
  res_t_pdf_gw_41_EM[i] = EM_L1_L2(n_i,t_grids[[4]],pdf_t,range_t)$L1
  res_t_pdf_gw_42_EM[i] = EM_L1_L2(n_i,t_grids[[4]],pdf_t,range_t)$L2
  
  # t Distribution - GW5
  
  n_i <- hist(x, breaks = t_grids[[5]], plot = F)$counts
  
  res_t_pdf_gw_51_method1[i] = L1_Distance_method1(n_i,t_grids[[5]],pdf_t,range_t)
  res_t_pdf_gw_52_method1[i] = L2_Distance_method1(n_i,t_grids[[5]],pdf_t,range_t)
  
  res_t_pdf_gw_51_method2[i] = L1_Distance_method2(n_i,t_grids[[5]],pdf_t,range_t)
  res_t_pdf_gw_52_method2[i] = L2_Distance_method2(n_i,t_grids[[5]],pdf_t,range_t)
  
  res_t_pdf_gw_51_method3[i] = L1_Distance_method3(n_i,t_grids[[5]],pdf_t,range_t)
  res_t_pdf_gw_52_method3[i] = L2_Distance_method3(n_i,t_grids[[5]],pdf_t,range_t)
  
  res_t_pdf_gw_51_EM[i] = EM_L1_L2(n_i,t_grids[[5]],pdf_t,range_t)$L1
  res_t_pdf_gw_52_EM[i] = EM_L1_L2(n_i,t_grids[[5]],pdf_t,range_t)$L2
  
}

# Method-1 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L1-GW/L1-PDF-1-t-gw-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pdf_gw_11_method1,res_t_pdf_gw_21_method1,res_t_pdf_gw_31_method1,res_t_pdf_gw_41_method1,res_t_pdf_gw_51_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,1.5))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-1/L2-GW/L2-PDF-1-t-gw-Method1.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pdf_gw_12_method1,res_t_pdf_gw_22_method1,res_t_pdf_gw_32_method1,res_t_pdf_gw_42_method1,res_t_pdf_gw_52_method1,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,0.7))
dev.off()

# Method-2 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L1-GW/L1-PDF-1-t-gw-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pdf_gw_11_method2,res_t_pdf_gw_21_method2,res_t_pdf_gw_31_method2,res_t_pdf_gw_41_method2,res_t_pdf_gw_51_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,1.5))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-2/L2-GW/L2-PDF-1-t-gw-Method2.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pdf_gw_12_method2,res_t_pdf_gw_22_method2,res_t_pdf_gw_32_method2,res_t_pdf_gw_42_method2,res_t_pdf_gw_52_method2,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,0.7))
dev.off()


# Method-3 # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L1-GW/L1-PDF-1-t-gw-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pdf_gw_11_method3,res_t_pdf_gw_21_method3,res_t_pdf_gw_31_method3,res_t_pdf_gw_41_method3,res_t_pdf_gw_51_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,1.5))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method-3/L2-GW/L2-PDF-1-t-gw-Method3.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pdf_gw_12_method3,res_t_pdf_gw_22_method3,res_t_pdf_gw_32_method3,res_t_pdf_gw_42_method3,res_t_pdf_gw_52_method3,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,0.7))
dev.off()

# EM # 

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L1-GW/L1-PDF-1-t-gw-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pdf_gw_11_EM,res_t_pdf_gw_21_EM,res_t_pdf_gw_31_EM,res_t_pdf_gw_41_EM,res_t_pdf_gw_51_EM,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,1.5))
dev.off()

grDevices::pdf("Desktop/RA_Documents/RA-DOCUMENTS/PDF/EM/L2-GW/L2-PDF-1-t-gw-EM.pdf",width = 12, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)
boxplot(res_t_pdf_gw_12_EM,res_t_pdf_gw_22_EM,res_t_pdf_gw_32_EM,res_t_pdf_gw_42_EM,res_t_pdf_gw_52_EM,
        names = c("2.5","2.0","1.5","1.0","0.5"),
        col = colors5,ylim=c(0,0.7))
dev.off()

