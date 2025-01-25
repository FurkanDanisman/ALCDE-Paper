# Recovery Plot Comparisons 

# Normal 

f = function(x) dnorm(x,mean = mu,sd = sigma)

n = 1000

set.seed(9)

x <- sort(rnorm(n, mean = mu,sd=sigma))

MM = qnorm(0.9999999,mu,sigma) + 3
norm_grid <- seq(-MM,MM, by = norm_deltas[1])

n_i <- hist(x, breaks = norm_grid, plot = F)$counts

e = EM_Algorithm_Format(n_i,norm_grid)
em_r = em(e$bl, e$bu, e$Freq, e$Theta)
g = function(x) dnorm(x,mean = em_r$mu_estimate,sd = 1)
h = function(x) dnorm(x,mean = em_r$mu_estimate,sd = em_r$sigma_estimate)

grDevices::pdf("Recovery-Comparison-Normal.pdf",width = 10, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)

plot(x = 0,y = -1,
     ylim = c(0,0.4),col="white"
     ,xlab="x",ylab = "Probability",bty="l",
     xlim=c(-4,4))

# Method-1 # 
plot1 = Log_Con_Discr_Algorithm_method1(n_i,norm_grid)
lcm = gcmlcm(x = plot1$mle_x,  y = plot1$mle_psi,type = "lcm")
lines(x = plot1$x_knots,y = exp(lcm$y.knots) * plot1$Total_Coefficient,col="brown3", lwd = 3)

# Method-2 # 
plot2 = Log_Con_Discr_Algorithm_method2(n_i,norm_grid)
lcm = gcmlcm(x = plot2$mle_x,  y = plot2$mle_psi,type = "lcm")
lines(x = plot2$x_knots,y = lcm$y.knots * plot2$Total_Coefficient,col="purple3", lwd = 3)

# Method-3 # 
plot3 = Log_Con_Discr_Algorithm_method3(n_i,norm_grid)
x = plot3$x
f_3 = exp(plot3$phi)

lines(x,f_3,type="l",col="darkorange3",lwd = 3)
  
# EM # 
# plot(g,xlim = c(6,14),add = T, col = "lightblue", lwd = 3)
plot(h,xlim = c(-4,4),add = T, col = "navy", lwd = 3)

# True Distribution # 
plot(f,xlim=c(-4,4),add=T, col = "darkgreen",lwd = 3)

# legend("topright",
#       legend = c("ALCDE1","ALCDE2","ALCDE3","EM","True"), 
#       col = c("brown3","purple","darkorange3","navy","darkgreen"), 
#       lty = c(1,1,1,1,1),
#       lwd = c(3,3,3,3,3),
#       bty = "n", 
#       pt.cex = 2, 
#       cex = 3, 
#       text.col = "black", 
#       horiz = F, 
#       inset = c(0, 0),x.intersp = 0.2, y.intersp = 1,seg.len = 0.75) # -0.3

dev.off()

# Gamma Distribution 

f = function(x) dgamma(x,alpha,beta)

set.seed(9)

x <- sort(rgamma(n,alpha,beta))

MM = qgamma(0.9999999,alpha,beta) + 3
gamma_grid <- seq(0,MM, by = gamma_deltas[1])

n_i <- hist(x, breaks = gamma_grid, plot = F)$counts

e = EM_Algorithm_Format(n_i,gamma_grid)
em_r = em(e$bl, e$bu, e$Freq, theta_init = e$Theta)
g = function(x) dnorm(x,mean = em_r$mu_estimate,sd = 1)
h = function(x) dnorm(x,mean = em_r$mu_estimate,sd = em_r$sigma_estimate)

grDevices::pdf("Recovery-Comparison-Gamma.pdf",width = 10, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)

plot(x = 0,y = -1,
     ylim = c(0,2),col="white"
     ,xlab="x",ylab = "Probability",bty="l",
     xlim=c(-0.15,2))

# Method-1 # 
plot1 = Log_Con_Discr_Algorithm_method1(n_i,gamma_grid)
lcm = gcmlcm(x = plot1$mle_x,  y = plot1$mle_psi,type = "lcm")
lines(x = plot1$x_knots,y = exp(lcm$y.knots) * plot1$Total_Coefficient,col="brown3", lwd = 3)

# Method-2 # 
plot2 = Log_Con_Discr_Algorithm_method2(n_i,gamma_grid)
lcm = gcmlcm(x = plot2$mle_x,  y = plot2$mle_psi,type = "lcm")
lines(x = plot2$x_knots,y = lcm$y.knots * plot2$Total_Coefficient,col="purple3", lwd = 3)

# Method-3 # 
plot3 = Log_Con_Discr_Algorithm_method3(n_i,gamma_grid)
x = plot3$x
f_3 = exp(plot3$phi)

lines(x,f_3,type="l",col="darkorange3",lwd = 3)

# EM # 
# plot(g,xlim = c(6,14),add = T, col = "lightblue",lwd =3)
plot(h,xlim = c(0,2),add = T, col = "navy", lwd = 3)

# True Distribution # 
plot(f,xlim=c(0,2),add=T, col = "darkgreen", lwd = 3)

#legend("topright", 
#       legend = c("ALCDE1","ALCDE2","ALCDE3","EM","True"), 
#       col = c("brown3","purple","darkorange3","navy","darkgreen"), 
#       lty = c(1,1,1,1,1),
#       lwd = c(3,3,3,3,3),
#       bty = "n", 
#       pt.cex = 2, 
#       cex = 3, 
#       text.col = "black", 
#       horiz = F, 
#       inset = c(0, 0),x.intersp = 0.2, y.intersp = 1,seg.len = 0.75) # -0.3
dev.off()

# chisq Distribution 

f = function(x) dchisq(x,df)

set.seed(9)
x <- sort(rchisq(n,df))

MM = qchisq(0.9999999,df) + 3
chisq_grid <- seq(0,MM, by = chisqr_deltas[1])

n_i <- hist(x, breaks = chisq_grid, plot = F)$counts

e = EM_Algorithm_Format(n_i,chisq_grid)
em_r = em(e$bl, e$bu, e$Freq, e$Theta)
g = function(x) dnorm(x,mean = em_r$mu_estimate,sd = 1)
h = function(x) dnorm(x,mean = em_r$mu_estimate,sd = em_r$sigma_estimate)

grDevices::pdf("Recovery-Comparison-Chisqr.pdf",width = 10, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)

plot(x = 0,y = -1,
     ylim = c(0,0.3),col="white"
     ,xlab="x",ylab = "Probability",bty="l",
     xlim=c(-2.3,20))

# Method-1 # 
plot1 = Log_Con_Discr_Algorithm_method1(n_i,chisq_grid)
lcm = gcmlcm(x = plot1$mle_x,  y = plot1$mle_psi,type = "lcm")
lines(x = plot1$x_knots,y = exp(lcm$y.knots) * plot1$Total_Coefficient,col="brown3", lwd = 3)

# Method-2 # 
plot2 = Log_Con_Discr_Algorithm_method2(n_i,chisq_grid)
lcm = gcmlcm(x = plot2$mle_x,  y = plot2$mle_psi,type = "lcm")
lines(x = plot2$x_knots,y = lcm$y.knots * plot2$Total_Coefficient,col="purple3", lwd = 3)

# Method-3 # 
plot3 = Log_Con_Discr_Algorithm_method3(n_i,chisq_grid)
x = plot3$x
f_3 = exp(plot3$phi)

lines(x,f_3,type="l",col="darkorange3",lwd = 3)

# EM # 
# plot(g,xlim = c(6,14),add = T, col = "lightblue",lwd =3)
plot(h,xlim = c(0,20),add = T, col = "navy", lwd = 3)

# True Distribution # 
plot(f,xlim=c(0,20),add=T, col = "darkgreen", lwd = 3)

#legend("topright", 
#       legend = c("ALCDE1","ALCDE2","ACSCi","EM","True"), 
#       col = c("brown3","purple","darkorange3","navy","darkgreen"), 
#       lty = c(1,1,1,1,1),
#       lwd = c(3,3,3,3,3),
#       bty = "n", 
#       pt.cex = 2, 
#       cex = 3, 
#       text.col = "black", 
#       horiz = F, 
#       inset = c(0, 0),x.intersp = 0.2, y.intersp = 1,seg.len = 0.75) # -0.3
dev.off()


# t Distribution 

f = function(x) dt(x,df)

set.seed(9)
x <- sort(rt(n,df))

MM = qt(0.999999999999,df) + 30
t_grid <- seq(-MM,MM, by = t_deltas[1])

n_i <- hist(x, breaks = t_grid, plot = F)$counts

e = EM_Algorithm_Format(n_i,t_grid)
em_r = em(e$bl, e$bu, e$Freq, e$Theta)
g = function(x) dnorm(x,mean = em_r$mu_estimate,sd = 1)
h = function(x) dnorm(x,mean = em_r$mu_estimate,sd = em_r$sigma_estimate)

grDevices::pdf("Recovery-Comparison-t.pdf",width = 10, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)

plot(x = 0,y = -1,
     ylim = c(0,0.4),col="white"
     ,xlab="x",ylab = "Probability",bty="l",
     xlim=c(-15,15))

# Method-1 # 
plot1 = Log_Con_Discr_Algorithm_method1(n_i,t_grid)
lcm = gcmlcm(x = plot1$mle_x,  y = plot1$mle_psi,type = "lcm")
lines(x = plot1$x_knots,y = exp(lcm$y.knots) * plot1$Total_Coefficient,col="brown3", lwd = 3)

# Method-2 # 
plot2 = Log_Con_Discr_Algorithm_method2(n_i,t_grid)
lcm = gcmlcm(x = plot2$mle_x,  y = plot2$mle_psi,type = "lcm")
lines(x = plot2$x_knots,y = lcm$y.knots * plot2$Total_Coefficient,col="purple3", lwd = 3)

# Method-3 # 
plot3 = Log_Con_Discr_Algorithm_method3(n_i,t_grid)
x = plot3$x
f_3 = exp(plot3$phi)

lines(x,f_3,type="l",col="darkorange3",lwd = 3)

# EM # 
# plot(g,xlim = c(6,14),add = T, col = "lightblue",lwd =3)
plot(h,xlim = c(-15,15),add = T, col = "navy", lwd = 3)

# True Distribution # 
plot(f,xlim=c(-15,15),add=T, col = "darkgreen", lwd = 3)

#legend("topright", 
#       legend = c("ALCDE1","ALCDE2","ACSCi","EM","True"), 
#       col = c("brown3","purple","darkorange3","navy","darkgreen"), 
#       lty = c(1,1,1,1,1),
#       lwd = c(3,3,3,3,3),
#       bty = "n", 
#       pt.cex = 2, 
#       cex = 3, 
#       text.col = "black", 
#       horiz = F, 
#       inset = c(0, 0),x.intersp = 0.2, y.intersp = 1,seg.len = 0.75) # -0.3
dev.off()

# Beta Distribution 

f = function(x) dbeta(x,alpha,beta)

set.seed(9)

x <- sort(rbeta(n,alpha,beta))

MM = qbeta(0.9999999,alpha,beta) + 3
beta_grid <- seq(0,1.1, by = beta_deltas[1])

n_i <- hist(x, breaks = beta_grid, plot = F)$counts

e = EM_Algorithm_Format(n_i,beta_grid)
em_r = em(e$bl, e$bu, e$Freq, theta_init = e$Theta)
g = function(x) dnorm(x,mean = em_r$mu_estimate,sd = 1)
h = function(x) dnorm(x,mean = em_r$mu_estimate,sd = em_r$sigma_estimate)

grDevices::pdf("Recovery-Comparison-Beta.pdf",width = 10, height = 10)
par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=5,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=2,cex.axis=2,
    font.axis=1,cex.main=2)

plot(x = 0,y = -1,
     ylim = c(0,3),col="white"
     ,xlab="x",ylab = "Probability",bty="l",
     xlim=c(0,1))

# Method-1 # 
plot1 = Log_Con_Discr_Algorithm_method1(n_i,beta_grid)
lcm = gcmlcm(x = plot1$mle_x,  y = plot1$mle_psi,type = "lcm")
lines(x = plot1$x_knots,y = exp(lcm$y.knots) * plot1$Total_Coefficient,col="brown3", lwd = 3)

# Method-2 # 
plot2 = Log_Con_Discr_Algorithm_method2(n_i,beta_grid)
lcm = gcmlcm(x = plot2$mle_x,  y = plot2$mle_psi,type = "lcm")
lines(x = plot2$x_knots,y = lcm$y.knots * plot2$Total_Coefficient,col="purple3", lwd = 3)

# Method-3 # 
plot3 = Log_Con_Discr_Algorithm_method3(n_i,beta_grid)
x = plot3$x
f_3 = exp(plot3$phi)

lines(x,f_3,type="l",col="darkorange3",lwd = 3)

# EM # 
# plot(g,xlim = c(6,14),add = T, col = "lightblue",lwd =3)
plot(h,xlim = c(0,1),add = T, col = "navy", lwd = 3)

# True Distribution # 
plot(f,xlim=c(0,1),add=T, col = "darkgreen", lwd = 3)

#legend("topright", 
#       legend = c("ALCDE1","ALCDE2","ALCDE3","EM","True"), 
#       col = c("brown3","purple","darkorange3","navy","darkgreen"), 
#       lty = c(1,1,1,1,1),
#       lwd = c(3,3,3,3,3),
#       bty = "n", 
#       pt.cex = 2, 
#       cex = 3, 
#       text.col = "black", 
#       horiz = F, 
#       inset = c(0, 0),x.intersp = 0.2, y.intersp = 1,seg.len = 0.75) # -0.3
dev.off()