# Lindsey data analysis

# Read in and set up data for analysis
Lindsey <- read.csv("Lindsey.csv")

n <- 31
y <- matrix(0,n,5)
colnames(y) <- c("Id", "star_time", "stop_time", "surtim_id", "cen_id")
y[,1] <- c(1:31)
y[,2] <- Lindsey[,1]
y[,3] <- Lindsey[,2]
y[,4] <- 1
y[,5] <- Lindsey[,3]
X <- matrix(0,31,4)
X <- Lindsey[,4:7]

# Fit AH model
control <- AH_MPL.control(n, smooth = 1000, n.obs_basis = 2, max.iter <- c(100, 5000, 50000), tol_1 = 1e-5, tol_2 = 1, tau = 1000, min.theta = 1e-10)
fit <- AH_MPL(y, X, control) 
fit$Beta

# Fit L&Y Model
library(addhazard)
head(Lindsey)
Lind.cs <- data.frame(y = Lindsey[,2],
                      delta = as.numeric(Lindsey[,3]>0),
                      x1 = Lindsey[,4], x2 = Lindsey[,5], 
                      x3 = Lindsey[,6], x4 = Lindsey[,7])


fit.Lind = coxph(Surv(y, delta) ~ tt(x1) + tt(x2) + tt(x3) + tt(x4), data = Lind.cs, tt = function(x, t,...) -t*x)
summary(fit.Lind)




