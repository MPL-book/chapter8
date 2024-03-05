## AH model with current status data

# Generate current status data
n <- 0
y <- t <- c <- delta <- x1 <- x2 <- NULL

while(n < 100){
  x1_i <- rbinom(1, 1, 0.2)
  x2_i <- runif(1, 0, 0.5)
  logU <- log(runif(1))
  t_i <- uniroot(target_func, c(0, 10), 
                 x1 = x1_i, x2 = x2_i, logU = logU)$root
  if(t_i > 0){
    c_i <- rexp(1,1) 
    x1 <- c(x1, x1_i)
    x2 <- c(x2, x2_i)
    n <- n + 1
    if(c_i < t_i){
      delta <- c(delta, 1)
      y <- c(y, c_i)
    }else if(t_i < c_i){
      delta <- c(delta, 0)
      y <- c(y, c_i)
    }
    t <- c(t, t_i)
    c <- c(c, c_i)
  }
}

# Fit model to current status data
cs.data <- data.frame(y = y, delta = delta, X1 = x1,
                      X2 = x2)
table(cs.data$delta)

library(survival)
cs.surv <- Surv(time = cs.data$y, event = cs.data$delta)
fit.tvc = coxph(cs.surv ~ tt(X1) + tt(X2), data = cs.data, tt = function(x, t,...) -t*x)
summary(fit.tvc)


