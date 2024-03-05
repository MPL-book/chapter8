## AH model simulation study with right censored data

library(RConics)
library(addhazard)
source("AdditiveHazardMPL.R") #load package for fitting AH MPL model


# function for generating data - right censoring only
generate_ch8_7 <- function(n){
  
  x_1 <- rbinom(n, 1, 0.6)
  x_2 <- runif(n,1,2)
  
  neg.logU <- -log(runif(n, 0.1, 1))
  
  y <- t <- delta <-  rep(0, n)
  for(i in 1:n){
    y[i] <- as.numeric(cubic(c(1, 0, -x_1[i]*0.5 + x_2[i]*1.5, -neg.logU[i]))[1])
  }
  
  c <- rexp(n,3)
  
  delta <- as.numeric(y < c)
  t <- c
  t[delta == 1] <- y[delta == 1]
  
  
  data = data.frame(id = rep(1:n), t, delta, x_1, x_2)
  return(data)
  
}

# loop for running simulations (Scenario 1)
for(s in 1:200){
  ah_data = generate_ch8_7(100)
  
  # fit L&Y model
  aalen.surv <- Surv(time = ah_data$t, event = ah_data$delta)
  ah.fit <- ah(aalen.surv ~ x_1 + x_2, data = ah_data,
               ties = FALSE)
  
  save_87[s,1] = ah.fit$coef[1]
  save_87[s,2] = ah.fit$coef[2]
  save_87[s,3] = ah.fit$se[1]
  save_87[s,4] = ah.fit$se[1]
  
  # fit MPL model
  ah_data_new <- cbind(rep(0,100), ah_data$t, ah_data$delta, ah_data$x_1, 
                       ah_data$x_2, ah_data$id)
  ah_data_new = data.frame(ah_data_new)
  colnames(ah_data_new) = c("start","stop", "event", "x_1", "x_2", "id")
  #ah_data_new$stop[which(ah_data_new$event == 0)] <-  Inf
  
  ctrl<-ah_mpl.control(basis = "uniform", smooth = NULL, 
                       max.iter = c(10, 1000, 1000), tol = 1e-05, 
                       n.knots = c(50,0), range.quant = c(0.1, 0.9), 
                       min.theta = 1e-10, penalty = 2L, order = 3L, 
                       epsilon = c(1e-16, 1e-10), ties = "epsilon", seed = NULL) 
  
  
  fit_mpl1 <- ah_survmpl(Surv(ah_data_new$start, ah_data_new$stop, ah_data_new$event) ~ x_1 +  x_2, 
                         data = ah_data_new, control = ctrl)
  
  save_87[s,5] = fit_mpl1$coef$Beta[1]
  save_87[s,6] = fit_mpl1$coef$Beta[2]
  save_87[s,7] = fit_mpl1$se$Beta[1,2]
  save_87[s,8] = fit_mpl1$se$Beta[2,2]
  
}


