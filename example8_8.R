## AH model simulation study with interval censoring

library(RConics)
library(addhazard)
source("AdditiveHazardMPL.R") #load package for fitting AH MPL model

# Function for generating interval censored data

generate_ch8_7_intcens <- function(n, pi_E){
  
  x_1 <- rbinom(n, 1, 0.6)
  x_2 <- runif(n,1,2)
  
  neg.logU <- -log(runif(n, 0.1, 1))
  
  y <-  rep(0, n)
  for(i in 1:n){
    y[i] <- as.numeric(cubic(c(1, 0, -x_1[i]*0.5 + x_2[i]*1.5, -neg.logU[i]))[1])
  }
  
  L <- runif(n)
  R <- L + runif(n)
  u_i <- runif(n)
  
  delta <- t_L <- t_R <- rep(0,n)
  
  for(i in 1:n){
    if(u_i[i] < pi_E){
      delta[i] <- 1
      t_L[i] <- 0
      t_R[i] <- y[i]
    }else if(y[i] < L[i]){
      delta[i] <- 2
      t_L[i] <- 0
      t_R[i] <- L[i]
    }else if(L[i] < y[i] & y[i] < R[i]){
      delta[i] <- 3
      t_L[i] <- L[i]
      t_R[i] <- R[i]
    }else if(R[i] < y[i]){
      delta[i] = 0
      t_L[i] <- 0
      t_R[i] <- R[i]
    }
  }
  
  data <- data.frame(id = rep(1:n), t_L, t_R, delta, x_1, x_2)
  return(data)
  
}


#set up simulation study

save <- matrix(0, ncol = 13, nrow = 300)
save_h0 <- NULL

for(s in 1:300){
  
  #generate and set up data
  dat <- generate_ch8_7_intcens(300, pi_E = 0.8)
  
  n <- 300
  y <- matrix(0,n,5)
  colnames(y) <- c("Id", "star_time", "stop_time", "surtim_id", "cen_id")
  y[,1] <- c(1:n)
  y[,2] <- dat$t_L
  y[,3] <- dat$t_R
  y[,4] <- 1
  y[,5] <- dat$delta
  
  X <- matrix(0,n,2)
  X <- dat[,5:6]
  
  # fit model
  control <- AH_MPL.control(n, smooth = 1000, n.obs_basis = 50, max.iter = c(100, 5000, 50000), 
                           tol_1 = 1e-5, tol_2 = 1, tau = 1000, min.theta = 1e-10)
  
  fit <- AH_MPL(y, X, control) 
  
  # regression coefficients
  save[s,1] <- fit$Beta[1]
  save[s,2] <- fit$Beta[2]
  
  save[s,3] <- fit$se.beta[1]
  save[s,4] <- fit$se.beta[2]
  
  # baseline hazard estimation
  h0.ind <- c(round(quantile(1:length(fit$bins$Alpha), c(0.25, 0.5, 0.75))))
  
  save[s,5] <- c(fit$Theta)[h0.ind[1]]
  save[s,6] <- c(fit$Theta)[h0.ind[2]]
  save[s,7] <- c(fit$Theta)[h0.ind[3]]
  
  save[s,8] <- (fit$bins$Alpha[-length(fit$bins$Alpha)])[h0.ind[1]]
  save[s,9] <- (fit$bins$Alpha[-length(fit$bins$Alpha)])[h0.ind[2]]
  save[s,10] <- (fit$bins$Alpha[-length(fit$bins$Alpha)])[h0.ind[3]]
  
  save[s,11] <- sqrt(diag(fit$cov.theta))[h0.ind[1]]
  save[s,12] <- sqrt(diag(fit$cov.theta))[h0.ind[2]]
  save[s,13] <- sqrt(diag(fit$cov.theta))[h0.ind[3]]
  
  save_h0 <- cbind(save_h0, fit$bins$Alpha[2:21],c(fit$Theta)[1:20])
  #select length of alpha and theta to save based on sample size & number of basis functions used
  
}

