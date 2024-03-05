## Random numbers with time-fixed and time-varying covariates

n = 50 #sample size
N_i <- 10  #maximum number of intermittent points
id_long <- x_long <- z1_long <- z2_long <- NULL 
start <- end <- event <- NULL

for(i in 1:n){
  
  # generate from Unif[0,1]
  u_i <- runif(1)
  neg.log.u <- -log(u_i)
  
  # generate time fixed covariate
  x <- rbinom(1, 1, 0.5)
  
  # generate change point in z2
  
  tau <- runif(1, 0.5, 2.5)
  H_tau <- tau^3 + 0.5 * tau * x - 1.5*(tau-1)^2
  
  # generate true event time
  if(neg.log.u < H_tau){
    y_i <- uniroot(function(t) t^3 + 0.5 * t * x - 1.5*(t-1)^2, c(0, 10))$root
  }else{
    y_i <- uniroot(function(t) t^3 + 0.5 * t * x - 1.5*(t-1)^2 - (t-tau), c(0, 10))$root
  }
  
  
  # generate censoring time
  c_i <- runif(1, 0.3, 0.8)
  
  # save event or censoring time
  t_i <- min(c_i, y_i)
  
  # generate observation times
  n_i <- sample(1:N_i, 1)
  t_ia <- quantile(c(0,t_i), c(seq(0,1,length.out = n_i+1)))[2:(n_i+1)]
  if (tau < t_i) {
    t_ia <- sort(c(t_ia, tau))
  }
  n_i = length(t_ia)
  
  # generate time varying covariates
  z_1 <- t_ia + 1
  z_2 <- rep(0, length(t_ia))
  z_2[which(t_ia > tau)] <- 1
  
  # create 'long' dataset
  id_long_i <- rep(i, n_i)
  x_long_i <- rep(x, n_i)
  event_i <- rep(0, n_i)
  if(y_i < c_i){
    event_i[n_i] <- 1
  }
  
  if(n_i > 1){
    end_i <- c(t_ia[2:n_i], t_i)
  }else{
    end_i <- (t_i)
  }
  
  id_long <- c(id_long, id_long_i)
  x_long <- c(x_long, x_long_i)
  z1_long <- c(z1_long, z_1)
  z2_long <- c(z2_long, z_2)
  start <- c(start, t_ia)
  end <- c(end, end_i)
  event <- c(event, event_i)
  
}

df <- data.frame(id_long, x_long, z1_long, z2_long, start, end, event)
