## Generate random numbers from an additive hazards model

# Function we need to solve for t
target_func <- function(t, x1, x2, logU){
  t^3 + (-0.5*x1 + 1.5*x2)*t + logU
}

# Loop to generate data
n <- 0
y <- delta <- x1 <- x2 <- NULL
while(n < 100){
  x1_i <- rbinom(1, 1, 0.6)
  x2_i <- runif(1, 1, 2)
  logU <- log(runif(1))
  t_i <- uniroot(target_func, c(0, 10), 
                 x1 = x1_i, x2 = x2_i, logU = logU)$root
  
  if(t_i > 0){
    c_i <- rexp(1,1.3)
    x1 <- c(x1, x1_i)
    x2 <- c(x2, x2_i)
    n <- n + 1
    if(c_i < t_i){
      delta <- c(delta, 0)
      y <- c(y, c_i)
    }else if(t_i < c_i){
      delta <- c(delta, 1)
      y <- c(y, t_i)
    }
  }
}

table(delta) #check censoring proportion
ah.data <- data.frame(y = y, delta = delta, X1 = x1, X2 = x2)

