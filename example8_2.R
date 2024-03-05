## R example for fitting AH model using Aalen's method

# Generate data using code from Example 8.1

# Load library and fit the model
library(timereg)
aalen.surv <- Surv(time = y, event = delta)
aalen.fit <- aalen(aalen.surv ~ X1 + X2, data = ah.data)

# Plot estimated cumulative beta functions (true in red)
plot(aalen.fit)
