## Lin & Yang method for an AH model

# Generate data using code in Example 8.1

# Load library and fit the model
library(addhazard)
aalen.surv <- Surv(time = y, event = delta)
ah.fit <- ah(aalen.surv ~ X1 + X2, data = ah.data, ties = FALSE)
summary(ah.fit)

# Predict baseline hazard function
newdata <- data.frame(X1=0, X2 =0)
time <- predict(ah.fit, newdata, newtime = seq(from=0.1,to=2, by=0.1))$time
baseline_haz <- predict(ah.fit, newdata, newtime = seq(from=0.1,to=2, by=0.1))$L
plot(baseline_haz ~ time, type = "s", ylab = "baseline hazard")

