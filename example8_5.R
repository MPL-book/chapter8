# Continuation of Example 8.4

#Fit the modified model to current status data
cs.data1 <- data.frame(y = y, delta = delta, X1 = x1,
                       X2 = x2)
cs.surv1 <- Surv(time = cs.data1$y, event = cs.data1$delta)

fit.tvc1 = coxph(cs.surv ~ tt(X1) + tt(X2) + X1 + X2, data = cs.data1, tt = function(x, t,...) -t*x)
summary(fit.tvc1)
