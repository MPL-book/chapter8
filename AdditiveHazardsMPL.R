# R code for additive hazards model
# AH_MPL.control: Function to control model parameters 
# n: sample size
# smooth: smooth parameter
# n.obs_basis: number of observed time per uniform basis
# max.iter: a vector with the first value defining the maximum number of iterations for the smooth parameter,
#           the second value the maximum number of iterations for the Theta, Beta and Gamma parameters to attempt to converge,
#           the third value is the total number of iterations allowed.
# tol_1: convergence tolerance value for the PD interior point algorithm.
# tol_2: convergence tolerance value for the smooth value selection algorithm
# tau: parameter in PD interior point algorithm
# min.theta: the minimal baseline hazard estimate value, below which the theta estimates will be considered to be 0.
AH_MPL.control = 
  function (n, smooth = 2000, n.obs_basis = 10, max.iter = c(100, 300, 5000), 
            tol_1 = 1e-10, tol_2 = 1e-5, tau = 1000, min.theta = 1e-10) 
  {
    if (!is.null(smooth)) {
      smooth = ifelse(smooth < 0, 1000, smooth)
    }
    else {
      smooth = 1000
    }
    if (!is.null(n.obs_basis)) {
      n.obs_basis = ifelse(n.obs_basis < 1 | n.obs_basis > 
                             floor(n/2), floor(n^(1/3) - 2), round(n.obs_basis))
    }
    else {
      n.obs_basis = floor(n^(1/3) - 2)
    }
    max.iter = c(ifelse(max.iter[1] > 0, as.integer(max.iter[1]), 100), 
                 ifelse(max.iter[2] > 0, as.integer(max.iter[2]), 300), 
                 ifelse(length(max.iter) == 2, 5000, 
                        ifelse(max.iter[3] > ifelse(max.iter[2] > 0, as.integer(max.iter[2]), 5000), as.integer(max.iter[3]), 5000)))
    tol_1 = ifelse(tol_1 > 0 & tol_1 < 1, tol_1, 1e-10)
    tol_2 = ifelse(tol_2 > 0, tol_2, 1e-3)
    min.theta = ifelse(min.theta > 0 & min.theta < 0.001, min.theta, 
                       1e-10)
    out = list(smooth = smooth, n.obs_basis = n.obs_basis, max.iter = max.iter, 
               tol_1 = tol_1, tol_2 = tol_2, tau = tau, min.theta = min.theta)
    out
  }
# bins_mpl: Function to determine bins and bin length  
# input:
#   control: model parameters in AH_MPL.control (control = AH_MPL.control())
#   t_i1: observed time for left censoring, uncensoring, right censoring and left-end point of finite interval censoring
#   t_i2: observed time for right-end point of finite interval censoring
#   ctype: a matrix for censoring indicators
#   example:              ctype
#          t_i1    t_i2     r    e    l     i
#           a        1      NA   NA   TRUE  NA   (left censored at a)
#           b        1      TRUE NA   NA    NA   (right censored at b)
#           c        d      NA   NA   NA    TRUE (finite interval censored between (c, d))
#           f        1      NA   TRUE NA    NA   (exactly observed at f)
# output:
#   m: number of bins
#   Alpha: edge points at each bin
#   Delta: bin length
bins_mpl = function (control, t_i1, t_i2, ctype) 
{
  event = c(t_i1[ctype[, "e"]], t_i1[ctype[, "i"]], t_i2[ctype[, "i"]], 
            t_i1[ctype[, "l"]])
  count = control$n.obs_basis
  order_event = sort(event)
  n.obs = length(order_event)
  m = ceiling(n.obs/count)
  Alpha = matrix(0, m + 1, 1)
  for (j in 2:m){
    Alpha[j] = order_event[(j-1) * count]
  }
  t.obs = c(t_i1[ctype[, "i"]], t_i2[ctype[, "i"]], 
            t_i1[ctype[, "e"]], t_i1[ctype[, "r"]], t_i1[ctype[, "l"]])
  Alpha[m + 1] = max(t.obs)
  Delta = diff(Alpha)
  list(m = m, Alpha = Alpha, Delta = Delta)
}


# penalty_mpl: R matrix in the penalty term 
# input:
#   bins: output of the function bins_mpl (bins = bins_mpl(control, t_i1, t_i2, ctype))
# output:
#   M_R_mm: R matrix in the penalty term
penalty_mpl = function (bins) 
{
  m = bins$m
  M_R_mm = matrix(0, m, m)
  for (u in 1:m){
    M_R_mm[u, u] = 6
  }
  for (u in 1: (m-1)){
    M_R_mm[u+1, u] = -4
    M_R_mm[u, u+1] = -4
  }
  for (u in 1: (m-2)){
    M_R_mm[u, u+2] = 1
    M_R_mm[u+2, u] = 1
  }
  M_R_mm[1, 1] = M_R_mm[m, m] = 1
  M_R_mm[2, 2] = M_R_mm[m-1, m-1] = 5
  M_R_mm[2, 1] = M_R_mm[1, 2] = M_R_mm[m, m-1] = M_R_mm[m-1, m] = -2
  return(M_R_mm)
}
# basis_mpl: Function to determine basis for approximating baseline hazards
# input: 
#   x: time points to be evaluated
#   bins: output of the function bins_mpl
#   which: 1 for baseine hazard and 2 for cumulative haseline hazard
# output:
#   M_psi_nm: basis for baseline hazard
#   M_Psi_nm: basis for cumulative haseline hazard
basis_mpl = function(x, bins, which = c(1, 2)) 
{
  which.matrix = rep(T, 2)
  which.matrix[-which] = FALSE
  n = length(x)
  Alpha = bins$Alpha
  Delta = bins$Delta
  m = bins$m
  M_Psi_nm = M_psi_nm = matrix(0, n, m)
  u_i = sapply(x, function(y, lim = Alpha[-1L]) sum(lim < y) + 1L)
  for (i in 1:n) {
    M_psi_nm[i, u_i[i]] = 1
    M_Psi_nm[i, 1:u_i[i]] = c(if (u_i[i] > 1) {
      Delta[1:(u_i[i] - 1)]
    }, x[i] - Alpha[u_i[i]])
  }
  if (all(which.matrix)) {
    list(psi = M_psi_nm, Psi = M_Psi_nm)
  }
  else {
    if (which.matrix[1]) {
      M_psi_nm
    }
    else {
      M_Psi_nm
    }
  }
}
# AH_MPL: Function to fit AH model 
# input: 
#   y: 1st column contains the ID number for subjects
#      2nd column contains starting time of monitoring time interval 
#      3rd column contains stopping time of monitoring time interval 
#      4th column contains indicators, value of 1 indicates the stopping time of the corresponding monitoring time interval
#          is the observed survival time and 0 the covariate changing time. (look at, particularly, the finite interval censored example below)
#      5th column contains indicator variable. For each subject, at the end of each monitor time interval, the indicator value is 0 until
#          at the last monitor time interval, at the end of which subject ceased to be investigated, 
#          the value is 0 for right censoring, 1 for event, 2 for left-censoring, or 3 for interval-censoring. 
#   X: covariates in AH regression
#   control: model parameters in function AH_MPL.control
#   Example for y, and X: 
#   y                                                      X (time varying)                
#   Id  2nd   3rd  4th  5th                                co-variate 1      co-variate 2   
#   1    0     2    0    0                                    x_11              x_12           
#   1    2     3    1    1 (uncensored at 3)                  y_11              y_12           
#   2    0     1    0    0                                    x_21              x_22
#   2    1     3    1    2 (left censored at 3)               y_21              y_22           
#   3    0     1    0    0                                    x_31              x_32
#   3    1     2    0    0                                    y_31              y_32              
#   3    2     4    1    0                                    y_31              y_32
#   3    4     5    0    0                                    z_31              z_32
#   3    5     6    1    3 (interval censored in (4, 6))      w_31              w_32           
#   4    0     2    0    0                                    x_41              x_42
#   4    2     4    1    0 (right censored at 4)              y_41              y_42           
# output:
#   cvg: matrix with the 1st column containing the number of internal iterations, 
#        the 2nd column containing mu values in each internal iteration,
#        the 3rd column containing penalized log-likelihood in each internal iteration.
#   Beta: coefficient estimates in AH regression
#   Theta: basis estimates
#   cov.theta: covariance matrix of basis estimates
#   se.beta: standard error of coefficient estimates in AH regression
#   t.obs: collection of observed data containing uncensored, left-censored, interval-censored and right censored data.
AH_MPL = 
  function (y, X, control) 
  {
    # Observed survival data with censoring type 
    n = length(unique(y[, 1]))
    t_i1 = t_i2  = status = matrix(1, n, 1)
    for (i in 1:n) {
      y.subject = y[which(y[, 1]==i), , drop=FALSE]
      n.monitor.time = nrow(y.subject)
      censoring_status = y[, 5][which(y[, 1] == i)]
      if (censoring_status[n.monitor.time] == 3) {
        interval_left_id = y.subject[, 4]
        t_i1[i] = y.subject[interval_left_id == 1, 3][1]
        t_i2[i] = y.subject[n.monitor.time, 3]
        status[i] = 3
      }
      else {
        t_i1[i] = y.subject[n.monitor.time, 3]
        status[i] = censoring_status[n.monitor.time]
      }
    }
    ctype = matrix(NA, nrow = n, ncol = 4)
    colnames(ctype) = c("r", "e", "l", "i")
    for (tw in 1:4) {
      ctype[, tw] = status == (tw - 1)
    }
    n.ctype = apply(ctype, 2, sum)
    ties_1 = duplicated(t_i1)
    if (any(ties_1)) {
      #      seed = 666
      t_i1[ties_1] = t_i1[ties_1] + runif(sum(ties_1), 0, 1e-5)
    }
    t_i2_r = t_i2[ctype[, 4]]
    ties_2 = duplicated(t_i2_r)
    if (any(ties_2)) {
      #      seed = 888
      t_i2_r[ties_2] = t_i2_r[ties_2] + runif(sum(ties_2), 0, 1e-5)
      t_i2[ctype[, 4]] = t_i2_r
    }
    monitor_time_covar = matrix(subset(y[, 3], y[, 4] == 0))
    n_monitor_covar = nrow(monitor_time_covar)
    # Covariate observations 
    beta_dim = ncol(X)
    xt1 = Xt1 = matrix(0, n, 2 * beta_dim)
    for (i in 1 : n ) {
      M.subject = X[which(y[, 1] == i), , drop = FALSE]
      y.subject = y[which(y[, 1]==i), , drop=FALSE]
      n.monitor.time = nrow(y.subject)
      if (status[i] == 3) {
        interval_left_id = y.subject[, 4]
        xt1[i, 1:beta_dim] = as.numeric(M.subject[interval_left_id == 1, ][1, ])
        xt1[i, (beta_dim + 1):(2 * beta_dim)] = as.numeric(M.subject[n.monitor.time, ])
        Xt1[i, 1:beta_dim] = colSums((y.subject[1 : which(interval_left_id == 1)[1], 3] - 
                                        y.subject[1 : which(interval_left_id == 1)[1], 2]) *
                                       M.subject[1 : which(interval_left_id == 1)[1], , drop = FALSE])
        Xt1[i, (beta_dim + 1):(2 * beta_dim)] = colSums((y.subject[, 3] - y.subject[, 2]) * M.subject)
      }
      else {
        xt1[i, 1:beta_dim] = as.numeric(M.subject[n.monitor.time, ]) 
        Xt1[i, 1:beta_dim] = colSums((y.subject[, 3] - y.subject[, 2]) * M.subject)
      }
    }
    M_xt1_nrp = xt1[ctype[, 1], (1:beta_dim), drop = F]
    M_xt1_nop = xt1[ctype[, 2], (1:beta_dim), drop = F]
    M_xt1_nlp = xt1[ctype[, 3], (1:beta_dim), drop = F]
    M_xt1_nip = xt1[ctype[, 4], (1:beta_dim), drop = F]
    M_xt2_nip = xt1[ctype[, 4], ((1 + beta_dim):(2 * beta_dim)), drop = F]
    M_Xt1_nrp = Xt1[ctype[, 1], (1:beta_dim), drop = F]
    M_Xt1_nop = Xt1[ctype[, 2], (1:beta_dim), drop = F]
    M_Xt1_nlp = Xt1[ctype[, 3], (1:beta_dim), drop = F]
    M_Xt1_nip = Xt1[ctype[, 4], (1:beta_dim), drop = F]
    M_Xt2_nip = Xt1[ctype[, 4], ((1 + beta_dim):(2 * beta_dim)), drop = F]
    M_txt1_nop = t(M_xt1_nop)
    M_tXt1_nrp = t(M_Xt1_nrp)
    M_tXt1_nop = t(M_Xt1_nop)
    M_tXt1_nlp = t(M_Xt1_nlp)
    M_tXt1_nip = t(M_Xt1_nip)
    M_tXt2_nip = t(M_Xt2_nip)
    M_x_m = subset(X, y[, 4] == 0)
    # Basis functions 
    bins = bins_mpl(control, t_i1, t_i2, ctype) 
    M_R_mm = penalty_mpl(bins)
    m = bins$m
    M_Rstar_ll = rbind(cbind(M_R_mm, matrix(0, m, beta_dim)),
                       matrix(0, beta_dim, m + beta_dim))
    M_psi_nom = basis_mpl(t_i1, bins, 
                          which = 1)[ctype[, 2], , drop = F]
    M_tpsi_nom = t(M_psi_nom)
    M_psi_nrm = basis_mpl(t_i1, bins,  
                          which = 1)[ctype[, 1], , drop = F]
    M_psi_nlm = basis_mpl(t_i1, bins, 
                          which = 1)[ctype[, 3], , drop = F]
    M_psi1_nim = basis_mpl(t_i1, bins, 
                           which = 1)[ctype[, 4], , drop = F]
    M_psi2_nim = basis_mpl(t_i2, bins,  
                           which = 1)[ctype[, 4], , drop = F]
    M_Psi_nom = basis_mpl(t_i1, bins,  
                          which = 2)[ctype[, 2], , drop = F]
    M_tPsi_nom = t(M_Psi_nom)
    M_Psi_nrm = basis_mpl(t_i1, bins,  
                          which = 2)[ctype[, 1], , drop = F]
    M_tPsi_nrm = t(M_Psi_nrm)
    M_Psi_nlm = basis_mpl(t_i1, bins,  
                          which = 2)[ctype[, 3], , drop = F]
    M_tPsi_nlm = t(M_Psi_nlm)
    M_Psi1_nim = basis_mpl(t_i1, bins, 
                           which = 2)[ctype[, 4], , drop = F]
    M_Psi2_nim = basis_mpl(t_i2, bins,  
                           which = 2)[ctype[, 4], , drop = F]
    M_tPsi1_nim = t(M_Psi1_nim)
    M_tPsi2_nim = t(M_Psi2_nim)
    if (n_monitor_covar == 0) {
      M_psi_m = NULL
    }
    else {
      M_psi_m = basis_mpl(monitor_time_covar, bins, which = 1)
    }
    # First derivatives of constraint function 
    C_M_h0 = - cbind(diag(1,m), matrix(0, m, beta_dim))
    C_M_ho = - cbind(M_psi_nom, M_xt1_nop)
    C_M_hr = - cbind(M_psi_nrm, M_xt1_nrp)
    C_M_hl = - cbind(M_psi_nlm, M_xt1_nlp)
    C_M_hi1 = - cbind(M_psi1_nim, M_xt1_nip)
    C_M_hi2 = - cbind(M_psi2_nim, M_xt2_nip)
    if (n_monitor_covar == 0) {
      C_M_m = NULL
    }
    else {
      C_M_m = - cbind(M_psi_m, M_x_m)
    }
    C_M = rbind(C_M_h0, C_M_ho, C_M_hr, C_M_hl, C_M_hi1, C_M_hi2, C_M_m)
    # initial smooth parameter 
    s_gamma = control$smooth
    # number of internal iterations
    K = control$max.iter[2]
    # parameters of the PD algorithm
    tau = control$tau
    xi = 0.5
    # tolerance in the internal algorithm 
    s_convlimit = control$tol_1
    # initial value of denominator in smooth parameter selection algorithm 
    s_df = -1
    # initial value of iteration steps 
    full.iter = 0
    # algorithm 
    for (iter in 1: control$max.iter[1]) {
      # initial values 
      M_theta_m1 = matrix(1, nrow = m, ncol = 1)
      M_beta_p1 = matrix(0, nrow = beta_dim, ncol = 1)
      M_lambda = M_s = matrix(1, nrow = m + n + n.ctype[4] + n_monitor_covar, ncol = 1)
      M_eta = rbind(M_theta_m1, M_beta_p1)
      mu = sum(M_lambda * M_s)/as.numeric(m + n + n.ctype[4] + n_monitor_covar)
      # Baseline and cumulative baseline hazard 
      M_h0_no1 = M_psi_nom %*% M_theta_m1
      M_H0_nl1 = M_Psi_nlm %*% M_theta_m1
      M_H01_ni1 = M_Psi1_nim %*% M_theta_m1
      M_H02_ni1 = M_Psi2_nim %*% M_theta_m1
      # hazard and cumulative Hazard 
      M_h_no1 = M_h0_no1 + M_xt1_nop %*% M_beta_p1
      M_H_nl1 = M_H0_nl1 + M_Xt1_nlp %*% M_beta_p1
      M_H1_ni1 = M_H01_ni1 + M_Xt1_nip %*% M_beta_p1
      M_H2_ni1 = M_H02_ni1 + M_Xt2_nip %*% M_beta_p1
      # Survival functions 
      M_S_nl1 = exp(-M_H_nl1)
      M_S1_ni1 = exp(-M_H1_ni1)
      M_S2_ni1 = exp(-M_H2_ni1)
      M_S_nl1[M_S_nl1 == 1] = 1 - 1e-16
      M_S1_ni1[M_S1_ni1 == 1] = 1 - 1e-16
      M_S2_ni1[M_S2_ni1 == 1] = 1 - 1e-16
      M_S1mS2_ni1 = M_S1_ni1 - M_S2_ni1
      M_S1mS2_ni1[M_S1mS2_ni1 < 1e-10] = 1e-10
      # First derivative of objective function
      M_gradtheta_m1 = M_tpsi_nom %*% (1/M_h_no1) - M_tPsi_nom %*% matrix(1, nrow = n.ctype[2], ncol = 1) -
        M_tPsi_nrm %*% matrix(1, nrow = n.ctype[1], ncol = 1) + M_tPsi_nlm %*% (M_S_nl1/(1 - M_S_nl1)) + 
        M_tPsi2_nim %*% (M_S2_ni1/M_S1mS2_ni1) - M_tPsi1_nim %*% (M_S1_ni1/M_S1mS2_ni1) 
      M_gradbeta_p1 = (M_txt1_nop %*% (1 / M_h_no1) - M_tXt1_nop %*% matrix(1, nrow = n.ctype[2], ncol = 1)) - M_tXt1_nrp %*% 
        matrix(1, nrow = n.ctype[1], ncol = 1) + M_tXt1_nlp %*% (M_S_nl1/(1 - M_S_nl1)) + (M_tXt2_nip %*% (M_S2_ni1/(M_S1mS2_ni1)) -
                                                                                             M_tXt1_nip %*% (M_S1_ni1/(M_S1mS2_ni1)))
      Rtheta = M_R_mm %*% M_theta_m1
      TwoLRtheta = s_gamma * 2 * Rtheta
      PL_FD =  -(rbind(M_gradtheta_m1, M_gradbeta_p1) - rbind(TwoLRtheta, matrix(0, beta_dim)))
      # collection of mu values and penalized log-likelihood
      cvg = matrix(0, K, 3)
      # internal algorithm
      for (k in 1:K) {
        # second derivative for objective function
        H = matrix(0, beta_dim + m, beta_dim + m)
        H[1:m, 1:m] = - M_tpsi_nom %*% diag(c(1/M_h_no1^2), n.ctype[2], n.ctype[2]) %*% M_psi_nom -
          M_tPsi_nlm %*% diag(c(M_S_nl1/(1-M_S_nl1)^2), n.ctype[3], n.ctype[3]) %*% M_Psi_nlm - 
          (M_tPsi2_nim - M_tPsi1_nim) %*% diag(c(M_S1_ni1 * M_S2_ni1/M_S1mS2_ni1^2), 
                                               n.ctype[4], n.ctype[4]) %*% (M_Psi2_nim - M_Psi1_nim)
        H[(m+1):(m+beta_dim), 1:m] = - M_txt1_nop %*% diag(c(1/M_h_no1^2), n.ctype[2], n.ctype[2]) %*% M_psi_nom -
          M_tXt1_nlp %*% diag(c(M_S_nl1/(1-M_S_nl1)^2), n.ctype[3], n.ctype[3]) %*% M_Psi_nlm -
          (M_tXt2_nip - M_tXt1_nip) %*% diag(
            c(M_S1_ni1 * M_S2_ni1/M_S1mS2_ni1^2), n.ctype[4], n.ctype[4]) %*% (M_Psi2_nim - M_Psi1_nim)
        H[(m+1):(m+beta_dim), (m+1):(m+beta_dim)] = - M_txt1_nop %*% diag(c(1/M_h_no1^2), n.ctype[2], n.ctype[2]) %*% M_xt1_nop -
          M_tXt1_nlp %*% diag(c(M_S_nl1/(1-M_S_nl1)^2), n.ctype[3], n.ctype[3]) %*% M_Xt1_nlp -
          (M_tXt2_nip - M_tXt1_nip) %*% diag(c(M_S1_ni1 * M_S2_ni1/M_S1mS2_ni1^2), n.ctype[4], n.ctype[4]) %*% (M_Xt2_nip - M_Xt1_nip)
        H[1:m, (m+1):(m+beta_dim)] = t(H[(m+1):(m+beta_dim), 1:m])
        SD_PL = -(H - s_gamma * 2 * M_Rstar_ll)
        # Newton direction
        D = solve(diag(c(M_s))) %*% diag(c(M_lambda))
        D_eta = solve(SD_PL + t(C_M) %*% D %*% C_M) %*% 
          (-PL_FD - t(C_M) %*% M_lambda - t(C_M) %*% D %*% (C_M %*% M_eta) - 
             xi * mu * t(C_M) %*% solve(diag(c(M_s))) %*% matrix(1, n + m + n.ctype[4] + n_monitor_covar, 1))
        D_lambda = solve(diag(c(M_s))) %*% (diag(c(M_lambda)) %*% C_M %*% D_eta + diag(c(M_lambda)) %*% (C_M %*% M_eta) +
                                              xi * mu * matrix(1, n + m + n.ctype[4] + n_monitor_covar, 1))
        D_s = -solve(diag(c(M_lambda))) %*% (diag(c(M_lambda)) %*% diag(c(M_s)) %*% matrix(1, n + m + n.ctype[4] + n_monitor_covar, 1) - 
                                               xi * mu * matrix(1, n + m + n.ctype[4] + n_monitor_covar, 1) + diag(c(M_s)) %*% D_lambda)
        # New estimates
        M_eta_OLD = M_eta
        M_lambda_OLD = M_lambda
        M_s_OLD = M_s
        mu_OLD = mu
        mu_bound = (1 - 0.01) * mu_OLD
        alpha = 1
        repeat{
          M_eta = M_eta_OLD + alpha * D_eta
          M_theta_m1 = matrix(M_eta[1 : m])
          M_beta_p1 = matrix(M_eta[(m + 1) : (m + beta_dim)])
          M_lambda = M_lambda_OLD + alpha * D_lambda
          M_s = M_s_OLD + alpha * D_s
          mu = sum(M_lambda * M_s)/(n + m + n.ctype[4] + n_monitor_covar)
          # new baseline hazard and cumulative hazard
          M_h0_no1 = M_psi_nom %*% M_theta_m1
          M_h0_nr1 = M_psi_nrm %*% M_theta_m1
          M_h0_nl1 = M_psi_nlm %*% M_theta_m1
          M_h01_ni1 = M_psi1_nim %*% M_theta_m1
          M_h02_ni1 = M_psi2_nim %*% M_theta_m1
          if (n_monitor_covar == 0) {
            M_h0_m = NULL
          }
          else {
            M_h0_m = M_psi_m %*% M_theta_m1
          }
          M_H0_no1 = M_Psi_nom %*% M_theta_m1
          M_H0_nr1 = M_Psi_nrm %*% M_theta_m1
          M_H0_nl1 = M_Psi_nlm %*% M_theta_m1
          M_H01_ni1 = M_Psi1_nim %*% M_theta_m1
          M_H02_ni1 = M_Psi2_nim %*% M_theta_m1
          # new hazard and cumulative Hazard
          M_h_no1 = M_h0_no1 + M_xt1_nop %*% M_beta_p1
          M_h_nr1 = M_h0_nr1 + M_xt1_nrp %*% M_beta_p1
          M_h_nl1 = M_h0_nl1 + M_xt1_nlp %*% M_beta_p1
          M_h1_ni1 = M_h01_ni1 + M_xt1_nip %*% M_beta_p1
          M_h2_ni1 = M_h02_ni1 + M_xt2_nip %*% M_beta_p1
          if (n_monitor_covar == 0) {
            M_h_m = NULL
          }
          else {
            M_h_m = M_h0_m + M_x_m %*% M_beta_p1
          }          
          M_H_no1 = M_H0_no1 + M_Xt1_nop %*% M_beta_p1
          M_H_nr1 = M_H0_nr1 + M_Xt1_nrp %*% M_beta_p1
          M_H_nl1 = M_H0_nl1 + M_Xt1_nlp %*% M_beta_p1
          M_H1_ni1 = M_H01_ni1 + M_Xt1_nip %*% M_beta_p1
          M_H2_ni1 = M_H02_ni1 + M_Xt2_nip %*% M_beta_p1
          # new survival
          M_S_nr1 = exp(-M_H_nr1)
          M_S_nl1 = exp(-M_H_nl1)
          M_S1_ni1 = exp(-M_H1_ni1)
          M_S2_ni1 = exp(-M_H2_ni1)
          M_S_nr1[M_S_nr1 == 1] = 1 - 1e-16
          M_S_nl1[M_S_nl1 == 1] = 1 - 1e-16
          M_S1_ni1[M_S1_ni1 == 1] = 1 - 1e-16
          M_S2_ni1[M_S2_ni1 == 1] = 1 - 1e-16
          M_S1mS2_ni1 = M_S1_ni1 - M_S2_ni1
          M_S1mS2_ni1[M_S1mS2_ni1 < 1e-10] = 1e-10
          # new first derivative of objective function
          M_gradtheta_m1 = M_tpsi_nom %*% (1/M_h_no1) - M_tPsi_nom %*% matrix(1, nrow = n.ctype[2], ncol = 1) -
            M_tPsi_nrm %*% matrix(1, nrow = n.ctype[1], ncol = 1) + M_tPsi_nlm %*% (M_S_nl1/(1 - M_S_nl1)) + 
            M_tPsi2_nim %*% (M_S2_ni1/M_S1mS2_ni1) - M_tPsi1_nim %*% (M_S1_ni1/M_S1mS2_ni1) 
          M_gradbeta_p1 = M_txt1_nop %*% (1/M_h_no1) - M_tXt1_nop %*% matrix(1, nrow = n.ctype[2], ncol = 1) - M_tXt1_nrp %*% 
            matrix(1, nrow = n.ctype[1], ncol = 1) + M_tXt1_nlp %*% (M_S_nl1/(1 - M_S_nl1)) + M_tXt2_nip %*% (M_S2_ni1/M_S1mS2_ni1) -
            M_tXt1_nip %*% (M_S1_ni1/M_S1mS2_ni1)
          Rtheta = M_R_mm %*% M_theta_m1
          TwoLRtheta = s_gamma * 2 * Rtheta
          PL_FD = -(rbind(M_gradtheta_m1, M_gradbeta_p1) - rbind(TwoLRtheta, matrix(0, beta_dim)))
          # new constraint functions
          f = - rbind(M_theta_m1, M_h_no1, M_h_nr1, M_h_nl1, M_h1_ni1, M_h2_ni1, M_h_m)
          # neighbour
          if (sqrt(t(PL_FD + t(C_M) %*% M_lambda) %*% (PL_FD + t(C_M) %*% M_lambda)) <= tau * mu_OLD & 
              sqrt(t(f + M_s) %*% (f + M_s)) <= tau * mu_OLD & 
              all(c(M_lambda * M_s) >= 0.5 *mu_OLD) &
              all(c(M_lambda, M_s) >= 0) & mu <= mu_bound)
            break
          # select the step length
          kapa = 0.6
          alpha = alpha * kapa
          # mu bound
          mu_bound = (1 - 0.01 * alpha) * mu_OLD
        }
        cvg[k, 1] = k
        cvg[k, 2] = mu
        cvg[k, 3] = sum(log(M_h_no1) - M_H_no1) - sum(M_H_nr1) + sum(log(1 - M_S_nl1)) + 
          sum(log(M_S1mS2_ni1)) - s_gamma * t(M_theta_m1) %*% M_R_mm %*% M_theta_m1
        print(k)
        if (cvg[k, 2] < s_convlimit){
          cvg = cvg[1:k, ]
          break
        }
      }
      # External Iteration 
      H = HRinv = matrix(0, beta_dim + m, beta_dim + m)
      H[1:m, 1:m] = - M_tpsi_nom %*% diag(c(1/M_h_no1^2), n.ctype[2], n.ctype[2]) %*% M_psi_nom -
        M_tPsi_nlm %*% diag(c(M_S_nl1/(1-M_S_nl1)^2), n.ctype[3], n.ctype[3]) %*% M_Psi_nlm - 
        (M_tPsi2_nim - M_tPsi1_nim) %*% diag(c(M_S1_ni1 * M_S2_ni1/M_S1mS2_ni1^2), 
                                             n.ctype[4], n.ctype[4]) %*% (M_Psi2_nim - M_Psi1_nim)
      H[(m+1):(m+beta_dim), 1:m] = - M_txt1_nop %*% diag(c(1/M_h_no1^2), n.ctype[2], n.ctype[2]) %*% M_psi_nom -
        M_tXt1_nlp %*% diag(c(M_S_nl1/(1-M_S_nl1)^2), n.ctype[3], n.ctype[3]) %*% M_Psi_nlm -
        (M_tXt2_nip - M_tXt1_nip) %*% diag(
          c(M_S1_ni1 * M_S2_ni1/M_S1mS2_ni1^2), n.ctype[4], n.ctype[4]) %*% (M_Psi2_nim - M_Psi1_nim)
      H[(m+1):(m+beta_dim), (m+1):(m+beta_dim)] = - M_txt1_nop %*% diag(c(1/M_h_no1^2), n.ctype[2], n.ctype[2]) %*% M_xt1_nop -
        M_tXt1_nlp %*% diag(c(M_S_nl1/(1-M_S_nl1)^2), n.ctype[3], n.ctype[3]) %*% M_Xt1_nlp -
        (M_tXt2_nip - M_tXt1_nip) %*% diag(c(M_S1_ni1 * M_S2_ni1/M_S1mS2_ni1^2), n.ctype[4], n.ctype[4]) %*% (M_Xt2_nip - M_Xt1_nip)
      H[1:m, (m+1):(m+beta_dim)] = t(H[(m+1):(m+beta_dim), 1:m])
      s_gamma_old = s_gamma
      s_df_old = s_df
      s_sigma2_old = 1 / (2 * s_gamma_old)
      pos = c((M_theta_m1 > control$min.theta), rep(TRUE, beta_dim))
      HRinv[pos, pos] = solve(-H[pos, pos] + (1/s_sigma2_old) * 
                                M_Rstar_ll[pos, pos])
      s_df = m - sum(diag(HRinv %*% M_Rstar_ll))/s_sigma2_old
      s_sigma2 = c(t(M_theta_m1) %*% M_R_mm %*% M_theta_m1/s_df)
      s_gamma = 1/(2 * s_sigma2)
      full.iter = full.iter + k
      if (full.iter > control$max.iter[3]) {
        break
      }
      if (abs(s_df - s_df_old) < control$tol_2) {
        break
      }
    }
    s_gamma  = s_gamma_old
    M_2 = -H + 2 * s_gamma * M_Rstar_ll
    HesbinvD = matrix(diag(solve(-H[(m + 1):(m + beta_dim), (m + 1):(m + beta_dim)])))
    Minv_2 = matrix(0, m + beta_dim, m + beta_dim)
    pos = c((M_theta_m1 > control$min.theta), rep(TRUE, beta_dim))
    temp = solve(M_2[pos, pos])
    Minv_2[pos, pos] = temp
    cov_NuNu_M2HM2 = Minv_2 %*% (-H) %*% Minv_2
    cov.theta = cov_NuNu_M2HM2[1:m, 1:m]
    var.beta = matrix(diag(cov_NuNu_M2HM2[(m + 1):(m + beta_dim), (m + 1):(m + beta_dim)]))
    var.beta[var.beta < 0] = HesbinvD[var.beta < 0]
    se.beta = sqrt(var.beta)
    colnames(se.beta) = c("se")
    rownames(se.beta) = paste("beta", 1:beta_dim, sep = "")
    Beta = M_beta_p1
    colnames(Beta) = c("estimate")
    rownames(Beta) = paste("beta", 1:beta_dim, sep = "")
    return(list(cvg = cvg, Beta = Beta, Theta = c(M_theta_m1) * (M_theta_m1 > control$min.theta),
                cov.theta = cov.theta, se.beta = se.beta, bins = bins,
                t_i1 = t_i1, t_i2 = t_i2, X = rbind(M_xt1_nrp, M_xt1_nop, M_xt1_nlp, M_xt1_nip, M_xt2_nip), 
                t_il = c(t_i1[ctype[, "i"]]), t_ir = c(t_i2[ctype[, "i"]]), t_r = c(t_i1[ctype[, "r"]]), 
                t_l = c(t_i1[ctype[, "l"]]), t_e = c(t_i1[ctype[, "e"]]), ctype = ctype,
                Gradbeta = round((PL_FD + t(C_M) %*% M_lambda)[(m + 1) : (m + beta_dim)], 5)))
  }

plot.baseline_hazard = function (x, upper.quantile = 0.95) 
{
  bins = x$bins
  V_x_X = bins$Alpha
  prob = upper.quantile
  cov_ThetaTheta = x$cov.theta
  M_Ppsi_Xm = basis_mpl(V_x_X, bins, which = 1)
  V_sd2.Hh0.X = diag(M_Ppsi_Xm %*% cov_ThetaTheta %*% t(M_Ppsi_Xm))
  pos.var = V_sd2.Hh0.X > 0
  V_Hh0_X = c(M_Ppsi_Xm %*% matrix(x$Theta, ncol = 1))[pos.var]
  V_x_X = V_x_X[pos.var]
  V_sd.Hh0.X = sqrt(V_sd2.Hh0.X[pos.var])
  upper = V_Hh0_X + prob * V_sd.Hh0.X
  lower = V_Hh0_X - prob * V_sd.Hh0.X
  lower[lower < 0] = 0
  plot(V_x_X, V_Hh0_X, ylim=c(min(V_Hh0_X, upper, lower), max(V_Hh0_X, upper, lower)), 
       type='l', col='red', xlab = 'Years', ylab='Baseline hazard estimate')
  lines(V_x_X, upper, lty = 2, col = 'blue')
  lines(V_x_X, lower, lty = 2, col = 'blue')
  legend("topleft", legend = c("Estimate", "95% conf. interval"), 
         lty = c(1, 2), col = c('red', 'blue'), ncol = 1, bty = "n", cex = 0.75)
}
