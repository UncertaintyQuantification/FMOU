KF_multidim <- function(Mt, Qt, Phi, Rt,m0, C0, output){
  
  n = ncol(output) 
  k = nrow(output)
  d = ncol(Mt)
  
  ## Kalman Filter
  a_record <- matrix(NA, nrow=d, ncol=n+1) 
  R_record <- as.list(1:(n+1)) 
  m_record <- matrix(NA, nrow=d, ncol=n+1) 
  C_record <- as.list(1:(n+1)) 
  
  a_record[,1] = NA
  R_record[[1]] = NA
  m_record[,1] = m0
  C_record[[1]] = C0
  
  for(t in 1:n){
    a_record[,(t+1)] = Phi %*% m_record[,t]
    R_record[[t+1]] = Phi %*% C_record[[(t)]] %*% t(Phi) + Qt
    Kt = R_record[[t+1]] %*% t(Mt) %*% solve(Mt %*% R_record[[t+1]] %*% t(Mt) + Rt)
    m_record[,(t+1)] = a_record[,(t+1)] + Kt %*% (output[,t]-Mt %*% a_record[,(t+1)])
    C_record[[t+1]] = R_record[[t+1]] - Kt %*% Mt %*% R_record[[t+1]]
  }
  
  # backward smoothing
  S_record <- as.list(1:(n+1)) # Var(v_t|y_{1:n})
  s_record <- matrix(NA, nrow=d, ncol=(n+1)) # E(v_t|y_{1:n})
  tilde_S_record <-as.list(1:n) # Cov(v_{t}, v_{t-1}|y_{1:n}), t=1,...,n
  
  S_record[[n+1]] = C_record[[n+1]]
  s_record[,(n+1)] = m_record[,(n+1)]
  #tilde_S_record[[n]] = (diag(d) - Kt%*%Mt) %*% Phi %*% C_record[[n-1]]
  
  for(t in n:1){
    # pre compute
    J_t = C_record[[t]] %*% t(Phi) %*% solve(R_record[[(t+1)]])
    
    s_record[,t] = m_record[,t] + J_t %*% (s_record[,t+1]-a_record[,t+1])
    S_record[[t]] = C_record[[t]] - J_t %*% (R_record[[t+1]]-S_record[[t+1]]) %*% t(J_t)
    tilde_S_record[[t]] = S_record[[t+1]] %*% t(J_t)
  }
  
  
  res <- list(
    m=m_record,
    C=C_record,
    s=s_record,
    S=S_record,
    tilde_S=tilde_S_record)
  return(res)
}

EM_SS <- function(output, Mt, n_iter=10){
  
  k = nrow(output)
  n = ncol(output)
  d = ncol(Mt)
  
  ### initialize the parameters
  mu = as.vector(rep(0, d))
  Sigma = diag(0.001,d,d)
  Q = diag(0.001,d,d)
  R = diag(0.001,k,k)
  Phi = diag(d)
  
  ### initialize KF and RTS smoother
  kf = KF_multidim(Mt, Q, Phi, R,mu, Sigma, output)
  
  A = Reduce("+", kf$S) - kf$S[[(n+1)]] + kf$s[,1:n] %*% t(kf$s[,1:n])
  B = Reduce("+", kf$tilde_S) + kf$s[,2:(n+1)] %*% t(kf$s[,1:n])
  C = Reduce("+", kf$S) - kf$S[[1]] + kf$s[,2:(n+1)] %*% t(kf$s[,2:(n+1)])
  
  m=1
  while(m<=n_iter){
    Phi = B %*% solve(A)
    Q = (C-B%*%solve(A)%*% t(B))/n
    Q = (Q + t(Q))/2
    R = diag(0, k,k)
    for(t in 1:n){
      R = R + (output[,t]-Mt%*%kf$s[,(t+1)]) %*% t(output[,t]-Mt%*%kf$s[,(t+1)]) + Mt %*% kf$S[[t+1]] %*% t(Mt)
    }
    R = R/n
    R = diag(diag(R), k,k)
    mu = as.vector(kf$s[,1])
    Sigma = kf$S[[1]]
    
    kf = KF_multidim(Mt, Q, Phi, R,mu, Sigma, output)
    
    A = Reduce("+", kf$S) - kf$S[[(n+1)]] + kf$s[,1:n] %*% t(kf$s[,1:n])
    B = Reduce("+", kf$tilde_S) + kf$s[,2:(n+1)] %*% t(kf$s[,1:n])
    C = Reduce("+", kf$S) - kf$S[[1]] + kf$s[,2:(n+1)] %*% t(kf$s[,2:(n+1)])
    
    m = m + 1
  }
  
  kf = KF_multidim(Mt, Q, Phi, R,mu, Sigma, output)
  pred_mean = Mt %*% kf$s[,2:(n+1)]
  
  res = list(Phi=Phi,
             Q=Q,
             R=R,
             mu=mu,
             Sigma=Sigma,
             pred_mean=pred_mean)
  return(res)
}

