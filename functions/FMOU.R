KF <- function(F_t, V_t, G_t, W_t, y){
  # this function uses Kalman Filter to calculate the diagonal and primary off-diagonal term of Cov(Z_l|Y,theta)
  # a_vec is diagonal, b_vec is primary off-diagonal
  n = length(y)
  a_record <- rep(NA, n)
  R_record <- rep(NA, n)
  f_record <- rep(NA, n)
  Q_record <- rep(NA, n)
  m_record <- rep(NA, n)
  C_record <- rep(NA, n)
  
  a_record[1] = 0
  #R_record[1] = W_t
  R_record[1] = W_t/(1-G_t^2)
  f_record[1] = F_t * a_record[1]
  Q_record[1] = F_t * R_record[1] * F_t + V_t
  m_record[1] = a_record[1] + R_record[1] * F_t * (1/Q_record[1]) * (y[1]-f_record[1])
  C_record[1] = R_record[1] - R_record[1] * F_t * (1/Q_record[1]) * F_t * R_record[1]
  
  for(i in 2:n){
    a_record[i] <- G_t * m_record[(i-1)]
    R_record[i] <- G_t * C_record[(i-1)] * G_t + W_t
    f_record[i] <- F_t * a_record[i]
    Q_record[i] <- F_t * R_record[i] * F_t + V_t
    m_record[i] <- a_record[i] + R_record[i] * F_t * (1/Q_record[i]) * (y[i]-f_record[i])
    C_record[i] <- R_record[i] - R_record[i] * F_t * (1/Q_record[i]) * F_t * R_record[i]
  }
  
  S_record <- rep(NA, n)
  s_record <- rep(NA, n)
  b_record <- rep(NA, (n-1)) # primary off-diagonal terms
  s_record[n] = m_record[n]
  S_record[n] = C_record[n]
  
  for(i in (n-1):1){
    s_cur <- m_record[i] + C_record[i] * G_t * (1/R_record[(i+1)]) * (s_record[(i+1)]-a_record[(i+1)])
    S_cur <- C_record[i] - C_record[i] * G_t * (1/R_record[(i+1)]) * (R_record[(i+1)] - S_record[(i+1)]) * (1/R_record[(i+1)]) * G_t * C_record[i]
    b_cur <- C_record[i] * G_t * (1/R_record[(i+1)]) * S_record[(i+1)]
    s_record[i] <- s_cur
    S_record[i] <- S_cur
    b_record[i] <- b_cur
  }
  res <- list(#a_vec = S_record,
              b_vec = b_record,
              #a_t = a_record,
              s_t = s_record,
              f_t=f_record,
              Q_t=Q_record,
              a_t = a_record,
              S_t = S_record,
              R_t = R_record,
              m_t = m_record,
              C_t = C_record)
  #m_t = m_record) 
  return(res)
}


EM_alg_FMOU <- function(output, d, M=50, threshold=10^{-4},
                        est_U0=T, U0=NA, est_sigma0_2=T, sigma0_2=NA, 
                        U_init = NA, rho_init=NA, sigma2_init=NA,
                        track_iterations=F,
                        neg_log_lik=F, 
                        track_Q = F){
  # This function is used to implement EM alg for multi dimensional HOU
  # @output: a k*n matrix, output=(y_1,...,y_n).
  # @d: number of latent factors
  # @M: number of iterations in EM algorithm
  # @threshold: for the predictive mean of observations
  # @est_U0: estimate U0 or not
  # @U0: a k*d matrix. If est_UO=F, then U0 must NOT be NA
  # @est_sigma0_2: estimate noise level or not
  # @sigma0_2: if est_sigma0_2=F, then sigma0_2 must NOT be NA. It is a positive scalar
  # @U_init, rho_init, sigma_init: user-specified U, rho and sigma
  # @neg_log_lik: the default is FALSE - not computing negative log likelihood
  # @track_Q: the default is FALSE - not tracking the values of Q-function
  
  n = ncol(output)
  k = nrow(output)
  
  # check if d<=k
  if(d>k){
    return("Error: d can not be greater than the dimension of each observation.")
  }
  
  # check if user specified initial values are suitable
  if(!any(is.na(U_init))){
    if(nrow(U_init)!=k){
      return("Error: U_init has k rows")
    } else if(ncol(U_init)!=d){
      return("Error: U_init has d columns")
    } else if(sqrt(mean((t(U_init)%*%U_init-diag(d))^2))>10^{-8}){
      return("Error: U_init should be an k*d orthogonal matrix")
    }
  }
  
  if(!any(is.na(rho_init))){
    if(length(rho_init)!=d){
      return("Error: rho has d elements")
    }
  }
  
  if(!any(is.na(sigma2_init))){
    if(length(sigma2_init)!=d){
      return("Error: sigma2 has d elements")
    }
  }
  
  
  
  ############ start initialization ############
  if(est_U0){
    if(sum(is.na(U_init))==0){ # U_init is given
      U_pre = U_init
    } else{
      if(d<=min(n,k)){
        U_pre = svd(output)$u[,1:d]
      } else{
        U_pre = diag(nrow=k,ncol=d)
      }
    }
  } else{
    U_pre = U0
  }
  
  if(sum(is.na(rho_init))==0){ # rho_init is given
    rho_pre = rho_init
  } else{
    rho_pre = seq(0.8,0.99,length.out=d)
  }
  
  if(sum(is.na(sigma2_init))==0){ # sigma2_init is given
    sigma2_pre <- sigma2_init
  } else{
    sigma2_pre <- seq(0.5,1,length.out=d)
  }
  
  
  tr_Y_Yt <- sum(output^2)
  output_tilde = t(U_pre) %*% output # d*n matrix
  #sigma2_pre <- rep(log(1.5),d)
  #rho_pre <- rep(log(1.5),d)
  
  if(est_sigma0_2){
    sigma2_0_pre <- log(1.5)
    #sigma2_0_pre <- runif(1,1,100)
  } else{
    sigma2_0_pre = sigma0_2
  }
  
  Z_hat_pre <- matrix(NA, d, n)
  diag_Sigma = as.list(1:d)
  off_diag_Sigma = as.list(1:d)
  
  for(l in 1:d){
    kf = KF(F_t=1, V_t=sigma2_0_pre, G_t=rho_pre[l], W_t=sigma2_pre[l], output_tilde[l,])
    Z_hat_pre[l,] = kf$s_t
    diag_Sigma[[l]] = kf$S_t
    off_diag_Sigma[[l]] = kf$b_vec
  }
  
  pred_cur <- U_pre %*% Z_hat_pre
  pred_pre <- pred_cur + 1
  
  record_rho <- matrix(NA,d, M)
  record_sigma2 <- matrix(NA, d, M)
  record_sigma2_0 <- rep(NA, M)
  record_neg_lik <- rep(NA, M)
  record_Q_func <- rep(NA, M)
  
  
  ######## start EM alg
  m=1
  while(sqrt(mean((pred_cur - pred_pre)^2)) > threshold & m<=M){
    Theta_old = list(S = diag_Sigma, S_tilde=off_diag_Sigma,Z_hat=Z_hat_pre)
    pred_pre = pred_cur
    
    ## update U
    if(est_U0){
      Z_hat_Y_t <- Z_hat_pre %*% t(output)
      svd_Z_hat_Y_t <- svd(Z_hat_Y_t)
      U_cur <- svd_Z_hat_Y_t$v %*% t(svd_Z_hat_Y_t$u)
    } else{
      U_cur = U0
    }
    output_tilde = t(U_cur) %*% output
    
    ## update sigma2_0
    if(est_sigma0_2){
      tr_Z_hat_Y_t_V <- sum(Z_hat_pre*output_tilde)
      tr_Z_t_Z <- sum(Z_hat_pre*Z_hat_pre)
      sigma2_0_cur <- (tr_Y_Yt - 2*tr_Z_hat_Y_t_V + sum(unlist(diag_Sigma)) + tr_Z_t_Z)/(n*k) 
    } else{
      sigma2_0_cur <- sigma0_2
    }
    record_sigma2_0[m] = sigma2_0_cur
    
    log_Q=rep(NA,d)
    for(l in 1:d){
      ## update rho
      p3 = n*sum(Z_hat_pre[l,2:(n-1)]^2) + n*sum(diag_Sigma[[l]][2:(n-1)])-sum(Z_hat_pre[l,2:(n-1)]^2) - sum(diag_Sigma[[l]][2:(n-1)])
      p2 = 2*sum(Z_hat_pre[l,2:n]*Z_hat_pre[l,1:(n-1)]) + 2*sum(off_diag_Sigma[[l]]) - n*sum(Z_hat_pre[l,2:n]*Z_hat_pre[l,1:(n-1)]) - n*sum(off_diag_Sigma[[l]])
      p1 = -sum(Z_hat_pre[l,]^2) - sum(diag_Sigma[[l]])-n*sum(Z_hat_pre[l,2:(n-1)]^2)-n*sum(diag_Sigma[[l]][2:(n-1)])
      p0 = n*(sum(Z_hat_pre[l,2:n]*Z_hat_pre[l,1:(n-1)]) + sum(off_diag_Sigma[[l]]))
      cubic_equation_solution = polyroot(c(p0,p1,p2,p3))
      #print(cubic_equation_solution)
      my_solution_index = which(abs(Im(cubic_equation_solution))<10^{-10} & Re(cubic_equation_solution)<=1 & Re(cubic_equation_solution)>=-1)
      #print(my_solution_index)
      rho_cur = Re(cubic_equation_solution[my_solution_index])
      record_rho[l,m] = rho_cur
      
      #print(cubic_equation_solution)
      
      ## update sigma2
      sigma2_cur = ((1-rho_cur^2)*Z_hat_pre[l,1]^2 + sum((Z_hat_pre[l,2:n]-rho_cur*Z_hat_pre[l,1:(n-1)])^2) + sum(diag_Sigma[[l]]) + rho_cur^2*sum(diag_Sigma[[l]][2:(n-1)]) - 2*rho_cur*sum(off_diag_Sigma[[l]]))/n
      record_sigma2[l,m] = sigma2_cur
      
      ## update Z_hat_pre
      kf = KF(F_t=1, V_t=sigma2_0_cur, G_t=rho_cur, W_t=sigma2_cur, output_tilde[l,])
      Z_hat_pre[l,] = kf$s_t
      diag_Sigma[[l]] = kf$S_t
      off_diag_Sigma[[l]] = kf$b_vec
      log_Q[l]=sum(log(kf$Q_t)) 
    }
    
    pred_cur = U_cur %*% Z_hat_pre
    
    # record neg_lik
    if(neg_log_lik){
      S_2_here = tr_Y_Yt - sum(output_tilde*Z_hat_pre)
      neg_log_det=-1/2*sum(log_Q)+n*d/2*log(sigma2_0_cur)
      record_neg_lik[m] = -(neg_log_det-(n*k)/2*log(S_2_here))
    }
    
    # record Q_function
    if(track_Q){
      if(est_U0){
        Q_function_part1 = -(n*k/2)*log(sigma2_0_cur) + sum(svd_Z_hat_Y_t$d)/sigma2_0_cur 
      }else{
        Q_function_part1 = -(n*k/2)*log(sigma2_0_cur) + sum(diag((t(output)%*%U0%*%Z_hat_pre)))/sigma2_0_cur 
      }
      Q_function_part2 = -0.5*sum(log(record_sigma2[,m]^{n}/(1-record_rho[,m]^2))) - sum(Theta_old$Z_hat*Theta_old$Z_hat)*0.5/sigma2_0_cur - sum(unlist(Theta_old$diag_Sigma))*0.5/sigma2_0_cur-tr_Y_Yt*0.5/sigma2_0_cur
      Q_function_part3 = 0
      for(l in 1:d){
        cur_val_num = (1-record_rho[l,m]^2)*Theta_old$Z_hat[l,1] + sum((Theta_old$Z_hat[l,2:n]-record_rho[l,m]*Theta_old$Z_hat[l,1:(n-1)])^2) + sum(Theta_old$diag_Sigma[[l]]) + record_rho[l,m]^2*sum(Theta_old$diag_Sigma[[l]][2:(n-1)])-2*record_rho[l,m]*sum(Theta_old$off_diag_Sigma[[l]])
        Q_function_part3 = Q_function_part3 - cur_val_num/(2*record_sigma2[l,m])
      }
      record_Q_func[m] = Q_function_part1 +Q_function_part2 + Q_function_part3
    }

    m = m+1
  }
  
  a_t = matrix(NA, d,n)
  S_t = matrix(NA, d, n)
  output_tilde = t(U_cur) %*% output
  for(l in 1:d){
    kf = KF(F_t=1, V_t=record_sigma2_0[(m-1)], G_t=record_rho[l,(m-1)], 
            W_t=record_sigma2[l, (m-1)], output_tilde[l,])
    a_t[l,] = kf$a_t
    S_t[l,] = kf$S_t
  }
  pred_mean_var = matrix(NA,dim(U_cur)[1],n)
  for(tt in 1:n){
    pred_mean_var[,tt] = rowSums(U_cur*t(matrix(S_t[,tt],d,dim(U_cur)[1]))*U_cur)
  }
  

  pred_mean_lower = U_cur %*% Z_hat_pre - 1.96 * sqrt(pred_mean_var)
  pred_mean_upper = U_cur %*% Z_hat_pre + 1.96 * sqrt(pred_mean_var)
  
  if(track_iterations){
    res <- list(Z_hat = Z_hat_pre,
                record_sigma2=record_sigma2[,1:(m-1)],
                record_rho=record_rho[,1:(m-1)],
                record_sigma2_0 = record_sigma2_0[1:(m-1)],
                pred_mean = pred_cur,
                pred_mean_95lower = pred_mean_lower,
                pred_mean_95upper = pred_mean_upper,
                U = U_cur,
                a_t = a_t,
                num_iterations = m-1,
                record_neg_lik=record_neg_lik,
                #posterior_var = S_t,
                record_Q_func=record_Q_func,
                Z_post_sd = sqrt(S_t))
  } else{
    res <- list(Z_hat = Z_hat_pre,
                sigma2=record_sigma2[,(m-1)],
                rho=record_rho[,(m-1)],
                sigma2_0 = record_sigma2_0[(m-1)],
                pred_mean = pred_cur,
                pred_mean_95lower = pred_mean_lower,
                pred_mean_95upper = pred_mean_upper,
                U = U_cur,
                a_t = a_t,
                num_iterations = m-1,
                record_neg_lik=record_neg_lik[m-1],
                #posterior_var = S_t,
                record_Q_func=record_Q_func,
                Z_post_sd = sqrt(S_t))
  }
  
  
  return(res)
} 





