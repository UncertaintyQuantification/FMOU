## This is the codes for Experiement 2: Correctly specified data-generating processes with estimated factor loading matrices

library(rstiefel)
library(rospca)
library(ggplot2)
library(gridExtra)

source('../functions/DMD.R') 
source('../functions/FMOU.R')

d = 5 
k = 20 # number of stations
n1 = 100 # 100 timesteps
n2 = 200 # 200 timesteps
n3 = 400 # 400 timesteps
sigma_0_2_list = c(1, 2)
N = 20 # number of replications 

######################################## n=100 ########################################
# largest principal angel between U0 and estimated U0, true d is given
angle_diff_FMOU_n1 = matrix(NA, N, length(sigma_0_2_list)) 
angle_diff_QW1_n1 = matrix(NA, N, length(sigma_0_2_list)) 
angle_diff_QW5_n1 = matrix(NA, N, length(sigma_0_2_list)) 
angle_diff_DMD_n1 = matrix(NA, N, length(sigma_0_2_list)) 

# largest principal angel between U0 and estimated U0, d is estimated
angle_diff_FMOU_n1_est_d = matrix(NA, N, length(sigma_0_2_list)) 
angle_diff_QW1_n1_est_d = matrix(NA, N, length(sigma_0_2_list))
angle_diff_QW5_n1_est_d = matrix(NA, N, length(sigma_0_2_list))
angle_diff_DMD_n1_est_d = matrix(NA, N, length(sigma_0_2_list))

# RMSE of mean of output, true d is given
rmse_mean_FMOU_n1 = matrix(NA, N, length(sigma_0_2_list)) 
rmse_mean_QW1_n1 = matrix(NA, N, length(sigma_0_2_list)) 
rmse_mean_QW5_n1 = matrix(NA, N, length(sigma_0_2_list)) 
rmse_mean_DMD_n1 = matrix(NA, N, length(sigma_0_2_list)) 

# RMSE of mean of output, d is estimated
rmse_mean_FMOU_n1_est_d = matrix(NA, N, length(sigma_0_2_list)) 
rmse_mean_QW1_n1_est_d = matrix(NA, N, length(sigma_0_2_list))
rmse_mean_QW5_n1_est_d = matrix(NA, N, length(sigma_0_2_list))
rmse_mean_DMD_n1_est_d = matrix(NA, N, length(sigma_0_2_list))

# 95% CI of predictive mean by FMOU, d is given
pred_mean_FMOU_lower_n1 = as.list(1:length(sigma_0_2_list))
pred_mean_FMOU_upper_n1 = as.list(1:length(sigma_0_2_list))

# 95% CI of predictive mean by FMOU, d is estimated
pred_mean_FMOU_lower_n1_est_d = as.list(1:length(sigma_0_2_list))
pred_mean_FMOU_upper_n1_est_d = as.list(1:length(sigma_0_2_list))

# length of 95% CI
L_95_n1 = matrix(NA, nrow=N, ncol=length(sigma_0_2_list)) # d is given
L_95_n1_est_d = matrix(NA, nrow=N, ncol=length(sigma_0_2_list)) # d is estimated

# coverage of 95% CI
P_95_n1 = matrix(NA, nrow=N, ncol=length(sigma_0_2_list)) # d is given
P_95_n1_est_d = matrix(NA, nrow=N, ncol=length(sigma_0_2_list)) # d is estimated

# estimated number of d by FMOU, QW1, QW5, DMD
est_d_FMOU_n1 = matrix(NA, N, length(sigma_0_2_list))
est_d_QW1_n1 = matrix(NA, N, length(sigma_0_2_list))
est_d_QW5_n1 = matrix(NA, N, length(sigma_0_2_list))
est_d_DMD_n1 = matrix(NA, N, length(sigma_0_2_list))

# record signals, observations and predictions from the first replication
example_obs_n1 = as.list(1:length(sigma_0_2_list))
example_mean_n1 = as.list(1:length(sigma_0_2_list))
example_pred_mean_FMOU_n1 = as.list(1:length(sigma_0_2_list))
example_pred_mean_QW1_n1 = as.list(1:length(sigma_0_2_list))
example_pred_mean_QW5_n1 = as.list(1:length(sigma_0_2_list))
example_pred_mean_DMD_n1 = as.list(1:length(sigma_0_2_list))

# computational time of FMOU
time_record_FMOU_n1 = matrix(NA, N, length(sigma_0_2_list))

# start simulation
for(idx_noise in 1: length(sigma_0_2_list)){
  # specify noise level
  noise_level= sigma_0_2_list[idx_noise]
  
  for(it in 1:N){
    set.seed(it)
    if(it%%10 == 0){print(it)}
  ## generate latent processes, factor loading matrix
    U = rustiefel(k, k)
    z = matrix(NA, d, n1)
    sigma_2 = runif(d, 0.5, 1)
    rho = runif(d, 0.95, 1)
    for(l in 1:d){
      R = matrix(NA, n1, n1)
      diag(R) = 1
      for(ir in 1:n1){
        for(ic in 1:n1){
          R[ir, ic] = rho[l]^(abs(ir-ic)) * R[ir, ir]
        }
      }
      R = (sigma_2[l]/(1-rho[l]^2) )* R
      z[l, ] = t(chol(R)) %*% rnorm(n1)
    }
  ## generate signal and observation
    signal = U[,1:d] %*% z
    y = signal + matrix(rnorm(n1*k,mean=0,sd=sqrt(noise_level)),k,n1)
    if(it==1){
      example_obs_n1[[idx_noise]] = y
      example_mean_n1[[idx_noise]] = signal
    }
  
  ## 1. FMOU
  ## select d
    svd_output=svd(y)
    FMOU_loss_score = NULL
    criteria_val_pre = log(mean((y - svd_output$u[,1]%*%t(svd_output$u[,1])%*%y)^2)) + (k+n1)/(k*n1)*log(k*n1/(k+n1))
    FMOU_loss_score = c(FMOU_loss_score,  criteria_val_pre)
    for(i_d in 2:ceiling(dim(y)[1]*2/3)){
      criteria_val_cur = log(mean((y - svd_output$u[,1:i_d]%*%t(svd_output$u[,1:i_d])%*%y)^2)) + i_d*(k+n1)/(k*n1)*log(k*n1/(k+n1))
      FMOU_loss_score = c(FMOU_loss_score,  criteria_val_cur)
    }
    est_d_FMOU_here = which.min(FMOU_loss_score) # i_d-1
    est_d_FMOU_n1[it, idx_noise] = which.min(FMOU_loss_score) # i_d-1
  ## fit FMOU
    time1 = Sys.time()
    tilde_y = t(U[, 1:d])%*%y
    em_alg <- EM_alg_FMOU(output=y, d=d, M=1000, threshold=10^{-6},
                                     est_U0=T, est_sigma0_2=T)
    U_est_FMOU = em_alg$U
    time2 = Sys.time()
    
    tilde_y = t(U[, 1:est_d_FMOU_here])%*%y
    em_alg_est_d <- EM_alg_FMOU(output=y, d=est_d_FMOU_here, M=1000, threshold=10^{-6},
                                           est_U0=T, est_sigma0_2=T)
    U_est_FMOU_est_d = em_alg_est_d$U
    
    time_record_FMOU_n1[it, idx_noise] = time2 - time1
    
    fit_latent_state_EM = U_est_FMOU %*% em_alg$Z_hat
    fit_latent_state_EM_est_d = U_est_FMOU_est_d %*% em_alg_est_d$Z_hat
    
    angle_diff_FMOU_n1[it, idx_noise] = rospca::angle(U[,1:d], U_est_FMOU)
    angle_diff_FMOU_n1_est_d[it, idx_noise] = rospca::angle(U[,1:d], U_est_FMOU_est_d)
    rmse_mean_FMOU_n1[it, idx_noise] = sqrt(mean(( fit_latent_state_EM - signal)^2))
    rmse_mean_FMOU_n1_est_d[it, idx_noise] = sqrt(mean(( fit_latent_state_EM_est_d - signal)^2))
   
    L_95_n1[it,idx_noise] = mean(em_alg$pred_mean_95upper-em_alg$pred_mean_95lower)
    L_95_n1_est_d[it,idx_noise] = mean(em_alg_est_d$pred_mean_95upper-em_alg_est_d$pred_mean_95lower)
    
    P_95_n1[it, idx_noise] = mean(em_alg$pred_mean_95upper>signal &em_alg$pred_mean_95lower< signal)
    P_95_n1_est_d[it, idx_noise] = mean(em_alg_est_d$pred_mean_95upper>signal &em_alg_est_d$pred_mean_95lower< signal)
    
      
    if(it==1){
      example_pred_mean_FMOU_n1[[idx_noise]] = list(with_true_d = fit_latent_state_EM,
                                                    with_est_d =fit_latent_state_EM_est_d)
      pred_mean_FMOU_lower_n1[[idx_noise]] = em_alg$pred_mean_95lower
      pred_mean_FMOU_upper_n1[[idx_noise]] = em_alg$pred_mean_95upper
      pred_mean_FMOU_lower_n1_est_d[[idx_noise]] = em_alg_est_d$pred_mean_95lower
      pred_mean_FMOU_upper_n1_est_d[[idx_noise]] = em_alg_est_d$pred_mean_95upper
    }
    
    ## 2. DMD
    DMD_fit = DMD_alg(y, fix_r=T, r=d)
    DMD_fit_est_d = DMD_alg(y, threshold=0.99)
    rmse_mean_DMD_n1[it, idx_noise] = sqrt(mean((cbind(y[,1], DMD_fit$in_sample_pred) - signal)^2))
    rmse_mean_DMD_n1_est_d[it, idx_noise] = sqrt(mean((cbind(y[,1], DMD_fit_est_d$in_sample_pred) - signal)^2))
    if(it==1){
      example_pred_mean_DMD_n1[[idx_noise]] = list(with_true_d = cbind(y[,1], DMD_fit$in_sample_pred),
                                                   with_est_d = cbind(y[,1], DMD_fit_est_d$in_sample_pred)
      )
    }
    
    angle_diff_DMD_n1[it, idx_noise] = rospca::angle(U[,1:d], DMD_fit$A_hat[,1:DMD_fit$rank])
    angle_diff_DMD_n1_est_d[it, idx_noise] = rospca::angle(U[,1:d], DMD_fit_est_d$A_hat[,1:DMD_fit_est_d$rank])
    est_d_DMD_n1[it, idx_noise] = DMD_fit_est_d$rank
    
    ## 3. LY1
    cov_lag_1 = matrix(0, k, k)
    for(lag in 1:1){
      cov_lag_1_here = matrix(NA, k, k)
      for (i in 1:k) {
        for (j in 1:k) {
          cov_lag_1_here[i, j] <- cov(y[i , 1: (ncol(y)-lag)], y[j, (1+lag):(ncol(y))])
        }
      }
      cov_lag_1 = cov_lag_1 + cov_lag_1_here %*% t(cov_lag_1_here)
    }
    eigen_QW_1 = eigen(cov_lag_1)
    lambda_hat = eigen_QW_1$values
    eigen_ratio = lambda_hat[2:k] / lambda_hat[1:(k-1)]
    est_d_QW1_here = which.min(eigen_ratio)
    est_d_QW1_n1[it, idx_noise] = est_d_QW1_here
    U_est_QW1 = eigen_QW_1$vectors[,1:d]
    U_est_QW1_est_d = eigen_QW_1$vectors[,1:est_d_QW1_here]
    rmse_mean_QW1_n1[it, idx_noise] = sqrt(mean(( U_est_QW1%*%t(U_est_QW1)%*%y - signal)^2))
    rmse_mean_QW1_n1_est_d[it, idx_noise] = sqrt(mean(( U_est_QW1_est_d%*%t(U_est_QW1_est_d)%*%y - signal)^2))
    if(it==1){
      example_pred_mean_QW1_n1[[idx_noise]] = list(with_true_d = U_est_QW1%*%t(U_est_QW1)%*%y,
                                                   with_est_d = U_est_QW1_est_d%*%t(U_est_QW1_est_d)%*%y)
    }
    angle_diff_QW1_n1[it, idx_noise] = rospca::angle(U[,1:d], U_est_QW1)
    angle_diff_QW1_n1_est_d[it, idx_noise] = rospca::angle(U[,1:d], U_est_QW1_est_d)
    
    ## 4. LY5
    cov_lag_5 = matrix(0, k, k)
    for(lag in 1:5){
      cov_lag_5_here = matrix(NA, k, k)
      for (i in 1:k) {
        for (j in 1:k) {
          cov_lag_5_here[i, j] <- cov(y[i , 1: (ncol(y)-lag)], y[j, (1+lag):(ncol(y))])
        }
      }
      cov_lag_5 = cov_lag_5 + cov_lag_5_here %*% t(cov_lag_5_here)
    }
    eigen_QW_5 = eigen(cov_lag_5)
    lambda_hat = eigen_QW_5$values
    eigen_ratio = lambda_hat[2:k] / lambda_hat[1:(k-1)]
    est_d_QW5_here = which.min(eigen_ratio)
    est_d_QW5_n1[it, idx_noise] = est_d_QW5_here
    U_est_QW5 = eigen_QW_5$vectors[,1:d]
    U_est_QW5_est_d = eigen_QW_5$vectors[,1:est_d_QW5_here]
    rmse_mean_QW5_n1[it, idx_noise] = sqrt(mean(( U_est_QW5%*%t(U_est_QW5)%*%y - signal)^2))
    rmse_mean_QW5_n1_est_d[it, idx_noise] = sqrt(mean(( U_est_QW5_est_d%*%t(U_est_QW5_est_d)%*%y - signal)^2))
    if(it==1){
      example_pred_mean_QW5_n1[[idx_noise]] = list(with_true_d = U_est_QW5%*%t(U_est_QW5)%*%y,
                                                   with_est_d = U_est_QW5_est_d%*%t(U_est_QW5_est_d)%*%y
      )
    }
    
    angle_diff_QW5_n1[it, idx_noise] = rospca::angle(U[,1:d], U_est_QW5)
    angle_diff_QW5_n1_est_d[it, idx_noise] = rospca::angle(U[,1:d], U_est_QW5_est_d)

  }
}

angle_diff_n1_noise_1 = cbind(angle_diff_FMOU_n1[,1], angle_diff_DMD_n1[,1], angle_diff_QW1_n1[,1], angle_diff_QW5_n1[,1])
angle_diff_n1_noise_2 = cbind(angle_diff_FMOU_n1[,2], angle_diff_DMD_n1[,2], angle_diff_QW1_n1[,2], angle_diff_QW5_n1[,2])
angle_diff_n1_est_d_noise_1 = cbind(angle_diff_FMOU_n1_est_d[,1], angle_diff_DMD_n1_est_d[,1], angle_diff_QW1_n1_est_d[,1], angle_diff_QW5_n1_est_d[,1])
angle_diff_n1_est_d_noise_2 = cbind(angle_diff_FMOU_n1_est_d[,2], angle_diff_DMD_n1_est_d[,2], angle_diff_QW1_n1_est_d[,2], angle_diff_QW5_n1_est_d[,2])
rmse_mean_n1_noise_1 = cbind(rmse_mean_FMOU_n1[,1], rmse_mean_DMD_n1[,1], rmse_mean_QW1_n1[,1], rmse_mean_QW5_n1[,1])
rmse_mean_n1_noise_2 = cbind(rmse_mean_FMOU_n1[,2], rmse_mean_DMD_n1[,2], rmse_mean_QW1_n1[,2], rmse_mean_QW5_n1[,2])
rmse_mean_n1_est_d_noise_1 = cbind(rmse_mean_FMOU_n1_est_d[,1], rmse_mean_DMD_n1_est_d[,1], rmse_mean_QW1_n1_est_d[,1], rmse_mean_QW5_n1_est_d[,1])
rmse_mean_n1_est_d_noise_2 = cbind(rmse_mean_FMOU_n1_est_d[,2], rmse_mean_DMD_n1_est_d[,2], rmse_mean_QW1_n1_est_d[,2], rmse_mean_QW5_n1_est_d[,2])

est_d_all_result_n1_noise_1 = cbind(est_d_FMOU_n1[,1], est_d_DMD_n1[,1], est_d_QW1_n1[,1], est_d_QW5_n1[,1])
est_d_all_result_n1_noise_2 = cbind(est_d_FMOU_n1[,2], est_d_DMD_n1[,2], est_d_QW1_n1[,2], est_d_QW5_n1[,2])


######################################## n=200 ########################################
angle_diff_FMOU_n2 = matrix(NA, N, length(sigma_0_2_list)) 
angle_diff_QW1_n2 = matrix(NA, N, length(sigma_0_2_list)) 
angle_diff_QW5_n2 = matrix(NA, N, length(sigma_0_2_list)) 
angle_diff_DMD_n2 = matrix(NA, N, length(sigma_0_2_list)) 

angle_diff_FMOU_n2_est_d = matrix(NA, N, length(sigma_0_2_list)) 
angle_diff_QW1_n2_est_d = matrix(NA, N, length(sigma_0_2_list))
angle_diff_QW5_n2_est_d = matrix(NA, N, length(sigma_0_2_list))
angle_diff_DMD_n2_est_d = matrix(NA, N, length(sigma_0_2_list))

rmse_mean_FMOU_n2 = matrix(NA, N, length(sigma_0_2_list)) 
rmse_mean_QW1_n2 = matrix(NA, N, length(sigma_0_2_list)) 
rmse_mean_QW5_n2 = matrix(NA, N, length(sigma_0_2_list)) 
rmse_mean_DMD_n2 = matrix(NA, N, length(sigma_0_2_list)) 

rmse_mean_FMOU_n2_est_d = matrix(NA, N, length(sigma_0_2_list)) 
rmse_mean_QW1_n2_est_d = matrix(NA, N, length(sigma_0_2_list))
rmse_mean_QW5_n2_est_d = matrix(NA, N, length(sigma_0_2_list))
rmse_mean_DMD_n2_est_d = matrix(NA, N, length(sigma_0_2_list))

L_95_n2 = matrix(NA, nrow=N, ncol=length(sigma_0_2_list)) # d is given
L_95_n2_est_d = matrix(NA, nrow=N, ncol=length(sigma_0_2_list)) # d is estimated

P_95_n2 = matrix(NA, nrow=N, ncol=length(sigma_0_2_list)) # d is given
P_95_n2_est_d = matrix(NA, nrow=N, ncol=length(sigma_0_2_list)) # d is estimated

est_d_FMOU_n2 = matrix(NA, N, length(sigma_0_2_list))
est_d_QW1_n2 = matrix(NA, N, length(sigma_0_2_list))
est_d_QW5_n2 = matrix(NA, N, length(sigma_0_2_list))
est_d_DMD_n2 = matrix(NA, N, length(sigma_0_2_list))

example_obs_n2 = as.list(1:length(sigma_0_2_list))
example_mean_n2 = as.list(1:length(sigma_0_2_list))
example_pred_mean_FMOU_n2 = as.list(1:length(sigma_0_2_list))
example_pred_mean_QW1_n2 = as.list(1:length(sigma_0_2_list))
example_pred_mean_QW5_n2 = as.list(1:length(sigma_0_2_list))
example_pred_mean_DMD_n2 = as.list(1:length(sigma_0_2_list))
example_pred_mean_uncertainty_FMOU_n2 = as.list(1:length(sigma_0_2_list))

time_record_FMOU_n2 = matrix(NA, N, length(sigma_0_2_list)) 

for(idx_noise in 1: length(sigma_0_2_list)){
  # specify noise level
  noise_level= sigma_0_2_list[idx_noise]
  for(it in 1:N){
    set.seed(2*it)
    if(it%%10 == 0){print(it)}
    ## generate latent processes and factor loading matrix
    U = rustiefel(k, k)
    z = matrix(NA, d, n2)
    sigma_2 = runif(d, 0.5, 1)
    rho = runif(d, 0.95, 1)
    for(l in 1:d){
      R = matrix(NA, n2, n2)
      diag(R) = 1
      for(ir in 1:n2){
        for(ic in 1:n2){
          R[ir, ic] = rho[l]^(abs(ir-ic)) * R[ir, ir]
        }
      }
      R = (sigma_2[l]/(1-rho[l]^2) )* R
      z[l, ] = t(chol(R)) %*% rnorm(n2)
    }
    
    ## generate signal and observations
    signal = U[,1:d] %*% z
    y = signal + matrix(rnorm(n2*k,mean=0,sd=sqrt(noise_level)),k,n2)
    if(it==1){
      example_obs_n2[[idx_noise]] = y
      example_mean_n2[[idx_noise]] = signal
    }
    
    ## 1. FMOU
    ## select d 
    svd_output=svd(y)
    FMOU_loss_score = NULL
    criteria_val_pre = log(mean((y - svd_output$u[,1]%*%t(svd_output$u[,1])%*%y)^2)) + (k+n2)/(k*n2)*log(k*n2/(k+n2))
    FMOU_loss_score = c(FMOU_loss_score,  criteria_val_pre)
    
    for(i_d in 2:ceiling(dim(y)[1]*2/3)){
      criteria_val_cur = log(mean((y - svd_output$u[,1:i_d]%*%t(svd_output$u[,1:i_d])%*%y)^2)) + i_d*(k+n2)/(k*n2)*log(k*n2/(k+n2))
      FMOU_loss_score = c(FMOU_loss_score,  criteria_val_cur)
    }
    est_d_FMOU_here = which.min(FMOU_loss_score) # i_d-1
    est_d_FMOU_n2[it, idx_noise] = which.min(FMOU_loss_score) # i_d-1
    
    ## fit FMOU
    time1 = Sys.time()
    tilde_y = t(U[, 1:d])%*%y
    em_alg <- EM_alg_FMOU(output=y, d=d, M=1000, threshold=10^{-6},
                                     est_U0=T, U0=NA, est_sigma0_2=T)
    U_est_FMOU = em_alg$U
    time2 = Sys.time()
    
    tilde_y = t(U[, 1:est_d_FMOU_here])%*%y
    em_alg_est_d <- EM_alg_FMOU(output=y, d=est_d_FMOU_here, M=1000, threshold=10^{-6},
                                           est_U0=T, U0=NA, est_sigma0_2=T)
    U_est_FMOU_est_d = em_alg_est_d$U
    
    time_record_FMOU_n2[it, idx_noise] = time2 - time1
    
    fit_latent_state_EM = U_est_FMOU %*% em_alg$Z_hat
    fit_latent_state_EM_est_d = U_est_FMOU_est_d %*% em_alg_est_d$Z_hat
    
    angle_diff_FMOU_n2[it, idx_noise] = rospca::angle(U[,1:d], U_est_FMOU)
    angle_diff_FMOU_n2_est_d[it, idx_noise] = rospca::angle(U[,1:d], U_est_FMOU_est_d)
    rmse_mean_FMOU_n2[it, idx_noise] = sqrt(mean(( fit_latent_state_EM - signal)^2))
    rmse_mean_FMOU_n2_est_d[it, idx_noise] = sqrt(mean(( fit_latent_state_EM_est_d - signal)^2))
    
    L_95_n2[it,idx_noise] = mean(em_alg$pred_mean_95upper-em_alg$pred_mean_95lower)
    L_95_n2_est_d[it,idx_noise] = mean(em_alg_est_d$pred_mean_95upper-em_alg_est_d$pred_mean_95lower)
    
    P_95_n2[it, idx_noise] = mean(em_alg$pred_mean_95upper>signal &em_alg$pred_mean_95lower< signal)
    P_95_n2_est_d[it, idx_noise] = mean(em_alg_est_d$pred_mean_95upper>signal &em_alg_est_d$pred_mean_95lower< signal)
    
    if(it==1){
      example_pred_mean_FMOU_n2[[idx_noise]] = list(with_true_d = fit_latent_state_EM,
                                                    with_est_d = fit_latent_state_EM_est_d)
    }
    
    ## 2. DMD
    DMD_fit = DMD_alg(y, fix_r=T, r=d)
    DMD_fit_est_d = DMD_alg(y, threshold=0.99)
    rmse_mean_DMD_n2[it, idx_noise] = sqrt(mean((cbind(y[,1], DMD_fit$in_sample_pred) - signal)^2))
    rmse_mean_DMD_n2_est_d[it, idx_noise] = sqrt(mean((cbind(y[,1], DMD_fit_est_d$in_sample_pred) - signal)^2))
    if(it==1){
      example_pred_mean_DMD_n2[[idx_noise]] = list(with_true_d = cbind(y[,1], DMD_fit$in_sample_pred),
                                                   with_est_d = cbind(y[,1], DMD_fit_est_d$in_sample_pred)
      )
    }
    
    angle_diff_DMD_n2[it, idx_noise] = rospca::angle(U[,1:d], DMD_fit$A_hat)
    angle_diff_DMD_n2_est_d[it, idx_noise] = rospca::angle(U[,1:d], DMD_fit_est_d$A_hat)
    est_d_DMD_n2[it, idx_noise] = DMD_fit_est_d$rank
    
    ## 3. LY1
    cov_lag_1 = matrix(0, k, k)
    for(lag in 1:1){
      cov_lag_1_here = matrix(NA, k, k)
      for (i in 1:k) {
        for (j in 1:k) {
          cov_lag_1_here[i, j] <- cov(y[i , 1: (ncol(y)-lag)], y[j, (1+lag):(ncol(y))])
        }
      }
      cov_lag_1 = cov_lag_1 + cov_lag_1_here %*% t(cov_lag_1_here)
    }
    eigen_QW_1 = eigen(cov_lag_1)
    lambda_hat = eigen_QW_1$values
    eigen_ratio = lambda_hat[2:k] / lambda_hat[1:(k-1)]
    est_d_QW1_here = which.min(eigen_ratio)
    est_d_QW1_n2[it, idx_noise] = est_d_QW1_here
    U_est_QW1 = eigen_QW_1$vectors[,1:d]
    U_est_QW1_est_d = eigen_QW_1$vectors[,1:est_d_QW1_here]
    rmse_mean_QW1_n2[it, idx_noise] = sqrt(mean(( U_est_QW1%*%t(U_est_QW1)%*%y - signal)^2))
    rmse_mean_QW1_n2_est_d[it, idx_noise] = sqrt(mean(( U_est_QW1_est_d%*%t(U_est_QW1_est_d)%*%y - signal)^2))
    if(it==1){
      example_pred_mean_QW1_n2[[idx_noise]] = list(with_true_d = U_est_QW1%*%t(U_est_QW1)%*%y,
                                                   with_est_d = U_est_QW1_est_d%*%t(U_est_QW1_est_d)%*%y
      )
    }
    
    angle_diff_QW1_n2[it, idx_noise] = rospca::angle(U[,1:d], U_est_QW1)
    angle_diff_QW1_n2_est_d[it, idx_noise] = rospca::angle(U[,1:d], U_est_QW1_est_d)
    
    
    ## 4. LY5
    cov_lag_5 = matrix(0, k, k)
    for(lag in 1:5){
      cov_lag_5_here = matrix(NA, k, k)
      for (i in 1:k) {
        for (j in 1:k) {
          cov_lag_5_here[i, j] <- cov(y[i , 1: (ncol(y)-lag)], y[j, (1+lag):(ncol(y))])
        }
      }
      cov_lag_5 = cov_lag_5 + cov_lag_5_here %*% t(cov_lag_5_here)
    }
    eigen_QW_5 = eigen(cov_lag_5)
    lambda_hat = eigen_QW_5$values
    eigen_ratio = lambda_hat[2:k] / lambda_hat[1:(k-1)]
    est_d_QW5_here = which.min(eigen_ratio)
    est_d_QW5_n2[it, idx_noise] = est_d_QW5_here
    U_est_QW5 = eigen_QW_5$vectors[,1:d]
    U_est_QW5_est_d = eigen_QW_5$vectors[,1:est_d_QW5_here]
    rmse_mean_QW5_n2[it, idx_noise] = sqrt(mean(( U_est_QW5%*%t(U_est_QW5)%*%y - signal)^2))
    rmse_mean_QW5_n2_est_d[it, idx_noise] = sqrt(mean(( U_est_QW5_est_d%*%t(U_est_QW5_est_d)%*%y - signal)^2))
    if(it==1){
      example_pred_mean_QW5_n2[[idx_noise]] = list(with_true_d = U_est_QW5%*%t(U_est_QW5)%*%y,
                                                   with_est_d = U_est_QW5_est_d%*%t(U_est_QW5_est_d)%*%y
      )
    }
    
    angle_diff_QW5_n2[it, idx_noise] = rospca::angle(U[,1:d], U_est_QW5)
    angle_diff_QW5_n2_est_d[it, idx_noise] = rospca::angle(U[,1:d], U_est_QW5_est_d)

  }
}

angle_diff_n2_noise_1 = cbind(angle_diff_FMOU_n2[,1], angle_diff_DMD_n2[,1], angle_diff_QW1_n2[,1], angle_diff_QW5_n2[,1])
angle_diff_n2_noise_2 = cbind(angle_diff_FMOU_n2[,2], angle_diff_DMD_n2[,2], angle_diff_QW1_n2[,2], angle_diff_QW5_n2[,2])
angle_diff_n2_est_d_noise_1 = cbind(angle_diff_FMOU_n2_est_d[,1], angle_diff_DMD_n2_est_d[,1], angle_diff_QW1_n2_est_d[,1], angle_diff_QW5_n2_est_d[,1])
angle_diff_n2_est_d_noise_2 = cbind(angle_diff_FMOU_n2_est_d[,2], angle_diff_DMD_n2_est_d[,2], angle_diff_QW1_n2_est_d[,2], angle_diff_QW5_n2_est_d[,2])
rmse_mean_n2_noise_1 = cbind(rmse_mean_FMOU_n2[,1], rmse_mean_DMD_n2[,1], rmse_mean_QW1_n2[,1], rmse_mean_QW5_n2[,1])
rmse_mean_n2_noise_2 = cbind(rmse_mean_FMOU_n2[,2], rmse_mean_DMD_n2[,2], rmse_mean_QW1_n2[,2], rmse_mean_QW5_n2[,2])
rmse_mean_n2_est_d_noise_1 = cbind(rmse_mean_FMOU_n2_est_d[,1], rmse_mean_DMD_n2_est_d[,1], rmse_mean_QW1_n2_est_d[,1], rmse_mean_QW5_n2_est_d[,1])
rmse_mean_n2_est_d_noise_2 = cbind(rmse_mean_FMOU_n2_est_d[,2], rmse_mean_DMD_n2_est_d[,2], rmse_mean_QW1_n2_est_d[,2], rmse_mean_QW5_n2_est_d[,2])

est_d_all_result_n2_noise_1 = cbind(est_d_FMOU_n2[,1], est_d_DMD_n2[,1], est_d_QW1_n2[,1], est_d_QW5_n2[,1])
est_d_all_result_n2_noise_2 = cbind(est_d_FMOU_n2[,2], est_d_DMD_n2[,2], est_d_QW1_n2[,2], est_d_QW5_n2[,2])

######################################## n=400 ########################################
angle_diff_FMOU_n3 = matrix(NA, N, length(sigma_0_2_list)) 
angle_diff_QW1_n3 = matrix(NA, N, length(sigma_0_2_list)) 
angle_diff_QW5_n3 = matrix(NA, N, length(sigma_0_2_list)) 
angle_diff_DMD_n3 = matrix(NA, N, length(sigma_0_2_list)) 

angle_diff_FMOU_n3_est_d = matrix(NA, N, length(sigma_0_2_list)) 
angle_diff_QW1_n3_est_d = matrix(NA, N, length(sigma_0_2_list))
angle_diff_QW5_n3_est_d = matrix(NA, N, length(sigma_0_2_list))
angle_diff_DMD_n3_est_d = matrix(NA, N, length(sigma_0_2_list))

rmse_mean_FMOU_n3 = matrix(NA, N, length(sigma_0_2_list)) 
rmse_mean_QW1_n3 = matrix(NA, N, length(sigma_0_2_list)) 
rmse_mean_QW5_n3 = matrix(NA, N, length(sigma_0_2_list)) 
rmse_mean_DMD_n3 = matrix(NA, N, length(sigma_0_2_list)) 

rmse_mean_FMOU_n3_est_d = matrix(NA, N, length(sigma_0_2_list)) 
rmse_mean_QW1_n3_est_d = matrix(NA, N, length(sigma_0_2_list))
rmse_mean_QW5_n3_est_d = matrix(NA, N, length(sigma_0_2_list))
rmse_mean_DMD_n3_est_d = matrix(NA, N, length(sigma_0_2_list))

L_95_n3 = matrix(NA, nrow=N, ncol=length(sigma_0_2_list)) # d is given
L_95_n3_est_d = matrix(NA, nrow=N, ncol=length(sigma_0_2_list)) # d is estimated

P_95_n3 = matrix(NA, nrow=N, ncol=length(sigma_0_2_list)) # d is given
P_95_n3_est_d = matrix(NA, nrow=N, ncol=length(sigma_0_2_list)) # d is estimated

est_d_FMOU_n3 = matrix(NA, N, length(sigma_0_2_list))
est_d_QW1_n3 = matrix(NA, N, length(sigma_0_2_list))
est_d_QW5_n3 = matrix(NA, N, length(sigma_0_2_list))
est_d_DMD_n3 = matrix(NA, N, length(sigma_0_2_list))

example_obs_n3 = as.list(1:length(sigma_0_2_list))
example_mean_n3 = as.list(1:length(sigma_0_2_list))
example_pred_mean_FMOU_n3 = as.list(1:length(sigma_0_2_list))
example_pred_mean_QW1_n3 = as.list(1:length(sigma_0_2_list))
example_pred_mean_QW5_n3 = as.list(1:length(sigma_0_2_list))
example_pred_mean_DMD_n3 = as.list(1:length(sigma_0_2_list))
example_pred_mean_uncertainty_FMOU_n3 = as.list(1:length(sigma_0_2_list))

time_record_FMOU_n3 = matrix(NA, N, length(sigma_0_2_list)) 


for(idx_noise in 1: length(sigma_0_2_list)){
  # specify noise level
  noise_level= sigma_0_2_list[idx_noise]
  for(it in 1:N){
    set.seed(3*it)
    if(it%%10 == 0){print(it)}
    ## generate latent processes and factor loading matrix
    U = rustiefel(k, k)
    z = matrix(NA, d, n3)
    sigma_2 = runif(d, 0.5, 1)
    rho = runif(d, 0.95, 1)
    for(l in 1:d){
      R = matrix(NA, n3, n3)
      diag(R) = 1
      for(ir in 1:n3){
        for(ic in 1:n3){
          R[ir, ic] = rho[l]^(abs(ir-ic)) * R[ir, ir]
        }
      }
      R = (sigma_2[l]/(1-rho[l]^2) )* R
      z[l, ] = t(chol(R)) %*% rnorm(n3)
    }
    ## record noise level
    signal = U[,1:d] %*% z
    y = signal + matrix(rnorm(n3*k,mean=0,sd=sqrt(noise_level)),k,n3)
    if(it==1){
      example_obs_n3[[idx_noise]] = y
      example_mean_n3[[idx_noise]] = signal
    }
    
    ## 1. FMOU
    # select d
    svd_output=svd(y)
    FMOU_loss_score = NULL
    criteria_val_pre = log(mean((y - svd_output$u[,1]%*%t(svd_output$u[,1])%*%y)^2)) + i_d*(k+n3)/(k*n3)*log(k*n3/(k+n3))
    FMOU_loss_score = c(FMOU_loss_score,  criteria_val_pre)
    
    for(i_d in 2:ceiling(dim(y)[1]*2/3)){
      criteria_val_cur = log(mean((y - svd_output$u[,1:i_d]%*%t(svd_output$u[,1:i_d])%*%y)^2)) + i_d*(k+n3)/(k*n3)*log(k*n3/(k+n3))
      FMOU_loss_score = c(FMOU_loss_score,  criteria_val_cur)
    }
    est_d_FMOU_here = which.min(FMOU_loss_score) # i_d-1
    est_d_FMOU_n3[it, idx_noise] = which.min(FMOU_loss_score) # i_d-1
    
    # fit FMOU
    time1 = Sys.time()
    tilde_y = t(U[, 1:d])%*%y
    em_alg <- EM_alg_FMOU(output=y, d=d, M=1000, threshold=10^{-6},
                          est_U0=T, U0=NA, est_sigma0_2=T)
    U_est_FMOU = em_alg$U
    time2 = Sys.time()
    
    tilde_y = t(U[, 1:est_d_FMOU_here])%*%y
    em_alg_est_d <- EM_alg_FMOU(output=y, d=est_d_FMOU_here, M=1000, threshold=10^{-6},
                          est_U0=T, U0=NA, est_sigma0_2=T)
    U_est_FMOU_est_d = em_alg_est_d$U
    
    time_record_FMOU_n3[it, idx_noise] = time2 - time1
    
    fit_latent_state_EM = U_est_FMOU %*% em_alg$Z_hat
    fit_latent_state_EM_est_d = U_est_FMOU_est_d %*% em_alg_est_d$Z_hat
    
    angle_diff_FMOU_n3[it, idx_noise] = rospca::angle(U[,1:d], U_est_FMOU)
    angle_diff_FMOU_n3_est_d[it, idx_noise] = rospca::angle(U[,1:d], U_est_FMOU_est_d)
    rmse_mean_FMOU_n3[it, idx_noise] = sqrt(mean(( fit_latent_state_EM - signal)^2))
    rmse_mean_FMOU_n3_est_d[it, idx_noise] = sqrt(mean(( fit_latent_state_EM_est_d - signal)^2))
    
    L_95_n3[it,idx_noise] = mean(em_alg$pred_mean_95upper-em_alg$pred_mean_95lower)
    L_95_n3_est_d[it,idx_noise] = mean(em_alg_est_d$pred_mean_95upper-em_alg_est_d$pred_mean_95lower)
    
    P_95_n3[it, idx_noise] = mean(em_alg$pred_mean_95upper>signal &em_alg$pred_mean_95lower< signal)
    P_95_n3_est_d[it, idx_noise] = mean(em_alg_est_d$pred_mean_95upper>signal &em_alg_est_d$pred_mean_95lower< signal)
    
    if(it==1){
      example_pred_mean_FMOU_n3[[idx_noise]] = list(with_true_d = fit_latent_state_EM,
                                                    with_est_d = fit_latent_state_EM_est_d)
      #example_pred_mean_uncertainty_FMOU_n3[[idx_noise]] = list(upper =em_alg$pred_mean_95upper,
      #                                                          lower =em_alg$pred_mean_95lower)
    }
    
    # 2. DMD
    DMD_fit = DMD_alg(y, fix_r=T, r=d)
    DMD_fit_est_d = DMD_alg(y, threshold=0.99)
    rmse_mean_DMD_n3[it, idx_noise] = sqrt(mean((cbind(y[,1], DMD_fit$in_sample_pred) - signal)^2))
    rmse_mean_DMD_n3_est_d[it, idx_noise] = sqrt(mean((cbind(y[,1], DMD_fit_est_d$in_sample_pred) - signal)^2))
    if(it==1){
      example_pred_mean_DMD_n3[[idx_noise]] = list(with_true_d = cbind(y[,1], DMD_fit$in_sample_pred),
                                                   with_est_d = cbind(y[,1], DMD_fit_est_d$in_sample_pred)
      )
    }
    
    angle_diff_DMD_n3[it, idx_noise] = rospca::angle(U[,1:d], DMD_fit$A_hat)
    angle_diff_DMD_n3_est_d[it, idx_noise] = rospca::angle(U[,1:d], DMD_fit_est_d$A_hat)
    est_d_DMD_n3[it, idx_noise] = DMD_fit_est_d$rank
    
    # 3. LY1
    cov_lag_1 = matrix(0, k, k)
    for(lag in 1:1){
      cov_lag_1_here = matrix(NA, k, k)
      for (i in 1:k) {
        for (j in 1:k) {
          cov_lag_1_here[i, j] <- cov(y[i , 1: (ncol(y)-lag)], y[j, (1+lag):(ncol(y))])
        }
      }
      cov_lag_1 = cov_lag_1 + cov_lag_1_here %*% t(cov_lag_1_here)
    }
    eigen_QW_1 = eigen(cov_lag_1)
    lambda_hat = eigen_QW_1$values
    eigen_ratio = lambda_hat[2:k] / lambda_hat[1:(k-1)]
    est_d_QW1_here = which.min(eigen_ratio)
    est_d_QW1_n3[it, idx_noise] = est_d_QW1_here
    U_est_QW1 = eigen_QW_1$vectors[,1:d]
    U_est_QW1_est_d = eigen_QW_1$vectors[,1:est_d_QW1_here]
    rmse_mean_QW1_n3[it, idx_noise] = sqrt(mean(( U_est_QW1%*%t(U_est_QW1)%*%y - signal)^2))
    rmse_mean_QW1_n3_est_d[it, idx_noise] = sqrt(mean(( U_est_QW1_est_d%*%t(U_est_QW1_est_d)%*%y - signal)^2))
    if(it==1){
      example_pred_mean_QW1_n3[[idx_noise]] = list(with_true_d = U_est_QW1%*%t(U_est_QW1)%*%y,
                                                   with_est_d = U_est_QW1_est_d%*%t(U_est_QW1_est_d)%*%y
      )
    }
    
    angle_diff_QW1_n3[it, idx_noise] = rospca::angle(U[,1:d], U_est_QW1)
    angle_diff_QW1_n3_est_d[it, idx_noise] = rospca::angle(U[,1:d], U_est_QW1_est_d)
    
    # 4. LY5
    cov_lag_5 = matrix(0, k, k)
    for(lag in 1:5){
      cov_lag_5_here = matrix(NA, k, k)
      for (i in 1:k) {
        for (j in 1:k) {
          cov_lag_5_here[i, j] <- cov(y[i , 1: (ncol(y)-lag)], y[j, (1+lag):(ncol(y))])
        }
      }
      cov_lag_5 = cov_lag_5 + cov_lag_5_here %*% t(cov_lag_5_here)
    }
    eigen_QW_5 = eigen(cov_lag_5)
    lambda_hat = eigen_QW_5$values
    eigen_ratio = lambda_hat[2:k] / lambda_hat[1:(k-1)]
    est_d_QW5_here = which.min(eigen_ratio)
    est_d_QW5_n3[it, idx_noise] = est_d_QW5_here
    U_est_QW5 = eigen_QW_5$vectors[,1:d]
    U_est_QW5_est_d = eigen_QW_5$vectors[,1:est_d_QW5_here]
    rmse_mean_QW5_n3[it, idx_noise] = sqrt(mean(( U_est_QW5%*%t(U_est_QW5)%*%y - signal)^2))
    rmse_mean_QW5_n3_est_d[it, idx_noise] = sqrt(mean(( U_est_QW5_est_d%*%t(U_est_QW5_est_d)%*%y - signal)^2))
    if(it==1){
      example_pred_mean_QW5_n3[[idx_noise]] = list(with_true_d = U_est_QW5%*%t(U_est_QW5)%*%y,
                                                   with_est_d = U_est_QW5_est_d%*%t(U_est_QW5_est_d)%*%y
      )
    }
    
   
    angle_diff_QW5_n3[it, idx_noise] = rospca::angle(U[,1:d], U_est_QW5)
    angle_diff_QW5_n3_est_d[it, idx_noise] = rospca::angle(U[,1:d], U_est_QW5_est_d)
  }
}


angle_diff_n3_noise_1 = cbind(angle_diff_FMOU_n3[,1], angle_diff_DMD_n3[,1], angle_diff_QW1_n3[,1], angle_diff_QW5_n3[,1])
angle_diff_noise_1 = cbind(angle_diff_n1_noise_1, angle_diff_n2_noise_1, angle_diff_n3_noise_1)

angle_diff_n3_noise_2 = cbind(angle_diff_FMOU_n3[,2], angle_diff_DMD_n3[,2], angle_diff_QW1_n3[,2], angle_diff_QW5_n3[,2])
angle_diff_noise_2 = cbind(angle_diff_n1_noise_2, angle_diff_n2_noise_2, angle_diff_n3_noise_2)

angle_diff_n3_est_d_noise_1 = cbind(angle_diff_FMOU_n3_est_d[,1], angle_diff_DMD_n3_est_d[,1], angle_diff_QW1_n3_est_d[,1], angle_diff_QW5_n3_est_d[,1])
angle_diff_noise_est_d_1 = cbind(angle_diff_n1_est_d_noise_1, angle_diff_n2_est_d_noise_1, angle_diff_n3_est_d_noise_1)

angle_diff_n3_est_d_noise_2 = cbind(angle_diff_FMOU_n3_est_d[,2], angle_diff_DMD_n3_est_d[,2], angle_diff_QW1_n3_est_d[,2], angle_diff_QW5_n3_est_d[,2])
angle_diff_noise_est_d_2 = cbind(angle_diff_n1_est_d_noise_2, angle_diff_n2_est_d_noise_2, angle_diff_n3_est_d_noise_2)

rmse_mean_n3_noise_1 = cbind(rmse_mean_FMOU_n3[,1], rmse_mean_DMD_n3[,1], rmse_mean_QW1_n3[,1], rmse_mean_QW5_n3[,1])
rmse_mean_noise_1 = cbind(rmse_mean_n1_noise_1, rmse_mean_n2_noise_1, rmse_mean_n3_noise_1)

rmse_mean_n3_noise_2 = cbind(rmse_mean_FMOU_n3[,2], rmse_mean_DMD_n3[,2], rmse_mean_QW1_n3[,2], rmse_mean_QW5_n3[,2])
rmse_mean_noise_2 = cbind(rmse_mean_n1_noise_2, rmse_mean_n2_noise_2, rmse_mean_n3_noise_2)

rmse_mean_n3_est_d_noise_1 = cbind(rmse_mean_FMOU_n3_est_d[,1], rmse_mean_DMD_n3_est_d[,1], rmse_mean_QW1_n3_est_d[,1], rmse_mean_QW5_n3_est_d[,1])
rmse_mean_noise_est_d_1 = cbind(rmse_mean_n1_est_d_noise_1, rmse_mean_n2_est_d_noise_1, rmse_mean_n3_est_d_noise_1)

rmse_mean_n3_est_d_noise_2 = cbind(rmse_mean_FMOU_n3_est_d[,2], rmse_mean_DMD_n3_est_d[,2], rmse_mean_QW1_n3_est_d[,2], rmse_mean_QW5_n3_est_d[,2])
rmse_mean_noise_est_d_2 = cbind(rmse_mean_n1_est_d_noise_2, rmse_mean_n2_est_d_noise_2, rmse_mean_n3_est_d_noise_2)

est_d_all_result_n3_noise_1 = cbind(est_d_FMOU_n3[,1], est_d_DMD_n3[,1], est_d_QW1_n3[,1], est_d_QW5_n3[,1])
est_d_all_result_noise_1 = cbind(est_d_all_result_n1_noise_1, est_d_all_result_n2_noise_1, est_d_all_result_n3_noise_1)

est_d_all_result_n3_noise_2 = cbind(est_d_FMOU_n3[,2], est_d_DMD_n3[,2], est_d_QW1_n3[,2], est_d_QW5_n3[,2])
est_d_all_result_noise_2 = cbind(est_d_all_result_n1_noise_2, est_d_all_result_n2_noise_2, est_d_all_result_n3_noise_2)

################################# Tables #################################
# 1. Table 2 in manuscript
rmse_mean_fix_d_noise1 = cbind(apply(rmse_mean_n1_noise_1,2,mean),
                               apply(rmse_mean_n2_noise_1,2,mean),
                               apply(rmse_mean_n3_noise_1,2,mean))
rownames(rmse_mean_fix_d_noise1) = c("FMOU","DMD","LY1","LY5")
colnames(rmse_mean_fix_d_noise1) = c("n=100", "n=200", "n=400")
rmse_mean_fix_d_noise1

rmse_mean_fix_d_noise2 = cbind(apply(rmse_mean_n1_noise_2,2,mean),
                               apply(rmse_mean_n2_noise_2,2,mean),
                               apply(rmse_mean_n3_noise_2,2,mean))
rownames(rmse_mean_fix_d_noise2) = c("FMOU","DMD","LY1","LY5")
colnames(rmse_mean_fix_d_noise2) = c("n=100", "n=200", "n=400")
rmse_mean_fix_d_noise2


# Table S1 in supplement
rmse_mean_est_d_noise1 = cbind(apply(rmse_mean_n1_est_d_noise_1,2,mean),
                               apply(rmse_mean_n2_est_d_noise_1,2,mean),
                               apply(rmse_mean_n3_est_d_noise_1,2,mean))
rownames(rmse_mean_est_d_noise1) = c("FMOU","DMD","LY1","LY5")
colnames(rmse_mean_est_d_noise1) = c("n=100", "n=200", "n=400")
rmse_mean_est_d_noise1

rmse_mean_est_d_noise2 = cbind(apply(rmse_mean_n1_est_d_noise_2,2,mean),
                               apply(rmse_mean_n2_est_d_noise_2,2,mean),
                               apply(rmse_mean_n3_est_d_noise_2,2,mean))
rownames(rmse_mean_est_d_noise2) = c("FMOU","DMD","LY1","LY5")
colnames(rmse_mean_est_d_noise2) = c("n=100", "n=200", "n=400")
rmse_mean_est_d_noise2


# Table S2 : coverage and length of 95% CI
coverage_95CI_true_d = rbind(apply(P_95_n1,2,mean), apply(P_95_n2,2,mean), apply(P_95_n3,2,mean))
rownames(coverage_95CI_true_d) <- c("n=100", "n=200", "n=400")
colnames(coverage_95CI_true_d) <- c("sigma0^2=1", "sigma0^2=2")
coverage_95CI_true_d

len_95CI_true_d = rbind(apply(L_95_n1,2,mean), apply(L_95_n2,2,mean), apply(L_95_n3,2,mean))
rownames(len_95CI_true_d) <- c("n=100", "n=200", "n=400")
colnames(len_95CI_true_d) <- c("sigma0^2=1", "sigma0^2=2")
len_95CI_true_d


coverage_95CI_est_d = rbind(apply(P_95_n1_est_d,2,mean), apply(P_95_n2_est_d,2,mean), apply(P_95_n3_est_d,2,mean))
rownames(coverage_95CI_est_d) <- c("n=100", "n=200", "n=400")
colnames(coverage_95CI_est_d) <- c("sigma0^2=1", "sigma0^2=2")
coverage_95CI_est_d

len_95CI_est_d = rbind(apply(L_95_n1_est_d,2,mean), apply(L_95_n2_est_d,2,mean), apply(L_95_n3_est_d,2,mean))
rownames(len_95CI_est_d) <- c("n=100", "n=200", "n=400")
colnames(len_95CI_est_d) <- c("sigma0^2=1", "sigma0^2=2")
len_95CI_est_d

####################################### Plots #####################################
# 1. Figure 1: Largest principal angle between U0 and hat_U0
par(mfrow=c(1,2),  mar=c(3.9+0.15, 3.8, 3.5, 1.8), oma=c(2, 2, 2, 2))
par(las=2) 
boxplot(angle_diff_noise_1*(pi/2), las=2, names=c('FMOU','DMD', 'LY1','LY5', 'FMOU','DMD', 'LY1','LY5', 'FMOU','DMD', 'LY1','LY5'),
        col=c("#00a6fb", "#8758FF", "#FF9100","#6C946F", "#00a6fb", "#8758FF", "#FF9100","#6C946F", "#00a6fb", "#8758FF", "#FF9100","#6C946F"),
        at=c(1:4, 6:9, 11:14), ylab="Largest principle angle",
        main=expression(d == 5 ~ "," ~ sigma[0]^2 == 1),
        ylim=c(0,1))

boxplot(angle_diff_noise_2*(pi/2), las=2, names=c('FMOU','DMD', 'LY1','LY5', 'FMOU','DMD', 'LY1','LY5', 'FMOU','DMD', 'LY1','LY5'),
        col=c("#00a6fb", "#8758FF", "#FF9100","#6C946F", "#00a6fb", "#8758FF", "#FF9100","#6C946F", "#00a6fb", "#8758FF", "#FF9100","#6C946F"), 
        at=c(1:4, 6:9, 11:14), ylab="Largest principle angle",
        main=expression(d == 5 ~ "," ~ sigma[0]^2 == 2),
        ylim=c(0,1.5))
par(mfrow=c(1,1))


# 2. Figure 3: Predictive signal by FMOU, DMD, LY1, LY5
row_index = 1
fix_d_plot_df = data.frame(method=factor(rep(c("Signal","FMOU","DMD","LY1","LY5"), each=n1),levels=c("Signal", "FMOU","DMD","LY1","LY5")),
                           val = c(example_mean_n1[[1]][row_index,], 
                                   example_pred_mean_FMOU_n1[[1]]$with_true_d[row_index,],
                                   example_pred_mean_DMD_n1[[1]]$with_true_d[row_index,],
                                   example_pred_mean_QW1_n1[[1]]$with_true_d[row_index,],
                                   example_pred_mean_QW5_n1[[1]]$with_true_d[row_index,]),
                           time = rep(1:n1, 5))

obs_df = data.frame(time=1:n1, obs =example_obs_n1[[1]][row_index,] )
em_pred_interval = data.frame(time = 1:n1,
                              lower = pred_mean_FMOU_lower_n1[[1]][row_index, ],
                              upper = pred_mean_FMOU_upper_n1[[1]][row_index, ]) 


fig3_left <- ggplot() + 
  geom_ribbon(data=em_pred_interval, aes(x=time, ymin=lower, ymax=upper), 
              fill="#00a6fb", alpha = 0.3) + 
  geom_path(data=fix_d_plot_df, aes(time, val, group=method,
                                    col=method, linetype=method, linewidth=method)) +
  geom_point(data=obs_df, aes(time, obs), col="black", size=1.2, shape=20) + 
  #geom_point(data=fix_d_plot_df, aes(time, val, group=method,
  #                                   col=method, shape=method)) + 
  scale_linetype_manual(values = c("solid","solid", "dashed", "solid", "dotdash")) + 
  scale_linewidth_manual(values=c(0.8,0.9,0.8,0.8,0.8)) + 
  scale_color_manual(values=c("black", "#00a6fb", "#8758FF", "#FF9100","#6C946F")) + 
  #scale_shape_manual(values=c(20, 1, 4, 5, 18, 15)) + 
  theme_classic() +
  ggtitle(expression(sigma[0]^2 == 1)) + 
  xlab('t') +
  ylab(expression(hat(y)[1])) +
  ylim(c(range(obs_df$obs)[1], range(obs_df$obs)[2])) + 
  theme_classic() + 
  theme(#legend.position=c(0.9,0.79), 
    legend.position="none",
    legend.text = element_text(size = 9 ),
    legend.title= element_blank(),
    axis.title = element_text(size = 12), 
    axis.text = element_text(size=12),
    panel.border = element_rect(colour = "black", fill=NA, linewidth =0.9),
    plot.title = element_text(hjust = 0.5, size = 18)
  )

fig3_left


row_index = 3
fix_d_plot_df = data.frame(method=factor(rep(c("Signal","FMOU","DMD","LY1","LY5"), each=n1),levels=c("Signal", "FMOU","DMD","LY1","LY5")),
                           val = c(example_mean_n1[[2]][row_index,], 
                                   example_pred_mean_FMOU_n1[[2]]$with_true_d[row_index,],
                                   example_pred_mean_DMD_n1[[2]]$with_true_d[row_index,],
                                   example_pred_mean_QW1_n1[[1]]$with_true_d[row_index,],
                                   example_pred_mean_QW5_n1[[2]]$with_true_d[row_index,]),
                           time = rep(1:n1, 5))

obs_df = data.frame(time=1:n1, obs =example_obs_n1[[2]][row_index,] )
em_pred_interval = data.frame(time = 1:n1,
                              lower = pred_mean_FMOU_lower_n1[[2]][row_index, ],
                              upper = pred_mean_FMOU_upper_n1[[2]][row_index, ]) 

fig3_right <- ggplot() + 
  geom_ribbon(data=em_pred_interval, aes(x=time, ymin=lower, ymax=upper), 
              fill="#00a6fb", alpha = 0.3) + 
  geom_path(data=fix_d_plot_df, aes(time, val, group=method,
                                    col=method, linetype=method, linewidth=method)) +
  geom_point(data=obs_df, aes(time, obs), col="black", size=1.2, shape=20) + 
  #geom_point(data=fix_d_plot_df, aes(time, val, group=method,
  #                                   col=method, shape=method)) + 
  scale_linetype_manual(values = c("solid","solid", "dashed", "solid", "dotdash")) + 
  scale_linewidth_manual(values=c(0.8,0.9,0.8,0.8,0.8)) + 
  scale_color_manual(values=c("black", "#00a6fb", "#8758FF", "#FF9100","#6C946F")) + 
  #scale_shape_manual(values=c(20, 1, 4, 5, 18, 15)) + 
  theme_classic() +
  xlab('t') +
  ylab(expression(hat(y)[3])) +
  ylim(c(range(obs_df$obs)[1]-0.8, range(obs_df$obs)[2]+1.5)) + 
  theme_classic() + 
  ggtitle(expression(sigma[0]^2 == 2)) + 
  theme(legend.position=c(0.1,0.2), 
        #legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title= element_blank(),
        axis.title = element_text(size = 12), 
        axis.text = element_text(size=12),
        panel.border = element_rect(colour = "black", fill=NA, linewidth =0.9),
        plot.title = element_text(hjust = 0.5, size = 18)
  )

fig3_right

#pdf("figures/Exp1_pred_n100_fix_d.pdf", height=4.5, width=11)
grid.arrange(fig3_left, fig3_right, nrow=1)
#dev.off()

## 3. Figure S1: Difference between selected d and true d by different methods
par(mfrow=c(1,2))
par(mar=c(5, 4.1, 4, 2) + 0.1) 
par(las=2) 
d_min = min(est_d_all_result_noise_1)
d_max = max(est_d_all_result_noise_1)
boxplot(est_d_all_result_noise_1 - d, las=2, names=c('FMOU','DMD', 'LY1','LY5', 'FMOU','DMD', 'LY1','LY5', 'FMOU','DMD', 'LY1','LY5'),
        main=expression(k==20 * " and " * "d = 5"), 
        col=c("#00a6fb", "#8758FF", "#FF9100","#6C946F", "#00a6fb", "#8758FF", "#FF9100","#6C946F", "#00a6fb", "#8758FF", "#FF9100","#6C946F"), 
        ylim=c(d_min-d-1, d_max-d), at=c(1:4, 6:9, 11:14), ylab=expression(hat(d) - d))

d_min = min(est_d_all_result_noise_2)
d_max = max(est_d_all_result_noise_2)
boxplot(est_d_all_result_noise_2 - d, las=2, names=c('FMOU','DMD', 'LY1','LY5', 'FMOU','DMD', 'LY1','LY5', 'FMOU','DMD', 'LY1','LY5'),
        main=expression(k==20 * " and " * "d = 5"), 
        col=c("#00a6fb", "#8758FF", "#FF9100","#6C946F", "#00a6fb", "#8758FF", "#FF9100","#6C946F", "#00a6fb", "#8758FF", "#FF9100","#6C946F"), 
        ylim=c(d_min-d-1, d_max-d), at=c(1:4, 6:9, 11:14), ylab=expression(hat(d) - d))

#save.image("saved_workspace/Exp1.RData")
