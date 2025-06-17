# this is the codes for Exp3 - Linear diffusion equation
library(rstiefel)
library(ReacTran)
library(deSolve)
library(plot3D)
library(MASS)
library(ggplot2)
library(Rcpp)
library(astsa)

source("../functions/FMOU.R")
source('../functions/DMD.R')
source('../functions/EM_VAR1.R')


N = 20 # repeat times
k = 300 # dimension of observations
n = 300
sigma_0_2_list = c(0.01, 0.05, 0.3)^2

# record RMSE
rmse_mean_FMOU = matrix(NA, N, length(sigma_0_2_list)) 
rmse_mean_QW1 = matrix(NA, N, length(sigma_0_2_list)) 
rmse_mean_QW5 = matrix(NA, N, length(sigma_0_2_list)) 
rmse_mean_DMD = matrix(NA, N, length(sigma_0_2_list)) 
rmse_mean_IC = matrix(NA, N, length(sigma_0_2_list)) 
rmse_mean_ar1_em = matrix(NA, N, length(sigma_0_2_list))

# record estimated number of parameters
est_d_FMOU = matrix(NA, N, length(sigma_0_2_list))
est_d_DMD = matrix(NA, N, length(sigma_0_2_list))
est_d_QW1 = matrix(NA, N, length(sigma_0_2_list))
est_d_QW5 = matrix(NA, N, length(sigma_0_2_list))
est_d_IC = matrix(NA, N, length(sigma_0_2_list))

# record signal and observations
signal = as.list(1:length(sigma_0_2_list))
observation = as.list(1:length(sigma_0_2_list))

# record computational time of selecting d
time_est_d_FMOU = matrix(NA, N, length(sigma_0_2_list))

# record computational time (exclude the time of selecting d)
time_FMOU = matrix(NA, N, length(sigma_0_2_list))
time_DMD = matrix(NA, N, length(sigma_0_2_list))
time_QW1 = matrix(NA, N, length(sigma_0_2_list))
time_QW5 = matrix(NA, N, length(sigma_0_2_list))
time_ar1_em = matrix(NA, N, length(sigma_0_2_list))

# record predictive signals
pred_mean_FMOU = as.list(1:length(sigma_0_2_list))
pred_mean_DMD = as.list(1:length(sigma_0_2_list))
pred_mean_QW1 = as.list(1:length(sigma_0_2_list))
pred_mean_QW5 = as.list(1:length(sigma_0_2_list))
pred_mean_IC = as.list(1:length(sigma_0_2_list))
pred_mean_ar1_em = as.list(1:length(sigma_0_2_list))

# length of 95% CI
L_95 = matrix(NA, nrow=N, ncol=length(sigma_0_2_list)) 

# coverage of 95% CI
P_95 = matrix(NA, nrow=N, ncol=length(sigma_0_2_list)) 


for(idx_noise in 1: length(sigma_0_2_list)){
  noise_level= sigma_0_2_list[idx_noise]
  
  for(it in 1:N){
    if(it%%10 == 0){print(it)}
    set.seed(it)
    
    #################### generate signal by numerical solver #######################
    Grid <- setup.grid.1D(N=k, L=1)
    pde1D <- function(t, C, parms){
      tran <- tran.1D(C=C, D=D, C.down=Cext, dx=Grid)$dC
      list(tran-Q)
    }
    # initial states
    D = 1
    Q = 0
    Cext = 1
    times <- seq(0, 0.2, length.out=n)
    reality <- ode.1D(y = rep(0, Grid$N), times = times, func = pde1D, parms = NULL, nspec = 1)
    reality <- t(reality[,2:(k+1)]) # values between 0 and 1
    if(it==1){
      signal[[idx_noise]] = reality
    }
    
    ################ generate noisy observations ################
    y = reality + matrix(rnorm(n*k,mean=0,sd=sqrt(noise_level)),k,n)
    if(it==1){
      observation[[idx_noise]] = y
    }
    
    ############################ FMOU  ############################
    # 1. select d by IC
    time = system.time({
      svd_output=svd(y)
      U = svd_output$u
      FMOU_loss_score = NULL
      for(i_d in 1:ceiling(dim(y)[1]*2/3)){
        criteria_val_cur = log(mean((y - U[,1:i_d]%*%t(U[,1:i_d])%*%y)^2)) + i_d*(k+n)/(k*n)*log(k*n/(k+n))
        FMOU_loss_score = c(FMOU_loss_score,  criteria_val_cur)
      }
      est_d_FMOU_here = which.min(FMOU_loss_score)
      est_d_FMOU[it, idx_noise] = est_d_FMOU_here
    })
    time_est_d_FMOU[it, idx_noise] = time[3]
    
    # fit FMOU
    time = system.time({
      em_alg <- EM_alg_FMOU(output=y, d=est_d_FMOU_here, M=1000, threshold=10^{-6},
                            est_U0=T, est_sigma0_2=T)
    })
    rmse_mean_FMOU[it, idx_noise] = sqrt(mean(( em_alg$pred_mean - reality)^2))
    L_95[it, idx_noise] = mean(em_alg$pred_mean_95upper-em_alg$pred_mean_95lower)
    P_95[it, idx_noise] = mean(em_alg$pred_mean_95upper>reality & em_alg$pred_mean_95lower<reality)
    time_FMOU[it,idx_noise] = time[3]
    
    if(it==1){
      pred_mean_FMOU[[idx_noise]] = em_alg$pred_mean
    }
    
    ###################### DMD ######################
    time = system.time({
      DMD_fit = DMD_alg(y, threshold=0.99)
    })
    rmse_mean_DMD[it, idx_noise] = sqrt(mean((cbind(y[,1], DMD_fit$in_sample_pred) - reality)^2))
    est_d_DMD[it, idx_noise] = DMD_fit$rank
    time_DMD[it, idx_noise] = time[3]
    if(it==1){
      pred_mean_DMD[[idx_noise]] = cbind(y[,1], DMD_fit$in_sample_pred)
    }

    # ####################### LY1 ######################
    time = system.time({
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
      U_est_QW1 = eigen_QW_1$vectors[,1:est_d_QW1_here]
    })
    time_QW1[it,idx_noise] = time[3]
    est_d_QW1[it, idx_noise] = est_d_QW1_here
   

    rmse_mean_QW1[it, idx_noise] = sqrt(mean(( U_est_QW1%*%t(U_est_QW1)%*%y - reality)^2))

    if(it==1){
      pred_mean_QW1[[idx_noise]] = U_est_QW1%*%t(U_est_QW1)%*%y
    }

    # ####################### LY5 ######################
    time = system.time({
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
      U_est_QW5 = eigen_QW_5$vectors[,1:est_d_QW5_here]
    })
    
    est_d_QW5[it, idx_noise] = est_d_QW5_here
    time_QW5[it, idx_noise] = time[3]

    rmse_mean_QW5[it, idx_noise] = sqrt(mean(( U_est_QW5%*%t(U_est_QW5)%*%y - reality)^2))
    if(it==1){
      pred_mean_QW5[[idx_noise]] = U_est_QW5%*%t(U_est_QW5)%*%y
    }
    
    # ####################### AR1_EM ######################
    # A = diag(k)
    # #### estimate parameters by EM
    # time = system.time({
    #   ar1_em_fit = EM_SS(y, A, n_iter=8)
    # })
    #   #### predict mean
    # kf = KF_multidim(A, ar1_em_fit$Q, ar1_em_fit$Phi, ar1_em_fit$R, ar1_em_fit$mu, ar1_em_fit$Sigma, y)
    # time_ar1_em[it, idx_noise] = time[3]
    # rmse_mean_ar1_em[it, idx_noise] = sqrt(mean((reality-kf$s[,2:(n+1)])^2))
    # if(it==1){
    # pred_mean_ar1_em[[idx_noise]] = kf$s[,2:(n+1)]
    # }
    
  }
}

rmse_mean_noise_1 = cbind(rmse_mean_FMOU[,1], rmse_mean_DMD[,1], rmse_mean_QW1[,1], rmse_mean_QW5[,1], rmse_mean_ar1_em[,1])
rmse_mean_noise_2 = cbind(rmse_mean_FMOU[,2], rmse_mean_DMD[,2], rmse_mean_QW1[,2], rmse_mean_QW5[,2], rmse_mean_ar1_em[,2])
rmse_mean_noise_3 = cbind(rmse_mean_FMOU[,3], rmse_mean_DMD[,3], rmse_mean_QW1[,3], rmse_mean_QW5[,3], rmse_mean_ar1_em[,3])

time_noise_1 = cbind(time_FMOU[,1], time_DMD[,1], time_QW1[,1], time_QW5[,1], time_ar1_em[,1])
time_noise_2 = cbind(time_FMOU[,2], time_DMD[,2], time_QW1[,2], time_QW5[,2], time_ar1_em[,2])
time_noise_3 = cbind(time_FMOU[,3], time_DMD[,3], time_QW1[,3], time_QW5[,3], time_ar1_em[,3])

######################### Table #########################
# 1.Table 3
rmse_mean = cbind(apply(rmse_mean_noise_1,2,mean),
                  apply(rmse_mean_noise_2,2,mean),
                  apply(rmse_mean_noise_3,2,mean))
rownames(rmse_mean) = c("FMOU", "DMD","LY1","LY5","AR1-EM")
colnames(rmse_mean) = sigma_0_2_list
rmse_mean


# 2. Table S3 
apply(L_95,2,mean)
apply(P_95,2,mean)

# 3. avg computational time
time_mean = cbind(apply(time_noise_1,2,mean),
                  apply(time_noise_2,2,mean),
                  apply(time_noise_3,2,mean))
rownames(time_mean) = c("FMOU", "DMD","LY1","LY5","AR1-EM")
colnames(time_mean) = sigma_0_2_list
rowMeans(time_mean)


######################### Figure 4 #########################
fmou_fit = EM_alg_FMOU(output=observation[[2]], d=est_d_FMOU[1,2],M=100, threshold=10^{-6},
                       est_U0=T, est_sigma0_2=T, track_iterations = T, neg_log_lik = T)
mycolor = colorRampPalette(c("blue", "white", "red"))(100)
#png("../figures/linear_diffusion_pred_fmou_dmd_June25.png", width=1000, height=250)
par(mfrow=c(1,4),mgp = c(1.25, 0.5, 0))
par(mar=c(3,3,3,1.5))
image2D(t(reality), main="(a) True mean", xlab='t',ylab='x',xaxt = "n", 
        yaxt = "n",col = mycolor,cex.main=2.2,cex.lab = 2, colkey = list(cex.axis=1.8,cex.clab=1.8))
par(mar=c(3,3,3,1.5))
image2D(t(pred_mean_FMOU[[2]]), main="(b) Pred mean (FMOU)", 
        xlab="t", ylab="",xaxt = "n", yaxt = "n",
        col = mycolor,cex.main=2.2,cex.lab = 2, colkey = list(cex.axis=1.8,cex.clab=1.8))
par(mar=c(3,3,3,2.2))
image2D(t(pred_mean_DMD[[2]]), main="(c) Pred mean (DMD)",
        xlab="t",ylab="", xaxt = "n", yaxt = "n",col = mycolor,
        cex.main=2.2,cex.lab =2, colkey = list(cex.axis=1.8,cex.clab=1.8))
par(mar=c(4,5,3,2.5))
par(mgp=c(2.5, 1, 0)) # increase the distance between lab and ticks
plot(-fmou_fit$record_neg_lik,type='l',
     ylab='log likelihood',xlab='E-M iteration',main='(d) Log lik, FMOU',cex.main=2.2, 
     cex.lab=1.8,cex.axis=1.4)
#dev.off()


