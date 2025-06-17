library(lhs)
library(RobustGaSP)
library(plot3D)
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyr)
library(latex2exp)
library(FastGaSP)

source("../functions/FMOU.R")
source('../functions/DMD.R')
source('../functions/EM_VAR1.R')
#### Generate mean of observations ---------------------------
branin <- function(xx, a=1, b=5.1/(4*pi^2), c=5/pi, r=6, s=10, t=1/(8*pi))
{
  
  x1 <- xx[1]
  x2 <- xx[2]
  
  term1 <- a * (x2 - b*x1^2 + c*x1 - r)^2
  term2 <- s*(1-t)*cos(x1)
  
  y <- term1 + term2 + s
  return(y)
}

n1=300
n2=300
N=n1*n2
p=2 ##2D input
set.seed(1)
LB=c(-5,0)
UB=c(10,15)
range=UB-LB
input1=LB[1]+range[1]*seq(0,1,1/(n1-1)) ##row
input2=LB[2]+range[2]*seq(0,1,1/(n2-1)) ##column
input=cbind(rep(input1,n2),as.vector(t(matrix(input2,n2,n1))))  

output=matrix(NA,N,1)
f=matrix(NA,N,1)
for(i in 1:N){
  f[i]=branin(input[i,])
}

f_mat = matrix(f,n1,n2)

num_rep = 20
sigma_0_2_list = c(1,25,400)

# record RMSE
rmse_mean_FMOU = matrix(NA, num_rep, length(sigma_0_2_list)) 
rmse_mean_QW1 = matrix(NA, num_rep, length(sigma_0_2_list)) 
rmse_mean_QW5 = matrix(NA, num_rep, length(sigma_0_2_list)) 
rmse_mean_DMD = matrix(NA, num_rep, length(sigma_0_2_list)) 
rmse_mean_ar1_em = matrix(NA, num_rep, length(sigma_0_2_list))

# record estimated number of d
est_d_FMOU = matrix(NA, num_rep, length(sigma_0_2_list))
est_d_DMD = matrix(NA, num_rep, length(sigma_0_2_list))
est_d_QW1 = matrix(NA, num_rep, length(sigma_0_2_list))
est_d_QW5 = matrix(NA, num_rep, length(sigma_0_2_list))

# record signal and observations
signal = f_mat
observation = as.list(1:length(sigma_0_2_list))

# record computational time of selecting d
time_est_d_FMOU = matrix(NA, num_rep, length(sigma_0_2_list))

# record computational time (exclude the time of selecting d)
time_FMOU = matrix(NA, num_rep, length(sigma_0_2_list))
time_DMD = matrix(NA, num_rep, length(sigma_0_2_list))
time_QW1 = matrix(NA, num_rep, length(sigma_0_2_list))
time_QW5 = matrix(NA, num_rep, length(sigma_0_2_list))
time_ar1_em = matrix(NA, num_rep, length(sigma_0_2_list))

# record predictive signals
pred_mean_FMOU = as.list(1:length(sigma_0_2_list))
pred_mean_DMD = as.list(1:length(sigma_0_2_list))
pred_mean_QW1 = as.list(1:length(sigma_0_2_list))
pred_mean_QW5 = as.list(1:length(sigma_0_2_list))
pred_mean_ar1_em = as.list(1:length(sigma_0_2_list))

# length of 95% CI
L_95 = matrix(NA, nrow=num_rep, ncol=length(sigma_0_2_list)) 
# coverage of 95% CI
P_95 = matrix(NA, nrow=num_rep, ncol=length(sigma_0_2_list))

n=n1
k=n2
for(idx_noise in 1: length(sigma_0_2_list)){
  noise_level= sigma_0_2_list[idx_noise]
  
  for(it in 1:num_rep){
    if(it%%10 == 0){print(it)}
    set.seed(it)
    
    ################ generate noisy observations ################
    y = f_mat + matrix(rnorm(n*k,mean=0,sd=sqrt(noise_level)),k,n)
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
      em_alg <- EM_alg_FMOU(output=y, d=est_d_FMOU_here, M=100, threshold=10^{-6},
                            est_U0=T, est_sigma0_2=T)
    })
    rmse_mean_FMOU[it, idx_noise] = sqrt(mean(( em_alg$pred_mean - f_mat)^2))
    L_95[it, idx_noise] = mean(em_alg$pred_mean_95upper-em_alg$pred_mean_95lower)
    P_95[it, idx_noise] = mean(em_alg$pred_mean_95upper>f_mat & em_alg$pred_mean_95lower<f_mat)
    time_FMOU[it,idx_noise] = time[3]

    if(it==1){
      pred_mean_FMOU[[idx_noise]] = em_alg$pred_mean
    }

    ###################### DMD ######################
    time = system.time({
      DMD_fit = DMD_alg(y, threshold=0.99)
    })
    rmse_mean_DMD[it, idx_noise] = sqrt(mean((cbind(y[,1], DMD_fit$in_sample_pred) - f_mat)^2))
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


    rmse_mean_QW1[it, idx_noise] = sqrt(mean(( U_est_QW1%*%t(U_est_QW1)%*%y - f_mat)^2))

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

    rmse_mean_QW5[it, idx_noise] = sqrt(mean(( U_est_QW5%*%t(U_est_QW5)%*%y -f_mat)^2))
    if(it==1){
      pred_mean_QW5[[idx_noise]] = U_est_QW5%*%t(U_est_QW5)%*%y
    }

    # ####################### AR1_EM ######################
    A = diag(k)
    #### estimate parameters by EM
    time = system.time({
      ar1_em_fit = EM_SS(y, A, n_iter=10)
    })
    #### predict mean
    kf = KF_multidim(A, ar1_em_fit$Q, ar1_em_fit$Phi, ar1_em_fit$R, ar1_em_fit$mu, ar1_em_fit$Sigma, y)
    time_ar1_em[it, idx_noise] = time[3]
    rmse_mean_ar1_em[it, idx_noise] = sqrt(mean((f_mat-kf$s[,2:(n+1)])^2))
    if(it==1){
      pred_mean_ar1_em[[idx_noise]] = kf$s[,2:(n+1)]
    }
  }
}

rmse_mean_noise_1 = cbind(rmse_mean_FMOU[,1], rmse_mean_DMD[,1], rmse_mean_QW1[,1], rmse_mean_QW5[,1], rmse_mean_ar1_em[,1])
rmse_mean_noise_2 = cbind(rmse_mean_FMOU[,2], rmse_mean_DMD[,2], rmse_mean_QW1[,2], rmse_mean_QW5[,2], rmse_mean_ar1_em[,2])
rmse_mean_noise_3 = cbind(rmse_mean_FMOU[,3], rmse_mean_DMD[,3], rmse_mean_QW1[,3], rmse_mean_QW5[,3], rmse_mean_ar1_em[,3])

time_noise_1 = cbind(time_FMOU[,1], time_DMD[,1], time_QW1[,1], time_QW5[,1], time_ar1_em[,1])
time_noise_2 = cbind(time_FMOU[,2], time_DMD[,2], time_QW1[,2], time_QW5[,2], time_ar1_em[,2])
time_noise_3 = cbind(time_FMOU[,3], time_DMD[,3], time_QW1[,3], time_QW5[,3], time_ar1_em[,3])

######################### Table -----------------------------------
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
rownames(time_mean) = c("FMOU", "DMD","LY1","LY5","EM-AR(1)")
colnames(time_mean) = sigma_0_2_list
rowMeans(time_mean)


######################### Figure 4 -----------------------------------
mycolor = colorRampPalette(c("blue", "white", "red"))(100)
fmou_fit = EM_alg_FMOU(output=observation[[2]], d=est_d_FMOU[1,2],M=100, threshold=10^{-6},
                       est_U0=T, est_sigma0_2=T, track_iterations = T, neg_log_lik = T)
#png("../figures/branin_pred_fmou_dmd_June25.png", width=1000, height=250)
par(mfrow=c(1,4),mgp = c(1.25, 0.5, 0))
par(mar=c(3,3,3,1.5))
image2D(t(f_mat), main="(e) True mean", xlab=expression(x[2]),ylab=expression(x[1]),xaxt = "n", 
        yaxt = "n",col = mycolor,cex.main=2.2,cex.lab = 2, colkey = list(cex.axis=1.8,cex.clab=1.8))
par(mar=c(3,3,3,1.5))
image2D(t(pred_mean_FMOU[[2]]), main="(f) Pred mean (FMOU)", 
        xlab=expression(x[2]), ylab="",xaxt = "n", yaxt = "n",
        col = mycolor,cex.main=2.2,cex.lab = 2, colkey = list(cex.axis=1.8,cex.clab=1.8))
par(mar=c(3,3,3,2.2))
image2D(t(pred_mean_DMD[[2]]), main="(g) Pred mean (DMD)",
        xlab=expression(x[2]),ylab="", xaxt = "n", yaxt = "n",col = mycolor,
        cex.main=2.2,cex.lab = 2, colkey = list(cex.axis=1.8,cex.clab=1.8))
par(mar=c(4,5,3,2.5))
par(mgp=c(2.5, 1, 0))  # increase the distance between lab and ticks
plot(-fmou_fit$record_neg_lik,type='l',
     ylab='log likelihood',xlab='E-M iteration',main='(h) Log lik, FMOU',cex.main=2.2, 
     cex.lab=1.8,cex.axis=1.4)
#dev.off()