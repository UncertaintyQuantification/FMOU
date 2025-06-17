library(rstiefel)
library(FastGaSP)
library(Rcpp)
library(RcppEigen)
library(rospca)
library(ggplot2)
library(gridExtra)
library(cowplot)

source('../functions/FMOU.R')
source('../functions/DMD.R')
source('../functions/functions_geophysics_gppca.R')
sourceCpp(file='../functions/functions.cpp') 


N = 20 # repeat times
d_vec = c(5, 10)
k_vec = c(20, 40) # number of stations
n1 = 100 # number of epochs
n2 = 200
n3 = 400 
noise_level = 0.2 


### n = n1
angle_diff_FMOU_n1 = matrix(NA, N, length(d_vec)) 
angle_diff_GPPCA_n1 = matrix(NA, N, length(d_vec)) 
angle_diff_FMOU_n1_est_d = matrix(NA, N, length(d_vec)) 

rmse_mean_FMOU_n1 = matrix(NA, N, length(d_vec)) 
rmse_mean_GPPCA_n1 = matrix(NA, N, length(d_vec)) 

time_record_FMOU_n1 = matrix(NA, N, length(d_vec)) 
time_record_GPPCA_n1 = matrix(NA, N, length(d_vec)) 

for(idx_d in 1: length(d_vec)){
  d = d_vec[idx_d]
  k = k_vec[idx_d]
  for(it in 1:N){
    set.seed(it)
    if(it%%10 == 0){print(it)}
    ## generate latent processes
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
    signal = U[,1:d] %*% z
    y = signal + matrix(rnorm(n1*k,mean=0,sd=sqrt(noise_level)),k,n1)
    
    
    ####################### Start FMOU with true d ######################
    time1 = Sys.time()
    em_alg <- EM_alg_FMOU(output=y, d=d, threshold=10^{-8}, M=100,
                          est_U0=T, U0=NA, est_sigma0_2=T)
    U_est_FMOU = em_alg$U
    time2 = Sys.time()
    time_record_FMOU_n1[it, idx_d] = time2 - time1

    fit_latent_state_EM = U_est_FMOU %*% em_alg$Z_hat
   
    angle_diff_FMOU_n1[it, idx_d] = rospca::angle(U[,1:d], U_est_FMOU)
    rmse_mean_FMOU_n1[it, idx_d] = sqrt(mean(( fit_latent_state_EM - signal)^2))
    
    
    ###################### GPPCA ######################
    G_hat=as.list(1:d)
    kernel_type='exp'
    ### these are required by GPPCA
    input = seq(1:ncol(y))
    output = y
    delta_x= rep(1, ncol(y)-1)
    n = n1
    ###
    time1 = Sys.time()
    output_2=y%*%t(y)
    trace_output_2=sum(y^2)
    svd_output=svd(y)
    A_ini=svd_output$u[,1:d]
    param_ini=c(rep(log(.1),d),rep(log(1),d))
    m=try(optim(param_ini,neg_log_lik_diff_cov_FFBS, A_ini=A_ini, kernel_type=kernel_type, method="L-BFGS-B",lower=c(rep(log(0.001),d),rep(log(.1),d)),
                upper=c(rep(log(7),d),rep(log(1000),d)), control=list(maxit=50)), silent=T)
    # m=try(optim(param_ini,neg_log_lik_diff_cov_FFBS, A_ini=A_ini, kernel_type=kernel_type, method="L-BFGS-B",
    #             control=list(maxit=100)), silent=T)
    beta_hat=exp(m$par[1:d])
    tau_hat=exp(m$par[(d+1):(2*d)])
    G_log_det_cov_hat=Get_G_log_det_cov(m$par, output, delta_x,d, kernel_type = kernel_type)
    for(i in 1:d){
      G_hat[[i]]=output_2-G_log_det_cov_hat[[i]]
    }
    A_hat=Optimization_Stiefel_Manifold(A_ini, G=G_hat, max_iter=100)
    sigma_2_0_hat=(trace_output_2+F_Funct(A_hat,G_hat))/(length(output))
    
    record_param = rep(NA, 2*d + 1)
    record_param[1:d]=beta_hat
    record_param[(d+1):(2*d)]=tau_hat*(sigma_2_0_hat)
    record_param[2*d+1]=sigma_2_0_hat
    
    AZ_posterior=pred_FFBS_FastGaSP(param=record_param,A_hat=A_hat,input=input,testing_input=input, output_here=output,d=d,kernel_type=kernel_type)
    time2<-Sys.time()
    time_record_GPPCA_n1[it, idx_d] = as.numeric(difftime(time2, time1, units="secs"))
    
    rmse_mean_GPPCA_n1[it, idx_d] = sqrt(mean(( AZ_posterior[[1]] - signal)^2))
    angle_diff_GPPCA_n1[it, idx_d] = rospca::angle(U[,1:d], A_hat)
    
  }
}


angle_diff_n1_setup_1 = cbind(angle_diff_FMOU_n1[,1], angle_diff_GPPCA_n1[,1])
angle_diff_n1_setup_2 = cbind(angle_diff_FMOU_n1[,2], angle_diff_GPPCA_n1[,2])
rmse_mean_n1_setup_1 = cbind(rmse_mean_FMOU_n1[,1], rmse_mean_GPPCA_n1[,1])
rmse_mean_n1_setup_2 = cbind(rmse_mean_FMOU_n1[,2], rmse_mean_GPPCA_n1[,2])

time_n1_setup_1 = cbind(time_record_FMOU_n1[,1], time_record_GPPCA_n1[,1])
time_n1_setup_2 = cbind(time_record_FMOU_n1[,2], time_record_GPPCA_n1[,2])




### n = n2
angle_diff_FMOU_n2 = matrix(NA, N, length(d_vec)) 
angle_diff_GPPCA_n2 = matrix(NA, N, length(d_vec)) 
angle_diff_FMOU_n2_est_d = matrix(NA, N, length(d_vec)) 

rmse_mean_FMOU_n2 = matrix(NA, N, length(d_vec)) 
rmse_mean_GPPCA_n2 = matrix(NA, N, length(d_vec)) 

time_record_FMOU_n2 = matrix(NA, N, length(d_vec)) 
time_record_GPPCA_n2 = matrix(NA, N, length(d_vec)) 

for(idx_d in 1: length(d_vec)){
  d = d_vec[idx_d]
  k = k_vec[idx_d]
  for(it in 1:N){
    set.seed(2*it)
    if(it%%10 == 0){print(it)}
    ## generate latent processes
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
    signal = U[,1:d] %*% z
    y = signal + matrix(rnorm(n2*k,mean=0,sd=sqrt(noise_level)),k,n2)
    
    
    ####################### Start FMOU with true d ######################
    time1 = Sys.time()
    em_alg <- EM_alg_FMOU(output=y, d=d, threshold=10^{-8}, M=100,
                          est_U0=T, U0=NA, est_sigma0_2=T)
    U_est_FMOU = em_alg$U
    time2 = Sys.time()
    time_record_FMOU_n2[it, idx_d] = time2 - time1
    
    fit_latent_state_EM = U_est_FMOU %*% em_alg$Z_hat
    
    angle_diff_FMOU_n2[it, idx_d] = rospca::angle(U[,1:d], U_est_FMOU)
    rmse_mean_FMOU_n2[it, idx_d] = sqrt(mean(( fit_latent_state_EM - signal)^2))
    
    
    ###################### GPPCA ######################
    G_hat=as.list(1:d)
    kernel_type='exp'
    ### these are required by GPPCA
    input = seq(1:ncol(y))
    output = y
    delta_x= rep(1, ncol(y)-1)
    n = n2
    ###
    time1 = Sys.time()
    output_2=y%*%t(y)
    trace_output_2=sum(y^2)
    svd_output=svd(y)
    A_ini=svd_output$u[,1:d]
    param_ini=c(rep(log(.1),d),rep(log(1),d))
    m=try(optim(param_ini,neg_log_lik_diff_cov_FFBS, A_ini=A_ini, kernel_type=kernel_type, method="L-BFGS-B",lower=c(rep(log(0.001),d),rep(log(.1),d)),
                upper=c(rep(log(7),d),rep(log(1000),d)), control=list(maxit=50)), silent=T)
    # m=try(optim(param_ini,neg_log_lik_diff_cov_FFBS, A_ini=A_ini, kernel_type=kernel_type, method="L-BFGS-B",
    #             control=list(maxit=100)), silent=T)
    beta_hat=exp(m$par[1:d])
    tau_hat=exp(m$par[(d+1):(2*d)])
    G_log_det_cov_hat=Get_G_log_det_cov(m$par, output, delta_x,d, kernel_type = kernel_type)
    for(i in 1:d){
      G_hat[[i]]=output_2-G_log_det_cov_hat[[i]]
    }
    A_hat=Optimization_Stiefel_Manifold(A_ini, G=G_hat, max_iter=100)
    sigma_2_0_hat=(trace_output_2+F_Funct(A_hat,G_hat))/(length(output))
    
    record_param = rep(NA, 2*d + 1)
    record_param[1:d]=beta_hat
    record_param[(d+1):(2*d)]=tau_hat*(sigma_2_0_hat)
    record_param[2*d+1]=sigma_2_0_hat
    
    AZ_posterior=pred_FFBS_FastGaSP(param=record_param,A_hat=A_hat,input=input,testing_input=input, output_here=output,d=d,kernel_type=kernel_type)
    time2<-Sys.time()
    time_record_GPPCA_n2[it, idx_d] = as.numeric(difftime(time2, time1, units="secs"))
    
    rmse_mean_GPPCA_n2[it, idx_d] = sqrt(mean(( AZ_posterior[[1]] - signal)^2))
    angle_diff_GPPCA_n2[it, idx_d] = rospca::angle(U[,1:d], A_hat)
    
  }
}



angle_diff_n2_setup_1 = cbind(angle_diff_FMOU_n2[,1], angle_diff_GPPCA_n2[,1])
angle_diff_n2_setup_2 = cbind(angle_diff_FMOU_n2[,2], angle_diff_GPPCA_n2[,2])
rmse_mean_n2_setup_1 = cbind(rmse_mean_FMOU_n2[,1], rmse_mean_GPPCA_n2[,1])
rmse_mean_n2_setup_2 = cbind(rmse_mean_FMOU_n2[,2], rmse_mean_GPPCA_n2[,2])

time_n2_setup_1 = cbind(time_record_FMOU_n2[,1], time_record_GPPCA_n2[,1])
time_n2_setup_2 = cbind(time_record_FMOU_n2[,2], time_record_GPPCA_n2[,2])



### n = n3
angle_diff_FMOU_n3 = matrix(NA, N, length(d_vec)) 
angle_diff_GPPCA_n3 = matrix(NA, N, length(d_vec)) 
angle_diff_FMOU_n3_est_d = matrix(NA, N, length(d_vec)) 

rmse_mean_FMOU_n3 = matrix(NA, N, length(d_vec)) 
rmse_mean_GPPCA_n3 = matrix(NA, N, length(d_vec)) 

time_record_FMOU_n3 = matrix(NA, N, length(d_vec)) 
time_record_GPPCA_n3 = matrix(NA, N, length(d_vec)) 

for(idx_d in 1: length(d_vec)){
  d = d_vec[idx_d]
  k = k_vec[idx_d]
  for(it in 1:N){
    set.seed(3*it)
    if(it%%10 == 0){print(it)}
    ## generate latent processes
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
    signal = U[,1:d] %*% z
    y = signal + matrix(rnorm(n3*k,mean=0,sd=sqrt(noise_level)),k,n3)
    
    
    ####################### Start FMOU with true d ######################
    time1 = Sys.time()
    em_alg <- EM_alg_FMOU(output=y, d=d, threshold=10^{-8}, M=100,
                          est_U0=T, U0=NA, est_sigma0_2=T)
    U_est_FMOU = em_alg$U
    time2 = Sys.time()
    time_record_FMOU_n3[it, idx_d] = time2 - time1
    
    fit_latent_state_EM = U_est_FMOU %*% em_alg$Z_hat
    
    angle_diff_FMOU_n3[it, idx_d] = rospca::angle(U[,1:d], U_est_FMOU)
    rmse_mean_FMOU_n3[it, idx_d] = sqrt(mean(( fit_latent_state_EM - signal)^2))
    
    
    ###################### GPPCA ######################
    G_hat=as.list(1:d)
    kernel_type='exp'
    ### these are required by GPPCA
    input = seq(1:ncol(y))
    output = y
    delta_x= rep(1, ncol(y)-1)
    n = n3
    ###
    time1 = Sys.time()
    output_2=y%*%t(y)
    trace_output_2=sum(y^2)
    svd_output=svd(y)
    A_ini=svd_output$u[,1:d]
    param_ini=c(rep(log(.1),d),rep(log(1),d))
    m=try(optim(param_ini,neg_log_lik_diff_cov_FFBS, A_ini=A_ini, kernel_type=kernel_type, method="L-BFGS-B",lower=c(rep(log(0.001),d),rep(log(.1),d)),
                upper=c(rep(log(7),d),rep(log(1000),d)), control=list(maxit=50)), silent=T)
    # m=try(optim(param_ini,neg_log_lik_diff_cov_FFBS, A_ini=A_ini, kernel_type=kernel_type, method="L-BFGS-B",
    #             control=list(maxit=100)), silent=T)
    beta_hat=exp(m$par[1:d])
    tau_hat=exp(m$par[(d+1):(2*d)])
    G_log_det_cov_hat=Get_G_log_det_cov(m$par, output, delta_x,d, kernel_type = kernel_type)
    for(i in 1:d){
      G_hat[[i]]=output_2-G_log_det_cov_hat[[i]]
    }
    A_hat=Optimization_Stiefel_Manifold(A_ini, G=G_hat, max_iter=100)
    sigma_2_0_hat=(trace_output_2+F_Funct(A_hat,G_hat))/(length(output))
    
    record_param = rep(NA, 2*d + 1)
    record_param[1:d]=beta_hat
    record_param[(d+1):(2*d)]=tau_hat*(sigma_2_0_hat)
    record_param[2*d+1]=sigma_2_0_hat
    
    AZ_posterior=pred_FFBS_FastGaSP(param=record_param,A_hat=A_hat,input=input,testing_input=input, output_here=output,d=d,kernel_type=kernel_type)
    time2<-Sys.time()
    time_record_GPPCA_n3[it, idx_d] = as.numeric(difftime(time2, time1, units="secs"))
    
    rmse_mean_GPPCA_n3[it, idx_d] = sqrt(mean(( AZ_posterior[[1]] - signal)^2))
    angle_diff_GPPCA_n3[it, idx_d] = rospca::angle(U[,1:d], A_hat)
    
  }
}


angle_diff_n3_setup_1 = cbind(angle_diff_FMOU_n3[,1], angle_diff_GPPCA_n3[,1])
angle_diff_setup_1 = cbind(angle_diff_n1_setup_1, angle_diff_n2_setup_1, angle_diff_n3_setup_1)

angle_diff_n3_setup_2 = cbind(angle_diff_FMOU_n3[,2], angle_diff_GPPCA_n3[,2])
angle_diff_setup_2 = cbind(angle_diff_n1_setup_2, angle_diff_n2_setup_2, angle_diff_n3_setup_2)

rmse_mean_n3_setup_1 = cbind(rmse_mean_FMOU_n3[,1], rmse_mean_GPPCA_n3[,1])
rmse_mean_setup_1 = cbind(rmse_mean_n1_setup_1, rmse_mean_n2_setup_1, rmse_mean_n3_setup_1)

rmse_mean_n3_setup_2 = cbind(rmse_mean_FMOU_n3[,2], rmse_mean_GPPCA_n3[,2])
rmse_mean_setup_2 = cbind(rmse_mean_n1_setup_2, rmse_mean_n2_setup_2, rmse_mean_n3_setup_2)

time_n3_setup_1 = cbind(time_record_FMOU_n3[,1], time_record_GPPCA_n3[,1])
time_all_result_setup_1 = cbind(time_n1_setup_1, time_n2_setup_1, time_n3_setup_1)

time_n3_setup_2 = cbind(time_record_FMOU_n3[,2], time_record_GPPCA_n3[,2])
time_all_result_setup_2 = cbind(time_n1_setup_2, time_n2_setup_2, time_n3_setup_2)



avg_time_setup_1 = matrix( apply(time_all_result_setup_1, 2, mean), ncol=3)
avg_time_setup_2 = matrix( apply(time_all_result_setup_2, 2, mean), ncol=3)

######################### Figure S1 #########################
# pdf('Exp1_supplement_results.pdf',height=6.5, width=10)
n_values <- c(n1, n2, n3) 
# 1st Row: 
data <- as.data.frame(angle_diff_setup_1) * (pi / 2)
data_long <- data.frame(
  value = as.vector(as.matrix(data)),
  method = rep(c("FMOU", "GPPCA"), times = 3, each = nrow(data)),
  position = rep(c(1, 2, 4, 5, 7, 8), each = nrow(data))
)
boxplot1 <- ggplot(data_long, aes(x = factor(position), y = value, fill = method)) +
  geom_boxplot(outlier.shape =1, outlier.size = 1) + theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  scale_fill_manual(values =c("#00a6fb", "#bf963c")) +
  ylim(0, 0.225) +
  labs(x = "", y = "Largest principal angle") +
  theme_minimal()+
  #theme(legend.position = "none")+
  scale_x_discrete(labels = c("FMOU", "GPPCA", "FMOU", "GPPCA", "FMOU", "GPPCA")) +
  theme(legend.position = c(0.68, 0.9), legend.text = element_text(size = 5),
        legend.background = element_rect(fill = "white", color = "black"), legend.direction="horizontal",
        legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )

data <- as.data.frame(rmse_mean_setup_1)
data_long <- data.frame(
  value = as.vector(as.matrix(data)),
  method = rep(c("FMOU", "GPPCA", "FMOU", "GPPCA", "FMOU", "GPPCA"), each = nrow(data)),
  position = rep(c(1, 2, 4, 5, 7, 8), each = nrow(data))
)
boxplot2 <- ggplot(data_long, aes(x =factor(position), y = value, fill = method)) +
  geom_boxplot(outlier.shape =1, outlier.size = 1) + theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  scale_fill_manual(values = c("#00a6fb", "#bf963c")) +
  ylim(0.15, 0.35) +
  labs(x = "", y = "RMSE of mean") +
  theme_minimal()+
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("FMOU", "GPPCA", "FMOU", "GPPCA", "FMOU", "GPPCA")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )

bar_data_setup_1 <- data.frame(
  n = rep(n_values, times = 2),
  time = as.vector(t(avg_time_setup_1[1:2, ])),
  method = rep(c("FMOU", "GPPCA"), each = length(n_values))
)
barplot1 <- ggplot(bar_data_setup_1, aes(x = factor(n), y = time, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") + theme(plot.margin = unit(c(1, 1, 2.5, 1), "cm")) +
  geom_text(aes(label = round(time, 2)), position = position_dodge(width = 0.9), vjust = -0.5) + 
  scale_fill_manual(values = c("FMOU" = "#00a6fb", "GPPCA" = "#bf963c")) +
  ylim(0, 26)  + 
  labs(x = "n", y = "Running time (seconds)", fill = "") +
  theme_minimal() + 
  theme(legend.position = c(0.32, 0.9), legend.text = element_text(size = 5),
        legend.background = element_rect(fill = "white", color = "black"), legend.direction="horizontal",
        legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )


# 2nd Row:
data <- as.data.frame(angle_diff_setup_2) * (pi / 2)
data_long <- data.frame(
  value = as.vector(as.matrix(data)),
  method = rep(c("FMOU", "GPPCA", "FMOU", "GPPCA", "FMOU", "GPPCA"), each = nrow(data)),
  position = rep(c(1, 2, 4, 5, 7, 8), each = nrow(data))
)
boxplot3 <- ggplot(data_long, aes(x =factor(position), y = value, fill = method)) +
  geom_boxplot(outlier.shape =1, outlier.size = 1) + theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  scale_fill_manual(values = c("#00a6fb", "#bf963c")) +
  ylim(0, 0.8) +
  labs(x = "", y = "Largest principal angle") +
  theme_minimal()+
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("FMOU", "GPPCA", "FMOU", "GPPCA", "FMOU", "GPPCA")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )
# theme(legend.position = c(0.1, 0.9), # Move legend inside the plot (80% right, 80% up)
#       legend.background = element_rect(fill = "white", color = "black"), # Optional: Add a background box
#       legend.title = element_blank())

data <- as.data.frame(rmse_mean_setup_2)
data_long <- data.frame(
  value = as.vector(as.matrix(data)),
  method = rep(c("FMOU", "GPPCA", "FMOU", "GPPCA", "FMOU", "GPPCA"), each = nrow(data)),
  position = rep(c(1, 2, 4, 5, 7, 8), each = nrow(data))
)
boxplot4 <- ggplot(data_long, aes(x = factor(position), y = value, fill = method)) +
  geom_boxplot(outlier.shape =1, outlier.size = 1) + theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))+
  scale_fill_manual(values = c("#00a6fb", "#bf963c")) +
  ylim(0.15, 0.5) +
  labs(x = "", y = "RMSE of mean") +
  theme_minimal()+
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("FMOU", "GPPCA", "FMOU", "GPPCA", "FMOU", "GPPCA")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )

bar_data_setup_2 <- data.frame(
  n = rep(n_values, times = 2),
  time = as.vector(t(avg_time_setup_2[1:2, ])),
  method = rep(c("FMOU", "GPPCA"), each = length(n_values))
)
barplot2 <- ggplot(bar_data_setup_2, aes(x = factor(n), y = time, fill = method)) +
  geom_text(aes(label = round(time, 2)), position = position_dodge(width = 0.9), vjust = -0.5) + 
  geom_bar(stat = "identity", position = "dodge") + theme(plot.margin = unit(c(1, 1, 2.5, 1), "cm")) +
  scale_fill_manual(values = c("FMOU" = "#00a6fb", "GPPCA" = "#bf963c")) +
  ylim(0, 165)  + 
  labs(x = "n", y = "Running time (seconds)", fill = "") +
  theme_minimal() + 
  theme(legend.position = "none") +
  # theme(legend.position = c(0.22, 0.85), legend.text = element_text(size = 3.8),
  #       legend.background = element_rect(fill = "white", color = "black"), 
  #       legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )

# Define a uniform theme to apply to all plots
uniform_theme <- theme(
  legend.text = element_text(size = 8.3),  # Legend text size
  legend.title = element_blank(),        # Remove legend title
  axis.text = element_text(size = 8.5),   # Axis tick labels size
  axis.title = element_text(size = 10.6),  # Axis labels size
  plot.title = element_text(size = 10.5)   # Plot title size (if applicable)
)

boxplot1 <- boxplot1 + uniform_theme
boxplot2 <- boxplot2 + uniform_theme
boxplot3 <- boxplot3 + uniform_theme
boxplot4 <- boxplot4 + uniform_theme
barplot1 <- barplot1 + uniform_theme
barplot2 <- barplot2 + uniform_theme

aligned_plots <- plot_grid(
  boxplot1, boxplot2, barplot1, boxplot3, boxplot4, barplot2,
  ncol = 3, align = "hv", axis = "tb" 
)

aligned_plots
# dev.off()

