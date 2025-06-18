library(geometry)
library(plot3D)
library(MASS)
library(FastGaSP)
library(RobustGaSP)
library(pracma)

# functions used for NIF and FMOU
source('../functions/FMOU.R')
source('../functions/functions_geophysics_gppca.R')

## Read data and results summarized from Noel's MATLAB source
# load estimated slips (after backward smoothing) by NIF
slips = read.csv("../real_data/Cascadia_data/est_slip_large_mesh.csv", header=FALSE)
# load Greens functions
GF = read.csv("../real_data/Cascadia_data/Green_fun_large_mesh.csv", header=FALSE)
# load observations in 3 directions over 88 days (in meters)
obs = read.csv("../real_data/Cascadia_data/obs_KF.csv", header=FALSE)
# load estimated BM by NIF
bench_BM = read.csv("../real_data/Cascadia_data/est_BM_large_mesh.csv", header=FALSE)
# load estimated frame motions by NIF
ref_frame = read.csv("../real_data/Cascadia_data/est_frame_large_mesh.csv", header=FALSE)
# reads date (epochs) of observations
doy = read.csv("../real_data/Cascadia_data/DOY.csv", header=FALSE)
# the forward results of KF (state vectors) in NIF
NIF_forward = read.csv("../real_data/Cascadia_data/NIF_forward_large_mesh.csv", header=FALSE)
# the backward results of KF (state vectors) in NIF
NIF_backward = read.csv("../real_data/Cascadia_data/NIF_backward_large_mesh.csv", header=FALSE)
# obsrvation uncertainties (in meters)
obs_uncertainty = read.csv("../real_data/Cascadia_data/obs_uncertainty.csv", header=FALSE)
# get the name of stations we used in analysis:
names = read.csv("../real_data/Cascadia_data/station_names.csv", header=FALSE)

slips = as.matrix(slips)
GF = GF_ori = as.matrix(GF)
obs = as.matrix(obs)
bench_BM = as.matrix(bench_BM)
ref_frame = as.matrix(ref_frame)
NIF_forward = as.matrix(NIF_forward)
NIF_backward = as.matrix(NIF_backward)
names = as.vector(as.matrix(names))
doy = as.numeric(doy)
obs_uncertainty_ori = obs_uncertainty = as.matrix(obs_uncertainty)

epoch = dim(slips)[2]
N_site = dim(obs)[1] / 3 # number of GPS stations
N_xi = dim(GF)[2] # number of fault patches

# Stations in the index list: (5  16  17  24  30  80  83  88 104 105 107 109 110 111 112 113 114) are missing
missing_location = c(5,  16,  17,  24,  30,  80,  83,  88, 104, 105, 107, 109, 110, 111, 112, 113, 114)
# We only model the measurments in East and North
N_site_used = N_site - length(missing_location)
GF = GF[-c(3*missing_location, 3*missing_location-1, 3*missing_location-2), ]
GF= GF[-seq(3, 3*N_site_used, by=3),]
bench_BM = bench_BM[-c(3*missing_location, 3*missing_location-1, 3*missing_location-2), ]
bench_BM = bench_BM[-seq(3, 3*N_site_used, by=3),]
ref_frame = ref_frame[-c(3*missing_location, 3*missing_location-1, 3*missing_location-2), ]
ref_frame = ref_frame[-seq(3, 3*N_site_used, by=3),]
observed_location = (1:N_site)[-missing_location]
obs_uncertainty = obs_uncertainty[c(mapply(c, 3*observed_location-2, 3*observed_location-1)), ]

sum(is.na(obs[-c(3*missing_location, 3*missing_location-1, 3*missing_location-2), ])) / (300*88)
# About 15% of the observations are missing


plot_data_idx = c(16, 33, 54, 67, 90, 119, 150, 193)
plot_slip_idx = c(50, 100, 350, 500, 756, 850, 921, 1100, 1250, 1500, 1700, 1800, 1900, 1950)


obs_used = obs
# fill the missing parts at each station:
t = c(1: epoch)
for(i in 1:length(observed_location)){
  index_here = c(3*observed_location[i]-2, 3*observed_location[i]-1, 3*observed_location[i])
  for(j in 1:3){
    if(sum(obs[index_here[j], ] == "NaN")>0){
      NA_index = as.numeric(which(obs[index_here[j], ] == "NaN"))
      input_here = t[-NA_index]
      output_here = as.numeric(obs[index_here[j], - NA_index])
      m_rgasp = rgasp(input_here, output_here, zero.mean="Yes", nugget.est=T, kernel_type = "pow_exp", alpha=c(1))
      
      ## fill the missing values:
      if (length(NA_index)>1){
        pred.model = predict(m_rgasp, testing_input=as.numeric(t[NA_index]))
        set.seed(i*j)
        obs_used[index_here[j], NA_index] = rnorm(length(pred.model$mean), mean=pred.model$mean, sd=pred.model$sd)
      }
      else{ 
        pred.model=predict(m_rgasp, testing_input=as.numeric(c(t[1], t[NA_index])))
        obs_used[index_here[j], NA_index] = rnorm(1, mean=pred.model$mean[2], sd=pred.model$sd[2])
      }
    }
  }
}

obs_cleaned = obs_used[c(mapply(c, 3*observed_location-2, 3*observed_location-1)), ] 
sum(is.na(obs_cleaned))
# Make the initial value of observations as zero (then cummulative displacements will be zero at t_1)
output = obs_cleaned  - obs_cleaned[, 1] 

# Try to add frame motions, where the displacements at all stations in each direction 
# is the same. It's equivalinet to model Y = G*s + (f_1, f_2)_{300 x 88} + \epsilon
G = cbind(GF, rep(c(1,0), 100), rep(c(0,1), 100)) # add 2 columns for East/North frame motions
dim(G)


eigen_Phi=eigen(G%*%t(G))
Lambda = eigen_Phi$values

# We can use the observation uncertainty provided to approximate the errors approximately
obs_uncertainty_E = mean(obs_uncertainty[seq(1,200,2), ], na.rm=TRUE) # 4.151058e-06
obs_uncertainty_N = mean(obs_uncertainty[seq(2,200,2), ], na.rm=TRUE) # 4.205995e-06
avg_obs_uncertainty = 0.5*mean(obs_uncertainty[seq(1,200,2), ], na.rm=TRUE) + 0.5*mean(obs_uncertainty[seq(2,200,2), ], na.rm=TRUE) # 4.178526e-06

obs_uncertainty_cleaned = obs_uncertainty
# fill the missing parts at each station:
t = c(1: epoch)
for(j in 1:nrow( obs_uncertainty)){
  if(sum(obs_uncertainty[j, ] == "NaN")>0){
    NA_index = as.numeric(which(obs_uncertainty[j, ] == "NaN"))
    if(j %% 2 == 1){
      obs_uncertainty_cleaned[j, NA_index] = obs_uncertainty_E
    }
    else{
      obs_uncertainty_cleaned[j, NA_index] = obs_uncertainty_N
    }
  }
}
sum(is.na(obs_uncertainty_cleaned))

U = eigen_Phi$vectors
Psi = t(U)%*%G    # compute basis functions

# And we select "d" based on observations' uncertainties
d_list = seq(10, 100)
#d_list = seq(1,dim(U)[1])
sd_slip = rep(NA, length(d_list))
noise_est = rep(NA, length(d_list))


time1 = Sys.time()
left = d_list[1]
right = rev(d_list)[1]
while(left <= right){
  mid = left + floor((right-left)/2)
  
  em_fit = EM_alg_FMOU(output, d=mid, est_U0=F, U0=U[,1:mid], est_sigma0_2=T, threshold=10^{-6}) 
  est_noise = em_fit$sigma2_0
  if(est_noise < avg_obs_uncertainty){
    right = mid-1
  }
  else{
    left = mid + 1
  }
}
Nbasis_EM  = mid
time2 = Sys.time()
time_FMOU_select_d = time2 - time1
time_FMOU_select_d 



time1 = Sys.time()
tilde_output = t(U[,1:Nbasis_EM])%*%output
system.time({
  em_alg = EM_alg_FMOU(output, d=Nbasis_EM, est_U0=F, est_sigma0_2=F, threshold = 10^{-6},
                       U0=U[,1:Nbasis_EM], M=100, sigma0_2 = avg_obs_uncertainty)
  Sp_latent_state_EM = t(Psi[1:Nbasis_EM,])%*%diag(1/Lambda[1:Nbasis_EM])%*%em_alg$Z_hat
  fit_latent_state_EM = U[,1:Nbasis_EM]%*%em_alg$Z_hat
})[3]
time2 = Sys.time()
fmou_time = time2-time1
fmou_time




################## Let's use the original NIF ########################
time1 = Sys.time()
Nbasis = 200
statedim = 3*Nbasis # dimension of latent states for Kalman Filter
G_NIF = matrix(U%*%diag(Lambda)[,1:Nbasis], ncol=Nbasis)
theta_1 =  logspace(2, 4, 5)
theta_2 = c(0.5, 1, 1.5, 2)
loglik_theta_combination = matrix(NA, length(theta_1), length(theta_2))
for(i1 in 1:length(theta_1)){
  for(i2 in 1:length(theta_2)){
    param_here = c(theta_1[i1], theta_2[i2])
    loglik_theta_combination[i1, i2] = neg_loglik_fix_noise(param_here, t=doy, Data=output, G_NIF, Lambda, avg_obs_uncertainty)
  }
}
optim_loglik_index = which.min(loglik_theta_combination)
theta_2_index = ceiling(optim_loglik_index / length(theta_1))
theta_1_index = optim_loglik_index - (theta_2_index - 1)*length(theta_1)
theta_hat = theta_1[theta_1_index]
theta2_hat = theta_2[theta_2_index]


tau = 0
res = netloglik3p(t=doy, Data=output, G_NIF, Lambda, theta_hat, theta2_hat)
sig2_hat = avg_obs_uncertainty
alpha_hat = sqrt(sig2_hat)*theta_hat
gamma_hat = sqrt(sig2_hat)*theta2_hat

# construct the initial states
x0 = matrix(0, statedim, 1)
var0 = 1.0^(-8)*Lambda[1]*diag(statedim)
for(j in 1:Nbasis)
{var0[(3*(j-1)+1),(3*(j-1)+1)] = gamma_hat^2 *Lambda[j]}

# Run Filter
res =  netfilt(t=doy, Data=output, G=matrix(G_NIF[,1:Nbasis], ncol=Nbasis), sigmas=c(tau, sqrt(sig2_hat), alpha_hat), x0, var0, 1)
x_kgk=res[[1]];  sigma_kgk = res[[2]]
x_kgn = res[[3]]; sigma_kgn= res[[4]]

c_NIF = matrix(0,Nbasis, length(doy))  # matrix of stochastic coefficients
Ts = matrix(0, Nbasis, 3*Nbasis)  # we can think it's "G" in KF 

for(h in 1:length(doy))
  for(i in 1:Nbasis)
  { Ts[i,(3*(i-1)+1):(3*(i-1)+3)] = c(doy[h]-doy[1],1,0) 
  c_NIF[,h] = Ts%*%(x_kgn[h,])
  }
NIF_slip = matrix(t(Psi[1:Nbasis,]), ncol=Nbasis)%*%c_NIF  # estimated slips
NIF_signal = G_NIF %*% c_NIF  
time2 = Sys.time()
NIF_time = time2-time1
NIF_time
#time2 - time1 

rm(sigma_kgk)
rm(sigma_kgn)
rm(GF_ori)
rm(res)
rm(NIF_backward)
rm(NIF_forward)
save.image(file="FMOU_real_data.RData")


################## Let's use the original NIF with reduced dimension ########################
# time1 = Sys.time()
# Nbasis = Nbasis_EM
# statedim = 3*Nbasis # dimension of latent states for Kalman Filter
# G_NIF = matrix(U%*%diag(Lambda)[,1:Nbasis], ncol=Nbasis)
# theta_1 =  logspace(2, 4, 5)
# theta_2 = c(0.5, 1, 1.5, 2)
# loglik_theta_combination = matrix(NA, length(theta_1), length(theta_2))
# for(i1 in 1:length(theta_1)){
#   for(i2 in 1:length(theta_2)){
#     param_here = c(theta_1[i1], theta_2[i2])
#     loglik_theta_combination[i1, i2] = neg_loglik_fix_noise(param_here, t=doy, Data=output, G_NIF, Lambda, avg_obs_uncertainty)
#   }
# }
# optim_loglik_index = which.min(loglik_theta_combination)
# theta_2_index = ceiling(optim_loglik_index / length(theta_1))
# theta_1_index = optim_loglik_index - (theta_2_index - 1)*length(theta_1)
# theta_hat = theta_1[theta_1_index]
# theta2_hat = theta_2[theta_2_index]
# 
# 
# tau = 0
# res = netloglik3p(t=doy, Data=output, G_NIF, Lambda, theta_hat, theta2_hat)
# sig2_hat = avg_obs_uncertainty
# alpha_hat = sqrt(sig2_hat)*theta_hat
# gamma_hat = sqrt(sig2_hat)*theta2_hat
# 
# # construct the initial states
# x0 = matrix(0, statedim, 1)
# var0 = 1.0^(-8)*Lambda[1]*diag(statedim)
# for(j in 1:Nbasis)
# {var0[(3*(j-1)+1),(3*(j-1)+1)] = gamma_hat^2 *Lambda[j]}
# 
# # Run Filter
# res_select_d =  netfilt(t=doy, Data=output, G=matrix(G_NIF[,1:Nbasis], ncol=Nbasis), sigmas=c(tau, sqrt(sig2_hat), alpha_hat), x0, var0, 1)
# x_kgk_select_d = res_select_d[[1]];  sigma_kgk_select_d = res_select_d[[2]]
# x_kgn_select_d = res_select_d[[3]]; sigma_kgn_select_d = res_select_d[[4]]
# 
# c_NIF_select_d = matrix(0,Nbasis, length(doy))  # matrix of stochastic coefficients
# Ts_select_d = matrix(0, Nbasis, 3*Nbasis)  # we can think it's "G" in KF 
# 
# for(h in 1:length(doy))
#   for(i in 1:Nbasis)
#   { Ts_select_d[i,(3*(i-1)+1):(3*(i-1)+3)] = c(doy[h]-doy[1],1,0) 
#   c_NIF_select_d[,h] = Ts_select_d%*%(x_kgn_select_d[h,])
#   }
# NIF_slip_select_d = matrix(t(Psi[1:Nbasis,]), ncol=Nbasis)%*%c_NIF_select_d  # estimated slips
# NIF_signal_select_d = G_NIF %*% c_NIF_select_d  
# time2 = Sys.time()
# time2 - time1 # 89.1813 seconds


######################################################################

#raw_obs_drop_up = obs[c(mapply(c, 3*observed_location-2, 3*observed_location-1)), ]
# compare RMSE of fitting data
#sqrt(mean((raw_obs_drop_up - (NIF_signal + obs_cleaned[,1]))^2, na.rm=T)) # NIF
#sqrt(mean((raw_obs_drop_up - (NIF_signal_select_d + obs_cleaned[,1]))^2, na.rm=T)) # NIF selected d
#sqrt(mean((raw_obs_drop_up - (fit_latent_state_EM + obs_cleaned[,1]))^2, na.rm=T)) # FMOU
# percentage of negative slips
#sum(NIF_slip < 0) / (1978*88) # NIF
#sum(NIF_slip_select_d < 0) / (1978*88) # NIF
#sum(Sp_latent_state_EM[1:1978,] < 0) / (1978*88) # FMOU



# ### output slip rates
# slip_rate_FMOU = matrix(0, nrow=1978, ncol = 88)
# slip_rate_NIF = matrix(0, nrow=1978, ncol = 88)
# slip_rate_NIF_select_d = matrix(0, nrow=1978, ncol = 88)
# for(tt in 2:88){
#   slip_rate_FMOU[,tt-1] = (Sp_latent_state_EM[1:1978,tt] - Sp_latent_state_EM[1:1978,tt-1]) / (doy[tt]-doy[tt-1])
#   #slip_rate_NIF[,tt-1] = (NIF_slip[1:1978,tt] - NIF_slip[1:1978,tt-1]) / (doy[tt]-doy[tt-1])
#   #slip_rate_NIF_select_d[,tt-1] = (NIF_slip[1:1978,tt] - NIF_slip[1:1978,tt-1]) / (doy[tt]-doy[tt-1])
# }
# 
# slip_rate_FMOU[slip_rate_FMOU<0] = 0
# slip_rate_NIF[slip_rate_NIF<0] = 0
# slip_rate_NIF_select_d[slip_rate_NIF_select_d<0] = 0
# 
# write.csv(Sp_latent_state_EM[1:1978, ], "saved_data/slip_FMOU_large_mesh.csv", row.names=FALSE)
# write.csv(NIF_slip[1:1978, ], "saved_data/slip_NIF_large_mesh.csv", row.names=FALSE)
# write.csv(NIF_slip_select_d[1:1978, ], "saved_data/slip_NIF_large_mesh_select_d.csv", row.names=FALSE)
# write.csv(slip_rate_NIF, "saved_data/slip_rate_NIF_large_mesh.csv", row.names=FALSE)
# write.csv(slip_rate_NIF_select_d, "saved_data/slip_rate_NIF_large_mesh_select_d.csv", row.names=FALSE)
# write.csv(slip_rate_FMOU, "saved_data/slip_rate_FMOU_large_mesh.csv", row.names=FALSE)
