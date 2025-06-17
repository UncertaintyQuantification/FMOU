library(geometry)
library(plot3D)
library(MASS)
library(pracma)


# functions used for NIF and FMOU
source('../functions/functions_geophysics_gppca.R')
source('../functions/FMOU.R')
source('../functions/functions.R')


GF = read.csv("../real_data/Cascadia_data/Green_fun_large_mesh.csv", header=FALSE)
# reads date (epochs) of observations
doy = read.csv("../real_data/Cascadia_data/DOY.csv", header=FALSE)
# obsrvation uncertainties (in meters)
obs_uncertainty = read.csv("../real_data/Cascadia_data/obs_uncertainty.csv", header=FALSE)

GF = GF_ori = as.matrix(GF)
doy = as.numeric(doy)
obs_uncertainty_ori = obs_uncertainty = as.matrix(obs_uncertainty)

epoch = length(doy)
N_site = 351 # number of GPS stations
N_xi = dim(GF)[2] # number of fault patches

########################### The Elliptical Slip example #####################################

# nd_ll: nodes of a triangular mesh. Each row is a node:  col 1 for longitudes, col 2 for latitudes
nd_ll = read.csv("../real_data/Cascadia_data/nd_ll_large_mesh.csv", header=FALSE)
#  indices to the nodes in nd_ll that define each triangle. each row is a triangle with 3 columns that correspond to row indices of nd_ll.
el = read.csv("../real_data/Cascadia_data/el_large_mesh.csv", header=FALSE)
# load latitude and longitude coordinates
lon = read.csv("../real_data/Cascadia_data/lon_G.csv", header=FALSE)
lat = read.csv("../real_data/Cascadia_data/lat_G.csv", header=FALSE)

nd_ll = as.matrix(nd_ll)
el = as.matrix(el)
lon = as.matrix(lon)
lat = as.matrix(lat)

### make the max slip-rate = 3
c= 3
r = 1

# compute the center of each triangle
tricenter = matrix(0, nrow=dim(GF)[2] , 2) # longitude and latitude
for(i in 1:dim(GF)[2]){
  tricenter[i,1]= mean(nd_ll[el[i, ], 1])
  tricenter[i,2]= mean(nd_ll[el[i, ], 2])
}

range(tricenter[,1])
range(tricenter[,2])
range(lon)
range(lat)

semi_minor = 0.7 
# the semi-major will grow linearly, starting from a length semi_major_ini
semi_major_ini = semi_minor/3
# the center of the ellipse
center = c(-123.2, 46.0)

ellipse_len_t = length(doy)
ellipse_t = doy
# the major-axis grows linearly with time, with a grow speed
major_grow_rate = 8


### The grid here and below is just used to show the history (propagation) of slips and slip rates, 
### for simulation we still need to use the coordinate of each patch
x_grid = seq(-126, -121, length.out = 100)
y_grid = seq(43, 49, length.out = 200)

slip_rate_record = as.list(1:ellipse_len_t)
slip_record = as.list(1:ellipse_len_t)

alpha = 2 # used in function grid_slip_update(), grid_slip_rate_update() and point_slip_update() and point_rate_slip_update()

grid_slip_rate_dist_mat = matrix(0, nrow=length(y_grid)*length(x_grid), ncol=ellipse_len_t)
grid_slip_dist_mat = matrix(0, nrow=length(y_grid)*length(x_grid), ncol=ellipse_len_t)
for(i_t in 1:ellipse_len_t){
  grid_slip_dist = grid_slip_update(t0=ellipse_t[1], current_time=ellipse_t[i_t], semi_minor=semi_minor,
                                    semi_major_ini=semi_major_ini, center=center)
  grid_slip_dist_mat[, i_t] = as.vector(grid_slip_dist)
}

for(i_t in 1:(ellipse_len_t-1)){
  grid_slip_rate_dist_mat[,i_t] = (grid_slip_dist_mat[,(i_t+1)] - grid_slip_dist_mat[,i_t])/(ellipse_t[i_t+1]-ellipse_t[i_t]) /365
}


t_seq = c(5, 15, 25, 35, 50, 60, 65, 80)
par(mfrow=c(2,4))
for(i_t in t_seq){
  # grid_slip_dist = grid_slip_update(t0=ellipse_t[1], current_time=ellipse_t[i_t], semi_minor=semi_minor,
  #                                   semi_major_ini=semi_major_ini, center=center)
  grid_slip_dist = matrix(grid_slip_dist_mat[, i_t], nrow=length(y_grid))
  #image2D(t(grid_slip_dist), x=x_grid, y= y_grid, NAcol = "white",
  #        xlab="Longitude", ylab="Latitude", main="slips")
  #points(tricenter[,1], tricenter[,2], cex=0.05, pch=4)
}
par(mfrow=c(1,1))


t_seq = c(5, 15, 25, 35, 50, 60, 65, 80)
par(mfrow=c(2,4))
for(i_t in t_seq){
  # grid_slip_dist = grid_slip_update(t0=ellipse_t[1], current_time=ellipse_t[i_t], semi_minor=semi_minor,
  #                                   semi_major_ini=semi_major_ini, center=center)
  grid_slip_rate_dist = matrix(grid_slip_rate_dist_mat[, i_t], nrow=length(y_grid))
  #image2D(t(grid_slip_rate_dist), x=x_grid, y= y_grid, NAcol = "white",
  #        xlab="Longitude", ylab="Latitude", main="slip rates")
  #points(tricenter[,1], tricenter[,2], cex=0.05, pch=4)
}
par(mfrow=c(1,1))




##### below is the simulation of slips on our triangular fault mesh with 1978 patches (or 691 pathces for rough grids).
# save the history of slips and slip rates in a 2D matrix
point_slip_rate_dist_mat = matrix(0, nrow=dim(el)[1], ncol=ellipse_len_t)
point_slip_dist_mat = matrix(0, nrow=dim(el)[1], ncol=ellipse_len_t)
for(idx in 1:ellipse_len_t){
  point_slip_dist = point_slip_update(t0=ellipse_t[1], current_time=ellipse_t[idx], semi_minor=semi_minor,
                                      semi_major_ini=semi_major_ini, center=center)
  point_slip_dist_mat[,idx] = as.vector(point_slip_dist)
  # point_slip_dist = point_slip_update(t0=ellipse_t[1], current_time=ellipse_t[idx], semi_minor=semi_minor,
  #                                     semi_major_ini=semi_major_ini, center=center)
}

for(idx in 1:(ellipse_len_t-1)){
  # compute the slip rates (per year)
  point_slip_rate_dist_mat[,idx] = (point_slip_dist_mat[,(idx+1)] - point_slip_dist_mat[,idx])/(ellipse_t[idx+1]-ellipse_t[idx])/365
}


missing_location = c(5,  16,  17,  24,  30,  80,  83,  88, 104, 105, 107, 109, 110, 111, 112, 113, 114)
# We only model the measurments in East and North
N_site_used = N_site - length(missing_location)
GF = GF[-c(3*missing_location, 3*missing_location-1, 3*missing_location-2), ]
GF= GF[-seq(3, 3*N_site_used, by=3),]

ellipse_signal = GF %*% point_slip_dist_mat
sigma0_list = c(0.01, 0.05, 0.2)

###### select d ######
eigen_GF=eigen(GF%*%t(GF))
Lambda = eigen_GF$values
U = eigen_GF$vectors
G = U%*%diag(Lambda)
Psi = t(U)%*%GF 

## Select the number of factors, with binary search
seed_list = c(2,0,1)
d_list = rep(0,length(sigma0_list))
ellipse_output_list = as.list(1:length(sigma0_list))
time_select_d = as.list(1:length(sigma0_list))

d_select = NULL

for(idx_noise in 1:length(sigma0_list)){
  sigma0 = sigma0_list[idx_noise]
  set.seed(seed_list[idx_noise])
  ellipse_output = ellipse_signal + matrix(rnorm(length(ellipse_signal), 0, sd=sigma0), nrow=nrow(ellipse_signal), ncol=ncol(ellipse_signal))
  ellipse_output_list[[idx_noise]] = ellipse_output
  
  # use variance matching to select d, binary search
  time = system.time({
    left = 1
    right = dim(U)[2]
    
    while(left <= right){
      mid = left + floor((right-left)/2)
      em_fit = EM_alg_FMOU(ellipse_output, d=mid,
                           est_U0=F, U0=U[,1:mid], est_sigma0_2=T, threshold=10^{-6}) # don't fix noise in this step
      est_noise = em_fit$sigma2_0
      if(est_noise < sigma0^2){
        # overfit, reduce d
        right = mid-1
      }
      else{
        # underfit, increase d
        left = mid + 1
      }
    }
  })
  
  
  d_here = mid
  d_list[idx_noise] = d_here
  time_select_d[idx_noise] = time[3]
  
}

d_list
time_select_d

##########################################################################################################
# now we get numerical results under noise 
FMOU_rmse_mean = rep(NA, length(d_list))
FMOU_rmse_slip = rep(NA, length(d_list))
FMOU_slip_est = as.list(length(d_list))
FMOU_mean_est = as.list(length(d_list))
FMOU_sd_slip = rep(NA, length(d_list))
FMOU_neg_percent = rep(NA, length(d_list))
FMOU_time_record = rep(NA, length(d_list))

L_95_FMOU = rep(NA, 3)
P_95_FMOU = rep(NA, 3)

L_95_slip_FMOU = rep(NA, 3)
P_95_slip_FMOU = rep(NA, 3)

NIF_rmse_mean = rep(NA, length(d_list))
NIF_rmse_slip = rep(NA, length(d_list))
NIF_slip_est = as.list(length(d_list))
NIF_mean_est = as.list(length(d_list))
NIF_sd_slip = rep(NA, length(d_list))
NIF_neg_percent = rep(NA, length(d_list))
NIF_time_record = rep(NA, length(d_list))

NIF_full_rmse_mean = rep(NA, length(d_list))
NIF_full_rmse_slip = rep(NA, length(d_list))
NIF_full_slip_est = as.list(length(d_list))
NIF_full_mean_est = as.list(length(d_list))
NIF_full_sd_slip = rep(NA, length(d_list))
NIF_full_neg_percent = rep(NA, length(d_list))
NIF_full_time_record = rep(NA, length(d_list))



for(nd in 1:length(d_list)){
  print(nd)
  d_here = d_list[nd]
  sigma0 = sigma0_list[nd]
  ellipse_output = ellipse_output_list[[nd]]

  time1 = Sys.time()
  em_alg_here = EM_alg_FMOU(ellipse_output, d=d_here, est_U0=F, est_sigma0_2=F, 
                            U0=U[,1:d_here], sigma0_2=sigma0^2, threshold=10^{-6})
  Sp_latent_state_EM_here = t(Psi[1:d_here,])%*%diag(1/Lambda[1:d_here])%*%em_alg_here$Z_hat
  fit_latent_state_EM_here = (U%*%diag(Lambda))[,1:d_here] %*% diag(1/Lambda[1:d_here])%*%em_alg_here$Z_hat  #estimated data
  time2 = Sys.time()
  FMOU_time_record[nd] = time2 - time1 
  FMOU_rmse_mean[nd] = sqrt(mean((ellipse_signal - fit_latent_state_EM_here)^2, na.rm=T))
  FMOU_rmse_slip[nd] = sqrt(mean((point_slip_dist_mat - Sp_latent_state_EM_here)^2, na.rm=T))
  FMOU_sd_slip[nd] = sd(Sp_latent_state_EM_here)
  FMOU_slip_est[[nd]] = Sp_latent_state_EM_here
  FMOU_mean_est[[nd]] = fit_latent_state_EM_here
  FMOU_neg_percent[nd]  = sum(Sp_latent_state_EM_here < 0) / (ncol(GF)*88)
  
  L_95_FMOU[nd] = mean(em_alg_here$pred_mean_95upper-em_alg_here$pred_mean_95lower)
  P_95_FMOU[nd] = mean(em_alg_here$pred_mean_95upper>ellipse_signal & em_alg_here$pred_mean_95lower<ellipse_signal)
  
  post_var = (em_alg_here$Z_post_sd)^2
  pred_slip_var = NULL
  for(tt in 1:ncol(ellipse_output)){
    pred_slip_var = cbind(pred_slip_var, 
                          diag(t(GF) %*% U[,1:d_here]%*%diag(1/Lambda[1:d_here]) %*% diag(post_var[,tt], nrow=d_here, ncol=d_here)%*%diag(1/Lambda[1:d_here])%*%t(U[,1:d_here])%*%GF))
  }
  pred_slip_lower = Sp_latent_state_EM_here - 1.96 * sqrt(pred_slip_var)
  pred_slip_upper = Sp_latent_state_EM_here + 1.96 * sqrt(pred_slip_var)
  L_95_slip_FMOU[nd] = mean(pred_slip_upper-pred_slip_lower)
  P_95_slip_FMOU[nd] = mean(pred_slip_upper>point_slip_dist_mat & pred_slip_lower <point_slip_dist_mat)
  
  ##### For NIF, Use grid-search instead of optim()
  # NIF - full rank, d=k
  time1 = Sys.time()
  G_NIF = U%*%diag(Lambda)
  Nbasis = nrow(U)
  statedim = 3*Nbasis # dimension of latent states for Kalman Filter
  G_NIF = matrix(G_NIF[,1:Nbasis], ncol=Nbasis)
  theta_1 = c(1000, 1500, 2000, 2500, 3000)
  theta_2 = c(100, 150, 200, 250, 300)
  loglik_theta_combination = matrix(NA, length(theta_1), length(theta_2))
  for(i1 in 1:length(theta_1)){
    for(i2 in 1:length(theta_2)){
      param_here = c(theta_1[i1], theta_2[i2])
      loglik_theta_combination[i1, i2] = neg_loglik_fix_noise(param_here, t=ellipse_t, Data=ellipse_output, G_NIF, Lambda, sigma0^2)
    }
  }
  optim_loglik_index = which.min(loglik_theta_combination)
  theta_2_index = ceiling(optim_loglik_index / length(theta_1))
  theta_1_index = optim_loglik_index - (theta_2_index - 1)*length(theta_1)
  theta_hat = theta_1[theta_1_index]
  theta2_hat = theta_2[theta_2_index]

  tau = 0
  res = netloglik3p(t=ellipse_t, Data=ellipse_output, G_NIF, Lambda, theta_hat, theta2_hat)
  sig2_hat = sigma0^2
  alpha_hat = sqrt(sig2_hat)*theta_hat
  gamma_hat = sqrt(sig2_hat)*theta2_hat

  # construct the initial states
  x0 = matrix(0, statedim, 1)
  var0 = 1.0^(-8)*Lambda[1]*diag(statedim)
  for(j in 1:Nbasis)
  {var0[(3*(j-1)+1),(3*(j-1)+1)] = gamma_hat^2 *Lambda[j]}

  # Run Filter
  res =  netfilt(t=ellipse_t, Data=ellipse_output, G=matrix(G_NIF[,1:Nbasis], ncol=Nbasis), sigmas=c(tau, sqrt(sig2_hat), alpha_hat), x0, var0, 1)
  x_kgk=res[[1]];  sigma_kgk = res[[2]]
  x_kgn = res[[3]]; sigma_kgn= res[[4]]

  c_NIF = matrix(0,Nbasis, ellipse_len_t)  # matrix of stochastic coefficients
  Ts = matrix(0, Nbasis, 3*Nbasis)  # we can think it's "G" in KF

  for(h in 1:ellipse_len_t)
    for(i in 1:Nbasis)
    { Ts[i,(3*(i-1)+1):(3*(i-1)+3)] = c(ellipse_t[h]-ellipse_t[1],1,0)
      c_NIF[,h] = Ts%*%(x_kgn[h,])
    }
  NIF_slip = matrix(t(Psi[1:Nbasis,]), ncol=Nbasis)%*%c_NIF  # estimated slips
  NIF_signal = G_NIF %*% c_NIF   # estimated data
  time2 = Sys.time()

  NIF_full_time_record[nd] = as.numeric(difftime(time2, time1, units="secs"))
  NIF_full_rmse_mean[nd] = sqrt(mean((ellipse_signal - NIF_signal)^2, na.rm=T))
  NIF_full_rmse_slip[nd] = sqrt(mean((point_slip_dist_mat - NIF_slip)^2, na.rm=T))
  NIF_full_sd_slip[nd] = sd(NIF_slip)
  NIF_full_slip_est[[nd]] = NIF_slip
  NIF_full_mean_est[[nd]] = NIF_signal
  NIF_full_neg_percent[nd] = sum(NIF_slip < 0) / (ncol(GF)*88)

  # # NIF - select d
  time1 = Sys.time()
  G_NIF = U%*%diag(Lambda)
  Nbasis = d_here
  statedim = 3*Nbasis # dimension of latent states for Kalman Filter
  G_NIF = matrix(G_NIF[,1:Nbasis], ncol=Nbasis)
  theta_1 = c(1000, 1500, 2000, 2500, 3000)
  theta_2 = c(100, 150, 200, 250, 300)
  loglik_theta_combination = matrix(NA, length(theta_1), length(theta_2))
  for(i1 in 1:length(theta_1)){
    for(i2 in 1:length(theta_2)){
      param_here = c(theta_1[i1], theta_2[i2])
      loglik_theta_combination[i1, i2] = neg_loglik_fix_noise(param_here, t=ellipse_t, Data=ellipse_output, G_NIF, Lambda, sigma0^2)
    }
  }
  optim_loglik_index = which.min(loglik_theta_combination)
  theta_2_index = ceiling(optim_loglik_index / length(theta_1))
  theta_1_index = optim_loglik_index - (theta_2_index - 1)*length(theta_1)
  theta_hat = theta_1[theta_1_index]
  theta2_hat = theta_2[theta_2_index]

  tau = 0
  res = netloglik3p(t=ellipse_t, Data=ellipse_output, G_NIF, Lambda, theta_hat, theta2_hat)
  sig2_hat = sigma0^2
  alpha_hat = sqrt(sig2_hat)*theta_hat
  gamma_hat = sqrt(sig2_hat)*theta2_hat

  # construct the initial states
  x0 = matrix(0, statedim, 1)
  var0 = 1.0^(-8)*Lambda[1]*diag(statedim)
  for(j in 1:Nbasis)
  {var0[(3*(j-1)+1),(3*(j-1)+1)] = gamma_hat^2 *Lambda[j]}

  # Run Filter
  res =  netfilt(t=ellipse_t, Data=ellipse_output, G=matrix(G_NIF[,1:Nbasis], ncol=Nbasis), sigmas=c(tau, sqrt(sig2_hat), alpha_hat), x0, var0, 1)
  x_kgk=res[[1]];  sigma_kgk = res[[2]]
  x_kgn = res[[3]]; sigma_kgn= res[[4]]

  c_NIF = matrix(0,Nbasis, ellipse_len_t)  # matrix of stochastic coefficients
  Ts = matrix(0, Nbasis, 3*Nbasis)  # we can think it's "G" in KF

  for(h in 1:ellipse_len_t)
    for(i in 1:Nbasis)
    { Ts[i,(3*(i-1)+1):(3*(i-1)+3)] = c(ellipse_t[h]-ellipse_t[1],1,0)
    c_NIF[,h] = Ts%*%(x_kgn[h,])
    }
  NIF_slip = matrix(t(Psi[1:Nbasis,]), ncol=Nbasis)%*%c_NIF  # estimated slips
  NIF_signal = G_NIF %*% c_NIF   # estimated data
  time2 = Sys.time()

  NIF_time_record[nd] = as.numeric(difftime(time2, time1, units="secs"))
  NIF_rmse_mean[nd] = sqrt(mean((ellipse_signal - NIF_signal)^2, na.rm=T))
  NIF_rmse_slip[nd] = sqrt(mean((point_slip_dist_mat - NIF_slip)^2, na.rm=T))
  NIF_sd_slip[nd] = sd(NIF_slip)
  NIF_slip_est[[nd]] = NIF_slip
  NIF_mean_est[[nd]] = NIF_signal
  NIF_neg_percent[nd] = sum(NIF_slip < 0) / (ncol(GF)*88)
}

#save.image(file="FMOU_Exp5_saved_results.RData")
##### output data for visualization. These three .csv files are needed for running Figure6.m
# Sp_latent_state_EM = FMOU_slip_est[[3]]
# NIF_slips = NIF_slip_est[[3]]
# write.csv(Sp_latent_state_EM, "saved_data/FMOU_ellipse_slip_large_mesh.csv", row.names=FALSE)
# write.csv(point_slip_dist_mat, "saved_data/true_ellipse_slip_large_mesh.csv", row.names=FALSE)
# write.csv(NIF_slips, "saved_data/NIF_ellipse_slip_large_mesh.csv", row.names=FALSE)

