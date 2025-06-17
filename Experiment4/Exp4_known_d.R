## Exp 4: in all approaches we assume true d is known.

library(rstiefel)
library(Rcpp)
source('../functions/functions_geophysics_gppca.R')
source('../functions/FMOU.R')


### Case 1: k=25, k'=150, d=6, sigma_0^2=1.5, n=100/200/300

N = 20
d_true = 6
k = 25 # number of stations
n_vec = c(100, 200, 300) # number of time points
n_xi= 150 # number of slips
sigma_0_2 = 1.5


rmse_mean_FMOU = matrix(NA, N, 3) 
rmse_slip_FMOU = matrix(NA, N, 3)
rmse_mean_NIF = matrix(NA, N, 3) # NIF is slower
rmse_slip_NIF = matrix(NA, N, 3) 


L_95_FMOU = matrix(NA, N, 3)
P_95_FMOU = matrix(NA, N, 3)

L_95_slip_FMOU = matrix(NA, N, 3)
P_95_slip_FMOU = matrix(NA, N, 3)

signal_list = as.list(1:3)
slip_list = as.list(1:3)

for(idx in 1: 3){
  # specify time steps n
  n = n_vec[idx]
  
  for(it in 1:N){
    set.seed(it + 100*(idx-1))
    if(it%%10 == 0){print(it)}
    ## generate G, slip and latent factors
    U = rustiefel(k, k)
    V = rustiefel(n_xi, k)
    lambda = sort(runif(k, 0, 1), decreasing=TRUE) # assume these are the singular values
    sv = lambda^{1/2}
    G = U %*% diag(sv) %*% t(V)
    tilde_z = z = matrix(NA, d_true, n)
    sigma_2 = runif(d_true, 1, 2)
    rho = runif(d_true, 0.95, 1)
    for(l in 1:d_true){
      R = matrix(NA, n, n)
      for(ir in 1:n){
        for(ic in 1:n){
          R[ir, ic] = rho[l]^(abs(ir-ic))
        }
      }
      R = (sigma_2[l]/(1-rho[l]^2) )* R
      tilde_z[l, ] = t(chol(R)) %*% rnorm(n)
    }
    z = diag(1/lambda[1:d_true]) %*% tilde_z
    slip = t(G) %*% U[, 1:d_true] %*% z
    
    ## generate signal and observation
    signal = G %*% slip  # G G^T U[,1:d_true] %*% z = U Lambda 
    # signal = U[, 1:d_true] %*% tilde_z 
    y = signal + matrix(rnorm(n*k,mean=0,sd=sqrt(sigma_0_2)),k,n)
    Psi = t(U) %*% G
    
    signal_list[[idx]] = signal
    slip_list[[idx]] = slip
    
    
    ### use true d 
    # 1. FMOU, where U and d are known, sigma_0 is unknown
    em_fit <- EM_alg_FMOU(y, d=d_true, 
                          est_U0=F, U0=U[,1:d_true], est_sigma0_2=T)
    Sp_latent_state_EM_fix_noise = t(G) %*% U[,1:d_true] %*% diag(1/lambda[1:d_true], nrow=d_true) %*% em_fit$Z_hat
    U_est_FMOU = as.matrix(U[,1:d_true], ncol=d_true)
    fit_latent_state_EM_fix_noise = U_est_FMOU %*% em_fit$Z_hat
    rmse_slip_FMOU[it, idx] = sqrt(mean((Sp_latent_state_EM_fix_noise - slip)^2))
    rmse_mean_FMOU[it, idx] = sqrt(mean((fit_latent_state_EM_fix_noise - signal)^2))
    
    L_95_FMOU[it, idx] = mean(em_fit$pred_mean_95upper-em_fit$pred_mean_95lower)
    P_95_FMOU[it, idx] = mean(em_fit$pred_mean_95upper>signal & em_fit$pred_mean_95lower<signal)
    
    post_var = (em_fit$Z_post_sd)^2
    pred_slip_var = NULL
    for(tt in 1:n){
      pred_slip_var = cbind(pred_slip_var, 
                            diag(t(G) %*% U[,1:d_true]%*%diag(1/lambda[1:d_true]) %*% diag(post_var[,tt], nrow=d_true, ncol=d_true)%*%diag(1/lambda[1:d_true])%*%t(U[,1:d_true])%*%G))
    }
    
    pred_slip_lower = Sp_latent_state_EM_fix_noise - 1.96 * sqrt(pred_slip_var)
    pred_slip_upper = Sp_latent_state_EM_fix_noise + 1.96 * sqrt(pred_slip_var)
    L_95_slip_FMOU[it,idx] = mean(pred_slip_upper-pred_slip_lower)
    P_95_slip_FMOU[it,idx] = mean(pred_slip_upper>slip & pred_slip_lower <slip)
    
    
    # 2. NIF, true d
    NIF_true_d = NIF(y, U, lambda, d=d_true, Phi=G, t=c(1:n), len_t=n, maxit=50)
    rmse_slip_NIF[it, idx] = sqrt(mean((NIF_true_d$est_slip - slip)^2))
    rmse_mean_NIF[it, idx] = sqrt(mean((NIF_true_d$est_mean - signal)^2))
  }
}


rmse_mean_all_result_n1 = cbind(rmse_mean_FMOU[,1], rmse_mean_NIF[,1])
rmse_mean_all_result_n2 = cbind(rmse_mean_FMOU[,2], rmse_mean_NIF[,2])
rmse_mean_all_result_n3 = cbind(rmse_mean_FMOU[,3], rmse_mean_NIF[,3])
rmse_mean_all_result_e1 = cbind(rmse_mean_all_result_n1, rmse_mean_all_result_n2, rmse_mean_all_result_n3)



rmse_slip_all_result_n1 = cbind(rmse_slip_FMOU[,1], rmse_slip_NIF[,1])
rmse_slip_all_result_n2 = cbind(rmse_slip_FMOU[,2], rmse_slip_NIF[,2])
rmse_slip_all_result_n3 = cbind(rmse_slip_FMOU[,3], rmse_slip_NIF[,3])
rmse_slip_all_result_e1 = cbind(rmse_slip_all_result_n1, rmse_slip_all_result_n2, rmse_slip_all_result_n3)

### Case 2: k=32, k'=100, d=8, sigma_0^2=1.5, n=100/200/300
N = 20
d_true = 8
k = 32 # number of stations
n_vec = c(100, 200, 300)
n_xi= 100 # number of slips
sigma_0_2 = 1.5


sim2_rmse_mean_FMOU = matrix(NA, N, 3) 
sim2_rmse_slip_FMOU = matrix(NA, N, 3) 

sim2_rmse_mean_NIF = matrix(NA, N, 3) 
sim2_rmse_slip_NIF = matrix(NA, N, 3) 

sim2_L_95_FMOU = matrix(NA, N, 3)
sim2_P_95_FMOU = matrix(NA, N, 3)

sim2_L_95_slip_FMOU = matrix(NA, N, 3)
sim2_P_95_slip_FMOU = matrix(NA, N, 3)

sim2_signal_list = as.list(1:3)
sim2_slip_list = as.list(1:3)

for(idx in 1: 3){
  # specify the number of time steps
  n = n_vec[idx]
  
  for(it in 1:N){
    set.seed(it+10*idx)
    if(it%%10 == 0){print(it)}
    ## generate latent processes
    U = rustiefel(k, k)
    V = rustiefel(n_xi, k)
    lambda = sort(runif(k, 0, 1), decreasing=TRUE) # assume these are the singular values
    sv = lambda^{1/2}
    G = U %*% diag(sv) %*% t(V)
    tilde_z = z = matrix(NA, d_true, n)
    sigma_2 = runif(d_true, 1, 2)
    rho = runif(d_true, 0.95, 1)
    for(l in 1:d_true){
      R = matrix(NA, n, n)
      for(ir in 1:n){
        for(ic in 1:n){
          R[ir, ic] = rho[l]^(abs(ir-ic))
        }
      }
      R = (sigma_2[l]/(1-rho[l]^2) )* R
      tilde_z[l, ] = t(chol(R)) %*% rnorm(n)
    }
    z = diag(1/lambda[1:d_true]) %*% tilde_z
    slip = t(G) %*% U[, 1:d_true] %*% z
    signal = G %*% slip
    y = signal + matrix(rnorm(n*k,mean=0,sd=sqrt(sigma_0_2)),k,n)
    Psi = t(U) %*% G
    sim2_signal_list[[idx]] = signal
    sim2_slip_list[[idx]] = slip
    
    # 1. FMOU, where U and d are known, sigma_0 is unknown
    em_fit <- EM_alg_FMOU(y, d=d_true, 
                          est_U0=F, U0=U[,1:d_true], est_sigma0_2=T)
    Sp_latent_state_EM_fix_noise = t(G) %*% U[,1:d_true] %*% diag(1/lambda[1:d_true], nrow=d_true) %*% em_fit$Z_hat
    U_est_FMOU = as.matrix(U[,1:d_true], ncol=d_true)
    fit_latent_state_EM_fix_noise = U_est_FMOU %*% em_fit$Z_hat
    sim2_rmse_slip_FMOU[it, idx] = sqrt(mean((Sp_latent_state_EM_fix_noise - slip)^2))
    sim2_rmse_mean_FMOU[it, idx] = sqrt(mean((fit_latent_state_EM_fix_noise - signal)^2))
    
    sim2_L_95_FMOU[it, idx] = mean(em_fit$pred_mean_95upper-em_fit$pred_mean_95lower)
    sim2_P_95_FMOU[it, idx] = mean(em_fit$pred_mean_95upper>signal & em_fit$pred_mean_95lower<signal)
    
    post_var = (em_fit$Z_post_sd)^2
    pred_slip_var = NULL
    for(tt in 1:n){
      pred_slip_var = cbind(pred_slip_var, 
                            diag(t(G) %*% U[,1:d_true]%*%diag(1/lambda[1:d_true]) %*% diag(post_var[,tt], nrow=d_true, ncol=d_true)%*%diag(1/lambda[1:d_true])%*%t(U[,1:d_true])%*%G))
    }
    
    pred_slip_lower = Sp_latent_state_EM_fix_noise - 1.96 * sqrt(pred_slip_var)
    pred_slip_upper = Sp_latent_state_EM_fix_noise + 1.96 * sqrt(pred_slip_var)
    sim2_L_95_slip_FMOU[it,idx] = mean(pred_slip_upper-pred_slip_lower)
    sim2_P_95_slip_FMOU[it,idx] = mean(pred_slip_upper>slip & pred_slip_lower <slip)
    
    
    # 2. NIF, true d 
    NIF_true_d = NIF(y, U, lambda, d=d_true, Phi=G, t=c(1:n), len_t=n, maxit=50)
    sim2_rmse_slip_NIF[it, idx] = sqrt(mean((NIF_true_d$est_slip - slip)^2))
    sim2_rmse_mean_NIF[it, idx] = sqrt(mean((NIF_true_d$est_mean - signal)^2))
  }
}


sim2_rmse_mean_all_result_n1 = cbind(sim2_rmse_mean_FMOU[,1], sim2_rmse_mean_NIF[,1])
sim2_rmse_mean_all_result_n2 = cbind(sim2_rmse_mean_FMOU[,2], sim2_rmse_mean_NIF[,2])
sim2_rmse_mean_all_result_n3 = cbind(sim2_rmse_mean_FMOU[,3], sim2_rmse_mean_NIF[,3])
sim2_rmse_mean_all_result_e2 = cbind(sim2_rmse_mean_all_result_n1, sim2_rmse_mean_all_result_n2, sim2_rmse_mean_all_result_n3)


sim2_rmse_slip_all_result_n1 = cbind(sim2_rmse_slip_FMOU[,1], sim2_rmse_slip_NIF[,1])
sim2_rmse_slip_all_result_n2 = cbind(sim2_rmse_slip_FMOU[,2], sim2_rmse_slip_NIF[,2])
sim2_rmse_slip_all_result_n3 = cbind(sim2_rmse_slip_FMOU[,3], sim2_rmse_slip_NIF[,3])
sim2_rmse_slip_all_result_e2 = cbind(sim2_rmse_slip_all_result_n1, sim2_rmse_slip_all_result_n2, sim2_rmse_slip_all_result_n3)


################################ Table S4 ----------------------------
# length of 95% CI
apply(L_95_FMOU,2,mean)
apply(sim2_L_95_FMOU,2,mean)

apply(L_95_slip_FMOU,2,mean)
apply(sim2_L_95_slip_FMOU,2,mean)

apply(P_95_FMOU,2,mean)
apply(sim2_P_95_FMOU,2,mean)


#################################### Figure 5 ####################
#pdf('../figures/exp3_true_d_rmse_mean_slip_cut.pdf',height=3.5, width=10)
old_par = par()
# par(mfrow=c(2,2))
par(mfrow=c(1,4),  mar=c(3.9+0.15, 3.8, 3.5, 1.8), oma=c(2, 2, 2, 2),
    mgp = c(1.8, 0.5, 0))
par(las=2) 
boxplot(rmse_mean_all_result_e1, las=2, names=c('FMOU','NIF', 'FMOU','NIF', 'FMOU','NIF'),
        main=expression(k==25 * ", k'=150" * " and " * "d = 6"), 
        col=c("#00a6fb", "#FFF5CD", "#00a6fb", "#FFF5CD", "#00a6fb", "#FFF5CD"),
        ylim=c(0.2,4), at=c(1:2, 4:5, 7:8), ylab=expression(RMSE[m]),
        cex.lab=1.4, cex.axis=1.4,cex.main=1.4)
boxplot(rmse_slip_all_result_e1, las=2, names=c('FMOU','NIF', 'FMOU','NIF', 'FMOU','NIF'),
        main=expression(k==25 * ", k'=150" * " and " * "d = 6"), 
        col=c("#00a6fb", "#FFF5CD", "#00a6fb", "#FFF5CD", "#00a6fb", "#FFF5CD"),
        ylim=c(0, 2), at=c(1:2, 4:5, 7:8), ylab=expression(RMSE[s]),
        cex.lab=1.4, cex.axis=1.4,cex.main=1.4)
boxplot(sim2_rmse_mean_all_result_e2, las=2, names=c('FMOU','NIF', 'FMOU','NIF', 'FMOU','NIF'),
        main=expression(k==32 * ", k'=100" * " and " * "d = 8"), 
        col=c("#00a6fb", "#FFF5CD", "#00a6fb", "#FFF5CD", "#00a6fb", "#FFF5CD"),
        ylim=c(0.2, 4), at=c(1:2, 4:5, 7:8),  ylab=expression(RMSE[m]),
        cex.lab=1.4, cex.axis=1.4,cex.main=1.4)
boxplot(sim2_rmse_slip_all_result_e2, las=2, names=c('FMOU','NIF', 'FMOU','NIF', 'FMOU','NIF'),
        main=expression(k==32 * ", k'=100" * " and " * "d = 8"), 
        col=c("#00a6fb", "#FFF5CD", "#00a6fb", "#FFF5CD", "#00a6fb", "#FFF5CD"),
        ylim=c(0, 1.8), at=c(1:2, 4:5, 7:8), ylab=expression(RMSE[s]),
        cex.lab=1.4, cex.axis=1.4,cex.main=1.4)
par(las=0) 
par(old_par)
par(mfrow=c(1,1))
#dev.off()


