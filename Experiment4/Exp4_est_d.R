# U is known, we select d by IC
library(rstiefel)
library(ggplot2)


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
est_d_IC = matrix(NA, N, 3)

rmse_mean_FMOU_knownU = matrix(NA, N, 3)
rmse_slip_FMOU_knownU = matrix(NA, N, 3)
est_d_IC_knownU = matrix(NA, N, 3)

rmse_mean_FMOU_VM = matrix(NA, N, 3)
rmse_slip_FMOU_VM = matrix(NA, N, 3)
est_d_VM = matrix(NA, N, 3)

L_95_FMOU = matrix(NA, N, 3)
P_95_FMOU = matrix(NA, N, 3)

L_95_slip_FMOU = matrix(NA, N, 3)
P_95_slip_FMOU = matrix(NA, N, 3)

rmse_mean_NIF_full = matrix(NA, N, 3) # NIF is slower
rmse_slip_NIF_full = matrix(NA, N, 3) 

signal_list = as.list(1:3)

for(idx in 1: 3){
	# specify n
  n = n_vec[idx]
  
  for(it in 1:N){
    set.seed(it + 100*(idx-1))
    if(it%%10 == 0){print(it)}
    ## generate latent processes and slips
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
    
    ## generate signal and observations
    signal = G %*% slip  # G G^T U[,1:d_true] %*% z = U Lambda  
    y = signal + matrix(rnorm(n*k,mean=0,sd=sqrt(sigma_0_2)),k,n)
    Psi = t(U) %*% G
    
    signal_list[[idx]] = signal
    
    ## select d by IC (PCA)
    U_pca = svd(y)$u
    loss_score = rep(NA, floor(k*(2/3)))    
    
    for(i_d in 1:length(loss_score)){
      criteria_val_cur = log(mean((y - U_pca[,1:i_d]%*%t(U_pca[,1:i_d])%*%y)^2)) + i_d*(k+n)/(k*n)*log(k*n/(k+n))
      loss_score[i_d] = criteria_val_cur
      
    }
    selected_d = which.min(loss_score)
    est_d_IC[it, idx] = which.min(loss_score)
    
      
    ## Perform FMOU with selected d (PCA)
    
    em_alg <- EM_alg_FMOU(y, d= selected_d, est_U0=F, U0=U[,1: selected_d], est_sigma0_2=T) # fix U, but estimated noise
    Sp_latent_state_EM_fix_noise = t(G) %*% U[,1:selected_d] %*% diag(1/lambda[1: selected_d], nrow= selected_d) %*% em_alg$Z_hat
    U_est_FMOU = as.matrix(U[,1:selected_d], ncol= selected_d)
    fit_latent_state_EM_fix_noise = U_est_FMOU %*% em_alg$Z_hat
    rmse_slip_FMOU[it, idx] = sqrt(mean((Sp_latent_state_EM_fix_noise - slip)^2))
    rmse_mean_FMOU[it, idx] = sqrt(mean((fit_latent_state_EM_fix_noise - signal)^2))
    
    L_95_FMOU[it, idx] = mean(em_alg$pred_mean_95upper-em_alg$pred_mean_95lower)
    P_95_FMOU[it, idx] = mean(em_alg$pred_mean_95upper>signal & em_alg$pred_mean_95lower<signal)
    
    post_var = (em_alg$Z_post_sd)^2
    pred_slip_var = NULL
    for(tt in 1:n){
      pred_slip_var = cbind(pred_slip_var, 
                            diag(t(G) %*% U[,1:selected_d]%*%diag(1/lambda[1:selected_d]) %*% diag(post_var[,tt], nrow=selected_d, ncol=selected_d)%*%diag(1/lambda[1:selected_d])%*%t(U[,1:selected_d])%*%G))
    }
    
    pred_slip_lower = Sp_latent_state_EM_fix_noise - 1.96 * sqrt(pred_slip_var)
    pred_slip_upper = Sp_latent_state_EM_fix_noise + 1.96 * sqrt(pred_slip_var)
    L_95_slip_FMOU[it,idx] = mean(pred_slip_upper-pred_slip_lower)
    P_95_slip_FMOU[it,idx] = mean(pred_slip_upper>slip & pred_slip_lower <slip)
    
    ## select d by IC (known U)
    loss_score_knownU = rep(NA, floor(k*(2/3)))

    for(i_d in 1:length(loss_score)){
      criteria_val_cur = log(mean((y - U[,1:i_d]%*%t(U[,1:i_d])%*%y)^2)) + i_d*(k+n)/(k*n)*log(k*n/(k+n))
      loss_score_knownU[i_d] = criteria_val_cur

    }
    selected_d_knownU = which.min(loss_score_knownU)
    est_d_IC_knownU[it, idx] = which.min(loss_score_knownU)

    ## Perform FMOU with selected d(known U)
    em_alg <- EM_alg_FMOU(y, d= selected_d_knownU, est_U0=F, M=500,
                          U0=U[,1: selected_d_knownU], est_sigma0_2=T) # fix U, but estimated noise
    Sp_latent_state_EM_fix_noise = t(G) %*% U[,1: selected_d_knownU] %*% diag(1/lambda[1: selected_d_knownU], nrow= selected_d_knownU) %*% em_alg$Z_hat
    U_est_FMOU = as.matrix(U[,1: selected_d_knownU], ncol= selected_d_knownU)
    fit_latent_state_EM_fix_noise = U_est_FMOU %*% em_alg$Z_hat
    rmse_slip_FMOU_knownU[it, idx] = sqrt(mean((Sp_latent_state_EM_fix_noise - slip)^2))
    rmse_mean_FMOU_knownU[it, idx] = sqrt(mean((fit_latent_state_EM_fix_noise - signal)^2))

    ## select d by VM
    est_noise_VM = rep(NA, k)
    for(i_d in 1:k){ # k is small so I didn't use binary search here
      em_alg <- EM_alg_FMOU(y, d=i_d, est_U0=F, U0=U[,1: i_d,drop=F],
                            est_sigma0_2=T) # estimate noise
      est_noise_VM[i_d] = em_alg$sigma2_0
    }
    selected_d_VM = which.min(abs(est_noise_VM-sigma_0_2))
    est_d_VM[it, idx] = selected_d_VM

    ## Perform FMOU with selected d by VM
    em_alg <- EM_alg_FMOU(y, d= selected_d_VM, est_U0=F, U0=U[,1: selected_d_VM],
                          est_sigma0_2=F, sigma0_2 = sigma_0_2) # fix U and noise
    Sp_latent_state_EM_fix_noise = t(G) %*% U[,1: selected_d_VM] %*% diag(1/lambda[1: selected_d_VM], nrow= selected_d_VM) %*% em_alg$Z_hat
    U_est_FMOU = as.matrix(U[,1: selected_d_VM], ncol= selected_d_VM)
    fit_latent_state_EM_fix_noise = U_est_FMOU %*% em_alg$Z_hat
    rmse_slip_FMOU_VM[it, idx] = sqrt(mean((Sp_latent_state_EM_fix_noise - slip)^2))
    rmse_mean_FMOU_VM[it, idx] = sqrt(mean((fit_latent_state_EM_fix_noise - signal)^2))

    ## NIF with a full dimension (d=k)
    NIF_full = NIF(y, U, lambda, d=k, Phi=G, t=c(1:n), len_t=n, maxit=50)
    rmse_slip_NIF_full[it, idx] =  sqrt(mean((NIF_full$est_slip - slip)^2))
    rmse_mean_NIF_full[it, idx] = sqrt(mean((NIF_full$est_mean - signal)^2))
    
  }
}



rmse_mean_all_result_n1 = cbind(rmse_mean_FMOU[,1], rmse_mean_FMOU_knownU[,1], 
                                rmse_mean_FMOU_VM[,1],rmse_mean_NIF_full[,1])
rmse_mean_all_result_n2 = cbind(rmse_mean_FMOU[,2], rmse_mean_FMOU_knownU[,2], 
                                rmse_mean_FMOU_VM[,2],rmse_mean_NIF_full[,2])
rmse_mean_all_result_n3 = cbind(rmse_mean_FMOU[,3], rmse_mean_FMOU_knownU[,2], 
                                rmse_mean_FMOU_VM[,3],rmse_mean_NIF_full[,3])
rmse_mean_all_result_e1 = cbind(rmse_mean_all_result_n1, rmse_mean_all_result_n2, 
                                rmse_mean_all_result_n3)



rmse_slip_all_result_n1 = cbind(rmse_slip_FMOU[,1], rmse_slip_FMOU_knownU[,1], 
                                rmse_slip_FMOU_VM[,1], rmse_slip_NIF_full[,1])
rmse_slip_all_result_n2 = cbind(rmse_slip_FMOU[,2], rmse_slip_FMOU_knownU[,2], 
                                rmse_slip_FMOU_VM[,2], rmse_slip_NIF_full[,2])
rmse_slip_all_result_n3 = cbind(rmse_slip_FMOU[,3], rmse_slip_FMOU_knownU[,3], 
                                rmse_slip_FMOU_VM[,3],rmse_slip_NIF_full[,3])
rmse_slip_all_result_e1 = cbind(rmse_slip_all_result_n1, rmse_slip_all_result_n2, rmse_slip_all_result_n3)




### Case 2: k=32, k'=100, d=8, sigma_0^2=1.5, n=100/200/300

N = 20
d_true = 8
k = 32 # number of stations
n_vec = c(100, 200, 300) # number of time points
n_xi= 100 # number of slips
sigma_0_2 = 1.5

sim2_est_d_IC = matrix(NA, N, 3)
sim2_rmse_mean_FMOU = matrix(NA, N, 3) 
sim2_rmse_slip_FMOU = matrix(NA, N, 3) 

sim2_rmse_mean_FMOU_knownU = matrix(NA, N, 3)
sim2_rmse_slip_FMOU_knownU = matrix(NA, N, 3)
sim2_est_d_IC_knownU = matrix(NA, N, 3)

sim2_rmse_mean_FMOU_VM = matrix(NA, N, 3)
sim2_rmse_slip_FMOU_VM = matrix(NA, N, 3)
sim2_est_d_VM = matrix(NA, N, 3)


sim2_rmse_mean_NIF_full = matrix(NA, N, 3)
sim2_rmse_slip_NIF_full = matrix(NA, N, 3)

sim2_L_95_FMOU = matrix(NA, N, 3)
sim2_P_95_FMOU = matrix(NA, N, 3)

sim2_L_95_slip_FMOU = matrix(NA, N, 3)
sim2_P_95_slip_FMOU = matrix(NA, N, 3)

sim2_sigma_list = as.list(1:3)

for(idx in 1: 3){
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
    # signal = U[, 1:d_true] %*% tilde_z
    y = signal + matrix(rnorm(n*k,mean=0,sd=sqrt(sigma_0_2)),k,n)
    Psi = t(U) %*% G
    sim2_sigma_list[[idx]] = signal
    
    
    ## select d by IC
    U_pca = svd(y)$u
    loss_score = rep(NA, floor(k*(2/3)))    
    
    for(i_d in 1:length(loss_score)){
      criteria_val_cur = log(mean((y - U_pca[,1:i_d]%*%t(U_pca[,1:i_d])%*%y)^2)) + i_d*(k+n)/(k*n)*log(k*n/(k+n))
      loss_score[i_d] = criteria_val_cur
      
    }
    selected_d = which.min(loss_score)
    sim2_est_d_IC[it, idx] = which.min(loss_score)
    
      
    ## Perform FMOU with selected d 
    
    em_alg <- EM_alg_FMOU(y, d= selected_d, M=500,
                          est_U0=F, U0=U[,1: selected_d], est_sigma0_2=T) # fix U, but estimated noise
    Sp_latent_state_EM_fix_noise = t(G) %*% U[,1:selected_d] %*% diag(1/lambda[1: selected_d], nrow= selected_d) %*% em_alg$Z_hat
    U_est_FMOU = as.matrix(U[,1:selected_d], ncol= selected_d)
    fit_latent_state_EM_fix_noise = U_est_FMOU %*% em_alg$Z_hat
    sim2_rmse_slip_FMOU[it, idx] = sqrt(mean((Sp_latent_state_EM_fix_noise - slip)^2))
    sim2_rmse_mean_FMOU[it, idx] = sqrt(mean((fit_latent_state_EM_fix_noise - signal)^2))
    
    sim2_L_95_FMOU[it, idx] = mean(em_alg$pred_mean_95upper-em_alg$pred_mean_95lower)
    sim2_P_95_FMOU[it, idx] = mean(em_alg$pred_mean_95upper>signal & em_alg$pred_mean_95lower<signal)
    
    post_var = (em_alg$Z_post_sd)^2
    pred_slip_var = NULL
    for(tt in 1:n){
      pred_slip_var = cbind(pred_slip_var, 
                            diag(t(G) %*% U[,1:selected_d]%*%diag(1/lambda[1:selected_d]) %*% diag(post_var[,tt], nrow=selected_d, ncol=selected_d)%*%diag(1/lambda[1:selected_d])%*%t(U[,1:selected_d])%*%G))
    }
    
    pred_slip_lower = Sp_latent_state_EM_fix_noise - 1.96 * sqrt(pred_slip_var)
    pred_slip_upper = Sp_latent_state_EM_fix_noise + 1.96 * sqrt(pred_slip_var)
    sim2_L_95_slip_FMOU[it,idx] = mean(pred_slip_upper-pred_slip_lower)
    sim2_P_95_slip_FMOU[it,idx] = mean(pred_slip_upper>slip & pred_slip_lower <slip)
    
    # select d by IC (known U)
    loss_score_knownU = rep(NA, floor(k*(2/3)))

    for(i_d in 1:length(loss_score)){
      criteria_val_cur = log(mean((y - U[,1:i_d]%*%t(U[,1:i_d])%*%y)^2)) + i_d*(k+n)/(k*n)*log(k*n/(k+n))
      loss_score_knownU[i_d] = criteria_val_cur

    }
    selected_d_knownU = which.min(loss_score_knownU)
    sim2_est_d_IC_knownU[it, idx] = which.min(loss_score_knownU)

    # Perform FMOU with selected d(known U)
    em_alg <- EM_alg_FMOU(y, d= selected_d_knownU, est_U0=F, U0=U[,1: selected_d_knownU], est_sigma0_2=T) # fix U, but estimated noise
    Sp_latent_state_EM_fix_noise = t(G) %*% U[,1: selected_d_knownU] %*% diag(1/lambda[1: selected_d_knownU], nrow= selected_d_knownU) %*% em_alg$Z_hat
    U_est_FMOU = as.matrix(U[,1: selected_d_knownU], ncol= selected_d_knownU)
    fit_latent_state_EM_fix_noise = U_est_FMOU %*% em_alg$Z_hat
    sim2_rmse_slip_FMOU_knownU[it, idx] = sqrt(mean((Sp_latent_state_EM_fix_noise - slip)^2))
    sim2_rmse_mean_FMOU_knownU[it, idx] = sqrt(mean((fit_latent_state_EM_fix_noise - signal)^2))


    ## select d by VM
    est_noise_VM = rep(NA, k)
    for(i_d in 1:k){ # k is small so I didn't use binary search
      em_alg <- EM_alg_FMOU(y, d= i_d, est_U0=F, U0=U[,1: i_d,drop=F],
                             est_sigma0_2=T) # estimate noise in VM
       est_noise_VM[i_d] = em_alg$sigma2_0
    }
     selected_d_VM = which.min(abs(est_noise_VM-sigma_0_2))
     sim2_est_d_VM[it, idx] = selected_d_VM

    ## Perform FMOU with selected d by VM
    em_alg <- EM_alg_FMOU(y, d= selected_d_VM, est_U0=F, U0=U[,1: selected_d_VM],
                         est_sigma0_2=F, sigma0_2 = sigma_0_2) # fix U and noise
    Sp_latent_state_EM_fix_noise = t(G) %*% U[,1: selected_d_VM] %*% diag(1/lambda[1: selected_d_VM], nrow= selected_d_VM) %*% em_alg$Z_hat
    U_est_FMOU = as.matrix(U[,1: selected_d_VM], ncol= selected_d_VM)
    fit_latent_state_EM_fix_noise = U_est_FMOU %*% em_alg$Z_hat
    sim2_rmse_slip_FMOU_VM[it, idx] = sqrt(mean((Sp_latent_state_EM_fix_noise - slip)^2))
    sim2_rmse_mean_FMOU_VM[it, idx] = sqrt(mean((fit_latent_state_EM_fix_noise - signal)^2))

    ## NIF with a full dimension (d=k)
     NIF_full = NIF(y, U, lambda, d=k, Phi=G, t=c(1:n), len_t=n, maxit=50)
     sim2_rmse_slip_NIF_full[it, idx] =  sqrt(mean((NIF_full$est_slip - slip)^2))
     sim2_rmse_mean_NIF_full[it, idx] = sqrt(mean((NIF_full$est_mean - signal)^2))
    
  }
}



sim2_rmse_mean_all_result_n1 = cbind(sim2_rmse_mean_FMOU[,1], sim2_rmse_mean_FMOU_knownU[,1], 
                                     sim2_rmse_mean_FMOU_VM[,1], sim2_rmse_mean_NIF_full[,1])
sim2_rmse_mean_all_result_n2 = cbind(sim2_rmse_mean_FMOU[,2], sim2_rmse_mean_FMOU_knownU[,2],
                                     sim2_rmse_mean_FMOU_VM[,2], sim2_rmse_mean_NIF_full[,2])
sim2_rmse_mean_all_result_n3 = cbind(sim2_rmse_mean_FMOU[,3], sim2_rmse_mean_FMOU_knownU[,3], 
                                     sim2_rmse_mean_FMOU_VM[,3], sim2_rmse_mean_NIF_full[,3])
sim2_rmse_mean_all_result_e2 = cbind(sim2_rmse_mean_all_result_n1, sim2_rmse_mean_all_result_n2, 		
									 sim2_rmse_mean_all_result_n3)


sim2_rmse_slip_all_result_n1 = cbind(sim2_rmse_slip_FMOU[,1], sim2_rmse_slip_FMOU_knownU[,1],
                                     sim2_rmse_slip_FMOU_VM[,1],sim2_rmse_slip_NIF_full[,1])
sim2_rmse_slip_all_result_n2 = cbind(sim2_rmse_slip_FMOU[,2], sim2_rmse_slip_FMOU_knownU[,2],
                                     sim2_rmse_slip_FMOU_VM[,2],sim2_rmse_slip_NIF_full[,2])
sim2_rmse_slip_all_result_n3 = cbind(sim2_rmse_slip_FMOU[,3], sim2_rmse_slip_FMOU_knownU[,3],
                                     sim2_rmse_slip_FMOU_VM[,3], sim2_rmse_slip_NIF_full[,3])
sim2_rmse_slip_all_result_e2 = cbind(sim2_rmse_slip_all_result_n1, sim2_rmse_slip_all_result_n2, sim2_rmse_slip_all_result_n3)


############################# Table S4 and S5
apply(L_95_FMOU, 2, mean)
apply(sim2_L_95_FMOU, 2, mean)
apply(L_95_slip_FMOU,2,mean)
apply(sim2_L_95_slip_FMOU,2,mean)
apply(P_95_FMOU,2,mean)
apply(sim2_P_95_FMOU, 2, mean)
apply(P_95_slip_FMOU,2,mean)
apply(sim2_P_95_slip_FMOU,2,mean)


## Figure S2: selected d by IC under different combinations of parameters. 

#pdf('../figures/exp3_est_d_IC.pdf', height=4, width=10)
par(mfrow=c(2,3), mar=c(2,5,2,2))
boxplot(est_d_IC-6, names=c("n=100", "n=200", "n=300"), col=c("#37bd79", "#37bd79","#37bd79"),
		ylab=expression(hat(d) - d), #main=expression(k==25 * ", k'=150" * " and " * "d = 6"))
		main="(a) IC", cex.main=1.5, cex.axis = 1.5, cex.lab=1.5)
boxplot(est_d_IC_knownU-6, names=c("n=100", "n=200", "n=300"), col=c("#a7e237", "#a7e237","#a7e237"),
		ylab=expression(hat(d) - d), #main=expression(k==25 * ", k'=150" * " and " * "d = 6"))
		main="(b) IC, known U", cex.main=1.5, cex.axis = 1.5, cex.lab=1.5)
boxplot(est_d_VM-6, names=c("n=100", "n=200", "n=300"), col=c("#f4e604", "#f4e604","#f4e604"),
        ylab=expression(hat(d) - d), #main=expression(k==25 * ", k'=150" * " and " * "d = 6") 
        main="(c) VM", cex.main=1.5, cex.axis = 1.5, cex.lab=1.5)
boxplot(sim2_est_d_IC-8, names=c("n=100", "n=200", "n=300"), col=c("#37bd79", "#37bd79","#37bd79"),
		ylab=expression(hat(d) - d), #main=expression(k==32 * ", k'=100" * " and " * "d = 8"))
		main="(d) IC", cex.main=1.5, cex.axis = 1.5, cex.lab=1.5)
boxplot(sim2_est_d_IC_knownU-8, names=c("n=100", "n=200", "n=300"), col=c("#a7e237", "#a7e237","#a7e237"),
		ylab=expression(hat(d) - d), #main=expression(k==32 * ", k'=100" * " and " * "d = 8"))	
		main="(e) IC, known U", cex.main=1.5, cex.axis = 1.5, cex.lab=1.5)
boxplot(sim2_est_d_VM-8, names=c("n=100", "n=200", "n=300"), col=c("#f4e604", "#f4e604","#f4e604"),
        ylab=expression(hat(d) - d), #main=expression(k==32 * ", k'=100" * " and " * "d = 8"))
        main="(f) VM", cex.main=1.5, cex.axis = 1.5, cex.lab=1.5)
#dev.off()		



## Figure S3: RMSE_m and RMSE_s by d_hat from FMOU, NIF and NIF_full.

#pdf('../figures/exp3_est_d_rmse_mean_slip.pdf',height=8.5, width=10.5)
old_par = par()
# par(mfrow=c(2,2))
par(mfrow=c(2,2),  mar=c(3.9+0.15, 4.8, 5, 1.8), oma=c(2, 2, 2, 2))
par(las=2) 
boxplot(rmse_mean_all_result_e1, las=2, names=c('FMOU','FMOU,U','FMOU,VM','NIF,full', 'FMOU','FMOU, U','FMOU, VM','NIF,full', 'FMOU','FMOU,U','FMOU,VM','NIF_full'),
        main=expression(k==25 * ", k'=150" * " and " * "d = 6"), 
        col=c("#00a6fb",'#B7E0FF','#EE4E4E', "#FFF5CD", "#00a6fb",'#B7E0FF','#EE4E4E', "#FFF5CD", "#00a6fb",'#B7E0FF','#EE4E4E', "#FFF5CD"),
        ylim=c(0.2, 4), at=c(1:4, 6:9, 11:14), ylab=expression(RMSE[m]),cex.main=1.5)
boxplot(sim2_rmse_mean_all_result_e2, las=2, names=c('FMOU','FMOU,U','FMOU,VM','NIF,full', 'FMOU','FMOU, U','FMOU, VM','NIF,full', 'FMOU','FMOU,U','FMOU,VM','NIF_full'),
        main=expression(k==32 * ", k'=100" * " and " * "d = 8"), 
        col=c("#00a6fb",'#B7E0FF','#EE4E4E', "#FFF5CD", "#00a6fb",'#B7E0FF','#EE4E4E', "#FFF5CD", "#00a6fb",'#B7E0FF','#EE4E4E', "#FFF5CD"),
        ylim=c(0.2, 4), at=c(1:4, 6:9, 11:14),  ylab=expression(RMSE[m]),cex.main=1.5)
boxplot(rmse_slip_all_result_e1, las=2, names=c('FMOU','FMOU,U','FMOU,VM','NIF,full', 'FMOU','FMOU, U','FMOU, VM','NIF,full', 'FMOU','FMOU,U','FMOU,VM','NIF_full'),
        main=expression(k==25 * ", k'=150" * " and " * "d = 6"), 
        col=c("#00a6fb",'#B7E0FF','#EE4E4E', "#FFF5CD", "#00a6fb",'#B7E0FF','#EE4E4E', "#FFF5CD", "#00a6fb",'#B7E0FF','#EE4E4E', "#FFF5CD"),
        ylim=c(0, 2), at=c(1:4, 6:9, 11:14), ylab=expression(RMSE[s]),cex.main=1.5)
boxplot(sim2_rmse_slip_all_result_e2, las=2, names=c('FMOU','FMOU,U','FMOU,VM','NIF,full', 'FMOU','FMOU, U','FMOU, VM','NIF,full', 'FMOU','FMOU,U','FMOU,VM','NIF_full'),
        main=expression(k==32 * ", k'=100" * " and " * "d = 8"), 
        col=c("#00a6fb",'#B7E0FF','#EE4E4E', "#FFF5CD", "#00a6fb",'#B7E0FF','#EE4E4E', "#FFF5CD", "#00a6fb",'#B7E0FF','#EE4E4E', "#FFF5CD"),
        ylim=c(0, 2), at=c(1:4, 6:9, 11:14), ylab=expression(RMSE[s]),cex.main=1.5)
par(las=0) 
par(old_par)
par(mfrow=c(1,1))
#dev.off()



