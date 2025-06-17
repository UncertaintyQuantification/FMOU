neg_log_lik_shared_cov_FFBS<-function(param,kernel_type){
  
  G_log_det_cov=Get_G_log_det_cov(param, output, delta_x,d=1,kernel_type=kernel_type)
  G=output_2-G_log_det_cov[[1]]
  
  eigen_G=eigen(G)
  n = ncol(output)
  k = nrow(output)
  -(-(sum(G_log_det_cov[[2]]))*d/2-(n*k)/2*log(trace_output_2-sum(eigen_G$values[1:d]) ))
  
}

neg_log_lik_shared_cov_FFBS_fixed_A<-function(param,fixed_A,kernel_type){
  
  G_log_det_cov=Get_G_log_det_cov(param, output, delta_x,d=1,kernel_type=kernel_type)
  G=output_2-G_log_det_cov[[1]]
  
  n = ncol(output)
  k = nrow(output)
  -(-(sum(G_log_det_cov[[2]]))*d/2-(n*k)/2*log(trace_output_2-sum( (fixed_A%*%t(fixed_A))*G)))
  
}

neg_log_lik_diff_cov_FFBS<-function(param,A_ini,kernel_type='matern_5_2'){
  
  G_log_det_cov=Get_G_log_det_cov(param,output, delta_x,d,kernel_type)
  # G_log_det_cov=Get_G_log_det_cov(c(rep(param[1],d),param[-1]), output_sub_mean, delta_x,d)
  
  G=as.list(1:d)
  
  for(i in 1:d){
    G[[i]]=output_2-G_log_det_cov[[i]]
  }
  n = ncol(output)
  k = nrow(output)
  ##my method
  A_hat_here=Optimization_Stiefel_Manifold(A_ini, G=G,max_iter=200)
  
  S_2=trace_output_2+F_Funct(A_hat_here,G)
  neg_log_lik=-(-(sum(G_log_det_cov[[d+1]]))/2-(n*k)/2*log(S_2))
  # print(param)
  # print(S_2/(n*k))
  # print(neg_log_lik)
  return(neg_log_lik)
  
  #-(-(sum(G_log_det_cov[[d+1]]))/2-(n*k)/2*log(trace_output_sub_mean_2+F_Funct(A_hat_here,G)))
}

neg_log_lik_diff_cov_FFBS_fixed_A<-function(param,fixed_A,kernel_type='matern_5_2'){
  #print(param)
  
  
  G_log_det_cov=Get_G_log_det_cov(param, output, delta_x,d,kernel_type)
  G=as.list(1:d)
  n = ncol(output)
  k = nrow(output)
  
  for(i in 1:d){
    G[[i]]=output_2-G_log_det_cov[[i]]
  }
  
  -(-(sum(G_log_det_cov[[d+1]]))/2-(n*k)/2*log(trace_output_2+F_Funct(fixed_A,G)))
  
}

pred_FFBS_FastGaSP<-function(param,A_hat,input,testing_input, output_here,d,var_data=F,kernel_type){
  # beta=param[1]
  # sigma_2=param[2]
  # sigma_2_0=param[3]
  
  beta=param[1:d]
  sigma_2=param[(d+1):(2*d)]
  sigma_2_0=param[2*d+1]
  num_testing=length(testing_input)
  num_obs=length(input)
  
  predict_all_dlm=matrix(0, num_obs,d)
  var_all_dlm=matrix(0, num_obs,d)
  output_t_A=t(output_here)%*%A_hat
  var_est=0
  
  #if(needcov==T){
  #  var_est=array(0,c(k,k,num_testing))
  #}
  
  for(i_d in 1:d){
    
    m.model=fgasp(input,output_t_A[,i_d],have_noise=TRUE,kernel_type=kernel_type)
    m.pred=predict(param=c(log(beta[i_d]),log(sigma_2_0/sigma_2[i_d])),object=m.model,
                   testing_input=testing_input,var_data=FALSE,sigma_2=sigma_2[i_d])
    predict_all_dlm[,i_d]=m.pred@mean
    
    var_all_dlm[,i_d]= m.pred@var  
    
    var_est= var_est+(A_hat[,i_d]^2)%*%t(as.matrix(var_all_dlm[,i_d],num_testing,1))
  }
  
  return.list=as.list(1:3)
  return.list[[1]]=A_hat%*%t(predict_all_dlm)
  ###here I include the noise in the data
  if(var_data==F){
    return.list[[2]]=var_est
    #t(var_all_dlm)
  }else{
    return.list[[2]]=var_est+rep(sigma_2_0,k)
  }
  
  return.list[[3]]=t(predict_all_dlm)
  
  return.list
  
}

get_chol<-function(x,beta){
  R0_00=abs(outer(x,x,'-'))
  R=matern_5_2_funct(R0_00,beta)
  rcppeigen_get_chol(R)
}

#############################################################################################
################################# slips and Green's function ##############################################

slipxt <- function(x, t_current, t0, tf, a=NA, L, slip_amp){
  # x = position on fault surface (\xi),  t_current = current time
  #	t0 = initial time,  tf = final time
  #	L = length of fault,  slip_amp = maximum slip amplitude
  if(is.na(a))
  { a = 0.25*L}
  else{a = a}

  # relt = (t-t0)/(tf-t0)
  temp = a + (L-2*a)*t_current/tf
  if(a^2 - (x - temp)^2 >=0 ) temp = (a^2 - (x - temp)^2)^(1.5) 
  else temp = 0
  return(temp*slip_amp/(a^3))
}

get_Phi <- function(D, x_loc, x){
  phi = D/(D^2 + (x_loc-x)^2)/pi
  #phi = D*atan(abs(x_loc))*sin(x)
  return(phi)
}

################################ Paul's NIF #################################################
makeH <- function(Nbasis, tk, G){
  # %%  Input:
  # %%	Nsites   = number of observation sites
  # %%	Nbasis	 = number of basis functions
  # %%	tk       = t_k time at present iteration
  # %%	G        = Matrix that maps slip to data
  sub = matrix(0, Nbasis, 3*Nbasis)
  for(i in 1:Nbasis)
  {sub[i,(3*(i-1)+1):(3*(i-1)+3)] =c(tk,1,0)}
  
  H = G%*%sub
  return(H)
}  # H maps state vectors to observations 

makeF <- function(statedim, Nbasis, delt){  # this is matrix "G" that update state vectors in KF
  F_mat = diag(statedim)
  
  for(i in 1:Nbasis){
    index = 3*(i-1)+2
    F_mat[index,index+1] = delt
  }
  return(F_mat)
} # thisi is matrix "G" that update state vectors in KF

makeQ <- function(statedim, Nbasis, delt, alpha, tau){
  Q = matrix(0, 3*Nbasis, 3*Nbasis)
  for(i in 1:Nbasis){
    index = 3*(i-1)+2
    Q[index,index] = alpha^2*delt^3/3
    Q[index,index+1] = alpha^2*delt^2/2
    Q[index+1,index] = Q[index,index+1]	
    Q[index+1,index+1] = alpha^2*delt		
  }
  return(Q)
} # Like the matrix "W" in KF

netloglik3p <- function(t, Data, G, Lambda, theta2, theta3){
  # %% Input :  
  # %%	t	= vector of observation times 
  # %%	Data 	= Nsites*Nepochs matrix of observations
  # %%	G 	= Nsites*Nbasis maps slip to Signal
  # %% 	Nsites	= number of observation sites
  # %%	Lambda	 = eigenvalues of the Gram Matrix (in vector form)
  # %%		theta2 = alpha/sigma;
  # %%		theta3 = gamma/sigma;
  # %% 		tau	= random walk standard deviation, mm/sqrt(yr)
  # %% 		sigma	= white noise standard deviation, mm
  # %% 		alpha 	= scale of Weiner process	
  # %% 		gamma 	= smoothness prior standard devation	
  #   %% Output:      val  = C - 2 * log(L)
  #   %% 		C  (= n - n log n) is a constant 
  #   %% 		L is the profile likelihood  
  #   %%		sig2 = estimated variance
  
  #     %%  State Vector: x = [v_1, W_1(t), dot{W_1(t)}, v_2, W_2(t)....]'
  
  #  Observation eq:  d_k = H_k*x_k + eps_k		eps_k ~ N(0,sig^2*R_k)
  #  Update equation:  x_k+1 = F_k*x_k + delta_k		delta_k ~ N(0,sig^2*Q_k)
  
  # Determine some dimensions
  Nsites = dim(Data)[1]
  Nepochs = dim(Data)[2]
  n = Nepochs
  ncheck = dim(G)[1]
  Nbasis  = dim(G)[2]
  
  statedim = 3*Nbasis
  t = t - min(t)
  
  tau = 0
  alpha=theta2
  vsum = 0
  rsum = 0
  
  # Matricies to be used in smoothing operation & for Square Root Filter
  s = matrix(0, statedim+Nsites,  statedim+Nsites)
  # Time invariant operators:
  R = diag(Nsites)
  
  # Starting (prior) values for state and covariance matrix:
  x = matrix(0, statedim,ncol=1)
  Small = 10^{-8}*Lambda[1]
  covx = Small*diag(statedim)
  for(j in 1:Nbasis)
  {covx[(3*(j-1)+1),(3*(j-1)+1)] = theta3^2*Lambda[j]}
  
  # First update step to get x_1/1 and covariance, using straight Kalman filter
  H = makeH(Nbasis, t[1], G)
  nu = Data[,1] - H%*%x           # "Inovation" or prediction error
  invH = solve(R + H%*%covx%*%t(H))
  
  g = covx%*%t(H)%*%invH    # Kalman gain
  x = x + g%*%nu
  covx = covx - g%*%H%*%covx
  
  # Components of Likelihood
  rsum = rsum + t(nu)%*%invH%*%nu
  vsum = vsum + log(  det(R + H%*%covx%*%t(H))  )
  
  # Run Filter
  for(k in 2:n){
    #  Prediction Step:
    delt = t[k] - t[k-1]
    F_mat = makeF(statedim, Nbasis, delt) # State Transition Matrix
    Q = makeQ(statedim, Nbasis, delt, alpha, tau)
    
    x = F_mat%*%x
    covx =  F_mat%*%covx%*%t(F_mat) + Q
    sqrsigma = t(chol(covx))
    
    # Update Step:
    H = makeH(Nbasis, t[k], G)
    nu = Data[,k] - H%*%x        # "Inovation" or prediction error
    u = H%*%sqrsigma
    invH = solve(R + u%*%t(u))
    g = covx%*%t(H)%*%invH		# Kalman gain
    x = x + g%*%nu
    
    s[1:Nsites, 1:Nsites] = t(chol(R))
    s[1:Nsites, ((Nsites+1):(statedim+Nsites))] = matrix(0, Nsites,statedim)
    s[(Nsites+1):(statedim+Nsites),1:Nsites] = t(u)
    s[(Nsites+1):(statedim+Nsites),((Nsites+1):(statedim+Nsites))] = t(sqrsigma)
    # s = [ t(chol(R)), matrix(0, Nsites,statedim); t(u), t(sqrsigma)]
    qrres = qr(s)
    q1 = qr.Q(qrres)
    r1 = qr.R(qrres)
    sqrsigma = t(r1[(Nsites+1):(Nsites+statedim), (Nsites+1):(Nsites+statedim)])
    covx = sqrsigma%*%t(sqrsigma)
    
    # Components of Likelihood
    rsum = rsum + t(nu)%*%invH%*%nu;
    vsum = vsum + log( det(R + u%*%t(u)) )
  }
  
  N = Nsites*Nepochs
  val = N*log(rsum) + vsum
  sig2 = rsum/N
  
  return(c(val, sig2))
} # run KF to compute profile likelihood

netfilt <- function(t, Data, G, sigmas, x0, var0, smooth){
  #	sigmas = [tau, sigma, alpha]
  #	x0 	= state vector at initial epoch (length statedim)
  #	var0 	= diagonal of state covariance at initial epoch 
  
  # %% Output : 
  # %%		x_kp1gk 	= predicted state vector, x_{k+1|k}
  # %%		sigma_kp1gk 	= corresponding covariance
  # %%		x_kgk		= filtered state vector, x_{k|k}
  # %%		sigma_kgk	= corresponding covariance
  # %%		x_kgn		= filtered state vector, x_{k|N}
  # %%		sigma_kgn	= corresponding covariance
  # %%
  
  # Determine some dimensions
  Nsites = dim(Data)[1]
  Nepochs = dim(Data)[2]
  n = Nepochs
  ncheck = dim(G)[1]
  Nbasis  = dim(G)[2]
  
  statedim = 3*Nbasis
  
  # Variance parameters
  tau = sigmas[1]/sigmas[2] 
  alpha=sigmas[3]/sigmas[2]
  
  # Data Covariance
  R = diag(Nsites)
  
  # Matricies to be used in smoothing operation & for Square Root Filter
  s = matrix(0, statedim+Nsites, statedim+Nsites)
  x_kp1gk = matrix(0, n,statedim); sigma_kp1gk =  matrix(0, n,statedim^2)
  x_kgk = matrix(0,n,statedim); sigma_kgk =  matrix(0,n,statedim^2)
  x_kgn =  matrix(0,n,statedim); sigma_kgn =  matrix(0,n,statedim^2)
  
  # Starting (prior) values for state and covariance matrix:
  x = x0
  covx = var0
  x_kp1gk[1,] = t(x)            	 # x_{k+1|k} for smoothing
  sigma_kp1gk[1,] = as.vector(covx)		# Sigma_{k+1|k} for smoothing
  
  
  # First update step to get x_1/1 and covariance, using straight Kalman filter
  H = makeH(Nbasis, t[1], G)
  nu = Data[,1] - H%*%x           # "Inovation" or prediction error
  g = covx%*%t(H)%*%solve(R + H%*%covx%*%t(H))    # Kalman gain
  x = x + g%*%nu
  covx = covx - g%*%H%*%covx
  
  x_kgk[1,] = t(x)                # k_{1|1} for smoothing operation
  sigma_kgk[1,] = as.vector(covx)     #Sigma_{1|1} for smoothing operation
  
  # Run Square Root Filter
  for(k in 2:n){
    # Prediction Step:
    delt = t[k] - t[k-1]
    F_mat = makeF(statedim, Nbasis, delt) # State Transition Matrix
    Q = makeQ(statedim, Nbasis, delt, alpha, tau)
    
    x = F_mat%*%x
    if(kappa(Q)>1e+8){
      covx =  F_mat%*%covx%*%t(F_mat) + Q + 10^{-8}*diag(statedim)
    }else{
      covx =  F_mat%*%covx%*%t(F_mat) + Q 
    }
    sqrsigma = t(chol(covx))
    
    x_kp1gk[k,] = t(x)		# k_{k+1|k} for smoothing 
    sigma_kp1gk[k,] = as.vector(covx)    # Sigma_{k+1|k} for smoothing
    
    # Update Step:
    H = makeH(Nbasis, t[k], G)
    nu = Data[,k] - H%*%x        # "Inovation" or prediction error
    u = H%*%sqrsigma
    invH = solve(R + u%*%t(u))
    g = covx%*%t(H)%*%invH		# Kalman gain
    x = x + g%*%nu
    
    s[1:Nsites, 1:Nsites] = t(chol(R))
    s[1:Nsites, ((Nsites+1):(statedim+Nsites))] = matrix(0, Nsites,statedim)
    s[(Nsites+1):(statedim+Nsites),1:Nsites] = t(u)
    s[(Nsites+1):(statedim+Nsites),((Nsites+1):(statedim+Nsites))] = t(sqrsigma)
    # s = [ t(chol(R)), matrix(0, Nsites,statedim); t(u), t(sqrsigma)]
    qrres = qr(s)
    q1 = qr.Q(qrres)
    r1 = qr.R(qrres)
    sqrsigma = t(r1[(Nsites+1):(Nsites+statedim), (Nsites+1):(Nsites+statedim)])
    covx = sqrsigma%*%t(sqrsigma)
    
    x_kgk[k,] = t(x)		# k_{k|k} for smoothing 
    sigma_kgk[k,] = as.vector(covx)	#Sigma_{k|k} for smoothing
    
  }
  
  if(smooth ==1){
    
    #Run Smoother 
    x_kgn[n,] = t(x)
    sigma_kgn[n,] = sigma_kgk[n,]
    
    s = matrix(0, statedim,statedim)
    sigkgk = diag(statedim)
    sigkp1gN = diag(statedim)
    
    for(k in (n-1):1){
      delt = t[k+1] - t[k]
      F_mat = makeF(statedim, Nbasis, delt) # State Transition Matrix
      Q = makeQ(statedim, Nbasis, delt, alpha, tau) # "Process" Covariance
      sigkgk = matrix(sigma_kgk[k,], statedim, statedim)
      sigkp1gk =  F_mat%*%sigkgk%*%t(F_mat) + Q
      s = sigkgk%*%t(F_mat)%*%solve(sigkp1gk + 10^{-8}*diag(statedim))
      
      x_k1gk = F_mat%*%(x_kgk[k,])
      x = x_kgk[k,] + s%*%(x_kgn[k+1,] -  x_k1gk)
      x_kgn[k,] = as.vector(x)
      
      sigkp1gN = matrix(sigma_kgn[k+1,],statedim,statedim)
      
      sigkgN = sigkgk + s%*%(sigkp1gN -  sigkp1gk)%*%t(s)
      sigma_kgn[k,] = as.vector(sigkgN)
    }
    
  }
  res = as.list(1:6)
  res[[1]] =x_kgk;   res[[2]] = sigma_kgk
  res[[3]] = x_kgn;   res[[4]] = sigma_kgn;
  res[[5]] = x_kp1gk;    res[[6]] = sigma_kp1gk;
  return(res)
} # the NIF model, given all parameters

neg_loglik <- function(param, t, Data, G, Lambda){
  # the matrix H where the linear trend of stochastic coeffcient c is dropped
  # %% Input :  
  # %%	t	= vector of observation times 
  # %%	Data 	= Nsites*Nepochs matrix of observations
  # %%	G 	= Nsites*Nbasis maps slip to Signal
  # %% 			Nsites	= number of observation sites
  # %% 			Nepochs	= number of observation epochs
  # %%	Lambda	 = eigenvalues of the Gram Matrix (in vector form)
  # %%		theta2 = alpha/sigma;
  # %%		theta3 = gamma/sigma;
  # %% 		tau	= random walk standard deviation, mm/sqrt(yr)
  # %% 		sigma	= white noise standard deviation, mm
  # %% 		alpha 	= scale of Weiner process	
  # %% 		gamma 	= smoothness prior standard devation	
  #   %% Output:      val  = C - 2 * log(L)
  #   %% 		C  (= n - n log n) is a constant 
  #   %% 		L is the profile likelihood  
  #   %%		sig2 = estimated variance
  
  #     %%  State Vector: x = [b_1, W_1(t), dot{W_1(t)}, b_2, W_2(t)....]'
  
  #  Observation eq:  d_k = H_k*x_k + eps_k		eps_k ~ N(0,sig^2*R_k)
  #  Update equation:  x_k+1 = F_k*x_k + delta_k		delta_k ~ N(0,sig^2*Q_k)
  theta2 = param[1]
  theta3=param[2]
  # Determine some dimensions
  Nsites = dim(Data)[1]
  Nepochs = dim(Data)[2]
  n = Nepochs
  ncheck = dim(G)[1]
  Nbasis  = dim(G)[2]
  
  statedim = 3*Nbasis
  t = t - min(t)
  
  Big = 10000
  Small = 1.0e-8
  #	tau = sigmas(1); sigma = sigmas(2); 
  tau = 0
  alpha=theta2
  vsum = 0
  rsum = 0
  
  # Matricies to be used in smoothing operation & for Square Root Filter
  s = matrix(0, statedim+Nsites,  statedim+Nsites)
  # Time invariant operators:
  R = diag(Nsites)
  
  # Starting (prior) values for state and covariance matrix:
  x = matrix(0, statedim,ncol=1)
  Small = 10^{-8}*Lambda[1]
  covx = Small*diag(statedim)
  for(j in 1:Nbasis)
  {covx[(3*(j-1)+1),(3*(j-1)+1)] = theta3^2 /Lambda[j]}
  
  
  # First update step to get x_1/1 and covariance, using straight Kalman filter
  H = makeH(Nbasis, t[1], G)
  nu = Data[,1] - H%*%x           # "Inovation" or prediction error
  invH = solve(R + H%*%covx%*%t(H))
  
  g = covx%*%t(H)%*%invH    # Kalman gain
  x = x + g%*%nu
  covx = covx - g%*%H%*%covx
  
  # Components of Likelihood
  rsum = rsum + t(nu)%*%invH%*%nu
  vsum = vsum + log(  det(R + H%*%covx%*%t(H))  )
  
  # Run Square Root Filter
  for(k in 2:n){
    #  Prediction Step:
    delt = t[k] - t[k-1]
    F_mat = makeF(statedim, Nbasis, delt) # State Transition Matrix
    Q = makeQ(statedim, Nbasis, delt, alpha, tau)
    
    x = F_mat%*%x
    covx =  F_mat%*%covx%*%t(F_mat) + Q
    sqrsigma = t(chol(covx))
    
    # Update Step:
    H = makeH(Nbasis, t[k], G)
    nu = Data[,k] - H%*%x        # "Inovation" or prediction error
    u = H%*%sqrsigma
    invH = solve(R + u%*%t(u))
    g = covx%*%t(H)%*%invH		# Kalman gain
    x = x + g%*%nu
    
    s[1:Nsites, 1:Nsites] = t(chol(R))
    s[1:Nsites, ((Nsites+1):(statedim+Nsites))] = matrix(0, Nsites,statedim)
    s[(Nsites+1):(statedim+Nsites),1:Nsites] = t(u)
    s[(Nsites+1):(statedim+Nsites),((Nsites+1):(statedim+Nsites))] = t(sqrsigma)
    # s = [ t(chol(R)), matrix(0, Nsites,statedim); t(u), t(sqrsigma)]
    qrres = qr(s)
    q1 = qr.Q(qrres)
    r1 = qr.R(qrres)
    sqrsigma = t(r1[(Nsites+1):(Nsites+statedim), (Nsites+1):(Nsites+statedim)])
    covx = sqrsigma%*%t(sqrsigma)
    
    # Components of Likelihood
    rsum = rsum + t(nu)%*%invH%*%nu;
    vsum = vsum + log( det(R + u%*%t(u)) )
  }
  
  N = Nsites*Nepochs
  val = N*log(rsum) + vsum
  #sig2 = rsum/N
  return(as.numeric(val))
}


neg_loglik_fix_noise <- function(param, t, Data, G, Lambda, sigma_0_2){
  # the matrix H where the linear trend of stochastic coeffcient c is dropped
  # %% Input :  
  # %%	t	= vector of observation times 
  # %%	Data 	= Nsites*Nepochs matrix of observations
  # %%	G 	= Nsites*Nbasis maps slip to Signal
  # %% 			Nsites	= number of observation sites
  # %% 			Nepochs	= number of observation epochs
  # %%	Lambda	 = eigenvalues of the Gram Matrix (in vector form)
  # %%		theta2 = alpha/sigma_0_2;
  # %%		theta3 = gamma/sigma_0_2;
  # %% 		tau	= random walk standard deviation
  # %% 		sigma	= white noise standard deviation
  # %% 		alpha 	= scale of Weiner process	
  # %% 		gamma 	= smoothness prior standard deviation	
  #   %% Output:      val  = C - 2 * log(L)
  #   %% 		C  (= n - n log n) is a constant 
  #   %% 		L is the profile likelihood  
  #   %%		sig2 = estimated variance
  
  #     %%  State Vector: x = [b_1, W_1(t), dot{W_1(t)}, b_2, W_2(t)....]'
  
  #  Observation eq:  d_k = H_k*x_k + eps_k		eps_k ~ N(0,sig^2*R_k)
  #  Update equation:  x_k+1 = F_k*x_k + delta_k		delta_k ~ N(0,sig^2*Q_k)
  theta2 = param[1]
  theta3=param[2]
  # Determine some dimensions
  Nsites = dim(Data)[1]
  Nepochs = dim(Data)[2]
  n = Nepochs
  ncheck = dim(G)[1]
  Nbasis  = dim(G)[2]
  
  statedim = 3*Nbasis
  t = t - min(t)
  
  Big = 10000
  Small = 1.0e-8
  #	tau = sigmas(1); sigma = sigmas(2); 
  tau = 0
  alpha=theta2
  vsum = 0
  rsum = 0
  
  # Matricies to be used in smoothing operation & for Square Root Filter
  s = matrix(0, statedim+Nsites,  statedim+Nsites)
  # Time invariant operators:
  R = diag(Nsites)
  
  # Starting (prior) values for state and covariance matrix:
  x = matrix(0, statedim,ncol=1)
  Small = 10^{-8}*Lambda[1]
  covx = Small*diag(statedim)
  for(j in 1:Nbasis)
  {covx[(3*(j-1)+1),(3*(j-1)+1)] = theta3^2 /Lambda[j]}
  
  # First update step to get x_1/1 and covariance, using straight Kalman filter
  H = makeH(Nbasis, t[1], G)
  nu = Data[,1] - H%*%x           # "Inovation" or prediction error
  invH = solve(R + H%*%covx%*%t(H))
  
  g = covx%*%t(H)%*%invH    # Kalman gain
  x = x + g%*%nu
  covx = covx - g%*%H%*%covx
  
  # Components of Likelihood
  rsum = rsum + t(nu)%*%invH%*%nu
  vsum = vsum + log(  det(R + H%*%covx%*%t(H))  )
  
  # Run Square Root Filter
  for(k in 2:n){
    #  Prediction Step:
    delt = t[k] - t[k-1]
    F_mat = makeF(statedim, Nbasis, delt) # State Transition Matrix
    Q = makeQ(statedim, Nbasis, delt, alpha, tau)
    
    x = F_mat%*%x
    covx =  F_mat%*%covx%*%t(F_mat) + Q
    sqrsigma = t(chol(covx))
    
    # Update Step:
    H = makeH(Nbasis, t[k], G)
    nu = Data[,k] - H%*%x        # "Inovation" or prediction error
    u = H%*%sqrsigma
    invH = solve(R + u%*%t(u))
    g = covx%*%t(H)%*%invH		# Kalman gain
    x = x + g%*%nu
    
    s[1:Nsites, 1:Nsites] = t(chol(R))
    s[1:Nsites, ((Nsites+1):(statedim+Nsites))] = matrix(0, Nsites,statedim)
    s[(Nsites+1):(statedim+Nsites),1:Nsites] = t(u)
    s[(Nsites+1):(statedim+Nsites),((Nsites+1):(statedim+Nsites))] = t(sqrsigma)
    # s = [ t(chol(R)), matrix(0, Nsites,statedim); t(u), t(sqrsigma)]
    qrres = qr(s)
    q1 = qr.Q(qrres)
    r1 = qr.R(qrres)
    sqrsigma = t(r1[(Nsites+1):(Nsites+statedim), (Nsites+1):(Nsites+statedim)])
    covx = sqrsigma%*%t(sqrsigma)
    
    # Components of Likelihood
    rsum = rsum + t(nu)%*%invH%*%nu;
    vsum = vsum + log( det(R + u%*%t(u)) )
  }
  
  N = Nsites*Nepochs
  val = rsum /sigma_0_2 + vsum + N*log(sigma_0_2)
  #sig2 = rsum/N
  return(as.numeric(val))
}


NIF_fix_noise <- function(y, U_NIF, lambda, d, Phi, t, len_t, noise_level, maxit=100){
  time1 = Sys.time()
  tau = 0 # set it to be 0 if we don't want to use a benchmark motion 
  G_NIF = U_NIF%*%diag(lambda)
  Nbasis = d
  G_NIF = matrix(G_NIF[,1:Nbasis], ncol=Nbasis)
  statedim = 3*Nbasis # dimension of latent states for Kalman Filter
  Psi = t(U_NIF)%*%Phi    # compute basis functions
  
  param_ini= c(1, 10)
  m=try(optim(param_ini, neg_loglik_fix_noise, t=t, Data=y, G=G_NIF, Lambda=lambda, sigma_0_2 = noise_level,
              lower = c(1e-6,1e-6), upper=c(1e+5, 1e+5),
              method="L-BFGS-B", control = list(maxit = maxit)),silent=T)
  # while(is.character(m)){
  #   param_ini=param_ini+runif(2)
  #   m=try(optim(param_ini,neg_loglik, t=ellipse_t, Data=ellipse_output, G=G_NIF, Lambda=lambda,lower = c(1e-6,1e-6),upper=c(1e+4,1e+4),
  #               method="L-BFGS-B", control = list(maxit = 100)),silent=T)
  # }
  theta_hat = m$par[1]
  theta2_hat = m$par[2]
  tau = 0
  
  res = netloglik3p(t, y, G_NIF, lambda, theta_hat, theta2_hat)
  # sig2_hat = res[2] 
  sig2_hat = noise_level
  alpha_hat = sqrt(sig2_hat)*theta_hat
  gamma_hat = sqrt(sig2_hat)*theta2_hat
  
  # construct the initial states
  x0 = matrix(0, statedim, 1)
  var0 = 1.0^(-10)*lambda[1]*diag(statedim)
  for(j in 1:Nbasis)
  {var0[(3*(j-1)+1),(3*(j-1)+1)] = gamma_hat^2 *lambda[j]}
  
  # Run Filter
  res =  netfilt(t=t, Data=y, G=matrix(G_NIF[,1:Nbasis], ncol=Nbasis), sigmas=c(tau, sqrt(sig2_hat), alpha_hat), x0, var0, 1)
  x_kgk=res[[1]];  sigma_kgk = res[[2]]
  x_kgn = res[[3]]; sigma_kgn= res[[4]]
  
  c_NIF = matrix(0,Nbasis, len_t)  # matrix of stochastic coefficients
  Ts = matrix(0, Nbasis, 3*Nbasis)  # we can think it's "G" in KF 
  
  for(h in 1:len_t)
    for(i in 1:Nbasis)
    { Ts[i,(3*(i-1)+1):(3*(i-1)+3)] = c(t[h]-t[1],1,0) 
    c_NIF[,h] = Ts%*%(x_kgn[h,])
    }
  NIF_slip = matrix(t(Psi[1:Nbasis,]), ncol=Nbasis)%*%c_NIF  # estimated slips
  NIF_signal = G_NIF %*% c_NIF   # estimated data 
  time2 = Sys.time()
  
  record_time = as.numeric(difftime(time2, time1, units="secs"))
  
  res = as.list(NULL)
  res[["time"]] = record_time
  res[["est_slip"]] = NIF_slip
  res[["est_mean"]] = NIF_signal
  
  return(res)
}


NIF <- function(y, U_NIF, lambda, d, Phi, t, len_t, maxit=100){
  time1 = Sys.time()
  tau = 0 # set it to be 0 if we don't want to use a benchmark motion 
  G_NIF = U_NIF%*%diag(lambda)
  Nbasis = d
  G_NIF = matrix(G_NIF[,1:Nbasis], ncol=Nbasis)
  statedim = 3*Nbasis # dimension of latent states for Kalman Filter
  Psi = t(U_NIF)%*%Phi    # compute basis functions
  
  param_ini= c(10, 10)
  m=try(optim(param_ini, neg_loglik, t=t, Data=y, G=G_NIF, Lambda=lambda,
              lower = c(1e-6,1e-6), upper=c(1e+5, 1e+5),
              method="L-BFGS-B", control = list(maxit = maxit)),silent=T)
  # while(is.character(m)){
  #   param_ini=param_ini+runif(2)
  #   m=try(optim(param_ini,neg_loglik, t=ellipse_t, Data=ellipse_output, G=G_NIF, Lambda=lambda,lower = c(1e-6,1e-6),upper=c(1e+4,1e+4),
  #               method="L-BFGS-B", control = list(maxit = 100)),silent=T)
  # }
  theta_hat = m$par[1]
  theta2_hat = m$par[2]
  tau = 0
  
  res = netloglik3p(t, y, G_NIF, lambda, theta_hat, theta2_hat)
  sig2_hat = res[2] 
  alpha_hat = sqrt(sig2_hat)*theta_hat
  gamma_hat = sqrt(sig2_hat)*theta2_hat
  
  # construct the initial states
  x0 = matrix(0, statedim, 1)
  var0 = 1.0^(-10)*lambda[1]*diag(statedim)
  for(j in 1:Nbasis)
  {var0[(3*(j-1)+1),(3*(j-1)+1)] = gamma_hat^2 *lambda[j]}
  
  # Run Filter
  res =  netfilt(t=t, Data=y, G=matrix(G_NIF[,1:Nbasis], ncol=Nbasis), sigmas=c(tau, sqrt(sig2_hat), alpha_hat), x0, var0, 1)
  x_kgk=res[[1]];  sigma_kgk = res[[2]]
  x_kgn = res[[3]]; sigma_kgn= res[[4]]
  
  c_NIF = matrix(0,Nbasis, len_t)  # matrix of stochastic coefficients
  Ts = matrix(0, Nbasis, 3*Nbasis)  # we can think it's "G" in KF 
  
  for(h in 1:len_t)
    for(i in 1:Nbasis)
    { Ts[i,(3*(i-1)+1):(3*(i-1)+3)] = c(t[h]-t[1],1,0) 
    c_NIF[,h] = Ts%*%(x_kgn[h,])
    }
  NIF_slip = matrix(t(Psi[1:Nbasis,]), ncol=Nbasis)%*%c_NIF  # estimated slips
  NIF_signal = G_NIF %*% c_NIF   # estimated data 
  time2 = Sys.time()
  
  record_time = as.numeric(difftime(time2, time1, units="secs"))
  
  res = as.list(NULL)
  res[["time"]] = record_time
  res[["est_slip"]] = NIF_slip
  res[["est_mean"]] = NIF_signal
  
  return(res)
}


################################ Latent factors ##################################################
join_diff_para_loglik <- function(param,input, output, have_noise=TRUE, kernel_type='matern_5_2', lambda){
  #%% lambda: the vector of sigular values of \Phi
  res = 0
  sum_S2 = 0 
  # parameters we need to include: c(beta_1,...,beta_M, tau_1,...,tau_M, sigma_0)
  # tau_i = sigma_i^2/sigma_0^2
  # beta_i: range parameter for the ith latent process
  for(i in 1:M){
    param_here = c(param[i], -param[M+i]-2*log(lambda[i]) )
    output_here = output[i,]/sqrt(exp(param[M+i]))
    log_det_S2_here = Get_log_det_S2(param_here, have_noise=have_noise, delta_x, output_here,
                                     kernel_type=kernel_type)
    
    res = res - log_det_S2_here[[1]]/2 - n*param[M+i]/2
    sum_S2 = sum_S2 + log_det_S2_here[[2]]
  }
  sigma_0_2_hat = sum_S2/(n*M)
  -( res - (n*M)/2*log(sigma_0_2_hat) )
}

individual_loglik <- function(param,input, output, have_noise=TRUE, kernel_type='matern_5_2', lambda, sigma0=sigma_0){
  param_here = c(param[1], log(sigma0^2)-2*log(lambda)- param[2] )
  log_det_S2_here = Get_log_det_S2(param_here, have_noise=have_noise, delta_x, output,
                                   kernel_type=kernel_type)
  
  log_det_S2_here[[1]]/2 + n*param[2]/2 + log_det_S2_here[[2]]/(2*exp(param[2]))
}
