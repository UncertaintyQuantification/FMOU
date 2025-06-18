function [x_f, covx_f] = backward(x_f, covx_f, Cascadia2, ...
				  G, delta_k, sig_frame, gparms)

  %-(1). Determine some dimensions ------------------------------

  Nepochs = size(Cascadia2,2);

  statedim = length(x_f{1,1});
  [Nsitecomps, Nbasis]  = size(G);

  %-(5). Running smoother

  %% Smoother 
  %%
  %%  x_k|N = x_k|k + S_k*(x_k+1|N - x_k+1|k)
  %%  C_k|N = C_k|k + S_k*(C_k+1|N - C_k+1|k)*S^T
  %%
  %%   where:
  %%   	S_k = C_k|k  *  F_k^T * inv(C_k+1|k)
  %%      
  %%   note that 
  %%		x_k+1|k = F_k * x_k|k 
  %%	& 	C_k+1|k = F_k * C_k|k * F_k^T  + Q_k

  disp('Running Smoother')

  %-(7.1). store the final state vector "x_n|n" and its covariance
  %        "C_n|n" into "x_kgn(n,:)" and "sigma_kgn(n,:)"

  kk_last = size(x_f);
  Mepochs = kk_last(1);

  %-(7.2). initialize state covariance matrix

  S_k = zeros(statedim,statedim);

  %-(7.3). backward smoothing

  k_last = delta_k*(kk_last-1) + 1;

  for k = (k_last(1)-delta_k):-delta_k:1
    kk = (k-1)/delta_k + 1;
    fprintf(1,'[k kk] = %d %d\n',k,kk);

    %-(7.3.1). define the epoch interval
    delt = abs( Cascadia2(k+1).year - Cascadia2(k).year );

    %-(7.3.2). construct the state transition matrix "F"
    F = makeF(statedim, Nbasis, Nsitecomps, delt, sig_frame);

    %-(7.3.4). smoothed state vector "x_k|n" and covariance "C_k|n"

    % AMB: I could use kf_qrsc_smooth here. This implementation is not a
    % square-root smoother and so can run into problems due to subtracting one
    % cov matrix from another. However, I haven't figured out if I can make
    % kf_qrsc_smooth faster, and since this code has been working so far and is
    % faster, I think I'll leave it for now.
    
    % sigkp1gk = zz'*zz
    zz = chol(covx_f(kk+1,kk));
    %% Smoothed state vector: x_k|N = x_k|k + S_k*(x_k+1|N - x_k+1|k)
    y = x_f{kk+1,Mepochs} - x_f{kk+1,kk};
    x_f{kk,Mepochs} = x_f{kk,kk} +...
	covx_f(kk,kk)*(F'*(zz \ (y' / zz)'));
    %% Smoothed covariance: C_k|N = C_k|k + S_k*(C_k+1|N - C_k+1|k)*S_k^T
    % AMB: This is where problems can arise: C is theoretically p.d. but may
    % be indefinite in finite precision arithmetic.
    C = covx_f{kk+1,Mepochs} - covx_f{kk+1,kk};
    BFt = covx_f(kk,kk)*F';
    covx_f{kk,Mepochs} = covx_f{kk,kk} +...
	BFt*(((zz \ (zz' \ C)) / zz) / zz')*BFt';
    
    if (gparms.proj_use_bwd)
      lb = [-inf*ones(1,Nbasis); zeros(1,Nbasis)];
      one = ones(length(x_f{kk,Mepochs}) - 2*Nbasis,1);
      lb = [lb(:); -inf*one];
      if (gparms.proj_weak)
	ub = inf*ones(1,2*Nbasis);
      else
	ub = [x_f{kk+1,Mepochs}(1:2:2*Nbasis)'; inf*ones(1,Nbasis)];
      end
      ub = [ub(:); inf*one];
      
      x_f{kk,Mepochs} = project(...
	  gparms.popts, x_f{kk,Mepochs}, covx_f{kk,Mepochs}, lb, ub);

    end
    
      
     if(mod(k, 50)==0)
        save backward_save *
    end
  

    %-(5.4.7). monitoring of fault slip rate
    
    mean_slip = mean( x_f{kk,Mepochs}(1:2:2*Nbasis-1) );
    mean_rate = mean( x_f{kk,Mepochs}(2:2:2*Nbasis) );
    
    fprintf(1,['[average slip [cm], average slip-rate [cm/yr]] = ',...
	       '%s %s\n'],num2str(100*mean_slip),num2str(100*mean_rate));
  end
