function [fit, val, sig_d2, rsum, vsum, N_all, ...
	  dchi, x_f, covx_f, t, dsave, dsave2] = ...
    forward(x0, var0, stationstruct, G_G, sites, ...
	    delta_k, sigmas, sig_frame, Lap, ...
	     smooth_mode,  gparms, ...
         staweight, w)
G=G_G;
 
 
  %-(1). Determine some dimensions ------------------------------

  Ntype = 2; % IRW for slip

  Nepochs = size(stationstruct,2);
  [Nsitecomps, Nbasis]  = size(G);
  Nsites = Nsitecomps/3;
  Npseudo = size(Lap,1);
  statedim = length(x0);

  dsave  = NaN*ones(Nsitecomps, Nepochs);
  dsave2 = NaN*ones(Nsitecomps, Nepochs);

  
  %-(2). Variance parameters ------------------------------------

  sig_d = sigmas(1);
  gamma_j = sigmas(4);

  %-(3). Matricies to be used in smoothing operation
  %      & for Square Root Filter -------------------------------

  %-(3.2). initialize variables

  x_f = state([zeros(statedim,1)]);
  covx_f = state(zeros(statedim,statedim), gparms.scratch_dir);
  x_f{1,0} = x0;
  covx_f{1,0} = var0;

  dchi = zeros(Nepochs);
  s = zeros(statedim+Nsitecomps);

  %-(4). x(0|0) and C(0|0) ------------------------

  N_all = 0; vsum = 0; rsum = 0;

  %-(5). Run Square Root Filter ----------------------------------
  %%
  %% Square Root filter for estimating network model
  %% Data is:     d(t) = G*s(t) + B(t) + white noise
  %% slip is: 	s(t) = b*t + W(t)

  %added for traingles near tremor restriction
  Lap_orig=Lap;
  
  for k = 1:delta_k:Nepochs
    tic_t = tic;
    kk = (k-1)/delta_k + 1;
    t(k) = stationstruct(k).year;
    t2(k)=t(k);
    
    if(k~=1)
        delta_t=t(k)-stationstruct(k-1).year;
    else
        delta_t=t(k);
    end
    
    fprintf(1,'[The epoch and the time] = %d %f\n',k,t(k));

    %-(5.1). Reading Data

    % set up data vector - pos (apriori coords);
    clear d_org ENU covENU IND

    sitecodes = stationstruct(k).names;
    % data and its uncertainities are converted from [mm] to [m].
    d_org = (stationstruct(k).data(:,1:3))'/1000;

    d_org = d_org(:);
    N_d = length(d_org);
    Nsta = N_d/3;
    ENU = NaN*ones(N_d,1);
    

    sigENU = (stationstruct(k).data(:,4:6))'/1000;
    sigENU = (sig_d^2)*sigENU(:).^2;
    %weighting matrix for stations to weight, weighting factor w (w=1 for
    %no weighting)
    W=eye(length(sigENU));
    for h=1:size(staweight, 1)
        hind=GetIndex(sitecodes, staweight(h, :));
        if(hind>0)
        W(3*hind-2:3*hind, 3*hind-2:3*hind)=w*W(3*hind-2:3*hind, 3*hind-2:3*hind);
        end
    end
    

sigENU=W*sigENU;
    
    
    
    covENU = sparse(diag(sigENU));

    sitecodes = stationstruct(k).names;
    IND = [];

    % secular velocities

    % contains Ev, Nv and Uv
    load(gparms.secular_mat); 
    % unit in [m/yr]
    tmp = [Ev, Nv, Uv; Ea1, Na1, Ua1; Ea2, Na2, Ua2; ...
	   Es1, Ns1, Us1; Es2, Ns2, Us2]'; 

    X_init = tmp(:);
    %X_init=[X_init; X_f(:, kk+k_init-1)];
    clear tmp

    for i = 1:Nsta   % loop over this epoch's stations

j = GetIndex(sites, sitecodes(i,:));
      IND = [IND; j];

      %THIS ASSUMES t in YEARS
     
      
      H_init = zeros(3, length(X_init));
      H_init(:, 3*j-2:3*j) = (t2(k)-t2(1))*eye(3);
      % Secular is [E,N,U, E,N,U, E,N,U, ...]

      %	H_init(3*i-2:3*i, 3*j-2:3*j) = t(k)*eye(3);
      %	    % Secular is [E,N,U, E,N,U, E,N,U, ...]
      p = Nsitecomps;
      H_init(:, p+3*j-2:p+3*j) ...
	  = sin(2*pi*t2(k))*eye(3);
      % Annual (sin)
      p = 2*Nsitecomps;
      H_init(:, p+3*j-2:p+3*j) ...
	  = cos(2*pi*t2(k))*eye(3);
      % Annual (cos)
      p = 3*Nsitecomps;
      H_init(:, p+3*j-2:p+3*j) ...
	  = sin(4*pi*t2(k))*eye(3);
      % Semi-annual (sin)
      p = 4*Nsitecomps;
      H_init(:, p+3*j-2:p+3*j) ...
	  = cos(4*pi*t2(k))*eye(3);
      % Semi-annual (cos)
%       if(length(sig_frame)>0)
%       	p = 5*Nsitecomps;
%       	H_init(:, p+3*j-2:p+3*j) = eye(3);
%       % Wobble is [E,N,U, E,N,U, E,N,U, ...]
%       	H_init(3*i-2:3*i, end-2:end) = zeros(3);
%       	    % frame 
%       end
    %data minus secular, annual, semiannual
    ENU(3*i-2:3*i) = d_org(3*i-2:3*i) - H_init*X_init;
      dsave(3*j-2:3*j,k) = ENU(3*i-2:3*i);
      dsave2(3*j-2:3*j,k) = d_org(3*i-2:3*i);
    fit(3*j-2:3*j, k)=H_init*X_init;
    end
    

    %-(5.2). Prediction step
    if (k ~= 1)
      delt = abs( t(k) - t(k-1) );

      %-(5.2.2). construct state transition matrix "F" and
      %          process noise matrix "Q"
      F = makeF(statedim, Nbasis, Nsitecomps, delt, sig_frame);
      Q = makeQ(statedim, Nbasis, Nsitecomps, sigmas, delt, sig_frame);

      %-(5.2.3). predicted state x(k+1|k) and its covariance C(k+1|k)
      x_f{kk,kk-1} = F*x_f{kk-1,kk-1};
      covx_f{kk,kk-1} = F*covx_f{kk-1,kk-1}*F' + Q;
    end 
    
    % Updates
    
    %-(5.4). PSEUDO-UPDATE by pseudo measurements
    %-(5.4.1). create data vector "d_t", its covariance "cov_t",
    %          and design matrix "H"
    
    %remove on tilt update steps ONLY, look for epochs of only 1 station
 %   if(length(IND)>1)
    
    [d_p, H_p] = smooth_slip_rate2(...
	statedim, Lap, Nbasis, smooth_mode, Npseudo, delta_t);
    N_p = length(d_p);
    diagR_p = (gamma_j^2)*ones(N_p,1);
    
    Rc = spdiags(sqrt(diagR_p),0,N_p,N_p); % because R_p is diag
    
 %   else
        
 %       H_p=[];
 %       Rc=[];
 %       d_p=[];
 %       N_p=0;
        
%    end
    
    %-(5.4). UPDATE by REAL DATA
    d_d = ENU;
    R_d = covENU;
    H_d = makeH(statedim, IND, G);
    


    % Decide whether to update in batch or sequentially. They are
    % mathematically equivalent; the only issue is computation time (and
    % memory, but I'm not sure how to assess the effect of that in a simple
    % way). - AMB
    if (N_p + N_d >= length(x_f{kk,kk-1}))
      % pseudo
      [x_pp sqrsigma_pp] = kf_qrsc_update(...
	  H_p,Rc,d_p,x_f{kk,kk-1},chol(full(covx_f{kk,kk-1})));    
      % real
      [x_hat sqrsigma_f nu_d V_dc] = kf_qrsc_update(...
	  H_d,chol(R_d),d_d,x_pp,sqrsigma_pp); 
      % Log likelihood of just the real data
      V_d_bs_nu_d = (V_dc \ (nu_d' / V_dc)');
      vsum_k = 2*sum(log(diag(V_dc)));
      rsum_k = nu_d'*V_d_bs_nu_d;
    else 
      % both
      [x_hat sqrsigma_f nu_d V_dc] = kf_qrsc_update(...
	  [H_p; H_d], blkdiag(Rc,chol(R_d)), [d_p; d_d],...
	  x_f{kk,kk-1}, chol(full(covx_f{kk,kk-1})));
      % Log likelihood of just the real data. This is
      %   p(yr_k | Y_k-1, yp_k) = p(y_k | Y_k-1) / p(yp_k | Y_k-1),
      % where yr is real data, yp is pseudo data, and y is both. The division
      % implies the following operations:
      rsum_k = dot(nu_d,(V_dc \ (nu_d' / V_dc)')) -...
	       dot(nu_d(1:N_p),...
		   (V_dc(1:N_p,1:N_p) \ (nu_d(1:N_p)' / V_dc(1:N_p,1:N_p))'));
      vsum_k = 2*sum(log(diag(V_dc(N_p+1:end,N_p+1:end))));
    end
    covx_f{k,k} = sqrsigma_f' * sqrsigma_f;
    fprintf(1,'rsum_k = %e  vsum_k = %e\n',rsum_k,vsum_k);
    
    %-(5.4.7). non-negative

    fprintf(1,'min, max slip = %s, %s\n',...
	    num2str(min(x_hat(1:2:2*Nbasis))),...
	    num2str(max(x_hat(1:2:2*Nbasis))));
    fprintf(1,'min, max rate = %s, %s\n',...
	    num2str(min(x_hat(2:2:2*Nbasis))),...
	    num2str(max(x_hat(2:2:2*Nbasis))));

    if (gparms.proj_use_fwd)
      if (gparms.proj_weak)
	lb = [-inf*ones(1,Nbasis); zeros(1,Nbasis)];
      else
	lb = [x_f{kk-1,kk-1}(1:2:2*Nbasis)'; zeros(1,Nbasis)];
      end
      lb = [lb(:); -inf*ones(length(x_hat)-2*Nbasis,1)];
      ub = inf*ones(size(lb));
      
      x_f{kk,kk} = project(...
	  gparms.popts, x_hat, covx_f{kk,kk}, lb, ub);
    else
      x_f{kk,kk} = x_hat;
    end

    fprintf(1,'min, max slip = %s, %s\n',...
	    num2str(min(x_hat(1:2:2*Nbasis))),...
	    num2str(max(x_hat(1:2:2*Nbasis))));
    fprintf(1,'min, max rate = %s, %s\n',...
	    num2str(min(x_hat(2:2:2*Nbasis))),...
	    num2str(max(x_hat(2:2:2*Nbasis))));

    %-(5.4.7). monitoring of fault slip rate
    mean_slip = mean( x_f{kk,kk}(1:2:2*Nbasis-1) );
    mean_rate = mean( x_f{kk,kk}(2:2:2*Nbasis) );

    fprintf(1,'[average slip [cm], average slip-rate [cm/yr]] = %s %s\n',...
	    num2str(100*mean_slip),num2str(100*mean_rate));

    %-(5.6). Components of Likelihood
    N_all = N_all + N_d;

    rsum = rsum + rsum_k;
    vsum = vsum + vsum_k;

    dchi(k) = sqrt(rsum_k/N_d);

    fprintf(1,'[chi-square value at each epoch] = %s\n',num2str(dchi(k)));
    fprintf(1,'loop wall clock time = %f\n',toc(tic_t));
    
    if(mod(k, 50)==0)
        save forward_save *
    end
    
  end

  %-(6). MLE value and its index
  val = N_all*log(rsum) + vsum;
  sig_d2 = rsum/N_all;
