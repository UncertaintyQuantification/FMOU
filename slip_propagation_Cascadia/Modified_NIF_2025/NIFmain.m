function NIFmain(gp)
% gp is a structure of filter setup info and parameters. Modify gparms.m to
% produce gp.
  
%%
%% Script to analyze Loma Prieta post seismic data
%% Including spatial variation in slip
%% Paul Segall July-August 1997
%%
%%  modifed for spatial smoothing at each step June 1998
%%  modified 10/15/98 to make first basisfunction uniform slip

tic; 
%-(1). clear all variables -------------------------------------
%      and loading the setup files -----------------------------

	
	% Use AMB changes? If so, set amb=1; else, amb=0. (This exists for
        % debugging only; nominally, amb=1 always.)
        %NEVER CHANGE THIS
	global amb
	amb = 1;
	
	% Save-file stuff
	% Form output mat name
	svfn = gp.svfn;
% 	% Error if it exists
% if (~isempty(dir([svfn '.mat'])))
% 	  error('File exists, my friend!');
%end
	% Write comment file
	fid = fopen([svfn '.cmt'],'w');
	flds = fieldnames(gp);
	for i = 1:length(flds)
	  v = gp.(flds{i});
	  if (~ischar(v)) v = num2str(v); end
	  fprintf(fid,'%15s = %s\n',flds{i},v);
	end
	fclose(fid);
	
	% Load setup file 
%	matpathcas
	format long
	el = [];
	load(gp.setup_mat);
    % rrake = (180-rake)/180*pi; 
     
     gp.popts = struct('Wid',gp.proj_Wid,'weak',gp.proj_weak);
	
%-(2). GPS time series -----------------------------------------

	% load stationstruct_data.mat
	load(gp.data_mat);
    %removed for resolution testing
    k_init=1;
    %k_init=2040; %starting epoch
    %time=50; %length of event
    %time=54;
	%stationstruct = stationstruct(k_init:k_init+time);
    

	Nepochs = size(stationstruct,2);
	for k = 1:Nepochs
		t(k) = stationstruct(k).year;
        %adjust GPS white noise level
        %stationstruct(k).data(:, 4:6)=ones(length(stationstruct(k).data), 1)*[6, 6, 18];
    end
    
    %change from years to hours
    %t=(t-2010)*365.24*24;
    
    

    
	delta_k = 1;

%-(3). state dimension -----------------------------------------
%
%      2*Nbasis1  =  Nbasis1 for slip  +  Nbasis1 for slip-rate
%      Ncomps     =  random walk
%      Nframe     =  reference frame term

	statedim = 2*Nbasis1 + Ncomps + gp.Nframe;



 
%minimum norm model instead of Laplacian smoothing
if(gp.minnorm==1)
    Lap=eye(size(Lap));
end
    
    %-(5). Set hyper-parameters
%      Those values are derived from Maximum Likelihood w/ full data set.

	%sigma_hat = 0.75;
	sigma_hat = gp.sigma_hat;
	tau_GPS_hat = gp.tau_GPS_hat;

	%alpha_hat = 8.0; 
	alpha_hat = gp.alpha_hat; 
	gamma_hat = gp.gamma_hat;


	%sig_frame = sigma_hat*[5, 5, 5]/1000;
    if(gp.Nframe>0)
        sig_frame = 5*ones(gp.Nframe, 1)/1000;
    else
        sig_frame=[];
    end

%-(6). initial state vector 'x0' ------------------------------------

	x0 = zeros(statedim,1);


%-(7). initial state covariance matrix 'var0' ----------------------------

	hundm = 1.0e+4;  tenm = 1.0e+2;   onem = 1.0e+0;
	tencm = 1.0e-2;  onecm = 1.0e-4;  onemm = 1.0e-6;
	submm = 1.0e-8;

	var0 = onemm*eye(statedim);

% 1[mm] and 1[cm/yr] initial constraint for transient fault slip

	for i = 1:Nbasis1
		var0(2*i-1,2*i-1) = gp.var0_slip;
		var0(2*i-0,2*i-0) = gp.var0_rate;
	end

% 1[cm] apriori constraint for site coordinates

	p = 2*Nbasis1;

	
	for i = 1:Ncomps
		var0(p+i,p+i) = gp.var0_coord;
	end

% for reference frame

	p = 2*Nbasis1 + Ncomps;

	for i = 1:gp.Nframe
		var0(p+i,p+i) = sig_frame(i)^2; 
	end


%-(6). run network inversion filter ------------------------------

%-(6.1).  SCALING by sigma_hat

        C0g0 = (sigma_hat^2)*var0;

        sigmas = sigma_hat*[1.0, tau_GPS_hat, alpha_hat, gamma_hat];
        %first comp of sigmas weights GPS errors
        
        
	%ID_T = [];
	%Kern_TILT = [];
	%sig_T = [];

	range = 90; % This is only for non-negative.
    %%
    

%-(6.2). forward filtering -------------------

w=gp.w;
staweight=gp.staweight;
ID=ID_G;

%save before forward
save([svfn '.mat']);

%check for mle
%changed for restricted trangles
if(gp.mle==1)
%for i=length(gp.alpha):-1:1
%    for j=length(gp.gamma):-1:1
for i=length(gp.alpha):-1:1
    for j=length(gp.gamma):-1:1
        
                
        disp(i)
        disp(j)
    
        sigmas = [sigma_hat, tau_GPS_hat sigma_hat*gp.alpha(i), sigma_hat*gp.gamma(j)]; 
        
               
   	svfn = sprintf('%s%sf%db%dW%d_%d%d',...
		       gp.save_dir,gp.savefile,gp.proj_use_fwd,gp.proj_use_bwd,...
		       gp.proj_Wid, i, j);

	[fit, minus2_logL, sig2, rsum, vsum, N_all, ...
		dchi, x_f, covx_f, t, dsave, dsave2] ...
	    = forward_nnls(x0, C0g0, stationstruct, Kern_GPS, ID, ...
		delta_k, sigmas, sig_frame, Lap, ...
		gp.smooth_mode, gp, staweight, w);

     m2LL(i, j)=minus2_logL;
     sig2_M(i, j)=sig2;
     
    
	%save([svfn '.mat']);
    end
end


else
    [fit, minus2_logL, sig2, rsum, vsum, N_all, ...
		dchi, x_f, covx_f, t, dsave, dsave2] ...
	    = forward_nnls(x0, C0g0, stationstruct, Kern_GPS, ID, ...
		delta_k, sigmas, sig_frame, Lap, ...
		gp.smooth_mode, gp, staweight, w);
    
   % save([svfn '.mat']);

end
    %%
    
    %find max likelihood
    %Changed for restricted traingles
    if(gp.mle==1)
        [ihat, jhat] = find(m2LL==min(min(m2LL)));
        alpha_hat=gp.alpha(ihat);
        gamma_hat=gp.gamma(jhat);
        
        sigma_hat=sqrt(sig2_M(ihat, jhat));
         sigmas = [sigma_hat, tau_GPS_hat, sigma_hat*gp.alpha_hat, sigma_hat*gp.gamma_hat];
        
%          [fit, minus2_logL, sig2, rsum, vsum, N_all, ...
%  		dchi, x_f, covx_f, t, dsave, dsave2] ...
%  	    = forward_nnls(x0, C0g0, stationstruct, Kern_GPS, ID, ...
%  		delta_k, sigmas, sig_frame, Lap, ...
%  		rrake, range, gp.smooth_mode, origin, gp, k_init, staweight, w);

        [fit, minus2_logL, sig2, rsum, vsum, N_all, ...
 		dchi, x_f, covx_f, t, dsave, dsave2] ...
 	    = forward_nnls(x0, C0g0, stationstruct, Kern_GPS, ID, ...
		delta_k, sigmas, sig_frame, Lap, ...
		gp.smooth_mode, gp, staweight, w);

        svfn = sprintf('%s%sf%db%dW%d_fin',...
		       gp.save_dir,gp.savefile,gp.proj_use_fwd,gp.proj_use_bwd,...
		       gp.proj_Wid);
    
    end
    
    
    %%

%-(6.3). backward smoothing -------------------

	[x_f, covx_f] = backward(x_f, covx_f, stationstruct, ...
		Kern_GPS, delta_k, sig_frame, gp);

	%save([svfn '.mat']);
    %%

%-(7). retrieving slip and slip-rate ------------------
disp('start get rate')
	[slip_f, rate_f, slip_b, rate_b, ...
		sig_slip_f, sig_rate_f, sig_slip_b, sig_rate_b, ...
		Cslip_b, sig_Cslip_b] ...
	= get_rate_slip(x_f, covx_f, Kern_GPS, Nepochs);
disp('done get rate')
save([svfn '.mat']);
%-(8). plot the slip ------------------
%%

	cum_slip = slip_f(:,end);

	figure
	trisurf(el,nd_ll(:,1),nd_ll(:,2),nd_ll(:,3),cum_slip);
	colorbar;
	view(2);
	daspect([1,1,1]);
	xlabel('Easting'); ylabel('Northing'); daspect([1,1,1]);
	hold on
	plot(lon_G, lat_G, 'kx', 'LineWidth', 2)
	%print -dpsc2 finalslip.ps
    %%

%-(9). estimates of coordinates with random walk ------------------
%  
%       "L" starts from 0 at the initial epoch
%       "pos_adj" includes offsets at the initial epoch

        [pos_adj, pos_sigma] = pos_sig(x_f, covx_f, ...
                stationstruct, delta_k, Nbasis1, Ncomps);

 	pos_0 = pos_adj(:,1)*ones(1,Nepochs); 
 	L = pos_adj - pos_0; 
    %%
    

%-(10). estimates of transformation parameters
    if(gp.Nframe>0)
	[TT, Trans] ...
		= get_frame(x_f, Nepochs, Ncomps/3, Nbasis1, Ncomps, gp.Nframe);
    end

svfn
%-(11). compute predicted data in ENU frame
%-       Random-walk component is included --------------------------

        % dsave = original GPS data
	% pos_0 = initial position of each station component
	% Trans = common mode (ref frame error)
    if(gp.Nframe>0)
        data = dsave - pos_0;
    else
        data = dsave - pos_0;
    end

        [E, N, U, Ep, Np, Up] = OP_GPS(t, data, Cslip_b);


        Zr = zeros(size(E));
        E_Sig = Zr; N_Sig = Zr; U_Sig = Zr; NE_Corr = Zr;

        %ZrT = zeros(size(ET));
        %ET_Sig = ZrT; NT_Sig = ZrT; NET_Corr = ZrT;

	save([svfn '.mat']);
%% save the running time
    total_time = toc;
    save([svfn '.mat']);


   


