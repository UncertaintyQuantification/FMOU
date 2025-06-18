 
function Q = makeQ(statedim, Nbasis, Nsitecomps,  ...
			sigmas, delt, sig_frame)

%%  function Q = makeQ(statedim, Nbasis, post_step, ...
%%			Ncomps, sigmas, delt, sig_frame)
%%
%%  subroutine to make "Process" Covariance for Kalman filter
%%  Q is (statedim x statedim)
%%
%%  Input:
%%	statedim  = dimension of state vector
%%	Nbasis	  = number of basis functions
%%	alpha     = scale parameter for integrated random walk
%%	tau       = scale parameter for random walk
%%	delt      = t_k - t_{k-1} at present epoch
%%	
%%  Paul Segall, July 1997

%-(0). Setup null matrix

	Q = zeros(statedim);

	tau_GPS = sigmas(2);
	alpha = sigmas(3);

%-(1). For Tokai slip

	for i = 1:Nbasis
            Q(2*i-1,2*i-1) = (alpha^2)*(delt^3)/3;
            Q(2*i-1,2*i-0) = (alpha^2)*(delt^2)/2;
            Q(2*i-0,2*i-1) = (alpha^2)*(delt^2)/2;
            Q(2*i-0,2*i-0) = (alpha^2)*delt;
	end

%-(4). For random walk

	p = 2*Nbasis;
%for GPS only
	for i = 1:(Nsitecomps/3)
		Q(p+(3*i-2),p+(3*i-2)) = tau_GPS^2*delt;
		Q(p+(3*i-1),p+(3*i-1)) = tau_GPS^2*delt;
		Q(p+(3*i-0),p+(3*i-0)) = tau_GPS^2*delt;
    end

%-(5). For Reference frame

	Nframe = length(sig_frame);

    if Nframe ~= 0

	p = 2*Nbasis + Nsitecomps;

	for i = 1:Nframe
		Q(p+i,p+i) = sig_frame(i)^2;
	end

   end

%-(6). Form a process noise matrix

	Q = sparse(Q);


