    function F = makeF(statedim, Nbasis, Ncomps, delt, sig_frame);
%%
%%    function F = makeF(statedim, Ntype, Nbasis, ...
%%			Ncomps, delt, sig_frame)
%%
%%  subroutine to make state transition matrix for Kalman filter
%%  F is (statedim x statedim)
%%
%%  Input:
%%	statedim = dimension of state vector
%%	Nbasis	 = number of basis functions
%%	deltt    = t_k - t_{k-1} at present iteration
%%

%-(0). Identity for random walk

	F = eye(statedim);
	
%-(1). For Tokachi slip

	for i = 1:Nbasis
		F(2*i-1,2*i-1) = 1;
		F(2*i-1,2*i-0) = delt;
		F(2*i-0,2*i-1) = 0;
		F(2*i-0,2*i-0) = 1;
	end

%-(4). benchmark

	% EYE

%-(5). For reference frame shift

	Nframe = length(sig_frame);
	p = 2*Nbasis + Ncomps;
	for i = 1:Nframe
		F(p+i,p+i) = 0;
	end  

%-(5). Form a transition matrix

	F = sparse(F);

