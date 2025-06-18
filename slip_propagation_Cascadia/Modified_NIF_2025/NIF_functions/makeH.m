	function H = makeH(statedim, IND, G)


%	function [d_t, cov_t, H] = ...
%	    makedandH(x, sites, ENU, covENU, ...
%		tk, G, ini_epoch, origin)
%
% input: 
%		sites	    = list of all stations
%		pos	    = positions observed at this epoch
%		pos_cov	    = position covariance at this epoch
%		epoch	    = this epoch
%		G	    = matrix relating displacement to slip
% output: 
%		d_t 	    = data vector for this epoch
%		cov_t	    = data covariance for this epoch
%		H	    = matrix relating data to state vector at
%				curent epoch
% dimension of state vector:
%               position    = m
%               velocoty    = mm/yr
%               translation = m
%               rotation    = mas
%               scale       = 10^(-8)



	[Nsitecomps, Nbasis] = size(G);
	Nsites = Nsitecomps/3;

	Nobs = length(IND);
	Nobscomps = 3*Nobs;

% for each station find matching station and extract 
% appropriate part of G into G_t. also form matrix Ir

	G_t = zeros(Nobscomps,Nbasis);
	Ir = zeros(Nobscomps,Nsitecomps);
    %include if frame only
    if((Nbasis*4+Nsitecomps)<statedim)
	Rtr = zeros(Nobscomps,statedim-Nbasis*2-Nsitecomps);
    end
	% subR = rotmatrix(origin(1), origin(2));


	for i = 1:Nobs
		Ir(3*i-2:3*i, 3*IND(i)-2:3*IND(i)) = eye(3);
		%G_t = [G_t; G(3*IND(i)-2:3*IND(i),:)];
		G_t(3*i-2:3*i,:) = G(3*IND(i)-2:3*IND(i),:);
		% Rtr = [Rtr; subR];
        if((Nbasis*2+Nsitecomps)<statedim)
		Rtr(3*i-2:3*i,:) = eye(3);
        end
	end


%% construct the submatrix [1,0]

	sub = zeros(Nbasis, 2*Nbasis);
	sub(:,1:2:2*Nbasis-1) = eye(Nbasis);

%% Form H_d
if((Nbasis*2+Nsitecomps)<statedim)
	H = [G_t*sub, Ir, Rtr];
else
    H = [G_t*sub, Ir];
end

