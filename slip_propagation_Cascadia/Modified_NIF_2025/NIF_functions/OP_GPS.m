function [E, N, U, Ep, Np, Up] = OP_GPS(t, Data, Dhat)

%% [RData, E, N, U, RHat, Ep, Np, Up, kk] = OPnew(t, Data, Dhat, X0_enu)
%% 
%% Function to pull out Observed "O" and predicted "P" data for plotting.
%% Used with plotOP

%% INPUT:
%%	t 	= vector of observation times
%%	Data	= matrix of observed data
%%	Dhat	= matrix of predicted data
%%	XO_enu	= Apriori coordinates in ENU coordinates
%%	sites   = vector of site names
%%	rel_site= fixed station to which all baselines are formed
%% OUTPUT:
%%	RData 	= Matrix of relative observed data
%%	E,N,U	= East observed, N...
%%	RHat	= Matrix of relative predicted data
%%	Ep,Np,Up= East predicted, N...
%%	kk	= index of rel_site

        [Nsitecomps,Nepochs] = size(Data);
	Nsites = Nsitecomps/3;

% pull out the east north and up components for plotting

for i = 1:Nsites
	E(i ,:) = Data(3*i-2,:);
	N(i ,:) = Data(3*i-1,:);
	U(i ,:) = Data(3*i,:);

	Ep(i ,:) = Dhat(3*i-2,:);
	Np(i ,:) = Dhat(3*i-1,:);
	Up(i ,:) = Dhat(3*i,:);
end



